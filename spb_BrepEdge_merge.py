"""
Alternative to _MergeEdge, this script:
    1. Calls code that will remove unneeded knots in curves used by BrepTrims.
    2. RebuildEdges with the new curves (#1).
    3. All contiguous edges that pass BrepEdgeList.MergeEdge threshold are merged.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221026: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Input as ri
import scriptcontext as sc

import spb_BrepTrim_removeInteriorKnots


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bDecimalTol'; keys.append(key)
    values[key] = False
    names[key] = 'ToleranceType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'AbsValueOfSrfParamSpace', 'DecimalOfTrimCrvDomain')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_Absolute'; keys.append(key)
    values[key] = 1e-9
    names[key] = 'Tol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol_Decimal'; keys.append(key)
    values[key] = 1e-9
    names[key] = 'Decimal'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)


    for key in keys:
        if key not in names:
            names[key] = key[1:]


    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fTol_Absolute':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fTol_Decimal':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get edges with optional input.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select edge to merge")

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.EdgeCurve

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('bDecimalTol')
        if Opts.values['bDecimalTol']:
            addOption('fTol_Decimal')
        else:
            addOption('fTol_Absolute')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()


        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            rdObjRef_Edge = go.Object(0)
            go.Dispose()
            return rdObjRef_Edge

        if res == ri.GetResult.Number:
            if Opts.values['bDecimalTol']:
                key = 'fTol_Decimal'
            else:
                key = 'fTol_Absolute'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def processBrep(rgBrep, idx_Edge, bDecimalTol=None, fTol=None, bEcho=True, bDebug=False):
    """
    rgBrep: Input and modified.
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if bDecimalTol is None: bDecimalTol = Opts.values['bDecimalTol']
    if fTol is None:
        if bDecimalTol:
            fTol = getOpt('fTol_Decimal')
        else:
            fTol = getOpt('fTol_Absolute')


    iTrims_Modified = []
    iFaces_Modified = []

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt

    iCt_Edges_Pre = rgBrep.Edges.Count

    iCt_Merged = rgBrep.Edges.MergeEdge(
        edgeIndex=idx_Edge,
        angleTolerance=sc.doc.ModelAngleToleranceRadians)

    if iCt_Merged == 0: return

    iCt_Edges_Post = rgBrep.Edges.Count

    #for edge in rgBrep.Edges:
    #    print(edge.EdgeIndex)

    #sc.doc.Objects.AddCurve(edge); sc.doc.Views.Redraw(); return

    return spb_BrepTrim_removeInteriorKnots.processBrep(
        rgBrep,
        idx_Edges=[rgBrep.Edges.Count-1],
        bDecimalTol=bDecimalTol,
        fTol = fTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )


def processBrepObject(rhEdge_In, bDecimalTol=None, fTol=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rhEdge_In: rd.ObjRef of BrepEdge
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if bDecimalTol is None: bDecimalTol = Opts.values['bDecimalTol']
    if fTol is None:
        if bDecimalTol:
            fTol = getOpt('fTol_Decimal')
        else:
            fTol = getOpt('fTol_Absolute')

    rdB = rhEdge_In.Object()
    edge = rhEdge_In.Edge()
    idx_Edge = edge.EdgeIndex

    rgB_In = rdB.Geometry


    bModified = processBrep(
        rgBrep=rdB.Geometry,
        idx_Edge=idx_Edge,
        bDecimalTol=bDecimalTol,
        fTol = fTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if not bModified:
        print("No knots were removed.  Brep was not modified.")
        return False

    if rdB.CommitChanges():
        if bEcho:
            print("Brep was modified.")
    else:
        print("Brep could not be modified.")
        return False

    return True


def main():

    objrefs_Edges_In = getInput()
    if objrefs_Edges_In is None: return

    bDecimalTol = Opts.values['bDecimalTol']
    if bDecimalTol:
        fTol = Opts.values['fTol_Decimal']
    else:
        fTol = Opts.values['fTol_Absolute']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    processBrepObject(
        objrefs_Edges_In,
        bDecimalTol=bDecimalTol,
        fTol=fTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()