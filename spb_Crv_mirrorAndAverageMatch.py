"""
Behavior notes:
    Only wire curves are modified.
    Symmetry planes are per World, not CPlane.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220727: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fDistTol'; keys.append(key)
    values[key] = 1.0
    names[key] = 'MaxDist'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iSymPlane'; keys.append(key)
    values[key] = 0
    listValues[key] = 'X', 'Y', 'Z'
    stickyKeys[key] = '{}({})'.format(key, __file__)

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
                # For OptionList.
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fDistTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            if cls.riOpts[key].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = cls.values[key] = Rhino.RhinoMath.ZeroTolerance
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
        else:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            elif key in cls.listValues:
                cls.values[key] = idxList
            else:
                return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def _getInput_CurveObjects():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select wire curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fDistTol')
        #addOption('bEcho')
        #addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fDistTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_SymmetryPlane():
    """
    Get option.
    """

    idxOpt_Default = Opts.values['iSymPlane']

    go = ri.Custom.GetOption()

    go.SetCommandPrompt("Pick mirror plane")
    go.SetCommandPromptDefault(Opts.listValues['iSymPlane'][idxOpt_Default])
    go.AcceptNothing(True)

    idxs_Opt = {}

    for sOpt in Opts.listValues['iSymPlane']:
        idxs_Opt[sOpt] = go.AddOption(sOpt)

    res = go.Get()


    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Nothing:
        idxOpt = idxOpt_Default
    elif res == ri.GetResult.Option:
        idxOpt = go.Option().Index - 1 # '- 1' because go.Option().Index is base 1.
        Opts.setValue('iSymPlane', idxOpt)

    go.Dispose()

    return Opts.listValues['iSymPlane'][idxOpt]


def createSymmetricPair(nc, plane_sym, bT1):
    
    nc_Out_A = nc.DuplicateCurve()
    nc_Out_B = nc.DuplicateCurve()
    
    xform = rg.Transform.Mirror(plane_sym)
    
    if not nc_Out_B.Transform(xform):
        print("Transformation failed.")
        nc_Out_A.Dispose()
        nc_Out_B.Dispose()
        return
    
    return rg.Curve.CreateMatchCurve(
        curve0=nc_Out_A,
        reverse0=bT1,
        continuity=rg.BlendContinuity.Tangency,
        curve1=nc_Out_B,
        reverse1=bT1,
        preserve=rg.PreserveEnd.Tangency,
        average=True)


def main():
    
    objrefs_In = _getInput_CurveObjects()
    if objrefs_In is None: return
    
    sSymPlane = _getInput_SymmetryPlane()
    if sSymPlane is None: return
    
    fDistTol = Opts.values['fDistTol']
    
    iCt_Results = 0
    
    for objref_In in objrefs_In:
        rgC_In = objref_In.Curve()
        if not isinstance(rgC_In, rg.NurbsCurve): continue
        
        nc_In = rgC_In
        
        if sSymPlane == 'X':
            plane_sym = rg.Plane.WorldZX
            fStartValue = nc_In.PointAtStart.Y
            fEndValue = nc_In.PointAtEnd.Y
        elif sSymPlane == 'Y':
            plane_sym = rg.Plane.WorldYZ
            fStartValue = nc_In.PointAtStart.X
            fEndValue = nc_In.PointAtEnd.X
        elif sSymPlane == 'Z':
            plane_sym = rg.Plane.WorldXY
            fStartValue = nc_In.PointAtStart.Z
            fEndValue = nc_In.PointAtEnd.Z
        else:
            raise ValueError("{} not a value symmetry plane description.".format(
                sSymPlane))

        if abs(fStartValue) <= fDistTol:
            ncs_WIP = createSymmetricPair(nc_In, plane_sym, False)
            if not ncs_WIP: continue
            if abs(fEndValue) <= fDistTol:
                ncs_WIP[1].Dispose()
                ncs_WIP = createSymmetricPair(ncs_WIP[0], plane_sym, True)
                if not ncs_WIP: continue
        elif abs(fEndValue) <= fDistTol:
            ncs_WIP = createSymmetricPair(nc_In, plane_sym, True)
        else:
            continue

        sc.doc.Objects.Replace(objref_In, ncs_WIP[0])

        sc.doc.Objects.AddCurve(ncs_WIP[1], objref_In.Object().Attributes)

        iCt_Results += 1


    sc.doc.Views.Redraw()


    print("{} out of {} curves made into symmetric pairs.".format(
        iCt_Results, len(objrefs_In)))


if __name__ == '__main__': main()