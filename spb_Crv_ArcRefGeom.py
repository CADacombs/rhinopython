"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250312: Created.

TODO:
    WIP: Only allow arc-shaped curves to be selected.
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
    riAddOpts = {}
    stickyKeys = {}


    #key = 'fScaleFactor'; keys.append(key)
    #values[key] = 2.0
    #names[key] = 'ScaleFactor'
    #riOpts[key] = ri.Custom.OptionDouble(
    #        initialValue=values[key],
    #        setLowerLimit=True,
    #        limit=Rhino.RhinoMath.SqrtEpsilon
    #)
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_IsArc'; keys.append(key)
    values[key] = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    names[key] = 'IsArcTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True,
        limit=Rhino.RhinoMath.ZeroTolerance)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    #key = 'bCopy'; keys.append(key)
    #values[key] = False
    #names[key] = 'Copy'
    #riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = 'Echo'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = 'Debug'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
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

        #if key == 'fSearchTol':
        #    if cls.riOpts[key].CurrentValue < 0.0:
        #        cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
        #        sc.sticky[cls.stickyKeys[key]] = cls.values[key]
        #        return

        #    sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
        #    return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def getInput(bDebug=False):
    """
    Get linear curves.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select arc-shaped curves")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    #go.GeometryAttributeFilter = (
    #    ri.Custom.GeometryAttributeFilter.

    #go.AcceptNumber(True, acceptZero=False)
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == ri.GetResult.Object 

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        #addOption('fScaleFactor')
        addOption('fTol_IsArc')
        #addOption('bCopy')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        #if res == ri.GetResult.Number:
        #    key = 'fScaleFactor'
        #    Opts.riOpts[key].CurrentValue = go.Number()
        #    Opts.setValue(key)
        #    continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def main():

    objrefs = getInput()
    if objrefs is None: return

    #fScaleFactor = Opts.values['fScaleFactor']
    fTol_IsArc = Opts.values['fTol_IsArc']
    #bCopy = Opts.values['bCopy']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    for objref in objrefs:
        rgC_In = objref.Curve()
        if rgC_In is None:
            continue

        bSuccess, arc = rgC_In.TryGetArc(tolerance=fTol_IsArc)
        if not bSuccess:
            continue

        if isinstance(rgC_In, rg.BrepEdge):
            bIsEdge = True
            rgC_In_NotProxy = rgC_In.DuplicateCurve()
        else:
            bIsEdge = False
            rgC_In_NotProxy = rgC_In

        if isinstance(rgC_In_NotProxy, rg.ArcCurve):
            arc = rgC_In_NotProxy.Arc

        radius = arc.Radius

        plane = arc.Plane

        ps = rg.PlaneSurface(
            plane,
            xExtents=rg.Interval(-radius, radius),
            yExtents=rg.Interval(-radius, radius))

        line = rg.Line(plane.Origin, span=arc.Radius*plane.Normal)

        sc.doc.Objects.AddSurface(ps)

        sc.doc.Objects.AddLine(line)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()