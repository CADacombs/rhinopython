"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250312-13: Created.

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
    listValues = {}
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

    key = 'bAddArcIfInputIsNotArcCrv'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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

        if key == 'fTol_IsArc':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < cls.riOpts[key].InitialValue:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

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
        addOption('bAddArcIfInputIsNotArcCrv')
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

        if res == ri.GetResult.Number:
            key = 'fTol_IsArc'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _formatDistance(fDistance, fPrecision=None):
    if fDistance is None:
        return "(None)"
    if fDistance == Rhino.RhinoMath.UnsetValue:
        return "(Infinite)"
    if fPrecision is None:
        fPrecision = sc.doc.ModelDistanceDisplayPrecision

    if fDistance < 10.0**(-(fPrecision-1)):
        # For example, if fDistance is 1e-5 and fPrecision == 5,
        # the end of this script would display only one digit.
        # Instead, this return displays 2 digits.
        return "{:.2e}".format(fDistance)

    return "{:.{}f}".format(fDistance, fPrecision)


def _formatRadius(radius):
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
        return _formatDistance(radius)

    return "{} {} [{} millimeters]".format(
        _formatDistance(radius),
        sc.doc.GetUnitSystemName(
            modelUnits=True,
            capitalize=False,
            singular=False,
            abbreviate=False),
        _formatDistance(radius*Rhino.RhinoMath.UnitScale(
            sc.doc.ModelUnitSystem, Rhino.UnitSystem.Millimeters),
            fPrecision=4)
        )


def main():

    objrefs = getInput()
    if objrefs is None: return

    #fScaleFactor = Opts.values['fScaleFactor']
    fTol_IsArc = Opts.values['fTol_IsArc']
    bAddArcIfInputIsNotArcCrv = Opts.values['bAddArcIfInputIsNotArcCrv']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    radii = []

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
            bIsArcCrv = True
            arc = rgC_In_NotProxy.Arc
        else:
            bIsArcCrv = False

        radius = arc.Radius
        radii.append(radius)

        plane = arc.Plane

        ps = rg.PlaneSurface(
            plane,
            xExtents=rg.Interval(-radius, radius),
            yExtents=rg.Interval(-radius, radius))

        line = rg.Line(plane.Origin, span=arc.Radius*plane.Normal)

        sc.doc.Objects.AddSurface(ps)

        sc.doc.Objects.AddLine(line)

        if bAddArcIfInputIsNotArcCrv and not bIsArcCrv:
            sc.doc.Objects.AddArc(arc)

    if bEcho:
        if len(objrefs) == 1:
            if len(radii) == 0:
                print("Curve is not arc-shaped.")
            else:
                print("Curve's radius is {}.".format(_formatRadius(radius)))
        else:
            if len(radii) == 0:
                print("None of the curves are arc-shaped.")
            elif len(objrefs) > 10:
                print("Arc-shaped curves found.")
            else:
                print("Curves' radii: {}".format(
                    ", ".join([_formatRadius(radius) for radius in radii])))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()