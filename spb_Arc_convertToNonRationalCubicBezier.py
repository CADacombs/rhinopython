"""
This script converts an arc up through 180 degrees into a non-rational, cubic, Bezier NurbsCurve.

One of the two approximation options uses my solution for converting an arc into a
non-rational, cubic Bezier while maintaining the same tangent vectors and
radius at both ends:
https://discourse.mcneel.com/t/approximating-arc-with-non-rational-cubic-bezier/133784
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
220128: Created.
220314: Bug fix in curve selection.
220501: Now can input multiple curves.  Refactored.
220824: Improved printed output included reporting curve deviations.
220909: Added support for 180 degree arcs.
220917: By default, will create a cubic Bezier that is coincident with the arc
at its midpoint and matches the arc's tangent vectors at its ends.
220918: Added Opts and getInput.
221005: Bug fix.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bMatchRadiusAtEnds'; keys.append(key)
    values[key] = False
    names[key] = 'Match'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'MidPt', 'RadiusAtEnds')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fInputArcTol'; keys.append(key)
    values[key] = 1e-9
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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

        if key == 'fInputArcTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

        if key == 'fOutputDevTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue < 1e-6:
                cls.values[key] = cls.riOpts[key].CurrentValue = 1e-6
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getPreselectedCurves():
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    if gObjs_Preselected:
        gCrvs_Preselected = []
        oes = rd.ObjectEnumeratorSettings()
        oes.NormalObjects = True
        oes.LockedObjects = False
        oes.IncludeLights = False
        oes.IncludeGrips = False
        for rdRhinoObject in sc.doc.Objects.GetObjectList(oes):
            if rdRhinoObject.Id in gObjs_Preselected:
                if rdRhinoObject.ObjectType == rd.ObjectType.Curve:
                    gCrvs_Preselected.append(rdRhinoObject.Id)
        if gCrvs_Preselected:
            if Opts.values['bEcho']:
                s  = "({} curves".format(len(gCrvs_Preselected))
                s += " were preselected.)"
                print(s)
            return tuple(gCrvs_Preselected)


def getInput():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select arcs with angle <= 180 deg")

    go.GeometryFilter = rd.ObjectType.Curve


    #def customGeometryFilter(rdObj, rgObj, compIdx):
    #    if not isinstance(rgObj, rg.Curve):
    #        return False
    #    if (
    #        isinstance(rgObj, rg.NurbsCurve) and
    #        rgObj.Degree == 3 and
    #        rgObj.Points.Count == 4 and
    #        not rgObj.IsRational
    #    ):
    #        return False
    #    if rgObj.IsArc(fInputArcTol):
    #        return True
    #    return False

    #go.SetCustomGeometryFilter(customGeometryFilter)


    go.SubObjectSelect = False

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)


    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bMatchRadiusAtEnds')
        addOption('fInputArcTol')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fDevTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createPoints_matchArcRadiusAtEnds(arc):

    if (arc.Angle - Rhino.RhinoMath.ZeroTolerance) > math.pi:
        return


    if abs(arc.Angle - math.pi) <= Rhino.RhinoMath.ZeroTolerance:
        print("180 degrees")
        return Rhino.Collections.Point3dList(
            arc.StartPoint,
            arc.StartPoint + arc.Radius * 0.75**(-0.5) * arc.TangentAt(arc.AngleDomain.T0),
            arc.EndPoint   - arc.Radius * 0.75**(-0.5) * arc.TangentAt(arc.AngleDomain.T1),
            arc.EndPoint)


    cos_ = math.cos((math.pi/2.0)-arc.Angle)
    tan_ = math.tan(arc.Angle/2.0)
    sqrt_ = (1.0 + 6.0*tan_/cos_)**0.5
    t = arc.Radius * cos_ * (sqrt_ - 1.0) / 3.0

    return Rhino.Collections.Point3dList(
        arc.StartPoint,
        arc.StartPoint + t * arc.TangentAt(arc.AngleDomain.T0),
        arc.EndPoint   - t * arc.TangentAt(arc.AngleDomain.T1),
        arc.EndPoint)


def createPoints_matchArcAtMidPoint(arc):

    if (arc.Angle - Rhino.RhinoMath.ZeroTolerance) > math.pi:
        return


    # Per https://pomax.github.io/bezierinfo/#circles_cubic
    k = 4.0 * math.tan(arc.Angle/4.0) / 3.0

    return Rhino.Collections.Point3dList(
        arc.StartPoint,
        arc.StartPoint + k * arc.TangentAt(arc.AngleDomain.T0) * arc.Radius,
        arc.EndPoint   - k * arc.TangentAt(arc.AngleDomain.T1) * arc.Radius,
        arc.EndPoint)


def processCurve(crv_In, bMatchRadiusAtEnds=False, fInputArcTol=None):

    if fInputArcTol is None: fInputArcTol = 1e-9

    bSuccess, arc = crv_In.TryGetArc(fInputArcTol)
    if not bSuccess:
        #raise ValueError(
        #    "Arc should have been obtained from Curve.TryGetArc "
        #    "since Curve.IsArc returned true.")
        return


    if arc.Angle > (math.pi + Rhino.RhinoMath.ZeroTolerance):
        print("Angle of arc is {:.2f}.".format(math.degrees(arc.Angle)),
            "Script only converts up to 180 degs.",
            "Curve skipped.",
            sep="  "
            )
        return

    if bMatchRadiusAtEnds:
        point3dList = createPoints_matchArcRadiusAtEnds(arc)
    else:
        point3dList = createPoints_matchArcAtMidPoint(arc)

    return rg.Curve.CreateControlPointCurve(point3dList, degree=3)


def formatDistance(fDistance):
    if fDistance is None:
        return "(None)"
    if fDistance == Rhino.RhinoMath.UnsetValue:
        return "(Infinite)"
    if fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-1)):
        return "{:.2e}".format(fDistance)
    return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def main():

    objrefs = getInput()
    if objrefs is None: return

    bMatchRadiusAtEnds = Opts.values['bMatchRadiusAtEnds']
    fInputArcTol = Opts.values['fInputArcTol']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if bDebug:
        sEval = "bMatchRadiusAtEnds"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "fInputArcTol"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "bReplace"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "bEcho"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "bDebug"; print("{}: {}".format(sEval, eval(sEval)))
    else:
        sc.doc.Views.RedrawEnabled = False


    gCrvs_Out = []
    sLogs = []
    devs = []

    for objref in objrefs:

        crv_In = objref.Curve()

        if not crv_In: return


        nc_Res = processCurve(
            crv_In,
            bMatchRadiusAtEnds=bMatchRadiusAtEnds,
            fInputArcTol=fInputArcTol,
            )

        if nc_Res is None: continue


        rc = rg.Curve.GetDistancesBetweenCurves(nc_Res, crv_In, 0.1*sc.doc.ModelAbsoluteTolerance)
        if not rc[0]:
            sLogs.append("GetDistancesBetweenCurves failed.")
        else:
            devs.append(rc[1])

        rdObj_In = objref.Object()

        if bReplace and isinstance(rdObj_In, rd.CurveObject):
            bReplaced = sc.doc.Objects.Replace(rdObj_In.Id, nc_Res)
            if bReplaced:
                sLogs.append("Replaced curve with non-rational, cubic Bezier.")
                gCrv_Out = rdObj_In.Id
                gCrvs_Out.append(gCrv_Out)
            else:
                sLogs.append("Could not replace curve.")
                continue
        else:
            gCrv_Out = sc.doc.Objects.AddCurve(nc_Res)
            if gCrv_Out != gCrv_Out.Empty:
                sLogs.append("Added non-rational, cubic Bezier curve.")
                gCrvs_Out.append(gCrv_Out)
            else:
                sLogs.append("Could not add curve.")
                continue

    for s in set(sLogs):
        print("{} of {}".format(sLogs.count(s), s))

    if not gCrvs_Out:
        print("No curves were modified/added.")
        return

    if len(devs) == 1:
        print("Deviation: {}".format(formatDistance(devs[0])))
    elif len(devs) > 1:
        if max(devs) - min(devs) < 10**(-(sc.doc.ModelDistanceDisplayPrecision+1)):
            print("Deviations: {}".format(formatDistance(devs[0])))
        else:
            print("Deviations: [{},{}]".format(
                formatDistance(min(devs)),
                formatDistance(max(devs)),
                ))


    sc.doc.Views.Redraw()

    return gCrvs_Out


if __name__ == '__main__': main()