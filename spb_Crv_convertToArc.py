"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180905: Created.
...
190609-11: Refactored.  Modified printed output.  Removed option for conversion between knots.
210316: Modified an option default value.  Bug fix.
220831: Added an option.  Added function, replacing an import.
230818: Reenabled allowance of modifying options when curves are preselected.
241204: Modified some option default values.
250313: Removed the restriction of selecting edges in blocks.

TODO: Add bMaxRadius option.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'bTolByRatio'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTolRatio'; keys.append(key)
    values[key] = 1000.0
    names[key] = 'Ratio'
    riOpts[key] = ri.Custom.OptionDouble(
            initialValue=values[key], setLowerLimit=True, limit=1.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fDevTol'; keys.append(key)
    values[key] = max((0.1*sc.doc.ModelAbsoluteTolerance, 1e-6))
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fTanTol'; keys.append(key)
    values[key] = sc.doc.ModelAngleToleranceDegrees
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fMinNewCrvLen'; keys.append(key)
    values[key] = 100.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'MinArcLength'
    riOpts[key] = ri.Custom.OptionDouble(
            initialValue=values[key],
            setLowerLimit=True,
            limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fMaxRadius'; keys.append(key)
    values[key] = 1e5
    riOpts[key] = ri.Custom.OptionDouble(
            values[key],
            setLowerLimit=True,
            limit=sc.doc.ModelAbsoluteTolerance)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bReplace'; keys.append(key)
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
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
    
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue


def getPreselectedCurves():
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    if gObjs_Preselected:
        gCrvs_Preselected = []
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
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
    go.SetCommandPrompt("Select curves")
    
    go.GeometryFilter = rd.ObjectType.Curve


    def customGeomFilter(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if isinstance(geom, rg.BrepEdge):
            # DuplicateCurve gets the edge as a curve, which may be a subset of the EdgeCurve.
            rgC = geom.DuplicateCurve()
        elif isinstance(geom, rg.Curve):
            rgC = geom
        else:
            return False

        if isinstance(rgC, rg.ArcCurve):
            return False

        return True

    go.SetCustomGeometryFilter(customGeomFilter)


    go.AcceptNumber(True, True)
    
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    bPreselectedObjsChecked = False
    
    idxs_Opts = {}
    
    while True:
        go.ClearCommandOptions()

        go.AddOptionToggle(Opts.names['bTolByRatio'], Opts.riOpts['bTolByRatio'])
        if Opts.values['bTolByRatio']:
            print("Tolerance Ratio = Minimum (length / curve) deviation")
            go.AddOptionDouble(Opts.names['fTolRatio'], Opts.riOpts['fTolRatio'])
        else:
            go.AddOptionDouble(Opts.names['fDevTol'], Opts.riOpts['fDevTol'])
            idxs_Opts['Rh0'] = go.AddOption('Rh0')
            idxs_Opts['1eN9'] = go.AddOption('1eN9')
            idxs_Opts['1eN6'] = go.AddOption('1eN6')
            idxs_Opts['Deci'] = go.AddOption('Deci')
            idxs_Opts['ModelTol'] = go.AddOption('ModelTol')
            idxs_Opts['Deca'] = go.AddOption('Deca')
        go.AddOptionDouble(Opts.names['fTanTol'], Opts.riOpts['fTanTol'])
        go.AddOptionDouble(Opts.names['fMinNewCrvLen'], Opts.riOpts['fMinNewCrvLen'])
        go.AddOptionDouble(Opts.names['fMaxRadius'], Opts.riOpts['fMaxRadius'])
        go.AddOptionToggle(Opts.names['bReplace'], Opts.riOpts['bReplace'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
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
            gCrvs0 = []; rgEdges = []
    
            for objref in go.Objects():
                if objref.GeometryComponentIndex.Index == -1:
                    gCrvs0.append(objref.ObjectId)
                else:
                    rgE = objref.Geometry()
                    if isinstance(rgE, rg.BrepEdge):
                        rgEdges.append(objref.Geometry())
                    else:
                        rdObj = objref.Object()
                        if rdObj.ObjectType == rd.ObjectType.InstanceReference:
                            print("Objects in block instances are not supported.")
                            continue

            go.Dispose()

            return tuple([gCrvs0 + rgEdges] + [Opts.values[key] for key in Opts.keys])

        # An option was selected or a number was entered.
        key = 'fDevTol'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        elif not Opts.values['bTolByRatio']:
            if go.OptionIndex() == idxs_Opts['Rh0']:
                Opts.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            elif go.OptionIndex() == idxs_Opts['1eN9']:
                Opts.riOpts[key].CurrentValue = 1e-9
            elif go.OptionIndex() == idxs_Opts['1eN6']:
                Opts.riOpts[key].CurrentValue = 1e-6
            elif go.OptionIndex() == idxs_Opts['Deci']:
                Opts.riOpts[key].CurrentValue = 0.1*sc.doc.ModelAbsoluteTolerance
            elif go.OptionIndex() == idxs_Opts['ModelTol']:
                Opts.riOpts[key].CurrentValue = 1.0*sc.doc.ModelAbsoluteTolerance
            elif go.OptionIndex() == idxs_Opts['Deca']:
                Opts.riOpts[key].CurrentValue = 10.*sc.doc.ModelAbsoluteTolerance
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        key = 'fTanTol'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        key = 'fMinNewCrvLen'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()


def getArcCurve(rgCrv0, bTolByRatio=False, fTolRatio=1000.0, fDevTol=1e-9, fTanTol=None, fMinNewCrvLen=None, fMaxRadius=None):
    """
    Uses a different approach of obtaining the Arc or Circle than Curve.TryGetArc.
    """
    
    sType = rgCrv0.GetType().Name

    if rgCrv0 is None:
        return None, "Geometry not found!  It will be skipped."
    if sType == 'ArcCurve':
        return rgCrv0.DuplicateCurve(), 0.0
    if sType == 'LineCurve':
        return None, "Skipped LineCurve."
    if rgCrv0.IsLinear(1e-9):
        return None, "Skipped linear {}.".format(sType)

    if fMinNewCrvLen is None: fMinNewCrvLen = 10.0*sc.doc.ModelAbsoluteTolerance
    if fTanTol is None: fTanTol = sc.doc.ModelAngleToleranceRadians

    if rgCrv0.IsClosed:
        t0, t1 = rgCrv0.DivideByCount(segmentCount=3, includeEnds=False)[:2]
        rgArcCrv1 = rg.ArcCurve(
            rg.Circle(
                rgCrv0.PointAtStart,
                rgCrv0.PointAt(t0),
                rgCrv0.PointAt(t1),
            )
        )
    else:
        rgArcCrv1 = rg.ArcCurve(
            rg.Arc(
                rgCrv0.PointAtStart,
                rgCrv0.PointAt(
                        rgCrv0.DivideByCount(segmentCount=2, includeEnds=False)[0]),
                rgCrv0.PointAtEnd
            )
        )

    s = ""
    length_crv1 = rgArcCrv1.GetLength()
    if length_crv1 < fMinNewCrvLen:
        s += "Curve is too short"
    
    if fMaxRadius is None:
        modelUnitPerMm = Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters,
            sc.doc.ModelUnitSystem)
        fMaxRadius = MAX_ACCEPTABLE_SIZE_MM * modelUnitPerMm
    
    if rgArcCrv1.Radius > fMaxRadius:
        if s: s += " and "
        s += "Radius is too large."
    else:
        if s: s += "."
    if s:
        return None, s
    
    rc = rg.Curve.GetDistancesBetweenCurves(
        rgCrv0,
        rgArcCrv1,
        sc.doc.ModelAbsoluteTolerance)
    if not rc[0]:
        return None, "Deviation was not returned from GetDistancesBetweenCurves."
    
    fDev = rc[1]
    
    # Check deviation to tolerance.
    if bTolByRatio:
        fRatio = length_crv1 / fDev
        if fRatio < fTolRatio:
            return (
                None,
                "Curve was not converted because its (length / deviation) ratio is too small.")
    else:
        # Check absolute deviation.
        if fDev > fDevTol:
            return None, fDev


    angle_Start = rg.Vector3d.VectorAngle(rgArcCrv1.TangentAtStart, rgCrv0.TangentAtStart)
    if Rhino.RhinoMath.ToDegrees(angle_Start) > fTanTol:
        return None, "Tangent vector difference at start was above tolerance."

    angle_End = rg.Vector3d.VectorAngle(rgArcCrv1.TangentAtEnd, rgCrv0.TangentAtEnd)
    if Rhino.RhinoMath.ToDegrees(angle_End) > fTanTol:
        return None, "Tangent vector difference at end was above tolerance."


    return rgArcCrv1, fDev


def processObjectWithCurve(curveOrEdge0, **kwargs):
    """
    Parameters:
        curveOrEdge0 = (GUID of CurveObjects) and/or BrepEdge
        bTolByRatio,
        fTolRatio,
        fDevTol,
        fMinNewCrvLen,
        fMaxRadius,
        bReplace,
        bEcho,
        bDebug,
    Returns on success:
    Returns on fail:
    """


    def setOpt(key, value=None):
        if key in kwargs:
            return kwargs[key]
        elif key in Opts.riOpts:
            return Opts.riOpts[key].InitialValue
        else:
            return value

    bTolByRatio = setOpt('bTolByRatio')
    fTolRatio = setOpt('fTolRatio')
    fDevTol = setOpt('fDevTol')
    fTanTol = setOpt('fTanTol')
    fMinNewCrvLen = setOpt('fMinNewCrvLen')
    fMaxRadius = setOpt('fMaxRadius')
    bReplace = setOpt('bReplace')
    bEcho = setOpt('bEcho')
    bDebug = setOpt('bDebug')


    type_rdCurve_or_rgEdge = curveOrEdge0.GetType()
    
    if type_rdCurve_or_rgEdge == Guid:
        gCrv0 = curveOrEdge0
        rgCrv0 = rs.coercecurve(gCrv0)
    elif type_rdCurve_or_rgEdge == rg.BrepEdge:
        # Do not use the EdgeCurve instead since it may have a longer domain
        # than the actual BrepEdge.
        rgCrv0 = curveOrEdge0.DuplicateCurve()
        gCrv0 = None
    else:
        s = "Not CurveObject or BrepEdge."
        if bEcho: print(s)
        return None, s
    
    if isinstance(rgCrv0, rg.ArcCurve):
        return None, "Curve is already an ArcCurve."
    
    rc = getArcCurve(
        rgCrv0=rgCrv0,
        bTolByRatio=bTolByRatio,
        fTolRatio=fTolRatio,
        fDevTol=fDevTol,
        fTanTol=fTanTol,
        fMinNewCrvLen=fMinNewCrvLen,
        fMaxRadius=fMaxRadius,
        )

    if rc is None: return
    
    if not rc[0]:
        if isinstance(rc[1], float):
            if bEcho:
                s = "Required distance deviation to convert: {}".format(rc[1])
                print(s)
            return None, rc[1]
        else:
            if bEcho:
                print(rc[1])
            return None, rc[1]

    # Success.
    
    rgArcCrv1, fDev = rc
    
    if bReplace and gCrv0 is not None:
        bReplaced = sc.doc.Objects.Replace(gCrv0, rgArcCrv1)
        if bReplaced:
            if bEcho:
                s = "Curve was replaced with another with ArcCurve geometry"
                s += " with a deviation of {0:{1}f}.".format(
                                fDev, sc.doc.ModelDistanceDisplayPrecision+1)
                print(s)
            return gCrv0, fDev
        else:
            s = "Curve could not be replaced."
            if bEcho: print(s)
            return None, s
    else:
        gArcCurve1 = sc.doc.Objects.AddCurve(rgArcCrv1)
        if gArcCurve1 != Guid.Empty:
            if bEcho:
                s  = "Added an ArcCurve"
                s += " with a deviation of {0:{1}f}.".format(
                                fDev, sc.doc.ModelDistanceDisplayPrecision+1)
            return gArcCurve1, fDev
        else:
            s  = "Curve could not be added."
            if bEcho: print(s)
            return None, s


def processObjectsWithCurves(curvesAndEdges0, **kwargs):
    """
    Parameters:
        curvesAndEdges0 = (GUIDs of CurveObjects) and/or BrepEdges
        bTolByRatio,
        fTolRatio,
        fDevTol,
        fTanTol,
        fMinNewCrvLen,
        fMaxRadius,
        bReplace,
        bEcho,
        bDebug,
    Returns on success:
    Returns on fail:
    """


    def setOpt(key, value=None):
        if key in kwargs:
            return kwargs[key]
        elif key in Opts.riOpts:
            return Opts.riOpts[key].InitialValue
        else:
            return value

    bTolByRatio = setOpt('bTolByRatio')
    fTolRatio = setOpt('fTolRatio')
    fDevTol = setOpt('fDevTol')
    fTanTol = setOpt('fTanTol')
    fMinNewCrvLen = setOpt('fMinNewCrvLen')
    fMaxRadius = setOpt('fMaxRadius')
    bReplace = setOpt('bReplace')
    bEcho = setOpt('bEcho')
    bDebug = setOpt('bDebug')


    gCrvs0_Replaced = []
    gArcCrvs1 = []
    devs_all = []
    
    fTols_needed = []
    sFails = []
    
    gCrvs0 = []
    for curveOrEdge0 in curvesAndEdges0:
        gCrv0 = rs.coerceguid(curveOrEdge0)
        if gCrv0:
            gCrvs0.append(gCrv0)
    
    idxs_AtTenths = [int(round(0.1*i*len(curvesAndEdges0),0)) for i in range(10)]
    
    for iC, curveOrEdge0 in enumerate(curvesAndEdges0):
        if iC in idxs_AtTenths:
            Rhino.RhinoApp.SetCommandPrompt(
                    "Processing curve {} ...".format(
                    "" if len(curvesAndEdges0) == 1 else "{} of {} ".format(
                        iC+1, len(curvesAndEdges0))))
        
        rc = processObjectWithCurve(
                curveOrEdge0=curveOrEdge0,
                bTolByRatio=bTolByRatio,
                fTolRatio=fTolRatio,
                fDevTol=fDevTol,
                fTanTol=fTanTol,
                fMinNewCrvLen=fMinNewCrvLen,
                fMaxRadius=fMaxRadius,
                bReplace=bReplace,
                bEcho=False if len(curvesAndEdges0) > 1 else bEcho,
                bDebug=bDebug,
        )
        if rc is None: continue
        if rc[0] is None:
            if isinstance(rc[1], float):
                fTols_needed.append(rc[1])
            else:
                sFails.append(rc[1])
            continue
        
        gArcCrvX, dev = rc
        if gArcCrvX in gCrvs0:
            gCrvs0_Replaced.append(gArcCrvX)
        else:
            gArcCrvs1.append(gArcCrvX)
        devs_all.append(dev)
    
    if not bEcho:
        return gCrvs0_Replaced, gArcCrvs1


    if len(curvesAndEdges0) > 1:
        s = "Out of {} total curves:\n".format(len(curvesAndEdges0))
    else:
        s = ""

    for sFail in set(sFails):
        s += "\n[{}] {}".format(sFails.count(sFail), sFail)

    if fTols_needed:
        s += "\n[{}] Arc not found within tolerance.".format(len(fTols_needed))
        s += "  Tolerances needed: [{0:.{2}e},{1:.{2}e}]".format(
                min(fTols_needed), max(fTols_needed), 3)

    if gCrvs0_Replaced:
        if len(gCrvs0_Replaced) == len(curvesAndEdges0):
            s += "\n[{}] Curve was replaced.".format(len(curvesAndEdges0))
        else:
            s += "\n[{}] Replaced curve.".format(len(gCrvs0_Replaced))
        s += "  Deviations: [{0:.{2}e},{1:.{2}e}]".format(
                min(devs_all), max(devs_all), 3)

    if gArcCrvs1:
        if len(gArcCrvs1) == len(curvesAndEdges0):
            s += "\n[{}] Curve was added.".format(len(curvesAndEdges0))
        else:
            s += "\n[{}] Curve was added.".format(len(gArcCrvs1))
        s += "  Deviations: [{0:.{2}e},{1:.{2}e}]".format(
                min(devs_all), max(devs_all), 3)

    if not (gCrvs0_Replaced or gArcCrvs1):
        s += "\nNo curves were added or replaced."

    print(s)


    return gCrvs0_Replaced, gArcCrvs1


def main():
    
    gCrvs0_Preselected = getPreselectedCurves()
    
    rc = getInput()
    if rc is None: return
    (
        curvesAndEdges0,
        bTolByRatio,
        fTolRatio,
        fDevTol,
        fTanTol,
        fMinNewCrvLen,
        fMaxRadius,
        bReplace,
        bEcho,
        bDebug,
        ) = rc
    
    if bEcho and not bTolByRatio:
        print("Converting curves to arcs within a tolerance of {} ...".format(fDevTol))
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll()
    
    gCrvs0_Replaced, gArcCrvs1 = None, None
    
    rc = processObjectsWithCurves(
        curvesAndEdges0=curvesAndEdges0,
        bTolByRatio=bTolByRatio,
        fTolRatio=fTolRatio,
        fDevTol=fDevTol,
        fTanTol=fTanTol,
        fMinNewCrvLen=fMinNewCrvLen,
        fMaxRadius=fMaxRadius,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if rc is not None:
        gCrvs0_Replaced, gArcCrvs1 = rc
    
    if gCrvs0_Preselected:
        [sc.doc.Objects.Select(objectId=_) for _ in gCrvs0_Preselected]
    else:
        if gCrvs0_Replaced:
            [sc.doc.Objects.Select(objectId=_) for _ in gCrvs0_Replaced]
        if gArcCrvs1:
            [sc.doc.Objects.Select(objectId=_) for _ in gArcCrvs1]
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()