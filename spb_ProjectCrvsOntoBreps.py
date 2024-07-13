"""
This script is an alternative to _Project.
Segments of polycurves can be projected individually, avoiding merges through joints.
Projections of curves onto planar surfaces of polyface breps are more similar to loose projections.
For increased accuracy, all curves projected to planar faces are projected loose
to the TryGetPlane plane, then split at the intersections with the monoface brep's edges.
With AttemptRebuild enabled, curve will be projected at half the tolerance,
and a rebuild (to uniform) of the projected curve will be attempted at half tolerance.
If not successful, the projection at full tolerance will be used.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
191124-25: Created.
...
220821: Modified an option default value.
230720: Now passes fTol to a method instead of ModelAbsoluteTolerance.
231106: Improved selection routine for faces.
240712: Refactored and modified behavior of simplification routines.

TODO:
    Determine solution for when a loose projected curve doesn't lie within
    the projection boundary of the brep.
    Possible solutions:
        1. Keep brep-missing Greville points at their non-projected locations.
        2. Trim the curve before projection.
    
    Add support for projections of curves in overlapping directions, e.g., projecting a circle parallel to its plane.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscript.utility
import rhinoscript.userinterface
import scriptcontext as sc

from System import Guid


sOpts_OutputLayer = [
    'Input',
    'Current',
    'TargetObject',
    ]

sOpts_Direction = [
    'CPlaneX',
    'CPlaneY',
    'CPlaneZ',
    'WorldX',
    'WorldY',
    'WorldZ',
    'View',
    'Custom',
    ]


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'bProjectCrvSegs'; keys.append(key)
    values[key] = True
    names[key] = "Project"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'WholeCrv', 'CrvSegs')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLoose'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPostProcess'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryOnlyUniformCubicForNonlinear'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryGetArcs'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAcceptRational'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryRebuildOthersUniform'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteInput'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iOutputLayer'; keys.append(key)
    values[key] = 1
    names[key] = key[1:]
    riAddOpts[key] = addOptionList(key, names, sOpts_OutputLayer, values)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDirection'; keys.append(key)
    values[key] = 2
    names[key] = key[1:]
    riAddOpts[key] = addOptionList(key, names, sOpts_Direction, values)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'vectCustom'; keys.append(key)
    values[key] = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(rdObjs_toHighlight, sPrompt, rdGeomFilter):
    """
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt(prompt=sPrompt)

    go.GeometryFilter = rdGeomFilter

    #go.AlreadySelectedObjectSelect = True # Default is False
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    idxs_Opts = {}

    bPreselectedObjsChecked = False
    go.EnablePreSelect(True, ignoreUnacceptablePreselectedObjects=True)

    while True:
        Opts.riAddOpts['bProjectCrvSegs'](go)
        Opts.riAddOpts['bLoose'](go)
        Opts.riAddOpts['fTol'](go)
        if not Opts.values['bLoose']: Opts.riAddOpts['bPostProcess'](go)
        if Opts.values['bPostProcess']:
            Opts.riAddOpts['bTryOnlyUniformCubicForNonlinear'](go)
            if not Opts.values['bTryOnlyUniformCubicForNonlinear']:
                Opts.riAddOpts['bTryGetArcs'](go)
                Opts.riAddOpts['bAcceptRational'](go)
                Opts.riAddOpts['bTryRebuildOthersUniform'](go)
        Opts.riAddOpts['bDeleteInput'](go)
        idxs_Opts['iOutputLayer'] = Opts.riAddOpts['iOutputLayer'](go)
        idxs_Opts['iDirection'] = Opts.riAddOpts['iDirection'](go)
        Opts.riAddOpts['bEcho'](go)
        Opts.riAddOpts['bDebug'](go)

        if Opts.values['bDebug']:
            print("Before GetMultiple")
            sEval = "  go.ObjectCount"; print("{}: {}".format(sEval, eval(sEval)))

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if Opts.values['bDebug']:
            print("After GetMultiple")
            sEval = "  go.ObjectCount"; print("{}: {}".format(sEval, eval(sEval)))

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if Opts.values['bDebug']:
            sEval = "  bPreselectedObjsChecked"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "  go.ObjectsWerePreselected"; print("{}: {}".format(sEval, eval(sEval)))

        if not bPreselectedObjsChecked:
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            bPreselectedObjsChecked = True
            if go.ObjectsWerePreselected:
                continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return [objrefs] + [Opts.values[key] for key in Opts.keys]

        # An option was selected or a number was entered.
        if res == ri.GetResult.Number:
            Opts.riOpts['fTol'].CurrentValue = go.Number()
        else:
            if go.Option().Index == idxs_Opts['iOutputLayer']:
                Opts.values['iOutputLayer'] = (
                        go.Option().CurrentListOptionIndex)
            elif go.Option().Index == idxs_Opts['iDirection']:
                Opts.values['iDirection'] = (
                        go.Option().CurrentListOptionIndex)

                if sOpts_Direction[go.Option().CurrentListOptionIndex] == 'Custom':
                    rc = rhinoscript.userinterface.GetLine(
                        mode=1, point=None,
                        message1='Projection direction',
                        message3='Second direction point',
                        )
                    if not rc:
                        Opts.values['iDirection'] = 0
                    else:
                        Opts.values['vectCustom'] = rg.Vector3d(rc[1] - rc[0])
                        Opts.values['vectCustom'].Unitize()

        key = 'fTol'
        if Opts.riOpts[key].CurrentValue <= (1.0/2**32):
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.Brep):
        return rhObj
    elif isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        geom = rhObj.Geometry()
        if rhObj.GeometryComponentIndex.ComponentIndexType == rg.ComponentIndexType.BrepFace:
            geom = geom.DuplicateFace(duplicateMeshes=False)
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    if isinstance(geom, rg.Brep):
        return geom


def formatDistance(fDistance):
    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def projectCurve_Loose(rgC0, rgB, vectDir, bProjectCrvSegs=True, fTol=None, bDebug=False):
    """
    Projecting to breps of individual faces results in simpler curves for any
    planar (not only PlaneSurface) faces.
    For more accurate results, the curves's control points should be projected to
    the Plane themselves.
    """

    if bProjectCrvSegs:
        cs_toProj = rgC0.DuplicateSegments()
        if not cs_toProj:
            cs_toProj = [rgC0.DuplicateCurve()]
    else:
        cs_toProj = [rgC0.DuplicateCurve()]

    if fTol is None: fTol = sc.doc.ModelAbsoluteTolerance


    rgCs_Proj_ThisC0 = []

    for c_toProj in cs_toProj:
        # Duplicate curve and translate Greville points to projected locations.
        nc_toProj = c_toProj.ToNurbsCurve()

        grPts1 = []
        for gr in nc_toProj.GrevillePoints(all=False):
            pts_Proj = rg.Intersect.Intersection.ProjectPointsToBreps(
                breps=[rgB],
                points=[gr],
                direction=vectDir,
                tolerance=0.1*fTol)
            if pts_Proj:
                grPts1.append(pts_Proj[0])
            else:
                grPts1.append(gr)
            
            #sc.doc.Objects.AddPoint(pts_Proj[0])
        #sc.doc.Views.Redraw(); 1/0

        nc_toProj.SetGrevillePoints(grPts1)

        rgCs_Proj_ThisC0.append(nc_toProj)

    rgCs_Proj_Joined_ThisC0 = rg.Curve.JoinCurves(rgCs_Proj_ThisC0)

    # This will clean up any short segments that had previous passed
    # the previous checks, but have been deformed to short segments
    # after JoinCurves.
    # _Project doesn't handle this problem.
    for rgC in rgCs_Proj_Joined_ThisC0:
        rgC.RemoveShortSegments(fTol)

    return rgCs_Proj_Joined_ThisC0


def projectCurve_Tight(rgC0, rgBs_1Face, rgPlanes, vectDir, **kwargs):
    """
    Parameters:
        bProjectCrvSegs
        fTol
        bPostProcess
        bTryOnlyUniformCubicForNonlinear
        bTryGetArcs
        bAcceptRational
        bTryRebuildOthersUniform


    Projecting to breps of individual faces results in simpler curves for any
    planar (not only PlaneSurface) faces.
    For more accurate results, the curves's control points should be projected to
    the Plane themselves.
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bProjectCrvSegs = getOpt('bProjectCrvSegs')
    fTol = getOpt('fTol')
    bPostProcess = getOpt('bPostProcess')
    bTryOnlyUniformCubicForNonlinear = getOpt('bTryOnlyUniformCubicForNonlinear')
    bTryGetArcs = getOpt('bTryGetArcs')
    bAcceptRational = getOpt('bAcceptRational')
    bTryRebuildOthersUniform = getOpt('bTryRebuildOthersUniform')
    bDebug = getOpt('bDebug')


    sc.doc.Objects.UnselectAll() # For debug.


    def projectCurveLooseToPlane(rgC_In, rgPlane):
        rgC_Working = rgC_In.ToNurbsCurve()

        grPts1 = []

        for gr0 in rgC_Working.GrevillePoints(all=False):

            rgLine = rg.Line(gr0, vectDir)
            bSuccess, tLine = Rhino.Geometry.Intersect.Intersection.LinePlane(
                rgLine, rgPlane);
            if not bSuccess:
                s  = "Point missed the plane when projecting loose to a plane."
                s += "  Check results for accuracy."
                print(s)
                rgC_Working.Dispose()
                return

            gr1 = rgLine.PointAt(tLine);

            grPts1.append(gr1)

        rgC_Working.SetGrevillePoints(grPts1)
        
        return rgC_Working


    def projectCrvToGeom(crv_In, tolerance):
        rgCs_Proj_ThisSeg = []
        
        rgPlanes_ = rgPlanes # For debugging.

        for rgB, rgPlane in zip(rgBs_1Face, rgPlanes_):
            rgCs_Proj_ThisSeg_1F = rg.Curve.ProjectToBrep(
                curve=crv_In,
                brep=rgB,
                direction=vectDir,
                tolerance=tolerance)

            if not rgCs_Proj_ThisSeg_1F: continue

            #for c in rgCs_Proj_ThisSeg_1F: sc.doc.Objects.AddCurve(c)
            #sc.doc.Views.Redraw()
            pass

            WipList = []

            for rgC in rgCs_Proj_ThisSeg_1F:
                if rgC.GetLength() > tolerance:
                    WipList.append(rgC)
                else:
                    rgC.Dispose()
                    print("Short curve created at {} tolerance ignored.".format(
                        formatDistance(tolerance)))

            rgCs_Proj_ThisSeg_1F = WipList

            if not rgCs_Proj_ThisSeg_1F: continue

            #            # Temp for debug. ############################
            #            rgCs_Proj_ThisSeg.extend(rgCs_Proj_ThisSeg_1F)
            #            continue

            if not rgPlane:
                rgCs_Proj_ThisSeg.extend(rgCs_Proj_ThisSeg_1F)
                continue


            # Curves are projected to a planar surface.


            # Results of ProjectToBrep:
            #   ArcCurves and NurbCurves with more than 2 points
            # are not projected accurately enough (sometimes only 3 decimal places).
            #   LineCurves, PolylineCurves, and 2-point NurbsCurves
            # are projected accurately.
            for c in rgCs_Proj_ThisSeg_1F:
                if isinstance(c, rg.NurbsCurve) and c.Points.Count > 2:
                    break
                if isinstance(c, rg.ArcCurve):
                    break
                elif isinstance(c, rg.NurbsCurve) and c.Points.Count == 2:
                    continue
                if isinstance(c, (rg.LineCurve, rg.PolylineCurve)):
                    continue
            else:
                # Yes, all of rgCs_Proj_ThisSeg_1F are goo.
                rgCs_Proj_ThisSeg.extend(rgCs_Proj_ThisSeg_1F)
                continue


            # Project loose to the Plane for higher accuracy than Curve.ProjectToBrep.

            seg1_PostProject_PreSplit = projectCurveLooseToPlane(crv_In, rgPlane)
            if seg1_PostProject_PreSplit is None:
                # Use tight projection instead.
                rgCs_Proj_ThisSeg.extend(rgCs_Proj_ThisSeg_1F)
                continue


            # Split.

            ts_atSplits = []

            for rgEdge in rgB.Edges:
                crvinters = rg.Intersect.Intersection.CurveCurve(
                    seg1_PostProject_PreSplit,
                    rgEdge,
                    tolerance=tolerance,
                    overlapTolerance=0.0)
                if crvinters.Count == 0: continue # to next curve.
                for crvinter in crvinters:
                    ts_atSplits.append(crvinter.ParameterA)

            if not ts_atSplits: 
                # Projected curve may lie completely within the face.
                rgCs_Proj_ThisSeg.append(seg1_PostProject_PreSplit)
                continue

            segs2_PostSplit = seg1_PostProject_PreSplit.Split(ts_atSplits)

            # Determine which of the split project-to-plane curves should be kept.
            pts_toDetermineTrims = []
            for rgC in rgCs_Proj_ThisSeg_1F:
                pts_toDetermineTrims.append(rgC.PointAt(rgC.Domain.Mid))
                rgC.Dispose()

            for rgC in segs2_PostSplit:
                for pt in pts_toDetermineTrims:
                    bSuccess, t = rgC.ClosestPoint(pt)
                    if not bSuccess: continue # To next point.
                    if rgC.PointAt(t).DistanceTo(pt) <= 0.1 * tolerance:
                        # Tight tolerance above is to avoid grabbing segments
                        # adjacent to correct results that also happen to be short.
                        rgCs_Proj_ThisSeg.append(rgC)
                        break # To next curve from split.


        #for rgC in rgCs_Proj_ThisSeg:
        #    sc.doc.Objects.AddCurve(rgC)
        #sc.doc.Views.Redraw(); 1/0

        for rgC in rgCs_Proj_ThisSeg:
            rgC.RemoveShortSegments(tolerance)

        return rgCs_Proj_ThisSeg


    def removeShortestCountDiffSegs():
        print("Counts of HalfTol vs. FullTol projected curves are different.")
        sEval='len(rgCs_Proj_ThisSeg_HalfTol)'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='len(rgCs_Proj_ThisSeg_FullTol)'; print("{}: {}".format(sEval, eval(sEval)))
        if len(rgCs_Proj_1Seg_HalfTol) > len(rgCs_Proj_1Seg_FullTol):
            # Remove n shortest segments of difference in counts.
            iCt_RemoveShortest = len(rgCs_Proj_1Seg_HalfTol) - len(rgCs_Proj_1Seg_FullTol)

            fLengths = []
            for rgC in rgCs_Proj_1Seg_HalfTol:
                fLength = rgC.GetLength()
                fLengths.append(fLength)

            for i in range(iCt_RemoveShortest):
                idx_ToRemove = fLengths.index(min(fLengths))
                del rgCs_Proj_1Seg_HalfTol[idx_ToRemove]
                del fLengths[idx_ToRemove]

            print("Removed {} segments from HalfTol result.".format(iCt_RemoveShortest))
        else:
            print("Check results!")


    def hasInternalPolyknots(nc):
        if not isinstance(nc, rg.NurbsCurve): return

        # Bezier?
        if nc.Points.Count == nc.Degree + 1:
            return False

        start = 0 if nc.IsPeriodic else nc.Degree
        end = nc.Knots.Count - (0 if nc.IsPeriodic else nc.Degree) - 1

        for i in range(start, end):
            if nc.Knots.KnotMultiplicity(i) > 1:
                return True

        return False


    def areEpsilonEqual(a, b, epsilon):
        # This is a relative comparison.
        delta = abs(a - b)
        fRelComp = delta / max(abs(a), abs(b))
        return fRelComp < epsilon


    def isUniformNurbsCurve(nc):
        if not isinstance(nc, rg.NurbsCurve): return


        # Bezier?
        if nc.Points.Count == nc.Degree + 1:
            return True


        # Any internal polyknots?
        start = 0 if nc.IsPeriodic else nc.Degree
        end = nc.Knots.Count - (0 if nc.IsPeriodic else nc.Degree) - 1

        for i in range(start, end):
            if nc.Knots.KnotMultiplicity(i) > 1:
                return False


        # Any non-uniform knot spans?
        start = 0 if nc.IsPeriodic else nc.Degree - 1
        end = nc.Knots.Count - (0 if nc.IsPeriodic else nc.Degree - 1) - 1

        span0 = nc.Knots[start+1] - nc.Knots[start]

        for i in range(start+1, end):
            if nc.Knots.KnotMultiplicity(i) > 1:
                return False
            if not areEpsilonEqual(
                span0, nc.Knots[i+1] - nc.Knots[i],
                epsilon=(1.0/2**32)):
                    return False

        return True


    def rebuild(nc, tolerance):

        for iDeg in 2, 3, 5:
            # Try rebuilding as a Bezier.
            rebuilt = nc.Rebuild(
                pointCount=iDeg + 1,
                degree=iDeg,
                preserveTangents=True)

            bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
                nc, rebuilt, tolerance=0.1*tolerance)[:2]

            if bSuccess and fDistMax <= tolerance:
                return rebuilt

            rebuilt.Dispose()


        iCt_MaxCp = round(nc.GetLength() / (100.0 * sc.doc.ModelAbsoluteTolerance))


        # Try rebuilding at Degree 5 and maximum allowed control point count.

        rebuilt = nc.Rebuild(
            pointCount=iCt_MaxCp,
            degree=5,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            nc, rebuilt, tolerance=0.1*tolerance)[:2]

        if not bSuccess or fDistMax > tolerance:
            # Fail at maximum control point count, so quit searching.
            rebuilt.Dispose()
            return


        # Binary search.
        rebuilt_LastSuccess = rebuilt
        iCts_Cps_Tried = [4, iCt_MaxCp]

        iCt_Cp_Try = (iCt_MaxCp + 4) // 2

        if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

        while iCt_Cp_Try not in iCts_Cps_Tried:
            sc.escape_test()


            for iDeg in 3, 5:

                if iCt_Cp_Try < iDeg+1:
                    pass

                rebuilt = nc.Rebuild(
                    pointCount=iCt_Cp_Try,
                    degree=iDeg,
                    preserveTangents=True)

                bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
                    nc, rebuilt, tolerance=0.1*tolerance)[:2]

                if bSuccess and fDistMax <= tolerance:
                    break


            iCts_Cps_Tried.append(iCt_Cp_Try)
            iCts_Cps_Tried.sort()

            if bSuccess and fDistMax <= tolerance:
                rebuilt_LastSuccess.Dispose()
                rebuilt_LastSuccess = rebuilt
                # Bisect left.
                iCt_Cp_Try = (
                    (iCt_Cp_Try +
                        iCts_Cps_Tried[iCts_Cps_Tried.index(iCt_Cp_Try)-1]) // 2)
            else:
                rebuilt.Dispose()
                # Bisect right.
                iCt_Cp_Try = (
                    (iCt_Cp_Try +
                        iCts_Cps_Tried[iCts_Cps_Tried.index(iCt_Cp_Try)+1]) // 2)

            if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

        return rebuilt_LastSuccess


    def rebuildUniformNonrationalCubic(rgCrv, tolerance):

        # Try rebuilding as a Bezier.
        rebuilt = rgCrv.Rebuild(
            pointCount=4,
            degree=3,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            rgCrv, rebuilt, tolerance=0.1*tolerance)[:2]

        if bSuccess and fDistMax <= tolerance:
            return rebuilt

        rebuilt.Dispose()


        iCt_MaxCp = round(rgCrv.GetLength() / (100.0 * sc.doc.ModelAbsoluteTolerance))


        # Try rebuilding at maximum allowed control point count.

        rebuilt = rgCrv.Rebuild(
            pointCount=iCt_MaxCp,
            degree=3,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            rgCrv, rebuilt, tolerance=0.1*tolerance)[:2]

        if not bSuccess or fDistMax > tolerance:
            # Fail at maximum control point count, so quit searching.
            rebuilt.Dispose()
            return


        # Binary search to find least control point count.
        rebuilt_LastSuccess = rebuilt
        iCts_Cps_Tried = [4, iCt_MaxCp]

        iCt_Cp_Try = (iCt_MaxCp + 4) // 2

        if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

        while iCt_Cp_Try not in iCts_Cps_Tried:
            sc.escape_test()


            # Degree 3.

            rebuilt = rgCrv.Rebuild(
                pointCount=iCt_Cp_Try,
                degree=3,
                preserveTangents=True)

            bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
                rgCrv, rebuilt, tolerance=0.1*tolerance)[:2]


            iCts_Cps_Tried.append(iCt_Cp_Try)
            iCts_Cps_Tried.sort()

            if bSuccess and fDistMax <= tolerance:
                rebuilt_LastSuccess.Dispose()
                rebuilt_LastSuccess = rebuilt
                # Bisect left.
                iCt_Cp_Try = (
                    (iCt_Cp_Try +
                        iCts_Cps_Tried[iCts_Cps_Tried.index(iCt_Cp_Try)-1]) // 2)
            else:
                rebuilt.Dispose()
                # Bisect right.
                iCt_Cp_Try = (
                    (iCt_Cp_Try +
                        iCts_Cps_Tried[iCts_Cps_Tried.index(iCt_Cp_Try)+1]) // 2)

            if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

        return rebuilt_LastSuccess


    def explodePolyCrvs(rgCrvs_In):
        # Explode PolyCurves.
        # A PolyCurve from a projection may result when a start/end
        # of a closed curve crosses over a face,
        # but the entire curve wasn't projected to the face.
        proj_NoPCs = []
        for proj in rgCrvs_In:
            if isinstance(proj, rg.PolyCurve):
                proj.RemoveNesting()
                fromExplode = proj.Explode()
                proj_NoPCs.extend(fromExplode)
                proj.Dispose()
            else:
                proj_NoPCs.append(proj)

        return proj_NoPCs[:]


    def convertPolyCrvsToNurbs(rgCrvs_toMod):
        for i, c in enumerate(rgCrvs_toMod):
            if isinstance(c, rg.PolyCurve):
                rgCrvs_toMod[i] = c.ToNurbsCurve()
                c.Dispose()


    def areAllNurbsCrvsUniform(rgCrvs):
        for c in rgCrvs:
            if not isinstance(c, (rg.NurbsCurve)):
                continue
            if not isUniformNurbsCurve(c):
                return False
        return True


    def areAllCrvsUniformNonrationalCubic_or_Lines(rgCrvs):
        for c in rgCrvs:
            if isinstance(c, rg.LineCurve):
                continue
            if not isinstance(c, (rg.NurbsCurve)):
                # Is Arc, Polyline?, PolyCurve?
                return False
            if c.IsRational:
                return False
            if not c.Degree == 3:
                return False
            if not isUniformNurbsCurve(c):
                return False
        return True


    def tryConvertNurbsCrvsToLinearCubicBeziers(rgCrvs_toMod, tol):
        """
        rgCrvs_toMod is modified.
        Returns: None
        """
        for i, c in enumerate(rgCrvs_toMod):
            if isinstance(c, rg.NurbsCurve):
                if not c.IsClosed:
                    if c.IsLinear(tol):
                        Lc = rg.LineCurve(
                            c.PointAtStart,
                            c.PointAtEnd)
                        rgCrvs_toMod[i] = Lc.ToNurbsCurve()
                        rgCrvs_toMod[i].IncreaseDegree(3)
                        c.Dispose()


    def tryConvertNurbsToLines(rgCrvs_toMod, tol):
        """
        rgCrvs_toMod is modified.
        Returns: None
        """
        for i, c in enumerate(rgCrvs_toMod):
            if isinstance(c, rg.NurbsCurve):
                if not c.IsClosed:
                    if c.IsLinear(tol):
                        rgCrvs_toMod[i] = rg.LineCurve(
                            c.PointAtStart,
                            c.PointAtEnd)
                        c.Dispose()


    def tryConvertNurbsToArcs(rgCrvs_toMod, tol):
        """
        rgCrvs_toMod is modified.
        Returns: None
        """
        fTol_MaxRad = 1e3
        
        for i, c in enumerate(rgCrvs_toMod):
            if isinstance(c, rg.NurbsCurve):
                if c.IsClosed:
                    b, circle = c.TryGetCircle(tol)
                    if b:
                        if circle.Radius <= fTol_MaxRad:
                            rgCrvs_toMod[i] = rg.ArcCurve(circle)
                            c.Dispose()
                else:
                    b, arc = c.TryGetArc(tol)
                    if b:
                        if arc.Radius <= fTol_MaxRad:
                            rgCrvs_toMod[i] = rg.ArcCurve(arc)
                            c.Dispose()


    def areAllCurvesSimplified(rgCrvs):
        for c in rgCrvs:
            if isinstance(c, rg.LineCurve):
                pass
            elif isinstance(c, rg.ArcCurve):
                if not bTryGetArcs:
                    return False
            elif isinstance(c, rg.NurbsCurve):
                if c.IsRational:
                    if not bAcceptRational:
                        return False
                if bTryRebuildOthersUniform and not isUniformNurbsCurve(c):
                    return False
            else:
                raise ValueError("{} in areAllCurvesSimplified.".format(c.GetType().Name))

        return True


    def reportCrvTypes(rgCrvs):
        sTypes = [c.GetType().Name for c in rgCrvs]
        for sType in set(sTypes):
            print("[{}] {}".format(sTypes.count(sType), sType))


    def projectThenProcessCrv_UniformCubicOnly(c_toProj, fTol_Proj_Total):
        projs_FullTol = projectCrvToGeom(c_toProj, tolerance=fTol_Proj_Total)
        if not projs_FullTol:
            return

        if bDebug:
            print('-'*20)
            reportCrvTypes(projs_FullTol)

        if rgC0.IsPeriodic:
            convertPolyCrvsToNurbs(projs_FullTol)
        else:
            projs_FullTol = explodePolyCrvs(projs_FullTol)

        tryConvertNurbsToLines(projs_FullTol, 1e-6)

        if areAllCrvsUniformNonrationalCubic_or_Lines(projs_FullTol):
            return projs_FullTol

        # Project to a different tolerance and rebuild.

        for fTol_Proj_WIP in (0.5*fTol_Proj_Total, 0.1*fTol_Proj_Total):

            projs_atTrialTol = projectCrvToGeom(c_toProj, tolerance=fTol_Proj_WIP)

            if rgC0.IsPeriodic:
                convertPolyCrvsToNurbs(projs_atTrialTol)
            else:
                projs_atTrialTol = explodePolyCrvs(projs_atTrialTol)

            if areAllCrvsUniformNonrationalCubic_or_Lines(projs_atTrialTol):
                return projs_atTrialTol

            for i, proj in enumerate(projs_atTrialTol):
                if areAllCrvsUniformNonrationalCubic_or_Lines([proj]):
                    continue
                rebuilt = rebuildUniformNonrationalCubic(
                    proj,
                    tolerance=fTol_Proj_Total-fTol_Proj_WIP)

                if rebuilt:
                    projs_atTrialTol[i] = rebuilt
                    proj.Dispose()
                    continue
                else:
                    for c in projs_atTrialTol: c.Dispose()
                    break # out of loop of projs_atTrialTol to next tolerance for rebuild.

            else:
                for proj in projs_FullTol: proj.Dispose()
                return projs_atTrialTol
        else:
            if bDebug:
                print("All rebuild tolerances failed for segment.")
            return projs_FullTol


    def projectThenProcessCrv_NotOnly_UniformCubic(c_toProj, fTol_Proj_Total):
        projs_FullTol = projectCrvToGeom(c_toProj, tolerance=fTol_Proj_Total)
        if not projs_FullTol:
            return

        if bDebug:
            print('-'*20)
            reportCrvTypes(projs_FullTol)

        tryConvertNurbsToLines(projs_FullTol, fTol_Proj_Total)

        if bTryGetArcs:
            tryConvertNurbsToArcs(projs_FullTol, fTol_Proj_Total)

        if rgC0.IsPeriodic:
            convertPolyCrvsToNurbs(projs_FullTol)
        else:
            projs_FullTol = explodePolyCrvs(projs_FullTol)

        if areAllCurvesSimplified(projs_FullTol):
            return projs_FullTol

        # Project to a different tolerance and rebuild.

        for fTol_Proj_WIP in (0.5*fTol_Proj_Total, 0.1*fTol_Proj_Total):

            projs_atTrialTol = projectCrvToGeom(c_toProj, tolerance=fTol_Proj_WIP)

            tryConvertNurbsToLines(projs_atTrialTol, fTol_Proj_Total-fTol_Proj_WIP)

            if bTryGetArcs:
                tryConvertNurbsToArcs(projs_atTrialTol, fTol_Proj_Total-fTol_Proj_WIP)

            if rgC0.IsPeriodic:
                convertPolyCrvsToNurbs(projs_atTrialTol)
            else:
                projs_atTrialTol = explodePolyCrvs(projs_atTrialTol)

            if areAllCurvesSimplified(projs_atTrialTol):
                for proj in projs_FullTol: proj.Dispose()
                return projs_atTrialTol

            potentials = []
            for proj in projs_atTrialTol:
                if isinstance(proj, (rg.LineCurve, rg.ArcCurve)):
                    potentials.append(proj)
                    continue
            
                if not isinstance(proj, rg.NurbsCurve):
                    raise ValueError("Not a NurbsCurve!")

                if proj.IsRational and rg.Curve.IsEllipse(proj, tolerance=fTol_Proj_Total-fTol_Proj_WIP):
                    potentials.append(proj)
                elif not bTryRebuildOthersUniform and proj.IsRational and isUniformNurbsCurve(proj):
                    potentials.append(proj)
                elif isUniformNurbsCurve(proj):
                    potentials.append(proj)
                else:
                    # Rebuild to complement of projection tolerance.
                    rebuilt = rebuild(
                        proj,
                        tolerance=fTol_Proj_Total-fTol_Proj_WIP)
                    if rebuilt:
                        potentials.append(rebuilt)
                        continue
                    else:
                        for rgC in potentials: rgC.Dispose()
                        break # out of for loop to next tolerance for rebuild.

            else:
                # Successful simplification.
                for proj in projs_FullTol: proj.Dispose()
                return potentials
        else:
            if bDebug:
                print("All rebuild tolerances failed for segment.")
            return projs_FullTol


    def getJoints(rgCrvs_In):
        # Use JoinCurves to quickly obtain joint locations.

        joined = rg.Curve.JoinCurves(rgCrvs_In)
        joints = []
        for rgC in joined:
            if not isinstance(rgC, rg.PolyCurve): continue

            if rgC.IsClosed:
                joints.append(rgC.PointAtStart)

            for iSeg in range(1, rgC.SegmentCount):
                seg = rgC.SegmentCurve(iSeg)
                joints.append(seg.PointAtStart)
                seg.Dispose()

            rgC.Dispose()

        return joints


    def joinCurves(rgCs_In):

        joints = getJoints(rgCs_In)

        rgCrvs_toJoin = [rgC.DuplicateCurve() for rgC in rgCs_In]

        iJ = 0

        # Before JoinCurves, adjust ends of NurbsCurves.

        for joint in joints:
            sc.escape_test()

            iCt_EndsProcessed = 0

            for rgCrv in rgCrvs_toJoin:
                
                # Only NurbsCurves should have their control points translated when needed.
                if not isinstance(rgCrv, rg.NurbsCurve):
                    #print("Curve is a {}.".format(rgCrv.GetType().Name)
                    #sc.doc.Objects.AddCurve(rgCrv)
                    #sc.doc.Views.Redraw(); 1/0
                    continue # to next curve to include in JoinCurves.

                if rgCrv.PointAtStart.DistanceTo(joint) < sc.doc.ModelAbsoluteTolerance:
                    rgCrv.Points.SetPoint(
                        index=0,
                        point=joint,
                        weight=rgCrv.Points[0].Weight)
                    #sc.doc.Objects.AddCurve(rgCrv); sc.doc.Views.Redraw(); #1/0
                    iCt_EndsProcessed += 1
                elif rgCrv.PointAtEnd.DistanceTo(joint) < sc.doc.ModelAbsoluteTolerance:
                    rgCrv.Points.SetPoint(
                        index=rgCrv.Points.Count-1,
                        point=joint,
                        weight=rgCrv.Points[0].Weight)
                    #sc.doc.Objects.AddCurve(rgCrv); sc.doc.Views.Redraw(); #1/0
                    iCt_EndsProcessed += 1

                if iCt_EndsProcessed == 2:
                    break # for loop to next joint.

        return list(rg.Curve.JoinCurves(rgCrvs_toJoin))


    if bProjectCrvSegs:
        cs_toProj = rgC0.DuplicateSegments()
        if not cs_toProj:
            cs_toProj = [rgC0.DuplicateCurve()]
    else:
        cs_toProj = [rgC0.DuplicateCurve()]

    fTol_Proj_Total = sc.doc.ModelAbsoluteTolerance if fTol is None else fTol


    if bPostProcess:
        rgCs_Proj_1C0 = []

        for i, c_toProj in enumerate(cs_toProj):
            if bTryOnlyUniformCubicForNonlinear:
                projected_segs = projectThenProcessCrv_UniformCubicOnly(
                    c_toProj, fTol_Proj_Total)
            else:
                projected_segs = projectThenProcessCrv_NotOnly_UniformCubic(
                    c_toProj, fTol_Proj_Total)
            if not projected_segs:
                raise ValueError("No projection result!")
            rgCs_Proj_1C0.extend(projected_segs)

        if bDebug:
            reportCrvTypes(rgCs_Proj_1C0)
    else:
        rgCs_Proj_1C0 = []
        for c in cs_toProj:
            rgCs_Proj_1C0.extend(projectCrvToGeom(c, tolerance=fTol_Proj_Total))


    if not rgCs_Proj_1C0: return []

    #for rgC in rgCs_Proj_1C0: sc.doc.Objects.AddCurve(rgC)
    #sc.doc.Views.Redraw(); 1/0


    if len(rgCs_Proj_1C0) == 1:
        rgCs_Proj_Joined_1C0 = rgCs_Proj_1C0
    elif len(rgCs_Proj_1C0) > 1:
        # Since JoinCurves, as well as _Join, can skew the knot vector
        # (See )
        # join the curves by only moving control points as needed.
        rgCs_Proj_Joined_1C0 = joinCurves(rgCs_Proj_1C0)


        # This will clean up any short segments that had previous passed
        # the previous checks, but have been deformed to short segments
        # after JoinCurves.
        # This problem also occurs using _Project.
        for rgC in rgCs_Proj_Joined_1C0:
            rgC.RemoveShortSegments(fTol_Proj_Total)

    return rgCs_Proj_Joined_1C0


def processDocObjects(rhObjs_toProj, rhBreps, vect, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bProjectCrvSegs = getOpt('bProjectCrvSegs')
    bLoose = getOpt('bLoose')
    fTol = getOpt('fTol')
    bPostProcess = getOpt('bPostProcess')
    bTryOnlyUniformCubicForNonlinear = getOpt('bTryOnlyUniformCubicForNonlinear')
    bTryGetArcs = getOpt('bTryGetArcs')
    bAcceptRational = getOpt('bAcceptRational')
    bTryRebuildOthersUniform = getOpt('bTryRebuildOthersUniform')
    bDeleteInput = getOpt('bDeleteInput')
    iOutputLayer = getOpt('iOutputLayer')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def coerceBrepObject(rhObj):
        rdObj = rhinoscript.utility.coercerhinoobject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def processCurveObject():

        gCrvs_Out = []

        if bLoose:
            for rdB0, rgB0 in zip(rdB0s, rgB0s):
                if iOutputLayer == sOpts_OutputLayer.index('TargetObject'):
                    attr_Out.LayerIndex = rdB0.Attributes.LayerIndex

                rgCs_Proj_Joined_ThisC0_ThisB = projectCurve_Loose(
                    rgC0,
                    rgB0,
                    vect,
                    bProjectCrvSegs,
                    fTol,
                    bDebug)

                for rgC in rgCs_Proj_Joined_ThisC0_ThisB:
                    gCrv_Out = sc.doc.Objects.AddCurve(rgC, attributes=attr_Out)
                    if gCrv_Out != Guid.Empty:
                        gCrvs_Out.append(gCrv_Out)

            return gCrvs_Out
        else:
            # Tight.
            for rdB0, rgBs_1Face, rgPlanes in zip(rdB0s, rgB1s_1F_PerB0, rgPlanes_PerB0):
                if iOutputLayer == sOpts_OutputLayer.index('TargetObject'):
                    attr_Out.LayerIndex = rdB0.Attributes.LayerIndex

                rgCs_Proj_Joined_ThisC0_ThisB = projectCurve_Tight(
                    rgC0,
                    rgBs_1Face,
                    rgPlanes,
                    vect,
                    bProjectCrvSegs=bProjectCrvSegs,
                    fTol=fTol,
                    bPostProcess=bPostProcess,
                    bTryOnlyUniformCubicForNonlinear=bTryOnlyUniformCubicForNonlinear,
                    bTryGetArcs=bTryGetArcs,
                    bAcceptRational=bAcceptRational,
                    bTryRebuildOthersUniform=bTryRebuildOthersUniform,
                    bDebug=bDebug)

                for rgC in rgCs_Proj_Joined_ThisC0_ThisB:
                    gCrv_Out = sc.doc.Objects.AddCurve(rgC, attributes=attr_Out)
                    if gCrv_Out != Guid.Empty:
                        gCrvs_Out.append(gCrv_Out)

            return gCrvs_Out


    def processPointObject():

        gPts_Out = []

        for rdB, rgB in zip(rdB0s, rgB0s):
            if iOutputLayer == sOpts_OutputLayer.index('TargetObject'):
                attr_Out.LayerIndex = rdB.Attributes.LayerIndex

            pts_Proj = rg.Intersect.Intersection.ProjectPointsToBreps(
                breps=[rgB],
                points=[rgPt0.Location],
                direction=vect,
                tolerance=0.1*fTol)
            if not pts_Proj: continue

            for rgPt in pts_Proj:
                gPt_Out = sc.doc.Objects.AddPoint(rgPt, attributes=attr_Out)
                if gPt_Out != Guid.Empty:
                    gPts_Out.append(gPt_Out)

        return gPts_Out



    rgB0s = []
    rgB1s_1F_PerB0 = []
    rgPlanes_PerB0 = []
    rdB0s = []

    for rhB in rhBreps:
        rgB = coerceBrep(rhB)
        if not rgB.IsValid:
            print("Invalid Brep in input, but projection to it" \
                " will still be attempted." \
                "  Check results!")
        rgB0s.append(rgB)

        if not bLoose:
            rgB1s_1F_PerB0.append([])
            rgPlanes_PerB0.append([])

            for rgF in rgB.Faces:
                rgB_1F = rgF.DuplicateFace(duplicateMeshes=False)

                # Shrink face, otherwise Curve.ProjectToBrep
                # for some toroidal RevSurfaces may return None.
                # If this doesn't work, conversion from RevSurface
                # to NurbsSurface may be a solution.
                rgB_1F.Faces.ShrinkFaces()

                rgB1s_1F_PerB0[-1].append(rgB_1F)
                rgSrf = rgF.UnderlyingSurface()
                bSuccess, rgPlane = rgSrf.TryGetPlane(tolerance=(1.0/2**32))
                if bSuccess:
                    rgPlanes_PerB0[-1].append(rgPlane)
                else:
                    rgPlanes_PerB0[-1].append(None)

        rdB0s.append(coerceBrepObject(rhB))

    rgCs_Proj_All = []

    attr_Out = rd.ObjectAttributes()

    g_Out_All = []

    if iOutputLayer == sOpts_OutputLayer.index('Current'):
        attr_Out.LayerIndex = sc.doc.Layers.CurrentLayerIndex

    len_rhObjs0 = len(rhObjs_toProj)
    idxs_AtTenths = [int(round(0.1*i*len(rhObjs_toProj),0)) for i in range(10)]

    for iO, rhObj_toProj in enumerate(rhObjs_toProj):
        if sc.escape_test(False):
            print("Projecting interrupted by user.")
            return

        if len_rhObjs0 > 10:
            if iO in idxs_AtTenths:
                Rhino.RhinoApp.SetCommandPrompt("Processing at {:d}% of {} object ...".format(
                    int(100.0 * (iO+1) / len_rhObjs0), len_rhObjs0))
        elif len_rhObjs0 > 1:
            Rhino.RhinoApp.SetCommandPrompt(
                "Projecting {} of {} objects ({})".format(
                    iO+1, len_rhObjs0, rhinoscript.utility.coerceguid(rhObj_toProj)))
        else:
            Rhino.RhinoApp.SetCommandPrompt("Projecting object ...")


        rdObj0 = rhinoscript.utility.coercerhinoobject(rhObj_toProj)

        if iOutputLayer == sOpts_OutputLayer.index('Input'):
            attr_Out.LayerIndex = rdObj0.Attributes.LayerIndex

        if rdObj0.ObjectType == rd.ObjectType.Curve:
            rgC0 = rdObj0.Geometry
            g_Out = processCurveObject()
            if g_Out:
                g_Out_All.extend(g_Out)
            if bDeleteInput:
                sc.doc.Objects.Delete(item=rdObj0)
        elif rdObj0.ObjectType == rd.ObjectType.Brep:
            rgC0 = rhinoscript.utility.coercecurve(rhObj_toProj)
            g_Out = processCurveObject()
            if g_Out:
                g_Out_All.extend(g_Out)
                # BrepObject is not deleted.
        elif rdObj0.ObjectType == rd.ObjectType.Point:
            rgPt0 = rdObj0.Geometry
            g_Out = processPointObject()
            if g_Out:
                g_Out_All.extend(g_Out)
            if bDeleteInput:
                sc.doc.Objects.Delete(item=rdObj0)
        else:
            continue

    for rgB in rgB0s: rgB.Dispose()

    if not bLoose:
        for rgB1s_1F in rgB1s_1F_PerB0:
            for rgB in rgB1s_1F:
                rgB.Dispose()


    return g_Out_All


def getDirectionVector(iDirection):
    if sOpts_Direction[iDirection] == 'CPlaneX':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().XAxis
    if sOpts_Direction[iDirection] == 'CPlaneY':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().YAxis
    if sOpts_Direction[iDirection] == 'CPlaneZ':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    elif sOpts_Direction[iDirection] == 'WorldX':
        return rg.Vector3d.XAxis
    elif sOpts_Direction[iDirection] == 'WorldY':
        return rg.Vector3d.YAxis
    elif sOpts_Direction[iDirection] == 'WorldZ':
        return rg.Vector3d.ZAxis
    elif sOpts_Direction[iDirection] == 'View':
        return sc.doc.Views.ActiveView.ActiveViewport.GetCameraFrame()[1].ZAxis
    elif sOpts_Direction[iDirection] == 'Custom':
        return Opts.values['vectCustom']


def main():
    
    rc = getInput(
        [],
        "Select curves and points to project",
        rd.ObjectType.Curve | rd.ObjectType.Point)
    if rc is None: return
    objrefs_Crvs_Pts = rc[0]

    #for o in objrefs_Crvs_Pts:
    #    geomCompIdx = o.GeometryComponentIndex
    #    rdCompIdxType = geomCompIdx.ComponentIndexType
    #    if  geomCompIdx.ComponentIndexType == rg.ComponentIndexType.BrepEdge:
    #        o.Object().HighlightSubObject(geomCompIdx, highlight=True)
    #sc.doc.Views.Redraw()

    sc.doc.Objects.UnselectAll()

    rdCrvs_toHighlight = [o.Object() for o in objrefs_Crvs_Pts]

    print("Now, select breps and faces ...")

    rc = getInput(
        rdCrvs_toHighlight,
        "Select surfaces and polysurfaces to project onto",
        rd.ObjectType.Brep)
    if rc is None: return
    objrefs_Breps = rc[0]
    iDirection = rc[6]

    # Determine vector before main routine.
    vect = getDirectionVector(Opts.values['iDirection'])
    #if sOpts_Direction[iDirection] == 'CPlaneZ':
    #    vect = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    #elif sOpts_Direction[iDirection] == 'WorldZ':
    #    vect = rg.Vector3d.ZAxis
    #elif sOpts_Direction[iDirection] == 'View':
    #    rc = sc.doc.Views.ActiveView.ActiveViewport.GetCameraFrame()
    #    if not rc[0]: return
    #    vect = rc[1].ZAxis
    #elif sOpts_Direction[iDirection] == 'Custom':
    #    rc = rhinoscript.userinterface.GetLine(
    #        mode=1, point=None,
    #        message1='Projection direction',
    #        message3='Second direction point',
    #        )
    #    if not rc: return
    #    vect = rg.Vector3d(rc[1] - rc[0])

    if Opts.values['bDebug']:
        pass
    else:
        pass

    g_Res = processDocObjects(objrefs_Crvs_Pts, objrefs_Breps, vect)

    if Opts.values['bEcho']:
        if len(objrefs_Crvs_Pts) == len(g_Res):
            print("{} objects projected to same number of objects.".format(
                len(objrefs_Crvs_Pts)))
        else:
            print("{} objects projected to {} objects.".format(
                len(objrefs_Crvs_Pts), len(g_Res)))

    if g_Res:
        sc.doc.Objects.UnselectAll()
        sc.doc.Objects.Select(objectIds=g_Res)
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()