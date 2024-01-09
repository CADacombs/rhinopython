"""
This script will create draft curves silhouettes as obtained via _DraftAngleAnalysis,
except with curve simplification options and a more convenient execution.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums:
https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
171019-21: Created.
200509: Updated for V7 using newly added Silhouette.ComputeDraftCurve method.
200511,13: Added routine that simplifies curves.  Added some options.
200805: Modified pullDirection sign to match change in ComputeDraftCurve in Rhino 7.0.20217.3575, 8/4/2020.  Bug fix in printed output.
210223: Added bArcOK.  Bug fixes.
210315-17: Modified to check accuracy of output curves against input tolerances.  Added bDeg5.
210515: Now by default, will try to rebuild both degrees 3 and 5 within tolerance.
        Now, supports sub-selected BrepFaces as input.
220821-22: Refactored main V5&6 from V7 routines into separate functions.  Bug fixes.
240106-08: Removed code that only works in V5 & V6.  Added more direction options.
        Modified simplification routine.  Refactored.

TODO:
    Investigate computed curves that are more than input tolerance off of their faces. Replace with new pulled curve?
    Similar to the no-simplify routines, split simplify function into draft and no-draft functions.
    Add non-G2 splits to simplification functions.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import math


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iPullDir'; keys.append(key)
    listValues[key] = (
        'CPlanePosX',
        'CPlaneNegX',
        'CPlanePosY',
        'CPlaneNegY',
        'CPlanePosZ',
        'CPlaneNegZ',
        'WorldPosX',
        'WorldNegX',
        'WorldPosY',
        'WorldNegY',
        'WorldPosZ',
        'WorldNegZ',
        'View',
        'Custom',
        )
    values[key] = 4
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'vectCustom'; keys.append(key)
    values[key] = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDraftAngle'; keys.append(key)
    values[key] = 0.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(
        values[key], setLowerLimit=True, limit=1e-6)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fAngleTol'; keys.append(key)
    values[key] = 0.5 * sc.doc.ModelAngleToleranceDegrees
    riOpts[key] = ri.Custom.OptionDouble(
        values[key], setLowerLimit=True, limit=1e-6)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSimplifyRes'; keys.append(key)
    values[key] = True
    names[key] = 'AttemptSimplifyResult'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bArcOK'; keys.append(key)
    values[key] = True
    names[key] = 'Arcs'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'OK')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg2'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg3'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg5'; keys.append(key)
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

        if key == 'fDraftAngle':
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Why is key, {}, here?  Value was not set or sticky-saved.".format(key))
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select brep(s/faces)")

    go.GeometryFilter = rd.ObjectType.Brep

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('iPullDir')
        addOption('fDraftAngle')
        addOption('fDistTol')
        addOption('fAngleTol')
        addOption('bSimplifyRes')
        if Opts.values['bSimplifyRes']:
            addOption('bArcOK')
            addOption('bDeg2')
            addOption('bDeg3')
            addOption('bDeg5')
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
            key = 'fDraftAngle'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        if go.OptionIndex() == idxs_Opt['iPullDir']:
            Opts.values['iPullDir'] = (
                    go.Option().CurrentListOptionIndex)

            if Opts.listValues['iPullDir'][go.Option().CurrentListOptionIndex] == 'Custom':
                rc = rs.GetLine(
                    mode=1, point=None,
                    message1="Point for pull direction from",
                    message3="Point for pull direction to",
                    )
                if not rc:
                    print("Failed to get pull direction.")
                    continue

                key = 'vectCustom'
                Opts.values[key] = rg.Vector3d(rc[1] - rc[0])
                Opts.values[key].Unitize()
                sc.sticky[Opts.stickyKeys[key]] = Opts.values[key]
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getDirectionVector(iPullDir):
    sDirs = Opts.listValues['iPullDir']
    if sDirs[iPullDir] == 'CPlanePosX':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().XAxis
    if sDirs[iPullDir] == 'CPlaneNegX':
        return -sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().XAxis
    if sDirs[iPullDir] == 'CPlanePosY':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().YAxis
    if sDirs[iPullDir] == 'CPlaneNegY':
        return -sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().YAxis
    if sDirs[iPullDir] == 'CPlanePosZ':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    if sDirs[iPullDir] == 'CPlaneNegZ':
        return -sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    elif sDirs[iPullDir] == 'WorldPosX':
        return rg.Vector3d.XAxis
    elif sDirs[iPullDir] == 'WorldNegX':
        return -rg.Vector3d.XAxis
    elif sDirs[iPullDir] == 'WorldPosY':
        return rg.Vector3d.YAxis
    elif sDirs[iPullDir] == 'WorldNegY':
        return -rg.Vector3d.YAxis
    elif sDirs[iPullDir] == 'WorldPosZ':
        return rg.Vector3d.ZAxis
    elif sDirs[iPullDir] == 'WorldNegZ':
        return -rg.Vector3d.ZAxis
    elif sDirs[iPullDir] == 'View':
        return sc.doc.Views.ActiveView.ActiveViewport.GetCameraFrame()[1].ZAxis
    elif sDirs[iPullDir] == 'Custom':
        return Opts.values['vectCustom']
    raise Exception("Error in getDirectionVector.")


def _tangentAngleDifference(cA, cB, bDebug=False):
    tanDiffAtStart = rg.Vector3d.VectorAngle(cA.TangentAtStart, cB.TangentAtStart)
    tanDiffAtEnd = rg.Vector3d.VectorAngle(cA.TangentAtEnd, cB.TangentAtEnd)
    if bDebug:
        sEval='Rhino.RhinoMath.ToDegrees(tanDiffAtStart)'; print(sEval+':',eval(sEval))
        sEval='Rhino.RhinoMath.ToDegrees(tanDiffAtEnd)'; print(sEval+':',eval(sEval))


def _areEpsilonEqual(a, b, epsilon):
    # This is a relative comparison.
    delta = abs(a - b)
    fRelComp = delta / max(abs(a), abs(b))
    return fRelComp < epsilon


def isUniformNurbsCurve(nc):
    if not isinstance(nc, rg.NurbsCurve):
        return

    if nc.SpanCount == 1:
        return True

    if nc.Knots.KnotStyle in (rg.KnotStyle.QuasiUniform, rg.KnotStyle.Uniform):
        return True

    return False

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
        if not _areEpsilonEqual(
            span0, nc.Knots[i+1] - nc.Knots[i],
            epsilon=(1.0/2**32)):
                return False

    return True


def is_NurbsCrv_internally_G2_continuous(nc, fAngleTol, bDebug=False):
    if not isinstance(nc, rg.NurbsCurve):
        raise ValueError("{} passed to is_NurbsCrv_internally_G2_continuous.".format(nc))

    if nc.Degree >= 3 and isUniformNurbsCurve(nc):
        return True

    bFoundG2_discont, t = nc.GetNextDiscontinuity(
        continuityType=rg.Continuity.G2_locus_continuous if nc.IsPeriodic else rg.Continuity.G2_continuous,
        t0=nc.Domain.T0,
        t1=nc.Domain.T1)
    if bDebug: sEval='bFoundG2_discont'; print(sEval+':',eval(sEval))
    if bFoundG2_discont:
        return False

    # Testing G1 separately to apply a custom tolerance.
    # If the discontinuity type is G2 instead, is the previous check necessary?
    bFoundG1_discont, t = nc.GetNextDiscontinuity(
        continuityType=rg.Continuity.G1_locus_continuous if nc.IsPeriodic else rg.Continuity.G1_continuous,
        t0=nc.Domain.T0,
        t1=nc.Domain.T1,
        cosAngleTolerance=math.cos(sc.doc.ModelAngleToleranceRadians),
        curvatureTolerance=Rhino.RhinoMath.SqrtEpsilon)
    if bFoundG1_discont:
        print(t)
        return False

    return True


def split_NurbsCrv_at_G2_discontinuities(nc_In, bDebug=False):
    if not isinstance(nc_In, rg.NurbsCurve):
        raise ValueError("{} passed to is_NurbsCrv_internally_G2_continuous.".format(nc_In))

    cosAngleTolerance = math.cos(sc.doc.ModelAngleToleranceRadians)

    t0 = nc_In.Domain.T0 # Start parameter of check.  Is itself ignored.

    ts_Found = []

    while True:
        sc.escape_test()

        bFoundG1_discont, t = nc_In.GetNextDiscontinuity(
            continuityType=rg.Continuity.G1_locus_continuous if nc_In.IsPeriodic else rg.Continuity.G1_continuous,
            t0=t0,
            t1=nc_In.Domain.T1,
            cosAngleTolerance=cosAngleTolerance,
            curvatureTolerance=Rhino.RhinoMath.SqrtEpsilon)

        if not bFoundG1_discont:
            break # out of while loop.

        ts_Found.append(t)

        t0 = t

    if len(ts_Found) == 0:
        return

    split = nc_In.Split(ts_Found)

    if len(split) < 2:
        raise Exception("Split produced {} curves. Should have been {}.".format(len(ts_Found)))

    joined = rg.Curve.JoinCurves(split)

    if len(joined) != 1:
        raise Exception("JoinCurves produced {} curves. Should have been 1.".format(len(joined)))

    for c in split: c.Dispose()

    return joined


def split_NurbsCrvs_at_G2_discontinuities(rgCrvs_In, bDebug=False):
    """
    Will return new list with splits and duplicates,
    so it is up to the calling code to Dispose the curves in the original list.
    """
    rgCs_Out = []
    for c in rgCrvs_In:
        if not isinstance(c, rg.NurbsCurve):
            rgCs_Out.append(c.DuplicateCurve())
            continue

        rc = split_NurbsCrv_at_G2_discontinuities(c, bDebug=bDebug)

        if rc is None:
            rgCs_Out.append(c.DuplicateCurve())
            continue

        rgCs_Out.extend(rc)

    return rgCs_Out


def rebuild_Bezier(nc_In, iDegs, tolerance, bDebug=False):
    """
    More strict tolerance used on degree 2 since that Bezier doesn't allow
    tangency matching of both ends when they are not already aligned.

    Returns on success: rg.NurbsCurve, float(deviation)
    Returns on fail: None
    """
    if bDebug: print("rebuild_Bezier with tolerance={}:".format(tolerance))

    for iDeg in iDegs:
        pointCount = iDeg + 1

        rebuilt = nc_In.Rebuild(
            pointCount=pointCount,
            degree=iDeg,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            nc_In, rebuilt, tolerance=0.1*tolerance)[:2]

        if bSuccess:
            if iDeg==2:
                if fDistMax < 1e-6:
                    if bDebug:
                        print("Rebuilt within {}  Deg:{}  PtCt:{}".format(
                            fDistMax, iDeg, pointCount))
                        _tangentAngleDifference(nc_In, rebuilt, bDebug)
                    return rebuilt, fDistMax
            elif fDistMax <= tolerance:
                if bDebug:
                    print("Rebuilt within {}  Deg:{}  PtCt:{}".format(
                        fDistMax, iDeg, pointCount))
                    _tangentAngleDifference(nc_In, rebuilt, bDebug)
                return rebuilt, fDistMax

        if bDebug:
            print("Rebuild requires {}  Deg:{}  PtCt:{}".format(
                fDistMax, iDeg, pointCount))
        rebuilt.Dispose()


def rebuild_MultiSpan(nc_In, iDegs, tolerance, bDebug=False):
    """
    Will not bother rebuilding to degree-2 since multiple spans are divided by
    G1-likely knots.
    """
    if bDebug: print("rebuild_MultiSpan with tolerance={}:".format(tolerance))

    iCt_MaxSpans = nc_In.SpanCount

    #iCt_MaxCp = int(round(nc.GetLength() / (100.0 * sc.doc.ModelAbsoluteTolerance)))
    #if bDebug: sEval='iCt_MaxCp'; print(sEval+':',eval(sEval))


    # Rebuild at maximum control point count.
    #pointCount = iCt_MaxCp

    iDegs_ForRebuildSearch = []
    rebuilts_LastSuccess = [] # per iDeg.
    fDevs = [] # per iDeg.

    for iDeg in iDegs:

        if iDeg < 3: continue

        pointCount = iDeg + iCt_MaxSpans

        rebuilt = nc_In.Rebuild(
            pointCount=pointCount,
            degree=iDeg,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            nc_In,
            rebuilt,
            tolerance=0.1*tolerance)[:2]

        if bSuccess and fDistMax <= tolerance:
            iDegs_ForRebuildSearch.append(iDeg)
            rebuilts_LastSuccess.append(rebuilt)
            fDevs.append(fDistMax)
            if bDebug:
                print("Rebuilt within {}  Deg:{}  PtCt:{}".format(
                    tolerance, iDeg, pointCount))
        else:
            rebuilt.Dispose()


    if not iDegs_ForRebuildSearch:
        if bDebug:
            print("Not rebuilt within {} at max. span ct. of {}, so quitting rebuilding.".format(
                tolerance, iCt_MaxSpans))
        return


    if bDebug: print("Binary search.")

    for i, iDeg in enumerate(iDegs_ForRebuildSearch):

        iCt_MaxCp = iDeg + iCt_MaxSpans

        iCts_Cps_Tried = [iDeg + 1, iCt_MaxCp]

        iCt_Cp_Try = (iCt_MaxCp + iDeg + 1) // 2

        while iCt_Cp_Try not in iCts_Cps_Tried:
            sc.escape_test()


            rebuilt = nc_In.Rebuild(
                pointCount=iCt_Cp_Try,
                degree=iDeg,
                preserveTangents=True)

            bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
                nc_In, rebuilt, tolerance=0.1*tolerance)[:2]

            if bDebug:
                print("Degree:{}  CPtCt:{}  fDistMax:{}  WithinTol:{}".format(
                    iDeg,
                    iCt_Cp_Try,
                    fDistMax if bSuccess else bSuccess,
                    fDistMax <= tolerance if bSuccess else bSuccess
                    ))


            iCts_Cps_Tried.append(iCt_Cp_Try)
            iCts_Cps_Tried.sort()

            if bSuccess and fDistMax <= tolerance:
                rebuilts_LastSuccess[i].Dispose()
                rebuilts_LastSuccess[i] = rebuilt
                fDevs[i] = fDistMax
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

    if len(rebuilts_LastSuccess) == 1:
        _tangentAngleDifference(nc_In, rebuilts_LastSuccess[0], bDebug)
        return rebuilts_LastSuccess[0], fDevs[0]

    # TODO: Change the winner to that with minimum deviation instead of least CPs?
    iCts_Pts = [nc_In.Points.Count for nc_In in rebuilts_LastSuccess]
    idx_Winner = iCts_Pts.index(min(iCts_Pts))
    _tangentAngleDifference(nc_In, rebuilts_LastSuccess[idx_Winner], bDebug)
    return rebuilts_LastSuccess[idx_Winner], fDevs[idx_Winner]


def removeMultiKnots(nc_In, iDegs, tolerance, bDebug):
    """
    Returns on success: rg.NurbsCurve, float(deviation)
    Returns on fail: None


    minimumMultiplicity: Remove knots with multiplicity > minimumKnotMultiplicity.
    maximumMultiplicity: Remove knots with multiplicity < maximumKnotMultiplicity.
    """

    if bDebug: print("removeMultiKnots:".format(tolerance))

    for iDeg in iDegs:
        if iDeg < nc_In.Degree: continue

        nc_WIP = nc_In.DuplicateCurve()

        if iDeg > nc_In.Degree:
            nc_WIP.IncreaseDegree(iDeg)

        for minimumMultiplicity in range(1, nc_WIP.Degree-2+1):

            iCt_KnotsRemoved = nc_WIP.Knots.RemoveMultipleKnots(
                minimumMultiplicity=minimumMultiplicity,
                maximumMultiplicity=nc_WIP.Degree+1,
                tolerance=1.0)
            if not iCt_KnotsRemoved:
                if bDebug: print("RemoveMultipleKnots failed to remove any knots.")
                continue # to next minimumMultiplicity.

            bSuccess, fDev = rg.Curve.GetDistancesBetweenCurves(
                nc_WIP, nc_In, tolerance=0.1*tolerance)[:2]
            if not bSuccess:
                continue # to next minimumMultiplicity.

            if bDebug: sEval='fDev'; print(sEval+':',eval(sEval))
            if fDev > tolerance:
                continue # to next minimumMultiplicity.

            _tangentAngleDifference(nc_In, nc_WIP, bDebug)
            return nc_WIP, fDev

        nc_WIP.Dispose()


def simplifyCurve(rgCrv_In, bArcOK=False, iDegs=(2,3,5), fTol=None, fAngleTol=None, bDebug=False):
    """
    Simplify includes:
        Uniform knot vectors.

    Returns on success: rg.ArcCurve, rgLineCurve, or simplified rg.NurbsCurve.
    Returns on fail: None
    """
    if bDebug: print("simplifyCurve with fTol={}:".format(fTol))


    rgCrvs_Out = []

    fTol_Min = 1e-6

    if fTol is None:
        fTol = fTol_Min

    if fAngleTol is None:
        fAngleTol = sc.doc.ModelAngleToleranceDegrees

    fTol_MaxRad = 1e3


    if isinstance(rgCrv_In, (rg.PolylineCurve, rg.PolyCurve)):
        raise Exception("Input is {}.".format(rgCrv_In.GetType().Name))

    if isinstance(rgCrv_In, (rg.ArcCurve, rg.LineCurve)):
        print("Curve is already a {}.  How did that happen?".format(rgCrv_In.GetType().Name))
        return rgCrv_In.DuplicateCurve()

    nc_In = rgCrv_In

    if not nc_In.IsClosed:
        if nc_In.IsLinear(fTol_Min):
            return rg.LineCurve(
                nc_In.PointAtStart,
                nc_In.PointAtEnd)

    if bArcOK:
        if nc_In.IsClosed:
            b, circle = nc_In.TryGetCircle(fTol_Min)
            if b:
                if circle.Radius <= fTol_MaxRad:
                    return rg.ArcCurve(circle)
        else:
            b, arc = nc_In.TryGetArc(fTol_Min)
            if b:
                if arc.Radius <= fTol_MaxRad:
                    return rg.ArcCurve(arc)

        ## Just keep rational result.
        ## TODO: Remove this?
        #if nc_In.IsRational:
        #    return nc_In.Duplicate()

    if isUniformNurbsCurve(nc_In):
        return nc_In.Duplicate()

    rc = rebuild_Bezier(
        nc_In,
        iDegs=iDegs,
        tolerance=fTol,
        bDebug=bDebug)

    if rc is not None:
        return rc[0]

    if is_NurbsCrv_internally_G2_continuous(nc_In, fAngleTol, bDebug):
        rc = removeMultiKnots(nc_In, iDegs, fTol, bDebug)
        if rc: return rc[0]

    rc = rebuild_MultiSpan(
        nc_In,
        iDegs=iDegs,
        tolerance=fTol,
        bDebug=bDebug)

    if rc:
        return rc[0]

    print("TODO: Split at non-G2 knots.")


def computeMaxDistFromFace(rgFace, rgC_on_face, tolerance):
    tolerance=1e-7
    pulled = rgC_on_face.PullToBrepFace(face=rgFace, tolerance=tolerance)
    if len(pulled) == 1:
        pulled = pulled[0]
    else:
        pullback = rgFace.Pullback(rgC_on_face, tolerance=tolerance)
        sEval='pullback.GetType().Name'; print(sEval+':',eval(sEval))
        print(pullback.GetType().Name)
        pushup = rgFace.Pushup(pullback, tolerance=tolerance)
        sEval='pushup.GetType().Name'; print(sEval+':',eval(sEval))
        pulled = pushup
#        sc.doc.Objects.AddCurve(rgC_on_face); sc.doc.Views.Redraw()
#        raise Exception("PullToBrepFace produced {} curves, but should have been 1.".format(
#            len(pulled)))
    bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
        rgC_on_face, pulled, tolerance=tolerance)[:2]
    if not bSuccess:
        raise Exception("GetDistancesBetweenCurves failed.")

    if fDistMax > tolerance:
        sc.doc.Objects.AddCurve(rgC_on_face)
        sc.doc.Objects.AddCurve(pulled)
        sc.doc.Views.Redraw()
        1/0

    pulled.Dispose()

    return fDistMax


def createSilhouetteCurves_Simplify(rgGeomForSilh, **kwargs):
    """
    V7+ function due to use of Silhouette.ComputeDraftCurve.

    Parameters:
        rgGeomForSilh
        kwargs:
            iPullDir
            fDraftAngle
            fDistTol
            fAngleTol
            bDeg3
            bDeg5
            bArcOK
            bEcho
            bDebug
    """

    if Rhino.RhinoApp.ExeVersion < 7:
        print("This function requires Rhino 7 or higher.")
        return

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPullDir = getOpt('iPullDir')
    fDraftAngle = getOpt('fDraftAngle')
    fDistTol = getOpt('fDistTol')
    fAngleTol = getOpt('fAngleTol')
    bArcOK = getOpt('bArcOK')
    bDeg2 = getOpt('bDeg2')
    bDeg3 = getOpt('bDeg3')
    bDeg5 = getOpt('bDeg5')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    iDegs = []
    if bDeg2:
        iDegs.append(2)
    if bDeg3:
        iDegs.append(3)
    if bDeg5:
        iDegs.append(5)
    if not iDegs and bEcho:
        print("No degrees for rebuild.")


    def computeDraftCurvesWithSimplify(rgGeomForSilh, bArcOK):
        if bDebug: print("computeDraftCurvesWithSimplify:")

        tolerance = 0.5 * fDistTol

        if bDebug:
            sEval='rgGeomForSilh'; print(sEval+':',eval(sEval))
            sEval='draftAngle'; print(sEval+':',eval(sEval))
            sEval='pullDirection'; print(sEval+':',eval(sEval))
            sEval='tolerance'; print(sEval+':',eval(sEval))
            sEval='angleToleranceRadians'; print(sEval+':',eval(sEval))

        silhouettes = rg.Silhouette.ComputeDraftCurve(
            geometry=rgGeomForSilh,
            draftAngle=draftAngle,
            pullDirection=pullDirection,
            tolerance=tolerance,
            angleToleranceRadians=angleToleranceRadians)

        if not silhouettes:
            if bDebug: print("Silhouette.ComputeDraftCurve produced no results.")
            return
        elif bDebug:
            print("Silhouette.ComputeDraftCurve produced {} results.".format(
                len(silhouettes)))


        rgCs_Res_FirstTol = []
        bNonG2KnotsFound = False
        fLengths_Totals = [0.0]
        for silh in silhouettes:
            if silh.Curve is None: continue
            if silh.SilhouetteType != rg.SilhouetteType.DraftCurve:
                if bDebug: sEval='silh.SilhouetteType'; print(sEval+':',eval(sEval))
            rgC_Silh = silh.Curve
            fLengths_Totals[-1] += rgC_Silh.GetLength()

            rgC_Simplified = simplifyCurve(
                rgC_Silh,
                bArcOK=bArcOK,
                iDegs=iDegs,
                fTol=max((fDistTol-tolerance, 1e-6)),
                fAngleTol=None,
                bDebug=bDebug)

            if rgC_Simplified is None:
                rgCs_Res_FirstTol.append(rgC_Silh)
                bNonG2KnotsFound = True
            else:
                rgCs_Res_FirstTol.append(rgC_Simplified)
                rgC_Silh.Dispose()

        if bDebug: sEval="fLengths_Totals[-1]"; print(sEval+':',eval(sEval))

        if not bNonG2KnotsFound:
            if bDebug:
                print("All curves are uniform at default tolerance, {}.".format(
                    fDistTol))
            return rgCs_Res_FirstTol


        # ComputeDraftCurve to a different tolerance.

        for pow in range(2, 4+1):
            fTol_forCDC_WIP = fDistTol / (2**pow)
            if bDebug: sEval='fTol_forCDC_WIP'; print(sEval+':',eval(sEval))

            silhouettes = rg.Silhouette.ComputeDraftCurve(
                geometry=rgGeomForSilh,
                draftAngle=draftAngle,
                pullDirection=pullDirection,
                tolerance=fTol_forCDC_WIP,
                angleToleranceRadians=angleToleranceRadians)

            if not silhouettes:
                if bDebug: print("Silhouette.ComputeDraftCurve produced no results.")
                return

            rgCs_Res_WIPTol = []
            fLengths_Totals.append(0.0)
            for silh in silhouettes:
                if silh.Curve is None: continue
                if silh.SilhouetteType != rg.SilhouetteType.DraftCurve:
                    if bDebug: sEval='silh.SilhouetteType'; print(sEval+':',eval(sEval))
                rgC_Silh = silh.Curve
                fLengths_Totals[-1] += rgC_Silh.GetLength()

                rgC_Simplified = simplifyCurve(
                    rgC_Silh,
                    bArcOK=bArcOK,
                    iDegs=iDegs,
                    fTol=fDistTol - fTol_forCDC_WIP,
                    fAngleTol=None,
                    bDebug=bDebug)

                if rgC_Simplified is None:
                    for c in rgCs_Res_WIPTol: c.Dispose()
                    break # to next tolerance.

                rgCs_Res_WIPTol.append(rgC_Simplified)
                rgC_Silh.Dispose()
            else:
                # Success.
                if bDebug: sEval="fLengths_Totals[-1]"; print(sEval+':',eval(sEval))
                for c in rgCs_Res_FirstTol: c.Dispose()
                return rgCs_Res_WIPTol

        if bDebug:
            sEval="fLengths_Totals[-1]"; print(sEval+':',eval(sEval))
            sEval="min(fLengths_Totals)"; print(sEval+':',eval(sEval))
            sEval="max(fLengths_Totals)"; print(sEval+':',eval(sEval))

        # Fail.
        return rgCs_Res_FirstTol


    draftAngle = Rhino.RhinoMath.ToRadians(fDraftAngle)
    pullDirection = _getDirectionVector(Opts.values['iPullDir'])
    angleToleranceRadians = Rhino.RhinoMath.ToRadians(fAngleTol)


    if bDebug: print("Using V7's Silhouette.ComputeDraftCurve.")


    if isinstance(rgGeomForSilh, rg.Brep):
        rgFs = rgGeomForSilh.Faces
    elif isinstance(rgGeomForSilh, rg.BrepFace):
        rgFs = [rgGeomForSilh]
    else:
        raise ValueError("{} is not supported for this script.".format(rgGeomForSilh.GetType().Name))




    rgCs_Res = []

    for rgF in rgFs:
        sc.escape_test()
        rgCs_Res_perFace = computeDraftCurvesWithSimplify(rgF, bArcOK)
        if rgCs_Res_perFace:
            rgCs_Res.extend(rgCs_Res_perFace)

            for c in rgCs_Res_perFace:
                fDev = computeMaxDistFromFace(rgF, c, fDistTol)

                if fDev > fDistTol:
                    print("Tolerance is {}, but curve is {} off of face.".format(
                        fDistTol, fDev))


    return rgCs_Res


def createSilhouetteCurves_NoDraft_NoSimplify(rgGeomForSilh, **kwargs):
    """
    V7+ function due to use of Silhouette.ComputeDraftCurve.

    ComputeDraftCurve is used with Compute since it has been noticed that they
    sometimes/often? produce different results.

    Parameters:
        rgGeomForSilh
        kwargs:
            iPullDir
            fDistTol
            fAngleTol
            bEcho
            bDebug
    """

    if Rhino.RhinoApp.ExeVersion < 7:
        print("This function requires Rhino 7 or higher.")
        return

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPullDir = getOpt('iPullDir')
    fDistTol = getOpt('fDistTol')
    fAngleTol = getOpt('fAngleTol')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    pullDirection = _getDirectionVector(Opts.values['iPullDir'])
    angleToleranceRadians = Rhino.RhinoMath.ToRadians(fAngleTol)


    if bDebug:
        sEval='pullDirection'; print(sEval+':',eval(sEval))
        sEval='fDistTol'; print(sEval+':',eval(sEval))
        sEval='fAngleTol'; print(sEval+':',eval(sEval))
        sEval='angleToleranceRadians'; print(sEval+':',eval(sEval))


    if isinstance(rgGeomForSilh, rg.Brep):
        rgFs = rgGeomForSilh.Faces
    elif isinstance(rgGeomForSilh, rg.BrepFace):
        rgFs = [rgGeomForSilh]
    else:
        raise ValueError("{} is not supported for this script.".format(rgGeomForSilh.GetType().Name))


    rgCs_Res_All = []

    for rgF in rgFs:
        silhouettes = rg.Silhouette.ComputeDraftCurve(
            geometry=rgF,
            draftAngle=0.0,
            pullDirection=pullDirection,
            tolerance=fDistTol,
            angleToleranceRadians=angleToleranceRadians)

        if not silhouettes:
            if bDebug: print("Silhouette.ComputeDraftCurve produced no results.")
            continue

        rgCs_Res_DraftCompute = []
        for silh in silhouettes:
            if silh.Curve is None:
                continue
            if silh.SilhouetteType != rg.SilhouetteType.Tangent:
                sEval="silh.SilhouetteType"; print(sEval+':',eval(sEval))
                continue
            rgCs_Res_DraftCompute.append(silh.Curve)

        silhouettes = rg.Silhouette.Compute(
            geometry=rgF,
            silhouetteType=rg.SilhouetteType.Tangent,
            parallelCameraDirection=pullDirection,
            tolerance=fDistTol,
            angleToleranceRadians=angleToleranceRadians)

        if not silhouettes:
            if bDebug: print("Silhouette.Compute produced no results.")
            continue

        rgCs_Res_NoDraftCompute = []
        for silh in silhouettes:
            if silh.Curve is None:
                continue
            if silh.SilhouetteType != rg.SilhouetteType.Tangent:
                continue
            rgCs_Res_NoDraftCompute.append(silh.Curve)

        if len(rgCs_Res_DraftCompute) != len(rgCs_Res_NoDraftCompute):
            if bDebug:
                sEval='len(rgCs_Res_DraftCompute)'; print(sEval+':',eval(sEval))
                sEval='len(rgCs_Res_NoDraftCompute)'; print(sEval+':',eval(sEval))

            if len(rgCs_Res_DraftCompute) == 0:
                rgCs_Res_All.extend(rgCs_Res_NoDraftCompute)
            elif len(rgCs_Res_NoDraftCompute) == 0:
                rgCs_Res_All.extend(rgCs_Res_DraftCompute)
            else:
                #for c in rgCs_Res_DraftCompute: sc.doc.Objects.AddCurve(c)
                #for c in rgCs_Res_NoDraftCompute: sc.doc.Objects.AddCurve(c)
                #sc.doc.Views.Redraw(); 1/0
                print("So adding both sets of curves.")
                rgCs_Res_All.extend(rgCs_Res_DraftCompute)
                rgCs_Res_All.extend(rgCs_Res_NoDraftCompute)
        else:
            for i in range(len(rgCs_Res_DraftCompute)):
                cDraft = rgCs_Res_DraftCompute[i]
                cNoDraft = rgCs_Res_NoDraftCompute[i]
                fDev_DraftCompute = computeMaxDistFromFace(rgF, cDraft, fDistTol)
                fDev_NoDraftCompute = computeMaxDistFromFace(rgF, cNoDraft, fDistTol)
                if bDebug:
                    sEval='fDev_DraftCompute'; print(sEval+':',eval(sEval))
                    sEval='fDev_NoDraftCompute'; print(sEval+':',eval(sEval))
                if fDev_DraftCompute == fDev_NoDraftCompute:
                    print("Identical curves?")
                    rgCs_Res_All.append(cNoDraft)
                elif fDev_DraftCompute > fDev_NoDraftCompute:
                    rgCs_Res_All.append(cNoDraft)
                else:
                    rgCs_Res_All.append(cDraft)

    rgCs_Out = split_NurbsCrvs_at_G2_discontinuities(rgCs_Res_All, bDebug=bDebug)
    for c in rgCs_Res_All: c.Dispose()

    return rgCs_Out


def createSilhouetteCurves_Draft_NoSimplify(rgGeomForSilh, **kwargs):
    """
    V7+ function due to use of Silhouette.ComputeDraftCurve.

    Parameters:
        rgGeomForSilh
        kwargs:
            iPullDir
            fDraftAngle
            fDistTol
            fAngleTol
            bEcho
            bDebug
    """

    if Rhino.RhinoApp.ExeVersion < 7:
        print("This function requires Rhino 7 or higher.")
        return

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPullDir = getOpt('iPullDir')
    fDraftAngle = getOpt('fDraftAngle')
    fDistTol = getOpt('fDistTol')
    fAngleTol = getOpt('fAngleTol')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    draftAngleRadians = Rhino.RhinoMath.ToRadians(fDraftAngle)
    pullDirection = _getDirectionVector(Opts.values['iPullDir'])
    angleToleranceRadians = Rhino.RhinoMath.ToRadians(fAngleTol)


    if bDebug:
        sEval='fDraftAngle'; print(sEval+':',eval(sEval))
        sEval='draftAngleRadians'; print(sEval+':',eval(sEval))
        sEval='pullDirection'; print(sEval+':',eval(sEval))
        sEval='fDistTol'; print(sEval+':',eval(sEval))
        sEval='fAngleTol'; print(sEval+':',eval(sEval))
        sEval='angleToleranceRadians'; print(sEval+':',eval(sEval))


    if isinstance(rgGeomForSilh, rg.Brep):
        rgFs = rgGeomForSilh.Faces
    elif isinstance(rgGeomForSilh, rg.BrepFace):
        rgFs = [rgGeomForSilh]
    else:
        raise ValueError("{} is not supported for this script.".format(rgGeomForSilh.GetType().Name))


    rgCs_Res_All = []

    for rgF in rgFs:
        silhouettes = rg.Silhouette.ComputeDraftCurve(
            geometry=rgF,
            draftAngle=draftAngleRadians,
            pullDirection=pullDirection,
            tolerance=fDistTol,
            angleToleranceRadians=angleToleranceRadians)

        if not silhouettes:
            if bDebug: print("Silhouette.ComputeDraftCurve produced no results.")
            continue

        for silh in silhouettes:
            if silh.Curve is None:
                continue
            if silh.SilhouetteType != rg.SilhouetteType.DraftCurve:
                continue
            rgCs_Res_All.append(silh.Curve)

    rgCs_Out = split_NurbsCrvs_at_G2_discontinuities(rgCs_Res_All, bDebug=bDebug)
    for c in rgCs_Res_All: c.Dispose()

    return rgCs_Out


def processDocObjects(rhBreps, **kwargs):
    """
    Parameters:
        rhBreps: list(objref or GUID)

    Returns: list(GUID) of CurveObjects
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPullDir = getOpt('iPullDir')
    fDraftAngle = getOpt('fDraftAngle')
    fDistTol = getOpt('fDistTol')
    fAngleTol = getOpt('fAngleTol')
    bSimplifyRes = getOpt('bSimplifyRes')
    bArcOK = getOpt('bArcOK')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getGeomForSilh(rhObj):
        if isinstance(rhObj, (rg.Brep, rg.BrepFace)):
            return rhObj

        if isinstance(rhObj, rd.ObjRef):
            #print(rhObj.GeometryComponentIndex.ComponentIndexType)
            rdObj = rhObj.Object()
            rgObj = rhObj.Geometry()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj)
            rgObj = rdObj.Geometry
        else:
            return

        if isinstance(rgObj, (rg.Brep, rg.BrepFace)):
            return rgObj

        print("{} {} passed to processBrepObjects.".format(rdObj, rgObj))


    rgBreps1 = []
    sLogs = []
    gBreps1 = []


    gCs_Out_All = []
    sCrvDescrs_All = []

    for rhBrep0 in rhBreps:
        rc = getGeomForSilh(rhBrep0)
        if not rc:
            raise ValueError(
                "{} passed to processBrepObjects." \
                    "  Only ObjRef, DocObject, Geometry, or an iterable combination" \
                    " are accepted input.")

        rdObj_In = rc

        if bSimplifyRes:
            rgCs_Res = createSilhouetteCurves_Simplify(
                rdObj_In,
                iPullDir=iPullDir,
                fDraftAngle=fDraftAngle,
                fDistTol=fDistTol,
                fAngleTol=fAngleTol,
                bArcOK=bArcOK,
                bEcho=bEcho,
                bDebug=bDebug)
        else:
            if fDraftAngle == 0.0:
                rgCs_Res = createSilhouetteCurves_NoDraft_NoSimplify(
                    rdObj_In,
                    iPullDir=iPullDir,
                    fDistTol=fDistTol,
                    fAngleTol=fAngleTol,
                    bEcho=bEcho,
                    bDebug=bDebug)
            else:
                rgCs_Res = createSilhouetteCurves_Draft_NoSimplify(
                    rdObj_In,
                    iPullDir=iPullDir,
                    fDraftAngle=fDraftAngle,
                    fDistTol=fDistTol,
                    fAngleTol=fAngleTol,
                    bEcho=bEcho,
                    bDebug=bDebug)

        if not rgCs_Res:
            continue

        for c in rgCs_Res:
            gC_Out = sc.doc.Objects.AddCurve(c)
            if gC_Out == gC_Out.Empty: continue

            gCs_Out_All.append(gC_Out)

            if isinstance(c, rg.NurbsCurve):
                if c.SpanCount == 1:
                    sCrvDescrs_All.append('Bezier')
                    continue
                sCrvDescrs_All.append("{} NURBS".format(c.Knots.KnotStyle))
            else:
                sCrvDescrs_All.append(c.GetType().Name)

    return gCs_Out_All, sCrvDescrs_All


def main():

    if Rhino.RhinoApp.ExeVersion < 7:
        print("This script requires Rhino 7 or higher.")
        return

    objrefs = getInput()
    if objrefs is None: return

    iPullDir = Opts.values['iPullDir']
    fDraftAngle = Opts.values['fDraftAngle']
    fDistTol = Opts.values['fDistTol']
    fAngleTol = Opts.values['fAngleTol']
    bSimplifyRes = Opts.values['bSimplifyRes']
    bArcOK = Opts.values['bArcOK']
    bDeg2 = Opts.values['bDeg2']
    bDeg3 = Opts.values['bDeg3']
    bDeg5 = Opts.values['bDeg5']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gCs_Res, sCrvDescrs = processDocObjects(rhBreps=objrefs)

    if not gCs_Res:
        if bEcho: print("No resultant curves.")
        sc.doc.Views.RedrawEnabled = True
        return

    if bEcho:
        print("Output curves:")
        for sCrvDescr in set(sCrvDescrs):
            print("{} {}".format(sCrvDescrs.count(sCrvDescr), sCrvDescr))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()