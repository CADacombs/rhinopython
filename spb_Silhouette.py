"""
This script will create draft curves silhouettes as obtained via _DraftAngleAnalysis,
except with curve simplification options and a more convenient execution.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
171019-21: Created.
200509: Updated for V7 using newly added Silhouette.ComputeDraftCurve method.
200511,13: Added routine that simplifies curves. Added some options.
200805: Modified pullDirection sign to match change in ComputeDraftCurve in Rhino 7.0.20217.3575, 8/4/2020. Bug fix in printed output.
210223: Added bArcOK. Bug fixes.
210315-17: Modified to check accuracy of output curves against input tolerances. Added bDeg5.
210515: Now by default, will try to rebuild both degrees 3 and 5 within tolerance.
        Now, supports sub-selected BrepFaces as input.
220821-22: Refactored main V5&6 from V7 routines into separate functions. Bug fixes.
240106-10: Removed code that only works in V5 & V6. Added more direction options.
        Modified simplification routine. Refactored.
        Now skips self-intersecting curves.
        Pulls curves not on face to the face. The error is not uncommon for ComputeDraftCurve.
        Added option and routine to split output curves at G2 discontinuities.
240725: Bug fix in finding G2 discontinuities.
        Bug fix in removing knots.
241229: Bug fixes.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

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

    key = 'fDraftAngle_Deg'; keys.append(key)
    values[key] = 0.0
    names[key] = 'DraftAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseDraftComputeFor0'; keys.append(key)
    values[key] = False
    names[key] = 'ComputeMethod'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'NoDraft', 'Draft')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistTol'; keys.append(key)
    values[key] = 0.5 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fAngle_Tol_Deg'; keys.append(key)
    values[key] = 0.5 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
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

    key = 'bSplitNonG2'; keys.append(key)
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

        if key == 'fDraftAngle_Deg':
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fDistTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fAngle_Tol_Deg':
            if cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

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
        addOption('fDraftAngle_Deg')
        if Opts.values['fDraftAngle_Deg'] == 0.0:
            addOption('bUseDraftComputeFor0')
        addOption('fDistTol')
        addOption('fAngle_Tol_Deg')
        addOption('bSimplifyRes')
        if Opts.values['bSimplifyRes']:
            addOption('bArcOK')
            addOption('bDeg2')
            addOption('bDeg3')
            addOption('bDeg5')
        addOption('bSplitNonG2')
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
            key = 'fDraftAngle_Deg'
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


def formatDistance(fDistance, iPrecision=15):
    if fDistance is None: return "(No deviation provided)"
    
    if fDistance < 0.01:
        return "{:.{}e}".format(fDistance, iPrecision)
    
    if fDistance < 0.1:
        return "{:.{}g}".format(fDistance, iPrecision+1)
    
    return "{:.{}g}".format(fDistance, iPrecision)


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


def doesCurveSelfIntersect(rgC, tolerance):
    if tolerance is None: tolerance = sc.doc.ModelAbsoluteTolerance
    rc = rg.Intersect.Intersection.CurveSelf(rgC, tolerance)
    return bool(rc)


def create_hi_def_pull_to_face(rgFace, rgC_on_face, tolerance_forRef=1e-7, bDebug=False):

    pulled = rgC_on_face.PullToBrepFace(face=rgFace, tolerance=tolerance_forRef)
    if len(pulled) == 1:
        return pulled[0]
    else:
        pullback = rgFace.Pullback(rgC_on_face, tolerance=tolerance_forRef)
        pushup = rgFace.Pushup(pullback, tolerance=tolerance_forRef)
        if bDebug:
            sEval='pullback.GetType().Name'; print(sEval+':',eval(sEval))
            sEval='pushup.GetType().Name'; print(sEval+':',eval(sEval))
#        sc.doc.Objects.AddCurve(rgC_on_face); sc.doc.Views.Redraw()
#        raise Exception("PullToBrepFace produced {} curves, but should have been 1.".format(
#            len(pulled)))
        return pushup


def create_normal_def_pull_to_face(rgFace, rgC_on_face, tolerance, bDebug=False):

    pulled = rgC_on_face.PullToBrepFace(face=rgFace, tolerance=tolerance)
    if len(pulled) == 1:
        return pulled[0]
    else:
        pullback = rgFace.Pullback(rgC_on_face, tolerance=tolerance)
        pushup = rgFace.Pushup(pullback, tolerance=tolerance)
        if bDebug:
            sEval='pullback.GetType().Name'; print(sEval+':',eval(sEval))
            sEval='pushup.GetType().Name'; print(sEval+':',eval(sEval))
#        sc.doc.Objects.AddCurve(rgC_on_face); sc.doc.Views.Redraw()
#        raise Exception("PullToBrepFace produced {} curves, but should have been 1.".format(
#            len(pulled)))
        return pushup


def isCurveOnFace(rgC_on_face, rgFace, tolerance, bDebug=False):

    hi_def_pull = create_hi_def_pull_to_face(rgFace, rgC_on_face)

    bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
        rgC_on_face, hi_def_pull, tolerance=0.1*tolerance)[:2]
    if not bSuccess:
        # Some results where this occurs:
        #   Closed curves that are collapse, appearing like an open segment.
        return
        sc.doc.Objects.AddCurve(rgC_on_face)
        sc.doc.Objects.AddCurve(hi_def_pull)
        sc.doc.Views.Redraw()
        raise Exception("GetDistancesBetweenCurves failed.")

    hi_def_pull.Dispose()

    # 1e-5 is an additional tolerance helper for when the tolerance is arguably negligibly missed.
    helper = 1e-5

    if bDebug:
        if fDistMax > tolerance:
            missed = fDistMax-tolerance
            if missed > helper:
                print("Tolerance was missed by {}.".format(missed))
                print("Recommended helper to deviation check: {:.0e}".format(missed))

    if fDistMax > (tolerance + helper):

        if bDebug:
            print("Tolerance is {}, but curve is {} off of face.".format(
                tolerance, fDistMax))

        #if isinstance(rgFace.UnderlyingSurface(), rg.RevSurface):
        #    pass
        #elif isinstance(rgFace.UnderlyingSurface(), rg.NurbsSurface) and rgFace.UnderlyingSurface().IsRational:
        #    pass
        #else:
        #    sc.doc.Objects.AddSurface(rgFace)
        #    sc.doc.Objects.AddCurve(rgC_on_face)
        #    sc.doc.Objects.AddCurve(hi_def_pull)
        #    sc.doc.Views.Redraw()
        #    1/0

        return False

    return True


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


def is_NurbsCrv_internally_G2_continuous(nc, fAngle_Tol_Deg, bDebug=False):
    if not isinstance(nc, rg.NurbsCurve):
        raise ValueError("{} passed to is_NurbsCrv_internally_G2_continuous.".format(nc))

    if nc.Degree >= 3 and isUniformNurbsCurve(nc):
        return True

    cosAngleTolerance = math.cos(math.radians(fAngle_Tol_Deg))

    bFoundG1_discont, t = nc.GetNextDiscontinuity(
        continuityType=rg.Continuity.G2_locus_continuous if nc.IsPeriodic else rg.Continuity.G2_continuous,
        t0=nc.Domain.T0,
        t1=nc.Domain.T1,
        cosAngleTolerance=cosAngleTolerance,
        curvatureTolerance=Rhino.RhinoMath.UnsetValue) # Value for curvatureTolerance doesn't seem to have an effect for finding G2 discontinuities.
    if bFoundG1_discont:
        if bDebug: print("First G2 discontinuity found at {}.".format(t))
        return False

    return True


def rebuild_to_Bezier(nc_In, iDegs, tolerance, bDebug=False):
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


def rebuild_to_MultiSpan(nc_In, iDegs, tolerance, bDebug=False):
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

    #iDegs = range(2, 7+1)

    for iDeg in iDegs:
        if iDeg < nc_In.Degree: continue

        for minimumMultiplicity in range(1, iDeg-2+1):
            nc_WIP = nc_In.DuplicateCurve()

            if iDeg > nc_In.Degree:
                nc_WIP.IncreaseDegree(iDeg)

            iCt_KnotsRemoved = nc_WIP.Knots.RemoveMultipleKnots(
                minimumMultiplicity=minimumMultiplicity,
                maximumMultiplicity=iDeg+1,
                tolerance=1.0)
            if not iCt_KnotsRemoved:
                if bDebug: print("RemoveMultipleKnots failed to remove any knots.")
                nc_WIP.Dispose()
                continue # to next minimumMultiplicity.

            bSuccess, fDev = rg.Curve.GetDistancesBetweenCurves(
                nc_WIP, nc_In, tolerance=0.1*tolerance)[:2]
            if not bSuccess:
                nc_WIP.Dispose()
                continue # to next minimumMultiplicity.

            if bDebug: sEval='fDev'; print(sEval+':',eval(sEval))
            if fDev > tolerance:
                nc_WIP.Dispose()
                continue # to next minimumMultiplicity.

            _tangentAngleDifference(nc_In, nc_WIP, bDebug)
            return nc_WIP, fDev



def simplifyCurve(rgCrv_In, bArcOK=False, iDegs=(2,3,5), fTol=None, fAngle_Tol_Deg=None, bDebug=False):
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
        fTol = max((fTol_Min, 0.1*sc.doc.ModelAbsoluteTolerance))

    if fAngle_Tol_Deg is None:
        fAngle_Tol_Deg = sc.doc.ModelAngleToleranceDegrees

    fTol_MaxRad = 1e3


    if isinstance(rgCrv_In, (rg.ArcCurve, rg.LineCurve)):
        print("Curve is already a {}. How did that happen?".format(
            rgCrv_In.GetType().Name))
        return rgCrv_In.DuplicateCurve()

    if not rgCrv_In.IsClosed:
        if rgCrv_In.IsLinear(fTol_Min):
            return rg.LineCurve(
                rgCrv_In.PointAtStart,
                rgCrv_In.PointAtEnd)

    if bArcOK:
        if rgCrv_In.IsClosed:
            b, circle = rgCrv_In.TryGetCircle(fTol_Min)
            if b:
                if circle.Radius <= fTol_MaxRad:
                    return rg.ArcCurve(circle)
        else:
            b, arc = rgCrv_In.TryGetArc(fTol_Min)
            if b:
                if arc.Radius <= fTol_MaxRad:
                    return rg.ArcCurve(arc)

        ## Just keep rational result.
        ## TODO: Remove this?
        #if nc_In.IsRational:
        #    return nc_In.Duplicate()

    if isinstance(rgCrv_In, rg.PolylineCurve):
        
        if sc.doc.Objects.AddCurve(rgCrv_In) == Guid.Empty:
            s = ""
        else:
            s = ", and was added to the document"
        raise Exception("Input is {}{}.".format(
            rgCrv_In.GetType().Name, s))

    if isinstance(rgCrv_In, rg.PolyCurve):
        nc_In = rgCrv_In.ToNurbsCurve()
        print("PolyCurve converted to a NurbsCurve for simplification processing.")
    else:
        nc_In = rgCrv_In.Duplicate()

    if isUniformNurbsCurve(nc_In):
        return nc_In.Duplicate()

    rc = rebuild_to_Bezier(
        nc_In,
        iDegs=iDegs,
        tolerance=fTol,
        bDebug=bDebug)
    if rc:
        return rc[0]

    rc = removeMultiKnots(nc_In, iDegs, fTol, bDebug)
    if rc: return rc[0]
    #if is_NurbsCrv_internally_G2_continuous(nc_In, fAngle_Tol_Deg, bDebug):
    #    rc = removeMultiKnots(nc_In, iDegs, fTol, bDebug)
    #    if rc: return rc[0]

    rc = rebuild_to_MultiSpan(
        nc_In,
        iDegs=iDegs,
        tolerance=fTol,
        bDebug=bDebug)

    if rc:
        return rc[0]

    nc_In.Dispose()


def split_NurbsCrv_at_G2_discontinuities(nc_In, angleToleranceRadians, bDebug=False):
    if not isinstance(nc_In, rg.NurbsCurve):
        raise ValueError("{} passed to split_NurbsCrv_at_G2_discontinuities.".format(nc_In))

    cosAngleTolerance = math.cos(angleToleranceRadians)

    t0 = nc_In.Domain.T0 # Start parameter of check. Is itself ignored.

    ts_Found = []

    while True:
        sc.escape_test()

        bFoundG1_discont, t = nc_In.GetNextDiscontinuity(
            continuityType=rg.Continuity.G2_locus_continuous if nc_In.IsPeriodic else rg.Continuity.G2_continuous,
            t0=t0,
            t1=nc_In.Domain.T1,
            cosAngleTolerance=cosAngleTolerance,
            curvatureTolerance=Rhino.RhinoMath.UnsetValue) # Value for curvatureTolerance doesn't seem to have an effect for finding G2 discontinuities.

        if not bFoundG1_discont:
            break # out of while loop.

        ts_Found.append(t)

        t0 = t

    if len(ts_Found) == 0:
        return

    split = nc_In.Split(ts_Found)

    if len(split) < 2:
        s = "Split produced {} curves. Should have been {}. Curves were added.".format(
            len(split), len(ts_Found)+1)
        print(s)
        sc.doc.Objects.AddCurve(nc_In)
        [sc.doc.Objects.AddCurve(c) for c in split]

        raise Exception(s)

    joined = rg.Curve.JoinCurves(split)

    if len(joined) != 1:
        raise Exception("JoinCurves produced {} curves. Should have been 1.".format(
            len(joined)))

    for c in split: c.Dispose()

    return joined


def split_NurbsCrvs_at_G2_discontinuities(rgCrvs_In, angleToleranceRadians, bDebug=False):
    """
    Return: list(Curves) with splits and duplicates,
            It is up to the calling code to Dispose the curves in the original list.
    """
    rgCs_Out = []
    for c in rgCrvs_In:
        if not isinstance(c, rg.NurbsCurve):
            rgCs_Out.append(c.DuplicateCurve())
            continue

        rc = split_NurbsCrv_at_G2_discontinuities(c, angleToleranceRadians, bDebug=bDebug)

        if rc is None:
            rgCs_Out.append(c.DuplicateCurve())
            continue

        rgCs_Out.extend(rc)

    return rgCs_Out


def createSilhouetteCurves(rgGeomForSilh, **kwargs):
    """
    Main function that calls the "createSilhouetteCurves..." functions,
    and post-processes the results.

    Parameters:
        rgGeomForSilh
        kwargs:
            iPullDir
            fDraftAngle_Deg
            fDistTol
            fAngle_Tol_Deg
            bSimplifyRes
            bArcOK
            bDeg3
            bDeg5
            bSplitNonG2
            bEcho
            bDebug
    """

    if Rhino.RhinoApp.ExeVersion < 7:
        print("This function requires Rhino 7 or higher.")
        return

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPullDir = getOpt('iPullDir')
    fDraftAngle_Deg = getOpt('fDraftAngle_Deg')
    bUseDraftComputeFor0 = getOpt('bUseDraftComputeFor0')
    fDistTol = getOpt('fDistTol')
    fAngle_Tol_Deg = getOpt('fAngle_Tol_Deg')
    bSimplifyRes = getOpt('bSimplifyRes')
    bArcOK = getOpt('bArcOK')
    bDeg2 = getOpt('bDeg2')
    bDeg3 = getOpt('bDeg3')
    bDeg5 = getOpt('bDeg5')
    bSplitNonG2 = getOpt('bSplitNonG2')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    # Prepare variables to be passed to Compute or ComputeDraftCurves.
    draftAngleRadians = math.radians(fDraftAngle_Deg)
    pullDirection = _getDirectionVector(Opts.values['iPullDir'])
    tol_ForCompute_or_CDC = 0.5*fDistTol if bSimplifyRes else fDistTol
    angleToleranceRadians = math.radians(fAngle_Tol_Deg)


    iDegs = []
    if bDeg2:
        iDegs.append(2)
    if bDeg3:
        iDegs.append(3)
    if bDeg5:
        iDegs.append(5)
    if not iDegs and bEcho:
        print("No degrees for rebuild.")




    if isinstance(rgGeomForSilh, rg.Brep):
        rgFs = rgGeomForSilh.Faces
    elif isinstance(rgGeomForSilh, rg.BrepFace):
        rgFs = [rgGeomForSilh]
    else:
        raise ValueError("{} is not supported for this script.".format(rgGeomForSilh.GetType().Name))

    if (fDraftAngle_Deg != 0.0) or bUseDraftComputeFor0:
        bUseDraftCompute = True
        silhTypes_Allowed = rg.SilhouetteType.DraftCurve, rg.SilhouetteType.Tangent
    else:
        # Use Silhouette.Compute.
        bUseDraftCompute = False
        silhTypes_Allowed = rg.SilhouetteType.Tangent,


    if bDebug:
        sEval='fDraftAngle_Deg'; print(sEval+':',eval(sEval))
        sEval='draftAngleRadians'; print(sEval+':',eval(sEval))
        sEval="Opts.values['iPullDir']"; print(sEval+':',eval(sEval))
        sEval='pullDirection'; print(sEval+':',eval(sEval))
        sEval='fDistTol'; print(sEval+':',eval(sEval))
        sEval='fAngle_Tol_Deg'; print(sEval+':',eval(sEval))
        sEval='angleToleranceRadians'; print(sEval+':',eval(sEval))
        sEval='iDegs'; print(sEval+':',eval(sEval))
        sEval='bUseDraftCompute'; print(sEval+':',eval(sEval))
        sEval='silhTypes_Allowed'; print(sEval+':',eval(sEval))



    rgCs_WIP = []
    #fLengths_Totals = [0.0]

    Rhino.RhinoApp.SetCommandPrompt("Collecting filtered silhouette curves ...")

    for rgF in rgFs:
        if bUseDraftCompute:
            silhouettes = rg.Silhouette.ComputeDraftCurve(
                geometry=rgF,
                draftAngle=draftAngleRadians,
                pullDirection=pullDirection,
                tolerance=tol_ForCompute_or_CDC,
                angleToleranceRadians=angleToleranceRadians)
        else:
            silhouettes = rg.Silhouette.Compute(
                geometry=rgF,
                silhouetteType=rg.SilhouetteType.Tangent,
                parallelCameraDirection=pullDirection,
                tolerance=tol_ForCompute_or_CDC,
                angleToleranceRadians=angleToleranceRadians)

        if not silhouettes:
            continue

        for silh in silhouettes:
            if silh.Curve is None:
                continue

            if silh.SilhouetteType not in silhTypes_Allowed:
                if bDebug: print("Skipped {}.".format(silh.SilhouetteType))
                continue

            fLength = silh.Curve.GetLength()
            if fLength < tol_ForCompute_or_CDC:
                if bDebug: print("Curve with length of {} skipped.".format(formatDistance(fLength)))
                continue
            #fLengths_Totals[-1] += fLength
            #if bDebug: sEval="fLengths_Totals[-1]"; print(sEval+':',eval(sEval))

            if doesCurveSelfIntersect(silh.Curve, tol_ForCompute_or_CDC):
                if bDebug: print("Self-intersecting curve skipped.")
                continue

            if isCurveOnFace(silh.Curve, rgF, tol_ForCompute_or_CDC, bDebug=bDebug):
                rgCs_WIP.append(silh.Curve)
            else:
                new_pull = create_normal_def_pull_to_face(rgF, silh.Curve, tol_ForCompute_or_CDC, bDebug=bDebug)
                silh.Curve.Dispose()
                rgCs_WIP.append(new_pull)

    if not rgCs_WIP:
        return rgCs_WIP


    if not bSimplifyRes and not bSplitNonG2:
        return rgCs_WIP


    if bSimplifyRes:
        Rhino.RhinoApp.SetCommandPrompt("Simplifying ...")
        rgCs_Post_simplified = []
        for c in rgCs_WIP:
            rgC_Simplified = simplifyCurve(
                c,
                bArcOK=bArcOK,
                iDegs=iDegs,
                fTol=tol_ForCompute_or_CDC,
                fAngle_Tol_Deg=fAngle_Tol_Deg,
                bDebug=bDebug)

            if rgC_Simplified is None:
                rgCs_Post_simplified.append(c)
            else:
                rgCs_Post_simplified.append(rgC_Simplified)
                c.Dispose()
        rgCs_WIP = rgCs_Post_simplified


    if bSplitNonG2:
        Rhino.RhinoApp.SetCommandPrompt("Splitting ...")
        rgCs_Post_split = split_NurbsCrvs_at_G2_discontinuities(
            rgCs_WIP,
            angleToleranceRadians,
            bDebug=bDebug,
            )
        for c in rgCs_WIP: c.Dispose()
        rgCs_WIP = rgCs_Post_split


    return rgCs_WIP


def processDocObjects(rhBreps, **kwargs):
    """
    Parameters:
        rhBreps: list(objref or GUID)

    Returns: list(GUID) of CurveObjects
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPullDir = getOpt('iPullDir')
    fDraftAngle_Deg = getOpt('fDraftAngle_Deg')
    bUseDraftComputeFor0 = getOpt('bUseDraftComputeFor0')
    fDistTol = getOpt('fDistTol')
    fAngle_Tol_Deg = getOpt('fAngle_Tol_Deg')
    bSimplifyRes = getOpt('bSimplifyRes')
    bArcOK = getOpt('bArcOK')
    bDeg2 = getOpt('bDeg2')
    bDeg3 = getOpt('bDeg3')
    bDeg5 = getOpt('bDeg5')
    bSplitNonG2 = getOpt('bSplitNonG2')
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

        rgCs_Res = createSilhouetteCurves(
            rdObj_In,
            iPullDir=iPullDir,
            fDraftAngle_Deg=fDraftAngle_Deg,
            bUseDraftComputeFor0=bUseDraftComputeFor0,
            fDistTol=fDistTol,
            fAngle_Tol_Deg=fAngle_Tol_Deg,
            bSimplifyRes=bSimplifyRes,
            bArcOK=bArcOK,
            bDeg2=bDeg2,
            bDeg3=bDeg3,
            bDeg5=bDeg5,
            bSplitNonG2=bSplitNonG2,
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
    fDraftAngle_Deg = Opts.values['fDraftAngle_Deg']
    bUseDraftComputeFor0 = Opts.values['bUseDraftComputeFor0']
    fDistTol = Opts.values['fDistTol']
    fAngle_Tol_Deg = Opts.values['fAngle_Tol_Deg']
    bSimplifyRes = Opts.values['bSimplifyRes']
    bArcOK = Opts.values['bArcOK']
    bDeg2 = Opts.values['bDeg2']
    bDeg3 = Opts.values['bDeg3']
    bDeg5 = Opts.values['bDeg5']
    bSplitNonG2 = Opts.values['bSplitNonG2']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    rc = processDocObjects(rhBreps=objrefs)

    if rc is None: return

    gCs_Res, sCrvDescrs = rc

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