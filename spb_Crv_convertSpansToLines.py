"""
This script will convert linear contiguous spans to LineCurves.
It can be used as an alternative to _SimplifyCrv and _Convert _Output=Lines.

Merge=No means not to combine spans.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180905: Created.
...
230824: Reenabled allowance of modifying options when curves are preselected.
        Fixed bug that continued to process whole curves that were found to be arc, etc., -shaped.
241220,21: Bug fix: Whole-curve conversion to a single line is no longer attempted on closed curves.
        Added tolerance per length deviation option.
241222: Now can optionally create a single line through end-start of closed curve when they are collinear.
        Now merge is optional.
241223: Refactored.

TODO:
    Add an angle tolerance?



From https://developer.rhino3d.com/api/rhinocommon :
    Curve.DuplicateSegments  Curve[] of segments  Duplicates curve segments. Explodes polylines, polycurves and G1 discontinuous NURBS curves. Single segment curves, such as lines, arcs, unkinked NURBS curves, are duplicated.
    PolyCurve.Explode  Curve[] of polycurve segments  Explodes this PolyCurve into a list of Curve segments. This willnot explodenested polycurves. Call RemoveNesting first if you need all individual segments.
    PolyCurve.RemoveNesting  bool  Explodes nested polycurve segments and reconstructs this curve from the shattered remains. The result will have not have any PolyCurves as segments but it will have identical locus and parameterization.
    PolyCurve.CleanUp  Curve or None  Removes any nesting of polycurves. If this polycurve has just a single segment, the segment is returned. If, after nest removal, there are adjacent segments which are polylines, they are combined into a single polyline. The new curve may have a different domain from this polycurve. If the start and end segments of a closed input are polylines, the result may have a different seam location since the start and end segments will be combined.
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


    key = 'fDevTol'; keys.append(key)
    values[key] = max((0.1 * sc.doc.ModelAbsoluteTolerance, 1e-6))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAlsoRatioTol'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTolPerLengthUnit'; keys.append(key)
    values[key] = 0.001
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fMinLineLength'; keys.append(key)
    values[key] = 100.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bProcessArcs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMerge'; keys.append(key)
    values[key] = True
    names[key] = 'MergeContiguousSegs'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMergeThruSeam'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

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

        if key in ('fDevTol', 'fMinLineLength', 'fTolPerLengthUnit'):
            if cls.riOpts[key].CurrentValue <= 1e-6:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def addGeoms(geoms, bRedraw=True):
    """ For debugging. """

    gOuts = []

    if not hasattr(geoms, '__iter__'):
        geoms = [geoms]

    for geom in geoms:
        try:
            rc = geom.IsValidWithLog()
            if not rc[0]:
                print("{} skipped.\n  {}".format(geom.GetType().Name, rc[1]))
                continue
        except:
            pass
        if isinstance(geom, tuple):
            gOut = sc.doc.Objects.AddSurface(geom[1])
        elif isinstance(geom, rg.Line):
            gOut = sc.doc.Objects.AddLine(geom)
        elif isinstance(geom, rg.Curve):
            gOut = sc.doc.Objects.AddCurve(geom)
        elif isinstance(geom, rg.Surface):
            gOut = sc.doc.Objects.AddSurface(geom)
        elif isinstance(geom, rg.Point3d):
            gOut = sc.doc.Objects.AddPoint(geom)
        elif isinstance(geom, rg.Plane):
            intrvl = rg.Interval(-sc.doc.ModelAbsoluteTolerance*1000.0, sc.doc.ModelAbsoluteTolerance*1000.0)
            psrf = rg.PlaneSurface(geom, intrvl, intrvl)
            gOut = sc.doc.Objects.AddSurface(psrf)
        else:
            raise ValueError("Method to add {} missing from addGeoms.".format(geom.GetType().Name))
        if gOut == gOut.Empty:
            print("{} could not be added to document.".format(rc.GetType().Name))
        else:
            gOuts.append(gOut)
    if bRedraw: sc.doc.Views.Redraw()
    return gOuts


def getFormattedDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def getInput():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")
    
    go.GeometryFilter = rd.ObjectType.Curve
    
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

        addOption('fDevTol')
        addOption('fMinLineLength')
        addOption('bAlsoRatioTol')
        if Opts.values['bAlsoRatioTol']:
            addOption('fTolPerLengthUnit')
        addOption('bProcessArcs')
        addOption('bMerge')
        if Opts.values['bMerge']:
            addOption('bMergeThruSeam')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

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
                    rdObj = objref.Object()
                    if rdObj.ObjectType == rd.ObjectType.InstanceReference:
                        print("Objects in block instances are not supported.")
                        continue
                    rgEdges.append(objref.Geometry())

            objrefs = go.Objects()

            go.Dispose()

            return objrefs

        if res == ri.GetResult.Number:
            key = 'fDevTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getPassingLineCurve(rgCrv_In, bProcessArcs=False, fDevTol=None, fTolPerLengthUnit=None, fMinLineLength=123, bDebug=False):
    """
    Returns on success: rg.LineCurve, float(Deviation), None
    Returns on deviation fail: None, float(Deviation), None
    Returns on other fails: None, None, str(Description of fail)
    """

    if rgCrv_In is None:
        return None, None, "Geometry not found!  It will be skipped."


    if not bProcessArcs:
        if isinstance(rgCrv_In, rg.ArcCurve):
            return None, None, "Skipped ArcCurve."
            
            
        def tryGetPracticalArc(rgCrv_In):
            tol_IsArc = Rhino.RhinoMath.ZeroTolerance #max((1e-6, 0.001*sc.doc.ModelAbsoluteTolerance))
            bSuccess, arc = rgCrv_In.TryGetArc(tol_IsArc)
            if not bSuccess: return
            maxRadius = sc.doc.ModelAbsoluteTolerance * 1e6
            if arc.Radius > maxRadius:
                if bDebug: print("Skipped arc shape with radius {}".format(arc.Radius))
                return
            minAngleDegrees = 0.1*sc.doc.ModelAngleToleranceDegrees
            if arc.AngleDegrees < minAngleDegrees:
                if bDebug: print("Skipped arc shape with AngleDegrees {}".format(arc.AngleDegrees))
                return
            if bDebug:
                sEval = "arc.AngleDegrees"; print(sEval,'=',eval(sEval))
                sEval = "arc.Circumference"; print(sEval,'=',eval(sEval))
                sEval = "arc.Radius"; print(sEval,'=',eval(sEval))
            return arc
        
        arc = tryGetPracticalArc(rgCrv_In)
        if arc:
            #sc.doc.Objects.AddArc(arc)
            return None, None, "Skipped arc-shaped non-ArcCurve."
        #if rgCrv_In.IsArc(tol_IsArc):
        #    return None, None, "Skipped arc-shaped non-ArcCurve."

    rgLc_Out = rg.LineCurve(rgCrv_In.PointAtStart, rgCrv_In.PointAtEnd)

    if fMinLineLength is not None:
        length_crv1 = rgLc_Out.GetLength()
        if length_crv1 < fMinLineLength:
            return None, None, "Curve is too short."

    if fDevTol is None:
        fDevTol = 0.1 * sc.doc.ModelAbsoluteTolerance

    if fTolPerLengthUnit:
        fLength = rgCrv_In.GetLength()
        fDevTol_Abs = min((fTolPerLengthUnit * fLength), fDevTol)
    else:
        fDevTol_Abs = fDevTol


    tol_GDBC = 0.1*min(sc.doc.ModelAbsoluteTolerance, fDevTol_Abs)

    rc = rg.Curve.GetDistancesBetweenCurves(
        rgCrv_In,
        rgLc_Out,
        tolerance=tol_GDBC)

    if not rc[0]:
        return None, None, "Deviation could not be determined by Curve.GetDistancesBetweenCurves."

    fDev = rc[1]

    if fDev <= fDevTol_Abs:
        return rgLc_Out, fDev, None

    return None, fDev, None


def _getUniqueKnots(nc, bDebug=False):
    # All unique knots along curve domain except allow
    # a single repeat for closed curves, but
    # no additional duplicates for periodic curves.
    # These start and stop work regardless of IsClosed and IsPeriodic.
    start = nc.Degree-1
    stop = nc.Knots.Count - nc.Degree + 1

    ts_knots = []
    for i in range(start, stop):
        t = nc.Knots[i]
        if t not in ts_knots:
            ts_knots.append(t)
            if bDebug: rs.AddTextDot('{0:.4f}'.format(t), nc.PointAt(t))

    return ts_knots


def _getArcEndsParameterPairs(pc):
    ts = []
    for i in range(pc.SegmentCount):
        seg = pc.SegmentCurve(i)
        if isinstance(seg, rg.ArcCurve):
            ts.append((
                rg.PolyCurve.PolyCurveParameter(pc, i, seg.Domain.T0),
                rg.PolyCurve.PolyCurveParameter(pc, i, seg.Domain.T1),
                ))
    return ts


def _removeInteriorArcParameters(rgC_In, ts_knots):
    ts_arcPairs = _getArcEndsParameterPairs(rgC_In)

    for t0, t1 in ts_arcPairs:
        for tK in reversed(ts_knots):
            if (t0 + Rhino.RhinoMath.ZeroTolerance) < tK < (t1 - Rhino.RhinoMath.ZeroTolerance):
                ts_knots.remove(tK)


def _doesPolyCrvContainArcCrvs(pc):
    for i in range(pc.SegmentCount):
        seg = pc.SegmentCurve(i)
        if isinstance(seg, rg.ArcCurve):
            return True
        if seg.IsArc():
            return True
    return False


def _doesPolyCrvContainNurbsCrvs(pc):
    for i in range(pc.SegmentCount):
        seg = pc.SegmentCurve(i)
        if isinstance(seg, rg.NurbsCurve):
            return True
    return False


def _splitPolyCrvByParams(pc_In, ts):
    pc_Segs_All = list(pc_In.Split(ts))
    for i in range(len(pc_Segs_All)):
        if pc_Segs_All[i].SpanCount != 1:
            raise
        # Overwrite PolyCurve-type segment with the
        # single non-PolyCurve-type segment in segment.
        pc_Segs_All[i] = pc_Segs_All[i].SegmentCurve(0)
    return pc_Segs_All


def convertSpans_of_PolylineCrv(plc_In, fDevTol=None, fTolPerLengthUnit=None, fMinLineLength=None, bProcessArcs=False, bMergeThruSeam=True, bDebug=False):
    """
    """

    if not isinstance(plc_In, rg.PolylineCurve):
        raise Exception("{} passed to convertSpans_of_PolylineCrv.".format(plc_In.GetType().Name))

    segs_All = plc_In.DuplicateSegments()
    pLine = plc_In.ToPolyline()
    pts_All = pLine.ToArray()

    idx_Pts_toKeep = [0]

    idx_t_linear_start = 0
    devs = []

    for i in range(1, pLine.Count-1):
        if sc.escape_test(throw_exception=False):
            raise Exception("Break at index {}.".format(i))

        if bDebug: sEval = "i"; print(sEval,'=',eval(sEval))

        pts_forSection = [pts_All[j] for j in range(idx_t_linear_start, i+1+1)]

        linearSection = rg.PolylineCurve(pts_forSection)

        rc = getPassingLineCurve(
            linearSection,
            bProcessArcs=False,
            fDevTol=fDevTol,
            fMinLineLength=None)

        linearSection.Dispose()

        if rc[0] is not None:
            bContiguousLinearFound = True
            devs.append(rc[1])

            continue # to the next index.

        # Failed, so record up to the current index.

        idx_Pts_toKeep.append(i)
        idx_t_linear_start = i

    # Record the last section.
    idx_Pts_toKeep.append(pLine.Count-1)


    if plc_In.IsClosed and bMergeThruSeam and len(idx_Pts_toKeep) > 3:
        pts_forSection = []
        pts_forSection = [pts_All[j] for j in
                         (idx_Pts_toKeep[-2], idx_Pts_toKeep[0], idx_Pts_toKeep[1])]
        linearSection = rg.PolylineCurve(pts_forSection)

        rc = getPassingLineCurve(
            linearSection,
            bProcessArcs=False,
            fDevTol=fDevTol,
            fMinLineLength=None)

        if rc[0] is not None:
            bContiguousLinearFound = True
            devs.append(rc[1])
        idx_Pts_toKeep[0] = idx_Pts_toKeep[-2]
        del idx_Pts_toKeep[-1]

    if not devs:
        return None, None, "No segments of PolylineCurve were combined."

    pLc_Out = rg.PolylineCurve([pts_All[i] for i in idx_Pts_toKeep])

    return (
        pLc_Out,
        max(devs),
        None,
        )


def convertSpans_withMerging(rgC_In, fDevTol=None, fTolPerLengthUnit=None, fMinLineLength=None, bProcessArcs=False, bMergeThruSeam=True, bDebug=False):
    """
    Processes NURBS-form of input curve even if it is a PolyCurve.
    """

    if isinstance(rgC_In, rg.PolylineCurve):
        return convertSpans_of_PolylineCrv(
            cleaned,
            fDevTol=fDevTol,
            fTolPerLengthUnit=fTolPerLengthUnit,
            fMinLineLength=fMinLineLength,
            bProcessArcs=bProcessArcs,
            bMergeThruSeam=bMergeThruSeam,
            bDebug=bDebug,
            )

    if isinstance(rgC_In, rg.NurbsCurve):
        nc_FromIn = rgC_In.Duplicate()
        ts_knots = _getUniqueKnots(nc_FromIn, bDebug=bDebug)
        segs_All = nc_FromIn.Split(ts_knots)
    elif isinstance(rgC_In, rg.PolyCurve):
        cleaned = rgC_In.CleanUp()

        if cleaned is None:
            pc_fromIn = rgC_In.Duplicate()
        else:
            if isinstance(cleaned, rg.LineCurve):
                return cleaned, 0.0, "PolyCurve contained just one LineCurve."
            if isinstance(cleaned, rg.PolylineCurve):
                return convertSpans_of_PolylineCrv(
                    cleaned,
                    fDevTol=fDevTol,
                    fTolPerLengthUnit=fTolPerLengthUnit,
                    fMinLineLength=fMinLineLength,
                    bProcessArcs=bProcessArcs,
                    bMergeThruSeam=bMergeThruSeam,
                    bDebug=bDebug,
                    )

            pc_fromIn = cleaned

        nc_FromIn = pc_fromIn.ToNurbsCurve()

        if not rgC_In.Domain.EpsilonEquals(nc_FromIn.Domain, 1e-9):
            raise Exception(
                "NurbsCurve's domain doesn't EpsilonEqual the input curve's.")

        ts_knots = _getUniqueKnots(nc_FromIn, bDebug=bDebug)

        if bDebug: print(ts_knots)
        _removeInteriorArcParameters(rgC_In, ts_knots)
        if bDebug: print(ts_knots)

        segs_All = _splitPolyCrvByParams(rgC_In, ts_knots)

    if bDebug:
        sEval = 'len(ts_knots)'; print(sEval + ':' , eval(sEval))
        sEval = 'ts_knots[0]'; print(sEval + ':' , eval(sEval))
        #rs.AddTextDot("t:{}".format(0), nc0.PointAt(ts_knots[0]))
        sEval = 'ts_knots[-1]'; print(sEval + ':' , eval(sEval))
        #rs.AddTextDot("t:{}".format(-1), nc0.PointAt(ts_knots[-1]))
        # Add segment index.
        for i in range(len(ts_knots)-1):
            t = (ts_knots[i] + ts_knots[i+1]) / 2.0
            rs.AddTextDot(str(i), nc_FromIn.PointAt(t))


    bools_IsSegLinear = []
    devs = []
    tols_Needed = []

    for i, seg in enumerate(segs_All):
        if i == 2:
            pass
        rc = getPassingLineCurve(
            seg,
            bProcessArcs=bProcessArcs,
            fDevTol=fDevTol,
            fTolPerLengthUnit=None, # Will later be checked for each tuples_splits.
            fMinLineLength=None, # Will later be checked for each tuples_splits.
            )
        if rc and rc[0]:
            bools_IsSegLinear.append(True)
            devs.append(rc[1])
        else:
            bools_IsSegLinear.append(False)
            tols_Needed.append(rc[1])

    if not any(bools_IsSegLinear):
        return None, None, "No linear segments in input."

    i = 0 # Segment and parameter index.
    idx_t_linear_start = None
    idx_t_nonlinear_start = None
    tuples_splits = [] # tuple(tuple(int(beginning nc segment index), int(ending nc segment index)), bool(True if indices are linear section))
    #bLinearFound = False
    bContiguousLinearFound = False

    while i < len(segs_All):
        if sc.escape_test(throw_exception=False):
            raise Exception("Break at index {}.".format(i))

        if bDebug: sEval = "i"; print(sEval,'=',eval(sEval))

        if not bools_IsSegLinear[i]:
            # Start of non-linear, contiguous group.

            if idx_t_nonlinear_start is None:
                idx_t_nonlinear_start = i

            if idx_t_linear_start is not None:
                # Record linear section.
                if bContiguousLinearFound:
                    tuples_splits.append(((idx_t_linear_start, i), True))
                else:
                    tuples_splits.append(((idx_t_linear_start, i), False))

                bContiguousLinearFound = False
                idx_t_linear_start = None

            i += 1
            continue

        # Start of linear, contiguous group.

        #bLinearFound = True

        if idx_t_linear_start is None:
            idx_t_linear_start = i

            if idx_t_nonlinear_start is not None:
                tuples_splits.append(((idx_t_nonlinear_start, i), False))
                idx_t_nonlinear_start = None

            #devs.append(devs_1stPass[i])

            i += 1
            continue

        # Create a temporary curve and test for its total linearity.
        seg = nc_FromIn.Trim(
            rg.Interval(
                ts_knots[idx_t_linear_start],
                ts_knots[i+1]))

        rc = getPassingLineCurve(
            seg,
            bProcessArcs=bProcessArcs,
            fDevTol=fDevTol,
            fTolPerLengthUnit=None,
            fMinLineLength=None)
        #addGeoms(seg)
        #print(rc[2])
        if rc[0]:

            # Now test with length checks.
            rc = getPassingLineCurve(
                seg,
                bProcessArcs=bProcessArcs,
                fDevTol=fDevTol,
                fTolPerLengthUnit=fTolPerLengthUnit,
                fMinLineLength=fMinLineLength)

            if rc[0]:
                bContiguousLinearFound = True
                devs.append(rc[1])

            seg.Dispose()

            # Try the next index.

            i += 1
            continue

        seg.Dispose()

        # MultiLinear failed.  Record up to the current index.

        if bContiguousLinearFound:
            tuples_splits.append(((idx_t_linear_start, i), True))
        else:
            tuples_splits.append(((idx_t_linear_start, i), False))

        bContiguousLinearFound = False
        idx_t_linear_start = None

    # Record the last section.
    if idx_t_linear_start is None:
        tuples_splits.append(((idx_t_nonlinear_start, i), False))
    else:
        if bContiguousLinearFound:
            tuples_splits.append(((idx_t_linear_start, i), True))
        else:
            tuples_splits.append(((idx_t_linear_start, i), False))


    def tryMergeThruStart(rgC_In, ts_knots, tuples_splits):
        #sEval = "rgC_In.IsClosed"; print(sEval,'=',eval(sEval))
        if not rgC_In.IsClosed:
            return

        #sEval = "tuples_splits[0]"; print(sEval,'=',eval(sEval))
        #sEval = "tuples_splits[-1]"; print(sEval,'=',eval(sEval))

        if not (tuples_splits[0] and tuples_splits[-1]):
            return

        seg_Start = nc_FromIn.Trim(
            rg.Interval(
                ts_knots[tuples_splits[0][0][0]],
                ts_knots[tuples_splits[0][0][1]]))

        seg_End = nc_FromIn.Trim(
            rg.Interval(
                ts_knots[tuples_splits[-1][0][0]],
                ts_knots[tuples_splits[-1][0][1]]))

        rv = rg.Curve.JoinCurves([seg_Start, seg_End])
        if rv.Count != 1:
            return

        joined = rv[0]

        rc = getPassingLineCurve(
            joined,
            fDevTol=fDevTol,
            fTolPerLengthUnit=fTolPerLengthUnit,
            fMinLineLength=fMinLineLength)
        #addGeoms(joined)
        #print(rc[2])
        joined.Dispose()
        if bool(rc[0]):
            return rc[1]

    if bMergeThruSeam:
        rv = tryMergeThruStart(rgC_In, ts_knots, tuples_splits)
        sEval = "rv"; print(sEval,'=',eval(sEval))
        if rv is not None:
            #sEval = "tuples_splits"; print(sEval,'=',eval(sEval))
            tuples_splits[0] = (tuples_splits[-1][0][0], tuples_splits[0][0][1]), True
            #sEval = "tuples_splits[0]"; print(sEval,'=',eval(sEval))
            del tuples_splits[-1]
            #sEval = "tuples_splits"; print(sEval,'=',eval(sEval))
            bContiguousLinearFound = True
            devs.append(rv)


    # Combine any contiguous non-linear tuples_splits.
    tuples_splits_before_combines = tuples_splits[:]
    i = len(tuples_splits) - 2
    while i >= 0:
        sc.escape_test()
        if not tuples_splits[i][1] and not tuples_splits[i+1][1]:
            tuples_splits[i] = (tuples_splits[i][0][0], tuples_splits[i+1][0][1]), False
            del tuples_splits[i+1]
        i -= 1


    #if not bLinearFound:
    #    return None, None, "No linear segments found."


    # Proceed replacing segments with LineCurves.
    bLineSubstitutions = False
    pc_Out = rg.PolyCurve()

    for (idxS, idxE), bLinear in tuples_splits:
        if bLinear:
            line = rg.Line(nc_FromIn.PointAt(ts_knots[idxS]),
                nc_FromIn.PointAt(ts_knots[idxE]))
            if line.Length >= fMinLineLength:
                pc_Out.Append(line)
                bLineSubstitutions = True
                #addGeoms(line)
                continue # to next tuples_splits.

        # Using rgC_In instead of nc_FromIn so as to not increase the degree of any NurbsCurves.
        seg_ForOut = rgC_In.Trim(
            rg.Interval(ts_knots[idxS], ts_knots[idxE]))
        pc_Out.Append(seg_ForOut)
        #addGeoms(seg_ForOut)

    if not bLineSubstitutions:
        pc_Out.Dispose()
        return None, None, "No lines replaced any part of curve."


    # Success:

    bSuccess, plc = pc_Out.TryGetPolyline()
    if bSuccess:
        return (
            plc,
            max(devs),
            None,
            )

    return (
        pc_Out,
        max(devs),
        None,
        )


def convertSpans_NoMerging(rgC_In, fDevTol=None, fTolPerLengthUnit=None, fMinLineLength=None, bProcessArcs=False, bDebug=False):
    """
    Processes NURBS-form of input curve even if it is a PolyCurve.
    """

    if isinstance(rgC_In, rg.PolylineCurve):
        return None, None, "PolylineCurve cannot be processed without allowing segment merging."


    # Process PolyCurve first in case other curve types are just hidden within nesting and a single segment.

    if isinstance(rgC_In, rg.NurbsCurve):
        nc_FromIn = rgC_In.Duplicate()
        ts_knots = _getUniqueKnots(nc_FromIn, bDebug=bDebug)
        segs_All = nc_FromIn.Split(ts_knots)
    elif isinstance(rgC_In, rg.PolyCurve):
        cleaned = rgC_In.CleanUp()

        if cleaned is None:
            pc_fromIn = rgC_In.Duplicate()
        else:
            if isinstance(cleaned, rg.LineCurve):
                return cleaned, 0.0, "PolyCurve contained just one LineCurve."
            if isinstance(cleaned, rg.PolylineCurve):
                return cleaned, 0.0, "PolyCurve simplified to PolylineCurve but cannot be further processed without allowing segment merging."

            pc_fromIn = cleaned

        if bProcessArcs:
            if not(_doesPolyCrvContainArcCrvs(pc_fromIn) or _doesPolyCrvContainNurbsCrvs(pc_fromIn)):
                pc_fromIn.Dispose()
                return None, None, "No ArcCurves or NurbsCurve segments found for further processing."
        else:
            if not _doesPolyCrvContainNurbsCrvs(pc_fromIn):
                pc_fromIn.Dispose()
                return None, None, "No NurbsCurve segments found for further processing."

        nc_FromIn = pc_fromIn.ToNurbsCurve()

        if not rgC_In.Domain.EpsilonEquals(nc_FromIn.Domain, 1e-9):
            raise Exception(
                "NurbsCurve's domain doesn't EpsilonEqual the input curve's.")

        ts_knots = _getUniqueKnots(nc_FromIn, bDebug=bDebug)

        if bDebug: print(ts_knots)
        _removeInteriorArcParameters(rgC_In, ts_knots)
        if bDebug: print(ts_knots)

        segs_All = _splitPolyCrvByParams(rgC_In, ts_knots)

    if bDebug:
        sEval = 'len(ts_knots)'; print(sEval + ':' , eval(sEval))
        sEval = 'ts_knots[0]'; print(sEval + ':' , eval(sEval))
        #rs.AddTextDot("t:{}".format(0), nc0.PointAt(ts_knots[0]))
        sEval = 'ts_knots[-1]'; print(sEval + ':' , eval(sEval))
        #rs.AddTextDot("t:{}".format(-1), nc0.PointAt(ts_knots[-1]))
        # Add segment index.
        for i in range(len(ts_knots)-1):
            t = (ts_knots[i] + ts_knots[i+1]) / 2.0
            rs.AddTextDot(str(i), nc_FromIn.PointAt(t))


    bools_IsLineCrv = []
    bools_Is_seg_linear_nonLine = []
    devs = []
    tols_Needed = []

    for i, seg in enumerate(segs_All):
        if isinstance(seg, rg.LineCurve):
            bools_IsLineCrv.append(True)
            bools_Is_seg_linear_nonLine.append(False)
            continue

        rc = getPassingLineCurve(
            seg,
            bProcessArcs=bProcessArcs,
            fDevTol=fDevTol,
            fTolPerLengthUnit=fTolPerLengthUnit,
            fMinLineLength=fMinLineLength,
            )
        if rc and rc[0]:
            bools_IsLineCrv.append(False)
            bools_Is_seg_linear_nonLine.append(True)
            devs.append(rc[1])
        else:
            bools_IsLineCrv.append(False)
            bools_Is_seg_linear_nonLine.append(False)
            if rc[1] is not None:
                tols_Needed.append(rc[1])

    if not any(bools_Is_seg_linear_nonLine) and not any(bools_IsLineCrv):
        return None, None, "No linear segments in input (with current settings)."

    if not any(bools_Is_seg_linear_nonLine):
        return None, None, "No segments can be converted to lines (with current settings)."


    i = 0 # Segment index.
    idx_t_notConverted_Start = None
    idxs_ts_combined = [] # tuple(tuple(int(beginning nc segment index), int(ending nc segment index)), bool(True if indices are linear section))
    #bLinearFound = False
    bContiguousLinearFound = False

    while i < len(segs_All):
        if sc.escape_test(throw_exception=False):
            raise Exception("Break at index {}.".format(i))

        if bDebug: sEval = "i"; print(sEval,'=',eval(sEval))

        if not bools_Is_seg_linear_nonLine[i]:
            # Start of non-linear, contiguous group.

            if idx_t_notConverted_Start is None:
                idx_t_notConverted_Start = i

            i += 1
            continue

        if idx_t_notConverted_Start is not None:
            # Record non-converted section.
            idxs_ts_combined.append(((idx_t_notConverted_Start, i), False))
            idx_t_notConverted_Start = None

        # Record segment to convert.
        idxs_ts_combined.append(((i, i+1), True))

        i += 1

    # Record the last section.
    if idx_t_notConverted_Start is not None:
        idxs_ts_combined.append(((idx_t_notConverted_Start, i), False))


    # Proceed replacing segments with LineCurves.
    pc_Out = rg.PolyCurve()

    for (idxS, idxE), bLinear in idxs_ts_combined:
        if bLinear:
            line = rg.Line(nc_FromIn.PointAt(ts_knots[idxS]),
                nc_FromIn.PointAt(ts_knots[idxE]))
            if line.Length >= fMinLineLength:
                pc_Out.Append(line)
                #addGeoms(line)
                continue # to next tuples_splits.

        # Using rgC_In instead of nc_FromIn so as to not increase the degree of any NurbsCurves.
        seg_ForOut = rgC_In.Trim(
            rg.Interval(ts_knots[idxS], ts_knots[idxE]))
        pc_Out.Append(seg_ForOut)
        #addGeoms(seg_ForOut)

    # Success.

    bSuccess, plc = pc_Out.TryGetPolyline()
    if bSuccess:
        return (
            plc,
            max(devs) if devs else Rhino.RhinoMath.ZeroTolerance,
            None,
            )

    return (
        pc_Out,
        max(devs) if devs else Rhino.RhinoMath.ZeroTolerance,
        None,
        )


def processCurve(rgCrv_In, fDevTol=None, fTolPerLengthUnit=None, fMinLineLength=None, bProcessArcs=False, bMerge=True, bMergeThruSeam=True, bDebug=False):
    """
    Parameters:
        rgCrv_In
        bAlsoRatioTol
        fDevTol (May be ratio or absoluate based on value of bAlsoRatioTol)
        fMinLineLength
        bProcessArcs
        bDebug
    Returns on success: rg.Curve, float(Deviation), None
    Returns on deviation fail: None, float(Deviation), str(Description of fail)
    Returns on other fails: None, None, str(Description of fail)
    """


    if isinstance(rgCrv_In, rg.LineCurve):
        return None, None, "Skipped LineCurve."

    if not bProcessArcs and isinstance(rgCrv_In, rg.ArcCurve):
        return None, None, "Skipped ArcCurve."

    if rgCrv_In.IsArc(1e-6):
        return None, None, "Skipped curve passing IsArc."

    if not bMerge and isinstance(rgCrv_In, rg.PolylineCurve):
        return None, None, "PolylineCurve cannot be processed without allowing segment merging."


    if not rgCrv_In.IsClosed:

        rvs = getPassingLineCurve(
            rgCrv_In,
            bProcessArcs=bProcessArcs,
            fDevTol=fDevTol,
            fTolPerLengthUnit=fTolPerLengthUnit,
            fMinLineLength=fMinLineLength)

        rgC_Res, fDev, sLog = rvs

        if rgC_Res is not None:
            return rvs

        # Conversion of entire curve to line failed.
        # Try converting segments between knots.


    if isinstance(rgCrv_In, rg.PolylineCurve):
        rvs = convertSpans_of_PolylineCrv(
            rgCrv_In,
            fDevTol=fDevTol,
            fTolPerLengthUnit=fTolPerLengthUnit,
            fMinLineLength=fMinLineLength,
            bProcessArcs=bProcessArcs,
            bMergeThruSeam=bMergeThruSeam,
            bDebug=bDebug,
            )
    else:
        if bMerge:
            rvs = convertSpans_withMerging(
                rgCrv_In,
                fDevTol=fDevTol,
                fTolPerLengthUnit=fTolPerLengthUnit,
                fMinLineLength=fMinLineLength,
                bProcessArcs=bProcessArcs,
                bMergeThruSeam=bMergeThruSeam,
                bDebug=bDebug,
                )
        else:
            rvs = convertSpans_NoMerging(
                rgCrv_In,
                fDevTol=fDevTol,
                fTolPerLengthUnit=fTolPerLengthUnit,
                fMinLineLength=fMinLineLength,
                bProcessArcs=bProcessArcs,
                bDebug=bDebug,
                )

    rgC_Res, fDev, sLog = rvs

    if rgC_Res is None:
        if sLog is None:
            return None, None, "No simplification found."
        return None, None, sLog

    return rvs



    # Optional code to ramp up tolerance.

    p = -6
    c_Out = None

    tol_Sum = 0.0
    bTolLimitReached = False

    while True:
        tol = 10.0**(p)

        if tol >= (fDevTol - tol_Sum):
            tol = (fDevTol - tol_Sum)
            bTolLimitReached = True

        print(p, tol)
        
        if p == -2:
            addGeoms(c_Out)
            pass

        if sc.escape_test(throw_exception=False):
            raise Exception("Break at tolerance {}.".format(tol))


        rvs = convertSpans_withMerging(
            rgCrv_In if c_Out is None else c_Out,
            fDevTol=tol,
            fTolPerLengthUnit=fTolPerLengthUnit,
            fMinLineLength=fMinLineLength,
            bProcessArcs=bProcessArcs,
            bDebug=bDebug)

        if rvs[0] is not None:
            #addGeoms(rc[0])
            c_Out = rvs[0]

        if bTolLimitReached:
            1/0
            if c_Out is None:
                return None, None, "No simplification found."
            return c_Out, None, None

        tol_Sum += tol

        p += 1


def processCurveObject(rhCrv_In, **kwargs):
    """
    curvesAndEdges0 = (GUIDs of CurveObjects) and/or BrepEdges
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAlsoRatioTol = getOpt('bAlsoRatioTol')
    fDevTol = getOpt('fDevTol')
    fTolPerLengthUnit = getOpt('fTolPerLengthUnit')
    fMinLineLength = getOpt('fMinLineLength')
    bProcessArcs = getOpt('bProcessArcs')
    bMerge = getOpt('bMerge')
    bMergeThruSeam = getOpt('bMergeThruSeam')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rgC_In = rs.coercecurve(rhCrv_In)

    gC_In = rs.coerceguid(rhCrv_In)


    rc = processCurve(
        rgC_In,
        fDevTol=fDevTol,
        fTolPerLengthUnit=fTolPerLengthUnit if bAlsoRatioTol else None,
        fMinLineLength=fMinLineLength,
        bProcessArcs=bProcessArcs,
        bMerge=bMerge,
        bMergeThruSeam=bMergeThruSeam,
        bDebug=bDebug,
        )

    rgC_Res, fDev, sLog = rc

    if rgC_Res is None:
        if bEcho and sLog is not None: print(sLog)
        return rc

    if bReplace:
        if sc.doc.Objects.Replace(gC_In, rgC_Res):
            return gC_In, fDev, "Replaced a curve."
        else:
            s = "Curve could not be replaced with new geometry."
            if bEcho: print(s)
            return None, None, s

    # Add curve.
    gC_Out = sc.doc.Objects.AddCurve(rgCrv_Res)
    if gC_Out != Guid.Empty:
        return gCs_Out, fDev, "Added a curve."
    else:
        s = "New curve could not be added."
        if bEcho: print(s)
        return None, None, s


def main():

    rgCrvsAndEdges_In = getInput()
    if rgCrvsAndEdges_In is None: return


    bAlsoRatioTol = Opts.values['bAlsoRatioTol']
    fDevTol = Opts.values['fDevTol']
    fTolPerLengthUnit = Opts.values['fTolPerLengthUnit']
    fMinLineLength = Opts.values['fMinLineLength']
    bProcessArcs = Opts.values['bProcessArcs']
    bMerge = Opts.values['bMerge']
    bMergeThruSeam = Opts.values['bMergeThruSeam']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    fDevs = []
    sLogs = []

    for iC, objref in enumerate(rgCrvsAndEdges_In):
        gCrv_In = rs.coerceguid(objref)
        rgCrv_In = rs.coercecurve(objref)
        #print(rs.coercerhinoobject(objref))

        Rhino.RhinoApp.SetCommandPrompt(
            "Processing curve {}...".format(
                "" if len(rgCrvsAndEdges_In) == 1 else "{} of {} ".format(iC+1, len(rgCrvsAndEdges_In))))

        if isinstance(rgCrv_In, rg.BrepEdge):
            gCrv_In = sc.doc.Objects.AddCurve(rgCrv_In)
            rc = processCurveObject(
                gCrv_In,
                bAlsoRatioTol=bAlsoRatioTol,
                fDevTol=fDevTol,
                fTolPerLengthUnit=fTolPerLengthUnit,
                fMinLineLength=fMinLineLength,
                bProcessArcs=bProcessArcs,
                bMerge=bMerge,
                bMergeThruSeam=bMergeThruSeam,
                bReplace=bReplace,
                bEcho=False if len(rgCrvsAndEdges_In) > 1 else bEcho,
                bDebug=bDebug,
                )
            gRes, fDev, sLog = rc
            if gRes is None:
                sc.doc.Objects.Delete(objectId=gCrv_In, quiet=False)
            elif fDev is not None:
                fDevs.append(fDev)
            if sLog is not None:
                sLogs.append(sLog)
        elif gCrv_In is not None:
            rc = processCurveObject(
                gCrv_In,
                bEcho=False if len(rgCrvsAndEdges_In) > 1 else bEcho)
            gRes, fDev, sLog = rc
            if gRes is None:
                pass
            elif fDev is not None:
                fDevs.append(fDev)
            if sLog is not None:
                sLogs.append(sLog)
        else:
            s = "Not CurveObject or BrepEdge."
            if bEcho: print(s)
            continue

    if bEcho and len(rgCrvsAndEdges_In) > 1:
        for sLog in set(sLogs):
            print("[{}] {}".format(sLogs.count(sLog), sLog))

    if len(fDevs) == 1:
        if bEcho: print("Deviation: {}".format(getFormattedDistance(fDevs[0])))
    elif len(fDevs) > 1:
        if bEcho: print("Deviations: [{},{}]".format(
            getFormattedDistance(min(fDevs)),
            getFormattedDistance(max(fDevs))))


    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()