"""
This script is an alternative to _OffsetNormal in that it:
    1. Has options on how to prepare the input curves.
    2. Will offset within an input tolerance (Loose=No).
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
231227: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from clr import StrongBox
from System import Array


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bUseFaceOfSelNakedEdge'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLoose'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExplodePolyCrv'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRebuild'; keys.append(key)
    values[key] = True
    names[key] = 'TryToRebuild'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitAtNonG2Knots'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistance'; keys.append(key)
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 1.0
    else:
        values[key] = 10.0 * Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol'; keys.append(key)
    # 1/3 for offset, 1/3 for pipe, 1/3 for pipe trim.
    values[key] = max((0.33*sc.doc.ModelAbsoluteTolerance, 1e-6))
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAlignEndDirs'; keys.append(key)
    values[key] = False
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fTol':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

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


def _addCommonOptions(go):
    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    addOption('bExplodePolyCrv')
    addOption('bRebuild')
    addOption('bSplitAtNonG2Knots')
    addOption('fDistance')
    addOption('bAlignEndDirs')
    addOption('bLoose')
    if not Opts.values['bLoose']:
        addOption('fTol')
    addOption('bEcho')
    addOption('bDebug')

    return idxs_Opt


def _getInput_Curve():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve on face")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    #bPreselectedObjsChecked = False

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()
        addOption('bUseFaceOfSelNakedEdge')
        idxs_Opt.update(_addCommonOptions(go))

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref_CrvOnFace = go.Object(0)
            go.Dispose()

            return objref_CrvOnFace

        if res == ri.GetResult.Number:
            key = 'fDistance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_Face():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select base face")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    go.AcceptNumber(True, acceptZero=True)
    
    idxs_Opt = {}

    bPreselectedObjsChecked = False

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()
        idxs_Opt.update(_addCommonOptions(go))

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
            objref_Face = go.Object(0)

            go.Dispose()
    
            return objref_Face


        if res == ri.GetResult.Number:
            key = 'fDistance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_Click():
    """
    Click to toggle angle and/or direction with optional input.

    Returns:
        True: To recalculate and reloop
        False: To not recalculate and break out of loop with current output.
        None: To not recalculate and return without output.
    """

    go = ri.Custom.GetPoint()

    go.SetCommandPrompt("Left click to flip direction")

    go.SetCommandPromptDefault("Accept result")

    go.AcceptNumber(True, acceptZero=True)
    go.AcceptNothing(True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    key = 'FlipDir'; idxs_Opt[key] = go.AddOption(key)

    idxs_Opt.update(_addCommonOptions(go))

    res = go.Get()

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Nothing:
        go.Dispose()
        return False

    if res == ri.GetResult.Point:
        Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
        Opts.setValue('fDistance')
        go.Dispose()
        return True

    if res == ri.GetResult.Number:
        key = 'fDistance'
        Opts.riOpts[key].CurrentValue = go.Number()
        Opts.setValue(key)
        go.Dispose()
        return True

    # An option was selected.

    if go.OptionIndex() == idxs_Opt['FlipDir']:
        Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
        Opts.setValue('fDistance')
        go.Dispose()
        return True

    for key in idxs_Opt:
        if go.Option().Index == idxs_Opt[key]:
            Opts.setValue(key, go.Option().CurrentListOptionIndex)
            break

    go.Dispose()
    return True


def _crvWithSpansCompletelyOnFace(rgCrv, rgFace, t_Crv_Pick, fTol, bDebug=False):
    """
    Only process spans of the curve whose spans start and ends are on the Face.
    """
    iSpans_OnFace = []
    for iSpan in xrange(rgCrv.SpanCount):
        t = rgCrv.SpanDomain(iSpan).T0
        pt_OnCrv = rgCrv.PointAt(t)
        bSuccess, u, v = rgFace.ClosestPoint(pt_OnCrv)
        if not bSuccess: continue
        pt_OnFace = rgFace.PointAt(u,v)
        dist = pt_OnFace.DistanceTo(pt_OnCrv)
        if dist > fTol:
            if bDebug:
                sEval='dist'; print(sEval+': ',eval(sEval))
                print("PointAtStart not on underlying surface.")
            continue
        
        t = rgCrv.SpanDomain(iSpan).T1
        pt_OnCrv = rgCrv.PointAt(t)
        bSuccess, u, v = rgFace.ClosestPoint(pt_OnCrv)
        if not bSuccess: continue
        pt_OnFace = rgFace.PointAt(u,v)
        dist = pt_OnFace.DistanceTo(pt_OnCrv)
        if dist > fTol:
            if bDebug:
                #sc.doc.Objects.AddPoint(pt_OnFace)
                sEval='dist'; print(sEval+': ',eval(sEval))
                print("PointAtEnd not on underlying surface.")
            continue
        
        iSpans_OnFace.append(iSpan)
    
    if not iSpans_OnFace:
        print("None of the spans of the curve are completely on the face.")
        rgFace.Brep.Dispose()
        return

    if len(iSpans_OnFace) == rgCrv.SpanCount:
        return rgCrv.Duplicate()

    if len(iSpans_OnFace) == 1:
        return rgCrv.Trim(
            rgCrv.SpanDomain(iSpans_OnFace[0]).T0,
            rgCrv.SpanDomain(iSpans_OnFace[0]).T1)

    # Create nested lists of contiguous spans so that curves that go off and
    # on the face can be correctly processed.
    iSpans_Contiguous_nests = [[iSpans_OnFace[0]]]
    for iSpan in iSpans_OnFace[1:]:
        if iSpan == iSpans_Contiguous_nests[-1][-1] + 1:
            iSpans_Contiguous_nests[-1].append(iSpan)
        else:
            iSpans_Contiguous_nests.append([iSpan])
        
    if rgCrv.IsClosed and len(iSpans_Contiguous_nests) > 1:
        if (
                iSpans_Contiguous_nests[0][0] == 0 and
                iSpans_Contiguous_nests[-1][-1] == rgCrv.SpanCount-1
        ):
            iSpans_Contiguous_nests[0] = iSpans_Contiguous_nests[-1] + iSpans_Contiguous_nests[0]
            iSpans_Contiguous_nests.pop()
    
    if len(iSpans_Contiguous_nests) == 1:
        rgC_Out = rgCrv.Trim(
                rgCrv.SpanDomain(iSpans_Contiguous_nests[0][0]).T0,
                rgCrv.SpanDomain(iSpans_Contiguous_nests[0][-1]).T1)
    elif len(iSpans_Contiguous_nests) > 1:
        for iSpan_NestIndex, iSpans_Contiguous in enumerate(iSpans_Contiguous_nests):
            for iSpan in iSpans_Contiguous:
                if rgCrv.SpanDomain(iSpan).T0 <= t_Crv_Pick <= rgCrv.SpanDomain(iSpan).T1:
                    rgC_Out = rgCrv.Trim(
                            rgCrv.SpanDomain(iSpans_Contiguous_nests[iSpan_NestIndex][0]).T0,
                            rgCrv.SpanDomain(iSpans_Contiguous_nests[iSpan_NestIndex][-1]).T1)
                else:
                    print("Curve was not picked within the face.")
                    rgCrv.Dispose()
                    return
        
    if bDebug:
        sc.doc.Objects.AddCurve(rgC_Out)

    return rgC_Out


def _getOffsetDeviation(crvA, crvB, fTarget):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None

    fDev_Max = rc[1]
    fDev_Min = rc[4]

    return max((abs(abs(fTarget)-fDev_Max), abs(abs(fTarget)-fDev_Min)))


def _do_curves_deviate_within_tolerance(crvA, crvB, fTarget, fTol, fSamplingDist=None, bDebug=False):
    """
    Uses Curve.ClosestPoint instead of Curve.GetDistancesBetweenCurves.
    """

    if fSamplingDist is None:
        fSamplingDist = 10.0 * fTol

    fLimit_Min = abs(fTarget) - fTol
    fLimit_Max = abs(fTarget) + fTol

    strongBox_points = StrongBox[Array[rg.Point3d]]()

    rc = crvA.DivideByLength(
        segmentLength=fSamplingDist,
        includeEnds=True,
        points=strongBox_points)

    pts = list(strongBox_points.Value)
    if bDebug: sEval = "len(pts)"; print("{}: {}".format(sEval, eval(sEval)))

    devs = []

    for ptA in pts:
        bSuccess, t = crvB.ClosestPoint(ptA)
        if not bSuccess:
            raise ValueError("Curve.ClosestPoint failed.")
        ptB = rg.Curve.PointAt(crvB, t)
        dev = ptA.DistanceTo(ptB)
        devs.append(dev)
        #if bDebug: sc.doc.Objects.AddLine(ptA, ptB)
        if dev < fLimit_Min:
            return False
        if dev > fLimit_Max:
            return False

    #sEval = "rg.Curve.GetDistancesBetweenCurves(crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)"
    #print("{}: {}".format(sEval, eval(sEval)))

    # For debugging.
    if bDebug:
        #fDev_Min = min(devs)
        sEval = "min(devs)"; print("{}: {}".format(sEval, eval(sEval)))
        #fDev_Max = max(devs)
        sEval = "max(devs)"; print("{}: {}".format(sEval, eval(sEval)))

    return True


def _getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None

    return rc[1]


def _matchCrvEndDirs(nc_ToMod, nc_Ref):
    """
    nc_ToMod is modified.
    """
    bSuccess = nc_ToMod.SetEndCondition(
        bSetEnd=False,
        continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
        point=nc_ToMod.PointAtStart,
        tangent=nc_Ref.TangentAtStart)
    if not bSuccess:
        print("SetEndCondition failed.")
        return False
    bSuccess = nc_ToMod.SetEndCondition(
        bSetEnd=True,
        continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
        point=nc_ToMod.PointAtEnd,
        tangent=nc_Ref.TangentAtEnd)
    if not bSuccess:
        print("SetEndCondition failed.")
        return False

    return True


def _rebuildCrv(rgCrv_In, fTol_Simplify, bDebug=False):
    """
    Returns:
        Curve, float(deviation)
        None for no Rebuild.
    """

    #if isinstance(rgCrv_In, rg.PolylineCurve):
    #    raise Exception("PolylineCurve is not supported.")
    if isinstance(rgCrv_In, (rg.LineCurve, rg.ArcCurve)):
        return

    nc_WIP = rgCrv_In.ToNurbsCurve()

    if nc_WIP.SpanCount == 1:
        if bDebug: print("Curve already has only 1 span.")
        return


    if nc_WIP.IsPeriodic and (nc_WIP.Knots.KnotStyle == rg.KnotStyle.Uniform):
        sEval = "nc_WIP.Knots.Count"; print("{}: {}".format(sEval, eval(sEval)))
        if bDebug: print("Curve is already periodic and uniform.")
        return


    # Try to rebuild as a Bezier.
    for degree in 2,3,5:
        nc_Rebuilt = nc_WIP.Rebuild(
            pointCount=degree+4,
            degree=degree,
            preserveTangents=True)

        dev = _getDistancesBetweenCurves(nc_WIP, nc_Rebuilt)
        if dev <= fTol_Simplify:
            nc_WIP.Dispose()
            return nc_Rebuilt, dev

        nc_Rebuilt.Dispose()


    if nc_WIP.Degree == 3:
        if (nc_WIP.Knots.KnotStyle in (
            rg.KnotStyle.QuasiUniform,
            rg.KnotStyle.Uniform)
        ):
            if bDebug: print("Curve is already {}.".format(nc_WIP.Knots.KnotStyle))
            return


    # Try to rebuild as degree-3 uniform.
    pointcount = 4
    while True:
        nc_Rebuilt = nc_WIP.Rebuild(
            pointCount=pointcount,
            degree=3,
            preserveTangents=True)

        dev = _getDistancesBetweenCurves(nc_WIP, nc_Rebuilt)
        if dev <= fTol_Simplify:
            nc_WIP.Dispose()
            return nc_Rebuilt, dev

        nc_Rebuilt.Dispose()

        pointcount += 1
        if pointcount > 103:
            return


def _split_NurbsCrv_at_nonG2_knots(nc_In):
    """
    Always returns a list of curves.
    """

    if not isinstance(nc_In, rg.NurbsCurve):
        return [nc_In]

    ts_polyknots = []

    if nc_In.IsPeriodic:
        iKs = range(nc_In.Knots.Count)
    elif nc_In.IsClosed:
        iKs = range(nc_In.Knots.Count - nc_In.Degree)
    else:
        iKs = range(nc_In.Degree, nc_In.Knots.Count - nc_In.Degree)

    for iK in iKs:
        if nc_In.Knots.KnotMultiplicity(iK) > (nc_In.Degree-2):
            ts_polyknots.append(nc_In.Knots[iK])

    if not ts_polyknots:
        return [nc_In.DuplicateCurve()]

    rc = nc_In.Split(ts_polyknots)

    if not rc:
        print("Splitting at non-G2 knots failed.  Check input.")
        return [nc_In]

    return rc


def _prepareCrvToOffset(rgCrv_In, bExplodePolyCrv, bRebuild, bSplitAtNonG2Knots, bMakeDeformable, fTol, bDebug=False):
    if bExplodePolyCrv and isinstance(rgCrv_In, rg.PolyCurve):
        ncs_WIP = [_.ToNurbsCurve() for _ in rgCrv_In.Explode()]
    else:
        ncs_WIP = [rgCrv_In.ToNurbsCurve()]

    fDev = 0.0

    if bRebuild:
        # Rebuilding to a smaller tolerance to reduce toleranc stackup.
        fTol_Simplify = max((0.5*fTol, 1e-4))
        if bDebug: sEval = "fTol_Simplify"; print("{}: {}".format(sEval, eval(sEval)))

        ncs_Simplified = []
        fDevs_thisRebuild = []
        for nc_WIP in ncs_WIP:
            rc = _rebuildCrv(nc_WIP, fTol_Simplify=fTol_Simplify, bDebug=bDebug)
            if rc:
                nc_WIP.Dispose()
                ncs_Simplified.append(rc[0])
                fDevs_thisRebuild.append(rc[1])
            else:
                ncs_Simplified.append(nc_WIP)
        ncs_WIP = ncs_Simplified
        if fDevs_thisRebuild:
            fDev = max(fDevs_thisRebuild)

    if bSplitAtNonG2Knots:
        ncs_Split = []
        for nc_WIP in ncs_WIP:
            rc = _split_NurbsCrv_at_nonG2_knots(nc_WIP)
            ncs_Split.extend(rc)
            nc_WIP.Dispose()
        ncs_WIP = ncs_Split

    if bRebuild:
        ncs_Simplified = []
        fDevs_thisRebuild = []
        for nc_WIP in ncs_WIP:
            rc = _rebuildCrv(nc_WIP, fTol_Simplify=fTol_Simplify, bDebug=bDebug)
            if rc:
                nc_WIP.Dispose()
                ncs_Simplified.append(rc)
                fDevs_thisRebuild.append(rc[1])
            else:
                ncs_Simplified.append(nc_WIP)
        ncs_WIP = ncs_Simplified
        if fDevs_thisRebuild:
            fDev += max(fDevs_thisRebuild)

    if bMakeDeformable:
        for i in range(len(ncs_WIP)):
            if ncs_WIP[i].Degree < 3 and ncs_WIP[i].Points.Count < 4:
                ncs_WIP[i].IncreaseDegree(3)

    return ncs_WIP, fDev


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.crvs = []
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        for crv in self.crvs:
            bbox = crv.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PreDrawObjects(self, drawEventArgs):

        color = sc.doc.Layers.CurrentLayer.Color

        for crv in self.crvs:
            drawEventArgs.Display.DrawCurve(
                curve=crv,
                color=color,
                thickness=self.crv_thk)


def _createNurbsCurve_with_more_knots(nc_In, iQtyToAddInEachSpan):
    """
    Knots are added uniformly with each span.

    Parameters:
        nc_In: rg.NurbsSurface,
        iQtyToAddInEachSpan: int

    Returns: rg.NurbsSurface
    """

    if iQtyToAddInEachSpan < 1:
        return

    nc_Out = nc_In.ToNurbsCurve()

    # Add knots from Domain end to beginning.

    for iK in range(nc_In.Knots.Count-1, 0, -1):

        k_R = nc_In.Knots[iK]
        k_L = nc_In.Knots[iK-1]

        for i in range(iQtyToAddInEachSpan):
            fraction_from_R = float(i+1) / float(iQtyToAddInEachSpan+1)
            k_M = fraction_from_R*k_R + (1.0-fraction_from_R)*k_L
            nc_Out.Knots.InsertKnot(k_M)

    return nc_Out


def createOffsetCurve(rgCrv_In, rgSrf, bLoose, bAlignEndDirs, fDistance, fTol, fSamplingDist, bDebug):

    sEval = "fTol"; print("{}: {}".format(sEval, eval(sEval)))

    nc_Offset = rgCrv_In.OffsetNormalToSurface(
        surface=rgSrf, height=fDistance)

    if (bAlignEndDirs and
        not _matchCrvEndDirs(nc_Offset, rgCrv_In)
    ):
        raise Exception("Alignment of end tangent failed.")

    if bLoose:
        return nc_Offset


    dev = _getOffsetDeviation(rgCrv_In, nc_Offset, fDistance)
    if bDebug: sEval = "dev"; print("{}: {}".format(sEval, eval(sEval)))
    if dev <= fTol:
        return nc_Offset

    #if _do_curves_deviate_within_tolerance(rgCrv_In, nc_Offset, fDistance, fTol, fSamplingDist, bDebug=bDebug):
    #    return nc_Offset

    nc_Offset.Dispose()

    # Add knots and try again.
    i = 0
    while True:
        sc.escape_test()

        i += 1
        if i == 51:
            return

        if bDebug: sEval = "i"; print("{}: {}".format(sEval, eval(sEval)))

        nc_ToOffset = _createNurbsCurve_with_more_knots(
            rgCrv_In,
            iQtyToAddInEachSpan=i)

        nc_Offset = nc_ToOffset.OffsetNormalToSurface(
            surface=rgSrf, height=fDistance)

        if (bAlignEndDirs and
            not _matchCrvEndDirs(nc_Offset, rgCrv_In)
        ):
            raise Exception("Alignment of end tangent failed.")

        dev = _getOffsetDeviation(rgCrv_In, nc_Offset, fDistance)
        if bDebug: sEval = "dev"; print("{}: {}".format(sEval, eval(sEval)))
        if dev <= fTol:
            print("Added {} knots between spans.".format(i))
            return nc_Offset

        #if _do_curves_deviate_within_tolerance(
        #    rgCrv_In, nc_Offset, fDistance, fTol, fSamplingDist,
        #    bDebug=bDebug
        #):
        #    return nc_Offset

        nc_ToOffset.Dispose()
        nc_Offset.Dispose()


def _createGeometryInteractively():
    """
    """

    objref_CrvToOffset = _getInput_Curve()
    if objref_CrvToOffset is None: return


    bUseFaceOfSelNakedEdge = Opts.values['bUseFaceOfSelNakedEdge']


    rgEdge = objref_CrvToOffset.Edge()

    if rgEdge and bUseFaceOfSelNakedEdge and rgEdge.Valence == rg.EdgeAdjacency.Naked:
        idxF = objref_CrvToOffset.Edge().AdjacentFaces()[0]
        rgF_In = rgEdge.Brep.Faces[idxF]
    else:
        sc.doc.Objects.UnselectAll()

        objref_Face = _getInput_Face()
        if objref_Face is None: return

        sc.doc.Objects.UnselectAll()


        rgF_In = objref_Face.Face()


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    bLoose = Opts.values['bLoose']
    bAlignEndDirs = Opts.values['bAlignEndDirs']
    bExplodePolyCrv = Opts.values['bExplodePolyCrv']
    bSplitAtNonG2Knots = Opts.values['bSplitAtNonG2Knots']
    bRebuild = Opts.values['bRebuild']
    fDistance = Opts.values['fDistance']
    fTol = Opts.values['fTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rgC_In, t_Crv0_Pick = objref_CrvToOffset.CurveParameter()


    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()


    rgC_In_TrimmedToFace = _crvWithSpansCompletelyOnFace(
        rgC_In, rgF_In, t_Crv0_Pick, fTol, bDebug)
    if rgC_In_TrimmedToFace is None: return

    if (
        bDebug and
        isinstance(rgC_In, rg.NurbsCurve) and
        isinstance(rgC_In_TrimmedToFace, rg.NurbsCurve)
    ):
        sEval = "rgC_In.EpsilonEquals(rgC_In_TrimmedToFace, 1e-6)"; print("{}: {}".format(sEval, eval(sEval)))

    if rgC_In_TrimmedToFace.IsClosed and fAngle_End_Deg:
        fAngle_End_Deg = Opts.values['fAngle_End_Deg'] = sc.sticky[Opts.stickyKeys['fAngle_End_Deg']] = None
        bVariableAngle = Opts.values['bVariableAngle'] = sc.sticky[Opts.stickyKeys['bVariableAngle']] = False


    sk_conduit = 'conduit({})'.format(__file__) # StickyKey
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit


    fSamplingDist = 100.0*sc.doc.ModelAbsoluteTolerance
    while True:
        sc.escape_test()

        ncs_toOffset, fDev_fromRebuilds = _prepareCrvToOffset(
            rgC_In_TrimmedToFace,
            bExplodePolyCrv=bExplodePolyCrv,
            bRebuild=bRebuild,
            bSplitAtNonG2Knots=bSplitAtNonG2Knots,
            bMakeDeformable=bAlignEndDirs,
            fTol=fTol,
            bDebug=bDebug)

        if bRebuild:
            sEval = "fDev_fromRebuilds"; print("{}: {}".format(sEval, eval(sEval)))


        ncs_FinEnd = []


        for nc_toOffset in ncs_toOffset:

            if bDebug: sEval = "nc_toOffset.SpanCount"; print("{}: {}".format(sEval, eval(sEval)))

            nc_FinEnd = createOffsetCurve(
                rgCrv_In=nc_toOffset,
                rgSrf=rgF_In,
                bLoose=bLoose,
                bAlignEndDirs=bAlignEndDirs,
                fDistance=fDistance,
                fTol=fTol-fDev_fromRebuilds,
                fSamplingDist=fSamplingDist,
                bDebug=bDebug)
            if nc_FinEnd is None:
                if bEcho: print("No solution found.")
                return

            ncs_FinEnd.append(nc_FinEnd)

            nc_toOffset.Dispose()



        conduit.crvs = ncs_FinEnd

        conduit.Enabled = True

        sc.doc.Views.Redraw()

        if bEcho:
            sOut = []
            if len(ncs_FinEnd) > 1: sOut.append("{} curves".format(len(ncs_FinEnd)))
            if sOut:
                print("Calculated {}.".format(", ".join(sOut)))


        rc = _getInput_Click()

        conduit.Enabled = False

        if rc is None:
            for _ in ncs_FinEnd: _.Dispose()
            return

        if not rc:
            return (
                ncs_FinEnd,
                bEcho)


        for _ in ncs_FinEnd: _.Dispose()



        bLoose = Opts.values['bLoose']
        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bExplodePolyCrv = Opts.values['bExplodePolyCrv']
        bRebuild = Opts.values['bRebuild']
        bSplitAtNonG2Knots = Opts.values['bSplitAtNonG2Knots']
        fDistance = Opts.values['fDistance']
        fTol = Opts.values['fTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


def main():

    rc = _createGeometryInteractively()
    if rc is None: return

    rgCs, bEcho = rc

    gCs = []
    for rgC in rgCs:
        gC = sc.doc.Objects.AddCurve(rgC)
        if gC != gC.Empty:
            gCs.append(gC)

    if bEcho:
        sOut = []
        if gCs: sOut.append("{} offset curves".format(len(gCs)))
        print("Added {}".format(", ".join(sOut)))


if __name__ == '__main__': main()
