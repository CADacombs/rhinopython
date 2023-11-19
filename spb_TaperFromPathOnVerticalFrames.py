"""
This script is an alternative to _ExrudeCrvTapered.

Although it doesn't have all _ExrudeCrvTapered's options, e.g., Corners,
it has a few additional features:
    It can create loose extrudes.
    It can create variable tapers, where the start and end angles are defined per extrusion segment.
    It can align the ends of the tapered-to curves with the curves to taper-from curves
        to aid in constructing contiguous breps with G1 continuity.

Limitations:
    The variable taper only works per segment, not through an entire polycurve.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190604: Created from a split from another script.  Added bAddTaperEndCrvs.
190624: Now, command waits for Enter or options when curves are preselected.
...
220828: Added an option.  Refactored.
231117-18: Added bAlignEndDirs.  Removed bRebuildPath.  Added bAddTaperStartCrvs.
        Replaced an import with some of its code, then simplified it.

TODO: Correctly create curve opposite path and brep when path curve is closed.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math

from System import Enum


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bNumResForDist_NotAngle'; keys.append(key)
    values[key] = True
    names[key] = 'NumberEntry'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'DraftAngle', 'Dist')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistance'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bProjDist'; keys.append(key)
    values[key] = True
    names[key] = 'DistType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'True', 'Projected')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTaper_Start_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'TaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bVariableTaper'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTaper_End_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'EndTaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTaperChangePerCrvParam'; keys.append(key)
    values[key] = False
    names[key] = 'TaperChangePerCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Length', 'Param')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCPlane'; keys.append(key)
    values[key] = True
    names[key] = 'PlanView'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitPathsAtG2PlusKnots'; keys.append(key)
    values[key] = False
    names[key] = 'OnlyBezierOut'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAlignEndDirs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtGrevilles'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtEqualDivisions'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDivisionCt'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(
            initialValue=values[key],
            setLowerLimit=True,
            limit=2)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddTaperStartCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddTaperedLines'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddTaperEndCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddBrep'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iBrepMethod'; keys.append(key)
    listValues[key] = 'Loft2Crvs', 'LoftSectionLines', 'Sweep2A', 'Sweep2B', 'Network' # All items must be strings.
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftType'; keys.append(key)
    listValues[key] = Enum.GetNames(rg.LoftType)
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fBrepTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=1e-6)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

        if key == 'fDistance':
            if -1e-6 < cls.riOpts[key].CurrentValue < 1e-6:
                print("fDistance input is too small.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fBrepTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
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


def getInput(bFirstGetObjects=False):
    """
    Get objects with optional input.
    
    Returns
        None to cancel.
        tuple of (list(ObjRefs), bool(Generate new geometry), bool(Accept results):
    """


    def wasThereChangeInObjectsSelected(gCrvs_PreSelctd, gBreps_PreSelctd, idxs_Edges_PerBrep):
        objrefs = go.Objects()
        for objref in objrefs:
            rgCrv = objref.Curve()
            if isinstance(rgCrv, rg.BrepEdge):
                rgCrv.Dispose()
                if not objref.ObjectId in gBreps_PreSelctd:
                    return True
                zipped = zip(gBreps_PreSelctd, idxs_Edges_PerBrep)
                for gBrep_PreSelctd, idx_Edge_PerBrep in zipped:
                    if gBrep_PreSelctd == objref.ObjectId:
                        if idx_Edge_PerBrep == idx_Edge_PerBrep:
                            break
                else:
                    # Edge not found.
                    return True
            else:
                # Curves other than BrepEdges.
                rgCrv.Dispose()
                if not objref.ObjectId in gCrvs_PreSelctd:
                    return True
        return False


    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select path curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    # Due to the object collection below, do not use AcceptNothing.
    go.AcceptNumber(True, acceptZero=False)


    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    addOption('bNumResForDist_NotAngle')
    addOption('fDistance')
    if not Opts.values['bVariableTaper']:
        addOption('bProjDist')
    addOption('fTaper_Start_Deg')
    addOption('bVariableTaper')
    if Opts.values['bVariableTaper']:
        addOption('fTaper_End_Deg')
        idxs_Opt['SwapAngles'] = go.AddOption('SwapAngles')
        addOption('bTaperChangePerCrvParam')
    idxs_Opt['FlipDir'] = go.AddOption('FlipDir')
    idxs_Opt['FlipAngle'] = go.AddOption('FlipAngle')
    addOption('bCPlane')
    addOption('bSplitPathsAtG2PlusKnots')
    addOption('bAlignEndDirs')
    addOption('bAtGrevilles')
    addOption('bAtKnots')
    addOption('bAtEqualDivisions')
    if Opts.values['bAtEqualDivisions']:
        addOption('iDivisionCt')
    addOption('bAddTaperStartCrvs')
    addOption('bAddTaperedLines')
    addOption('bAddTaperEndCrvs')
    addOption('bAddBrep')
    if Opts.values['bAddBrep']:
        addOption('iBrepMethod')
        if Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'LoftSectionLines':
            addOption('iLoftType')
        addOption('fBrepTol')
    addOption('bEcho')
    addOption('bDebug')


    go.AlreadySelectedObjectSelect = True # So objects can be reselected after being unselected in same go.GetMultiple.
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

    if not go.ObjectsWerePreselected:
        objrefs = None
        iCt_Crvs_PreSelctd = 0
    else:
        if bFirstGetObjects:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs, True, False

        iCt_Crvs_PreSelctd = go.ObjectCount
        objrefs = go.Objects()
        gCrvs_PreSelctd = []
        gBreps_PreSelctd = []
        idxs_Edges_PerBrep = []
        for objref in objrefs:
            rgCrv = objref.Curve()
            if isinstance(rgCrv, rg.BrepEdge):
                gBreps_PreSelctd.append(objref.ObjectId)
                idxs_Edges_PerBrep.append(rgCrv.EdgeIndex)
            else:
                # Curves other than BrepEdges.
                gCrvs_PreSelctd.append(objref.ObjectId)
            rgCrv.Dispose()

        go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Object:
        objrefs = go.Objects()
        if iCt_Crvs_PreSelctd == go.ObjectCount:
            if not wasThereChangeInObjectsSelected(
                    gCrvs_PreSelctd=gCrvs_PreSelctd,
                    gBreps_PreSelctd=gBreps_PreSelctd,
                    idxs_Edges_PerBrep=idxs_Edges_PerBrep
            ):
                go.Dispose()
                return objrefs, False, True
        go.Dispose()
        return objrefs, True, False

    if res == ri.GetResult.Number:
        key = 'fDistance' if Opts.values['bNumResForDist_NotAngle'] else 'fTaper_Start_Deg'
        Opts.riOpts[key].CurrentValue = go.Number()
        Opts.setValue(key)
    elif Opts.values['bVariableTaper'] and go.OptionIndex() == idxs_Opt['SwapAngles']:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue, Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                Opts.riOpts['fTaper_End_Deg'].CurrentValue, Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
        Opts.setValue('fTaper_Start_Deg')
        Opts.setValue('fTaper_End_Deg')
    elif go.OptionIndex() == idxs_Opt['FlipAngle']:
        for key in 'fTaper_Start_Deg', 'fTaper_End_Deg':
            Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
            Opts.setValue(key)
    elif go.OptionIndex() == idxs_Opt['FlipDir']:
        for key in 'fDistance', 'fTaper_Start_Deg', 'fTaper_End_Deg':
            Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
            Opts.setValue(key)
    elif Opts.values['bAddBrep'] and go.OptionIndex() == idxs_Opt['iBrepMethod']:
        Opts.setValue('iBrepMethod', go.Option().CurrentListOptionIndex)
        #Opts.values['iBrepMethod'] = go.Option().CurrentListOptionIndex
    elif Opts.values['bAddBrep'] and Opts.values['iBrepMethod'] == 1 and go.OptionIndex() == idxs_Opt['iLoftType']:
        #Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
        Opts.setValue('iLoftType', go.Option().CurrentListOptionIndex)
    else:
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break

    return objrefs, True, False


def continuityVectorsAt(nc, t, side=rg.CurveEvaluationSide.Default):
    """
    Returns: Tuple of these 4 items:
        3D point as a vector,
        Unit tangent vector,
        Curvature vector,
        Vector for comparing G3 continuity, not the G3 vector itself
        
        None for any of aforementioned vectors on fail.
    """
    
    if not isinstance(nc, rg.NurbsCurve): return
    
    vs = nc.DerivativeAt(
            t,
            derivativeCount=3,
            side=side)
    
    # Not using rg.Curve.TangentAt since it doesn't take into account CurveEvaluationSide.
    vTangency = vs[1]/vs[1].Length

    cross = rg.Vector3d.CrossProduct

    # For R3
    vCurvature = cross(cross(vs[1], vs[2]), vs[1]) / vs[1].Length**4

    return vs[0], vTangency, vCurvature


def getG2Discontinuities(rgCrv_In):
    """
    """

    fAngleTol_Rad = sc.doc.ModelAngleToleranceRadians
    fRadiusTol = sc.doc.ModelAbsoluteTolerance

    nc = rgCrv_In.ToNurbsCurve()

    t0 = nc.Domain.T0
    t1 = nc.Domain.T1

    ts_discontinuities = []

    bG3_discontinuousFound = False

    # This will also skip simple knot overlaps of Periodic curves.
    if nc.IsClosed and not nc.IsPeriodic:
        iK = 0
    else:
        iK = nc.Degree
    iK_Stop = nc.Knots.Count - nc.Degree

    while iK < iK_Stop:
        sc.escape_test()

        m = nc.Knots.KnotMultiplicity(iK)

        if m <= nc.Degree - 3:
            # Continuity at this knot is at least G3.
            iK += m
            continue

        if iK == 0:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                ) = continuityVectorsAt(
                    nc,
                    nc.Knots[nc.Knots.Count - 1],
                    rg.CurveEvaluationSide.Below)
        else:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                ) = continuityVectorsAt(
                    nc,
                    nc.Knots[iK],
                    rg.CurveEvaluationSide.Below)


        (
            vG0_Above,
            vG1_Above,
            vG2_Above,
            ) = continuityVectorsAt(
                nc,
                nc.Knots[iK],
                rg.CurveEvaluationSide.Above)


        # G1.

        if m > nc.Degree - 1:

            iParallel = vG1_Below.IsParallelTo(
                other=vG1_Above, angleTolerance=fAngleTol_Rad)

            if not iParallel == 1:
                ts_discontinuities.append(nc.Knots[iK])
                iK += m
                continue


        # G2.

        if m > nc.Degree - 2:

            if vG2_Below.IsTiny() and vG2_Above.IsTiny():
                # Linear.
                pass
            else:

                iParallel = vG2_Below.IsParallelTo(
                    other=vG2_Above, angleTolerance=fAngleTol_Rad)

                fCrvDelta_Degrees = Rhino.RhinoMath.ToDegrees(
                    rg.Vector3d.VectorAngle(vG2_Below, vG2_Above))

                if iParallel != 1:
                    ts_discontinuities.append(nc.Knots[iK])
                    iK += m
                    continue
                else:
                    kappa_Below = vG2_Below.Length
                    kappa_Above = vG2_Above.Length
                    ratio_of_curvature = (abs(kappa_Below-kappa_Above) /
                                        max(kappa_Below, kappa_Above))
                    delta_radius = abs(1.0/kappa_Below-1.0/kappa_Above)
                    if delta_radius > fRadiusTol:
                        ts_discontinuities.append(nc.Knots[iK])
                        iK += m
                        continue

        iK += m

    nc.Dispose()

    return ts_discontinuities


def _prepareCurves(rgCrvs_In, bSplitPathsAtG2PlusKnots=False, bMakeDeformable=True):
    """
    Prepare curves, including splitting per options.

    Returns: list of (profile) lists of new curves
    """

    # JoinCurves groups curves for further processing
    # and also aligns direction within each PolyCurve.
    rgCrvs_Joined = rg.Curve.JoinCurves(
        rgCrvs_In,
        joinTolerance=sc.doc.ModelAbsoluteTolerance)
    if not rgCrvs_Joined: return

    rgCs_Out_Nested = [] # Per connected path.

    for rgCrv_Joined in rgCrvs_Joined:

        if isinstance(rgCrv_Joined, rg.PolyCurve):
            rgCrv_Joined.RemoveNesting()

        if isinstance(rgCrv_Joined, rg.PolyCurve):
            rc = rgCrv_Joined.Explode()
            if not rc:
                print("Could not Explode PolyCurve.")
            rgCrvs_ExplodedPoly = []
            for c in rc:
                rgCrvs_ExplodedPoly.append(c.ToNurbsCurve())
                c.Dispose()
        else:
            # not PolyCurve.
            rgCrvs_ExplodedPoly = [rgCrv_Joined.ToNurbsCurve()]

        rgCrv_Joined.Dispose()

        # Split each curve at each G2 discontinuity.
        rgCs_SplitAtNonG2_AllPolycrvs = []
        for rgCrv_SplitPoly in rgCrvs_ExplodedPoly:
            ts = getG2Discontinuities(rgCrv_SplitPoly)
            rgCs_SplitAtNonG2_ThisPolycrv = rgCrv_SplitPoly.Split(
                [rgCrv_SplitPoly.Domain.T0] + ts + [rgCrv_SplitPoly.Domain.T1])
            rgCs_SplitAtNonG2_AllPolycrvs.extend(rgCs_SplitAtNonG2_ThisPolycrv)

        if not bSplitPathsAtG2PlusKnots:
            rgCs_AfterG2PlusSplits = rgCs_SplitAtNonG2_AllPolycrvs[:] # Referencing same Geometry, so do not Dispose yet.
        else:
            for c in rgCs_SplitAtNonG2_AllPolycrvs:
                ts_SpanBoundaries = [c.Domain.T0]
                for iSpan in xrange(c.SpanCount):
                    ts_SpanBoundaries.append(c.SpanDomain(iSpan).T1)
                rgCrvs_SplitAtKnots = c.Split(ts_SpanBoundaries)
                rgCs_AfterG2PlusSplits.extend(rgCrvs_SplitAtKnots)
                c.Dispose()

        rgCs_Out_1Profile = []

        if not bMakeDeformable:
            rgCs_Out_1Profile.extend(rgCs_SplitAtNonG2_AllPolycrvs)
        else:
            for c in rgCs_SplitAtNonG2_AllPolycrvs:
                if c.Points.Count != 3:
                    rgCs_Out_1Profile.append(c)
                else:
                    if not c.IncreaseDegree(3):
                        print("Could not raise degree of {}.".format(c.GetType().Name))
                    rgCs_Out_1Profile.append(c)

        rgCs_Out_Nested.append(rgCs_Out_1Profile)

    return rgCs_Out_Nested


def createArrayedGeometry(rgCrv_Path, rgObjs_ToArray, plane_Proj, fTaper_Start_Deg, fTaper_End_Deg, bTaperChangePerCrvParam, bAtGrevilles, bAtKnots, bAtEqualDivisions, iDivisionCt=None, bDebug=False):
    """
    """

    xform_Proj = rg.Transform.PlanarProjection(plane_Proj)

    rgObjs1_Arrayed_PerSection = []

    nc2_Path = rgCrv_Path.ToNurbsCurve()
    if nc2_Path is None:
        print("NurbsCurve could not be calculated from curve.")
        return
    #sc.doc.Objects.AddCurve(nc2_Path)
        
    ts = []

    if bAtGrevilles:
        rc = nc2_Path.GrevilleParameters()
        if rc: ts.extend(rc)
        if nc2_Path.IsClosed:
            ts.pop()

    if bAtEqualDivisions:
        rc = nc2_Path.DivideByCount(
                segmentCount=iDivisionCt,
                includeEnds=True)
        if rc: ts.extend(rc)
        if nc2_Path.IsClosed:
            ts.append(nc2_Path.Domain.T1)
    
    if bAtKnots:
        rc = set(nc2_Path.Knots)
        if rc: ts.extend(rc)
        
    if ts is None:
        print("No parameters were obtained.")
        return

    ts = sorted(set(ts)) # Remove duplicates and sort.

    # Remove overlaps for closed (including periodic) curves.
    if nc2_Path.IsClosed:
        if bDebug:
            print(nc2_Path.Domain)
            print(ts)
        ts_WIP = []
        for t in ts:
            if t >= nc2_Path.Domain.T0 and t < nc2_Path.Domain.T1:
                ts_WIP.append(t)
        ts = ts_WIP
        if bDebug: print(ts)

    if fTaper_Start_Deg == fTaper_End_Deg:
        angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
    elif not bTaperChangePerCrvParam:
        length_Full = rgCrv_Path.GetLength()

    rgCrv_Path_Flattened = nc2_Path.Duplicate()
    rgCrv_Path_Flattened.Transform(xform_Proj)

    #sc.doc.Objects.AddCurve(rgCrv_Path_Flattened); return
    rgObjs1_Arrayed_PerSection = []
    for iT, t in enumerate(ts):
        bSuccess, frame = rgCrv_Path_Flattened.PerpendicularFrameAt(t=t)
        if not bSuccess:
            print("Perpendicular frame could not be calculated.")
            continue
                
        angle_StraightenFrame_Rad = rg.Vector3d.VectorAngle(frame.YAxis, plane_Proj.ZAxis, frame)

        if fTaper_Start_Deg != fTaper_End_Deg:
            if t == rgCrv_Path.Domain.T0:
                angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
            elif t == rgCrv_Path.Domain.T1:
                angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_End_Deg)
            #elif len(ts) >= 4 and iT == 1:
            #    # Same angle so that the brep starts at a tangent.
            #    angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
            #elif len(ts) >= 4 and iT == len(ts) - 2:
            #    # Same angle so that the brep end at a tangent.
            #    angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_End_Deg)

            else:
                if bTaperChangePerCrvParam:
                    t_Normalized = rgCrv_Path.Domain.NormalizedParameterAt(t)
                    angle_Rad = Rhino.RhinoMath.ToRadians(
                            fTaper_Start_Deg * (1.0 - t_Normalized) +
                            fTaper_End_Deg * t_Normalized)
                else:
                    length_to_t = rgCrv_Path.GetLength(
                            subdomain=rg.Interval(rgCrv_Path.Domain.T0, t))
                    angle_Rad = Rhino.RhinoMath.ToRadians(
                            fTaper_Start_Deg * (1.0 - length_to_t/length_Full) +
                            fTaper_End_Deg * length_to_t/length_Full)


            if bDebug: print(angle_Rad)
        frame.Rotate(angle=angle_StraightenFrame_Rad+angle_Rad, axis=frame.ZAxis)

        # Debug frame orientation.
        #                sc.doc.Objects.AddPoint(frame.PointAt(0.0,0.0,0.0))
        #                attr_Red = rd.ObjectAttributes()
        #                attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
        #                attr_Red.ObjectColor = Color.Red
        #                sc.doc.Objects.AddPoint(frame.PointAt(1.0,0.0,0.0), attr_Red)
        #                attr_Green = rd.ObjectAttributes()
        #                attr_Green.ColorSource = rd.ObjectColorSource.ColorFromObject
        #                attr_Green.ObjectColor = Color.Lime
        #                sc.doc.Objects.AddPoint(frame.PointAt(0.0,1.0,0.0), attr_Green)
                
        xform1 = rg.Transform.PlaneToPlane(plane_Proj, frame)
        xform2 = rg.Transform.Translation(
                nc2_Path.PointAt(t)-rgCrv_Path_Flattened.PointAt(t))
        xform3 = xform2 * xform1
        rgObjs1_Arrayed_PerSection.append([])

        for rgObj0_ToArray in rgObjs_ToArray:
            rgObj1_Arrayed_WIP = rgObj0_ToArray.Duplicate()
            rgObj1_Arrayed_WIP.Transform(xform3)
            rgObjs1_Arrayed_PerSection[-1].append(rgObj1_Arrayed_WIP)
        if not rgObjs1_Arrayed_PerSection[-1]:
            del rgObjs1_Arrayed_PerSection[-1]

    nc2_Path.Dispose()
    rgCrv_Path_Flattened.Dispose()

    return rgObjs1_Arrayed_PerSection


def _matchCrvEndDirs(nc_ToMod, nc_Ref):
    """
    nc_ToMod is modified.
    """
    
    if nc_ToMod.Points.Count < 4:
        print("_matchCrvEndDirs skipped due to only {} control points along surface side.".format(nc_ToMod.Points.Count))
        return False
    
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


def createBrep(iBrepMethod, iLoftType, fBrepTol, rgCrv_Path, rgNurbsCrv_TaperEnd_1Seg, rgLineCrvs_Arrayed):
    """
    """

    rgNurbsCrv1_PathSeg = rgCrv_Path.ToNurbsCurve()

    rgBs_Out = []
    
    if Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Loft2Crvs':
        rgBs_Out = rg.Brep.CreateFromLoft(
                curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg],
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=rg.LoftType.Straight,
                closed=False)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'LoftSectionLines':
        for L in rgLineCrvs_Arrayed:
            print(L.PointAtStart, L.PointAtEnd)
        print(rgLineCrvs_Arrayed[0].PointAtEnd.EpsilonEquals(rgLineCrvs_Arrayed[-1].PointAtEnd, epsilon=1e-12))
        rgBs_Out = rg.Brep.CreateFromLoft(
                curves=rgLineCrvs_Arrayed,
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=Enum.ToObject(rg.LoftType, iLoftType),
                closed=rgNurbsCrv1_PathSeg.IsClosed)
        print(rgBs_Out)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Sweep2A':
        rgBs_Out = rg.Brep.CreateFromSweep(
                rail1=rgNurbsCrv1_PathSeg,
                rail2=rgNurbsCrv_TaperEnd_1Seg,
                shapes=rgLineCrvs_Arrayed,
                closed=rgNurbsCrv1_PathSeg.IsClosed,
                tolerance=fBrepTol)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Sweep2B':
        rgSweep2 = rg.SweepTwoRail()
        #rgSweep2.AngleToleranceRadians
        rgSweep2.ClosedSweep = rgNurbsCrv1_PathSeg.IsClosed
        rgSweep2.MaintainHeight = False
        rgSweep2.SweepTolerance = fBrepTol
        rgBs_Out = rgSweep2.PerformSweep(
                rail1=rgNurbsCrv1_PathSeg,
                rail2=rgNurbsCrv_TaperEnd_1Seg,
                crossSections=rgLineCrvs_Arrayed)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Network':
        rgNurbsSrf, iError = rg.NurbsSurface.CreateNetworkSurface(
                curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg]+rgLineCrvs_Arrayed,
                continuity=1,
                edgeTolerance=fBrepTol,
                interiorTolerance=fBrepTol,
                angleTolerance=0.1*sc.doc.ModelAngleToleranceDegrees)
        if iError:
            print("CreateNetworkSurface error code: {}".format(iError))
        else:
            rgBs_Out = [rgNurbsSrf.ToBrep()]
            rgNurbsSrf.Dispose()

    rgNurbsCrv1_PathSeg.Dispose()

    return rgBs_Out


def main():

    while True:
        sc.escape_test()

        rc = getInput(bFirstGetObjects=True)
        if rc is None: return

        objrefs_Paths, bGenerateNew, bAcceptResults = rc
        if not objrefs_Paths: continue

        fDistance = Opts.values['fDistance']
        bProjDist = Opts.values['bProjDist']
        fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
        bVariableTaper = Opts.values['bVariableTaper']
        fTaper_End_Deg = Opts.values['fTaper_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bCPlane = Opts.values['bCPlane']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPathsAtG2PlusKnots = Opts.values['bSplitPathsAtG2PlusKnots']
        bAddTaperStartCrvs = Opts.values['bAddTaperStartCrvs']
        bAddTaperedLines = Opts.values['bAddTaperedLines']
        bAddTaperEndCrvs = Opts.values['bAddTaperEndCrvs']
        bAddBrep = Opts.values['bAddBrep']
        iBrepMethod = Opts.values['iBrepMethod']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        break


    rgLine_ToArray = None
    rgCrvs0_Path = []
    rgCs_Paths_PreparedForSampling_AllProfiles = []

    while True:

        while not (bAtGrevilles or bAtKnots or bAtEqualDivisions):
            print("No path point sampling is enabled.")
            sc.doc.Views.Redraw()
            rc = getInput(bFirstGetObjects=False)
            if rc is None:
                if rgLine_ToArray: rgLine_ToArray.Dispose()
                for c in rgCrvs0_Path: c.Dispose()
                for c in rgCs_Paths_PreparedForSampling_AllProfiles: c.Dispose()
                return

            objrefs_Paths, bGenerateNew, bAcceptResults = rc

            if bAcceptResults:
                if rgLine_ToArray: rgLine_ToArray.Dispose()
                for c in rgCrvs0_Path: c.Dispose()
                for c in rgCs_Paths_PreparedForSampling_AllProfiles: c.Dispose()
                return

            fDistance = Opts.values['fDistance']
            bProjDist = Opts.values['bProjDist']
            fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
            bVariableTaper = Opts.values['bVariableTaper']
            fTaper_End_Deg = Opts.values['fTaper_End_Deg']
            bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
            bAlignEndDirs = Opts.values['bAlignEndDirs']
            bCPlane = Opts.values['bCPlane']
            bAtGrevilles = Opts.values['bAtGrevilles']
            bAtKnots = Opts.values['bAtKnots']
            bAtEqualDivisions = Opts.values['bAtEqualDivisions']
            iDivisionCt = Opts.values['iDivisionCt']
            bSplitPathsAtG2PlusKnots = Opts.values['bSplitPathsAtG2PlusKnots']
            bAddTaperStartCrvs = Opts.values['bAddTaperStartCrvs']
            bAddTaperedLines = Opts.values['bAddTaperedLines']
            bAddTaperEndCrvs = Opts.values['bAddTaperEndCrvs']
            bAddBrep = Opts.values['bAddBrep']
            iBrepMethod = Opts.values['iBrepMethod']
            iLoftType = Opts.values['iLoftType']
            fBrepTol = Opts.values['fBrepTol']
            bEcho = Opts.values['bEcho']
            bDebug = Opts.values['bDebug']


        rgCrvs0_Path = []
        for objref_Path in objrefs_Paths:
            c = objref_Path.Curve()
            rgCrvs0_Path.append(c)

        bMakeDeformable = bAlignEndDirs or (fAngle_End_Deg and fAngle_End_Deg != fAngle_Start_Deg)

        Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")

        rgCs_Paths_PreparedForSampling_AllProfiles = _prepareCurves(
            rgCrvs_In=rgCrvs0_Path,
            bSplitPathsAtG2PlusKnots=bSplitPathsAtG2PlusKnots,
            bMakeDeformable=bMakeDeformable)
        if not rgCs_Paths_PreparedForSampling_AllProfiles:
            print("Path curve was not prepared.")
            return

        if bProjDist and not bVariableTaper:
            rgLine_ToArray = rg.Line(
                rg.Point3d(0.0, 0.0, 0.0),
                rg.Point3d(
                    0.0,
                    fDistance/math.cos(math.radians(fTaper_Start_Deg)),
                    0.0))
        else:
            rgLine_ToArray = rg.Line(
                rg.Point3d(0.0, 0.0, 0.0),
                rg.Point3d(0.0, fDistance, 0.0))

        if bCPlane:
            view_Active = sc.doc.Views.ActiveView
            plane_Proj = view_Active.ActiveViewport.ConstructionPlane()
        
            xform1 = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, plane_Proj)
            rgLine_ToArray.Transform(xform1)
        else:
            plane_Proj = rg.Plane.WorldXY

        rgLineCrv_ToArray = rg.LineCurve(rgLine_ToArray)

        gCrvs_TaperStart_All = []
        gLineCrvs1_Arrayed = []
        gCrvs_TaperEnd_All = []
        gBs_Out = []

        rgCrvs_PathSegs = None

        Rhino.RhinoApp.SetCommandPrompt("Creating geometry ...")

        for i_Path_Profile, rgCs_Path_Profile in enumerate(rgCs_Paths_PreparedForSampling_AllProfiles):

            rgBs_1Profile = []

            for i_Path_1Seg, rgCrv1_Path_1Seg in enumerate(rgCs_Path_Profile):
                rc = createArrayedGeometry(
                    rgCrv_Path=rgCrv1_Path_1Seg,
                    rgObjs_ToArray=[rgLineCrv_ToArray],
                    plane_Proj=plane_Proj,
                    fTaper_Start_Deg=fTaper_Start_Deg,
                    fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                    bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                    bAtGrevilles=bAtGrevilles,
                    bAtKnots=bAtKnots,
                    bAtEqualDivisions=bAtEqualDivisions,
                    iDivisionCt=iDivisionCt,
                    bDebug=bDebug)
                if rc is None: return
                print(rc)
                # Flatten list.
                rgLineCrvs_Arrayed_1PathSeg = [rgLs[0] for rgLs in rc]
                print(rgLineCrvs_Arrayed_1PathSeg)

                if bAddTaperedLines:
                    pts_Start = []
                    pts_Ends = []
                    for i, rgLc in enumerate(rgLineCrvs_Arrayed_1PathSeg):
                        if i_Path_1Seg > 0:
                            if rgLc.PointAtStart.EpsilonEquals(
                                rgLc_LastPathSeg.PointAtStart,
                                epsilon=sc.doc.ModelAbsoluteTolerance
                            ):
                                if rgLc.PointAtEnd.EpsilonEquals(
                                    rgLc_LastPathSeg.PointAtEnd,
                                    epsilon=sc.doc.ModelAbsoluteTolerance
                                ):
                                    print("Skipped duplicate line.")
                                    continue
                        gL = sc.doc.Objects.AddCurve(rgLc)
                        if gL != gL.Empty:
                            gLineCrvs1_Arrayed.append(gL)
                            #rgDot_ = rg.TextDot(str(i), rgObj1_Arrayed_1Seg.PointAtStart)
                            #sc.doc.Objects.AddTextDot(rgDot_)
                    rgLc_LastPathSeg = rgLc
                else:
                    # bAddTaperedLines == False.
                    if bAtKnots or bAtEqualDivisions:
                        if bAtKnots and iBrepMethod == 0:
                            s  = "AtKnots"
                            s += " only affects added arrayed lines"
                            s += " when BrepMethod == {},".format(Opts.listValues['iBrepMethod'][iBrepMethod])
                            s += " but AddArrayed option is disabled."
                            print(s)
                        if bAtEqualDivisions and iBrepMethod == 0:
                            s  = "AtEqualDivisions"
                            s += " only affects added arrayed lines"
                            s += " when BrepMethod == {},".format(Opts.listValues['iBrepMethod'][iBrepMethod])
                            s += " but AddArrayed option is disabled."
                            print(s)

                if bAddTaperStartCrvs or bAddTaperEndCrvs or bAddBrep:

                    if not bAtGrevilles or bAtKnots or bAtEqualDivisions:
                        rc = createArrayedGeometry(
                            rgCrv_Path=rgCrv1_Path_1Seg,
                            rgObjs_ToArray=[rgLineCrv_ToArray],
                            plane_Proj=plane_Proj,
                            fTaper_Start_Deg=fTaper_Start_Deg,
                            fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                            bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                            bAtGrevilles=True,
                            bAtKnots=False,
                            bAtEqualDivisions=False,
                            iDivisionCt=0,
                            bDebug=bDebug)
                        if rc is None: continue
                        # Flatten list.
                        rgLineCrvs_Arrayed_1PathSeg_GrevsOnly = [L[0] for L in rc]
                    else:
                        rgLineCrvs_Arrayed_1PathSeg_GrevsOnly = [L.Duplicate() for L in rgLineCrvs_Arrayed_1PathSeg]
                
                    pts_EndOf_LineCrvs_Arrayed = []
                    for rgLineCrv in rgLineCrvs_Arrayed_1PathSeg_GrevsOnly:
                        pts_EndOf_LineCrvs_Arrayed.append(rgLineCrv.PointAtEnd)

                        nc_OppOfPath = rgCrv1_Path_1Seg.ToNurbsCurve()
                        nc_OppOfPath.SetGrevillePoints(pts_EndOf_LineCrvs_Arrayed)

                    if bAlignEndDirs:
                        _matchCrvEndDirs(nc_OppOfPath, rgCrv1_Path_1Seg)

                    if bAddTaperStartCrvs:
                        gCrv_TaperStart_1Seg = sc.doc.Objects.AddCurve(rgCrv1_Path_1Seg)
                        if gCrv_TaperStart_1Seg != gCrv_TaperStart_1Seg.Empty:
                            gCrvs_TaperStart_All.append(gCrv_TaperStart_1Seg)

                    if bAddTaperEndCrvs:
                        gCrv_TaperEnd_1Seg = sc.doc.Objects.AddCurve(nc_OppOfPath)
                        if gCrv_TaperEnd_1Seg != gCrv_TaperEnd_1Seg.Empty:
                            gCrvs_TaperEnd_All.append(gCrv_TaperEnd_1Seg)
                    if bAddBrep:
                        rc = createBrep(
                                iBrepMethod=iBrepMethod,
                                iLoftType=iLoftType,
                                fBrepTol=fBrepTol,
                                rgCrv_Path=rgCrv1_Path_1Seg,
                                rgNurbsCrv_TaperEnd_1Seg=nc_OppOfPath,
                                rgLineCrvs_Arrayed=rgLineCrvs_Arrayed_1PathSeg_GrevsOnly)
                        if rc is None:
                            print("Cannot create brep(s).  Check input.")
                        else:
                            rgBs_1Profile.extend(rc)
                        rgCrv1_Path_1Seg.Dispose()
                
                    for c in rgLineCrvs_Arrayed_1PathSeg_GrevsOnly:
                        c.Dispose()


            if rgBs_1Profile:
                rgBreps1_Joined = rg.Brep.JoinBreps(
                        rgBs_1Profile,
                        tolerance=0.5*sc.doc.ModelAbsoluteTolerance)
                for b in rgBs_1Profile: b.Dispose()

                if rgBreps1_Joined:
                    #gBs_Out = []
                    iCt_Faces_Added = 0
                    for rgBrep1_Joined in rgBreps1_Joined:
                        gBrep1 = sc.doc.Objects.AddBrep(rgBrep1_Joined)
                        if gBrep1 != gBrep1.Empty:
                            gBs_Out.append(gBrep1)
                            rdB_Out = sc.doc.Objects.FindId(gBrep1)
                            iCt_Faces_Added += rdB_Out.BrepGeometry.Faces.Count
                    if bEcho:
                        print("{} breps with {} faces created.".format(
                            len(gBs_Out), iCt_Faces_Added))

        sc.doc.Views.Redraw()

        rc = getInput(bFirstGetObjects=False)
        if rc is None:
            for g_ in gCrvs_TaperStart_All:
                sc.doc.Objects.Delete(g_, True)
            for g_ in gLineCrvs1_Arrayed:
                sc.doc.Objects.Delete(g_, True)
            for g_ in gCrvs_TaperEnd_All:
                sc.doc.Objects.Delete(g_, True)
            for g_ in gBs_Out:
                sc.doc.Objects.Delete(g_, True)
            break

        objrefs_Paths, bGenerateNew, bAcceptResults = rc

        if bAcceptResults:
            break

        if not bGenerateNew:
            continue

        fDistance = Opts.values['fDistance']
        bProjDist = Opts.values['bProjDist']
        fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
        bVariableTaper = Opts.values['bVariableTaper']
        fTaper_End_Deg = Opts.values['fTaper_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bCPlane = Opts.values['bCPlane']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPathsAtG2PlusKnots = Opts.values['bSplitPathsAtG2PlusKnots']
        bAddTaperStartCrvs = Opts.values['bAddTaperStartCrvs']
        bAddTaperedLines = Opts.values['bAddTaperedLines']
        bAddTaperEndCrvs = Opts.values['bAddTaperEndCrvs']
        bAddBrep = Opts.values['bAddBrep']
        iBrepMethod = Opts.values['iBrepMethod']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


        for g_ in gCrvs_TaperStart_All:
            sc.doc.Objects.Delete(g_, True)
        for g_ in gLineCrvs1_Arrayed:
            sc.doc.Objects.Delete(g_, True)
        for g_ in gCrvs_TaperEnd_All:
            sc.doc.Objects.Delete(g_, True)
        for g_ in gBs_Out:
            sc.doc.Objects.Delete(g_, True)

    for c in rgCrvs0_Path: c.Dispose()
    for cs in rgCs_Paths_PreparedForSampling_AllProfiles:
        for c in cs:
            c.Dispose()
    rgLineCrv_ToArray.Dispose()
    for c in rgLineCrvs_Arrayed_1PathSeg:
        c.Dispose()


if __name__ == '__main__': main()