"""
200619: Created.
200629: Bug fix.
200729: Import-related update.
200810: Added a couple of functions.
210114: Removed raise for when createPoint3dOnInterior doesn't find a point.
210209: Now won't raise an error when 2 curves do not close in splitFace.
210407: Added more support for new-in-V7's BrepFace.PerFaceColor.
210412: Added more times fEdgeLen_Min is used in splitFace.
210503: Now, curves on SENW are ignored for retrim when underlying surface is not being replaced.
210929: is3dPointOnFace now converts non-NurbsSurfaces into NurbsSurfaces on ClosestPoint fails.
220317: Updated to use a new overload.  Import-related update.
230830: Disabled a routine in splitFace for testing purposes.
230901: Added code for debugging.
250324: Added routines that adjust curves to meet with natural corners and/or sides of surfaces.
        Added check for surfaces that create bad breps.
        Changed 2.0 to 1.8 to match Rhino (8)'s tolerances for joining curves.
        Modified what is used for curve joining tolerance.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc

from System import Guid
from System import Random

import xCurve
import xSurface


def createPoint3dOnInterior(rgFace, fMinDistFromBorder=None, bDebug=False):
    """
    Parameters:
        rgFace: BrepFace
        fMinDistFromBorder: float or (None to loop through various)
    Returns:
        Point3d on success, None on failure
    """
    
    fMinDistFromBorder_In = fMinDistFromBorder
    
    rgB_1F = rgFace.DuplicateFace(duplicateMeshes=False)
    rgB_1F.Faces.ShrinkFaces()
    rgB_1F.Compact() # Is this cecessary to remove duplicate surfaces?
    rgSrf_Shrunk = rgB_1F.Surfaces[0]
    
    domU_Shrunk = rgSrf_Shrunk.Domain(0)
    domV_Shrunk = rgSrf_Shrunk.Domain(1)

    #areaMassProp = rg.AreaMassProperties.Compute(rgFace.DuplicateFace(True)) # Otherwise, Compute uses underlying Surface for BrepFace.
    #if areaMassProp is None:
    #    print "Face[{}]'s AreaMassProperties cannot be calculated.".format(
    #        rgFace.FaceIndex)
    #    u = domU_Shrunk.Mid
    #    v = domV_Shrunk.Mid
    #else:
    #    ptCentrdW = areaMassProp.Centroid
    #    bSuccess, u, v = rgFace.ClosestPoint(ptCentrdW)

    u_Shrunk = domU_Shrunk.Mid
    v_Shrunk = domV_Shrunk.Mid

    # To3dCurve will include curves of seams.  This is advantageous in this circumstance.
    cs_All_loops_for_split = [loop.To3dCurve() for loop in rgFace.Loops]

    rand = Random()

    if fMinDistFromBorder_In is None:
        rangeTol = (
            10.0*sc.doc.ModelAbsoluteTolerance,
            sc.doc.ModelAbsoluteTolerance,
            0.5*sc.doc.ModelAbsoluteTolerance
            )
    else:
        rangeTol = fMinDistFromBorder_In,

    for fMinDistFromBorder in rangeTol:
        for i in xrange(1000):
            pt = rgSrf_Shrunk.PointAt(u_Shrunk, v_Shrunk)
            
            b, u, v = rgFace.ClosestPoint(pt)
            if not b: raise ValueError("ClosestPoint failed in createPointOnFace.")
            
            ptFaceRel = rgFace.IsPointOnFace(u, v)
            #sc.doc.Objects.AddPoint(pt); sc.doc.Views.Redraw()

            if bDebug: sEval = "ptFaceRel"; print sEval+':',eval(sEval)

            if ptFaceRel == rg.PointFaceRelation.Interior:
                # If point is not at least fMinDistFromBorder from border (all edges),
                # continue searching.
                for c in cs_All_loops_for_split:
                    b, t = c.ClosestPoint(pt)
                    if not b: raise ValueError("ClosestPoint failed in createPointOnFace.")
    
                    dist = pt.DistanceTo(c.PointAt(t))
                    if dist < fMinDistFromBorder:
                        break # to get another u and v.
                else: # Good point
                    for c in cs_All_loops_for_split: c.Dispose()
                    rgB_1F.Dispose()
                    if bDebug: sEval = "i"; print sEval+':',eval(sEval)
                    return pt
    
            # Get new parameters for point.
            u_Shrunk = rand.NextDouble() * domU_Shrunk.Length + domU_Shrunk.Min
            v_Shrunk = rand.NextDouble() * domV_Shrunk.Length + domV_Shrunk.Min

    from System.Drawing import Color
    attr = rd.ObjectAttributes()
    attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr.ObjectColor = Color.Red
    sc.doc.Objects.AddBrep(rgFace.DuplicateFace(duplicateMeshes=True), attributes=attr)

    print "Failed to find an interior point on face." \
        "  Its monoface has been added to the document."

    #sc.doc.Views.Redraw()
    #raise ValueError(
    #    "Failed to find an interior point on face[{}].".format(
    #        rgFace.FaceIndex))


def is3dPointOnFace(rgFace, pt, fTolerance=None):
    """
    Returns:
        False when the point is at least fTolerance from the face's underlying surface.
        Otherwise, rg.PointFaceRelation.Interior or Boundary.
    """
    
    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance
    
    b, u, v = rgFace.ClosestPoint(pt) # ClosestPoint is a Surface method and ignores edges of BrepFaces.
    if not b:
        s = "ClosestPoint on {} failed ".format(rgFace.UnderlyingSurface().GetType().Name)
        s += "at {:.15g},{:.15g},{:.15g}.".format(pt.X, pt.Y, pt.Z)
        if isinstance(rgFace.UnderlyingSurface(), rg.NurbsSurface):
            print s
            return

        ns = rgFace.UnderlyingSurface().ToNurbsSurface()
        b, u, v = ns.ClosestPoint(pt)
        ns.Dispose()
        if not b:
            s += "  It also failed for the NurbsSurface equivalent."
            print s
            return
        s += "  But NurbsSurface equivalent passed."
        print s

    dist = rgFace.PointAt(u,v).DistanceTo(pt)
    if dist > fTolerance:
        # Point isn't even on underlying surface.
        return False

    try:
        ptFaceRel = rgFace.IsPointOnFace(u, v, fTolerance) # Overload added to RhinoCommon 7.0.
    except:
        ptFaceRel = isParameterPointOnFace_beforeV7(
            rgFace, u, v, fTolerance)

    if ptFaceRel == rg.PointFaceRelation.Interior:
        pass
    elif ptFaceRel == rg.PointFaceRelation.Exterior:
        pass
    elif ptFaceRel == rg.PointFaceRelation.Boundary:
        pass
    else:
        pass

    return ptFaceRel


def isParameterPointOnFace_beforeV7(rgFace, u, v, f3dTolerance=None):
    """
    Hard-coded (3D) tolerance of IsPointOnFace up through at least V6.28 (7/17/2020)
    is Rhino.RhinoMath.SqrtEpsilon.
    Request to add tolerance parameter:
    https://discourse.mcneel.com/t/wish-add-tolerance-parameter-to-ispointonface/106801
    Until tolerance parameter is added to IsPointOnFace,
    check False values against distance from closest edge.
    """
    
    if f3dTolerance is None:
        f3dTolerance = sc.doc.ModelAbsoluteTolerance
    
    ptAt = rgFace.PointAt(u, v)
    # Note that this Point3d will already be pulled to the face.
    # For true 3D point check, use is3dPointOnFace instead.
    
    ptFaceRel = rgFace.IsPointOnFace(u, v)
    
    if ptFaceRel == rg.PointFaceRelation.Boundary:
        return rg.PointFaceRelation.Boundary
    
    for iE in rgFace.AdjacentEdges():
        bSuccess, t = rgFace.Brep.Edges[iE].ClosestPoint(
            testPoint=ptAt,
            maximumDistance=f3dTolerance)
        if not bSuccess:
            continue
        pt_onEdge = rgFace.Brep.Edges[iE].PointAt(t)
        dist_toCrv = ptAt.DistanceTo(other=pt_onEdge)
        if dist_toCrv <= f3dTolerance:
            return rg.PointFaceRelation.Boundary
    else:
        # No distances smaller than f3dTolerance found.
        return ptFaceRel # Exterior or boundary.


def separate_short_curves_from_others_in_loop(cs_ofLoop, fTol_MinDistBetweenOpenEnds=None):
    """
    """

    if fTol_MinDistBetweenOpenEnds is None:
        fTol_MinDistBetweenOpenEnds = 1.8 * sc.doc.ModelAbsoluteTolerance

    cs_WIP = []
    cs_Shorts = []
    for c in cs_ofLoop:
        
        # 250324: Commented out and replaced by below this comment block.
        #fLength = c.GetLength()
        #if fLength > max(fSplitTol, fEdgeLen_Min):
        #    cs_WIP.append(c)
        #else:
        #    c.Dispose()
        
        if c.IsClosed:
            if len(cs_ofLoop) > 1:
                print("Closed curve found in loop of {}.".format(len(cs_ofLoop)))
                return

            fLength = c.GetLength()
            if fLength <= (2.0 * fTol_MinDistBetweenOpenEnds):
                import rhinoscriptsyntax as rs
                point = c.PointAt(c.Domain.Mid)
                rs.AddTextDot(
                    text="Tiny trim loop?",
                    point=point)
                break
        
            cs_WIP.append(c)
            continue
        
        fDist_BtwnEndPts = c.PointAtStart.DistanceTo(c.PointAtEnd)
        if fDist_BtwnEndPts <= fTol_MinDistBetweenOpenEnds:
            fLength = c.GetLength()
            if fLength >= (2.0 * fDist_BtwnEndPts): # 2.0 is a guess. See print below.
                import rhinoscriptsyntax as rs
                #_rgB = rgFace.DuplicateFace(duplicateMeshes=False)
                #_rgB.Faces.ShrinkFaces()
                #_rgF = _rgB.Faces[0]
                #point = _rgF.PointAt(u=_rgF.Domain(0).Mid, v=_rgF.Domain(1).Mid)
                point = c.PointAt(c.Domain.Mid)
                rs.AddTextDot(
                    text="Large loop in trim?",
                    point=point)
                #_rgB.Dispose()
                print("End points are close but length is at least double the distance."
                    "This face will not be retrimmed."
                    "Dot added.")
                return
            #sEval="c.IsDocumentControlled"; print(sEval,'=',eval(sEval))
            #c.Dispose()
            cs_Shorts.append(c)
            continue

        cs_WIP.append(c)

    return cs_WIP, cs_Shorts


def rebuild_short_contiguous_crvs(cs_Short, fTol_MinDistBetweenOpenEnds=None):
    """
    """

    if fTol_MinDistBetweenOpenEnds is None:
        fTol_MinDistBetweenOpenEnds = 1.8 * sc.doc.ModelAbsoluteTolerance

    joinTolerance = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)

    cs_joined = rg.Curve.JoinCurves(cs_Short, joinTolerance=joinTolerance)

    cs_Out = []

    for cs in cs_joined:
        if not isinstance(cs, (rg.PolyCurve, rg.PolylineCurve)):
            continue

        fDist = cs.PointAtStart.DistanceTo(cs.PointAtEnd)

        if fDist <= fTol_MinDistBetweenOpenEnds:
            continue

        rebuilt = cs.Rebuild(
            pointCount=4,
            degree=3,
            preserveTangents=False)

        rvs = rg.Curve.GetDistancesBetweenCurves(cs, rebuilt, tolerance=0.1*joinTolerance)
        if not rvs[0]:
            continue

        if rvs[1] > 0.1 * sc.doc.ModelAbsoluteTolerance:
            continue

        cs_Out.append(rebuilt)

    return cs_Out


def _getNaturalCornerPts(srf):
    return (
        srf.PointAt(srf.Domain(0).Min, srf.Domain(1).Min),
        srf.PointAt(srf.Domain(0).Min, srf.Domain(1).Max),
        srf.PointAt(srf.Domain(0).Max, srf.Domain(1).Min),
        srf.PointAt(srf.Domain(0).Max, srf.Domain(1).Max)
        )


def _getCrvsOfIsoCrvs(srf):
    return (
        srf.IsoCurve(direction=1, constantParameter=srf.Domain(0).T0),
        srf.IsoCurve(direction=0, constantParameter=srf.Domain(1).T0),
        srf.IsoCurve(direction=1, constantParameter=srf.Domain(0).T1),
        srf.IsoCurve(direction=0, constantParameter=srf.Domain(1).T1)
        )


def tryToExtendCrvToNearPtOfSet(crv_In, side, pts_Target, tolerance):
    dists = []
    if side == rg.CurveEnd.Start:
        pt_c = crv_In.PointAtStart
    elif side == rg.CurveEnd.End:
        pt_c = crv_In.PointAtEnd
    else:
        raise Exception("{} passed as the side.".format(side))

    for pt_Target in pts_Target:
        dist = pt_c.DistanceTo(pt_Target)
        dists.append(dist)

    if min(dists) > tolerance:
        return

    idx_Winner = dists.index(min(dists))

    return crv_In.Extend(
        side=side,
        style=rg.CurveExtensionStyle.Smooth,
        endPoint=pts_Target[idx_Winner])


def tryToSetCrvEndToPtOfSet(crv_In, side, pts_Target, tolerance):
    dists = []
    if side == rg.CurveEnd.Start:
        pt_c = crv_In.PointAtStart
        setPoint = rg.Curve.SetStartPoint
    elif side == rg.CurveEnd.End:
        pt_c = crv_In.PointAtEnd
        setPoint = rg.Curve.SetEndPoint
    else:
        raise Exception("{} passed as the side.".format(side))

    for pt_Target in pts_Target:
        dist = pt_c.DistanceTo(pt_Target)
        dists.append(dist)

    if min(dists) > tolerance:
        return

    idx_Winner = dists.index(min(dists))

    crv_Out = crv_In.DuplicateCurve()

    bSuccess = setPoint(crv_Out, point=pts_Target[idx_Winner])

    if bSuccess:
        return crv_Out

    crv_Out.Dispose()


def tryToExtendCrvToNearCrvOfSet(crv_In, side, crvs_Target, tolerance):
    dists = []
    if side == rg.CurveEnd.Start:
        pt_c = crv_In.PointAtStart
    elif side == rg.CurveEnd.End:
        pt_c = crv_In.PointAtEnd
    else:
        raise Exception("{} passed as the side.".format(side))

    for crv_Target in crvs_Target:
        bSuccess, t = crv_Target.ClosestPoint(pt_c)
        dist = pt_c.DistanceTo(crv_Target.PointAt(t))
        dists.append(dist)

    if min(dists) > tolerance:
        return

    idx_Winner = dists.index(min(dists))

    return crv_In.Extend(
        side=side,
        style=rg.CurveExtensionStyle.Smooth,
        geometry=[crvs_Target[idx_Winner]])


def is_pt_near_a_pt_in_set(pt_toTest, pts_Ref, tolerance):
    for pt_Ref in pts_Ref:
        dist = pt_toTest.DistanceTo(pt_Ref)
        if dist <= tolerance:
            return True
    return False


def is_pt_near_a_crv_of_set(pt_toTest, crvs_Ref, tolerance):
    for crv_Ref in crvs_Ref:
        bSuccess, t = crv_Ref.ClosestPoint(pt_toTest)
        dist = pt_toTest.DistanceTo(crv_Ref.PointAt(t))
        if dist <= tolerance:
            return True
    return False


def tryToAdjustCrvEndsToSENW(crv_In, rgSrf_toSplit, tolerance):
    zeroTol = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)

    pts_Corner = _getNaturalCornerPts(rgSrf_toSplit)

    bStartIsAtCorner = is_pt_near_a_pt_in_set(
        crv_In.PointAtStart, pts_Corner, tolerance=zeroTol)

    bEndIsAtCorner = is_pt_near_a_pt_in_set(
        crv_In.PointAtEnd, pts_Corner, tolerance=zeroTol)

    if bStartIsAtCorner and bEndIsAtCorner:
        return

    bProcessStart = not bStartIsAtCorner and is_pt_near_a_pt_in_set(
        crv_In.PointAtStart, pts_Corner, tolerance=tolerance)
    bProcessEnd = not bEndIsAtCorner and is_pt_near_a_pt_in_set(
        crv_In.PointAtEnd, pts_Corner, tolerance=tolerance)

    crv_Out = None
    crv_WIP = crv_In

    if bProcessStart:
        rv = tryToExtendCrvToNearPtOfSet(
            crv_WIP,
            side=rg.CurveEnd.Start,
            pts_Target=pts_Corner,
            tolerance=tolerance)
        if rv:
            crv_Out = crv_WIP = rv
        else:
            rv = tryToSetCrvEndToPtOfSet(
                crv_WIP,
                side=rg.CurveEnd.Start,
                pts_Target=pts_Corner,
                tolerance=tolerance)
            if rv:
                crv_Out = crv_WIP = rv

    if bProcessEnd:
        rv = tryToExtendCrvToNearPtOfSet(
            crv_WIP,
            side=rg.CurveEnd.End,
            pts_Target=pts_Corner,
            tolerance=tolerance)
        if rv:
            crv_Out = crv_WIP = rv
        else:
            rv = tryToSetCrvEndToPtOfSet(
                crv_WIP,
                side=rg.CurveEnd.End,
                pts_Target=pts_Corner,
                tolerance=tolerance)
            if rv:
                crv_Out = crv_WIP = rv

    bStartIsAtCorner = is_pt_near_a_pt_in_set(
        crv_In.PointAtStart, pts_Corner, tolerance=zeroTol)

    bEndIsAtCorner = is_pt_near_a_pt_in_set(
        crv_In.PointAtEnd, pts_Corner, tolerance=zeroTol)

    if bStartIsAtCorner and bEndIsAtCorner:
        return crv_Out

    isoCrvs = _getCrvsOfIsoCrvs(rgSrf_toSplit)

    bStartIsOnSide = is_pt_near_a_crv_of_set(
        crv_In.PointAtStart, isoCrvs, tolerance=zeroTol)

    bEndIsAtCorner = is_pt_near_a_crv_of_set(
        crv_In.PointAtEnd, isoCrvs, tolerance=zeroTol)

    if bStartIsOnSide and bEndIsAtCorner:
        return crv_Out

    bProcessStart = not bStartIsOnSide and is_pt_near_a_crv_of_set(
        crv_In.PointAtStart, isoCrvs, tolerance=tolerance)
    bProcessEnd = not bEndIsAtCorner and is_pt_near_a_crv_of_set(
        crv_In.PointAtEnd, isoCrvs, tolerance=tolerance)

    if (not bProcessStart) and (not bProcessEnd):
        return crv_Out

    if bProcessStart:
        rv = tryToExtendCrvToNearCrvOfSet(
            crv_WIP,
            side=rg.CurveEnd.Start,
            crvs_Target=isoCrvs,
            tolerance=tolerance)
        if rv:
            crv_Out = crv_WIP = rv

    if bProcessEnd:
        rv = tryToExtendCrvToNearCrvOfSet(
            crv_WIP,
            side=rg.CurveEnd.End,
            crvs_Target=isoCrvs,
            tolerance=tolerance)

        if rv:
            crv_Out = crv_WIP = rv

    return crv_Out


def _isTrim_SENW(rgTrim):
    return (rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.South or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.East or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.North or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.West)


def splitFace(rgFace, **kwargs):
    """
    Parameters:
        rgFace
        rgSrf_Replacement: If None, use UnderlyingSurface of BrepFace.
        rgCrvs_Splitters_Replacement: If None, use 3D curves of rgFace loops.
        fSplitTol
        bTrimToSegs
        fEdgeLen_Min
        bOutputTrimmingCrvs
        bEcho
        bDebug
    Returns on success:
        rg.Brep
    Returns on fail:
        None
    """


    def setOpt(key, value): return kwargs[key] if key in kwargs else value

    rgSrf_Replacement = setOpt('rgSrf_Replacement', None)
    rgCrvs_Splitters_Replacement = setOpt('rgCrvs_Splitters_Replacement', None)
    fSplitTol = setOpt('fSplitTol', sc.doc.ModelAbsoluteTolerance)
    bTrimToSegs = setOpt('bTrimToSegs', True)
    fEdgeLen_Min = setOpt('fEdgeLen_Min', 1.8*sc.doc.ModelAbsoluteTolerance)
    bOutputTrimmingCrvs = setOpt('bOutputTrimmingCrvs', False)
    bEcho = setOpt('bEcho', True)
    bDebug = setOpt('bDebug', True)


    if bDebug: print "Function splitFace():"
    
    rgF_In = rgFace

    rgSrf_toSplit = rgF_In.UnderlyingSurface() if rgSrf_Replacement is None else rgSrf_Replacement

    if not xSurface.is_srf_valid_for_brep(rgSrf_toSplit):
        if isinstance(rgSrf_toSplit, rg.NurbsSurface):
            return
        rv = xSurface.repair_srf_invalid_for_brep(rgSrf_toSplit)
        if rv is None: return
        rgSrf_toSplit = rv

    # Reject surfaces that contain natural edges < 2 * fSplitTol.
    # 230830: Disabled this to enable surfaces with poles and to find examples where this routine is useful.
    #    for iDir, t in zip(
    #        (0, 0, 1, 1),
    #        (rgSrf_toSplit.Domain(1).T0, rgSrf_toSplit.Domain(1).T1,
    #            rgSrf_toSplit.Domain(0).T0, rgSrf_toSplit.Domain(1).T1,),
    #    ):
    #        isoCrv = rgSrf_toSplit.IsoCurve(iDir, t)
    #        fLength = isoCrv.GetLength()
    #        isoCrv.Dispose()
    #        if 0.0 < fLength <= 1.8 * sc.doc.ModelAbsoluteTolerance:
    #            return

    # 250324: Commented out.
    #fEdgeTols = []
    #for idxE in rgF_In.AdjacentEdges():
    #    fEdgeDev = rgF_In.Brep.Edges[idxE].Tolerance
    #    if fEdgeDev is not None:
    #        fEdgeTols.append(fEdgeDev)
    #fMaxEdgeTol = max(fEdgeTols)



    # Although splits may be successful if trims include seams,
    # they are avoided for cleaner input.  Otherwise, this could be used:
    # cs_All_loops_for_split = [loop.To3dCurve() for loop in rgFace_In.Loops]

    isoCrvs = _getCrvsOfIsoCrvs(rgSrf_toSplit)

    cs_All_loops_for_split = []

    for iL in range(rgF_In.Loops.Count):
        rgL = rgF_In.Loops[iL]
        cs_ofLoop = []

        bNaturalEdgesSkipped = False
        bNearNaturalEdgeSkipped = False

        for iT in range(rgL.Trims.Count):
            rgT = rgL.Trims[iT]
            if rgT.TrimType == rg.BrepTrimType.Seam: continue
            if rgT.TrimType == rg.BrepTrimType.Singular: continue
            if rgT.TrimType not in (
                rg.BrepTrimType.Boundary, rg.BrepTrimType.Mated
            ):
                print rgT.TrimType

            if rgSrf_Replacement is None:
                # Ignore curves on SENW.
                if _isTrim_SENW(rgT):
                    if bDebug: print "SENW curve ignored."
                    bNaturalEdgesSkipped = True
                    continue

            crv_EdgeDup = rgT.Edge.DuplicateCurve()

            for isoCrv in isoCrvs:
                rvs = rg.Curve.GetDistancesBetweenCurves(crv_EdgeDup, isoCrv, tolerance=0.1*fSplitTol)
                if not rvs[0]:
                    continue
                if rvs[1] <= fSplitTol:
                    bNearNaturalEdgeSkipped = True
                    break
            else:
                cs_ofLoop.append(crv_EdgeDup)

        #for c in cs_ofLoop: sc.doc.Objects.AddCurve(c)
        #sc.doc.Views.Redraw(); 1/0
        
        rvs = separate_short_curves_from_others_in_loop(
            cs_ofLoop,
            fTol_MinDistBetweenOpenEnds=fEdgeLen_Min)
        if rvs is None:
            return
        
        cs_toJoin, cs_Short = rvs
        
        if cs_Short:
            rvs = rebuild_short_contiguous_crvs(
                cs_Short, 
                fTol_MinDistBetweenOpenEnds=fEdgeLen_Min)
            if rvs:
                cs_toJoin.extend(rvs)

        if len(cs_toJoin) == 0:
            continue

        #cs_toJoin = cs_WIP[:]
        
        #250324: Replaced this line with below: fJoinTol = max(2.0*fMaxEdgeTol, 2.0*fSplitTol, fEdgeLen_Min)
        fJoinTol = fEdgeLen_Min
        
        rgCrvs_Joined = rg.Curve.JoinCurves(
            cs_toJoin,
            joinTolerance=fJoinTol)

        #for c in rgCrvs_Joined: sc.doc.Objects.AddCurve(c)
        #sc.doc.Views.Redraw(); 1/0


        for c in rgCrvs_Joined:
            # Open loop is OK if natural edges were skipped.
            if not c.IsClosed:
                if not bNaturalEdgesSkipped and not bNearNaturalEdgeSkipped:
                    if len(cs_toJoin) == 2:
                        print "Loop may be thin and didn't close."
                        return

                    sc.doc.Objects.AddCurve(c); sc.doc.Views.Redraw()
                    raise ValueError("Curve is not closed.")

                rv = tryToAdjustCrvEndsToSENW(c, rgSrf_toSplit, tolerance=fJoinTol)
                if rv:
                    c = rv

            cs_All_loops_for_split.append(c)

    #for c in cs_All_loops_for_split: sc.doc.Objects.AddCurve(c)
    #sc.doc.Views.Redraw()
    #return

    rgCrvs_Splitters_WIP = [c.Duplicate() for c in cs_All_loops_for_split]

    rc = xCurve.getCurvesToSplitSurface(
        rgCrvs_Splitters_WIP, rgSrf_toSplit, fSplitTol)
    if not rc:
        rgB_Out = rgSrf_toSplit.ToBrep()
        if Rhino.RhinoApp.ExeVersion >= 7:
            rgB_Out.Faces[0].PerFaceColor = rgF_In.PerFaceColor
        return rgB_Out

    for c in rgCrvs_Splitters_WIP: c.Dispose()
    rgCrvs_Splitters_WIP = rc

    iCt_CrvsWithShortEdgesRemoved = 0
    if fEdgeLen_Min:
        for c in rgCrvs_Splitters_WIP:
            if c.RemoveShortSegments(fEdgeLen_Min):
                iCt_CrvsWithShortEdgesRemoved += 1

    if bTrimToSegs:
        segs = xCurve.duplicateSegments(rgCrvs_Splitters_WIP)
        if segs:
            for c in rgCrvs_Splitters_WIP: c.Dispose()
            rgCrvs_Splitters_WIP = segs

    rgCrvs_Splitters = rgCrvs_Splitters_WIP

    if bOutputTrimmingCrvs: map(sc.doc.Objects.AddCurve, rgCrvs_Splitters)
    #if bDebug and not bOutputTrimmingCrvs: map(sc.doc.Objects.AddCurve, rgCrvs_Splitters)
    #if bDebug: sc.doc.Objects.AddSurface(rgSrf_toSplit)

    rc = xSurface.splitSurfaceIntoBrep(
        rgSrf_toSplit,
        rgCrvs_Splitters,
        fTolerance=fSplitTol,
        bTryOtherTolsOnFail=True,
        bDebug=bDebug)
    if rc is None:
        if bEcho: print "Split failed."
        return
    rgB_Out = rc
    if Rhino.RhinoApp.ExeVersion >= 7:
        for rgF in rgB_Out.Faces:
            rgF.PerFaceColor = rgF_In.PerFaceColor
    return rgB_Out


def retrimFace(rgFace, **kwargs):
    """
    Parameters:
        rgFace,
        rgSrf_Replacement=None, # If None, UnderlyingSurface of BrepFace is used.
        rgCrvs_Splitters_Replacement=None, # If None, 3D curves of rgFace loops are used.
        fSplitTol=sc.doc.ModelAbsoluteTolerance,
        bTrimToSegs=True,
        fEdgeLen_Min,
        bOutputSplitOnFail=False,
        bOutputTrimmingCrvs=False,
        bEcho=True,
        bDebug=False)
    Returns on success:
        rg.Brep
    Returns on fail:
        None
    """


    def setOpt(key, value): return kwargs[key] if key in kwargs else value

    rgSrf_Replacement = setOpt('rgSrf_Replacement', None)
    rgCrvs_Splitters_Replacement = setOpt('rgCrvs_Splitters_Replacement', None)
    fSplitTol = setOpt('fSplitTol', sc.doc.ModelAbsoluteTolerance)
    bTrimToSegs = setOpt('bTrimToSegs', True)
    fEdgeLen_Min = setOpt('fEdgeLen_Min', 1.8*sc.doc.ModelAbsoluteTolerance)
    bOutputSplitOnFail = setOpt('bOutputSplitOnFail', False)
    bOutputTrimmingCrvs = setOpt('bOutputTrimmingCrvs', False)
    bEcho = setOpt('bEcho', True)
    bDebug = setOpt('bDebug', False)


    def faceAtPoint(rgBrep, pt_onFace, bDebug=False):
        """
        Returns face index.
        """
        for idxF, rgFace in enumerate(rgBrep.Faces):
            b, u, v = rgFace.ClosestPoint(pt_onFace)
            if not b: raise ValueError("ClosestPoint returned None.")
            ptFaceRel = rgFace.IsPointOnFace(u, v)
            if ptFaceRel == rg.PointFaceRelation.Interior:
                return idxF
        sc.doc.Objects.AddBrep(rgBrep); sc.doc.Views.Redraw()
        raise ValueError("Face with point on its interior not found.")


    if bDebug: print "Function retrimFace():"
    
    rgFace_In = rgFace

    rgB_Split = splitFace(
        rgFace=rgFace_In,
        rgSrf_Replacement=rgSrf_Replacement,
        rgCrvs_Splitters_Replacement=rgCrvs_Splitters_Replacement,
        fSplitTol=fSplitTol,
        bTrimToSegs=bTrimToSegs,
        fEdgeLen_Min=fEdgeLen_Min,
        bOutputSplitOnFail=bOutputSplitOnFail,
        bOutputTrimmingCrvs=bOutputTrimmingCrvs,
        bEcho=bEcho,
        bDebug=bDebug)

    if rgB_Split is None:
        if bEcho: print "Split rejected/failed."
        return


    pt_onFace_In = createPoint3dOnInterior(rgFace_In)
    if not pt_onFace_In:
        if bOutputSplitOnFail: sc.doc.Objects.AddBrep(rgB_Split)
        rgB_Split.Dispose()
        return


    idxF_Pos = faceAtPoint(rgB_Split, pt_onFace_In)
    if idxF_Pos is None:
        if bOutputSplitOnFail: sc.doc.Objects.AddBrep(rgB_Split)
        rgB_Split.Dispose()
        return


    rgBrep_Retrmd = rgB_Split.Faces[idxF_Pos].DuplicateFace(duplicateMeshes=False)

    rgB_Split.Dispose()

    return rgBrep_Retrmd


def createModifiedFace(rgFace_In, surfaceFunc, bDebug=False):
    """
    Parameters:
        rgFace_In,
        surfaceFunc,
        bDebug,
    Returns on success:
        tuple(
            rg.Brep,
            float(Surface deviation from input))
        str(Feedback)
    Returns on fail:
        None,
        str(Feedback)
    """


    rgB_1F_In = rgFace_In.DuplicateFace(duplicateMeshes=True)
    area_In = rgB_1F_In.GetArea()
    rgB_1F_In.Dispose()
    if area_In <= (10.0*sc.doc.ModelAbsoluteTolerance)**2:
        return None, "Area of face is too small to process."

    rgSrf_In = rgFace_In.UnderlyingSurface()


    res, sLog = surfaceFunc(rgSrf_In)

    if not res:
        return None, sLog

    ns_Res, srf_dev = res

    if not ns_Res:
        return None, "NurbsSurface not created."
    if not ns_Res.IsValid:
        return None, "Invalid Surface geometry after Fit."

    if rgFace_In.IsSurface:
        rgBrep1_BeforeTrim = ns_Res.ToBrep()
        if not rgBrep1_BeforeTrim.IsValid:
            return None, "Invalid brep geometry after ToBrep."

        # Success.
        return (rgBrep1_BeforeTrim, srf_dev), "Trim skipped because face IsSurface."

    rgBrep1_AfterTrim = retrimFace(
        rgFace_In,
        rgSrf_Replacement=ns_Res, # If None, UnderlyingSurface of BrepFace is used.
        rgCrvs_Splitters_Replacement=None, # If None, 3D curves of rgFace loops are used.
        fSplitTol=sc.doc.ModelAbsoluteTolerance,
        bTrimToSegs=True,
        fEdgeLen_Min=1.8*sc.doc.ModelAbsoluteTolerance,
        bOutputSplitOnFail=False,
        bOutputTrimmingCrvs=False,
        bDebug=bDebug)

    if rgBrep1_AfterTrim is None:
        # Check areas.
        rgBrep1_BeforeTrim = ns_Res.ToBrep()

        area_Out = rgBrep1_BeforeTrim.GetArea()
        area_Diff = abs(area_In-area_Out)
        if area_Diff <= (10.0*sc.doc.ModelAbsoluteTolerance)**2:
            return (rgBrep1_BeforeTrim, srf_dev), "Face not trimmed because it is almost IsSurface."

        if bDebug:
            if not rgBrep1_BeforeTrim.IsValid:
                return None, "Invalid brep geometry after ToBrep.  (Split fail)."
            attrRed = rd.ObjectAttributes()
            attrRed.LayerIndex = sc.doc.Layers.CurrentLayerIndex
            attrRed.ColorSource = rd.ObjectColorSource.ColorFromObject
            attrRed.ObjectColor = Color.Red
            gBrep1 = sc.doc.Objects.AddBrep(rgBrep1_BeforeTrim, attrRed)
            if gBrep1 == Guid.Empty:
                return None, "Untrimmed monoface brep could not be added to document."
            return None, "Trim fail.  Split brep added."
        return None, "Trim fail."

    if not rgBrep1_AfterTrim.IsValid:
        if bDebug:
            rgBrep1_BeforeTrim = ns_Res.ToBrep()
            if not rgBrep1_BeforeTrim.IsValid:
                return None, "Invalid brep geometry after ToBrep.  (Split fail)."
            attrRed = rd.ObjectAttributes()
            attrRed.LayerIndex = sc.doc.Layers.CurrentLayerIndex
            attrRed.ColorSource = rd.ObjectColorSource.ColorFromObject
            attrRed.ObjectColor = Color.Red
            gBrep1 = sc.doc.Objects.AddBrep(rgBrep1_BeforeTrim, attrRed)
            if gBrep1 == Guid.Empty:
                return None, "Untrimmed monoface brep could not be added to document."
        return None, "Trim fail."

    # Success.
    return (rgBrep1_AfterTrim, srf_dev), None


