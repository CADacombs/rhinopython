"""
200619: Created.
200629: Bug fix.
200729: Import-related update.
200810: Added a couple of functions.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc

from System import Guid
from System import Random

import xBrep_splitSurfaceWithCurves
import xCurve


def formatDistance(fDistance):
    if fDistance is None:
        return "(No deviation was provided.)"
    if fDistance == 0.0:
        return "Exactly zero"
    if fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-3)):
        return "{:.1e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


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
    rgCrvs_fromLoops = [loop.To3dCurve() for loop in rgFace.Loops]

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
                for c in rgCrvs_fromLoops:
                    b, t = c.ClosestPoint(pt)
                    if not b: raise ValueError("ClosestPoint failed in createPointOnFace.")
    
                    dist = pt.DistanceTo(c.PointAt(t))
                    if dist < fMinDistFromBorder:
                        break # to get another u and v.
                else: # Good point
                    for c in rgCrvs_fromLoops: c.Dispose()
                    rgB_1F.Dispose()
                    if bDebug: sEval = "i"; print sEval+':',eval(sEval)
                    return pt
    
            # Get new parameters for point.
            u_Shrunk = rand.NextDouble() * domU_Shrunk.Length + domU_Shrunk.Min
            v_Shrunk = rand.NextDouble() * domV_Shrunk.Length + domV_Shrunk.Min

#    sc.doc.Objects.AddBrep(rgFace.DuplicateFace(duplicateMeshes=True))
#    sc.doc.Views.Redraw()

    raise ValueError(
        "Failed to find an interior point on face[{}].".format(
            rgFace.FaceIndex))


def is3dPointOnFace(rgFace, pt, fTolerance=None):
    """
    Returns:
        False when the point is at least fTolerance from the face's underlying surface.
        Otherwise, rg.PointFaceRelation.Interior or Boundary.
    """
    
    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance
    
    b, u, v = rgFace.ClosestPoint(pt)
    if not b:
        print "ClosestPoint return was not a success."
        return

    dist = rgFace.PointAt(u,v).DistanceTo(pt)
    if dist > fTolerance:
        # Point isn't even on underlying surface.
        return False

    ptFaceRel = isParameterPointOnFace(
        rgFace, u, v, fTolerance)

    return ptFaceRel


def isParameterPointOnFace(rgFace, u, v, f3dTolerance=None):
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
    
    pt_3d = rgFace.PointAt(u, v)
    # Note that this Point3d will already be pulled to the face.
    # For true 3D point check, use is3dPointOnFace instead.
    
    dist = rgFace.PointAt(u,v).DistanceTo(pt_3d)
    if dist > f3dTolerance:
        # Point must have pulled outside of face's boundary.
        return rg.PointFaceRelation.Exterior

    ptFaceRel = rgFace.IsPointOnFace(u, v)
    
    if ptFaceRel == rg.PointFaceRelation.Boundary:
        return rg.PointFaceRelation.Boundary
    
    for iE in rgFace.AdjacentEdges():
        bSuccess, t = rgFace.Brep.Edges[iE].ClosestPoint(
            testPoint=pt_3d,
            maximumDistance=f3dTolerance)
        if not bSuccess:
            continue
        pt_onEdge = rgFace.Brep.Edges[iE].PointAt(t)
        dist_toCrv = pt_3d.DistanceTo(other=pt_onEdge)
        if dist_toCrv <= f3dTolerance:
            return rg.PointFaceRelation.Boundary
    else:
        # No distances smaller than f3dTolerance found.
        return rg.PointFaceRelation.Interior


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
    fEdgeLen_Min = setOpt('fEdgeLen_Min', 2.0*sc.doc.ModelAbsoluteTolerance)
    bOutputTrimmingCrvs = setOpt('bOutputTrimmingCrvs', False)
    bEcho = setOpt('bEcho', True)
    bDebug = setOpt('bDebug', True)


    if bDebug: print "Function splitFace():"
    
    rgF_In = rgFace

    rgSrf_toSplit = rgF_In.UnderlyingSurface() if rgSrf_Replacement is None else rgSrf_Replacement
    
    # Reject surfaces that contain natural edges < 2 * fSplitTol.
    for iDir, t in zip(
        (0, 0, 1, 1),
        (rgSrf_toSplit.Domain(1).T0, rgSrf_toSplit.Domain(1).T1,
            rgSrf_toSplit.Domain(0).T0, rgSrf_toSplit.Domain(1).T1,),
    ):
        isoCrv = rgSrf_toSplit.IsoCurve(iDir, t)
        fLength = isoCrv.GetLength()
        isoCrv.Dispose()
        if 0.0 < fLength <= 2.0 * sc.doc.ModelAbsoluteTolerance:
            return


    fEdgeDevs = []
    for idxE in rgF_In.AdjacentEdges():
        fEdgeDev = rgF_In.Brep.Edges[idxE].Tolerance
        if fEdgeDev is not None:
            fEdgeDevs.append(fEdgeDev)
    fMaxEdgeDev = max(fEdgeDevs)


    # Although splits may be successful if trims include seams,
    # they are avoided for cleaner input.  Otherwise, this could be used:
    # rgCrvs_fromLoops = [loop.To3dCurve() for loop in rgFace_In.Loops]

    rgCrvs_fromLoops = []
    for loop in rgF_In.Loops:
        cs_ofLoop = []
        for trim in loop.Trims:
            # TODO: Also remove curves on surface extents?
            #       Only if new surface is not provided.
            if trim.TrimType == rg.BrepTrimType.Seam: continue
            if trim.TrimType == rg.BrepTrimType.Singular: continue
            if trim.TrimType not in (
                rg.BrepTrimType.Boundary, rg.BrepTrimType.Mated
            ):
                print trim.TrimType
            cs_ofLoop.append(trim.Edge.DuplicateCurve())
        
        #for c in cs_ofLoop: sc.doc.Objects.AddCurve(c)
        #sc.doc.Views.Redraw(); 1/0
        
        cs_WIP = []
        for c in cs_ofLoop:
            fLength = c.GetLength()
            if fLength > fSplitTol:
                cs_WIP.append(c)
            else:
                c.Dispose()
        
        cs_ofLoop = cs_WIP[:]
        
        fJoinTol = 2.01*max(fMaxEdgeDev, fSplitTol)
        
        rgCrvs_Joined = rg.Curve.JoinCurves(
            cs_ofLoop,
            joinTolerance=fJoinTol)

        #for c in rgCrvs_Joined: sc.doc.Objects.AddCurve(c)
        #    sc.doc.Views.Redraw(); 1/0

        for c in rgCrvs_Joined:
            if not c.IsClosed:
                sc.doc.Objects.AddCurve(c); sc.doc.Views.Redraw()
                raise ValueError("Curve is not closed.")
            rgCrvs_fromLoops.append(c)

    #for c in rgCrvs_fromLoops: sc.doc.Objects.AddCurve(c)
    #sc.doc.Views.Redraw()
    #return

    rgCrvs_Splitters_WIP = [c.Duplicate() for c in rgCrvs_fromLoops]

    rc = xBrep_splitSurfaceWithCurves.getCurvesToSplitSurface(
        rgCrvs_Splitters_WIP, rgSrf_toSplit, fSplitTol)
    if not rc:
        return rgSrf_toSplit.ToBrep()

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

    if bDebug or bOutputTrimmingCrvs: map(sc.doc.Objects.AddCurve, rgCrvs_Splitters)
    if bDebug: sc.doc.Objects.AddSurface(rgSrf_toSplit)

    rc = xBrep_splitSurfaceWithCurves.splitSurfaceIntoBrep(
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
        fEdgeLen_Min=2.0*sc.doc.ModelAbsoluteTolerance,
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
    fEdgeLen_Min = setOpt('fEdgeLen_Min', 2.0*sc.doc.ModelAbsoluteTolerance)
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
        fEdgeLen_Min=2.0*sc.doc.ModelAbsoluteTolerance,
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


