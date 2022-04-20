"""
191118-19: Created by moving main function from another module.
191215: Corrected debug printed output.
200401: Import-related bug fix.
200611: Implemented more correct use of IsPointOnFace.
200715: Joined 
210423: replaceShape now uses RemoveShortSegments for splitting curves.
        WIP: Make moveSurfaceSeamsToAvoidCurves avoid endpoints of curves as well.
220420: Repaired a RhinoCommon 7.17 script-breaking change.

TODO:
    Rotate RevSurface's before moving seams to avoid conversion to NurbsSurface.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

import xBrep_findMatchingFace
import xPlaneSurface

import random

maxExponentFor2 = 4


def rotateSphericalBrepToAvoidCurves(rgBrep0, rgCrvs_ToAvoid, rotationCenter, fTolerance=None, bDebug=False):
    """
    Parameters:
        rgBrep0
        rgCrvs_ToAvoid
        rotationCenter
    Returns: New brep.
    Moves surface (from sphere) seam away from rgCrvs_ToAvoid.
    Do not use this for spherical surfaces.
    """


    if fTolerance is None:
        fTolerance = max(
            2.0*sc.doc.ModelAbsoluteTolerance,
            max([e.Tolerance for e in rgBrep0.Edges])
            )

    # Join curves to be avoided.  Joining them now will reduce the number of tests.
    rgCrvs_Joined = rg.Curve.JoinCurves(rgCrvs_ToAvoid, fTolerance)
    if rgCrvs_Joined is None:
        if bDebug: print "JoinCurves failed!"
        return
    if len(rgCrvs_Joined) > 1:
        if bDebug: print "JoinCurves resulted in {} curves.".format(len(rgCrvs_Joined))
        return
    rgCrv_ToAvoid_Joined = rgCrvs_Joined[0]
    if not rgCrv_ToAvoid_Joined.IsClosed:
        if bDebug: print "Joined curve is not closed."
        return
    
    #sc.doc.Objects.AddCurve(rgCrv_ToAvoid_Joined); sc.doc.Views.Redraw()
    
    rgCrvs_rgBrep0Edges = rgBrep0.Curves3D
    #sc.doc.Objects.AddBrep(rgBrep0)

    rgCrvs_Joined = rg.Curve.JoinCurves(rgCrvs_rgBrep0Edges, fTolerance)
    if rgCrvs_Joined is None:
        if bDebug: print "JoinCurves failed!"
        return
    if len(rgCrvs_Joined) > 1:
        s  = "JoinCurves resulted in"
        s += " {} curves.".format(len(rgCrvs_Joined))
        s += "  Quitting rotation of this brep."
        print s
        return
    rgCrv_rgBrep0Edges_Joined = rgCrvs_Joined[0]
    
    #sc.doc.Objects.AddCurve(rgCrv_rgBrep0Edges_Joined); sc.doc.Views.Redraw()

    # First try starting baseball sphere.
    rc = rgCrv_ToAvoid_Joined.ClosestPoints(otherCurve=rgCrv_rgBrep0Edges_Joined)
    b, ptOnCrv_ToAvoid_Joined, ptOnBaseballEdgeCrv_Joined = rc
    if not b:
        print "ClosestPoints not found!"
        return
    dist = ptOnCrv_ToAvoid_Joined.DistanceTo(ptOnBaseballEdgeCrv_Joined)
    if dist > 10.0*max((fTolerance, sc.doc.ModelAbsoluteTolerance)):
        # TODO: Is 10.0* a good multiplier?
        return rgBrep0
        
    angles = [0.0,0.0,0.0]
    for N in xrange(maxExponentFor2):
        sc.escape_test()
        denominator = 2.0**N
        for numerator in xrange(1, denominator, 2):
            sc.escape_test()
            perunus = float(numerator) / float(denominator)
            angles[0] = Rhino.RhinoMath.ToRadians(360.0)*perunus
            for N in xrange(maxExponentFor2):
                sc.escape_test()
                denominator = 2.0**N
                for numerator in xrange(1, denominator, 2):
                    sc.escape_test()
                    perunus = float(numerator) / float(denominator)
                    angles[1] = Rhino.RhinoMath.ToRadians(360.0)*perunus
                    for N in xrange(maxExponentFor2):
                        sc.escape_test()
                        denominator = 2.0**N
                        for numerator in xrange(1, denominator, 2):
                            perunus = float(numerator) / float(denominator)
                            angles[2] = Rhino.RhinoMath.ToRadians(360.0)*perunus
                            rgBaseballEdgeCrv_Joined_WIP = rgCrv_rgBrep0Edges_Joined.Duplicate()
                            rgBaseballEdgeCrv_Joined_WIP.Rotate(
                                    angleRadians=angles[0],
                                    rotationAxis=rg.Vector3d.XAxis,
                                    rotationCenter=rotationCenter)
                            rgBaseballEdgeCrv_Joined_WIP.Rotate(
                                    angleRadians=angles[1],
                                    rotationAxis=rg.Vector3d.YAxis,
                                    rotationCenter=rotationCenter)
                            rgBaseballEdgeCrv_Joined_WIP.Rotate(
                                    angleRadians=angles[2],
                                    rotationAxis=rg.Vector3d.ZAxis,
                                    rotationCenter=rotationCenter)
                                
                            rc = rgCrv_ToAvoid_Joined.ClosestPoints(otherCurve=rgBaseballEdgeCrv_Joined_WIP)
                            b, ptOnCrv_ToAvoid_Joined, ptOnBaseballEdgeCrv_Joined = rc
                            if not b:
                                print "ClosestPoints not found!"
                                return
                            dist = ptOnCrv_ToAvoid_Joined.DistanceTo(ptOnBaseballEdgeCrv_Joined)
                            if dist <= 10.0*max((fTolerance, sc.doc.ModelAbsoluteTolerance)):
                                # TODO: Is 10.0* a good multiplier?
                                #sc.doc.Objects.AddSphere(sphereRotX)
                                continue # to next axis rotation.
                                
                            # Found good orientation.
                            #sc.doc.Objects.AddCurve(rgBaseballEdgeCrv_Joined_WIP)
                            rgBrep_Rotated = rgBrep0.Duplicate()
                            if not rgBrep_Rotated.Rotate(
                                angleRadians=angles[0],
                                rotationAxis=rg.Vector3d.XAxis,
                                rotationCenter=rotationCenter
                            ):
                                return
                            if not rgBrep_Rotated.Rotate(
                                angleRadians=angles[1],
                                rotationAxis=rg.Vector3d.YAxis,
                                rotationCenter=rotationCenter
                            ):
                                return
                            if not rgBrep_Rotated.Rotate(
                                angleRadians=angles[2],
                                rotationAxis=rg.Vector3d.ZAxis,
                                rotationCenter=rotationCenter
                            ):
                                return
                            return rgBrep_Rotated


def moveSurfaceSeamsToAvoidCurves(rgSrf0, rgCrvs_ToAvoid, iDirections, tolerance=None, bDebug=False):
    """
    Parameters:
        rgCrvs_ToAvoid: Curves from which to move the seam.
        
    Returns: New Surface with modified seam, if necessary.
        
    Do not use this for spherical surfaces.
    """
        
    rgBrep_1F_ToMoveSrfSeam = rgSrf0.ToBrep()
    if rgBrep_1F_ToMoveSrfSeam is None: return
        
    if tolerance is None:
        tolerance = 4.0*sc.doc.ModelAbsoluteTolerance
    
    # Test at existing seam location(s) in case moving the seam is not necessary.
    # This also gets a base count for the number of faces from split.
    rgFace_ToMoveSeam = rgBrep_1F_ToMoveSrfSeam.Faces[0]
        
    rgCrvs_Joined = rg.Curve.JoinCurves(rgCrvs_ToAvoid)
    if rgCrvs_Joined is None:
        if bDebug: print "JoinCurves failed!"
        return
    #for c in rgCrvs_Joined: sc.doc.Objects.AddCurve(c)
        
    ct_rgCrvs_Joined = len(rgCrvs_Joined)
    if ct_rgCrvs_Joined > 1:
        if bDebug: s = "JoinCurves resulted in {} curves.".format(len(rgCrvs_Joined))
        bOpen = [True for c in rgCrvs_Joined if not c.IsClosed]
        if bOpen:
            if bDebug: 
                s += "  {} are open.".format(len(bOpen))
                s += "  This merge will be skipped."
                print s
            return
        # Check whether they probably should be joined.
        ct_TooClose = 0
        for i in range(len(rgCrvs_Joined)):
            for j in range(i+1, len(rgCrvs_Joined)):
                cI = rgCrvs_Joined[i]
                cJ = rgCrvs_Joined[j]
                b, ptI, ptJ = cI.ClosestPoints(cJ)
                if not b:
                    if bDebug:
                        s += "  Error in checking distances between curves."
                        print s
                dist = ptI.DistanceTo(ptJ)
                if dist <= tolerance:
                    if bDebug:
                        s += "  Distance between 2 curves is only {}.".\
                                format(dist)
                        s += "  This merge will be skipped."
                        print s
                    return
    elif not rg.Curve.IsClosed:
        if bDebug: print "Joined curve is not closed."
        return
        
    rgBrep_NewSeamLoc_WIP = rgSrf0.ToBrep()
    rgFace_ToMoveSeam = rgBrep_NewSeamLoc_WIP.Faces[0]
        
    for iDir in iDirections:
        if not rgFace_ToMoveSeam.IsClosed(iDir): continue
        def findLeastIntersectingSeamParameter(iDir):
            parameter_LeastInters = None
            iCt_Inters_Least = None
            for N in xrange(maxExponentFor2):
                denominator = 2.0**N
                numerator_range = (0,) if denominator == 1 else xrange(1, denominator, 2)
                for numerator in numerator_range:
                    perunus = float(numerator) / float(denominator)
                    T0 = rgFace_ToMoveSeam.Domain(direction=iDir).T0
                    T1 = rgFace_ToMoveSeam.Domain(direction=iDir).T1
                    parameter = T0*(1-perunus)+T1*perunus
                    rgSeamCrv = rgFace_ToMoveSeam.IsoCurve(direction=(iDir+1)%2, constantParameter=parameter) # direction is opposite of constantParameter's direction.
                    #sc.doc.Objects.AddCurve(rgSeamCrv); sc.doc.Views.Redraw()
                    iCt_Inters = 0
                    for rgCrv_Joined in rgCrvs_Joined:
                        #sc.doc.Objects.AddCurve(rgCrv_Joined); sc.doc.Views.Redraw()
                        rc = rgCrv_Joined.ClosestPoints(otherCurve=rgSeamCrv)
                        b, ptOnBoundaryCrv, pointOnSeamCurve = rc
                        if not b:
                            if bDebug: print "ClosestPoints not found!"
                            return
                        dist = ptOnBoundaryCrv.DistanceTo(other=pointOnSeamCurve)
                        if dist <= tolerance:
                            iCt_Inters += 1
                    if iCt_Inters == 0:
                        return parameter
                    elif iCt_Inters_Least is None or iCt_Inters < iCt_Inters_Least:
                            parameter_LeastInters = parameter
                            iCt_Inters_Least = iCt_Inters
            if iCt_Inters_Least == len(rgCrvs_Joined):
                return
            else:
                return parameter
        parameter = findLeastIntersectingSeamParameter(iDir)
        if parameter is None:
            if len(iDirections) == 1:
                # Seam may have to be present as in the case of a
                # cylindrical or conical face with a full circumference.
                rgSrf1 = rgBrep_NewSeamLoc_WIP.Surfaces[0].Duplicate()
                rgBrep_NewSeamLoc_WIP.Dispose()
                return rgSrf1
            if bDebug: print "Surface is closed in both directions, so trying to move the other seam first."
            for iDir in reversed(iDirections):
                parameter = findLeastIntersectingSeamParameter(iDir)
                if parameter is None:
                    if bDebug: print "Parameter not found!"
                    return
                # Check whether seam is different than current.
                if parameter == rgFace_ToMoveSeam.Domain(direction=iDir).T0:
                    rgBrep_NewSeamLoc_WIP
                else:
                    rgBrep_NewSeamLoc_WIP = rg.Brep.ChangeSeam(
                            face=rgFace_ToMoveSeam,
                            direction=iDir,
                            parameter=parameter,
                            tolerance=sc.doc.ModelAbsoluteTolerance)
                    if rgBrep_NewSeamLoc_WIP is None:
                        if bDebug: print "ChangeSeam resulted in None!"
                        return
                # For a Torus, move second seam after moving first seam.
                if len(iDirections) == 2:
                    rgFace_ToMoveSeam = rgBrep_NewSeamLoc_WIP.Faces[0]
            return rgBrep_NewSeamLoc_WIP
        # Move seam if the good parameter is different than current seam.
        if parameter == rgFace_ToMoveSeam.Domain(direction=iDir).T0:
            rgBrep_NewSeamLoc_WIP
        else:
            rgBrep_NewSeamLoc_WIP = rg.Brep.ChangeSeam(
                    face=rgFace_ToMoveSeam,
                    direction=iDir,
                    parameter=parameter,
                    tolerance=sc.doc.ModelAbsoluteTolerance)
            if rgBrep_NewSeamLoc_WIP is None:
                if bDebug: print "ChangeSeam resulted in None!"
                return
        # For a Torus, move second seam after moving first seam.
        if len(iDirections) == 2:
            rgFace_ToMoveSeam = rgBrep_NewSeamLoc_WIP.Faces[0]
    
    rgSrf1 = rgBrep_NewSeamLoc_WIP.Surfaces[0].Duplicate()
    rgBrep_NewSeamLoc_WIP.Dispose()
    return rgSrf1


def getRandomPointOnFace(rgBrep, idxFace, fDistMin=None):
    """
    Parameters:
        rgBrep: Used to get the edges from the edge indices.
        idxFace
        fDistMin: Minimum distance from the face border.
    Returns:
        Point3d on success, None on failure
    """
    
    rgFace = rgBrep.Faces[idxFace]
    
    if fDistMin is None:
        fDistMin = 10.0*sc.doc.ModelAbsoluteTolerance
    
    areaMassProp = Rhino.Geometry.AreaMassProperties.Compute(rgFace)
    if areaMassProp is None:
        print "Face index {} skipped because its AreaMassProperties, " \
        "and thus its centroid, cannot be calculated.".format(
                idxFace)
        return
    
    ptCentrdW = areaMassProp.Centroid
    getrc, u, v = rgFace.ClosestPoint(ptCentrdW)
    domnU, domnV = rgFace.Domain(0), rgFace.Domain(1)
    fIE = 0. # Ratio from surface borders
    
    rgEdges = [rgBrep.Edges[idxEdge] for idxEdge in rgFace.AdjacentEdges()]
    
    for i in xrange(100000): # If point is not on face, continue searching.
        ptFaceRel = rgFace.IsPointOnFace(u, v)
        if ptFaceRel == rg.PointFaceRelation.Interior:
            # If point is not at least fDistMin from border, continue searching.
            pt = rgFace.PointAt(u, v)
            for rgEdge in rgEdges:
                b, t = rgEdge.ClosestPoint(pt)
                if b:
                    if pt.DistanceTo(rgEdge.PointAt(t)) < fDistMin: # Too close
                        break
            else: # Good point
                map(lambda x: x.Dispose(), rgEdges)
                rgFace.Dispose()
                return pt
        
        # Get new parameters for point.
        u = random.uniform(domnU.T0 + fIE*domnU.Length,
                domnU.T1 - fIE*domnU.Length)
        v = random.uniform(domnV.T0 + fIE*domnV.Length,
                domnV.T1 - fIE*domnV.Length)
    
    map(lambda x: x.Dispose(), rgEdges)
    rgFace.Dispose()


def getFaceAtPoint(rgBrep_Split, ptOnFace, bEcho=False):
    """
    """

    if bEcho: print 'getFaceAtPoint()...'
    
    # Test the point distance to each face.
    for iF, rgFace in enumerate(rgBrep_Split.Faces):
        rgBrep1_1F = rgFace.DuplicateFace(False)
        if rgBrep1_1F is None:
            if bEcho: sPrint = 'rgBrep1_1F'; print sPrint + ':', eval(sPrint)
            return
        
        b, u, v = rgFace.ClosestPoint(ptOnFace)
        if b:
            # Shouldn't need to implement this if IsPointOnFace correctly
            # returns rg.PointFaceRelation.Interior.
            #pt_Closest = rgFace.PointAt(u, v)
            #fDist = pt_Closest.DistanceTo(ptOnFace)

            ptFaceRel = rgFace.IsPointOnFace(u, v)
            if ptFaceRel == rg.PointFaceRelation.Interior:
                rgFace.Dispose()
                rgBrep1_1F.Dispose()
                return iF
        
        rgFace.Dispose()
        rgBrep1_1F.Dispose()


def extendCylinderToObjectSize(rgCyl, rgObjForSizeRef, bDebug=False):
    rgBbox_ForSizeRef = rgObjForSizeRef.GetBoundingBox(accurate=False)
    fAddLen = 1.1 * rgBbox_ForSizeRef.Min.DistanceTo(rgBbox_ForSizeRef.Max)

    if not rgCyl.IsFinite: # Fix cylinder that is infinite (has 0 height).
        if bDebug: print "Making cylinder finite...",
        rgCyl.Height1 = -fAddLen
        rgCyl.Height2 = fAddLen
        if bDebug: print "Height of cylinder after: {}".format(rgCyl.TotalHeight)
    else: # Increase length of cylinder so that it extends the reference object.
        if bDebug: print "Height of cylinder before: {}".format(rgCyl.TotalHeight),
        if rgCyl.Height1 < rgCyl.Height2:
            rgCyl.Height1 -= fAddLen
            rgCyl.Height2 += fAddLen
        else:
            rgCyl.Height1 += fAddLen
            rgCyl.Height2 -= fAddLen
        if bDebug: print "Height of cylinder after: {}".format(rgCyl.TotalHeight)


def replaceShape(rgBrep0, shape, fTolerance=None, bDebug=False):
    """
    shape_ForMerge can be NurbsSurface, Cone, Cylinder, Plane, Sphere, or Torus.
    
    Returns: New Brep.
    """

    if bDebug: print '-'*80 + '\n' + 'replaceShape()'
    
    rgBrep_In = rgBrep0

    shape_In = shape
    
    if not shape_In.IsValid:
        print "{} is NOT valid in replaceShape!".format(shape_In)

    if fTolerance is None:
        fTolerance = max(
            2.0*sc.doc.ModelAbsoluteTolerance,
            max([e.Tolerance for e in rgBrep_In.Edges])
            )

    # 200715: Disabled rgBrep_0RebuiltEs.
    #    rgBrep_0RebuiltEs = rgBrep_In.DuplicateBrep()
    #    for f in rgBrep_0RebuiltEs.Faces:
    #        f.RebuildEdges(
    #                tolerance=fTolerance,
    #                rebuildSharedEdges=True,
    #                rebuildVertices=True)

    # Prepare trimming curves.
    rgCrvs_NEs_Otr = rgBrep_In.DuplicateNakedEdgeCurves(True, False) # nakedOuter, nakedInner
    
    #map(sc.doc.Objects.AddCurve, rgCrvs_ToJoin_Otr); sc.doc.Views.Redraw()
    
    #    rgCrvs_ToJoin_Otr = []
    #    for c in rgCrvs_NEs_Otr:
    #        if c.PointAtStart.DistanceTo(c.PointAtEnd) > fTolerance:
    #            rgCrvs_ToJoin_Otr.append(c)
    #        else:
    #            c.Dispose()
    
    rgCrvs_ToJoin_Otr = rgCrvs_NEs_Otr[:]
    
    #map(sc.doc.Objects.AddCurve, rgCrvs_ToJoin_Otr); sc.doc.Views.Redraw()
    
    rgCrvs_NEs_Inr = rgBrep_In.DuplicateNakedEdgeCurves(False, True) # nakedOuter, nakedInner
    
    #map(sc.doc.Objects.AddCurve, rgCrvs_Inr); sc.doc.Views.Redraw()
    
    #    rgCrvs_ToJoin_Inr = []
    #    for c in rgCrvs_NEs_Inr:
    #        if c.PointAtStart.DistanceTo(c.PointAtEnd) > fTolerance:
    #            rgCrvs_ToJoin_Otr.append(c)
    #        else:
    #            c.Dispose()
    
    rgCrvs_ToJoin_Inr = rgCrvs_NEs_Inr[:]
    
    #map(sc.doc.Objects.AddCurve, rgCrvs_ToJoin_Inr); sc.doc.Views.Redraw()
    
    # Join curves because it has been found that small segments can cause
    # BrepFace.Split to fail.
    # A side effect of doing this is that some edges will be merged.
    
    rgCrvs_Joined_Otr = rg.Curve.JoinCurves(rgCrvs_ToJoin_Otr, fTolerance)

    #map(sc.doc.Objects.AddCurve, rgCrvs_Joined_Otr); sc.doc.Views.Redraw()
    
    rgCrvs_Joined_Inr = rg.Curve.JoinCurves(rgCrvs_ToJoin_Inr, fTolerance)

    #map(sc.doc.Objects.AddCurve, rgCrvs_Joined_Inr); sc.doc.Views.Redraw()
    
    rgCrvs_Joined_All = rg.Curve.JoinCurves(
        list(rgCrvs_Joined_Otr) + list(rgCrvs_Joined_Inr),
        fTolerance)
    
    for c in rgCrvs_Joined_All:
        if not c.IsClosed:
            #sc.doc.Objects.AddCurve(c); sc.doc.Views.Redraw(); 1/0
            print "Curve is not closed in xBrep.replaceShape."
            return
    
    #rgCrvs_Joined_All = rgCrvs_Joined_Otr + rgCrvs_Joined_Inr
    
    if bDebug: sEval = 'len(rgCrvs_Joined_All)'; print sEval + ':', eval(sEval)
    #map(sc.doc.Objects.AddCurve, rgCrvs_Joined_All); sc.doc.Views.Redraw(); 1/0
    
    # Determine how many faces into which the brep should split.
    # As in the case of a complete cylinder, the outer loop can have more than 1 trim set.
    iCt_Crvs_Joined_Otr = rgCrvs_Joined_Otr.Count
    if bDebug: sEval = 'iCt_Crvs_Joined_Otr'; print sEval + ':', eval(sEval)
    numTargetSplits_Max = (iCt_Crvs_Joined_Otr + rgBrep_In.Loops.Count -
            rgBrep_In.Faces.Count + 1)
    if bDebug: sEval = 'numTargetSplits_Max'; print sEval + ':', eval(sEval)
    
    rgCrvs_ForSplit = []
    for c in rgCrvs_Joined_All:
        print c.RemoveShortSegments(tolerance=fTolerance)
        segs = c.DuplicateSegments()
        if segs:
            rgCrvs_ForSplit.extend(c.DuplicateSegments())
        else:
            rgCrvs_ForSplit.append(c)
            print "DuplicateSegments returned None."
            # Because DuplicateSegments returns None when curve cannot be segmented.
    if bDebug: sEval = 'len(rgCrvs_ForSplit)'; print sEval + ':', eval(sEval)
    
    # Get final brep.
    
    sShape = shape_In.GetType().Name
    
    if sShape == 'Plane':
        rgPlaneSrf1 = xPlaneSurface.createFromPlaneAndObjectSize(
                rgPlane=shape_In, obj_ForSize=rgBrep_In)
        rgBrep_FullSrf = rgPlaneSrf1.ToBrep()
        rgFace_FullSrf = rgBrep_FullSrf.Faces[0]
        
        rgBrep_Split = rgFace_FullSrf.Split(
                curves=rgCrvs_ForSplit,
                tolerance=sc.doc.ModelAbsoluteTolerance)

    elif sShape == 'NurbsSurface':
        rgNurbsSrf1 = shape_In.Duplicate()
        rgBrep_FromNurbsSrf = rgNurbsSrf1.ToBrep()
        rgFace_FullSrf = rgBrep_FromNurbsSrf.Faces[0]
        
        directions = []
        if rgNurbsSrf1.IsClosed(0):
            directions.append(0)
        if rgNurbsSrf1.IsClosed(1):
            directions.append(1)
        if directions:
            rgSrf1_MovedSeam = moveSurfaceSeamsToAvoidCurves(
                    rgNurbsSrf1,
                    rgCrvs_Joined_Otr,
                    directions,
                    bDebug=bDebug)
            if rgSrf1_MovedSeam is None: return
            rgFace_FullSrf = rgSrf1_MovedSeam.ToBrep().Faces[0]
        
        rgBrep_Split = rgFace_FullSrf.Split(
                curves=rgCrvs_ForSplit,
                tolerance=sc.doc.ModelAbsoluteTolerance)
    
    else:
        # Primitive of RevSurface.
        
        # Check for case of sphere or torus not requiring any trim.
        if len(rgCrvs_Joined_All) == 0:
            if sShape == 'Sphere' or sShape == 'Torus':
                rgRevSrf1_ForSolid = shape_In.ToRevSurface()
                rgBrep_FromRevSrf = rgRevSrf1_ForSolid.ToBrep()
                rgRevSrf1_ForSolid.Dispose()
                return rgBrep_FromRevSrf

        # Extend rgShapeX of cylinder or cone.
        if sShape == 'Cylinder':
            extendCylinderToObjectSize(shape_In, rgBrep_In, bDebug=bDebug)
            #sc.doc.Objects.AddBrep(shape_In.ToBrep(False, False))
        elif sShape == 'Cone': # Extend length of cone so that surface is larger than trim.
            fScale = 1.1
            shape_In.Radius *= fScale
            shape_In.Height *= fScale

        rgRevSrf1_FromShape = shape_In.ToRevSurface() # Shape is primitive of RevSurface

        if isinstance(shape_In, rg.Cylinder) or isinstance(shape_In, rg.Cone):
            rgBrep_FullCyl_FromRevSrf = rgRevSrf1_FromShape.ToBrep()
            rgFace_FullSrf = rgBrep_FullCyl_FromRevSrf.Faces[0]
            
            if rgFace_FullSrf.IsClosed(0):
                directions = 0,
            elif rgFace_FullSrf.IsClosed(1):
                directions = 1,
            else:
                print "No seam on {}!".format(shape_In)
                return
            rgSrf1_MovedSeam = moveSurfaceSeamsToAvoidCurves(
                    rgRevSrf1_FromShape,
                    rgCrvs_Joined_Otr,
                    directions,
                    bDebug=bDebug)
            if rgSrf1_MovedSeam is None: return
            #sc.doc.Objects.AddSurface(rgSrf1_MovedSeam)
            
            rgB_FullSrf = rgSrf1_MovedSeam.ToBrep()
            
            rgFace_FullSrf = rgB_FullSrf.Faces[0]
            
            rgBrep_Split = rgFace_FullSrf.Split(
                    curves=rgCrvs_ForSplit,
                    tolerance=sc.doc.ModelAbsoluteTolerance)

            if not rgBrep_Split.IsValid:
                if bDebug:
                    print "rgBreps3_1F is invalid."
                print rgBrep_Split.IsValidWithLog()
                rgBrep_Split.Dispose()
                sc.doc.Objects.AddBrep(rgB_FullSrf)
                for crv in rgCrvs_ForSplit:
                    sc.doc.Objects.AddCurve(crv)
                return

        elif isinstance(shape_In, rg.Torus):
            rgBrep_FullTorus_FromRevSrf = rgRevSrf1_FromShape.ToBrep()
            rgFace_FullSrf = rgBrep_FullTorus_FromRevSrf.Faces[0]
            
            directions = 0,1
            rgSrf1_MovedSeam = moveSurfaceSeamsToAvoidCurves(
                    rgRevSrf1_FromShape,
                    rgCrvs_Joined_Otr,
                    directions,
                    bDebug=bDebug)
            if rgSrf1_MovedSeam is None: return
            #            sc.doc.Objects.AddBrep(rgBrep_FullTorus_FromRevSrf)
            #            for rgCrvX in rgCrvs_ForSplit:
            #                sc.doc.Objects.AddCurve(rgCrvX)
            
            if isinstance(rgSrf1_MovedSeam, rg.Brep):
                pass
            
            rgFace_FullSrf = rgSrf1_MovedSeam.ToBrep().Faces[0]
            
            rgBrep_Split = rgFace_FullSrf.Split(
                    curves=rgCrvs_ForSplit,
                    tolerance=sc.doc.ModelAbsoluteTolerance)
    
        elif isinstance(shape_In, rg.Sphere):
            sphere = shape_In
            
            rgBrep_Split = None
            
            # createBrepWithAdjustedSurfaceSeam doesn't
            # work as well with spheres due to poles, so try baseball first.

            rgBrep_Baseball0 = rg.Brep.CreateBaseballSphere(
                    center=sphere.Center,
                    radius=sphere.Radius,
                    tolerance=Rhino.RhinoMath.ZeroTolerance)
            
            #sc.doc.Objects.AddBrep(rgBrep_Baseball0); sc.doc.Views.Redraw()

            rgBrep_BaseballSphere_Rotated = rotateSphericalBrepToAvoidCurves(
                    rgBrep0=rgBrep_Baseball0,
                    rgCrvs_ToAvoid=rgCrvs_Joined_Otr,
                    rotationCenter=sphere.Center,
                    fTolerance=fTolerance,
                    bDebug=bDebug)
            if rgBrep_BaseballSphere_Rotated is not None:
                # Find correct Face of baseball sphere brep.
                pt_ToFindFace = rgCrvs_Joined_Otr[0].PointAtStart
                for iF, rgFace1 in enumerate(rgBrep_BaseballSphere_Rotated.Faces):
                    b, u,v = rgFace1.ClosestPoint(pt_ToFindFace)
                    if not b: return

                    fDist = pt_ToFindFace.DistanceTo(rgFace1.PointAt(u,v))
                    if fDist > sc.doc.ModelAbsoluteTolerance:
                        continue

                    ptFaceRel = rgFace1.IsPointOnFace(u,v)
                    if ptFaceRel == rg.PointFaceRelation.Interior:
                        rgBrep_Split = rgFace1.Split(
                                curves=rgCrvs_ForSplit,
                                tolerance=sc.doc.ModelAbsoluteTolerance)
                        break
            if rgBrep_Split is None:
                rgBrep_FullSphere_FromRevSrf = rgRevSrf1_FromShape.ToBrep()
                rgFace_FullSrf = rgBrep_FullSphere_FromRevSrf.Faces[0]
                
                if rgFace_FullSrf.IsClosed(0):
                    directions = 0,
                elif rgFace_FullSrf.IsClosed(1):
                    directions = 1,
                
                rgSrf1_MovedSeam = moveSurfaceSeamsToAvoidCurves(
                        rgRevSrf1_FromShape,
                        rgCrvs_Joined_Otr,
                        directions,
                        bDebug=bDebug)
                if rgSrf1_MovedSeam is not None:
                    rgFace_FullSrf = rgSrf1_MovedSeam.ToBrep().Faces[0]
                    rgBrep_Split = rgFace_FullSrf.Split(
                            curves=rgCrvs_ForSplit,
                            tolerance=sc.doc.ModelAbsoluteTolerance)
                
                if rgBrep_Split is None:
                    rgBrep_FullSphere_FromRevSrf_Rotated = rotateSphericalBrepToAvoidCurves(
                            rgBrep0=rgBrep_FullSphere_FromRevSrf,
                            rgCrvs_ToAvoid=rgCrvs_Joined_Otr,
                            rotationCenter=sphere.Center,
                            fTolerance=fTolerance,
                            bDebug=bDebug)
                    if rgBrep_FullSphere_FromRevSrf_Rotated is not None:
                        rgFace_FullSrf = rgBrep_FullSphere_FromRevSrf_Rotated.Faces[0]
                        rgBrep_Split = rgFace_FullSrf.Split(
                                curves=rgCrvs_ForSplit,
                                tolerance=sc.doc.ModelAbsoluteTolerance)
                        if rgBrep_Split is None:
                            print "Split for non-seam crossing trim surfaces for sphere cannot be obtained!"
                            sc.doc.Objects.AddBrep(rgBrep_FullSphere_FromRevSrf)
                            rgBrep_FullSphere_FromRevSrf.Dispose()
                            return
                    
                    if rgBrep_Split is None:
                        print "Non-seam crossing trim surfaces for sphere cannot be obtained!"
                        sc.doc.Objects.AddBrep(rgBrep_FullSphere_FromRevSrf)
                        rgBrep_FullSphere_FromRevSrf.Dispose()
                        return
                rgBrep_FullSphere_FromRevSrf.Dispose()
        else:
            print "Skipping this shape:"
            sEval = 'sShape'; print sEval + ':', eval(sEval)
            return
        
    if rgBrep_Split is None:
        sEval = 'rgBrep_Split'; print sEval + ':', eval(sEval)
        return
    if bDebug:
        print "Find newly trimmed face:"
        sEval = 'rgBrep_Split'; print sEval + ':', eval(sEval)
        sEval = 'rgBrep_Split.Faces.Count'; print sEval + ':', eval(sEval)
    
    idx_rgFace_Pos = None
    ptOnFace = getRandomPointOnFace(
            rgBrep=rgBrep_In,
            idxFace=0,
            fDistMin=2.0*sc.doc.ModelAbsoluteTolerance)
    if not ptOnFace:
        ptOnFace = getRandomPointOnFace(
                rgBrep=rgBrep_In,
                idxFace=0,
                fDistMin=1.0*sc.doc.ModelAbsoluteTolerance)
    if ptOnFace:
        #sc.doc.Objects.AddPoint(ptOnFace)
        idx_rgFace_Pos = getFaceAtPoint(
                rgBrep_Split, ptOnFace, bDebug)
    
    if idx_rgFace_Pos is None:
        idx_rgFace_Pos = xBrep_findMatchingFace.usingBoundingBoxOfBrep(
                rgBrep_Split, rgBrep_In)
        if idx_rgFace_Pos is None:
            idx_rgFace_Pos = xBrep_findMatchingFace.usingBoundingBoxOfEdges(
                    rgBrep_Split, rgBrep_In)
            if idx_rgFace_Pos is None:
                return
    
    rgBrep_New1F = rgBrep_Split.Faces[idx_rgFace_Pos].DuplicateFace(False)
    #sc.doc.Objects.AddBrep(rgBrep_New1F); 1/0
    rgBrep_Split.Dispose()
    
    if bDebug: sEval = 'rgBrep_New1F'; print sEval + ':', eval(sEval)
    
    if rgBrep_New1F is None:
        if bDebug: print "rgBreps3_1F is None."
        return
    if not rgBrep_New1F.IsValid:
        if bDebug:
            print "rgBreps3_1F is invalid."
            print rgBrep_New1F.IsValidWithLog()
        rgBrep_New1F.Dispose()
        return
    
    return rgBrep_New1F


