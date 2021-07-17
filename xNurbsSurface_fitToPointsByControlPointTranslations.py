"""
Started with https://github.com/mcneel/rhinopython101manual/raw/master/Scripts/DistributedSurfaceFitter_neat.py
190130: Corrected and improved speed.
190205-21: -
190322: Added fDivLen
190323: -
190409: Changed some printed output.
191010: VS normalized line ending.
191011: Removed a function.
191017: Added DisplyConduit for drawing iteration.  Some function and variable name changes.
191017-18: Added a function for ZAxis-only translations.

TODO:
    Further test determining magnitudes using U and V separately versues composite UV distances.
"""

from clr import StrongBox
from System import Array
from System import Object
from System.Drawing import Color

import Rhino
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc


PI = Rhino.RhinoMath.ToRadians(180.0)

sDebugRange = ':20' # '6:ct_cp-6'


class DrawSurfaceIterationDisplayConduit(Rhino.Display.DisplayConduit):

    def __init__(self, srfForBbox, ptsForBBox):
        self.srf = None
        self.color = Color.Blue
        self.material = Rhino.Display.DisplayMaterial()
        self.material.Diffuse = self.color
        self.bbox = srfForBbox.GetBoundingBox(accurate=True)
        for pt in ptsForBBox:
            self.bbox.Union(point=pt)

        # this is called every frame inside the drawing code so try to do as little as possible
        # in order to not degrade display speed. Don't create new objects if you don't have to as this
        # will incur an overhead on the heap and garbage collection.

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        # Since we are dynamically drawing geometry, we needed to override
        # CalculateBoundingBox. Otherwise, there is a good chance that our
        # dynamically drawing geometry would get clipped.
    
        # Include the mesh's bbox with the scene's bounding box
        calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):
        drawEventArgs.Display.DrawSurface(
                surface=self.srf,
                wireColor=self.color,
                wireDensity=1)


def fit_ZAxisProjected_Distributed(rgNurbsSrf, pts_Target, fDistTol=0.1*sc.doc.ModelAbsoluteTolerance, sSamplingMode='Distributed', sWeightCurve='Linear', fPercentOfDomain=50.0, bShowIterations=False, bDebug=False):
    """
    Translate control points along Z axis only.
    """
    
    if not rgNurbsSrf.IsValid:
        print "Geometry of NurbsSurface to morph is invalid!"
        return
    
    ns_WIP = rgNurbsSrf.Duplicate() # Keep copy of original because ns_WIP will be modified per iteration.
    
    # Make domains of both directions [0,1].
    interval01 = rg.Interval(t0=0.0, t1=1.0)
    ns_WIP.SetDomain(direction=0, domain=interval01)
    ns_WIP.SetDomain(direction=1, domain=interval01)
    
    ct_cp = ns_WIP.Points.CountU * ns_WIP.Points.CountV
    
    nss_Morphed = []
    deviations_PostMorph_AllIters_Max = []
    transLengths_AllIters_Max = []
    
    if bDebug: print "Using {} sampling mode with {} weight curve.".format(sSamplingMode, sWeightCurve)


    def create_pts_projd_to_srf(ns, pts0):
        """Try to pull all target points to surface."""
        brep = ns.ToBrep()
        b_Project_Successes = []
        pts_Projd = []
        pts0_Hits = []
        for pt0 in pts0:
            rc = rg.Intersect.Intersection.ProjectPointsToBreps(
                    breps=[brep],
                    points=[pt0],
                    direction=rg.Vector3d.ZAxis,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
            if not rc:
                #
                # TODO: Test angle of closest point with that of the normal at UV.
                # For angles within ?, add ClosestPoints to pt_Pulled.
                #
                b_Project_Successes.append(False)
            else:
                b_Project_Successes.append(True)
                pts_Projd.append(rc[0])
                pts0_Hits.append(pt0)
        brep.Dispose()
        return b_Project_Successes, pts_Projd, pts0_Hits

    rc = create_pts_projd_to_srf(ns_WIP, pts_Target)
    b_Projd_successes, pts_Projd, pts_Target_with_proj_hit = rc
    
    if not pts_Projd:
        print "No target points can be projected to surface."
        return
    if bDebug: print "{} target points were projected to surface.".format(len(pts_Projd))
    
    
    def roundedList(i_list):
        o_list = []
        for i in i_list:
            if isinstance(i, (tuple, list, rg.Point2d, rg.Point3d, rg.Vector3d)):
                o_list.append([])
                for j in i:
                    o_list[-1].append(str(round(j,9)))
            else:
                o_list.append(str(round(i,9)))
        return o_list
    
    
    def get_uvs_of_Point3ds(srf, pts):
        """Apply Surface.ClosestPoint to a list of Point3ds."""
        uvs = []
        for pt in pts:
            rc, u, v = srf.ClosestPoint(pt)
            if not rc:
                print "ClosestPoint returned None"
                return
            #            if Rhino.RhinoMath.EpsilonEquals(u, 0.0, 1e-12): u = 0.0
            #            if Rhino.RhinoMath.EpsilonEquals(v, 0.0, 1e-12): v = 0.0
            uvs.append((u, v))
        return uvs
    uvs_Target_projd_to_ns = get_uvs_of_Point3ds(
            ns_WIP, pts_Projd) # Only targets that pulled to surface.
    if bDebug: sEval = 'uvs_Target_projd_to_ns[:12]'; print '{}:{}'.format(sEval,roundedList(eval(sEval)))
    
    def get_ws_TargetsFromProjected(ns, pts_ProjdFrom, pts_ProjdHit, uvs_ProjHits):
        """
        For every target point that can be pulled to NurbsSurface,
        w, the signed Euclidian distance to surface, will be generated.
        """
        dists_z = []
        for pt_ProjdFrom, pt_ProjHit, uv_atHit in zip(pts_ProjdFrom, pts_ProjdHit, uvs_ProjHits):
            normal_at_T = ns.NormalAt(*uv_atHit) # is directed at target if target can be pulled to surface.
            dist_Srf_to_target = pt_ProjHit.DistanceTo(other=pt_ProjdFrom)
            vect_Projection_to_target = pt_ProjdFrom - pt_ProjHit
            # Determine direction for w in uvw.
            angle = rg.Vector3d.VectorAngle(
                    normal_at_T, vect_Projection_to_target)
            if angle > PI / 4.0:
                dist_Srf_to_target = -dist_Srf_to_target
            dists_z.append((dist_Srf_to_target))
        return dists_z
    dists_z = get_ws_TargetsFromProjected(
            ns_WIP, pts_Target_with_proj_hit, pts_Projd, uvs_Target_projd_to_ns)
    if bDebug: sEval = 'dists_z[:20]'; print '{}:{}'.format(sEval,roundedList(eval(sEval)))


    # Check deviation before morphing.
    deviations_PreMorph_ThisIter_Max = max(abs(z) for z in dists_z)
    if deviations_PreMorph_ThisIter_Max <= fDistTol:
        if bDebug: print "Target is already within tolerance."
        return ns_WIP, deviations_PreMorph_ThisIter_Max
    
    #
    def getGrevillePointParameters(ns):
        """Create flat list of Greville UV's."""
        uvs_Greville = []
        for iU in xrange(ns.Points.CountU):
            for iV in xrange(ns.Points.CountV):
                uv_Greville = ns.Points.GetGrevillePoint(iU,iV)
                uvs_Greville.append(uv_Greville)
        return uvs_Greville
    uvs_Greville = getGrevillePointParameters(ns_WIP)
    if bDebug: sEval = 'uvs_Greville[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
    
    b_FreezeCps = [False]*ct_cp
    
    def get_uvDistTol_Upr(ns, fPercentOfDomain):
        domLen_U = abs(ns.Domain(0).Length)
        domLen_V = abs(ns.Domain(1).Length)
        uvDist_SrfHypotenuse = ( domLen_U**2.0 + domLen_V**2.0 )**0.5
        cp_span = 2
        pt_ct_U = ns.Points.CountU
        pt_ct_V = ns.Points.CountV
        dist_GrevilleSpan_U = domLen_U / ((pt_ct_U-1) / float(cp_span))
        dist_GrevilleSpan_V = domLen_V / ((pt_ct_V-1) / float(cp_span))
        
        
        uvDistTol_Upr = min(
                domLen_U*fPercentOfDomain/100.0,
                domLen_V*fPercentOfDomain/100.0)
        return uvDistTol_Upr
    uvDistTol_Upr = get_uvDistTol_Upr(ns_WIP, fPercentOfDomain)
    
    def get_uDistTol_Upr___uDistTol_Upr(ns):
        percent = 100.0
        rc = [abs(ns.Domain(d).Length) * percent / 100.0 for d in 0,1]
        return tuple(rc)
    u__vDistTol_Upr = get_uDistTol_Upr___uDistTol_Upr(ns_WIP)
    
    uvDistTol_Lwr = 0.01 * uvDistTol_Upr
    u__vDistTol_Lwr = 0.01 * u__vDistTol_Upr[0], 0.01 * u__vDistTol_Upr[1]
    
    # Corner points may be treated differently than other control points.
    def get_iG_Corners(ns):
        iG_Corners = (
                0,
                ns.Points.CountV - 1,
                (ns.Points.CountU - 1) * ns.Points.CountV,
                (ns.Points.CountU * ns.Points.CountV) - 1,
        )
        return iG_Corners
    iG_Corners = get_iG_Corners(ns_WIP)
    
    weightPerIter = 1.0
    weight_ForNonCornerCtrlPts = 2.0**0.5
    
    gSrfB = None
    
    if sWeightCurve == 'Linear':
        def getWeightPerDistance(distance, L, U):
            weight = 1.0 - (distance - L) / (U - L)
            return weight
    elif sWeightCurve == 'Hyperbolic':
        def getWeightPerDistance(distance, L, U):
            """At distance of U, weight = 0."""
            weight = 1.0 - (1.0 - L/distance) / (1.0 - L/U)
            return weight
    else:
        print "Weight curve, {}, not supported.".format(sWeightCurve)
        return
    
    if sSamplingMode == 'Distributed':
        def get_weightedDists(dists_params_TsInGs, dists_z, uvDistTol_Lwr, uvDistTol_Upr):
            weights_TsPerGs = []
            dists_Weighted_TsPerGs = []
            for iG, dists_uv_TsInG in enumerate(dists_params_TsInGs):
                weights_TsPerGs.append([])
                dists_Weighted_TsPerGs.append([])
                zipped = zip(dists_uv_TsInG, dists_z)
                for iT, (dist_uv_TargetToGrev, distW) in enumerate(zipped):
                    if 1.001*dist_uv_TargetToGrev >= uvDistTol_Upr:
                        # 1.001 multiple will accept more borderline cases.
                        weight_PerDist = 0.0
                    elif dist_uv_TargetToGrev <= uvDistTol_Lwr:
                        weight_PerDist = 1.0
                    else:
                        weight_PerDist = getWeightPerDistance(
                                distance=dist_uv_TargetToGrev,
                                L=uvDistTol_Lwr,
                                U=uvDistTol_Upr)
                    # Notice that W distance will affect the target choice below.
                    weights_TsPerGs[-1].append(weight_PerDist)
                    dists_Weighted_TsPerGs[-1].append(weight_PerDist*distW)
            dists_Weighted = []
            zipped = zip(dists_Weighted_TsPerGs, weights_TsPerGs)
            for iG, (dists_Weighted_TsPerG, weights_TsPerG) in enumerate(zipped):
                if weights_TsPerG.count(1.0) == 1:
                    # When there is a single, 'perfect' target pull to
                    # Greville point match, ignore other targets
                    # for this Greville.
                    dist_Weighted_Normalized = dists_Weighted_TsPerG[
                            weights_TsPerG.index(1.0)]
                else:
                    sum_weights_TsPerG = sum(weights_TsPerG)
                    if sum_weights_TsPerG == 0.0:
                        dist_Weighted_Normalized = 0.0
                    else:
                        dist_Weighted_Normalized = (
                                sum(dists_Weighted_TsPerG) / sum_weights_TsPerG)
                dists_Weighted.append(dist_Weighted_Normalized)
            return dists_Weighted
        
        def WIP_get_weightsPerDists():
            weightsPerDist_u__v_TsPerGs = [] # [G0[T0[u,v],T1[u,v],...],G1[T0[u,v],T1[u,v],...],...]
            for iG, dists_u__v_TsInG in enumerate(dists_u__v_TsInGs):
                weightsPerDist_u__v_TsPerGs.append([])
                for iT, dists_u__v_TInG in enumerate(dists_u__v_TsInG):
                    weightsPerDist_u__v_TsPerGs[-1].append([])
                    for iDir in 0,1:
                        uOrvDistTol_Lwr = u__vDistTol_Lwr[iDir]
                        uOrvDistTol_Upr = u__vDistTol_Upr[iDir]
                        dist_uOrv_TargetToGrev = dists_u__v_TInG[iDir]
                        if dist_uOrv_TargetToGrev >= uOrvDistTol_Upr:
                            weight_PerDist = 0.0
                        elif dist_uOrv_TargetToGrev <= uOrvDistTol_Lwr:
                            weight_PerDist = 1.0
                        else:
                            weight_PerDist = getWeightPerDistance(
                                    distance=dist_uOrv_TargetToGrev,
                                    L=uOrvDistTol_Lwr,
                                    U=uOrvDistTol_Upr)
                            weight_PerDist = round(
                                    number=weight_PerDist, ndigits=4)
                        #                            if weight_PerDist == 0.0:
                        #                                factor = 0
                        #                            elif iG in iG_Corners:
                        #                                # Don't multiply factors of corner points with
                        #                                # index+1 since points lie on the surface.
                        #                                factor = weight_PerDist * distW
                        #                            else:
                        #                                factor = weight_PerDist * distW * float(iMorph+1)
                        weightsPerDist_u__v_TsPerGs[-1][-1].append(weight_PerDist)
            print weightsPerDist_u__v_TsPerGs
            magnitudes = []
            for iG, weightsPerDist_u__v_TsPerG in enumerate(weightsPerDist_u__v_TsPerGs):
                magnitudes_TsPerG = []
                for iT, weightsPerDist_u__v_TPerG in enumerate(weightsPerDist_u__v_TsPerG):
                    distW = uvws_Target[iT][2]
                    
                    # Is min better than max or average?
                    weight_PerDist = min(weightsPerDist_u__v_TPerG)
                    
                    if weight_PerDist == 0.0:
                        mag = 0
                    
                    magnitudes_TsPerG.append(mag)
                magnitudes.append(sum(magnitudes_TsPerG) / len(magnitudes_TsPerG))
            return magnitudes
    elif sSamplingMode == 'NearestPoint':
        def get_weightedDists(dists_params_TsInGs, dists_z, uvDistTol_Lwr, uvDistTol_Upr):
            """Use each Greville closest target to calculate the Greville's control point translation magnitude."""
            dists_Weighted = []
            for iG, dists_params_TsInG in enumerate(dists_params_TsInGs):
                # weight_PerDist == 1.0 when uv of target is at uv of Greville.
                dist_NearestT = min(dists_params_TsInG)
                if 1.0001*dist_NearestT >= uvDistTol_Upr:
                    # 1.0001 multiplier will accept more borderline cases.
                    weight_PerDist = 0.0
                    dist_Weighted = 0.0
                elif dist_NearestT <= uvDistTol_Lwr:
                    weight_PerDist = 1.0
                else:
                    weight_PerDist = getWeightPerDistance(
                            distance=dist_NearestT,
                            L=uvDistTol_Lwr,
                            U=uvDistTol_Upr)
                    weight_PerDist = round(weight_PerDist, 6)
                    if weight_PerDist <= 0.0:
                        weight_PerDist = 0.0
                        dist_Weighted = 0.0
                    if weight_PerDist >= 1.0:
                        weight_PerDist = 1.0
                if weight_PerDist > 0.0:
                    dists_Weighted_Min = []
                    for iT, dist_params_TInG in enumerate(dists_params_TsInG):
                        if round(dist_params_TInG, 6) == round(dist_NearestT, 6):
                            dists_Weighted_Min.append(weight_PerDist * dists_z[iT])
                    dist_Weighted = sum(dists_Weighted_Min) / float(len(dists_Weighted_Min))
                
                # Limit magnitude to prevent wild deformations.
                #        if magnitude > 1000.0 * uvDistTol_Lwr:
                #            magnitude = 1000.0 * uvDistTol_Lwr
                
                dists_Weighted.append(dist_Weighted)
            return dists_Weighted
    else:
        print "Sampling mode, {}, not supported.".format(sSamplingMode)
        return
    
    
    def get_pts_cp(ns):
        # Create flat list of control points' points.
        pts_cp_ns0 = []
        for iU in xrange(ns.Points.CountU):
            for iV in xrange(ns.Points.CountV):
                pt_cp_ns0 = ns.Points.GetControlPoint(iU,iV).Location
                pts_cp_ns0.append(pt_cp_ns0)
        return pts_cp_ns0
    
    
    #    def get_dists_separate_u_and_v_components_TsInGs(uvs_Greville, uvs_Target_projd_to_ns):
    #        # Compile all UV distances of projected 
    #        dists_TsInGs = [] # Nested list of target tuples per Greville: [G0[T0,T1,...],G1[T0,T1,...],...]
    #        for uv_Greville in uvs_Greville:
    #            dists_TsInGs.append([])
    #            for uv_Target in uvs_Target_projd_to_ns:
    #                dists_TsInGs[-1].append((
    #                        abs(uv_Greville[0] - uv_Target[0]),
    #                        abs(uv_Greville[1] - uv_Target[1])
    #                ))
    #        return dists_TsInGs
    #    
    #    
    #    def get_dists_param_TsInGs(dists_u__v_TsInGs):
    #        # Compile all UV distances of projected 
    #        dists_params_TsInGs = [] # Nested list of target per Greville: [G0[T0,T1,...],G1[T0,T1,...],...]
    #        for iG, dists_u__v_TsInG in enumerate(dists_u__v_TsInGs):
    #            dists_params_TsInGs.append([])
    #            for dist_u, dist_v in dists_u__v_TsInG:
    #                dist_uv_TargetToGrev = ((dist_u)**2.0 + (dist_v)**2.0)**0.5
    #                dists_params_TsInGs[-1].append(dist_uv_TargetToGrev)
    #        return dists_params_TsInGs
    
    
    def get_dists_param_TsInGs(uvs_Greville, uvs_Target_projd_to_ns):
        # Compile all UV distances of projected 
        dists_params_TsInGs = [] # Nested list of target tuples per Greville: [G0[T0,T1,...],G1[T0,T1,...],...]
        for uv_Greville in uvs_Greville:
            dists_params_TsInGs.append([])
            for uv_Target in uvs_Target_projd_to_ns:
                dist = (
                        ((abs(uv_Greville[0] - uv_Target[0]))**2.0 +
                        (abs(uv_Greville[1] - uv_Target[1]))**2.0)**0.5
                )
                #dist = round(dist, 9)
                dists_params_TsInGs[-1].append(dist)
        return dists_params_TsInGs
    
    
    def get_magnitudes(dists_Weighted, weight_ForNonCornerCtrlPts, weightPerIter):
        magnitudes = []
        for iG, dist_Weighted in enumerate(dists_Weighted):
            mag = dist_Weighted
            #                if iG in iG_Corners:
            #                    # Don't multiply factors of corner points, which lie on surface,
            #                    # with any factors besides W, the distance from surface.
            #                    mag = dist_Weighted
            #                else:
            #                    mag = dist_Weighted * weight_ForNonCornerCtrlPts
            #                    #mag = dist_Weighted * weight_ForNonCornerCtrlPts * weightPerIter
            magnitudes.append(mag)
        return magnitudes
    
    
    def createZVectorsPerNormalsAtGrevillePoints(ns):
        """Create flat list of normal unit vectors at Grevilles."""
        rads_90Degs = Rhino.RhinoMath.ToRadians(90.0)
        vects_GrevNormal = []
        for iU in xrange(ns.Points.CountU):
            for iV in xrange(ns.Points.CountV):
                vect_GrevNormal = ns.NormalAt(*tuple(ns.Points.GetGrevillePoint(iU,iV)))
                if (
                        rg.Vector3d.VectorAngle(vect_GrevNormal, rg.Vector3d.ZAxis) <
                        rads_90Degs
                ):
                    vects_GrevNormal.append(rg.Vector3d.ZAxis)
                else:
                    vects_GrevNormal.append(-rg.Vector3d.ZAxis)
        return vects_GrevNormal
    
    
    def set_b_FreezeCps(b_FreezeCps, uvs_Greville, uvs_Target_projd_to_ns, dists_z):
        for iG, uv_Greville in enumerate(uvs_Greville):
            if b_FreezeCps[iG]:
                continue # to next iG.
            for iT, uv_Target in enumerate(uvs_Target_projd_to_ns):
                if (not Rhino.RhinoMath.EpsilonEquals(
                    uv_Greville[0], uv_Target[0], epsilon=1e-6)
                ):
                    continue
                if (not Rhino.RhinoMath.EpsilonEquals(
                    uv_Greville[1], uv_Target[1], epsilon=1e-6)
                ):
                    continue
                if (not Rhino.RhinoMath.EpsilonEquals(
                    dists_z[iT], 0.0, epsilon=1e-6)
                ):
                    continue
                b_FreezeCps[iG] = True
                break # to next iG.
    
    
    def get_translations(vects, magnitudes, b_FreezeCps):
        translations = []
        for iGrev, (v, m) in enumerate(zip(vects, magnitudes)):
            if b_FreezeCps[iGrev]:
                translations.append(rg.Vector3d.Zero)
            else:
                translations.append(v*m)
        return translations
    
    
    def get_translatedPoints(pts_cp_ns0, translations):
        pts1 = []
        for iG, (pt0, t) in enumerate(zip(pts_cp_ns0, translations)):
            pt1 = pt0 + t
            #            pt1.X = round(pt1.X, 6)
            #            pt1.Y = round(pt1.Y, 6)
            #            pt1.Z = round(pt1.Z, 6)
            pts1.append(pt1)
        return pts1
    
    
    def setNurbsSurfacePointLocations(ns, pts1):
        for iV in xrange(ns.Points.CountV):
            for iU in xrange(ns.Points.CountU):
                cp_ns1 = ns.Points.GetControlPoint(u=iU, v=iV)
                cp_ns1.Location = pts1[iU*ns.Points.CountV+iV]
                ns.Points.SetControlPoint(u=iU, v=iV, cp=cp_ns1)


    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt


    if bShowIterations:
        conduit = DrawSurfaceIterationDisplayConduit(rgNurbsSrf, pts_Target)


    #
    # Main morphing iterations.
    for iMorph in xrange(100):
        if sc.escape_test(throw_exception=False):
            print "Break in main morphing iterations ..."
            break
        
        if bDebug: print "{}\nIterationIndex:{}".format('#'*16, iMorph)
        
        pts_cp_ns0 = get_pts_cp(ns_WIP)
        if bDebug: sEval = 'pts_cp_ns0[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        #        dists_u__v_TsInGs = get_dists_separate_u_and_v_components_TsInGs(
        #                uvs_Greville, uvs_Target_projd_to_ns) # Tuples of u,v.
        
        #        dists_params_TsInGs = get_dists_param_TsInGs(dists_u__v_TsInGs)
        dists_params_TsInGs = get_dists_param_TsInGs(uvs_Greville, uvs_Target_projd_to_ns)
        if bDebug: sEval = 'dists_params_TsInGs[0][{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        dists_Weighted = get_weightedDists(dists_params_TsInGs, dists_z, uvDistTol_Lwr, uvDistTol_Upr)
        if bDebug: sEval = 'dists_Weighted[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        magnitudes = get_magnitudes(dists_Weighted, weight_ForNonCornerCtrlPts, weightPerIter)
        if bDebug: sEval = 'magnitudes[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        vects_atGrev = createZVectorsPerNormalsAtGrevillePoints(ns_WIP)
        if bDebug: sEval = 'vects_atGrev[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        set_b_FreezeCps(b_FreezeCps, uvs_Greville, uvs_Target_projd_to_ns, dists_z)
        if bDebug: sEval = 'b_FreezeCps'; print '{}:{}'.format(sEval,eval(sEval))
        
        translations = get_translations(vects_atGrev, magnitudes, b_FreezeCps)
        if bDebug: sEval = 'translations[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        pts_cv_ns1 = get_translatedPoints(pts_cp_ns0, translations)
        if bDebug: sEval = 'pts_cv_ns1[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        
        ns_Morphed = ns_WIP.Duplicate()
        setNurbsSurfacePointLocations(ns_Morphed, pts_cv_ns1)
        
        if bDebug: sc.doc.Views.Redraw()
        
        transLengths_ThisIter = [translation.Length for translation in translations]
        if bDebug: sEval = 'transLengths_ThisIter[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,eval(sEval))
        
        b_Pulled_Successes_Prev = b_Projd_successes
        
        b_Projd_successes, pts_Projd, pts_Target_with_proj_hit = create_pts_projd_to_srf(ns_Morphed, pts_Target)
        if not pts_Projd:
            print "No target points can be pulled to surface."
            break
        
        if b_Projd_successes == b_Pulled_Successes_Prev:
            weightPerIter += 1.0
        else:
            weightPerIter = 1.0
        if bDebug: sEval = 'weightPerIter'; print '{}:{}'.format(sEval,eval(sEval))
        
        uvs_Target_projd_to_ns = get_uvs_of_Point3ds(ns_Morphed, pts_Projd) # Only targets that pulled to surface.
        if bDebug: sEval = 'uvs_Target_projd_to_ns[:12]'; print '{}:{}'.format(sEval,eval(sEval))
        
        dists_z = get_ws_TargetsFromProjected(ns_Morphed, pts_Target_with_proj_hit, pts_Projd, uvs_Target_projd_to_ns)
        if bDebug: sEval = 'dists_z[:12]'; print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        nss_Morphed.append(ns_Morphed)
        deviations_PostMorph_ThisIter = [abs(w) for w in dists_z]
        deviations_PostMorph_ThisIter_Max = max(deviations_PostMorph_ThisIter)
        deviations_PostMorph_AllIters_Max.append(deviations_PostMorph_ThisIter_Max)
        transLength_ThisIter_Max = max(transLengths_ThisIter)
        transLengths_AllIters_Max.append(transLength_ThisIter_Max)
        #
        

        if bShowIterations:
            conduit.srf = ns_Morphed
            conduit.Enabled = True
            sc.doc.Views.Redraw()

        s  = "Iteration:{0}".format(iMorph+1)
        s += " MaxDeviation:{0:.{1}f}".format(deviations_PostMorph_ThisIter_Max, sc.doc.ModelDistanceDisplayPrecision+1)
        s += " MaxTranslation:{0:.{1}f}".format(transLength_ThisIter_Max, sc.doc.ModelDistanceDisplayPrecision+1)
        if sCmdPrompt0:
            s += "  "
        Rhino.RhinoApp.CommandPrompt = sCmdPrompt0 + s
        
        if bDebug: print s
        
        if deviations_PostMorph_ThisIter_Max <= fDistTol:
            if bDebug:
                sEvals = 'deviations_PostMorph_ThisIter_Max', 'fDistTol'; print 'Break from iterations due to {} <= {} ({} <= {})'.format(sEvals[0],sEvals[1],eval(sEvals[0]),eval(sEvals[1]))
            break
        
        if iMorph > 0:
            if deviations_PostMorph_ThisIter_Max >= min(deviations_PostMorph_AllIters_Max[:-1]):
                # Deviations are not decreasing.
                if bDebug:
                    sEvals = 'deviations_PostMorph_ThisIter_Max', 'min(deviations_PostMorph_AllIters_Max[:-1])'; print 'Break from iterations due to {} >= {} ({} >= {})'.format(sEvals[0],sEvals[1],eval(sEvals[0]),eval(sEvals[1]))
                break
            #        elif transLength_ThisIter_Max <= (0.1 * fDistTol):
            #            break
            elif (
                    (abs(
                            deviations_PostMorph_AllIters_Max[-2] -
                            deviations_PostMorph_AllIters_Max[-1]) <
                    0.1 * fDistTol)
            ):
                # Difference in maximum deviation of last 2 iterations is less than
                # 0.1*fDistTol,
                # so don't bother continuing morphing the current starting surface.
                #print "  Deviation difference: {}".format(abs(deviations_PostMorph_AllIters_Max[-2] - deviations_PostMorph_AllIters_Max[1]))
                if bDebug: print "Break from iterations due to (abs(deviations_PostMorph_AllIters_Max[-2] - deviations_PostMorph_AllIters_Max[-1]) < (0.1*fDistTol))."
                break
        
        ns_WIP = ns_Morphed
        ns_Morphed = None
        
        if deviations_PostMorph_ThisIter_Max is None:
            conduit.Enabled = False
            return
    
    conduit.Enabled = False

    if not nss_Morphed:
        return
    else:
        if bDebug:
            print '\n' + '#'*40
        
        deviation_MinOfMax = min(deviations_PostMorph_AllIters_Max)
        idx_deviation_MinOfMax = deviations_PostMorph_AllIters_Max.index(deviation_MinOfMax)
        
        transLength_MinOfMax = min(transLengths_AllIters_Max)
        idx_transLength_MinOfMax = transLengths_AllIters_Max.index(transLength_MinOfMax)
        
        if bDebug:
            s  = "After {} iteration(s),".format(iMorph+1 if ns_Morphed else iMorph)
            s += "  MaximumDeviationFromTargetPoints:{0:.{1}f}".format(
                    deviation_MinOfMax,
                    sc.doc.ModelDistanceDisplayPrecision+2)
            s += " @Iteration:{}".format(idx_deviation_MinOfMax+1)
            s += "  MaximumControlPointTranslation:{0:.{1}f}".format(
                    transLength_MinOfMax,
                    sc.doc.ModelDistanceDisplayPrecision+2)
            s += " @Iteration:{}".format(idx_transLength_MinOfMax+1)
            print s
        
        ns_MinDev = nss_Morphed[idx_deviation_MinOfMax]

        for i, ns in enumerate(nss_Morphed):
            if i == idx_deviation_MinOfMax: continue
            ns.Dispose()
        
        ns_MinDev.SetDomain(direction=0, domain=rgNurbsSrf.Domain(direction=0))
        ns_MinDev.SetDomain(direction=1, domain=rgNurbsSrf.Domain(direction=1))
        
        return ns_MinDev, deviation_MinOfMax


def fit_Pulled(rgNurbsSrf, pts_Target, bRecalcNormals=True, fDistTol=0.1*sc.doc.ModelAbsoluteTolerance, sSamplingMode='Distributed', sWeightCurve='Linear', fPercentOfDomain=50.0, bShowIterations=False, bDebug=False):
    """
    bRecalcNormals: True to obtain the normals for each surface iteration
                    False to use the normals of the starting surface only.
    """
    
    if not rgNurbsSrf.IsValid:
        print "Geometry of NurbsSurface to morph is invalid!"
        return
    
    ns_WIP = rgNurbsSrf.Duplicate() # Keep copy of original because ns_WIP will be modified per iteration.
    
    # Make domains of both directions [0,1].
    interval01 = rg.Interval(t0=0.0, t1=1.0)
    ns_WIP.SetDomain(direction=0, domain=interval01)
    ns_WIP.SetDomain(direction=1, domain=interval01)
    
    ct_cp = ns_WIP.Points.CountU * ns_WIP.Points.CountV
    
    nss_Morphed = []
    deviations_PostMorph_AllIters_Max = []
    transLengths_AllIters_Max = []
    
    if bDebug: print "Using {} sampling mode with {} weight curve.".format(sSamplingMode, sWeightCurve)


    def create_pts_pulled_to_srf(ns, pts0):
        """Try to pull all target points to surface."""
        brep = ns.ToBrep()
        face = brep.Faces[0]
        b_Pulled_Successes = []
        pts_Pulled = []
        pts0_Hits = []
        for pt0 in pts0:
            rc = face.PullPointsToFace(
                    points=[pt0],
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
            )
            if not rc:
                #
                # TODO: Test angle of closest point with that of the normal at UV.
                # For angles within ?, add ClosestPoints to pt_Pulled.
                #
                b_Pulled_Successes.append(False)
            else:
                b_Pulled_Successes.append(True)
                pt_Pulled = rc[0]
                pts_Pulled.append(pt_Pulled)
                pts0_Hits.append(pt0)
        brep.Dispose()
        return b_Pulled_Successes, pts_Pulled, pts0_Hits
    rc = create_pts_pulled_to_srf(ns_WIP, pts_Target)
    b_Pulled_Successes, pts_Target_Pulled, pts_Target_Hits = rc
    
    if not pts_Target_Pulled:
        print "No target points can be pulled to surface."
        return
    if bDebug: print "{} target points were pulled to surface.".format(len(pts_Target_Pulled))
    
    
    def roundedList(i_list):
        o_list = []
        for i in i_list:
            if isinstance(i, (tuple, list, rg.Point2d, rg.Point3d, rg.Vector3d)):
                o_list.append([])
                for j in i:
                    o_list[-1].append(str(round(j,12)))
            else:
                o_list.append(str(round(i,12)))
        return o_list
    
    
    def get_uvs_of_Point3ds(srf, pts):
        """Apply Surface.ClosestPoint to a list of Point3ds."""
        uvs_Target_projd_to_ns = []
        for pt in pts:
            rc, u, v = srf.ClosestPoint(pt)
            if not rc:
                print "ClosestPoint returned None"
                return
            #            if Rhino.RhinoMath.EpsilonEquals(u, 0.0, 1e-12): u = 0.0
            #            if Rhino.RhinoMath.EpsilonEquals(v, 0.0, 1e-12): v = 0.0
            uvs_Target_projd_to_ns.append((u, v))
        return uvs_Target_projd_to_ns
    uvs_Target_projd_to_ns = get_uvs_of_Point3ds(
            ns_WIP, pts_Target_Pulled) # Only targets that pulled to surface.
    if bDebug: sEval = 'uvs_Target_projd_to_ns[:12]'; print '{}:{}'.format(sEval,roundedList(eval(sEval)))
    
    def get_ws_TargetsFromPulled(ns, pts_ToPull, pts_Pulled, uvs):
        """
        For every target point that can be pulled to NurbsSurface,
        w, the signed Euclidian distance to surface, will be generated.
        """
        dists_w = []
        for pt_ToPull, pt_Pulled, uv in zip(pts_ToPull, pts_Pulled, uvs):
            normal_at_T = ns.NormalAt(*uv) # is directed at target if target can be pulled to surface.
            dist_Srf_to_target = pt_Pulled.DistanceTo(other=pt_ToPull)
            vect_Projection_to_target = pt_ToPull - pt_Pulled
            # Determine direction for w in uvw.
            angle = rg.Vector3d.VectorAngle(
                    normal_at_T, vect_Projection_to_target)
            if angle > PI / 4.0:
                dist_Srf_to_target = -dist_Srf_to_target
            dists_w.append((dist_Srf_to_target))
        return dists_w
    dists_w = get_ws_TargetsFromPulled(
            ns_WIP, pts_Target_Hits, pts_Target_Pulled, uvs_Target_projd_to_ns)
    if bDebug: sEval = 'dists_w[:20]'; print '{}:{}'.format(sEval,roundedList(eval(sEval)))


    # Check deviation before morphing.
    deviations_PreMorph_ThisIter_Max = max(abs(w) for w in dists_w)
    if deviations_PreMorph_ThisIter_Max <= fDistTol:
        if bDebug: print "Target is already within tolerance."
        return ns_WIP, deviations_PreMorph_ThisIter_Max
    
    #
    def getGrevillePointParameters(ns):
        """Create flat list of Greville UV's."""
        uvs_Greville = []
        for iU in xrange(ns.Points.CountU):
            for iV in xrange(ns.Points.CountV):
                uv_Greville = ns.Points.GetGrevillePoint(iU,iV)
                uvs_Greville.append(uv_Greville)
        return uvs_Greville
    uvs_Greville = getGrevillePointParameters(ns_WIP)
    if bDebug: sEval = 'uvs_Greville[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
    
    b_FreezeCps = [False]*ct_cp
    
    def get_uvDistTol_Upr(ns, fPercentOfDomain):
        domLen_U = abs(ns.Domain(0).Length)
        domLen_V = abs(ns.Domain(1).Length)
        uvDist_SrfHypotenuse = ( domLen_U**2.0 + domLen_V**2.0 )**0.5
        cp_span = 2
        pt_ct_U = ns.Points.CountU
        pt_ct_V = ns.Points.CountV
        dist_GrevilleSpan_U = domLen_U / ((pt_ct_U-1) / float(cp_span))
        dist_GrevilleSpan_V = domLen_V / ((pt_ct_V-1) / float(cp_span))
        
        
        uvDistTol_Upr = min(
                domLen_U*fPercentOfDomain/100.0,
                domLen_V*fPercentOfDomain/100.0)
        #        if sWeightCurve == 'Hyperbolic':
        #            uvDistTol_Upr = uvDist_SrfHypotenuse
        #        elif sWeightCurve == 'Linear':
        #            # max, min, average?
        #            iUpperLimit = 2
        #            if iUpperLimit == 0:
        #                uvDistTol_Upr = min(dist_GrevilleSpan_U, dist_GrevilleSpan_V)
        #            elif iUpperLimit == 1:
        #                uvDistTol_Upr = max(dist_GrevilleSpan_U, dist_GrevilleSpan_V)
        #            else:
        #                uvDistTol_Upr = (dist_GrevilleSpan_U + dist_GrevilleSpan_V) / 2.0
        #        else:
        #            print "Weight curve, {}, not supported.".format(sWeightCurve)
        #            return
        return uvDistTol_Upr
    uvDistTol_Upr = get_uvDistTol_Upr(ns_WIP, fPercentOfDomain)
    
    def get_uDistTol_Upr___uDistTol_Upr(ns):
        percent = 100.0
        rc = [abs(ns.Domain(d).Length) * percent / 100.0 for d in 0,1]
        return tuple(rc)
    u__vDistTol_Upr = get_uDistTol_Upr___uDistTol_Upr(ns_WIP)
    
    uvDistTol_Lwr = 0.01 * uvDistTol_Upr
    u__vDistTol_Lwr = 0.01 * u__vDistTol_Upr[0], 0.01 * u__vDistTol_Upr[1]
    
    # Corner points may be treated differently than other control points.
    def get_iG_Corners(ns):
        iG_Corners = (
                0,
                ns.Points.CountV - 1,
                (ns.Points.CountU - 1) * ns.Points.CountV,
                (ns.Points.CountU * ns.Points.CountV) - 1,
        )
        return iG_Corners
    iG_Corners = get_iG_Corners(ns_WIP)
    
    weightPerIter = 1.0
    weight_ForNonCornerCtrlPts = 2.0**0.5
    
    gSrfB = None
    
    if sWeightCurve == 'Linear':
        def getWeightPerDistance(distance, L, U):
            weight = 1.0 - (distance - L) / (U - L)
            return weight
    elif sWeightCurve == 'Hyperbolic':
        def ALT_getWeightPerDistance(distance, L, U):
            """At distance of U, weight = L."""
            weight = ((1.0 - L/distance) / (1.0 - L/U)) * (L - 1.0) + 1.0
            return weight
        
        
        def getWeightPerDistance(distance, L, U):
            """At distance of U, weight = 0."""
            weight = 1.0 - (1.0 - L/distance) / (1.0 - L/U)
            return weight
    else:
        print "Weight curve, {}, not supported.".format(sWeightCurve)
        return
    
    if sSamplingMode == 'Distributed':
        def get_weightedDists(dists_params_TsInGs, dists_w, uvDistTol_Lwr, uvDistTol_Upr):
            weights_TsPerGs = []
            dists_Weighted_TsPerGs = []
            for iG, dists_uv_TsInG in enumerate(dists_params_TsInGs):
                weights_TsPerGs.append([])
                dists_Weighted_TsPerGs.append([])
                zipped = zip(dists_uv_TsInG, dists_w)
                for iT, (dist_uv_TargetToGrev, distW) in enumerate(zipped):
                    if 1.001*dist_uv_TargetToGrev >= uvDistTol_Upr:
                        # 1.001 multiple will accept more borderline cases.
                        weight_PerDist = 0.0
                    elif dist_uv_TargetToGrev <= uvDistTol_Lwr:
                        weight_PerDist = 1.0
                    else:
                        weight_PerDist = getWeightPerDistance(
                                distance=dist_uv_TargetToGrev,
                                L=uvDistTol_Lwr,
                                U=uvDistTol_Upr)
                    # Notice that W distance will affect the target choice below.
                    weights_TsPerGs[-1].append(weight_PerDist)
                    dists_Weighted_TsPerGs[-1].append(weight_PerDist*distW)
            dists_Weighted = []
            zipped = zip(dists_Weighted_TsPerGs, weights_TsPerGs)
            for iG, (dists_Weighted_TsPerG, weights_TsPerG) in enumerate(zipped):
                if weights_TsPerG.count(1.0) == 1:
                    # When there is a single, 'perfect' target pull to
                    # Greville point match, ignore other targets
                    # for this Greville.
                    dist_Weighted_Normalized = dists_Weighted_TsPerG[
                            weights_TsPerG.index(1.0)]
                else:
                    sum_weights_TsPerG = sum(weights_TsPerG)
                    if sum_weights_TsPerG == 0.0:
                        dist_Weighted_Normalized = 0.0
                    else:
                        dist_Weighted_Normalized = (
                                sum(dists_Weighted_TsPerG) / sum_weights_TsPerG)
                dists_Weighted.append(dist_Weighted_Normalized)
            return dists_Weighted
        
        def WIP_get_weightsPerDists():
            weightsPerDist_u__v_TsPerGs = [] # [G0[T0[u,v],T1[u,v],...],G1[T0[u,v],T1[u,v],...],...]
            for iG, dists_u__v_TsInG in enumerate(dists_u__v_TsInGs):
                weightsPerDist_u__v_TsPerGs.append([])
                for iT, dists_u__v_TInG in enumerate(dists_u__v_TsInG):
                    weightsPerDist_u__v_TsPerGs[-1].append([])
                    for iDir in 0,1:
                        uOrvDistTol_Lwr = u__vDistTol_Lwr[iDir]
                        uOrvDistTol_Upr = u__vDistTol_Upr[iDir]
                        dist_uOrv_TargetToGrev = dists_u__v_TInG[iDir]
                        if dist_uOrv_TargetToGrev >= uOrvDistTol_Upr:
                            weight_PerDist = 0.0
                        elif dist_uOrv_TargetToGrev <= uOrvDistTol_Lwr:
                            weight_PerDist = 1.0
                        else:
                            weight_PerDist = getWeightPerDistance(
                                    distance=dist_uOrv_TargetToGrev,
                                    L=uOrvDistTol_Lwr,
                                    U=uOrvDistTol_Upr)
                            weight_PerDist = round(
                                    number=weight_PerDist, ndigits=4)
                        #                            if weight_PerDist == 0.0:
                        #                                factor = 0
                        #                            elif iG in iG_Corners:
                        #                                # Don't multiply factors of corner points with
                        #                                # index+1 since points lie on the surface.
                        #                                factor = weight_PerDist * distW
                        #                            else:
                        #                                factor = weight_PerDist * distW * float(iMorph+1)
                        weightsPerDist_u__v_TsPerGs[-1][-1].append(weight_PerDist)
            print weightsPerDist_u__v_TsPerGs
            magnitudes = []
            for iG, weightsPerDist_u__v_TsPerG in enumerate(weightsPerDist_u__v_TsPerGs):
                magnitudes_TsPerG = []
                for iT, weightsPerDist_u__v_TPerG in enumerate(weightsPerDist_u__v_TsPerG):
                    distW = uvws_Target[iT][2]
                    
                    # Is min better than max or average?
                    weight_PerDist = min(weightsPerDist_u__v_TPerG)
                    
                    if weight_PerDist == 0.0:
                        mag = 0
                    
                    magnitudes_TsPerG.append(mag)
                magnitudes.append(sum(magnitudes_TsPerG) / len(magnitudes_TsPerG))
            return magnitudes
    elif sSamplingMode == 'NearestPoint':
        def get_weightedDists(dists_params_TsInGs, dists_w, uvDistTol_Lwr, uvDistTol_Upr):
            """Use each Greville closest target to calculate the Greville's control point translation magnitude."""
            dists_Weighted = []
            for iG, dists_params_TsInG in enumerate(dists_params_TsInGs):
                # weight_PerDist == 1.0 when uv of target is at uv of Greville.
                dist_NearestT = min(dists_params_TsInG)
                if 1.0001*dist_NearestT >= uvDistTol_Upr:
                    # 1.0001 multiplier will accept more borderline cases.
                    weight_PerDist = 0.0
                    dist_Weighted = 0.0
                elif dist_NearestT <= uvDistTol_Lwr:
                    weight_PerDist = 1.0
                else:
                    weight_PerDist = getWeightPerDistance(
                            distance=dist_NearestT,
                            L=uvDistTol_Lwr,
                            U=uvDistTol_Upr)
                    weight_PerDist = round(weight_PerDist, 6)
                    if weight_PerDist <= 0.0:
                        weight_PerDist = 0.0
                        dist_Weighted = 0.0
                    if weight_PerDist >= 1.0:
                        weight_PerDist = 1.0
                if weight_PerDist > 0.0:
                    dists_Weighted_Min = []
                    for iT, dist_params_TInG in enumerate(dists_params_TsInG):
                        if round(dist_params_TInG, 6) == round(dist_NearestT, 6):
                            dists_Weighted_Min.append(weight_PerDist * dists_w[iT])
                    dist_Weighted = sum(dists_Weighted_Min) / float(len(dists_Weighted_Min))
                
                # Limit magnitude to prevent wild deformations.
                #        if magnitude > 1000.0 * uvDistTol_Lwr:
                #            magnitude = 1000.0 * uvDistTol_Lwr
                
                dists_Weighted.append(dist_Weighted)
            return dists_Weighted
    else:
        print "Sampling mode, {}, not supported.".format(sSamplingMode)
        return
    
    
    def get_pts_cp(ns):
        # Create flat list of control points' points.
        pts_cp_ns0 = []
        for iU in xrange(ns.Points.CountU):
            for iV in xrange(ns.Points.CountV):
                pt_cp_ns0 = ns.Points.GetControlPoint(iU,iV).Location
                pts_cp_ns0.append(pt_cp_ns0)
        return pts_cp_ns0
    
    
    #    def get_dists_separate_u_and_v_components_TsInGs(uvs_Greville, uvs_Target_projd_to_ns):
    #        # Compile all UV distances of projected 
    #        dists_TsInGs = [] # Nested list of target tuples per Greville: [G0[T0,T1,...],G1[T0,T1,...],...]
    #        for uv_Greville in uvs_Greville:
    #            dists_TsInGs.append([])
    #            for uv_Target in uvs_Target_projd_to_ns:
    #                dists_TsInGs[-1].append((
    #                        abs(uv_Greville[0] - uv_Target[0]),
    #                        abs(uv_Greville[1] - uv_Target[1])
    #                ))
    #        return dists_TsInGs
    #    
    #    
    #    def get_dists_param_TsInGs(dists_u__v_TsInGs):
    #        # Compile all UV distances of projected 
    #        dists_params_TsInGs = [] # Nested list of target per Greville: [G0[T0,T1,...],G1[T0,T1,...],...]
    #        for iG, dists_u__v_TsInG in enumerate(dists_u__v_TsInGs):
    #            dists_params_TsInGs.append([])
    #            for dist_u, dist_v in dists_u__v_TsInG:
    #                dist_uv_TargetToGrev = ((dist_u)**2.0 + (dist_v)**2.0)**0.5
    #                dists_params_TsInGs[-1].append(dist_uv_TargetToGrev)
    #        return dists_params_TsInGs
    
    
    def get_dists_param_TsInGs(uvs_Greville, uvs_Target_projd_to_ns):
        # Compile all UV distances of projected 
        dists_params_TsInGs = [] # Nested list of target tuples per Greville: [G0[T0,T1,...],G1[T0,T1,...],...]
        for uv_Greville in uvs_Greville:
            dists_params_TsInGs.append([])
            for uv_Target in uvs_Target_projd_to_ns:
                dist = (
                        ((abs(uv_Greville[0] - uv_Target[0]))**2.0 +
                        (abs(uv_Greville[1] - uv_Target[1]))**2.0)**0.5
                )
                #dist = round(dist, 9)
                dists_params_TsInGs[-1].append(dist)
        return dists_params_TsInGs
    
    
    def get_magnitudes(dists_Weighted, weight_ForNonCornerCtrlPts, weightPerIter):
        magnitudes = []
        for iG, dist_Weighted in enumerate(dists_Weighted):
            mag = dist_Weighted
            #                if iG in iG_Corners:
            #                    # Don't multiply factors of corner points, which lie on surface,
            #                    # with any factors besides W, the distance from surface.
            #                    mag = dist_Weighted
            #                else:
            #                    mag = dist_Weighted * weight_ForNonCornerCtrlPts
            #                    #mag = dist_Weighted * weight_ForNonCornerCtrlPts * weightPerIter
            magnitudes.append(mag)
        return magnitudes
    
    
    def getNormalsAtGrevillePoints(ns):
        """Create flat list of normal unit vectors at Grevilles."""
        vects_GrevNormal = []
        for iU in xrange(ns.Points.CountU):
            for iV in xrange(ns.Points.CountV):
                vect_GrevNormal = ns.NormalAt(*tuple(ns.Points.GetGrevillePoint(iU,iV)))
                vects_GrevNormal.append(vect_GrevNormal)
        return vects_GrevNormal
    
    
    def set_b_FreezeCps(b_FreezeCps, uvs_Greville, uvs_Target_projd_to_ns, dists_w):
        for iG, uv_Greville in enumerate(uvs_Greville):
            if b_FreezeCps[iG]:
                continue # to next iG.
            for iT, uv_Target in enumerate(uvs_Target_projd_to_ns):
                if (not Rhino.RhinoMath.EpsilonEquals(
                    uv_Greville[0], uv_Target[0], epsilon=1e-6)
                ):
                    continue
                if (not Rhino.RhinoMath.EpsilonEquals(
                    uv_Greville[1], uv_Target[1], epsilon=1e-6)
                ):
                    continue
                if (not Rhino.RhinoMath.EpsilonEquals(
                    dists_w[iT], 0.0, epsilon=1e-6)
                ):
                    continue
                b_FreezeCps[iG] = True
                break # to next iG.
    
    
    def get_translations(vects_GrevNormal, magnitudes, b_FreezeCps):
        translations = []
        for iG, (v, m) in enumerate(zip(vects_GrevNormal, magnitudes)):
            if b_FreezeCps[iG]:
                translations.append(rg.Vector3d.Zero)
            else:
                translations.append(v*m)
        return translations
    
    
    def get_translatedPoints(pts_cp_ns0, translations):
        pts1 = []
        for iG, (pt0, t) in enumerate(zip(pts_cp_ns0, translations)):
            pt1 = pt0 + t
            #            pt1.X = round(pt1.X, 6)
            #            pt1.Y = round(pt1.Y, 6)
            #            pt1.Z = round(pt1.Z, 6)
            pts1.append(pt1)
        return pts1
    
    
    def setNurbsSurfacePointLocations(ns, pts1):
        for iV in xrange(ns.Points.CountV):
            for iU in xrange(ns.Points.CountU):
                cp_ns1 = ns.Points.GetControlPoint(u=iU, v=iV)
                cp_ns1.Location = pts1[iU*ns.Points.CountV+iV]
                ns.Points.SetControlPoint(u=iU, v=iV, cp=cp_ns1)


    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt


    if bShowIterations:
        conduit = DrawSurfaceIterationDisplayConduit(rgNurbsSrf, pts_Target)


    #
    # Main morphing iterations.
    for iMorph in xrange(100):
        if sc.escape_test(throw_exception=False):
            print "Break in main morphing iterations ..."
            break
        
        if bDebug: print "{}\nIterationIndex:{}".format('#'*16, iMorph)
        
        pts_cp_ns0 = get_pts_cp(ns_WIP)
        if bDebug: sEval = 'pts_cp_ns0[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        #        dists_u__v_TsInGs = get_dists_separate_u_and_v_components_TsInGs(
        #                uvs_Greville, uvs_Target_projd_to_ns) # Tuples of u,v.
        
        #        dists_params_TsInGs = get_dists_param_TsInGs(dists_u__v_TsInGs)
        dists_params_TsInGs = get_dists_param_TsInGs(uvs_Greville, uvs_Target_projd_to_ns)
        if bDebug: sEval = 'dists_params_TsInGs[0][{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        dists_Weighted = get_weightedDists(dists_params_TsInGs, dists_w, uvDistTol_Lwr, uvDistTol_Upr)
        if bDebug: sEval = 'dists_Weighted[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        magnitudes = get_magnitudes(dists_Weighted, weight_ForNonCornerCtrlPts, weightPerIter)
        if bDebug: sEval = 'magnitudes[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        if bRecalcNormals or iMorph == 0:
            vects_GrevNormal = getNormalsAtGrevillePoints(ns_WIP)
            if bDebug: sEval = 'vects_GrevNormal[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        set_b_FreezeCps(b_FreezeCps, uvs_Greville, uvs_Target_projd_to_ns, dists_w)
        if bDebug: sEval = 'b_FreezeCps'; print '{}:{}'.format(sEval,eval(sEval))
        
        translations = get_translations(vects_GrevNormal, magnitudes, b_FreezeCps)
        if bDebug: sEval = 'translations[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        pts_cv_ns1 = get_translatedPoints(pts_cp_ns0, translations)
        if bDebug: sEval = 'pts_cv_ns1[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        #    map(sc.doc.Objects.AddPoint, pts_cv_ns1)
        #    sc.doc.Views.Redraw()
        
        # Check for symmetry.
        #        for iP in xrange(len(pts_cv_ns1)//2):
        #            if not (
        #                Rhino.RhinoMath.EpsilonEquals(
        #                        abs(pts_cv_ns1[iP].X),
        #                        abs(pts_cv_ns1[-1-iP].X),
        #                        1e-6)
        #                and
        #                Rhino.RhinoMath.EpsilonEquals(
        #                        abs(pts_cv_ns1[iP].Y),
        #                        abs(pts_cv_ns1[-1-iP].Y),
        #                        1e-6)
        #                and
        #                Rhino.RhinoMath.EpsilonEquals(
        #                        abs(pts_cv_ns1[iP].Z),
        #                        abs(pts_cv_ns1[-1-iP].Z),
        #                        1e-6)
        #            ):
        #                print "Not symmetric at point {} in iteration index {}!".format(iP, iMorph)
        #                print pts_cv_ns1[iP], pts_cv_ns1[-1-iP]
        #                sc.doc.Objects.AddSurface(ns_WIP)
        #                ns_Morphed = setNurbsSurfacePointLocations(ns_WIP.Duplicate(), pts_cv_ns1)
        #                sc.doc.Objects.AddSurface(ns_Morphed)
        #                1/0
        
        ns_Morphed = ns_WIP.Duplicate()
        setNurbsSurfacePointLocations(ns_Morphed, pts_cv_ns1)
        
        if bDebug: sc.doc.Views.Redraw()
        
        transLengths_ThisIter = [translation.Length for translation in translations]
        if bDebug: sEval = 'transLengths_ThisIter[{}]'.format(sDebugRange); print '{}:{}'.format(sEval,eval(sEval))
        
        b_Pulled_Successes_Prev = b_Pulled_Successes
        
        b_Pulled_Successes, pts_Target_Pulled, pts_Target_Hits = create_pts_pulled_to_srf(ns_Morphed, pts_Target)
        if not pts_Target_Pulled:
            print "No target points can be pulled to surface."
            break
        #if bDebug: sEval = 'map(Object.ToString, pts_Target_Hits)'; print '{}:{}'.format(sEval,eval(sEval))
        #if bDebug: sEval = 'map(Object.ToString, pts_Target_Pulled)'; print '{}:{}'.format(sEval,eval(sEval))
        
        if b_Pulled_Successes == b_Pulled_Successes_Prev:
            weightPerIter += 1.0
        else:
            weightPerIter = 1.0
        if bDebug: sEval = 'weightPerIter'; print '{}:{}'.format(sEval,eval(sEval))
        
        uvs_Target_projd_to_ns = get_uvs_of_Point3ds(ns_Morphed, pts_Target_Pulled) # Only targets that pulled to surface.
        if bDebug: sEval = 'uvs_Target_projd_to_ns[:12]'; print '{}:{}'.format(sEval,eval(sEval))
        
        dists_w = get_ws_TargetsFromPulled(ns_Morphed, pts_Target_Hits, pts_Target_Pulled, uvs_Target_projd_to_ns)
        if bDebug: sEval = 'dists_w[:12]'; print '{}:{}'.format(sEval,roundedList(eval(sEval)))
        
        nss_Morphed.append(ns_Morphed)
        deviations_PostMorph_ThisIter = [abs(w) for w in dists_w]
        deviations_PostMorph_ThisIter_Max = max(deviations_PostMorph_ThisIter)
        deviations_PostMorph_AllIters_Max.append(deviations_PostMorph_ThisIter_Max)
        transLength_ThisIter_Max = max(transLengths_ThisIter)
        transLengths_AllIters_Max.append(transLength_ThisIter_Max)
        #
        

        if bShowIterations:
            conduit.srf = ns_Morphed
            conduit.Enabled = True
            sc.doc.Views.Redraw()

            #if gSrfB:
            #    rs.DeleteObject(gSrfB)
            #gSrfB = sc.doc.Objects.AddSurface(ns_Morphed)
            #sc.doc.Views.Redraw()
        
        s  = "Iteration:{0}".format(iMorph+1)
        s += " MaxDeviation:{0:.{1}f}".format(deviations_PostMorph_ThisIter_Max, sc.doc.ModelDistanceDisplayPrecision+1)
        s += " MaxTranslation:{0:.{1}f}".format(transLength_ThisIter_Max, sc.doc.ModelDistanceDisplayPrecision+1)
        if sCmdPrompt0:
            s += "  "
        Rhino.RhinoApp.CommandPrompt = sCmdPrompt0 + s
        
        if bDebug: print s
        
        if deviations_PostMorph_ThisIter_Max <= fDistTol:
            if bDebug:
                sEvals = 'deviations_PostMorph_ThisIter_Max', 'fDistTol'; print 'Break from iterations due to {} <= {} ({} <= {})'.format(sEvals[0],sEvals[1],eval(sEvals[0]),eval(sEvals[1]))
            break
        
        if iMorph > 0:
            if deviations_PostMorph_ThisIter_Max >= min(deviations_PostMorph_AllIters_Max[:-1]):
                # Deviations are not decreasing.
                if bDebug:
                    sEvals = 'deviations_PostMorph_ThisIter_Max', 'min(deviations_PostMorph_AllIters_Max[:-1])'; print 'Break from iterations due to {} >= {} ({} >= {})'.format(sEvals[0],sEvals[1],eval(sEvals[0]),eval(sEvals[1]))
                break
            #        elif transLength_ThisIter_Max <= (0.1 * fDistTol):
            #            break
            elif (
                    (abs(
                            deviations_PostMorph_AllIters_Max[-2] -
                            deviations_PostMorph_AllIters_Max[-1]) <
                    0.1 * fDistTol)
            ):
                # Difference in maximum deviation of last 2 iterations is less than
                # 0.1*fDistTol,
                # so don't bother continuing morphing the current starting surface.
                #print "  Deviation difference: {}".format(abs(deviations_PostMorph_AllIters_Max[-2] - deviations_PostMorph_AllIters_Max[1]))
                if bDebug: print "Break from iterations due to (abs(deviations_PostMorph_AllIters_Max[-2] - deviations_PostMorph_AllIters_Max[-1]) < (0.1*fDistTol))."
                break
        
        ns_WIP = ns_Morphed
        ns_Morphed = None
        
        if deviations_PostMorph_ThisIter_Max is None:
            conduit.Enabled = False
            return
    
    conduit.Enabled = False

    if not nss_Morphed:
        return
    else:
        if bDebug:
            print '\n' + '#'*40
        
        deviation_MinOfMax = min(deviations_PostMorph_AllIters_Max)
        idx_deviation_MinOfMax = deviations_PostMorph_AllIters_Max.index(deviation_MinOfMax)
        
        transLength_MinOfMax = min(transLengths_AllIters_Max)
        idx_transLength_MinOfMax = transLengths_AllIters_Max.index(transLength_MinOfMax)
        
        if bDebug:
            s  = "After {} iteration(s),".format(iMorph+1 if ns_Morphed else iMorph)
            s += "  MaximumDeviationFromTargetPoints:{0:.{1}f}".format(
                    deviation_MinOfMax,
                    sc.doc.ModelDistanceDisplayPrecision+2)
            s += " @Iteration:{}".format(idx_deviation_MinOfMax+1)
            s += "  MaximumControlPointTranslation:{0:.{1}f}".format(
                    transLength_MinOfMax,
                    sc.doc.ModelDistanceDisplayPrecision+2)
            s += " @Iteration:{}".format(idx_transLength_MinOfMax+1)
            print s
        
        ns_MinDev = nss_Morphed[idx_deviation_MinOfMax]

        for i, ns in enumerate(nss_Morphed):
            if i == idx_deviation_MinOfMax: continue
            ns.Dispose()
        
        #if bShowIterations:
        #    if gSrfB:
        #        rs.DeleteObject(gSrfB)
        #    sc.doc.Views.Redraw()
        
        ns_MinDev.SetDomain(direction=0, domain=rgNurbsSrf.Domain(direction=0))
        ns_MinDev.SetDomain(direction=1, domain=rgNurbsSrf.Domain(direction=1))
        
        return ns_MinDev, deviation_MinOfMax