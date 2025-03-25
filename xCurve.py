"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190629: Created.
190825: Fixed bug in getEllipticalNurbsCurve.
190912: Added some default parameter values.
200119: Added Nurbs class.
        Modified return for some curve types and shapes.
        Modified some default parameter values.
200223: Bug fix in creating circular Arc in getArcCurve.
200619, 220317: Added functions.
230108: Bug fix in filterCurvesOnSurface for ignoring curves with only endpoints on face.
230814: Replaced a couple references to ModelAbsoluteTolerance with tolerance parameter.
230901: Added code for debugging.
250324: Added code for Rhino V5.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

from clr import StrongBox
from System import Array

import xBrepFace

MAX_ACCEPTABLE_SIZE_MM = 1E4


class Nurbs():

    @staticmethod
    def hasInternalPolyknots(nc0):
        if nc0.SpanCount == 1: return False
        if nc0.Degree == 1: return False
        knots = nc0.Knots
        if nc0.IsPeriodic:
            iK = 0
            nK = knots.Count - 2 # 1 is unnecessary.
        else:
            iK = nc0.Degree
            nK = -nc0.Degree + knots.Count - 2 # 1 is unnecessary.
        while True:
            m = knots.KnotMultiplicity(iK)
            if m > 1:
                return True
            iK += m
            if iK > nK:
                return False


    @staticmethod
    def hasSomeFullyMultiplePolyknots(nc0):
        if nc0.SpanCount == 1: return False
        if nc0.Degree == 1: return False
        knots = nc0.Knots
        if nc0.IsPeriodic:
            iK = 0
            nK = knots.Count - 2 # 1 is unnecessary.
        else:
            iK = nc0.Degree
            nK = -nc0.Degree + knots.Count - 2 # 1 is unnecessary.
        while True:
            m = knots.KnotMultiplicity(iK)
            if m == nc0.Degree:
                return True
            iK += m
            if iK > nK:
                return False


    @staticmethod
    def isUniform(nc0):
        if nc0.SpanCount == 1: return True
        if nc0.Degree == 1: return False
        knots = nc0.Knots
        if nc0.IsPeriodic:
            iK = 0
            nK = knots.Count - 1
            tDelta0 = knots[1] - knots[0]
        else:
            iK = nc0.Degree
            nK = knots.Count - nc0.Degree
            tDelta0 = knots[iK] - knots[iK-1]
        while True:
            m = knots.KnotMultiplicity(iK)
            if m > 1:
                return False
            iK += m
            tDelta = knots[iK] - knots[iK-1]
            if abs(tDelta - tDelta0) > 1e-9:
                return False
            if iK >= nK:
                return True


def duplicateSegments(rgCrvs, bExplodePolyCrvs=True):
    """
    Uses Curve.DuplicateSegments, but also:
        Returns single segments when only 1 segment exists in input.
        Includes option to explode PolyCurves.
    """


    def dupSegs_1Crv(crv_In):
        cs_Out = []
        
        if isinstance(crv_In, rg.BrepEdge):
            edge = crv_In
            if crv_In.Domain.CompareTo(edge.EdgeCurve.Domain) == 0:
                # Domains are the same.
                crv_In = edge.EdgeCurve
            else:
                print("Domains are different.",
                    "TODO: Remove this print after some time to study when this occurs.")
                crv_In = edge.EdgeCurve.Trim(domain=crv_In.Domain)

        if isinstance(crv_In, rg.PolyCurve):
            if bExplodePolyCrvs:
                pc_Dup = crv_In.DuplicatePolyCurve()
                pc_Dup.RemoveNesting() # Returns boolean.
                cs_Out = pc_Dup.Explode()
                pc_Dup.Dispose()
            else:
                cs_Out = crv_In.DuplicateSegments()
        elif isinstance(crv_In, rg.PolylineCurve):
            if crv_In.PointCount > 2:
                cs_Out = crv_In.DuplicateSegments()
            else:
                # Because DuplicateSegments produces no output for single span PolylineCurves.
                cs_Out = [rg.LineCurve(crv_In.PointAtStart, to=crv_In.PointAtEnd)]
        else:
            cs_Out = [crv_In.DuplicateCurve()]

        return cs_Out

    cs_In = rgCrvs

    try:
        len(cs_In)
        cs_Out = []
        for crv_In in cs_In:
            cs_Ret = dupSegs_1Crv(crv_In)
            cs_Out.extend(cs_Ret)
        return cs_Out
    except:
        # Possibly 1 curve.
        return dupSegs_1Crv(cs_In)


def getArcCurve(rgCrv0, bTolByRatio=False, fTolRatio=1000.0, fDevTol=1e-9, fMinNewCrvLen=None, fMaxRadius=None):
    """
    Uses a different approach of obtaining the Arc or Circle than Curve.TryGetArc.
    """
    
    sType = rgCrv0.GetType().Name

    if rgCrv0 is None:
        return None, "Geometry not found!  It will be skipped."
    if sType == 'ArcCurve':
        return rgCrv0.DuplicateCurve(), 0.0
    if sType == 'LineCurve':
        return None, "Skipped LineCurve."
    if rgCrv0.IsLinear(1e-9):
        return None, "Skipped linear {}.".format(sType)

    if fMinNewCrvLen is None: fMinNewCrvLen = 10.0*sc.doc.ModelAbsoluteTolerance

    if rgCrv0.IsClosed:
        t0, t1 = rgCrv0.DivideByCount(segmentCount=3, includeEnds=False)[:2]
        rgArcCrv1 = rg.ArcCurve(
            rg.Circle(
                rgCrv0.PointAtStart,
                rgCrv0.PointAt(t0),
                rgCrv0.PointAt(t1),
            )
        )
    else:
        rgArcCrv1 = rg.ArcCurve(
            rg.Arc(
                rgCrv0.PointAtStart,
                rgCrv0.PointAt(
                        rgCrv0.DivideByCount(segmentCount=2, includeEnds=False)[0]),
                rgCrv0.PointAtEnd
            )
        )

    s = ""
    length_crv1 = rgArcCrv1.GetLength()
    if length_crv1 < fMinNewCrvLen:
        s += "Curve is too short"
    
    if fMaxRadius is None:
        modelUnitPerMm = Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters,
            sc.doc.ModelUnitSystem)
        fMaxRadius = MAX_ACCEPTABLE_SIZE_MM * modelUnitPerMm
    
    if rgArcCrv1.Radius > fMaxRadius:
        if s: s += " and "
        s += "Radius is too large."
    else:
        if s: s += "."
    if s:
        return None, s
    
    rc = rg.Curve.GetDistancesBetweenCurves(
        rgCrv0,
        rgArcCrv1,
        sc.doc.ModelAbsoluteTolerance)
    if not rc[0]:
        return None, "Deviation was not returned from GetDistancesBetweenCurves."
    
    fDev = rc[1]
    
    # Check deviation to tolerance.
    if bTolByRatio:
        fRatio = length_crv1 / fDev
        if fRatio < fTolRatio:
            s  = "Curve was not converted because its (length / deviation) ratio is too small."
            return None, s
    else:
        # Check absolute deviation.
        if fDev > fDevTol:
            return None, fDev
    
    return rgArcCrv1, fDev


def getCurvesToSplitSurface(rgCrvs_tryForSplitters, rgSrf, fTolerance, bDebug=False):
    """
    Duplicates input curves that are on surface and splits to edges.
    """
    rc = filterCurvesOnSurface(
        rgCrvs_tryForSplitters,
        rgSrf,
        fTolerance=fTolerance,
        bDebug=bDebug)
    if not rc: return
    (
        crvs_CompletelyOnSrf,
        crvs_PartiallyOnSrf,
        ) = rc

    if not crvs_PartiallyOnSrf:
        return crvs_CompletelyOnSrf

    # Trim the curves that lie partially on the surface mostly to get rid of the overlapping (Boundary) portions.
    rc = splitCurvesWithSurfaceEdges(
        rgCrvs_In=crvs_PartiallyOnSrf,
        rgSrf_In=rgSrf,
        fTolerance=fTolerance)
    (
        crvs_Trimmed_Interior,
        crvs_Trimmed_Boundary,
        crvs_Trimmed_Exterior
        ) = rc

    for c in crvs_Trimmed_Boundary + crvs_Trimmed_Exterior: c.Dispose()
    for c in crvs_PartiallyOnSrf: c.Dispose()

    return crvs_CompletelyOnSrf + crvs_Trimmed_Interior


def getEllipticalNurbsCurve(rgCrv0, bTolByRatio=False, fTolRatio=1000.0, fDevTol=1e-9, fMinNewCrvLen=None, fMaxRadius=None):
    """
    Uses a different approach of obtaining the Arc or Circle than Curve.TryGetArc.
    """
    
    if rgCrv0 is None:
        return None, "Geometry not found!  It will be skipped."
    if isinstance(rgCrv0, rg.ArcCurve):
        return rgCrv0.ToNurbsCurve(), 0.0
    if isinstance(rgCrv0, rg.LineCurve):
        return None, "Skipped LineCurve."
    if rgCrv0.IsLinear(1e-9):
        return None, "Skipped linear {}.".format(sType)

    if fMinNewCrvLen is None: fMinNewCrvLen = 10.0*sc.doc.ModelAbsoluteTolerance

    rc = rgCrv0.TryGetEllipse(fDevTol)
    
    if not rc[0]:
        return None, "Ellipse within tolerance not returned from TryGetEllipse."
     
    rgEllipse = rc[1]
    
    rgNurbsCurve1 = rgEllipse.ToNurbsCurve()
    
    if fMaxRadius is None:
        modelUnitPerMm = Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters,
            sc.doc.ModelUnitSystem)
        fMaxRadius = MAX_ACCEPTABLE_SIZE_MM * modelUnitPerMm
    
    s = ""
    length_crv1 = rgNurbsCurve1.GetLength()
    if length_crv1 < fMinNewCrvLen:
        s += "Curve is too short"
    if rgEllipse.Radius1 > fMaxRadius:
        if s: s += " and "
        s += "Radius is too large."
    if rgEllipse.Radius2 > fMaxRadius:
        if s: s += " and "
        s += "Radius is too large."
    else:
        if s: s += "."
    if s:
        return None, s
    
    rc = rg.Curve.GetDistancesBetweenCurves(
        rgCrv0,
        rgNurbsCurve1,
        sc.doc.ModelAbsoluteTolerance)
    if not rc[0]:
        return None, "Deviation was not returned from GetDistancesBetweenCurves."
    
    fDev = rc[1]
    
    # Check deviation to tolerance.
    if bTolByRatio:
        fRatio = length_crv1 / fDev
        if fRatio < fTolRatio:
            s  = "Curve was not converted because its (length / deviation) ratio is too small."
            return None, s
    else:
        # Check absolute deviation.
        if fDev > fDevTol:
            return None, fDev
    
    return rgNurbsCurve1, fDev


def filterCurvesOnSurface(rgCrvs_In, rgSrf, fSamplingResolution=None, fTolerance=None, bDebug=False):
    """
    Parameteters:
        curves
        rgSrf: If a BrepFace, its BrepEdges are used instead of the underlying surface natural border.
        fSamplingResolution: float of curve division length.
        fTolerance: float of distance from curve to surface.
        bDebug
    Returns:
        list(Curves of rgCrvs_In completely on rgSrf),
        list(Curves of rgCrvs_In only partially on rgSrf and intersects at least one (sur)face border)
    """


    try: rgCrvs_In = list(rgCrvs_In)
    except: rgCrvs_In = [rgCrvs_In]

    if fSamplingResolution is None:
        fSamplingResolution = 100.0*sc.doc.ModelAbsoluteTolerance
    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance

    if isinstance(rgSrf, rg.BrepFace):
        rgFace = rgSrf
        rgB_Temp = None
    else:
        rgB_Temp = rgSrf.ToBrep()
        rgFace = rgB_Temp.Faces[0]


    def overlapCurvesOnFace(curve, face, tolerance):
        """
        Unfortunately, at least up through 7.25, Intersection.CurveBrepFace
        sometimes doesn't produce overlapCurves as expected, resulting in false negatives.
        """
        rc = rg.Intersect.Intersection.CurveBrepFace(
            curve, face, tolerance)

        bSuccess, overlapCurves, intersectionPoints = rc
        if not bSuccess: return []

        return overlapCurves


    def getCrvDivisionPts(crv):

        strongBox_points = StrongBox[Array[rg.Point3d]]()

        if Rhino.RhinoApp.ExeVersion == 5:
            rc = crv.DivideByLength(
                fSamplingResolution,
                includeStart=True,
                points=strongBox_points)
        else:
            rc = crv.DivideByLength(
                fSamplingResolution,
                includeEnds=True,
                points=strongBox_points)

        if rc:
            return list(strongBox_points.Value)

        crv_GetLength = crv.GetLength()

        if crv_GetLength <= fTolerance:
            if bDebug:
                print("No points for curve that is {} long.".format(
                    crv_GetLength))
                #sc.doc.Objects.AddCurve(crv)
            return []
        else:
            rc = crv.DivideByCount(
                segmentCount=2,
                includeEnds=True,
                points=strongBox_points)
            pts = list(strongBox_points.Value)
            #for pt in pts: sc.doc.Objects.AddPoint(pt)
            #sc.doc.Views.Redraw(); 1/0

        return pts


    def count_crvPtsOnFace(pts_Crv, face, fTolerance):

        iCt_CrvPts_On_Face = 0

        for pt in pts_Crv:
            rc = xBrepFace.is3dPointOnFace(
                face, pt, fTolerance)
            
            if rc:
                iCt_CrvPts_On_Face += 1
            else:
                pass
                #sc.doc.Objects.AddPoint(pt); 1/0

        return iCt_CrvPts_On_Face


    def isCrvCompletelyOnFace(crv, pts_Crv):

        # Quick check for False.
        rc = xBrepFace.is3dPointOnFace(
            rgFace, crv.PointAtStart, fTolerance)

        if not rc:
            # This includes None, False, and rg.PointFaceRelation.Exterior.
            return False


        pts = getCrvDivisionPts(crv)
        if len(pts) == 0: return False

        for pt in pts_Crv:
            rc = xBrepFace.is3dPointOnFace(
                rgFace, pt, fTolerance)

            if not rc:
                # This includes None, False, and rg.PointFaceRelation.Exterior.
                return False

        return True


    def isCurveCompletelyOnFaceBorder(crv):

        strongBox_points = StrongBox[Array[rg.Point3d]]()

        rc = crv.DivideByLength(
            fSamplingResolution,
            includeEnds=True,
            points=strongBox_points)

        if rc:
            pts = list(strongBox_points.Value)
            #for pt in pts: sc.doc.Objects.AddPoint(pt)
            #sc.doc.Views.Redraw(); 1/0
            if len(pts) == 2:
                rc = crv.DivideByCount(
                    segmentCount=2,
                    includeEnds=True,
                    points=strongBox_points)
                pts = list(strongBox_points.Value)
        else:
            crv_GetLength = crv.GetLength()
            
            if crv_GetLength <= fTolerance:
                if bDebug:
                    print("No points for curve that is {} long.".format(
                        crv_GetLength))
                    #sc.doc.Objects.AddCurve(crv)
                return
            else:
                rc = crv.DivideByCount(
                    segmentCount=2,
                    includeEnds=True,
                    points=strongBox_points)
                pts = list(strongBox_points.Value)
                #for pt in pts: sc.doc.Objects.AddPoint(pt)
                #sc.doc.Views.Redraw(); 1/0

        for pt in pts:
            # TODO: Review the results of the following using more tolerance for natural surfaces.
            rc = xBrepFace.is3dPointOnFace(
                rgFace, pt, 1.5*fTolerance if rgB_Temp else fTolerance)

            if rc != rg.PointFaceRelation.Boundary:
                return False

        return True


    def doesCrvIntersectFaceBorder(crv, rgFace):
        for edge in [rgFace.Brep.Edges[iE] for iE in rgFace.AdjacentEdges()]:
            intrscts = rg.Intersect.Intersection.CurveCurve(
                curveA=crv,
                curveB=edge,
                tolerance=fTolerance,
                overlapTolerance=0.1*fTolerance)

            if bDebug: print("Intersection count: {}".format(intrscts.Count))

            if intrscts.Count > 0:
                return True

        return False


    #sc.doc.Objects.AddBrep(rgFace.Brep); sc.doc.Views.Redraw(); 1/0

    cs_CompletelyOn = []
    cs_PartiallyOn_BorderIntersect = []

    for c in rgCrvs_In:
        pts_Crv = getCrvDivisionPts(c)
        if len(pts_Crv) == 0:
            continue
        ct_ptsOnFace = count_crvPtsOnFace(pts_Crv, rgFace, fTolerance)
        #if len(overlapCurvesOnFace(c, rgFace, fTolerance)) == 0:
        #    continue
        if ct_ptsOnFace == 0:
            continue
        if ct_ptsOnFace == len(pts_Crv):
            if isCurveCompletelyOnFaceBorder(c):
                continue
            cs_CompletelyOn.append(c.DuplicateCurve())
            continue

        # Partially on face
        if doesCrvIntersectFaceBorder(c, rgFace):
            cs_PartiallyOn_BorderIntersect.append(c.DuplicateCurve())
        else:
            print("Curve is partially on face but does not intersect any border.")
    #map(sc.doc.Objects.AddCurve, cs_CompletelyOn); 1/0; return


    if rgB_Temp: rgB_Temp.Dispose()

    return cs_CompletelyOn, cs_PartiallyOn_BorderIntersect


def splitCurvesWithSurfaceEdges(rgCrvs_In, rgSrf_In, fTolerance=None, bDebug=False):
    """
    Parameters:
        rgCrvs_In
        rgSrf_In: If rg.BrepFace, then adjacent edges are used, otherwise surface domain isocurves are used.
        fTolerance
        bDebug
    Returns:
        rgCrvs_Out_SrfInterior
        rgCrvs_Out_SrfBoundary
        rgCrvs_Out_SrfExterior
        
    Intersect.Intersection.CurveCurve sometimes returns intersections at the exact ends of curves.
    In this case, the entire curve will be the result of Curve.Split.
    """
    
    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance

    rgEdges = [] # Face edges or Surface domain extent curves.

    if isinstance(rgSrf_In, rg.BrepFace):
        idxs_Edges = rgSrf_In.AdjacentEdges()
        for idxE in idxs_Edges:
            rgE = rgSrf_In.Brep.Edges[idxE]
            rgEdges.append(rgE)
    elif isinstance(rgSrf_In, rg.Surface):
        # SENW == 0,1,2,3 for IsSingular.
        
        # South.
        if not rgSrf_In.IsSingular(0):
            rgC = rgSrf_In.IsoCurve(0, rgSrf_In.Domain(1).T0)
            rgEdges.append(rgC)
        
        # East.
        if not rgSrf_In.IsSingular(1):
            rgC = rgSrf_In.IsoCurve(1, rgSrf_In.Domain(0).T1)
            rgEdges.append(rgC)
        
        # North.
        if not rgSrf_In.IsSingular(2):
            rgC = rgSrf_In.IsoCurve(0, rgSrf_In.Domain(1).T1)
            rgEdges.append(rgC)
        
        # West.
        if not rgSrf_In.IsSingular(3):
            rgC = rgSrf_In.IsoCurve(1, rgSrf_In.Domain(1).T0)
            rgEdges.append(rgC)
    else:
        raise ValueError("{} passed to splitCurvesWithSurfaceEdges".format(
            rgSrf_In.GetType().Name))


    rgCrvs_Out_SrfInterior = []
    rgCrvs_Out_SrfBoundary = []
    rgCrvs_Out_SrfExterior = []

    for rgC_toSplit in rgCrvs_In:
        t_Pts = []
        doms_Overlaps = []
        
        for rgEdge in rgEdges:
            intrscts = rg.Intersect.Intersection.CurveCurve(
                curveA=rgC_toSplit,
                curveB=rgEdge,
                tolerance=fTolerance,
                overlapTolerance=0.1*fTolerance)

            if bDebug: print("Intersection count: {}".format(intrscts.Count))

            for intsct in intrscts:
                if intsct.IsPoint:
                    #print(intsct.ParameterA)
                    t_Pts.append(intsct.ParameterA)
                if intsct.IsOverlap:
                    doms_Overlaps.append(intsct.OverlapA)
                    #print(intsct.OverlapA.T0, intsct.OverlapA.T1)
                    t_Pts.append(intsct.OverlapA.T0)
                    t_Pts.append(intsct.OverlapA.T1)

        t_Pts = sorted(set(t_Pts))
        
        # TODO: Remove near-duplicates.
        
        rgCs_fromSplit = rgC_toSplit.Split(t=t_Pts)


        # For debugging.
        if len(rgCs_fromSplit) == 1:
            pass
            #print("1 curve resulted from split.")
            #c = rgCs_fromSplit[0]
            #sc.doc.Objects.AddCurve(c)
            #continue


        for c in rgCs_fromSplit:

            fLength = c.GetLength()

            if fLength == 0.0:
                if bDebug: print("Zero length curve ignored.")
                c.Dispose()
                continue
            elif fLength < fTolerance:
                if bDebug:
                    print("Short ({}) curve ignored.  Domain: [{},{}]".format(
                        formatDistance(fLength),
                        c.Domain.T0,
                        c.Domain.T1))
                c.Dispose()
                continue

            for dom in doms_Overlaps:
                if c.Domain.IncludesParameter(dom.Mid):
                    rgCrvs_Out_SrfBoundary.append(c)
                    break # To next curve from split.
            else:
                # Curve is not part of Intersection overlaps but may still be along boundary.

                pt_MidCrv = c.PointAt(c.Domain.Mid)
                #sc.doc.Objects.AddPoint(pt_MidCrv)

                b, u, v = rgSrf_In.ClosestPoint(pt_MidCrv) # if rgSrf_In is BrepFace, UnderlyingSurface is still used.
                if not b:
                    raise ValueError("ClosestPoint failed.")

                ptAtUV = rgSrf_In.PointAt(u,v)
                #sc.doc.Objects.AddPoint(ptAtUV)


                # Check whether Closest point on surfaces is beyond tolerance from the test point.
                dist = ptAtUV.DistanceTo(pt_MidCrv)
                if dist > fTolerance:
                    rgCrvs_Out_SrfExterior.append(c)
                    continue


                # Point is either interior or along an edge.

                if isinstance(rgSrf_In, rg.BrepFace):

                    # TODO: Check distance from edge?

                    ptFaceRel = rgSrf_In.IsPointOnFace(u,v)
                    if ptFaceRel == rg.PointFaceRelation.Interior:
                        rgCrvs_Out_SrfInterior.append(c)
                    else:
                        rgCrvs_Out_SrfBoundary.append(c)
                else:
                    # rgSrf_In is a Surface but not a BrepFace.

                    # TODO: Check distance from surface side?

                    if (
                        rgSrf_In.Domain(0).IncludesParameter(u) and
                        rgSrf_In.Domain(1).IncludesParameter(v)
                    ):
                        rgCrvs_Out_SrfInterior.append(c)
                    else:
                        rgCrvs_Out_SrfBoundary.append(c)


    return rgCrvs_Out_SrfInterior, rgCrvs_Out_SrfBoundary, rgCrvs_Out_SrfExterior


