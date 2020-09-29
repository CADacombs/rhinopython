"""
190629: Created.
190825: Fixed bug in getEllipticalNurbsCurve.
190912: Added some default parameter values.
200119: Added Nurbs class.
        Modified return for some curve types and shapes.
        Modified some default parameter values.
200223: Bug fix in creating circular Arc in getArcCurve.
200619: Added a function.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

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
                print "Domains are different." # TODO: Remove this print after some time to study when this occurs.
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
