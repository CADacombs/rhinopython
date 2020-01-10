"""
160618: Created.
...
190226: Added import.
190505: Added some functions.
190513: Replaced an import with a new internal function.
190521-22: Added checks for cylinders, cones, spheres, and tori of large size.
        Merged 2 functions into one.
        Changed a default value.
190523-24: Added a function.
190617: Using NurbsSurface in place of other surfaces due to bugs in RhinoCommon stated below.
190720: A Cone-related function now check whether test point is at the Apex and returns accordingly.
190722: Split some functions.  Added a function.
190727: Corrected some output and notes.
190803: Modifications for when a value < 1e-12 is passed in as a deviation tolerance to some functions.
190810: getDescription now reports the construction curves for Rev and Sum surfaces. 
191118: Moved some functions to other modules.
191124: Import-related bug fix.
191126: Modified some printed output.

Known Issues with RhinoCommon:
    
    RevSurface.IsPlanar and
    SumSurface.IsPlanar (5 - 6.11+?):
        Some results are erroneously False.  Use NurbsSurface instead (ToNurbsSurface).
    
    RevSurface.TryGetPlane:
        Some erroneous Plane results.  Use NurbsSurface instead (ToNurbsSurface).
    
    RevSurface.TryGetSphere (5 - 6.11+?):
        Resultant sphere IsValid == False.  Use BrepFace or NurbsSurface conversion instead.
    
    RevSurface.TryGetCylinder (5 - 6.14+?):
        Cylinder not found when it should.  Convert to NurbsSurface and analyze that instead.
    
    NurbsSurface.IsSphere(tolerance) and
    NurbsSurface.TryGetSphere(tolerance) (5.? - 6.11+?):
        Sphere not found for some sizes at tolerance < 1e-9.
    
    NurbsSurface.IsTorus(tolerance) and
    NurbsSurface.TryGetTorus(tolerance) (5.? - 6.11+?):
        Torus not found for some sizes at tolerance < 1e-9.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

import System


maxExponentFor2 = 4


def getDescription(rgSrf0):
    """
    Returns: String
    """

    typeSrf = rgSrf0.GetType()
    
    tolerance = 1e-9

    # Check whether face's UnderlyingSurface is already a primitive shape.
    if typeSrf == rg.PlaneSurface:
        return "PlaneSurface"
    elif rgSrf0.IsPlanar(tolerance):
        return "Planar ({}) {}".format(tolerance, typeSrf.Name)
    #    elif (
    #            (typeSrf==rg.RevSurface or typeSrf==rg.SumSurface)
    #            and rgSrf0.IsPlanar(tolerance)
    #    ):
    #        return "Planar ({}) {}".format(tolerance, typeSrf.Name)
    
    if rgSrf0.IsSolid: s = "Solid"
    elif (rgSrf0.IsClosed(0) or rgSrf0.IsClosed(1)): s = "Closed"
    else: s = "Open"
    
    if rgSrf0.IsCylinder(tolerance):
        s += ", cylindrical ({})".format(tolerance)
    elif rgSrf0.IsCone(tolerance):
        s += ", conical ({})".format(tolerance)
    elif rgSrf0.IsSphere(tolerance):
        s += ", spherical ({})".format(tolerance)
    elif rgSrf0.IsTorus(tolerance):
        s += ", toric ({})".format(tolerance)
    
    s += " {}".format(typeSrf.Name)
    
    if typeSrf==rg.RevSurface:
        s += " with {} revolute".format(rgSrf0.Curve.GetType().Name)
    elif typeSrf==rg.SumSurface:
        s += " composed of {} and {}".format(
                rgSrf0.IsoCurve(0, rgSrf0.Domain(1).T0).GetType().Name,
                rgSrf0.IsoCurve(1, rgSrf0.Domain(0).T0).GetType().Name,
                )
    elif typeSrf == rg.NurbsSurface:
        s += "  Rational," if rgSrf0.IsRational else "  Non-rational,"
        s += "  Degrees {} x {}".format(rgSrf0.Degree(0), rgSrf0.Degree(1))
        s += "  CP counts {} x {}".format(rgSrf0.Points.CountU, rgSrf0.Points.CountV)
    
    return s


def isIsoStatusAtSeam(rgSrf, isoStatus):
    if isoStatus == Rhino.Geometry.IsoStatus.South:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Mid, rgSrf.Domain(1).Min)
    elif isoStatus == Rhino.Geometry.IsoStatus.East:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Max, rgSrf.Domain(1).Mid)
    elif isoStatus == Rhino.Geometry.IsoStatus.North:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Mid, rgSrf.Domain(1).Max)
    elif isoStatus == Rhino.Geometry.IsoStatus.West:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Min, rgSrf.Domain(1).Mid)


def isSrfSideSingularPerIsoStatus(rgSrf, isoStatus):
    if isoStatus == Rhino.Geometry.IsoStatus.South: return rgSrf.IsSingular(0)
    if isoStatus == Rhino.Geometry.IsoStatus.East: return rgSrf.IsSingular(1)
    if isoStatus == Rhino.Geometry.IsoStatus.North: return rgSrf.IsSingular(2)
    if isoStatus == Rhino.Geometry.IsoStatus.West: return rgSrf.IsSingular(3)


def pointsOfSenw(rgSrf, side):
    """
    side can be an Rhino.Geometry.IsoStatus or
    integers 0 (South), 1 (East), 2 (North), & 3 (West).
    """
    if side == 0 or side == Rhino.Geometry.IsoStatus.South:
        ptList = rgSrf.Points
        numU = ptList.CountU
        pts = [ptList.GetControlPoint(u, 0).Location for u in range(numU)]
    elif side == 1 or side == Rhino.Geometry.IsoStatus.East:
        ptList = rgSrf.Points
        numU = ptList.CountU
        numV = ptList.CountV
        pts = [ptList.GetControlPoint(numU-1, v).Location for v in range(numV)]
    elif side == 2 or side == Rhino.Geometry.IsoStatus.North:
        ptList = rgSrf.Points
        numU = ptList.CountU
        numV = ptList.CountV
        pts = [ptList.GetControlPoint(u, numV-1).Location for u in range(numU)]
    elif side == 3 or side == Rhino.Geometry.IsoStatus.West:
        ptList = rgSrf.Points
        numV = ptList.CountV
        pts = [ptList.GetControlPoint(0, v).Location for v in range(numV)]
    else: return
    
    return pts


def senwCurvePerSide(rgSrf, side):
    """
    side can be an Rhino.Geometry.IsoStatus or
    integers 0 (South), 1 (East), 2 (North), & 3 (West).
    Curve direction matches that of surface, which means that
    North and West curves are in opposition to the trim loop.
    """
    if side == 0 or side == Rhino.Geometry.IsoStatus.South:
        return rgSrf.IsoCurve(0, rgSrf.Domain(1).Min)
    elif side == 1 or side == Rhino.Geometry.IsoStatus.East:
        return rgSrf.IsoCurve(1, rgSrf.Domain(0).Max)
    elif side == 2 or side == Rhino.Geometry.IsoStatus.North:
        return rgSrf.IsoCurve(0, rgSrf.Domain(1).Max)
    elif side == 3 or side == Rhino.Geometry.IsoStatus.West:
        return rgSrf.IsoCurve(1, rgSrf.Domain(0).Min)


def shortSenws(rgSrf, fSenwLen_MinAllowed, bEcho=False):
    def maximumPointSpread(pts):
        distMax = 0.
        for i in range(len(pts)):
            ptA = pts[i]
            for j in range(i+1, len(pts)):
                ptB = pts[j]
                dist = ptA.DistanceTo(ptB)
                if dist > distMax: distMax = dist
        return distMax
    
    
    iSenwsWithinPtSpread = []
    for side in (0,1,2,3):
        if not rgSrf.IsSingular(side):
            crv = senwCurvePerSide(rgSrf, side)
            length = crv.GetLength()
            if length < fSenwLen_MinAllowed:
                iSenwsWithinPtSpread.append(side)
            else:
                pts_CP = pointsOfSenw(rgSrf, side)
                if maximumPointSpread(pts_CP) < fSenwLen_MinAllowed:
                    if bEcho: print "SENW found using maximumPointSpread instead of GetLength()."
                    iSenwsWithinPtSpread.append(side)
    return iSenwsWithinPtSpread


