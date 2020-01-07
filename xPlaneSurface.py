"""
160702: Added planeSurfaceFromPlaneAndFace() and planeSurfaceOfPlaneSurfaceFace().
160711: Redesigned createFromPlaneSurfaceFace().
160723: createFromPlaneSurfaceFace(): Surface is now oversized.
190128: Started this script as a branch from another.
190202: Renamed script.
190204: Updated bounding box in createFromPlaneAndObjectSize.
190207: Removed a couple import references.
190219: Modified a variable name and some comments.
190318: These notes were moved to their own file.
190319: These notes were returned to this file.
190407: Added a function.  Removed the function.
190507: Modifed a function.  Commented out another function.
190519: Fixed typo in an argument name.
190619: BoundingBox in createFromPlaneAndObjectSize is now inflated 10.0*sc.doc.ModelAbsoluteTolerance.
        Added bDebug to a function.
191126: Bug fixes in createFromFace.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import xPrimitiveShape


def createFromFace(rgFace0, bMatchToShrunkFace=True, fTolerance=0.1*sc.doc.ModelAbsoluteTolerance, bDebug=False):
    """
    Return on success: tuple(
        tuple(
            PlaneSurface(That at least covers the entire BrepFace)
            float(Tolerance used)
            bool(Whether shrunk surface was used),
        None)
    Return on fail: tuple(
        None
        str(LogStatement))
    """
    
    rgSrf_Face0 = rgFace0.UnderlyingSurface()
    
    if isinstance(rgSrf_Face0, rg.PlaneSurface):
        return None, "Surface is already a PlaneSurface."

    rc = xPrimitiveShape.BrepFace.tryGetPlane(
            rgFace0,
            bMatchToShrunkFace=bMatchToShrunkFace,
            fTolerance=fTolerance,
            bDebug=bDebug)
    if rc[0] is None:
        return rc
    else:
        rgPlane1, fTol_Used, bShrunkUsed = rc[0]
    
    rgSrf1 = createFromPlaneAndObjectSize(
        rgPlane=rgPlane1,
        obj_ForSize=rgFace0,
        bDebug=bDebug)
    if rgSrf1 is None:
        return None, "createFromPlaneAndObjectSize fail."
    
    return (rgSrf1, fTol_Used, bShrunkUsed), None


def createFromPlaneAndObjectSize(rgPlane, obj_ForSize, bDebug=False):
    if not isinstance(obj_ForSize, (list, tuple, set)):
        obj_ForSize = [obj_ForSize]
    
    bbox_Use = rg.BoundingBox.Empty
    
    plane_to_world = rg.Transform.ChangeBasis(
            plane0=rgPlane, plane1=rg.Plane.WorldXY)
    
    for obj in obj_ForSize:
        geom = rs.coercegeometry(obj)
        if isinstance(geom, rg.BrepFace):
            rgBrep_1Face = geom.DuplicateFace(duplicateMeshes=False)
            bBoxX = rgBrep_1Face.GetBoundingBox(plane=rgPlane)
            rgBrep_1Face.Dispose()
        else:
            bBoxX = geom.GetBoundingBox(plane=rgPlane)
        
        if not bBoxX.IsValid: return
        
        bbox_Use.Union(other=bBoxX)
    
    bbox_Use.Inflate(10.0*sc.doc.ModelAbsoluteTolerance)
    
    bbox_Use.Transform(xform=plane_to_world)
    
    if bDebug:
        box_Use = rg.Box(bbox_Use)
        sc.doc.Objects.AddBox(box_Use)
        sc.doc.Views.Redraw()
    
    rgSrf1 = rg.PlaneSurface.CreateThroughBox(plane=rgPlane, box=bbox_Use)
    return rgSrf1


# 190507: Retire attempt.
#def createFromPlaneAndSurface(rgPlane, rgSrf0):
#    bBox_rgSrf0 = rgSrf0.GetBoundingBox(True)
#    if not bBox_rgSrf0.IsValid: return
#    diagLen = bBox_rgSrf0.Min.DistanceTo(bBox_rgSrf0.Max)
#    bBox_rgSrf0.Inflate(0.05*diagLen)
#    rgSrf1 = rg.PlaneSurface.CreateThroughBox(rgPlane, bBox_rgSrf0)
#    return rgSrf1


def createFromPlaneSurfaceFace(rgFace):
    domainU = rgFace.Domain(0)
    domainV = rgFace.Domain(1)
    
    rgSrf = rgFace.UnderlyingSurface()
    if type(rgSrf) == Rhino.Geometry.PlaneSurface:
        b, rgPlane = rgFace.TryGetPlane()
        if b:
            return rg.PlaneSurface(rgPlane, domainU, domainV)
