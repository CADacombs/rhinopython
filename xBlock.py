"""
160618: Created by importing functions from SrfRadius.py.
160625: Fixed bugs in brepOrExtrusionOfBlock_WireframePick() and tryPickedFaceOfBlock().
160702: Changed EpsilonEquals tolerance in tryPickedFaceOfBlock from Rhino.RhinoMath.ZeroTolerance to trying from 1.e-11 to sc.doc.ModelAbsoluteTolerance.
160713: Fixed bug in tryPickedFaceOfBlock() and made more efficient by looping through the faces of the brep only once.
        Modified tolerance value in brepOrExtrusionOfBlock_WireframePick() from ...ZeroTolerance to (.001 * sc.doc.ModelAbsoluteTolerance).
160719: In tryPickedFaceOfBlock(), changed "return" to "return None, None".
190314: Updated an import name.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc

import xLayer


def brepOrExtrusionOfBlock_ShadedPick(rdInstRef, ptPicked, planeCam):
    """Create intersecting line that passes through the actual picked point.
    Its vector is that of the camera view.
    The line extends to the bounding box of the block instance."""
    
    vector = rg.Vector3d.Negate(planeCam.ZAxis)
    ptOnCameraPlane = planeCam.ClosestPoint(ptPicked)
    pts_bbox = rs.BoundingBox(rdInstRef.Id)
    bbox = rg.BoundingBox(pts_bbox[0], pts_bbox[6])
    
    rgLine = rg.Line(ptOnCameraPlane, vector, 1.)
    rgLine.ExtendThroughBox(bbox)
    rgLineCrv = rg.LineCurve(rgLine)
    
    rgBreps_Intrsct, pts_Intrsct_All = (
            brepsOrExtrusionsAndPtsInBlockInstThatIntrsctLine(rdInstRef, rgLineCrv))
    
    rgObj_Nearest = pt_Nearest = None
    
    for i, pts_Intrsct in enumerate(pts_Intrsct_All):
        for pt_Intrsct in pts_Intrsct:
            fDistPtToPlane = ptOnCameraPlane.DistanceTo(pt_Intrsct)
            if pt_Nearest is None or fDistPtToPlane < fDistPtToPlane_Min:
                rgObj_Nearest = rgBreps_Intrsct[i]
                pt_Nearest = pt_Intrsct
                fDistPtToPlane_Min = fDistPtToPlane
    
    return rgObj_Nearest, pt_Nearest


def brepOrExtrusionOfBlock_WireframePick(rdInstRef, pt, xform=rg.Transform.Identity):
    
    xform *= rdInstRef.InstanceXform
    rdInstDef = rdInstRef.InstanceDefinition
    rdObjs = rdInstDef.GetObjects()
    
    for rdObj in rdObjs:
        if not xLayer.areLayerAndAllAncestorsVisible(rdObj.Attributes.LayerIndex):
            continue
        elif rdObj.ObjectType == rd.ObjectType.InstanceReference:
            rgObj = brepOrExtrusionOfBlock_WireframePick(rdObj, pt, xform)
            if rgObj: return rgObj # Otherwise, continue with loop.
        else: # Object is not a block instance
            rgObj = rdObj.Geometry
            rgObj.Transform(xform)
            if rgObj.ObjectType == rd.ObjectType.Brep:
                ptOn = rgObj.ClosestPoint(pt)
            elif rgObj.ObjectType == rd.ObjectType.Extrusion:
                bPt, u, v = rgObj.ClosestPoint(pt)
                ptOn = rgObj.PointAt(u, v)
            else: continue # Not brep or extrusion.
            
            if pt.EpsilonEquals(ptOn, .001 * sc.doc.ModelAbsoluteTolerance):
                return rgObj


def brepsOrExtrusionsAndPtsInBlockInstThatIntrsctLine(rdInstRef, rgLineCrv, xform=rg.Transform.Identity, rgBreps_Intrsct=None, pts_Intrsct_All=None):
    
    if rgBreps_Intrsct is None: rgBreps_Intrsct = []
    if pts_Intrsct_All is None: pts_Intrsct_All = []
    
    fTol = .1 * sc.doc.ModelAbsoluteTolerance # Tolerances < 0  in Intersection.CurveBrep() seem to produce results within 1e-14 of each other.
    xform *= rdInstRef.InstanceXform
    rdInstDef = rdInstRef.InstanceDefinition
    rdObjs = rdInstDef.GetObjects()
    
    for rdObj in rdObjs:
        if not xLayer.areLayerAndAllAncestorsVisible(rdObj.Attributes.LayerIndex):
            continue
        elif rdObj.ObjectType == (
                rd.ObjectType.InstanceReference):
            rgBreps_Intrsct, pts_Intrsct_All = (
                    brepsOrExtrusionsAndPtsInBlockInstThatIntrsctLine(
                    rdObj, rgLineCrv, xform, rgBreps_Intrsct, pts_Intrsct_All))
        else: # Object is not a block instance
            rgObj = rdObj.Geometry
            if rgObj.ObjectType == rd.ObjectType.Brep:
                pass
            elif rgObj.ObjectType == rd.ObjectType.Extrusion:
                rgObj = rdObj.Geometry.ToBrep()
            else: continue
            
            rgObj.Transform(xform)
            bIntrsct, rgCrvs_Intrsct, pts_Intrsct = (
                    rg.Intersect.Intersection.CurveBrep(
                            rgLineCrv, rgObj, fTol))
            if len(pts_Intrsct) > 0: # bIntrsct = True even when no intersection
                rgBreps_Intrsct.append(rgObj)
                pts_Intrsct_All.append(tuple(pts_Intrsct))
            
            #return rgBreps_Intrsct, pts_Intrsct_All
    
    return rgBreps_Intrsct, pts_Intrsct_All


def tryPickedFaceOfBlock(rdObj, ptPicked):
    """
    Parameters:
        rdObj
        ptPicked
    Returns:
        rgFace
        ptPicked
    
    Although it would be better to capture the view properties
    immediately after the point is picked, the code is located here so
    that planeCam doesn't need to be passed to the current function.
    planeCam will be used when a shaded brep face is picked.
    """
    
    viewport = sc.doc.Views.ActiveView.MainViewport
    bCamPlane, planeCam = viewport.GetCameraFrame()
    if bCamPlane is None: return None, None
    
    bWireFrameView = not viewport.DisplayMode.SupportsShading
    if bWireFrameView:
        rgObj = brepOrExtrusionOfBlock_WireframePick(
                rdObj, ptPicked) # Get brep
    else: # Shaded brep face was picked.
        rgObj, ptPicked = brepOrExtrusionOfBlock_ShadedPick(
                rdObj, ptPicked, planeCam) # Note that ptPicked is reassigned.
    
    if rgObj is not None and (rgObj.ObjectType ==
            rd.ObjectType.Brep or rd.ObjectType.Extrusion): # Get face
        if rgObj.ObjectType == rd.ObjectType.Extrusion:
            rgObj = rgObj.ToBrep()
        fDist_Min = 1./Rhino.RhinoMath.ZeroTolerance
        for f, rgFace in enumerate(rgObj.Faces):
            bPt, u, v = rgFace.ClosestPoint(ptPicked)
            if rgFace.IsPointOnFace(u, v): # Works because values are 0 for Exterior, 1 for Interior, 2 for Boundary
                ptOn = rgFace.PointAt(u, v)
                fDist = ptPicked.DistanceTo(ptOn)
                if fDist < fDist_Min:
                    idxFace_Closest = f
                    fDist_Min = fDist
        #sPrint = 'fDist_Min'; print sPrint + ':', eval(sPrint)
        return rgObj.Faces[idxFace_Closest], ptPicked
    else: rgFace = None
    
    return None, None