"""
160531-0602: Created by importing functions from other scripts.
160607: Addition.
160609: Mod.
160611: Major additions for removeTrim.
160612: Bug fix in addOuterLoop.
160623: Added retrimOneFaceBrep() and trimOneFaceBrepWithCurves().
160630: Added some functions.
160706: Moved some functions to another script.
160707: Some functions now try to matching bounding boxes at ModelAbsoluteTolerance then, if necessary, 10*ModelAbsoluteTolerance.
        Added a function.
160708: Added a function.
160725: Added a function.
160819: Added indicesOfInvalidBrepFaces().
160824: Fixed bug in orderedLists_AdjFLT_perSelFBorders and renamed to orderedLists_AdjacentFLT_perSelFBorders.
160826: trimOneFaceBrepWithCurves(): Default value of bShrinkFaces changed from True to False.
160829: removeTrim: Fixed initializing list bug.
        Renamed ptOnFaceMinDistFromBorder() to ptOnFaceNotOnBorder and added default value for fDistMin.
160831: Moved addOuterLoop, indicesOfTrimsNotToRemove, removeTrim to untrimLocal.py
160901: Updated debugging.
160902: tryExistingEdgeMatchingCurvx: Replaced length comparison with curve deviation check.
160911: Moved orderedLists_AdjacentFLT_perSelFBorders to another script.
170630: ptOnFaceNotOnBorder: Now skips faces whose AreaMassProperties.Compute's output is None.
181020: Removed indicesOfInvalidBrepFaces.
181029: addSenwEdgeTrimToLooX: Added code with intrvl.Swap() to avoid
190128: Moved indicesOfContiguousFilletFaces to spbBrep.py.
190129: Updated a comment.
190201: Copied some code from spb_rc5\brep.py.
190505, 200107, 200401: Modified this history log.
200505: Simplified a function.
200522: Bug fix.
"""

import Rhino
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc


def epsilonEqualsBoundingBox(bBoxA, bBoxB, fTol=sc.doc.ModelAbsoluteTolerance):
    """
    Created because BoundingBox structure does not contain an EpsilonEquals method.
    """
    boxA = Rhino.Geometry.Box(bBoxA)
    if not boxA: return
    boxB = Rhino.Geometry.Box(bBoxB)
    if not boxB: return
    return boxA.EpsilonEquals(boxB, fTol)


def usingBoundingBoxOfBrep(rgBrep_toSearch, rgBrep_ForBBox, bDebug=False):
    """
    rgBrep_ForBBox will most likely be a brep or a face.
    """
    if bDebug: print 'xBrep_findMatchingFace.usingBoundingBoxOfBrep()...'
    
    rgBbox_B0 = rgBrep_ForBBox.GetBoundingBox(True)
    
    for e in 0,1:
        fEpsEqTol = sc.doc.ModelAbsoluteTolerance * 10.**e
        idx_rgFace_Pos = None
        
        # Test one bounding box of each face against bounding box of original brep.
        for f, rgFace in enumerate(rgBrep_toSearch.Faces):
            rgBrep1_1F = rgFace.DuplicateFace(False)
            if rgBrep1_1F is None:
                if bDebug: sPrint = 'rgBrep1_1F'; print sPrint + ':', eval(sPrint)
                return
            
            bBox_SplitF = rgBrep1_1F.GetBoundingBox(True)
            rgBrep1_1F.Dispose()
            if bDebug:
                sc.doc.Objects.AddBrep(rg.Brep.CreateFromBox(bBox_SplitF))
                sc.doc.Views.Redraw()
            
            if epsilonEqualsBoundingBox(rgBbox_B0, bBox_SplitF, fEpsEqTol): # From geometry.py.
                if idx_rgFace_Pos is None: idx_rgFace_Pos = f
                else:
                    idx_rgFace_Pos = None
                    break
        
        if idx_rgFace_Pos is not None: return idx_rgFace_Pos


def usingBoundingBoxOfEdges(rgBrep_toSearch, rgBrep0_1F, bDebug=False):
    if bDebug: print 'xBrep_findMatchingFace.usingBoundingBoxOfEdges()...'
    
    rgCrvs_B0Edges = rgBrep0_1F.DuplicateEdgeCurves() # Includes seams
    
    rgBbox_B0Crvs = None
    for rgCrv in rgCrvs_B0Edges:
        rgBbox = rgCrv.GetBoundingBox(True)
        rgCrv.Dispose()
        if rgBbox_B0Crvs is None: rgBbox_B0Crvs = rgBbox
        else: rgBbox_B0Crvs.Union(rgBbox)
    
    for e in (0,1):
        fEpsEqTol = sc.doc.ModelAbsoluteTolerance * 10.**e
        idx_rgFace_Pos = None
        
        for f, rgFace in enumerate(rgBrep_toSearch.Faces):
            rgBrep1_1F = rgFace.DuplicateFace(False)
            rgFace.Dispose()
            if rgBrep1_1F is None:
                if bDebug: sPrint = 'rgBrep1_1F'; print sPrint + ':', eval(sPrint)
                return
            
            rgCrvs_B1Edges = rgBrep1_1F.DuplicateEdgeCurves()
            
            rgBrep1_1F.Dispose()
            
            """ Test one bounding box for each complete set of edges of each face
            against the bounding box of the splitting curves."""
            bBox_SplitEs = None
            for rgCrv in rgCrvs_B1Edges:
                rgBbox = rgCrv.GetBoundingBox(True)
                rgCrv.Dispose()
                if bBox_SplitEs is None: bBox_SplitEs = rgBbox
                else: bBox_SplitEs.Union(rgBbox)
            
            if epsilonEqualsBoundingBox(rgBbox_B0Crvs, bBox_SplitEs, fEpsEqTol): # From geometry.py.
                if idx_rgFace_Pos is None: idx_rgFace_Pos = f
                else:
                    idx_rgFace_Pos = None
                    break
        
        if idx_rgFace_Pos is not None: return idx_rgFace_Pos


def usingPointOnFace(rgBrep_toSearch, ptOnFace, bDebug=False):
    """
    Returns face index.
    """
    for rgFace in rgBrep_toSearch.Faces:
        b, u, v = rgFace.ClosestPoint(ptOnFace)
        if b:
            print rgFace.IsPointOnFace(u, v)
            if rgFace.IsPointOnFace(u, v):
                return rgFace.FaceIndex
