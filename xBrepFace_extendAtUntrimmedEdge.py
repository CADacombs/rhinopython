"""
160827-29: Created.
160907-08: Now outer and inner trims are obtained separately and the former are
        extended before the extended surface is trimmed.
        Now will trim non-selected edges of same SENW.
        Modularizations.
        Although only one trim is allowed to be extended when running this script
        directly, multiple extended Trims can be passed to trimExtendedSrf.
160915: trimExtendedSrf: Fixed bug in IsoStatus matching.
190413,0423, 200109, 0619: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import xBrep_findMatchingFace
import xBrepFace
import xBrepTrim


def dispose(listDispose): map(lambda x: x.Dispose, listDispose)


def getBrepIdAndTrimIdx():
    
    disposeUs = []
    
    # Load sticky.
    stickyKeys = ['bSmooth({})'.format(__file__)]
    stickyValues = [True]
    for i, stickyKey in enumerate(stickyKeys):
        if sc.sticky.has_key(stickyKey): stickyValues[i] = sc.sticky[stickyKey]
    bSmooth, = stickyValues
    
    # Get untrimmed brep edge with optional input.
    go = ri.Custom.GetObject(); disposeUs.append(go)
    go.SetCommandPrompt("Select untrimmed edge to extend")
    go.GeometryFilter = rd.ObjectType.EdgeFilter
    # SurfaceBoundaryEdge doesn't include seams.
    go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.SurfaceBoundaryEdge)
    
    def is1FaceBrep (rdBrep, rgTrim, compIdx):
        return rdBrep.BrepGeometry.Faces.Count == 1
    go.SetCustomGeometryFilter(is1FaceBrep)
    
    optT_Smooth = ri.Custom.OptionToggle(bSmooth, 'RuledExtension', 'Smooth')
    go.AddOptionToggle('Type', optT_Smooth)
    
    while go.Get() != ri.GetResult.Object:
        sc.escape_test()
        if go.CommandResult() != Rhino.Commands.Result.Success: # Canceled.
            dispose(disposeUs); return
    
    objref = go.Object(0); disposeUs.append(objref)
    idBrep = objref.ObjectId
    compIdx = objref.GeometryComponentIndex.Index
    
    bSmooth = optT_Smooth.CurrentValue
    
    # Save sticky.
    stickyValues = bSmooth,
    for i, stickyKey in enumerate(stickyKeys):
        sc.sticky[stickyKey] = stickyValues[i]
    
    dispose(disposeUs)
    return idBrep, compIdx, bSmooth


def getExtensionFactor():
    
    # Load sticky.
    stickyKeys = ['fExtLength({})'.format(__file__)]
    stickyValues = [1.0]
    for i, stickyKey in enumerate(stickyKeys):
        if sc.sticky.has_key(stickyKey): stickyValues[i] = sc.sticky[stickyKey]
    fExtLength, = stickyValues
    
    rc, fExtLength = ri.RhinoGet.GetNumber("ExtensionFactor",
            False, fExtLength, 0., 1./Rhino.RhinoMath.ZeroTolerance)
    
    if rc == Rhino.Commands.Result.Success:
        # Save sticky.
        stickyValues = fExtLength,
        for i, stickyKey in enumerate(stickyKeys):
            sc.sticky[stickyKey] = stickyValues[i]
        
        return fExtLength


def isBrepReadyForExtend(rgBrep, bEcho=False, bDebug=False):
    
    if not rgBrep.IsValid:
        if bEcho: print "Surface is not valid.  Exiting..."
        return False
    if not rgBrep.IsManifold:
        if bEcho: print "Surface is non-manifold.  Exiting..."
        return False
    if rgBrep.IsSurface and rgBrep.IsSolid:
        if bEcho: print "Surface is a solid.  Exiting..."
        return False
    
    disposeUs = []
    
    rgFace = rgBrep.Faces[0]; disposeUs.append(rgFace)
    if not rgFace.IsValid:
        if bEcho:
            print "Warning: Face is invalid, but script will attempt to continue."
    
    rgSrf = rgFace.UnderlyingSurface(); disposeUs.append(rgSrf)
    if not rgSrf.IsValid:
        if bEcho: print "Underlying surface is invalid.  Exiting.."
        dispose(disposeUs); return 
    
    rgBrep_Extended_TrySplit = rgSrf.ToBrep(); disposeUs.append(rgBrep_Extended_TrySplit)
    
    fKinkTol = sc.doc.ModelAngleToleranceRadians
    bSplitKinkyFaces = rgBrep_Extended_TrySplit.Faces.SplitKinkyFaces(
            fKinkTol, True) # True doesn't mean that any splits had occurred.
    if not bSplitKinkyFaces:
        if bDebug: sPrint = 'bSplitKinkyFaces'; print sPrint + ':', eval(sPrint)
        if bEcho: print "Split kiny face check failed.  Exiting..."
        dispose(disposeUs); return
    if rgBrep_Extended_TrySplit.Faces.Count != 1:
        if bEcho: print "At {}{} tolerance, the surface has kinks.  " \
                "Repair surface before extending it.".format(
                sc.doc.ModelAngleToleranceDegrees, chr(176))
        dispose(disposeUs); return
    
    dispose(disposeUs)
    
    return True


def getAreaOfSurface(rgSrf, bEcho=False):
    """Tries brep if surface fails."""
    areaMassProp = rg.AreaMassProperties.Compute(rgSrf)
    if areaMassProp is None:
        if bEcho:
            print "AreaMassProperties cannot be computed for surface."
            print "Trying to compute AreaMassProperties for brep of surface..."
        rgBrepTemp = rgSrf.ToBrep()
        areaMassProp = rg.AreaMassProperties.Compute(rgBrepTemp)
        rgBrepTemp.Dispose()
        if areaMassProp is None: return
        elif bEcho: print "...Success!"
    area = areaMassProp.Area
    areaMassProp.Dispose()
    return area


def isExtendedSrfOk(rgSrf_Ext, rgSrf0, bEcho=False, bDebug=False):
    """
    Includes verification that area has changed.
    """
    
    if not rgSrf_Ext.IsValid:
        return False
    
    # Get AreaMassProperties of original and new surface.
    areaSrf0 = getAreaOfSurface(rgSrf0)
    if areaSrf0 is None:
        if bEcho: print "Area not computed for original surface."
        return False
    areaSrf_Extended = getAreaOfSurface(rgSrf_Ext)
    if areaSrf_Extended is None:
        if bEcho: print "Area of original surface not computed."
        return False
    
    # Compare original and new areas and, if the same, return.
    if Rhino.RhinoMath.EpsilonEquals(areaSrf_Extended, areaSrf0,
            sc.doc.ModelAbsoluteTolerance**2):
        return False
    
    return rgSrf_Ext


def trimExtendedSrf(rgSrf_Ext, rgBrep0, idxTrims_toSkip, bEcho=False, bDebug=False):
    """
    Parameters:
        rgSrf_Ext: Extended surface.
        rgBrep0: Original 1-face Brep.
        idxTrims_toSkip: List of SENW Trims not to trim.
    Returns:
        Trimmed 1-face Brep or None.
    """
    
    disposeUs = []
    
    fPtMatchTol = 1.e-6
    
    # Make IsoStatus list per idxTrims_toSkip.
    isoStats_Ext = []
    for i in idxTrims_toSkip:
        rgT_Ext = rgBrep0.Trims[i]; disposeUs.append(rgT_Ext)
        if rgT_Ext.IsoStatus not in isoStats_Ext:
            isoStats_Ext.append(rgT_Ext.IsoStatus)
    
    # Make outer Edge list and inner Edge list of only non-SENW Edges or
    # SENW Trim not extended for Trim.
    rgTrims_B0 = rgBrep0.Trims; disposeUs.extend(rgTrims_B0)
    rgEs_Otr = []; rgEs_Inr = []
    for i, rgTrim in enumerate(rgTrims_B0):
        if rgTrim.Loop.LoopType == rg.BrepLoopType.Outer: # Outer loop.
            # Only use other edges with selected edge SENW.
            if xBrepTrim.isSenw(rgTrim):
                if any(i == idx for idx in idxTrims_toSkip):
                    continue # Trims of idxTrims_toSkip will not be used.
                if any(rgTrim.IsoStatus == iso for iso in isoStats_Ext):
                    continue # This is a SENW different than that of any of the extended Trims.
            rgEs_Otr.append(rgTrim.Edge)
        else: # Inner loop
            rgEs_Inr.append(rgTrim.Edge)
    
    #map(sc.doc.Objects.AddCurve, rgEs_Otr)
    
    rgBrep_Extended = Rhino.Geometry.Brep.CreateFromSurface(rgSrf_Ext)
    
    # How could this happen?
    if rgEs_Otr.Count + rgEs_Inr.Count == 0:
        if bEcho: print "Original brep was trimmed, but " \
                "no edges found with which to trim.  Using full surface."
        return rgBrep_Extended
    
    # Trim face (surface) with both outer and inner edges.
    disposeUs.append(rgBrep_Extended)
    disposeUs.extend(rgEs_Otr + rgEs_Inr)
    
    # TO (Possibly) DO: If extending the curves when they don't need to
    # be becomes a problem, first check whether extended edge is a full
    # surface edge and create conditional based on that condition.
    
    # Join curves (edges).
    rgCs_Otr_Joined_All = rg.Curve.JoinCurves(rgEs_Otr)
    disposeUs.extend(rgCs_Otr_Joined_All)
    
    # Extend joined curves on surface, using SimplifyEnd at actual extensions
    # before the curves are exploded.
    rgFace_Ext = rgBrep_Extended.Faces[0]; disposeUs.append(rgFace_Ext)
    rgCs_Otr_Joined_Ext = []
    for rgC_Otr_Joined in rgCs_Otr_Joined_All:
        rgC_Otr_Joined_Ext = rgC_Otr_Joined.ExtendOnSurface(
                rg.CurveEnd.Both, rgFace_Ext)
        if rgC_Otr_Joined_Ext is None: continue
        disposeUs.append(rgC_Otr_Joined_Ext)
        
        ## Find extension "side" ends by comparing CurveEnd points, before vs. after ExtendOnSurface.
        # At extended "side" ends, use SimplifyEnd to simplify only that "side" end and adjacent segments.
        ptS0 = rgC_Otr_Joined.PointAtStart
        ptE0 = rgC_Otr_Joined.PointAtEnd
        ptS1 = rgC_Otr_Joined_Ext.PointAtStart
        ptE1 = rgC_Otr_Joined_Ext.PointAtEnd
        
        # CurveEnd.Start check and simplification.
        if not ptS1.EpsilonEquals(ptS0, fPtMatchTol):
            rgC_Otr_Joined_Ext_StartSimpl = rgC_Otr_Joined_Ext.SimplifyEnd(
                    rg.CurveEnd.Start, rg.CurveSimplifyOptions.All,
                    fPtMatchTol, sc.doc.ModelAngleToleranceRadians)
            if rgC_Otr_Joined_Ext_StartSimpl is None:
                rgC_Otr_Joined_Ext_StartSimpl = rgC_Otr_Joined_Ext # For CurveEnd.End processing.
            else:
                disposeUs.append(rgC_Otr_Joined_Ext)
        else: rgC_Otr_Joined_Ext_StartSimpl = rgC_Otr_Joined_Ext # For CurveEnd.End processing.
        
        # CurveEnd.End check and simplification.
        if not ptE1.EpsilonEquals(ptE0, fPtMatchTol):
            rgC_Otr_Joined_Ext = rgC_Otr_Joined_Ext_StartSimpl.SimplifyEnd(
                    rg.CurveEnd.End, rg.CurveSimplifyOptions.All,
                    fPtMatchTol, sc.doc.ModelAngleToleranceRadians)
            if rgC_Otr_Joined_Ext is None:
                rgC_Otr_Joined_Ext = rgC_Otr_Joined_Ext_StartSimpl
            else:
                disposeUs.append(rgC_Otr_Joined_Ext)
        rgCs_Otr_Joined_Ext.append(rgC_Otr_Joined_Ext)
    disposeUs.extend(rgCs_Otr_Joined_Ext)
    
    # Explode curves so that the proper vertices are replaced on the brep.
    rgCs_Otr_Ext_Segs_All = []
    for rgCrv_Otr_Ext in rgCs_Otr_Joined_Ext:
        rgCrvs_Otr_Ext_Segs = rgCrv_Otr_Ext.DuplicateSegments()
        rgCs_Otr_Ext_Segs_All.extend(rgCrvs_Otr_Ext_Segs)
    disposeUs.extend(rgCs_Otr_Ext_Segs_All)
    
    if bDebug: map(sc.doc.Objects.AddCurve, rgCs_Otr_Ext_Segs_All + rgEs_Inr)
    
    # Split full surface brep with curves.
    rgBrep_Split = rgFace_Ext.Split(rgCs_Otr_Ext_Segs_All + rgEs_Inr,
            sc.doc.ModelAbsoluteTolerance)
    if rgBrep_Split is None: dispose(disposeUs); return
    disposeUs.append(rgBrep_Split)
    
    # Check whether brep has more than one face before attempting to get correct face.
    if rgBrep_Split.Faces.Count > 1:
        # Get point on face of original brep for face matching.
        ptOnFace = xBrepFace.createPoint3dOnInterior(
            rgBrep0.Faces[0],
            fMinDistFromBorder=10.0*sc.doc.ModelAbsoluteTolerance)
        if ptOnFace is None:
            dispose(disposeUs); return
        
        # Get correct brep face.
        idx_rgFace_Pos = xBrep_findMatchingFace.usingPointOnFace(
                rgBrep_Split, ptOnFace)
        if idx_rgFace_Pos is None:
            dispose(disposeUs); return
        
        rgBrep_ForReplace = rgBrep_Split.Faces[
                idx_rgFace_Pos].DuplicateFace(False)
        disposeUs.append(rgBrep_ForReplace)
    else:
        if bEcho:
            print "Brep was not split.  Replacing original with full surface..."
        rgBrep_ForReplace = rgBrep_Split
    
    return rgBrep_ForReplace


def main(bEcho=False, bDebug=False):
    
    disposeUs = []
    
    sTitle = "Extend Base Surface"
    
    ret = getBrepIdAndTrimIdx()
    if ret is None: return
    idBrep0, idxTrim, bSmooth = ret
    
    rgBrep0 = sc.doc.Objects.Find(idBrep0).BrepGeometry
    disposeUs.append(rgBrep0)
    if not isBrepReadyForExtend(rgBrep0, bEcho=bEcho, bDebug=bDebug):
        dispose(disposeUs); return Rhino.Commands.Result.Failure
    
    rgSrf0 = rgBrep0.Faces[0].UnderlyingSurface(); disposeUs.append(rgSrf0)
    if rgSrf0 is None:
        dispose(disposeUs); return Rhino.Commands.Result.Failure
    
    fExtLength = getExtensionFactor()
    if fExtLength is None: dispose(disposeUs); return
    
    rgTrim = rgBrep0.Trims[idxTrim]; disposeUs.append(rgTrim)
    
    # Create extended surface.
    rgSrf_Ext = rgSrf0.Extend(rgTrim.IsoStatus, fExtLength, bSmooth)
    if rgSrf_Ext is None: dispose(disposeUs); return
    disposeUs.append(rgSrf_Ext)
    if not isExtendedSrfOk(rgSrf_Ext, rgSrf0, bEcho=bEcho, bDebug=bDebug):
        dispose(disposeUs); return
    #sc.doc.Objects.AddSurface(rgSrf_Ext)
    
    # Was brep a full surface?
    if rgBrep0.IsSurface:
        if bEcho: print "Original brep was not trimmed.  Using full surface."
        rgBrep_ForReplace = Rhino.Geometry.Brep.CreateFromSurface(rgSrf_Ext)
    else:
        rgBrep_ForReplace = trimExtendedSrf(rgSrf_Ext, rgBrep0, [idxTrim], bEcho=bEcho, bDebug=bDebug)
        if rgBrep_ForReplace is None: dispose(disposeUs); return
    
    if not sc.doc.Objects.Replace(idBrep0, rgBrep_ForReplace):
        dispose(disposeUs); return Rhino.Commands.Result.Failure
    
    sc.doc.Views.Redraw()
    
    dispose(disposeUs); return Rhino.Commands.Result.Success


if __name__ == '__main__': main(bEcho=bool(1), bDebug=bool(0))