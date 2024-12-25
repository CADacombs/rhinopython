"""
160601: Created by importing functions from other scripts.
160611-12: _
160618: Import-related update.
160721: Shortened sticky code.
160822: Minor changes.
190129-0605: Import-related update.
200107-08: Import-related update and various minor changes.
241224: Import-related update.
"""

import Rhino
import rhinoscriptsyntax as rs
import scriptcontext as sc

import spb_Brep_contiguousFilletFaces


def getChainFilletFaces(bDebug=False):
    
    # Load sticky.
    stickyKeys = [
        'bSameRad({})'.format(__file__),
        'fRadTol({})({})'.format(__file__, sc.doc.Name),
        'bTan({})'.format(__file__),
        'fAngleTol({})({})'.format(__file__, sc.doc.Name),
    ]
    stickyValues = [
        True,
        10.0*rs.UnitAbsoluteTolerance(),
        True,
        10.0*rs.UnitAngleTolerance(),
    ]
    for i, stickyKey in enumerate(stickyKeys):
        if sc.sticky.has_key(stickyKey): stickyValues[i] = sc.sticky[stickyKey]
    (
        bSameRad,
        fRadTol,
        bTan,
        fAngleTol,
    ) = stickyValues
    
    #
    # Get constant radius (rolling ball) fillet face of brep.
    go = Rhino.Input.Custom.GetObject()
    go.SetCommandPrompt("Select starting face")
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface
    go.GeometryAttributeFilter = (
            Rhino.Input.Custom.GeometryAttributeFilter.SubSurface) # Doesn't allow single surfaces to be selected.
    
    optT_SameRad = Rhino.Input.Custom.OptionToggle(bSameRad, 'No', 'Yes')
    optD_RadTol = Rhino.Input.Custom.OptionDouble(
            fRadTol, True, Rhino.RhinoMath.ZeroTolerance)
    optT_Tan = Rhino.Input.Custom.OptionToggle(bTan, 'No', 'Yes')
    optD_AngleTol = Rhino.Input.Custom.OptionDouble(
            fAngleTol, True, Rhino.RhinoMath.ZeroTolerance)
    
    go.AddOptionToggle('SameRad', optT_SameRad)
    if optT_SameRad.CurrentValue:
        go.AddOptionDouble('RadTol', optD_RadTol)
    go.AddOptionToggle('TanContinuity', optT_Tan)
    if optT_Tan.CurrentValue:
        go.AddOptionDouble('AngleTol', optD_AngleTol)
    
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == Rhino.Input.GetResult.Object 
    
    while go.Get() != Rhino.Input.GetResult.Object:
        if sc.escape_test(False): break
        
        res = go.CommandResult()
        
        if res != Rhino.Commands.Result.Success:
            go.Dispose()
            return
        
        bSameRad = optT_SameRad.CurrentValue
        bTan = optT_Tan.CurrentValue
        go.ClearCommandOptions()
        go.AddOptionToggle('SameRadius', optT_SameRad)
        if optT_SameRad.CurrentValue:
            fRadTol = optD_RadTol.CurrentValue
            go.AddOptionDouble('RadiusTolerance', optD_RadTol)
        go.AddOptionToggle('Tangent', optT_Tan)
        if bTan:
            fAngleTol = optD_AngleTol.CurrentValue
            go.AddOptionDouble('AngleTolerance', optD_AngleTol)
    
    bSameRad = optT_SameRad.CurrentValue
    bTan = optT_Tan.CurrentValue
    
    # Save sticky
    stickyValues = bSameRad, fRadTol, bTan, fAngleTol
    for i, stickyKey in enumerate(stickyKeys):
        sc.sticky[stickyKey] = stickyValues[i]
    
    rdORef = go.Object(0)
    
    rgFace = rdORef.Face()
    
    fPlanarTol = rs.UnitAbsoluteTolerance()
    if rgFace.IsPlanar(fPlanarTol):
        print "Planar face selected."
        go.Dispose()
        return
    
    rdBrep0 = rdORef.Object() # Parent
    
    rgBrep0 = rdORef.Brep()
    idxFace_Sel = rdORef.GeometryComponentIndex.Index
    
    rc = spb_Brep_contiguousFilletFaces.indicesOfContiguousFilletFaces(
        rgBrep0,
        idxFace_Sel,
        bSameRad,
        fRadTol,
        bTan,
        fAngleTol)

    if rc is None:
        print "Non-fillet face selected."
        go.Dispose()
        return

    idxFaces_Pass = sorted(rc)
    
    idBrep0 = rdORef.ObjectId
    
    go.Dispose()
    
    return idBrep0, idxFaces_Pass


# Tester
if __name__ == '__main__': print getChainFilletFaces(bDebug=True)