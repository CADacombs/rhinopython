"""
160406: Created.
161212: Commented out alias maker.
170501: Added sticky for fAngleTol_Deg.  Removed alias maker.
170507: Modified command option text.
171104: Split from 1 function into 3.
        Now waits for input even when objects are preselected.
        Can now escape out of face loop.
171108: Modified printed output.
180317: Added DotFontHt option and if executed in V6 or greater, default dot font height will be 3* that in V5.
180427: Renamed from findSharpEdges.py to markSharpEdges.py.
180725: Added concavity option.
        Renamed from markSharpEdges.py to sharpEdges.py.
181025: Bug fix in printed output.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


def getInput():
    
    sRhFileName = '{}'.format(sc.doc.Name) # If Name doesn't exist, an empty string is created.
    
    fAngleTol_Deg_Default = 50.0 * sc.doc.ModelAngleToleranceDegrees
    idxConcavityType_Default = idxConcavityType = 2 # 2 = Either concave or convex.
    iDotHeight_Default = 33 if Rhino.RhinoApp.ExeVersion >= 6 else 11
    
    # Create options with default values before any are changed with those in sticky.
    optD_fAngleTol_Deg = ri.Custom.OptionDouble(
            fAngleTol_Deg_Default, True, Rhino.RhinoMath.ZeroTolerance)
    sConcavityTypes = 'ConcaveOnly', 'ConvexOnly', 'EitherConcaveOrConvex'
    optT_bDupEdge = ri.Custom.OptionToggle(False, 'No', 'Yes')
    optT_bAddDot = ri.Custom.OptionToggle(True, 'No', 'Yes')
    optI_iDotFontHt = ri.Custom.OptionInteger(iDotHeight_Default, True, 3)
    
    # Load sticky.
    stickyKeys = [
            'fAngleTol_Deg({})'.format(__file__),
            'idxConcavityType({})'.format(__file__),
            'bDupEdge({})'.format(__file__),
            'bAddDot({})({})'.format(__file__, sRhFileName),
            'iDotFontHt({})'.format(__file__),
    ]
    stickyValues = [
            optD_fAngleTol_Deg.CurrentValue,
            idxConcavityType_Default,
            optT_bDupEdge.CurrentValue,
            optT_bAddDot.CurrentValue,
            optI_iDotFontHt.CurrentValue,
    ]
    for i, stickyKey in enumerate(stickyKeys):
        if sc.sticky.has_key(stickyKey): stickyValues[i] = sc.sticky[stickyKey]
    (
            optD_fAngleTol_Deg.CurrentValue,
            idxConcavityType,
            optT_bDupEdge.CurrentValue,
            optT_bAddDot.CurrentValue,
            optI_iDotFontHt.CurrentValue,
    ) = stickyValues
    
    # Get objects with optional input
    
    dictOptListIdxs = {'idxConcavityType': None}
    
    def addOptions():
        go.AddOptionDouble("AngleTolerance" + chr(176), optD_fAngleTol_Deg)
        dictOptListIdxs['idxConcavityType'] = go.AddOptionList(
                englishOptionName='ConcavityType',
                listValues=sConcavityTypes,
                listCurrentIndex=idxConcavityType)
        go.AddOptionToggle('DuplicateEdges', optT_bDupEdge)
        go.AddOptionToggle('AddDots', optT_bAddDot)
        if optT_bAddDot.CurrentValue:
            go.AddOptionInteger('DotFontHt', optI_iDotFontHt)
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select breps")
    
    go.GeometryFilter = rd.ObjectType.PolysrfFilter
    
    addOptions()
    
    go.AcceptNumber(True, True)
    
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    bPreselectedObjsChecked = False
    
    while True:
        if sc.escape_test(False): return
        
        bDupEdge_Old = optT_bDupEdge.CurrentValue
        bAddDot_Old = optT_bAddDot.CurrentValue
        
        res = go.GetMultiple(1, 0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return # Esc key was pressed.
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            optD_fAngleTol_Deg.CurrentValue = abs(go.Number())
        
        elif res == ri.GetResult.Option and go.Option().Index == 2:
            idxConcavityType = go.Option().CurrentListOptionIndex
        
        if optD_fAngleTol_Deg.CurrentValue == 0.0:
                optD_fAngleTol_Deg.CurrentValue = fAngleTol_Deg_Default
        
        # Do not allow both bDupEdge and bAddDot to be False.
        if optT_bDupEdge.CurrentValue != bDupEdge_Old:
            optT_bAddDot.CurrentValue = True
        elif optT_bAddDot.CurrentValue != bAddDot_Old:
            optT_bDupEdge.CurrentValue = True
        
        # Save sticky
        stickyValues = (
                optD_fAngleTol_Deg.CurrentValue,
                idxConcavityType,
                optT_bDupEdge.CurrentValue,
                optT_bAddDot.CurrentValue,
                optI_iDotFontHt.CurrentValue,
        )
        for i, stickyKey in enumerate(stickyKeys):
            sc.sticky[stickyKey] = stickyValues[i]
        
        go.ClearCommandOptions()
        addOptions()
    
    return (
            [objref.ObjectId for objref in go.Objects()],
            optD_fAngleTol_Deg.CurrentValue,
            idxConcavityType,
            optT_bDupEdge.CurrentValue,
            optT_bAddDot.CurrentValue,
            optI_iDotFontHt.CurrentValue,
    )


def sharpEdges(brep, fAngleTol_Deg=sc.doc.ModelAngleToleranceDegrees, idxConcavityType=2):
    """
    Parameters:
        brep: Can be Object document or geometry.
        fAngleTol_Deg: Used to determine how many iterations to use to refine fAngle.
        idxConcavityType: 0=Concave only, 1=Convex only , 2=Either concave or convex
    Returns:
        List of indices of sharp BrepEdges of desired concavity in brep.
    """
    
    rgBrep = brep if isinstance(brep, rg.Brep) else rs.coercebrep(brep)
    
    idx_rgEdges_Sharp = []
    gCrvs1 = [] # Accumulation of duplicated edges (curves)
    
    fAngleTol_Rad = Rhino.RhinoMath.ToRadians(fAngleTol_Deg)
    
    decimalFractions = 0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.0, 1.0
    
    # Accumulate indices of sharp edges.
    for rgEdge in rgBrep.Edges:
        if rgEdge.Valence != rg.EdgeAdjacency.Interior:
            # Skip edge because it doesn't have 2 faces.
            continue
        
        if rgEdge.IsSmoothManifoldEdge(fAngleTol_Rad): continue
        
        if idxConcavityType != 2:
            # Test for concavity vs. convexity against desired state.
            
            tMax = rgEdge.Domain.Max
            tMin = rgEdge.Domain.Min
            
            for d in decimalFractions:
                t = (tMax-tMin) * d + tMin
                
                concavity = rgEdge.ConcavityAt(
                        t=t,
                        tolerance=fAngleTol_Rad)
                if idxConcavityType == 0 and concavity == rg.Concavity.Concave:
                    # Concave found.
                    break
                elif idxConcavityType == 1 and concavity == rg.Concavity.Convex:
                    # Convex found.
                    break
            else:
                # Desired concavity not found, so continue onto next edge.
                continue
        
        idx_rgEdges_Sharp.append(rgEdge.EdgeIndex)
    
    return idx_rgEdges_Sharp


def main(bDebug=False):
    
    rc = getInput()
    if rc is None: return
    (
            gBreps,
            fAngleTol_Deg,
            idxConcavityType,
            bDupEdge,
            bAddDot,
            iDotFontHt,
    ) = rc
    
    nSharp_All = iBrepWithoutSharps_ct = 0
    
    if bDupEdge: sc.doc.Objects.UnselectAll()
    
    if idxConcavityType == 0:
        sConcavityPrint = 'some concavity'
    elif idxConcavityType == 1:
        sConcavityPrint = 'some convexity'
    else:
        sConcavityPrint = 'both concavity and convexity'
    
    print "Searching for sharp edges with {} within {}{}...".format(
            sConcavityPrint,
            fAngleTol_Deg,
            chr(176)
    )
    
    # Process breps.
    for gBrep in gBreps:
        rgBrep = rs.coercebrep(gBrep)
        
        rgEdges_All = rgBrep.Edges
        
        rc = sharpEdges(rgBrep, fAngleTol_Deg, idxConcavityType)
        if rc is None: return
        idx_rgEdges_Sharp = rc
        
        if len(idx_rgEdges_Sharp) == 0:
            # No sharps of desired concavity found in this brep.
            iBrepWithoutSharps_ct += 1
            continue
        
        for idx in idx_rgEdges_Sharp:
            rgEdge = rgBrep.Edges[idx]
            if bDupEdge:
                sc.doc.Objects.Select(sc.doc.Objects.AddCurve(rgEdge))
            if bAddDot:
                t = rgEdge.DivideByCount(segmentCount=2, includeEnds=False)[0]
                rgDot = rg.TextDot(str(''), rgEdge.PointAt(t))
                rgDot.FontHeight = iDotFontHt
                sc.doc.Objects.AddTextDot(rgDot)
        
        print "{} sharp, interior edges with {} found in {}.".format(
                len(idx_rgEdges_Sharp),
                sConcavityPrint,
                gBrep)
    
    if iBrepWithoutSharps_ct == len(gBreps):
        print "No sharps with {} found in {} polyface breps.".format(
                sConcavityPrint,
                iBrepWithoutSharps_ct)
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main(bDebug=bool(0))