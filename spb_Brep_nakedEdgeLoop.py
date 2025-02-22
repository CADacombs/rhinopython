"""
160924: Created.
161003: Created curves are now selected when script ends.
171227: Name changed from dupBorderSelect to duplicateBorderSelect.
        Modularizations.
180101: Now all edges are returned from rgCurvesOfBorders when Curve.JoinCurves doesn't create a single, closed curve.
180915-16: Refactored and externalized some functions to another function.
        Name changed from duplicateBorderSelect.py to dupBorderSelect.py.
191110-12: Returned some functions from another module.  Added more functions.
200109: Modified a function name.
200205: Refactored 2 functions into 1.
200602: Modified some variable names.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

#
# Functions only used externally.


def getNakedEdgeIndices_2EdgesPerVertexOnly(brep):


    def getBorderEdgeInfo_OnlyVertexCountsOf2(brep):

        idxs_NE_All = []
        idxs_NEVs_All = []

        for edge in brep.Edges:
            if edge.Valence == rg.EdgeAdjacency.Naked:
                idxs_NE_All.append(edge.EdgeIndex)
                idxs_NEVs_All.append(
                    (edge.StartVertex.VertexIndex,
                        edge.EndVertex.VertexIndex))

        # Edges of naked borders normally have a total of 2 count for each vertex.
        # Edges with vertex counts != 2 will be ignored.
        flatList = [item for pair in idxs_NEVs_All for item in pair]
        idxVs_not2Count = [item for item in set(flatList) if flatList.count(item) != 2]

        if not idxVs_not2Count:
            return idxs_NE_All, idxs_NEVs_All
        else:
            idxs_NE_2CtVs, idxs_NEVs_2CtVs = [], []
            for idxE, idxVs in zip(idxs_NE_All, idxs_NEVs_All):
                edge = brep.Edges[idxE]
                if (
                    (edge.PointAtStart not in idxVs_not2Count) and
                    (edge.PointAtEnd not in idxVs_not2Count)
                ):
                    idxs_NE_2CtVs.append(idxE)
                    idxs_NEVs_2CtVs.append(idxVs)
            return idxs_NE_2CtVs, idxs_NEVs_2CtVs


    idx_NEs, idxs_Vs  = getBorderEdgeInfo_OnlyVertexCountsOf2(brep)

    idxs_nakedEdges_perBorder = [[idx_NEs[0]]]

    while True:
        sc.escape_test()
            
        if not idxs_nakedEdges_perBorder[-1]:
            # Add first items in nested lists.
            for i in range(1, len(idx_NEs)):
                idxE_A = idx_NEs[i]
                if not any(idxE_A in idxs_Border for idxs_Border in idxs_nakedEdges_perBorder):
                    idxs_nakedEdges_perBorder[-1].append(idxE_A)
                    break
        else:
            idxE_A = idxs_nakedEdges_perBorder[-1][-1]
            
        for i in range(1, len(idx_NEs)):
            idxE_B = idx_NEs[i]
            if any(idxE_B in idxs_Border for idxs_Border in idxs_nakedEdges_perBorder):
                continue
            if not idxs_nakedEdges_perBorder[-1]:
                # This only occurs after the 'else' in 'else' below.
                idxs_nakedEdges_perBorder[-1].append(idxE_B)
                continue
            if (
                    idxs_Vs[idx_NEs.index(idxE_A)][0] in idxs_Vs[i] or
                    idxs_Vs[idx_NEs.index(idxE_A)][1] in idxs_Vs[i]
            ):
                idxs_nakedEdges_perBorder[-1].append(idxE_B)
                break # to restart 'for i' but with the last added Edge index as Edge A.
        else:
            if sum([len(border) for border in idxs_nakedEdges_perBorder]) == len(idx_NEs):
                # All naked edge indices are in idxs_nakedEdges_perBorder.
                break # out of 'while'.
            else:
                # Start a new set of border edges.
                idxs_nakedEdges_perBorder.append([])

    #print idxs_nakedEdges_perBorder

    return idxs_nakedEdges_perBorder


def createClosedCrvs_2EdgesPerVertexOnly(brep, bDebug=False):
    """
    Returns closed curves.
    """



    def joinEdgeCrvsToEdgeTols(brep, idxs_Border, bDebug=False):

        # Do not get EdgeCurves since they sometimes extend beyond the Edge Vertices.
        crvs_toJoin = []
        tols_Border = []
        for iE in idxs_Border:
            edge = brep.Edges[iE]
            crvs_toJoin.append(edge)
            tols_Border.append(edge.Tolerance)

        rgCrvs_Joined = rg.Curve.JoinCurves(
                crvs_toJoin,
                joinTolerance=2.01*max(tols_Border))
        if not rgCrvs_Joined:
            1/0
        elif len(rgCrvs_Joined) > 1:
            if bDebug:
                s  = "More than 1 closed curve joined for a single border."
                s += "  This brep will not be used for trimming."
                print s
            return

        rgCrv_Joined = rgCrvs_Joined[0]

        #if not rgCrv_Joined.IsClosed:
        #    sc.doc.Objects.AddCurve(rgCrv_Joined)
        #    raise ValueError("Warning... curves for JoinCurves did not close.")

        return rgCrv_Joined



    idxs_NEs_perBorder = getNakedEdgeIndices_2EdgesPerVertexOnly(brep)
    crvs_Border = []
    for idxs_Border in idxs_NEs_perBorder:
        crv_Border = joinEdgeCrvsToEdgeTols(brep, idxs_Border, bDebug=bDebug)
        if crv_Border and crv_Border.IsClosed: crvs_Border.append(crv_Border)
    return crvs_Border


#
#


def convertObjrefsToSortedBrepsAndEdges(objrefs):
    """
    Parameters:
        objrefs of Edges
    Returns on success:
        list(GUIDs of breps)
        list(list(Edge indices) per brep)
    """
    
    gBreps = []
    idxs_Edges = []
    
   # Organize input into synchronized lists of brep ids and edge indices.
    for objref in objrefs:
        
        gBrep = objref.ObjectId
        idx_Edge = objref.Edge().EdgeIndex
        
        if gBrep not in gBreps:
            # Brep not in list, so add it as well as objref's edge.
            gBreps.append(gBrep)
            idxs_Edges.append([idx_Edge])
        else:
            # Brep already in list, so add the objref's edge to the relative index.
            iB = gBreps.index(gBrep)
            if idx_Edge not in idxs_Edges[iB]:
                idxs_Edges[iB].append(idx_Edge)
    
    return gBreps, idxs_Edges


def getInput():
    
    # Get trimmed edges of breps.
    go = ri.Custom.GetObject()

    go.GeometryFilter = Rhino.DocObjects.ObjectType.EdgeFilter
    go.GeometryAttributeFilter = (ri.Custom.GeometryAttributeFilter.BoundaryEdge)

    go.SetCommandPrompt("Select at least one edge of each naked edge loop")

    go.EnablePreSelect(False, True)

    res = go.GetMultiple(1, 0)

    if res != ri.GetResult.Object: return
    
    objrefs = go.Objects()
    go.Dispose()

    return objrefs


def getEdgeIndicesOfBorders(rhBrep, idxs_Edges, bDebug=False):
    """
    Parameter:
        rhBrep
        idxs_Edges: List of BrepEdge indices or None if all borders (naked edges) are to be returned.
    Loop through each target edge of the brep:
        Find adjacent naked edges using vertices.
        Accumulate edges in a list of a naked edge loop.
    Join edges of naked edge loop into a closed curve.
    
    Returns: 1-level deep list: Indices of BrepEdge in list per border.
        The edges are in order in border loop.
    """

    rgBrep = rs.coercebrep(rhBrep)
    if not rgBrep.IsValid:
        print "Warning!  Brep is not valid.  Check results."
    
    idxs_Es_perBorder = [] # 1-level deep list.
    
    if not idxs_Edges:
        idxs_Edges = [edge.EdgeIndex for edge in rgBrep.Edges if edge.Valence == rg.EdgeAdjacency.Naked]
    
    for idx_Edge in idxs_Edges:
        if bDebug:
            print '-'*80
            sEval = 'idx_Edge'; print sEval+':', eval(sEval)
        
        # Skip edges that are already in join list.
        if idx_Edge in [i for b in idxs_Es_perBorder for i in b]:
            continue # To next edge.
        
        rgEdge = rgBrep.Edges[idx_Edge]
        idx_rgVertex_Start0 = rgEdge.StartVertex.VertexIndex
        if bDebug: sEval = 'idx_rgVertex_Start0'; print sEval+':', eval(sEval)
        idx_rgVertex_Last = None
        if bDebug: sEval = 'idx_rgVertex_Last'; print sEval+':', eval(sEval)
        idxs_Edges_in_border = []
        
        # Accumulate indices of a border of edges.
        while True:
            if bDebug: print '-'*40
            sc.escape_test()
            idxs_Edges_in_border.append(idx_Edge)
            rgVertex_Now = rgEdge.EndVertex
            idx_Vertex_Now = rgVertex_Now.VertexIndex
            if bDebug: sEval = 'idx_Vertex_Now'; print sEval+':', eval(sEval)
            
            # When does this happen?
            if bDebug: sEval = 'idx_Vertex_Now == idx_rgVertex_Last'; print sEval+':', eval(sEval)
            if idx_Vertex_Now == idx_rgVertex_Last:
                if bDebug:
                    print "End vertex matches the last end vertex, so using the start vertex instead."
                rgVertex_Now = rgEdge.StartVertex
                if bDebug: sEval = 'rgVertex_Now'; print sEval+':', eval(sEval)
                idx_Vertex_Now = rgVertex_Now.VertexIndex
                if bDebug: sEval = 'idx_Vertex_Now'; print sEval+':', eval(sEval)
            
            # Break out of while loop if this is the last edge of naked edge border.
            if bDebug: sEval = 'idx_Vertex_Now == idx_rgVertex_Start0'; print sEval+':', eval(sEval)
            if idx_Vertex_Now == idx_rgVertex_Start0:
                break
            
            if bDebug: sEval = 'rgVertex_Now.EdgeIndices()'; print sEval+':', eval(sEval)
            # Find next adjacent naked edge.
            for idx_Edge_at_Vertext in rgVertex_Now.EdgeIndices():
                if bDebug: sEval = 'idx_Edge_at_Vertext'; print sEval+':', eval(sEval)
                rgEdge = rgBrep.Edges[idx_Edge_at_Vertext]
                if bDebug: sEval = 'idx_Edge_at_Vertext != idx_Edge'; print sEval+':', eval(sEval)
                if bDebug: sEval = 'idx_Edge_at_Vertext not in idxs_Edges_in_border'; print sEval+':', eval(sEval)
                if bDebug: sEval = 'rgEdge.Valence == rg.EdgeAdjacency.Naked'; print sEval+':', eval(sEval)
                if (
                        idx_Edge_at_Vertext != idx_Edge and
                        idx_Edge_at_Vertext not in idxs_Edges_in_border and
                        rgEdge.Valence == rg.EdgeAdjacency.Naked
                ):
                    idx_Edge = idx_Edge_at_Vertext
                    if bDebug: sEval = 'idx_Edge'; print sEval+':', eval(sEval)
                    idx_rgVertex_Last = idx_Vertex_Now
                    if bDebug: sEval = 'idx_rgVertex_Last'; print sEval+':', eval(sEval)
                    break
            else: print "Next adjacent naked edge not found!"
        
        idxs_Es_perBorder.append(idxs_Edges_in_border)
    
    return idxs_Es_perBorder


def get_joined_curves(rhBrep, idxs_Edges_per_border, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBrep: rg.Brep, rd.BrepObject, or (GUID of brep)
        idxs_Edges_per_border: These are ordered to border loop.
        bEcho,
        bDebug
        
    Returns on success:
        list(rg.Curve (0-level deep))
    """
    
    rgBrep = rs.coercebrep(rhBrep)
    
    rgCrvs_Borders = []
    
    # Notice that joinTolerance is not provided to JoinCurves.
    
    for edgesOfBorder in idxs_Edges_per_border:
        rgEdges_in_border = [rgBrep.Edges[idx] for idx in edgesOfBorder]
        rgCrvs = rg.Curve.JoinCurves(rgEdges_in_border)
        rgCrvs_Borders.extend(rgCrvs)
    
    return rgCrvs_Borders


def add_CurvesOfBorders(rgCrvs, bEcho=True, bDebug=False):
    
    gCurves_borders = []
    
    for rgCrv in rgCrvs:
        if not rgCrv.IsClosed: bAllCurvesClosed = False
        gCrv = sc.doc.Objects.AddCurve(rgCrv)
        if gCrv == Guid.Empty: continue
        sc.doc.Objects.Select(gCrv)
        gCurves_borders.append(gCrv)
    
    return gCurves_borders


def processBrep(rgBrep, idxs_Trims, bEcho=True, bDebug=False):
    
    idxs_Edges_per_border = getEdgeIndicesOfBorders(rgBrep, idxs_Trims, bDebug)
    if idxs_Edges_per_border is None: return
    
    rgCrvs_borders = get_joined_curves(rgBrep, idxs_Edges_per_border, bEcho, bDebug)
    if rgCrvs_borders is None: return
    
    gCurves_borders = add_CurvesOfBorders(rgCrvs_borders, bEcho, bDebug)
    return gCurves_borders


def processBrepObjects(rhBreps, idx_Trims_PerBrep, bEcho=True, bDebug=False):
    """
    Parameters:
        rdBreps, # list(rd.BrepObject's or GUID's of them)
        idx_Trims_PerBrep,
        bEcho,
        bDebug
    Returns on success:
        list(list(GUIDs of curves) per brep))
    """
    
    gCurves_Borders_PerBrep = [] # 1-level deep list nesting.
    
    for rhBrep, idxs_Trims in zip(rhBreps, idx_Trims_PerBrep):
        rgBrep_In = rs.coercebrep(rhBrep)
        gCurves_Borders_PerBrep.append(
            processBrep(rgBrep_In, idxs_Trims, bEcho, bDebug))
    
    return gCurves_Borders_PerBrep


def main(bEcho=True, bDebug=False):
    
    objrefs = getInput()
    if not objrefs: return

    rc = convertObjrefsToSortedBrepsAndEdges(objrefs)
    if rc is None: return
    
    gBreps, idxs_Edges_per_Brep = rc
    
    sc.doc.Objects.UnselectAll()
    
    gCurves_Borders_PerBrep = processBrepObjects(
        gBreps, idxs_Edges_per_Brep, bEcho, bDebug)
    if gCurves_Borders_PerBrep is None: return
    
    sc.doc.Views.Redraw()
    
    gCurve_ct = 0
    iOpenCurve_ct = 0
    for gCurves_borders in gCurves_Borders_PerBrep:
        for gCurve in gCurves_borders:
            if not rs.IsCurveClosed(gCurve): iOpenCurve_ct += iOpenCurve_ct
            gCurve_ct += 1
    
    if not iOpenCurve_ct:
        print "{} curve(s) added to document, and all are closed.".format(gCurve_ct)
    else:
        print "{} curve(s) added to document, and {} are open.".format(
                gCurve_ct, iOpenCurve_ct)


if __name__ == '__main__': main(bEcho=bool(1), bDebug=bool(0))
