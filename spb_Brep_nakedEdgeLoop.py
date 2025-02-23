"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160924: Created.
...
250221,22: Added command options. Bug fix in printed output.

TODO:
    Move createClosedCrvs_2EdgesPerVertexOnly and getNakedEdgeIndices_2EdgesPerVertexOnly
    into their own script with a UI?
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fJoinTol'; keys.append(key)
    #values[key] = sc.doc.ModelAbsoluteTolerance
    value = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    values[key] = float(format(value, '.0e'))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fMinSegLengthToKeep'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    #value = 1e-3 * Rhino.RhinoMath.UnitScale(
    #    Rhino.UnitSystem.Millimeters,
    #    sc.doc.ModelUnitSystem)
    #values[key] = float(format(value, '.0e'))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)


    for key in keys:
        if key not in names:
            names[key] = key[1:]


    # Load sticky.
    for key in stickyKeys:
        _debug = sc.sticky
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key in ('fJoinTol', 'fMinSegLengthToKeep'):
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


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

    #print(idxs_nakedEdges_perBorder

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
                print(s)
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
    """
    Get naked edges of breps.
    """

    go = ri.Custom.GetObject()

    go.GeometryFilter = Rhino.DocObjects.ObjectType.EdgeFilter
    go.GeometryAttributeFilter = (ri.Custom.GeometryAttributeFilter.BoundaryEdge)

    go.SetCommandPrompt("Select one edge of each naked edge loop")

    go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fJoinTol')
        addOption('fMinSegLengthToKeep')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getEdgeIndicesOfBorders(rhBrep, idxs_Edges, bDebug=False):
    """
    Parameter:
        rhBrep
        idxs_Edges: List of BrepEdge indices or None if all borders (naked edges) are to be returned.
        bDebug: bool
    Loop through each target edge of the brep:
        Find adjacent naked edges using vertices.
        Accumulate edges in a list of a naked edge loop.
    Join edges of naked edge loop into a closed curve.
    
    Returns: 1-level deep list: Indices of BrepEdge in list per border.
        The edges are in order in border loop.
    """

    rgBrep = rhBrep if isinstance(rhBrep, rg.Brep) else rs.coercebrep(rhBrep)
    if not rgBrep.IsValid:
        print("Warning!  Brep is not valid.  Check results.")

    idxs_Es_perBorder = [] # 1-level deep list.

    if not idxs_Edges:
        idxs_Edges = [edge.EdgeIndex for edge in rgBrep.Edges if edge.Valence == rg.EdgeAdjacency.Naked]

    for idx_Edge in idxs_Edges:
        if bDebug:
            print('-'*80)
            sEval = 'idx_Edge'; print(sEval+':', eval(sEval))

        # Skip edges that are already in join list.
        if idx_Edge in [i for b in idxs_Es_perBorder for i in b]:
            continue # To next edge.

        rgEdge = rgBrep.Edges[idx_Edge]
        idx_rgVertex_Start0 = rgEdge.StartVertex.VertexIndex
        if bDebug: sEval = 'idx_rgVertex_Start0'; print(sEval+':', eval(sEval))
        idx_rgVertex_Last = None
        if bDebug: sEval = 'idx_rgVertex_Last'; print(sEval+':', eval(sEval))
        idxs_Edges_in_border = []

        # Accumulate indices of a border of edges.
        while True:
            if bDebug: print('-'*40)
            sc.escape_test()
            idxs_Edges_in_border.append(idx_Edge)
            rgVertex_Now = rgEdge.EndVertex
            idx_Vertex_Now = rgVertex_Now.VertexIndex
            if bDebug: sEval = 'idx_Vertex_Now'; print(sEval+':', eval(sEval))

            # When does this happen?
            if bDebug: sEval = 'idx_Vertex_Now == idx_rgVertex_Last'; print(sEval+':', eval(sEval))
            if idx_Vertex_Now == idx_rgVertex_Last:
                if bDebug:
                    print("End vertex matches the last end vertex, so using the start vertex instead.")
                rgVertex_Now = rgEdge.StartVertex
                if bDebug: sEval = 'rgVertex_Now'; print(sEval+':', eval(sEval))
                idx_Vertex_Now = rgVertex_Now.VertexIndex
                if bDebug: sEval = 'idx_Vertex_Now'; print(sEval+':', eval(sEval))
            
            # Break out of while loop if this is the last edge of naked edge border.
            if bDebug: sEval = 'idx_Vertex_Now == idx_rgVertex_Start0'; print(sEval+':', eval(sEval))
            if idx_Vertex_Now == idx_rgVertex_Start0:
                break
            
            if bDebug: sEval = 'rgVertex_Now.EdgeIndices()'; print(sEval+':', eval(sEval))
            # Find next adjacent naked edge.
            
            idx_Es_V_Now = rgVertex_Now.EdgeIndices()
            
            for idx_Edge_at_Vertext in idx_Es_V_Now:
                if bDebug: sEval = 'idx_Edge_at_Vertext'; print(sEval+':', eval(sEval))
                rgEdge = rgBrep.Edges[idx_Edge_at_Vertext]
                
                if bDebug:
                    sEval = 'idx_Edge_at_Vertext != idx_Edge'; print(sEval+':', eval(sEval))
                    sEval = 'idx_Edge_at_Vertext not in idxs_Edges_in_border'; print(sEval+':', eval(sEval))
                    sEval = 'rgEdge.Valence == rg.EdgeAdjacency.Naked'; print(sEval+':', eval(sEval))
                
                if (
                        idx_Edge_at_Vertext != idx_Edge and
                        idx_Edge_at_Vertext not in idxs_Edges_in_border and
                        rgEdge.Valence == rg.EdgeAdjacency.Naked
                ):
                    idx_Edge = idx_Edge_at_Vertext
                    if bDebug: sEval = 'idx_Edge'; print(sEval+':', eval(sEval))
                    idx_rgVertex_Last = idx_Vertex_Now
                    if bDebug: sEval = 'idx_rgVertex_Last'; print(sEval+':', eval(sEval))
                    break
            else: print("Next adjacent naked edge not found!")

        idxs_Es_perBorder.append(idxs_Edges_in_border)

    return idxs_Es_perBorder


def removeIndicesOfShortEdges(rhBrep, idxs_Edges_in_border, fMinSegLengthToKeep=None):
    rgBrep = rhBrep if isinstance(rhBrep, rg.Brep) else rs.coercebrep(rhBrep)
    for idx in reversed(idxs_Edges_in_border):
        fLength = rgBrep.Edges[idx].GetLength()
        if fLength < fMinSegLengthToKeep:
            idxs_Edges_in_border.remove(idx)


def get_joined_curves(rhBrep, idxs_Edges_per_border, fJoinTol=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBrep: rg.Brep, rd.BrepObject, or (GUID of brep)
        idxs_Edges_per_border: These are ordered to border loop.
        fJoinTol
        bEcho
        bDebug
        
    Returns on success:
        list(rg.Curve (0-level deep))
    """

    rgBrep = rs.coercebrep(rhBrep)

    rgCrvs_Borders = []

    for idxEs_of_border in idxs_Edges_per_border:
        rgEdges_in_border = [rgBrep.Edges[idx] for idx in idxEs_of_border]
        rgCrvs = rg.Curve.JoinCurves(rgEdges_in_border, joinTolerance=fJoinTol)
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


def processBrep(rgBrep, idxs_Trims, fMinSegLengthToKeep=None, fJoinTol=None, bEcho=True, bDebug=False):

    idxs_Edges_per_border = getEdgeIndicesOfBorders(
        rgBrep,
        idxs_Trims,
        bDebug=bDebug)
    if idxs_Edges_per_border is None: return

    if fMinSegLengthToKeep:
        for idxEs_of_border in idxs_Edges_per_border:
            removeIndicesOfShortEdges(
                rgBrep,
                idxEs_of_border,
                fMinSegLengthToKeep=fMinSegLengthToKeep)

    rgCrvs_Borders = get_joined_curves(
        rgBrep,
        idxs_Edges_per_border,
        fJoinTol=fJoinTol,
        bEcho=bEcho,
        bDebug=bDebug)
    return rgCrvs_Borders


def processBrepObjects(rhBreps, idx_Trims_PerBrep, fMinSegLengthToKeep=None, fJoinTol=None, bEcho=True, bDebug=False):
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

        rgCrvs_Borders = processBrep(
            rgBrep_In,
            idxs_Trims,
            fMinSegLengthToKeep=fMinSegLengthToKeep,
            fJoinTol=fJoinTol,
            bEcho=bEcho,
            bDebug=bDebug)
        if rgCrvs_Borders is None:
            continue

        gCurves_Borders = add_CurvesOfBorders(
            rgCrvs_Borders,
            bEcho,
            bDebug)

        if gCurves_Borders is None:
            continue

        gCurves_Borders_PerBrep.append(gCurves_Borders)
    
    return gCurves_Borders_PerBrep


def main():

    objrefs = getInput()
    if not objrefs: return

    fJoinTol = Opts.values['fJoinTol']
    fMinSegLengthToKeep = Opts.values['fMinSegLengthToKeep']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    #if not bDebug: sc.doc.Views.RedrawEnabled=False

    rc = convertObjrefsToSortedBrepsAndEdges(objrefs)
    if rc is None: return

    gBreps, idxs_Edges_per_Brep = rc

    sc.doc.Objects.UnselectAll()
    
    gCurves_Borders_PerBrep = processBrepObjects(
        gBreps,
        idxs_Edges_per_Brep,
        fMinSegLengthToKeep=fMinSegLengthToKeep,
        fJoinTol=fJoinTol,
        bEcho=bEcho,
        bDebug=bDebug)
    if gCurves_Borders_PerBrep is None: return
    
    sc.doc.Views.Redraw()
    
    iCt_Crvs_added = 0
    iCt_Open = 0
    for gCurves_borders in gCurves_Borders_PerBrep:
        for gCurve in gCurves_borders:
            if not rs.IsCurveClosed(gCurve):
                iCt_Open += 1
            iCt_Crvs_added += 1
    
    if iCt_Crvs_added == 1:
        if not iCt_Open:
            print("1 closed curve added to document.")
        else:
            print("1 open curve added to document.")
    else:
        if iCt_Open == 0:
            print("{} closed curves added to document.".format(iCt_Crvs_added))
        elif iCt_Open == iCt_Crvs_added:
            print("{} open curves added to document.".format(iCt_Open))
        else:
            print("{} open and {} closed curves added to document.".format(
                iCt_Open, iCt_Crvs_added-iCt_Open))


if __name__ == '__main__': main()