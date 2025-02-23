"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
181009-10: Created.
181026: Fixed bug in surface output.
190413: Added reference to replace local function.
200618: Added PerFaceColor support for V7.
250222: Added command options.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bDetachAll'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAllSimilar'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDupCrvsOnSrfBorders'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fMinSegLengthToKeep'; keys.append(key)
    values[key] = 2.0 * sc.doc.ModelAbsoluteTolerance
    #value = 1e-3 * Rhino.RhinoMath.UnitScale(
    #    Rhino.UnitSystem.Millimeters,
    #    sc.doc.ModelUnitSystem)
    #values[key] = float(format(value, '.0e'))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fJoinTol'; keys.append(key)
    value = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    values[key] = float(format(value, '.0e'))
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


def getInput():
    """
    Get faces or edges.
    """

    go = ri.Custom.GetObject()

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('bDetachAll')
        if not Opts.values['bDetachAll']:
            addOption('bAllSimilar')
        addOption('bDupCrvsOnSrfBorders')
        addOption('fMinSegLengthToKeep')
        addOption('fJoinTol')
        addOption('bEcho')
        addOption('bDebug')

        if Opts.values['bDetachAll']:
            go.SubObjectSelect = False
            go.SetCommandPrompt("Select face to untrim")
            go.GeometryFilter = rd.ObjectType.Surface
            go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.TrimmedSurface
        else:
            go.SubObjectSelect = True
            go.SetCommandPrompt("Select edge of loop to detach")
            go.GeometryFilter = rd.ObjectType.EdgeFilter

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


def brep_GUIDs_and_indices_of_trims_from_objrefs(objrefs):
    gBreps = []
    idx_Trims_per_brep = []

    for objref in objrefs:
        idx_rgTrim = objref.GeometryComponentIndex.Index
        gBrep = objref.ObjectId

        if not gBrep in gBreps:
            gBreps.append(gBrep)
            idx_Trims_per_brep.append([idx_rgTrim])
        else:
            idx_Trims_per_brep[gBreps.index(gBrep)].append(idx_rgTrim)

    return (
        gBreps,
        idx_Trims_per_brep,
        )


def printBrepElementCounts(rgBrep):
    
    print('-'*10)
    sEval = 'rgBrep.Surfaces.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Faces.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Curves2D.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Curves3D.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Loops.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Trims.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Edges.Count'; print(sEval+':', eval(sEval))
    sEval = 'rgBrep.Vertices.Count'; print(sEval+':', eval(sEval))
    print('-'*30)


def duplicateBrepPerLoops(rgBrep_In, idxLoops_ToDetach, bEcho=None, bDebug=None):
    
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    if bDebug: print("Make blank Brep.")
    
    rgB_Out = rg.Brep()
    
    #printBrepElementCounts(rgBrep1)
    
    # Outer loop.
    
    for idxLoop_ToDetach in idxLoops_ToDetach:
        rgF_In = rgBrep_In.Faces[0]
        if rgBrep_In.Loops[idxLoop_ToDetach].LoopType == rg.BrepLoopType.Outer:
            if bDebug: print("Add Face directly from Surface, including its outer loop.")
            rgF_B_Out = rgB_Out.Faces.Add(rgF_In.DuplicateSurface())
            break
    else:
        if bDebug: print("Add Surface.")
        idxSurface = rgB_Out.AddSurface(rgF_In.DuplicateSurface())
        if bDebug: print("Add Face.")
        rgF_B_Out = rgB_Out.Faces.Add(idxSurface)

    if Rhino.RhinoApp.ExeVersion >= 7:
        rgF_B_Out.PerFaceColor = rgF_In.PerFaceColor


    if bDebug: printBrepElementCounts(rgB_Out)
    
    if bDebug: print("Loop through loops, building Brep.")
    
    # Order of B0 lists will correspond those of B1.
    idx_Crvs2d_B0 = []
    idx_Crvs2d_B1 = []
    
    idx_Crvs3d_B0 = []
    idx_Crvs3d_B1 = []
    
    idxVertices_B0 = []
    idxVertices_B1 = []
    
    idx_Edges_B0 = []
    idx_Edges_B1 = []
    
    idx_Trims_B0 = []
    idx_Trims_B1 = []
    
    for iL, rgLoop_B0 in enumerate(rgBrep_In.Loops):
        if bDebug: print("Loop {}:".format(iL))
        if rgLoop_B0.LoopIndex in idxLoops_ToDetach:
            if bDebug: print("This loop will no be copied.  Continue to next loop.")
            continue
        
        if (
                rgLoop_B0.LoopType == rg.BrepLoopType.Outer
                and
                rgB_Out.Loops.Count > 0
        ):
            if bDebug: print("Outer loop is already in rgBrep1.  Continue to next loop.")
            continue
        
        if bDebug: print("Copy all related Brep elements of this loop to Brep1.")
        
        if bDebug: print("Add Loop.")
        rgLoop_B1 = rgB_Out.Loops.Add(rgLoop_B0.LoopType, rgB_Out.Faces[0])
        
        if bDebug: print("Loop through Trims in this Loop.")
        for rgTrim_B0 in rgLoop_B0.Trims:
            idxTrim_B0 = rgTrim_B0.TrimIndex
            idx_Trims_B0.append(idxTrim_B0)
            rgEdge_B0 = rgTrim_B0.Edge
            idxEdge_B0 = rgEdge_B0.EdgeIndex
            
            if idxEdge_B0 not in idx_Edges_B0:
                idx_Edges_B0.append(idxEdge_B0)
                rgEdge_B1 = idxEdge_B1 = None # to show that Edge needs to be added to Brep1.
            else:
                if bDebug: print("Using existing Edge.")
                idxEdge_B1 = idx_Edges_B1[idx_Edges_B0.index(idxEdge_B0)]
                rgEdge_B1 = rgB_Out.Edges[idxEdge_B1]
            
            if bDebug: print("Trim{}  Edge{}:".format(idxTrim_B0, idxEdge_B0))
            
            if bDebug: print("Duplicate and add Curve for Trim.")
            rgCrv2d_B0 = rgTrim_B0.TrimCurve
            idx_Crvs2d_B0.append(rgTrim_B0.TrimCurveIndex)
            rgCrv2d_B1 = rgCrv2d_B0.Duplicate()
            idxCrv2d_B1 = rgB_Out.Curves2D.Add(rgCrv2d_B1)
            idx_Crvs2d_B1.append(idxCrv2d_B1)
            
            
            if bDebug: print("Duplicate and add Curve for Edge.")
            idxCrv3d_B0 = rgEdge_B0.EdgeCurveIndex
            if idxCrv3d_B0 not in idx_Crvs3d_B0:
                idx_Crvs3d_B0.append(idxCrv3d_B0)
                rgCrv3d_B0 = rgEdge_B0.EdgeCurve
                rgCrv3d_B1 = rgCrv3d_B0.Duplicate()
                idxCrv3d_B1 = rgB_Out.Curves3D.Add(rgCrv3d_B1)
                idx_Crvs3d_B1.append(idxCrv3d_B1)
            else:
                if bDebug: print("Using existing Curve for Edge.")
            
            if bDebug: print("Process Start Vertex of Trim.")
            rgVertex_TrimStart_B0 = rgTrim_B0.StartVertex
            idxVertex_TrimStart_B0 = rgVertex_TrimStart_B0.VertexIndex
            if idxVertex_TrimStart_B0 not in idxVertices_B0:
                if bDebug: print("Add Vertex.")
                idxVertices_B0.append(idxVertex_TrimStart_B0)
                rgVertex_TrimStart_B1 = rgB_Out.Vertices.Add(rgVertex_TrimStart_B0.Location, vertexTolerance=sc.doc.ModelAbsoluteTolerance)
                idxVertex_TrimStart_B1 = rgVertex_TrimStart_B1.VertexIndex
                idxVertices_B1.append(idxVertex_TrimStart_B1)
            else:
                if bDebug: print("Using existing Vertex.")
                idxVertex_TrimStart_B1 = idxVertices_B1[idxVertices_B0.index(idxVertex_TrimStart_B0)]
            
            if bDebug: print("Process End Vertex of Trim.")
            rgVertex_TrimEnd_B0 = rgTrim_B0.EndVertex
            idxVertex_TrimEnd_B0 = rgVertex_TrimEnd_B0.VertexIndex
            if idxVertex_TrimEnd_B0 not in idxVertices_B0:
                if bDebug: print("Add Vertex.")
                idxVertices_B0.append(idxVertex_TrimEnd_B0)
                rgVertex_TrimEnd_B1 = rgB_Out.Vertices.Add(rgVertex_TrimEnd_B0.Location, vertexTolerance=sc.doc.ModelAbsoluteTolerance)
                idxVertex_TrimEnd_B1 = rgVertex_TrimEnd_B1.VertexIndex
                idxVertices_B1.append(idxVertex_TrimEnd_B1)
            else:
                if bDebug: print("Using existing Vertex.")
                idxVertex_TrimEnd_B1 = idxVertices_B1[idxVertices_B0.index(idxVertex_TrimEnd_B0)]
            
            if bDebug: print("Determine Start and End Vertices of Edge.")
            if rgTrim_B0.IsReversed():
                if bDebug: print("Trim IS reversed to Edge.")
                idxVertex_EdgeStart_B1 = idxVertex_TrimEnd_B1
                idxVertex_EdgeEnd_B1 = idxVertex_TrimStart_B1
                rev3d = True
            else:
                if bDebug: print("Trim is NOT reversed to Edge.")
                idxVertex_EdgeStart_B1 = idxVertex_TrimStart_B1
                idxVertex_EdgeEnd_B1 = idxVertex_TrimEnd_B1
                rev3d = False
            
            if idxEdge_B1 is None:
                if bDebug: print("Add Edge.")
                rgEdge_B1 = rgB_Out.Edges.Add(
                        rgB_Out.Vertices[idxVertex_EdgeStart_B1],
                        rgB_Out.Vertices[idxVertex_EdgeEnd_B1],
                        idxCrv3d_B1,
                        rgEdge_B0.Tolerance)
                idxEdge_B1 = rgEdge_B1.EdgeIndex
                idx_Edges_B1.append(idxEdge_B1)
            else:
                if bDebug: print("Using existing edge for Trim.")
            
            if bDebug: print("Add Trim.")
            rgTrim_B1 = rgB_Out.Trims.Add(
                    rgEdge_B1,
                    rev3d=rev3d,
                    loop=rgLoop_B1,
                    curve2dIndex=idxCrv2d_B1)
            idxTrim_B1 = rgTrim_B1.TrimIndex
            idx_Trims_B1.append(idxTrim_B1)
            
            if bDebug: print("Set Trim tolerance.")
            tolU, tolV = rgTrim_B0.GetTolerances()
            #            if bDebug:
            #                sEval = 'tolU'; print(sEval+':', eval(sEval))
            #                sEval = 'tolV'; print(sEval+':', eval(sEval))
            rgTrim_B1.SetTolerances(tolU, tolV)
            #            if bDebug:
            #                tolU, tolV = rgTrim_B1.GetTolerances()
            #                sEval = 'tolU'; print(sEval+':', eval(sEval))
            #                sEval = 'tolV'; print(sEval+':', eval(sEval))
        
            if bDebug: print("Set Trim IsoStatus.")
            rgTrim_B1.IsoStatus = rgTrim_B0.IsoStatus
            if bDebug: sEval = 'rgTrim_B1.IsoStatus'; print(sEval+':', eval(sEval))
    
    if bDebug: printBrepElementCounts(rgB_Out)
    
    return rgB_Out


def isSenw(rgTrim):
    return (rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.South or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.East or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.North or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.West)


def processBrep(gB_In, idx_rgTrims=None, bDetachAll=None, bAllSimilar=None, bDupCrvsOnSrfBorders=None, fMinSegLengthToKeep=None, fJoinTol=None, bEcho=True, bDebug=False):

    if bDetachAll is None: bDetachAll = Opts.values['bDetachAll']
    if bAllSimilar is None: bAllSimilar = Opts.values['bAllSimilar']
    if bDupCrvsOnSrfBorders is None: bDupCrvsOnSrfBorders = Opts.values['bDupCrvsOnSrfBorders']
    if fMinSegLengthToKeep is None: bDebug = Opts.values['fMinSegLengthToKeep']
    if fJoinTol is None: bDebug = Opts.values['fJoinTol']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']

    def outputUnderlyingSurface():
        """
        """
        

        rgF_In = rgBrep_In.Faces[0]

        srf = rgF_In.UnderlyingSurface()
        if srf is None:
            if bEcho: print("Underlying surface cannot be obtained.")
            return False
        
        if all(bAddCurve_successes):
            rgB_Out = srf.ToBrep()
            if Rhino.RhinoApp.ExeVersion >= 7:
                rgB_Out.Faces[0].PerFaceColor = rgF_In.PerFaceColor
            bAdded = sc.doc.Objects.Replace(gB_In, rgB_Out)
            if not bAdded:
                if bEcho: print("Brep could not be replaced with underlying surface.")
            return bAdded
        else:
            # Not all curves could be added.
            rgB_Out = srf.ToBrep()
            gB_Out = sc.doc.Objects.AddBrep(rgB_Out)
            if gB_Out == Guid.Empty:
                if bEcho: print("Underlying surface could not be added to document.")
                return False
            if bEcho:
                print("Underlying surface was added instead of untrimming the face was not untrimmed because not all trim edge curves could be added.")
            return gB_Out

    def addCurves():
        """
        """
        
        def edgesOfLoop():
            """
            """
            
            rgEdges_ThisLoop = []
            for rgTrim in rgLoop.Trims:
                if rgTrim.TrimType == rg.BrepTrimType.Seam:
                    continue
                elif not bDupCrvsOnSrfBorders and isSenw(rgTrim):
                    continue
                edge = rgTrim.Edge
                if edge is None:
                    raise Exception("Edge is None.")
                    #bAddCurve_successes.append(False)
                    #continue
                if fMinSegLengthToKeep:
                    fLength = edge.GetLength()
                    if fLength < fMinSegLengthToKeep:
                        continue
                bAddCurve_successes.append(True)
                rgEdges_ThisLoop.append(edge)
            return rgEdges_ThisLoop
        
        rgEdges_ThisLoop = edgesOfLoop()
        if rgEdges_ThisLoop is None:
            return

        i_Crv_cts[0] += len(rgEdges_ThisLoop)

        rgCrvs_Joined = rg.Curve.JoinCurves(
            rgEdges_ThisLoop,
            joinTolerance=fJoinTol)

        i_Crv_cts[1] += len(rgCrvs_Joined)

        for rgCrv in rgCrvs_Joined:
            gCrv = sc.doc.Objects.AddCurve(rgCrv)
            if gCrv == Guid.Empty:
                bAddCurve_successes.append(False)
            else:
                bAddCurve_successes.append(True)


    rgBrep_In = rs.coercebrep(gB_In)
    
    if rgBrep_In.IsSurface:
        if bEcho: print("Monoface brep is not trimmed.  Its borders are the same as those of its underlying surface.")
        return
    
    bAddCurve_successes = []
    
    i_Crv_cts = [0, 0] # Before join, After join
    
    if not bDetachAll:
        # Check for bAllSimilar with trims from both Outer and Inner loops selected.
        if bAllSimilar:
            bOuter = any(rgTrim.TrimIndex in idx_rgTrims for rgTrim in rgBrep_In.Faces[0].OuterLoop.Trims)
            if bDebug: sEval = 'bOuter'; print(sEval+':', eval(sEval))
            bInner = any(rgTrim.TrimIndex in idx_rgTrims for rgLoop in rgBrep_In.Loops if rgLoop.LoopType == rg.BrepLoopType.Inner for rgTrim in rgLoop.Trims)
            if bDebug: sEval = 'bInner'; print(sEval+':', eval(sEval))
            if bOuter and bInner:
                bDetachAll = True
    
    if bDetachAll:
        for rgLoop in rgBrep_In.Loops:
            addCurves()
        
        print("{} joined curves from {} curves added to document.".format(
                i_Crv_cts[1], i_Crv_cts[0]))
        
        if i_Crv_cts[1]:
            rgB_Res = outputUnderlyingSurface()
    else:
        # Not bDetachAll.
        
        idxLoops_ToDetach = []
        
        if bAllSimilar:
            if bOuter:
                rgLoop = rgBrep_In.Faces[0].OuterLoop
                
                addCurves()
                
                if all(bAddCurve_successes):
                    idxLoops_ToDetach.append(rgLoop.LoopIndex)
            
            if bInner:
                for rgLoop in rgBrep_In.Loops:
                    
                    if rgLoop.LoopType != rg.BrepLoopType.Inner:
                        continue
                    
                    addCurves()
                    
                    if all(bAddCurve_successes):
                        idxLoops_ToDetach.append(rgLoop.LoopIndex)
            
        else:
            # Not bAllSimilar.
            for rgLoop in rgBrep_In.Loops:
                
                for rgTrim in rgLoop.Trims:
                    if rgTrim.TrimIndex in idx_rgTrims:
                        # Detach this loop.
                        break
                else:
                    # No given trims found in this loop.
                    continue # to next loop.
                
                addCurves()
                
                if all(bAddCurve_successes):
                    idxLoops_ToDetach.append(rgLoop.LoopIndex)
        
        print("{} joined curves from {} curves added to document.".format(
            i_Crv_cts[1], i_Crv_cts[0]))
        
        if not i_Crv_cts[1]:
            if bEcho: print("Brep was not modified.")
            return
        
        # If all loops were selected, then output the entire underlying surface.
        if len(idxLoops_ToDetach) == rgBrep_In.Loops.Count:
            # Detach all trim.
            return outputUnderlyingSurface()
        
        if bDebug: print("Create brep with loops removed.")
        
        rgB_Res = duplicateBrepPerLoops(rgBrep_In, idxLoops_ToDetach)
        
        if rgB_Res is None:
            sEval = 'rgBrep1'; print(sEval+':', eval(sEval))
            return False
        if not rgB_Res.IsValid:
            sEval = 'rgBrep1.IsValid'; print(sEval+':', eval(sEval))
            b, sLog = rgB_Res.IsValidWithLog()
            if not b: print(sLog)
            return False
        
        if bDebug: sEval = 'rgBrep1.Loops.Count'; print(sEval+':', eval(sEval))
        
        if all(bAddCurve_successes):
            bReplaced = sc.doc.Objects.Replace(gB_In, rgB_Res)
        else:
            gBrep1 = sc.doc.Objects.AddBrep(rgB_Res)
            if gBrep1 == Guid.Empty:
                print("Modified brep could not be added to document.")
                return False
            print("Brep was not replaced because all trim edge curves could not be added.")
            return False
    
    return rgB_Res


def processBreps(gBreps0, idx_Trims_per_brep, bDetachAll=None, bAllSimilar=None, bDupCrvsOnSrfBorders=None, fMinSegLengthToKeep=None, fJoinTol=None, bEcho=True, bDebug=False):

    i_Success_ct = 0

    for gB_In, idx_rgTrims in zip(gBreps0, idx_Trims_per_brep):
        if processBrep(
            gB_In,
            idx_rgTrims,
            bDetachAll=bDetachAll,
            bAllSimilar=bAllSimilar,
            bDupCrvsOnSrfBorders=bDupCrvsOnSrfBorders,
            fMinSegLengthToKeep=fMinSegLengthToKeep,
            fJoinTol=fJoinTol,
            bEcho=bEcho,
            bDebug=bDebug
        ):
            i_Success_ct += 1

    print("{} out of {} monoface breps have been untrimmed.".format(
            i_Success_ct, len(gBreps0)))


def main():
    
    objrefs = getInput()
    if objrefs is None: return

    rv = brep_GUIDs_and_indices_of_trims_from_objrefs(objrefs)
    if not rv: return

    (
        gBreps,
        idx_Trims_per_brep,
        ) = rv

    bDetachAll = Opts.values['bDetachAll']
    bAllSimilar = Opts.values['bAllSimilar']
    bDupCrvsOnSrfBorders = Opts.values['bDupCrvsOnSrfBorders']
    fMinSegLengthToKeep = Opts.values['fMinSegLengthToKeep']
    fJoinTol = Opts.values['fJoinTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    processBreps(
        gBreps,
        idx_Trims_per_brep,
        bDetachAll=bDetachAll,
        bAllSimilar=bAllSimilar,
        bDupCrvsOnSrfBorders=bDupCrvsOnSrfBorders,
        fMinSegLengthToKeep=fMinSegLengthToKeep,
        fJoinTol=fJoinTol,
        bEcho=bEcho,
        bDebug=bDebug
        )

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()