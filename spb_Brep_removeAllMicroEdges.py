"""
Alternative to _RemoveAllNakedMicroEdges, this script will remove all micro-edges,
both naked and interior, with no limiation of that the edges
"fold or loop back on themselves and have no matching edge to which they can be joined"
as stated in https://docs.mcneel.com/rhino/7/help/en-us/index.htm#commands/removeallnakedmicroedges.htm

Both naked and interior micro-BrepEdges, as defined by the MaxEdgeLengthToRemove
option, will be removed along with the relative BrepTrims.
For each edge removed, both end vertices are replaced with a single vertex at
the average location.

If this script fails, run _RemoveAllNakedMicroEdges then this script.
If that also fails, _MergeAllEdges then this script.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums ( https://discourse.mcneel.com/ ).
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221114-18: Created.

WIP:
    Having trouble matching the ends of some Loop-consecutive Trims using MatchEnds method
    when matches occur loop by loop as the trims are added.  This may just be the case
    for mated edges and MatchEnds for the entire Brep is the best solution.
    
    For problems occuring at consecutive micro-edges,
    First, merge all consecutive micro-edges in each loop?
    Maybe the consecutive micro-edge Curve2Ds should be replaced with a single, cubic Bezier.

TODO:
    Merge edges first to handle some consecutive micro-edges.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fMaxMicroLength'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    names[key] = 'MaxEdgeLengthToRemove'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bDupEdgesRemoved'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

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

        if key == 'fMaxMicroLength':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal breps when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fMaxMicroLength')
        addOption('bDupEdgesRemoved')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            rdBrepObjs = [o.Object() for o in go.Objects()]
            go.Dispose()
            return rdBrepObjs

        if res == ri.GetResult.Nothing:
            oes = rd.ObjectEnumeratorSettings()
            oes.NormalObjects = True
            oes.LockedObjects = False
            oes.IncludeLights = False
            oes.IncludeGrips = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep
            rdBrepObjs = list(sc.doc.Objects.GetObjectList(oes))
            go.Dispose()
            if len(rdBrepObjs) == 0: return
            return rdBrepObjs

        if res == ri.GetResult.Number:
            key = 'fMaxMicroLength'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getTrimEndDeviation(rgLoop, idxLT):
    return rgLoop.Trims[idxLT].PointAtEnd.DistanceTo(
        rgLoop.Trims[(idxLT + 1) % rgLoop.Trims.Count].PointAtStart)


def _maxTrimEndDeviation(rgLoop, bDebug):
    fMax = 0.0
    iT_A_Max = None
    iT_B_Max = None
    for iT_A in range(rgLoop.Trims.Count):
        dev = _getTrimEndDeviation(rgLoop, iT_A)
        if dev > fMax:
            fMax = dev
            iT_A_Max = iT_A
            iT_B_Max = iT_A + 1
    if bDebug and (fMax > 0.0): print(iT_A_Max, iT_B_Max)
    return fMax


def createBrep_RemoveEdges(rgBrep_In, idxEdges_In, bDebug=False):
    """
    Parameters:
        rgBrep_In
        idxEdges_In: list(int)
        bDebug: bool

    Returns on success:
        Brep
        None

    Returns on fail:
        None
        str("Explanation of fail or that brep doesn't contain micro-edges.")
    """


    rgB_In = rgBrep_In


    idxEs_In_ToRemove = idxEdges_In

    rgB_Out = rg.Brep()


    idxVs_Out_Added = [] # 1:1 relationship with idxEs_In_ToRemove, so they share a list index.

    idxTs_In_ToRemove = []

    idxVs_In_ToRemove = []
    idxVs_Out_Added_Per_idxVs_In_ToRemove = []

    for idxE_In in idxEs_In_ToRemove:
        rgE_In = rgB_In.Edges[idxE_In]

        idxTs_In_ToRemove.extend(rgE_In.TrimIndices())

        pt = 0.5*(rgE_In.StartVertex.Location + rgE_In.EndVertex.Location)
        rgV = rgB_Out.Vertices.Add(
            point=pt,
            vertexTolerance=sc.doc.ModelAbsoluteTolerance) # SetTolerancesBoxesAndFlags below should correct this.
        idxVs_Out_Added.append(rgV.VertexIndex)

        idxV_In = rgE_In.StartVertex.VertexIndex
        if idxV_In not in idxVs_In_ToRemove:
            idxVs_In_ToRemove.append(idxV_In)
            idxVs_Out_Added_Per_idxVs_In_ToRemove.append(rgV.VertexIndex)
        idxV_In = rgE_In.EndVertex.VertexIndex
        if idxV_In not in idxVs_In_ToRemove:
            idxVs_In_ToRemove.append(idxV_In)
            idxVs_Out_Added_Per_idxVs_In_ToRemove.append(rgV.VertexIndex)


    if bDebug:
        sEval = "idxEs_In_ToRemove"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "idxTs_In_ToRemove"; print("{}: {}".format(sEval, eval(sEval)))


    # Add surfaces.
    for srf_In in rgB_In.Surfaces:
        srf_Out = srf_In.Duplicate()
        rgB_Out.AddSurface(srf_Out)


    # Add edge curves (Curve3D).
    idxC3s_In = []
    idxC3s_Out_Per_idxC3s_In = []
    for idxE_In, rgE_In in enumerate(rgB_In.Edges):
        if rgE_In.EdgeIndex in idxEs_In_ToRemove:
            continue
        rgC3_Out = rgE_In.EdgeCurve.Duplicate()
        idxC3_Out = rgB_Out.AddEdgeCurve(rgC3_Out)
        idxC3s_Out_Per_idxC3s_In.append(idxC3_Out)
        idxC3s_In.append(rgE_In.EdgeCurveIndex)


    # Add trim curves (Curve2D).
    idxC2s_In = []
    idxC2s_Out_Per_idxC2s_In = []
    for idxT_In, rgT_In in enumerate(rgB_In.Trims):
        if rgT_In.TrimIndex in idxTs_In_ToRemove:
            continue
        rgC2_Out = rgT_In.TrimCurve.Duplicate()
        idxC2_Out = rgB_Out.AddTrimCurve(rgC2_Out)
        idxC2s_Out_Per_idxC2s_In.append(idxC2_Out)
        idxC2s_In.append(rgT_In.TrimCurveIndex)


    # Add faces.
    for iF, face_In in enumerate(rgB_In.Faces):
        rgB_Out.Faces.Add(face_In.SurfaceIndex)
        rgB_Out.Faces[iF].OrientationIsReversed = face_In.OrientationIsReversed
        rgB_Out.Faces[iF].PerFaceColor = face_In.PerFaceColor


    # Add kept vertices.  New ones were already added.
    idxVs_Out_Kept = []
    idxVs_In_Per_idxVs_Out_Kept = []
    for idxV_In, vertex_In in enumerate(rgB_In.Vertices):
        if idxV_In in idxVs_In_ToRemove:
            continue

        v_Out = rgB_Out.Vertices.Add(
            point=vertex_In.Location,
            vertexTolerance=sc.doc.ModelAbsoluteTolerance)  # SetTolerancesBoxesAndFlags below should correct this.
        idxVs_Out_Kept.append(v_Out.VertexIndex)
        idxVs_In_Per_idxVs_Out_Kept.append(idxV_In)


    # Add loops.
    for iL, loop_In in enumerate(rgB_In.Loops):
        rgB_Out.Loops.Add(
            loopType=loop_In.LoopType,
            face=rgB_Out.Faces[loop_In.Face.FaceIndex])


    # Add edges.
    idxTs_In_ToKeep = []
    idxEs_Out_Per_idxTs_In_ToKeep = []

    for idxE_In, rgE_In in enumerate(rgB_In.Edges):
        if idxE_In in idxEs_In_ToRemove:
            continue

        idxV_In = rgE_In.StartVertex.VertexIndex
        if idxV_In in idxVs_In_ToRemove:
            idxV_S = idxVs_Out_Added_Per_idxVs_In_ToRemove[idxVs_In_ToRemove.index(idxV_In)]
        else:
            idxV_S = idxVs_Out_Kept[idxVs_In_Per_idxVs_Out_Kept.index(idxV_In)]

        idxV_In = rgE_In.EndVertex.VertexIndex
        if idxV_In in idxVs_In_ToRemove:
            idxV_E = idxVs_Out_Added_Per_idxVs_In_ToRemove[idxVs_In_ToRemove.index(idxV_In)]
        else:
            idxV_E = idxVs_Out_Kept[idxVs_In_Per_idxVs_Out_Kept.index(idxV_In)]

        idxC3_Out = idxC3s_Out_Per_idxC3s_In[idxC3s_In.index(rgE_In.EdgeCurveIndex)]
        rgC3_Out = rgB_Out.Curves3D[idxC3_Out]

        b_rgE_and_rgC3_domains_match = (
            rgE_In.Domain.T0 == rgC3_Out.Domain.T0 and
            rgE_In.Domain.T1 == rgC3_Out.Domain.T1)

        if idxV_S == idxV_E:
            if not rgC3_Out.IsClosed:
                rgC3_Out.MakeClosed(sc.doc.ModelAbsoluteTolerance)
            b_rgE_and_rgC3_domains_match = True

        if b_rgE_and_rgC3_domains_match:
            rgE_Out = rgB_Out.Edges.Add(
                startVertexIndex=idxV_S,
                endVertexIndex=idxV_E,
                curve3dIndex=idxC3_Out,
                edgeTolerance=sc.doc.ModelAbsoluteTolerance) # SetTolerancesBoxesAndFlags below should correct this.
        else:
            rgE_Out = rgB_Out.Edges.Add(
                startVertexIndex=idxV_S,
                endVertexIndex=idxV_E,
                curve3dIndex=idxC3_Out,
                subDomain=rgE_In.Domain,
                edgeTolerance=sc.doc.ModelAbsoluteTolerance) # SetTolerancesBoxesAndFlags below should correct this.

        if idxV_S == idxV_E:
            if not rgE_Out.IsClosed:
                rgE_Out.MakeClosed(sc.doc.ModelAbsoluteTolerance)


        #if idxV_S == idxV_E:
        #    print(rgE_Out.IsClosed, rgE_Out.EdgeIndex)

        #if rgE_Out.EdgeIndex == 741:
        #    sc.doc.Objects.AddCurve(rgE_Out); sc.doc.Views.Redraw()
        #    print(rgE_Out.IsClosed, rgE_Out.IsClosable(sc.doc.ModelAbsoluteTolerance), rgE_Out.EdgeIndex, idxV_S, idxV_E)

        idxE_Out = rgE_Out.EdgeIndex
        for idxT_In in rgE_In.TrimIndices():
            idxTs_In_ToKeep.append(idxT_In)
            idxEs_Out_Per_idxTs_In_ToKeep.append(idxE_Out)


    # Add trims NOT per loop.
    for idxT_In, rgT_In in enumerate(rgB_In.Trims):
        if idxT_In in idxTs_In_ToRemove:
            continue

        idxC2_Out = idxC2s_Out_Per_idxC2s_In[idxC2s_In.index(rgT_In.TrimCurveIndex)]

        rgL_Out = rgB_Out.Loops[rgT_In.Loop.LoopIndex]

        if rgT_In.TrimType == rg.BrepTrimType.Singular:

            idxV_In = rgT_In.StartVertex.VertexIndex

            if idxV_In in idxVs_In_ToRemove:
                idxV_Out = idxVs_Out_Added_Per_idxVs_In_ToRemove[idxVs_In_ToRemove.index(idxV_In)]
            else:
                idxV_Out = idxVs_Out_Kept[idxVs_In_Per_idxVs_Out_Kept.index(idxV_In)]

            rgT_Out = rgB_Out.Trims.AddSingularTrim(
                vertex=rgB_Out.Vertices[idxV_Out],
                loop=rgL_Out,
                iso=rgT_In.IsoStatus,
                curve2dIndex=idxC2_Out)
        else:
            idxE_Out = idxEs_Out_Per_idxTs_In_ToKeep[idxTs_In_ToKeep.index(idxT_In)]

            rgT_Out = rgB_Out.Trims.Add(
                edge=rgB_Out.Edges[idxE_Out],
                rev3d=rgT_In.IsReversed(),
                loop=rgL_Out,
                curve2dIndex=idxC2_Out)

            rgT_Out.IsoStatus = rgT_In.IsoStatus

        rgT_Out.SetTolerances(
            toleranceU=Rhino.RhinoMath.ZeroTolerance,
            toleranceV=Rhino.RhinoMath.ZeroTolerance) # SetTolerancesBoxesAndFlags below should correct this.




    # Add trims PER loop.
    #for idxL_In, rgL_In in enumerate(rgB_In.Loops):
    #    rgL_Out = rgB_Out.Loops[rgT_In.Loop.LoopIndex]
    #    for rgT_In in rgL_In.Trims:
    #        idxT_In = rgT_In.TrimIndex
    #    #for idxT_In, rgT_In in enumerate(rgB_In.Trims):
    #        if idxT_In in idxTs_In_ToRemove:
    #            continue

    #        idxC2_Out = idxC2s_Out_Per_idxC2s_In[idxC2s_In.index(rgT_In.TrimCurveIndex)]

    #        if rgT_In.TrimType == rg.BrepTrimType.Singular:

    #            idxV_In = rgT_In.StartVertex.VertexIndex

    #            if idxV_In in idxVs_In_ToRemove:
    #                idxV_Out = idxVs_Out_Added_Per_idxVs_In_ToRemove[idxVs_In_ToRemove.index(idxV_In)]
    #            else:
    #                idxV_Out = idxVs_Out_Kept[idxVs_In_Per_idxVs_Out_Kept.index(idxV_In)]

    #            rgT_Out = rgB_Out.Trims.AddSingularTrim(
    #                vertex=rgB_Out.Vertices[idxV_Out],
    #                loop=rgL_Out,
    #                iso=rgT_In.IsoStatus,
    #                curve2dIndex=idxC2_Out)
    #        else:
    #            idxE_Out = idxEs_Out_Per_idxTs_In_ToKeep[idxTs_In_ToKeep.index(idxT_In)]

    #            rgT_Out = rgB_Out.Trims.Add(
    #                edge=rgB_Out.Edges[idxE_Out],
    #                rev3d=rgT_In.IsReversed(),
    #                loop=rgL_Out,
    #                curve2dIndex=idxC2_Out)

    #            rgT_Out.IsoStatus = rgT_In.IsoStatus

    #        rgT_Out.SetTolerances(
    #            toleranceU=Rhino.RhinoMath.ZeroTolerance,
    #            toleranceV=Rhino.RhinoMath.ZeroTolerance) # SetTolerancesBoxesAndFlags below should correct this.

            #if rgT_Out.TrimIndex in (715,716):
            #    sc.doc.Objects.AddCurve(rgT_Out.Edge); sc.doc.Views.Redraw()

        #fMaxDev = _maxTrimEndDeviation(rgL_Out, bDebug)
        #if fMaxDev == 0.0:
        #    continue

        #if bDebug: print("L[{}], Dev:{}".format(idxL_In, fMaxDev))

        #bMatchedTrimEnds = rgB_Out.Trims.MatchEnds(rgL_Out)
        ##if bDebug:
        #if not bMatchedTrimEnds:
        #    for idxLT, rgT in enumerate(rgL_Out.Trims):
        #        if _getTrimEndDeviation(rgL_Out, idxLT) == 0.0:
        #            continue
        #        bMatchedTrimEnds_2 = rgB_Out.Trims.MatchEnds(rgT, rgL_Out.Trims[(idxLT + 1) % rgL_Out.Trims.Count])
        #        if bMatchedTrimEnds_2:
        #            continue

        #        if bDebug:
        #            for idxLT, rgT in enumerate(rgL_Out.Trims):
        #                sc.doc.Objects.AddCurve(rgT)
        #                rgDot = rg.TextDot(
        #                    "idxLT[{}], idxT[{}]".format(idxLT, rgT.TrimIndex),
        #                    rgT.PointAt(rgT.Domain.Mid))
        #                sc.doc.Objects.AddTextDot(rgDot)

        #            #sc.doc.Objects.AddBrep(rgL_In.Face.DuplicateFace(duplicateMeshes=True))

        #        rgB_Out.Dispose()
        #        #1/0
        #        return None, "Failed to match BrepTrim ends in Loop[{}].".format(idxL_In)

        #print("Matched BrepTrim ends in Loop[{}]: {}".format(idxL_In, bMatchedTrimEnds))

        #fMaxDev = _maxTrimEndDeviation(rgL_Out)
        #if fMaxDev > 0.0:
        #print(idxL_In, fMaxDev)




    bMatchedTrimEnds_All = rgB_Out.Trims.MatchEnds()

    if bDebug:
        sEval = "rgB_In.Loops.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "rgB_Out.Loops.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "rgB_Out.Edges.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "rgB_Out.Vertices.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "bMatchedTrimEnds_All"; print("{}: {}".format(sEval, eval(sEval)))


    bIsValid, sLog = rgB_Out.IsValidWithLog()

    if bIsValid:
        if bDebug: print("rgB_Out is initially valid.")

        rgB_Out.SetTolerancesBoxesAndFlags(
            bLazy=False,
            bSetVertexTolerances=True,
            bSetEdgeTolerances=True,
            bSetTrimTolerances=True,
            bSetTrimIsoFlags=True,
            bSetTrimTypeFlags=True,
            bSetLoopTypeFlags=True,
            bSetTrimBoxes=True)

        return rgB_Out, None


    if bDebug: print("Initial brep:", sLog)


    rgB_Out.Repair(sc.doc.ModelAbsoluteTolerance)

    bIsValid, sLog = rgB_Out.IsValidWithLog()

    if not bIsValid:
        print(sLog)
        #sc.doc.Objects.AddCurve(rgB_Out.Trims[715].Edge)
        #sc.doc.Objects.AddCurve(rgB_Out.Trims[716].Edge)
        #sc.doc.Views.Redraw(); 1/0
        rgB_Out.Dispose()
        return None, "After Brep.Repair brep: {}".format(sLog)

    if bDebug: print("Brep.Repair was successful.")

    rgB_Out.SetTolerancesBoxesAndFlags(
        bLazy=False,
        bSetVertexTolerances=True,
        bSetEdgeTolerances=True,
        bSetTrimTolerances=True,
        bSetTrimIsoFlags=True,
        bSetTrimTypeFlags=True,
        bSetLoopTypeFlags=True,
        bSetTrimBoxes=True)


    return rgB_Out, None


def createBrep_RemoveMicroEdges(rgBrep_In, fMaxMicroLength, bDupEdgesRemoved=False, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep_In
        fMaxMicroLength: float
        bDebug: bool

    Returns on success:
        Brep
        None

    Returns on fail:
        None
        str("Explanation of fail or that brep doesn't contain micro-edges.")
    """


    rgB_In = rgBrep_In


    def getMicroEdges(rgB):
        idxEs_Micro = []
        for rgE_In in rgB.Edges:
            if rgE_In.GetLength() <= fMaxMicroLength:
                idxEs_Micro.append(rgE_In.EdgeIndex)
        return idxEs_Micro


    idxEs_Micro_In = idxEs_Micro = getMicroEdges(rgB_In)


    if not idxEs_Micro:
        return None, "Brep does not have any micro-edges."

    #print("{} micro-edges found.".format(len(idxEs_Micro)))


    # For debugging
    #if bDebug:
    #    rgB_Res, sLog = createBrep_RemoveEdges(
    #        rgB_In,
    #        [idxEs_Micro[0]],
    #        bDebug=True)
    #    if rgB_Res:
    #        idxEs_Micro_End = getMicroEdges(rgB_Res)

    #        sLog = "Micro edge count change: {} -> {}".format(
    #            len(idxEs_Micro_In), len(idxEs_Micro_End))
    #        #if bEcho: print(sLog)
    #        return rgB_Res, sLog

    #    return None, sLog



    # Try to create brep with all micro-edges removed.

    rgB_Res, sLog = createBrep_RemoveEdges(
        rgB_In,
        idxEs_Micro,
        bDebug=bDebug)

    if rgB_Res:
        if bDupEdgesRemoved:
            for idxE in idxEs_Micro_In:
                sc.doc.Objects.AddCurve(rgB_In.Edges[idxE])

        idxEs_Micro_End = getMicroEdges(rgB_Res)

        sLog = "Micro edge count change: {} -> {}".format(
            len(idxEs_Micro_In), len(idxEs_Micro_End))
        #if bEcho: print(sLog)
        return rgB_Res, sLog


    if bEcho:
        print("Removing all micro-edges at once failed because {}".format(sLog))

    if len(idxEs_Micro_In) == 1:
        return None, sLog


    print("Now, trying to remove one at a time ...")

    rgB_LastGood = rgB_In.DuplicateBrep()

    i = 0

    while i < len(idxEs_Micro):
        sc.escape_test()

        rgB_Res, sLog = createBrep_RemoveEdges(
            rgB_LastGood,
            [idxEs_Micro[i]],
            bDebug=bDebug)

        if rgB_Res:
            rgB_LastGood.Dispose()
            rgB_LastGood = rgB_Res
            idxEs_Micro = getMicroEdges(rgB_LastGood)
            print("Removed edge {} of {} remaining.".format(
                i+1, len(idxEs_Micro)))
            if bDupEdgesRemoved:
                sc.doc.Objects.AddCurve(rgB_LastGood.Edges[idxEs_Micro[i]])
            i = 0
        else:
            print("Failed to remove edge {} of {} remaining because {}".format(
                i+1, len(idxEs_Micro), sLog))
            i += 1

    idxEs_Micro_End = getMicroEdges(rgB_LastGood)

    if len(idxEs_Micro_End) < len(idxEs_Micro_In):
        sLog = "Micro edge count change: {} -> {}".format(
            len(idxEs_Micro_In), len(idxEs_Micro_End))
        #if bEcho: print(sLog)
        return rgB_LastGood, sLog

    rgB_LastGood.Dispose()

    return None, sLog


def processBrepObjects(rdBs_In, fMaxMicroLength, bDupEdgesRemoved=False, bEcho=True, bDebug=False):
    """
    Parameters:
        rdBs_In: list(DocObjects.BrepObject)
    """

    gBs_Modified = []

    sLogs = []

    for iB, rdB in enumerate(rdBs_In):

        Rhino.RhinoApp.SetCommandPrompt(
                prompt="Processing brep {} of {} ...".format(iB+1, len(rdBs_In)))

        rgB_Res, sLog = createBrep_RemoveMicroEdges(
            rdB.Geometry,
            fMaxMicroLength=fMaxMicroLength,
            bDupEdgesRemoved=bDupEdgesRemoved,
            bEcho=bEcho,
            bDebug=bDebug,
            )

        if sLog is not None:
            sLogs.append(sLog)

        if not rgB_Res:
            continue

        if sc.doc.Objects.Replace(rdB.Id, brep=rgB_Res):
            gBs_Modified.append(rdB.Id)
            sLogs.append("Brep was modified.".format(rdB.Id))

    if len(rdBs_In) == 1:
        print(sLogs[0])
    else:
        for sLog in sorted(set(sLogs)):
            print("[{}] {}".format(sLogs.count(sLog), sLog))

    if not gBs_Modified:
        print("No breps were modified.")


    return gBs_Modified


def main():

    rdBs_In = getInput()
    if rdBs_In is None: return

    fMaxMicroLength = Opts.values['fMaxMicroLength']
    bDupEdgesRemoved = Opts.values['bDupEdgesRemoved']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    processBrepObjects(
        rdBs_In,
        fMaxMicroLength=fMaxMicroLength,
        bDupEdgesRemoved=bDupEdgesRemoved,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()