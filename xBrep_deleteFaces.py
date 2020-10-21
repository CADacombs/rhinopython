"""
160531-0601: Created.
160824: Now selection of faces can continue after chain fillet selection.
160829: Adjacent faces are now trimmed to one another and attempted to be joined back into a solid.
        Added bAddCrvs option.
160831: Import-related update.
        joinBrepsWithOneIndividually now allows main breps to start with no faces.
160907-08: Modularizations and now adjacent faces will extend if the solid doesn't heal using other techniques.
160911: Removed brep layer of nesting in orderedLists_AdjacentFLT_perSelFBorders.
        Moved orderedLists_AdjacentFLT_perSelFBorders from brep.py.
160914-15: Now curvesForTrim uses Curve.ClosestPoints before attempting to pull naked edges of polyface breps to face to split.
180612, 190423, 200107: Import-related update.
200107-10: Restarted development of this script.
           Added bHeal and bAddModFaces.
           Now tries Brep.JoinBreps instead of Brep.Join due to forced edge joining of the latter.
           Now determines whether resultant brep is final based on naked edge border counts before versus after processing.
200619: Import-related update.

TODO: 
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System.Diagnostics import Stopwatch

import xBrep_findMatchingFace
import xBrep_nakedEdgeLoop
import xBrepFace
import xBrepFace_extendAtUntrimmedEdge
import xBrepFace_removeTrim
import xBrepObject
import xInput


stopwatch = Stopwatch() # One instance will be used for all tests.
timeLimit = 0.1 # Occurrences of some methods requiring more time (seconds) than this value will be stated.


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'bHeal'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddModFaces'; keys.append(key)
    values[key] = True
    names[key] = 'AddFacesOnNoSolid'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddCrvs'; keys.append(key)
    values[key] = True
    names[key] = 'AddCrvsOnNoSolid'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """
    
    # Prepare lists for collecting.
    gBreps = []; idx_rgFs_perB = []
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select faces")
    #go.SetCommandPromptDefault("Enter when done")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface
    #go.GeometryAttributeFilter = (
    #        ri.Custom.GeometryAttributeFilter.SubSurface)

    go.AcceptNothing(True)

    idxs_Opts = {}

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


    def collectSelection():
        for rdObjRef in go.Objects():
            idx_rgFace = rdObjRef.GeometryComponentIndex.Index
            gBrep = rdObjRef.ObjectId
                
            if not gBrep in gBreps:
                gBreps.append(gBrep)
                idx_rgFs_perB.append(set([idx_rgFace]))
            else:
                if idx_rgFace in idx_rgFs_perB[gBreps.index(gBrep)]:
                    idx_rgFs_perB[gBreps.index(gBrep)].remove(idx_rgFace)
                else:
                    idx_rgFs_perB[gBreps.index(gBrep)].add(idx_rgFace)

    def unhighlightAllFacesInAllSelectedBreps():
        for gB in gBreps:
            rdB = sc.doc.Objects.FindId(gB) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gB)
            rdB.UnhighlightAllSubObjects()
        sc.doc.Views.Redraw()


    while True:
        idxs_Opts['ChainSelectFillets'] = go.AddOption('ChainSelectFillets')
        Opts.riAddOpts['bHeal'](go)
        if Opts.values['bHeal']:
            Opts.riAddOpts['bAddModFaces'](go)
            Opts.riAddOpts['bAddCrvs'](go)
        Opts.riAddOpts['bEcho'](go)
        Opts.riAddOpts['bDebug'](go)

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            unhighlightAllFacesInAllSelectedBreps()
            gBreps = []; idx_rgFs_perB = []
            break
        elif res == ri.GetResult.Nothing:
            break
        elif res == ri.GetResult.Object:
            collectSelection()
            break
        elif res == ri.GetResult.Option:

            # Collect any faces already selected.
            collectSelection()

            if go.Option().Index == idxs_Opts['ChainSelectFillets']:
                sc.doc.Objects.UnselectAll()
                sc.doc.Views.Redraw()

                rc = xInput.getChainFilletFaces()
                if rc is not None:
                    gBrep, idx_rgFaces_Chain = rc
                    if not gBrep in gBreps:
                        gBreps.append(gBrep)
                        idx_rgFs_perB.append(set(idx_rgFaces_Chain))
                    else:
                        idx_rgFs_perB[gBreps.index(gBrep)] |= set(idx_rgFaces_Chain)

            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()
        
        sc.doc.Objects.UnselectAll() # Necessary when go.Get() is repeated.
        
        # Highlight all faces in selections.
        for b, gB in enumerate(gBreps):
            rdB = sc.doc.Objects.FindId(gB) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gB)
            rgFs = rdB.BrepGeometry.Faces
            for idxF in idx_rgFs_perB[b]:
                rdB.HighlightSubObject(rgFs[idxF].ComponentIndex(), True)
            rgFs[idxF].Dispose()
        sc.doc.Views.Redraw()
    
    go.Dispose()
    
    unhighlightAllFacesInAllSelectedBreps()
    
    if len(gBreps) == 0:
        return
    
    # Convert nested sets into lists.
    for i, idx_rgFs in enumerate(idx_rgFs_perB):
        idx_rgFs_perB[i] = list(idx_rgFs)
    
    return tuple(
            [gBreps] +
            [idx_rgFs_perB] +
            [Opts.values[key] for key in Opts.keys])


def processBrep(rgBrep0, idx_rgFs, bEcho=False, bDebug=False):
    """
    Parameters:
        rgBrep0
        idx_rgFs: Face indices to delete.
    Returns:
        Success:
            list(rgBreps to replace original geometry)
            list(Original rgBrepFaces adjacent to faces deleted)
            list(Modified rgBrepFaces adjacent to faces deleted)
            list(rgCurves of unaffected edges of adjacent faces)
        Failure: None
    """


    rgB0 = rgBrep0
    idx_Fs0_toDel = idx_rgFs



    def createOrderedListsOfBrepComponents(rgBrep0, idx_rgFaces0_Sel):
        """
        Per edges that border input BrepFace indices,
        create ordered lists of adjacent face, loop, and trim indices.

        Returns list(BrepFace indices), list(BrepLoop indices), list(BrepTrim indices)
        """
    
        idx_Fs0_Sel = idx_rgFaces0_Sel

        idx_AdjFs = []
        idx_Ls = []
        idx_Ts = []
        idx_Es0_All = []
        idx_Es0_NotSharedBySel = []
    
        # Loop through each face of selected trims.
        for iF, idx_F0_Sel in enumerate(idx_Fs0_Sel):
            rgF_Sel = rgBrep0.Faces[idx_F0_Sel]
            for idx_E in rgF_Sel.AdjacentEdges():
                idx_Es0_All.append(idx_E)

        for idx_E in set(idx_Es0_All):
            if idx_Es0_All.count(idx_E) == 1:
                idx_Es0_NotSharedBySel.append(idx_E)
            else:
                # Shared by 2 select faces or is a seam.
                pass

        for idx_rgEdge0_NotShared in idx_Es0_NotSharedBySel:
            idx_rgTrims_ToCheck = rgBrep0.Edges[idx_rgEdge0_NotShared].TrimIndices()
            if len(idx_rgTrims_ToCheck) > 1:
                for idx_rgTrim in idx_rgTrims_ToCheck:
                    idx_rgFace = rgBrep0.Trims[idx_rgTrim].Face.FaceIndex
                    if idx_rgFace not in idx_Fs0_Sel:
                        rgTrim = rgBrep0.Trims[idx_rgTrim]
                        idx_rgFace = rgTrim.Face.FaceIndex
                        idx_rgLoop = rgTrim.Loop.LoopIndex
                        if not idx_rgFace in idx_AdjFs:
                            idx_AdjFs.append(idx_rgFace)
                            idx_Ls.append([idx_rgLoop])
                            idx_Ts.append([[idx_rgTrim]])
                        else:
                            iF = idx_AdjFs.index(idx_rgFace)
                            if not idx_rgLoop in idx_Ls[iF]:
                                idx_Ls[iF].append(idx_rgLoop)
                                idx_Ts[iF].append([idx_rgTrim])
                            else:
                                iL = idx_Ls[iF].index(idx_rgLoop)
                                idx_Ts[iF][iL].append(idx_rgTrim)
    
        return idx_AdjFs, idx_Ls, idx_Ts


    def curvesForTrim(rgB_ToSplit, rgBs_Splitters_Now, fIntrsctTol=0.1*sc.doc.ModelAbsoluteTolerance):
        rgFace = rgB_ToSplit.Faces[0]
        rgCrvs_ForTrim = []
        fDevTol = sc.doc.ModelAbsoluteTolerance
        for rgB_Splitter in rgBs_Splitters_Now:
            if rgB_Splitter.Faces.Count == 1:
                # Get brep intersections with single-face splitters.
                stopwatch.Restart()
                b, rgCrvs, rgPts = rg.Intersect.Intersection.BrepBrep(
                        rgB_ToSplit, rgB_Splitter, fIntrsctTol)
                stopwatch.Stop()
                timeElapsed = stopwatch.Elapsed.TotalSeconds
                if timeElapsed > timeLimit:
                    print "{:.2f} seconds for " \
                        "Intersect.Intersection.BrepBrep".format(
                        timeElapsed)
                if b and rgCrvs.Count > 0:
                    rgCrvs_ForTrim.extend(rgCrvs)
            elif rgB_Splitter.Faces.Count > 1:
                rgCrvs = rgB_Splitter.DuplicateNakedEdgeCurves(True, True)
                if rgCrvs.Count > 0:
                    # Get intersections of brep with polyface splitter edges.
                    # Only accept curves that lie on brep.
                    for cA in rgCrvs:
                        # Due to Curve.PullToBrepFace possibly hanging Rhino,
                        # first test with Curve.ClosestPoints.
                        stopwatch.Restart()
                        b, ptOnCrv, ptOnObj, iObj = cA.ClosestPoints([rgB_ToSplit])
                        stopwatch.Stop()
                        timeElapsed = stopwatch.Elapsed.TotalSeconds
                        if timeElapsed > timeLimit:
                            print "{:.2f} seconds for " \
                                "Curve.ClosestPoints".format(
                                timeElapsed)
                            #if b:
                                #sc.doc.Objects.AddCurve(cA)
                                #sc.doc.Objects.AddBrep(rgB_ToSplit)
                        if not b: continue
                        rgCrvs_Pulled = cA.PullToBrepFace(rgFace, fIntrsctTol)
                        if rgCrvs_Pulled.Count == 1: # Can > 1 curves be useful?
                            cB = rgCrvs_Pulled[0]
                            rc = rg.Curve.GetDistancesBetweenCurves(cA, cB, fDevTol)
                            if rc[0] and rc[1] <= fDevTol: rgCrvs_ForTrim.append(cB)
                else:
                    print "Not!"
        return rgCrvs_ForTrim


    def joinBrepsSimultaneously(rgBreps_toJoin, iCt_NakedBorders):
        """
        Returns:
            list of Breps
            None on failure
        """

        iCt_NB0 = iCt_NakedBorders

        fTol_Join = max(
            1.001 * max([edge.Tolerance for rgB in rgBreps_toJoin for edge in rgB.Edges]),
            2.1 * sc.doc.ModelAbsoluteTolerance,
        )


        rgBs_Joined = rg.Brep.JoinBreps(
                rgBreps_toJoin,
                tolerance=fTol_Join)

        if not rgBs_Joined:
            if bEcho: print "Join failed."
            return
        
        for rgB in rgBs_Joined:
            if not rgB.IsValid:
                print "Invalid brep from JoinBreps"
                for rgB in rgBs_Joined: rgB.Dispose()
                return

        # Check if all are solid for quick return from this function.
        for rgB in rgBs_Joined:
            if not rgB.IsSolid: break
        else:
            # All breps are solids.
            return rgBs_Joined


        # In case some naked edges need some coercing ...
        for iB in range(len(rgBs_Joined)):
            if not rgBs_Joined[iB].IsSolid:
                rgB_JNE = rgBs_Joined[iB].DuplicateBrep()

                iCt_Joins = rgB_JNE.JoinNakedEdges(fTol_Join)

                if iCt_Joins:
                    print "JoinNakedEdges produced {} joins.".format(
                        iCt_Joins)

                if rgB_JNE.IsSolid:
                    print "Brep is now a solid."

                if not rgB_JNE.IsValid:
                    print "But Brep is no longer valid." \
                        "  Previous will be used instead."
                    rgB_JNE.Dispose()
                else:
                    rgBs_Joined[iB].Dispose()
                    rgBs_Joined[iB] = rgB_JNE


        # Check if all are solid.
        for rgB in rgBs_Joined:
            if not rgB.IsSolid: break
        else:
            # All breps are solids.
            return rgBs_Joined


        # Check whether the resultant brep(s) has no more naked borders than the input count.
        iCt_NB1_All = 0

        for rgB in rgBs_Joined:

            if rgB.IsSolid: continue

            idx_Es_NB1 = xBrep_nakedEdgeLoop.getEdgeIndicesOfBorders(rgB, None)
            iCt_NB = len(idx_Es_NB1) # 0 == Solid.
            iCt_NB1_All += iCt_NB

        if iCt_NB1_All > iCt_NB0:
            for rgB in rgBs_Joined: rgB.Dispose()
            return


        return rgBs_Joined


    def joinBrepsSingly(rgBrep_Main, rgBreps_Others, bDebug=False):
        """
        Parameters:
        Returns: New brep on success.
        """

        fTol_Join = max(
            1.001 * max([edge.Tolerance for rgB in [rgBrep_Main] + rgBreps_Others for edge in rgB.Edges]),
            2.1 * sc.doc.ModelAbsoluteTolerance,
        )


        # Make working copies of the original input.
        rgB0_Joined_WIP = rgBrep_Main.DuplicateBrep()
        rgBreps_Others_Copy = rgBreps_Others[:]
    
        # If rgB0_Joined_WIP doesn't contain any faces,
        # Append the first brep of rgBreps_Others_Copy to it.
        if rgB0_Joined_WIP.Faces.Count == 0:
            rgB0_Joined_WIP.Append(rgBreps_Others_Copy[0])
            del rgBreps_Others_Copy[0]

        while True:
            for i, rgBrep_Try in enumerate(rgBreps_Others_Copy):
                bSuccess = rgB0_Joined_WIP.Join(
                    rgBrep_Try, fTol_Join,
                    compact=True)
                # compact=False sometimes creates an invalid brep
                # that doesn't become valid with a later Compact() call.
                if bSuccess:
                    # Successful join, so remove that face from the ID and RG lists
                    # and restart joining with the modified RG list.
                    del rgBreps_Others_Copy[i]
                
                    #sc.doc.Objects.AddBrep(rgB0_Joined_WIP); sc.doc.Views.Redraw(); #1/0
                
                    break # to restart for loop.
            else:
              # Join for each remaining brep to the main has been attempted
              # in an uninterrupted for loop.
              break # out of while loop.
    
        b, sLog = rgB0_Joined_WIP.IsValidWithLog()
        if not b:
            rgB0_Joined_WIP.Dispose()
            print sLog
            print '-'*40
            b, sLog = rgB1_1F.IsValidTopology()
            print sLog
            return

        #sc.doc.Objects.AddBrep(rgB0_Joined_WIP); sc.doc.Views.Redraw()#; 1/0

        if rgB0_Joined_WIP is None:
            if bEcho: print "Join failed."
            return

        if not rgB0_Joined_WIP.IsSolid:
            if bEcho:
                print "Failed to heal into a solid after trimming faces."
            rgB0_Joined_WIP.JoinNakedEdges(2.1*sc.doc.ModelAbsoluteTolerance)
            if rgB0_Joined_WIP.IsSolid:
                if bEcho: print "JoinNakedEdges resulted in a solid."
            else:
                if bEcho: print "JoinNakedEdges didn't result in a solid."
                rgB0_Joined_WIP.Dispose()
                return

        return rgB0_Joined_WIP


    def trim1FaceBrepsWithEachOther(rgBs_1F_ToTrim, ptsOnFaces, rgB_Base, bEcho=False, bDebug=False):
        """
        Parameters:
        Returns:
        """

        rgBs_1F_Out = []
        iCtIntersectFails = 0
    
        for iB in range(len(rgBs_1F_ToTrim)):
            rgB_ToSplit = rgBs_1F_ToTrim[iB].DuplicateBrep()
            rgBs_Splitters_Now = rgBs_1F_ToTrim[0:iB] + rgBs_1F_ToTrim[iB:] + rgB_Base
        
            # Create curves from intersections of face with splitters.
            stopwatch.Restart()
            rgCrvs_ForTrim = curvesForTrim(rgB_ToSplit, rgBs_Splitters_Now)
            stopwatch.Stop()
            timeElapsed = stopwatch.Elapsed.TotalSeconds
            if timeElapsed > timeLimit:
                print "{:.2f} seconds for " \
                    "curvesForTrim" \
                    " for brep {} of {}".format(
                        timeElapsed, iB, len(rgBs_1F_ToTrim))
        #            if len(rgCrvs_ForTrim) > 1:
        #                sc.doc.Objects.AddBrep(rgB_ToSplit)
        #                map(sc.doc.Objects.AddCurve, rgCrvs_ForTrim)
        
            if len(rgCrvs_ForTrim) == 0:
                iCtIntersectFails += 1
                if bEcho:
                    print "No intersections found for 1-face brep {} of {}.".format(
                            iB+1, len(rgBs_1F_ToTrim))
                rgBs_1F_Out.append(rgB_ToSplit)
                continue
        
            #if bDebug: map(sc.doc.Objects.AddCurve, rgCrvs_ForTrim)
        
            # Split face using intersection curves.
            rgF_ToSplit = rgB_ToSplit.Faces[0]
            rgB_Split = rgF_ToSplit.Split(rgCrvs_ForTrim,
                    sc.doc.ModelAbsoluteTolerance)
            if rgB_Split is None:
                rgBs_1F_Out.append(rgB_ToSplit)
                continue
        
            if rgB_Split.Faces.Count == 1:
                rgB_Split.Dispose()
                rgBs_1F_Out.append(rgB_ToSplit)
                continue
        
            #  Find correct face in split brep.
            idx_rgFace_Pos = xBrep_findMatchingFace.usingPointOnFace(
                    rgB_Split, ptsOnFaces[iB])
            if idx_rgFace_Pos is None:
                rgBs_1F_Out.append(rgB_ToSplit)
                continue
            rgB_1F_ForReplace = rgB_Split.Faces[
                    idx_rgFace_Pos].DuplicateFace(False)
            rgBs_1F_Out.append(rgB_1F_ForReplace)
    
        return rgBs_1F_Out, iCtIntersectFails


    def isoStatusesToExtend(idx_AdjFs, idx_Ls_perAdjF, idx_Ts_ToRemove_perL_perAdjF):
        """
        Returns: List of lists of isoStatuses per adjacent face.
        """
        # Use Surface.ClosestSide to determine SENW Trim(s) of each surface to extend.
        rgB0_Shrunk = rgB0.DuplicateBrep()
        rgB0_Shrunk.Faces.ShrinkFaces()
        isoStats_Closest_OnFaces = []

        for iF, idx_AdjF in enumerate(idx_AdjFs):
            rgFace = rgB0_Shrunk.Faces[idx_AdjF]
            isoStats_Closest_OnFaces.append([])
            for iL, idx_rgLoop in enumerate(idx_Ls_perAdjF[iF]):
                for idx_rgTrim_ToRemove in idx_Ts_ToRemove_perL_perAdjF[iF][iL]:
                    rgTrim = rgB0_Shrunk.Trims[idx_rgTrim_ToRemove]
                    t = rgTrim.DivideByCount(2, False)[0]
                    ptTrim = rgTrim.PointAt(t)
                    #pt3d = rgFace.PointAt(ptTrim.X, ptTrim.Y)
                    #sc.doc.Objects.AddPoint(pt3d)
                    isoStat_Closest_OnFace = rgFace.ClosestSide(ptTrim.X, ptTrim.Y)
                    if isoStat_Closest_OnFace not in isoStats_Closest_OnFaces[iF]:
                        isoStats_Closest_OnFaces[iF].append(isoStat_Closest_OnFace)

        rgB0_Shrunk.Dispose()

        return isoStats_Closest_OnFaces


    def extend1FaceBrep(rgB1_1F, isoStats_Closest_OnFace):
        
        rgSrf_Ext = rgB1_1F.Faces[0].UnderlyingSurface()
        
        # Extend all isostats of interest for each 1-face brep.
        for isostat in isoStats_Closest_OnFace:
            # Create extended surface.
            rgSrf_Ext = rgSrf_Ext.Extend(isostat, fExtLength, True) # True for smooth.
        
        # Create indices of Trims with which not to trim the extended surface.
        idxTrims_toSkip = []
        for idx, rgTrim in enumerate(rgB1_1F.Trims):
            if any(rgTrim.IsoStatus == iso for iso in
                    isoStats_Closest_OnFace):
                idxTrims_toSkip.append(idx)
        
        stopwatch.Restart()
        rgB1_1F_Adj_Ext = xBrepFace_extendAtUntrimmedEdge.trimExtendedSrf(
                rgSrf_Ext,
                rgB1_1F,
                idxTrims_toSkip,
                bEcho=bDebug,
                bDebug=bDebug)
        stopwatch.Stop()
        timeElapsed = stopwatch.Elapsed.TotalSeconds
        if timeElapsed > timeLimit:
            print "{:.2f} seconds for " \
                "xBrepFace_extendAtUntrimmedEdge.trimExtendedSrf".format(
                timeElapsed)
            #sc.doc.Objects.AddBrep(rgB1_1F_Adj_Ext)
        
        return rgB1_1F_Adj_Ext


    def extend1FaceBreps(rgBs1_1F_Adj_ForExt, isoStats_Closest_OnFaces):
        
        rgBs_Adj_Extended = []
        
        for iB in range(len(rgBs1_1F_Adj_ForExt)):
            rgB_In = rgBs1_1F_Adj_ForExt[iB]

            #sEval = 'rgB_In.IsSolid'; print sEval+':',eval(sEval)
            if rgB_In.IsSolid:
                rgBs_Adj_Extended.append(rgB_In.DuplicateBrep())
            elif not xBrepFace_extendAtUntrimmedEdge.isBrepReadyForExtend(
                    rgB_In, bEcho=bEcho, bDebug=bDebug):
                if bEcho:
                    s = "Brep {} of rgBs1_1F_Adj_ForExt is not ready to be " \
                          "extended and will be skipped.".format(iB)
                    print s
                #sc.doc.Objects.AddBrep(rgB1_1F)
                rgBs_Adj_Extended.append(rgB_In.DuplicateBrep())
            else:
                rgB_Out = extend1FaceBrep(rgB_In, isoStats_Closest_OnFaces[iB])
                
                if rgB_Out is not None:
                    rgBs_Adj_Extended.append(rgB_Out)
                    #sc.doc.Objects.AddBrep(rgB_Out); sc.doc.Views.Redraw()
        
        return rgBs_Adj_Extended



    # Collect into lists the faces, loops, and trims
    # of adjacent faces that share edges of faces to remove.
    rc = createOrderedListsOfBrepComponents(
        rgB0,
        idx_Fs0_toDel)
    (
        idx_AdjFs,
        idx_Ls_perAdjF,
        idx_Ts_ToRemove_perL_perAdjF
    ) = rc

    if bDebug:
        sEval = 'idx_AdjFs'; print sEval+':',eval(sEval)
        sEval = 'idx_Ls_perAdjF'; print sEval+':',eval(sEval)
        sEval = 'idx_Ts_ToRemove_perL_perAdjF'; print sEval+':',eval(sEval)


    # Remove selected and adjacent faces from brep.
    # It is OK if brep has no faces left.
    idx_rgFaces_ToRemove = idx_Fs0_toDel + idx_AdjFs
    idx_rgFaces_ToRemove.sort(reverse=True)
    rgB_Unaffected = rgB0.DuplicateBrep()
    for idx_rgFace_ToRemove in idx_rgFaces_ToRemove:
        rgB_Unaffected.Faces.RemoveAt(idx_rgFace_ToRemove)

    rgBs_Adj_0 = [rgB0.Faces[idx].DuplicateFace(True) for idx in idx_AdjFs]


    print "Phase 1: Untrim adjacent faces and join all."


    # From adjacent faces, remove trim common with those of faces to remove.
    rgBs_Adj_RemovedSharedT = []
    rgCrvs1_Joined = []

    for iF, idx_AdjF in enumerate(idx_AdjFs):
        
        if iF == 3:
            pass
        
        rc = xBrepFace_removeTrim.removeTrim(
                rgB0,
                idx_AdjF,
                idx_Ls_perAdjF[iF],
                idx_Ts_ToRemove_perL_perAdjF[iF],
                bDebug=bDebug)

        if rc is None or rc[0] is None:
            print "f[{}] not untrimmed!".format(idx_AdjF)
            for rgB in rgBs_Adj_RemovedSharedT: rgB.Dispose()
            return [rgB_Unaffected], rgBs_Adj_0, [], []

        rgB1_1F, rgEdges_NotSelButInLoopsWithSel = rc
        
        rgBs_Adj_RemovedSharedT.append(rgB1_1F)
        
        # Join curves.
        if len(rgEdges_NotSelButInLoopsWithSel) > 0:
            rgCrvs1_Joined.extend(Rhino.Geometry.Curve.JoinCurves(
                    rgEdges_NotSelButInLoopsWithSel))


    # Get point on each original adjacent face to later use to find correct faces of split brep.
    ptsOnFaces = []
    for idx_AdjF in idx_AdjFs:
        if idx_AdjF == 36:
            pass
        
        pt = xBrepFace.createPoint3dOnInterior(rgB0.Faces[idx_AdjF])
        if pt is None:
            print "Point on face not found!" \
                "Only non-modified adjacent faces will be added."
            for rgB in rgBs_Adj_RemovedSharedT: rgB.Dispose()

            sc.doc.Objects.AddBrep(rgB0.Faces[idx_AdjF].DuplicateFace(False))
            sc.doc.Views.Redraw(); return

            return [rgB_Unaffected], rgBs_Adj_0, [], []

        ptsOnFaces.append(pt)
    


    # Track number of naked borders of initial brep to determine
    # expectations of naked borders of final breps.
    idx_Es_nakedBorders0 = xBrep_nakedEdgeLoop.getEdgeIndicesOfBorders(rgB0, None)
    iCt_NakedBorders = len(idx_Es_nakedBorders0) # 0 == Solid.


    rgBs_Joined = joinBrepsSimultaneously(
        [rgB_Unaffected] + rgBs_Adj_RemovedSharedT,
        iCt_NakedBorders)
    if rgBs_Joined is not None:
        if bEcho:
            print "Brep healed after " \
                "untrimming adjacent faces."

        rgB_Unaffected.Dispose()
        for rgB in rgBs_Adj_0: rgB.Dispose()
        for rgB in rgBs_Adj_RemovedSharedT: rgB.Dispose()

        return rgBs_Joined, [], [], rgCrvs1_Joined

    #rgB_Joined = joinBrepsSingly(
    #    rgB_Unaffected,
    #    rgBs_Adj_RemovedSharedT,
    #    bDebug=bDebug)
    #if rgB_Joined is not None:
    #    if bEcho: print "Brep healed into a solid after untrimming adjacent faces."
    #    rgB_Unaffected.Dispose()
    #    for rgB in rgBs_Adj_RemovedSharedT: rgB.Dispose()
    #    return rgB_Joined, rgBs_Adj_0, [], rgCrvs1_Joined

    if bEcho: print "Failed to heal after untrimming adjacent faces."


    print "Phase 2: Trim adjacent faces and join all."

    if bEcho: print "Trimming adjacent faces with one another..."
    rc = trim1FaceBrepsWithEachOther(
            rgBs_Adj_RemovedSharedT,
            ptsOnFaces=ptsOnFaces,
            rgB_Base=[],
            bEcho=bEcho,
            bDebug=bDebug)
    rgBs_Adj_NewTrim, iCtIntersectFails = rc
    
    if iCtIntersectFails > 0:
        if bEcho:
            print "Intersects for trimming not found for {} out of {} faces.".format(
                    iCtIntersectFails, len(rgBs_Adj_NewTrim))


    rgBs_Joined = joinBrepsSimultaneously(
        [rgB_Unaffected] + rgBs_Adj_NewTrim,
        iCt_NakedBorders)
    if rgBs_Joined is not None:
        if bEcho:
            print "Brep healed after " \
                "trimming adjacent faces."

        rgB_Unaffected.Dispose()
        for rgB in rgBs_Adj_0: rgB.Dispose()
        for rgB in rgBs_Adj_NewTrim: rgB.Dispose()

        return rgBs_Joined, [], [], rgCrvs1_Joined

    #rgB_Joined = joinBrepsSingly(
    #    rgB_Unaffected,
    #    rgBs_Adj_NewTrim,
    #    bDebug=bDebug)
    #if rgB_Joined is not None:
    #    if bEcho: print "Brep healed into a solid after trimming adjacent faces."
    #    rgB_Unaffected.Dispose()
    #    for rgB in rgBs_Adj_NewTrim: rgB.Dispose()
    #    return rgB_Joined, rgBs_Adj_0, [], rgCrvs1_Joined

    for rgB in rgBs_Adj_NewTrim: rgB.Dispose()

    if bEcho: print "Failed to heal brep after trimming adjacent faces."


    print "Phase 3: Extend the adjacent faces and trim again."

    ## Shrink planar, cylindrical, and conical surfaces before extend.
    #rgBs1_1F_Adj_ForExt = []
    #tol = sc.doc.ModelAbsoluteTolerance
    #for rgB1_1F_Adj in rgBs1_1F_Adj:
    #    rgB1_1F_Adj_ForExt = rgB1_1F_Adj.DuplicateBrep()
    #    rgF = rgB1_1F_Adj_ForExt.Faces[0]
    #    if rgF.IsPlanar(tol) or rgF.IsCylinder(tol) or rgF.IsCone(tol):
    #        print rgB1_1F_Adj_ForExt.Faces.ShrinkFaces()
    #    rgBs1_1F_Adj_ForExt.append(rgB1_1F_Adj_ForExt)


    isoStats_Closest_OnFaces = isoStatusesToExtend(
        idx_AdjFs,
        idx_Ls_perAdjF,
        idx_Ts_ToRemove_perL_perAdjF)
    
    # For surface extension length,
    # use diagonal dimension of bounding box of faces to delete.
    rgBs_Sel = map(lambda x: rgB0.Faces[x].DuplicateFace(True), idx_Fs0_toDel)
    bbox = rgBs_Sel[0].GetBoundingBox(True)
    rgBs_Sel[0].Dispose()
    if rgBs_Sel.Count > 1:
        for i in range(1, len(rgBs_Sel)):
            bbox.Union(rgBs_Sel[i].GetBoundingBox(True))
            rgBs_Sel[i].Dispose()
    fExtLength = bbox.Min.DistanceTo(bbox.Max)

    rgBs_Adj_Extended = extend1FaceBreps(rgBs_Adj_RemovedSharedT, isoStats_Closest_OnFaces)

    rc = trim1FaceBrepsWithEachOther(
            rgBs_Adj_Extended,
            ptsOnFaces,[rgB_Unaffected],
            bEcho=bEcho,
            bDebug=bDebug)
    rgBs_Adj_NewTrim, iCtIntersectFails = rc

    for rgB in rgBs_Adj_Extended: rgB.Dispose()

    #for rgB in rgBs_Adj_NewTrim: sc.doc.Objects.AddBrep(rgB)
    #sc.doc.Views.Redraw(); return

    rgBs_Joined = joinBrepsSimultaneously(
        [rgB_Unaffected] + rgBs_Adj_NewTrim,
        iCt_NakedBorders)
    if rgBs_Joined is not None:
        if bEcho:
            print "Brep healed after " \
                "extending and trimming adjacent faces."

        rgB_Unaffected.Dispose()
        for rgB in rgBs_Adj_0: rgB.Dispose()
        for rgB in rgBs_Adj_NewTrim: rgB.Dispose()

        return rgBs_Joined, [], [], rgCrvs1_Joined

    #rgB_Joined = joinBrepsSingly(
    #    rgB_Unaffected,
    #    rgBs_Adj_NewTrim,
    #    bDebug=bDebug)
    #if rgB_Joined is not None:
    #    if bEcho: print "Brep healed into a solid after extending and trimming adjacent faces."
    #    rgB_Unaffected.Dispose()
    #    for rgB in rgBs_Adj_NewTrim: rgB.Dispose()
    #    return rgB_Joined, rgBs_Adj_0, [], rgCrvs1_Joined

    for rgB in rgBs_Adj_NewTrim: rgB.Dispose()

    if bEcho: print "Failed to heal after extending and trimming adjacent faces."

    return [rgB_Unaffected], [], rgBs_Adj_RemovedSharedT, rgCrvs1_Joined


def processBrepObjects(gBreps0, idx_rgFs_perB0, bHeal=True, bAddModFaces=True, bAddCrvs=True, bEcho=False, bDebug=False):
    """
    Parameters:
        idBreps0 (optional)
        idx_rgFs_perB (optional): Faces to be deleted
    Returns: None
    """
    
    Rhino.RhinoApp.SetCommandPrompt("Working ...")
    
    # Loop through breps.
    for iB0 in range(len(gBreps0)):
        gBrep0 = gBreps0[iB0]
        rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
        rgB0 = rdBrep0.BrepGeometry # Do not modify this brep.

        if not rgB0.IsValid:
            if bEcho: print "Invalid brep skipped."
            rgB0.Dispose()
            continue

        idx_rgFs0_thisB = idx_rgFs_perB0[iB0] # Faces of 1 brep, so no nesting.
        
        if rgB0.Faces.Count == len(idx_rgFs0_thisB):
            sc.doc.Objects.Delete(gBrep0, False)
            rgB0.Dispose()
            if bEcho:
                print "All faces of brep were selected," \
                    " so the entire brep was deleted."
            continue
        
        if not bHeal:
            xBrepObject.removeFaces(rdBrep0, idx_rgFs0_thisB)
            continue

        rc = processBrep(
            rgBrep0=rgB0,
            idx_rgFs=idx_rgFs0_thisB,
            bEcho=bEcho,
            bDebug=bDebug)
        rgB0.Dispose()
        
        if not rc: continue

        (
            rgB1s_forReplacement,
            rgBs_Adj_0,
            rgBs_Adj_Mod,
            rgCrvs_RemainingTsOfAdjFs,
        ) = rc

        if len(rgB1s_forReplacement) == 1 and not rgBs_Adj_0 and not rgBs_Adj_Mod:
            if not xBrepObject.replaceGeometry(rdBrep0, rgB1s_forReplacement):
                if bEcho: print "Replacing original brep failed."
        else:
            if bAddModFaces and rgBs_Adj_Mod:
                for iB in range(len(rgBs_Adj_Mod)):
                    idB1_1F = sc.doc.Objects.AddBrep(
                        rgBs_Adj_Mod[iB], rdBrep0.Attributes)
                if bEcho:
                    print "Modified " \
                        "faces adjacent with those removed have been extracted."
            elif rgBs_Adj_0:
                for iB in range(len(rgBs_Adj_0)):
                    idB1_1F = sc.doc.Objects.AddBrep(
                        rgBs_Adj_0[iB], rdBrep0.Attributes)
                if bEcho:
                    print "Unmodified " \
                        "faces adjacent with those removed have been extracted."

            # Replace original brep with one (or more) with removed faces.
            xBrepObject.replaceGeometry(rdBrep0, rgB1s_forReplacement)

            if bAddCrvs and (rgBs_Adj_0 or rgBs_Adj_Mod):
                for rgCrv in rgCrvs_RemainingTsOfAdjFs:
                    # Using default attributes so curves are more noticeable
                    # than the untrimmed faces.
                    sc.doc.Objects.AddCurve(rgCrv)
                sc.doc.Views.Redraw()

        for rgB in rgB1s_forReplacement: rgB.Dispose()
        for rgB in rgBs_Adj_0: rgB.Dispose()
        for rgB in rgBs_Adj_Mod: rgB.Dispose()

    sc.doc.Views.Redraw()


def main():
    
    rc = getInput()
    if rc is None: return

    (
        gBreps0,
        idx_rgFs_perB,
        bHeal,
        bAddModFaces,
        bAddCrvs,
        bEcho,
        bDebug,
    ) = rc

    if gBreps0 is None or idx_rgFs_perB is None: return

    if Opts.values['bDebug']:
        print "Running with Debug mode on."
        import sys
        for sModule in list(sys.modules):
            if sModule[0] == 'x':
                try:
                    reload(sys.modules[sModule])
                    print "{} reloaded.".format(sModule)
                except:
                    s  = "{} NOT reloaded.".format(sModule)
                    s += "  Does the module contain a bug"
                    s += " or was its name changed?"
                    print s

    processBrepObjects(
        gBreps0=gBreps0,
        idx_rgFs_perB0=idx_rgFs_perB,
        bHeal=bHeal,
        bAddModFaces=bAddModFaces,
        bAddCrvs=bAddCrvs,
        bEcho=bEcho,
        bDebug=bDebug)


if __name__ == '__main__': main()