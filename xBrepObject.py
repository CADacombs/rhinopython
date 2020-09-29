"""
181025: Created.
...
190810: replaceFaces now accepts surfaces as input for the new geometry.
        Added fTolerance_Join parameter for replaceFaces.  If one is not provided, the current brep's maximum edge tolerance is used.
190827: CANCELED.
190828: Maximum of 1.01 * brep's edge tolerance and input tolerance is now used for JoinBreps.
191010: Fixed bug where duplicate values in idxFaces caused some faces to be missing from brep with faces removed.
        Changed some parameter names.  Changed a printed output to bDebug.
191015: Corrected a variable name in 1 function.
191028: replaceFaces can now accept more types of input for breps and face indices.
191029: Refactored replaceGeometry.
191029-30: Created replaceFaces_WIP, a replacement for replaceFaces that allows the new geometry to be plural.
191030: Bug fix in extractFaces.
191031: Refactored replaceFaces.
191101: removeFaces now returns [] instead of True when all faces are removed.
        extractFaces now returns a tuple of (list of extracted), (list of remaining) instead of just a list of extracted Breps.
191117: CreateBooleanUnion for separating "shells" now accepts non-manifold breps.
191118: Modified determination of joining tolerances in 2 functions.
191208: Replaced an xrange with range so that the reversed result can be passed to length
191213: Bug fixed by changing JoinBrep tolerance value in addFromSubsetOfFaces.
200107: replaceGeometry now returns None when some of the new geometry is invalid.
200108: Modified how addFromSubsetOfFaces handles invalid breps created by BrepFace.DuplicateFace.
200129: Corrected a docstring.
200527: Modified conditionals for some printed output.
200619: Added more functions.  Purged some of this history.
200701: Now no face indices passed to extractFaces creates an _Explode.  Various minor changes.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Diagnostics import Stopwatch

import xBrepFace


def coerceRhinoObject(rhObj):
    """
    'Deleted objects cannot be found by id.'
    (https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_DocObjects_Tables_ObjectTable_FindId.htm)
    """
    rdObj = None
    if isinstance(rhObj, rd.RhinoObject):
        rdObj = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        rdObj = rhObj.Object()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
    return rdObj


def coerceBrepObject(rhObj):
    rdObj = coerceRhinoObject(rhObj)
    if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
        return rdObj


def coerceBrepGuid(rhObj):
    if isinstance(rhObj, Guid): return rhObj
    rdBrep = coerceBrepObject(rhObj)
    if rdBrep:
        return rdBrep.Id


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.Brep):
        return rhObj
    elif isinstance(rhObj, rg.Surface):
        return rhObj.ToBrep()
    elif isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        #print rhObj.GeometryComponentIndex.ComponentIndexType
        geom = rhObj.Geometry()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    if isinstance(geom, rg.Brep):
        return geom


def addFromSubsetOfFaces(rhBrep, idxFaces, bAddOnlyMonofaces=True, bRetainLayer=True, bRetainColor=True, bDebug=False):
    """
    """
    
    rdBrep_In = coerceBrepObject(rhBrep)
    if rdBrep_In is None: return
    rgBrep_In = rdBrep_In.BrepGeometry
    if not rgBrep_In.IsValid: return
    
    
    def addBrepOfSubsetOfFaces_JoinBreps():
        
        # Duplicate faces to their own breps to be joined.
        rgBreps1 = [] # Faces (breps) to be duplicated.
        
        for i in idxFaces:
            rgFace = rgBrep_In.Faces[i]
            rgBrep_1Face = rgFace.DuplicateFace(True)
            if rgBrep_1Face is None:
                if bDebug: print "Face {} could not be duplicated as a brep!".format(i)
                return None
            rgBreps1.append(rgBrep_1Face)
        
        # Join monoface breps.

        # Using a tight tolerance to rejoin only existing shared edges.
        fTol_Join = 1e-9

        rgBreps_Joined = rg.Brep.JoinBreps(rgBreps1, tolerance=fTol_Join)
        if rgBreps_Joined is None:
                if bDebug: print "Joining breps failed!"
                return
        for rgB in rgBreps_Joined:
            if not rgB.IsValid:
                if bDebug: print "Joined brep not valid.  Exiting..."
                return
        
        if any(b.Faces.Count > 1 for b in rgBreps_Joined):
            # Separate any brep shells of modified brep and act based on shell quantity.
            rgBreps_per_shell = rg.Brep.CreateBooleanUnion(
                rgBreps_Joined, tolerance=0.0, manifoldOnly=False)
            if rgBreps_per_shell is None:
                if bDebug: print "Error in attempting to separate brep shells.  No objects have been modified."
                return
        else:
            # Skipped attempting to Boolean union monoface breps in case any contact are in contact with one another.
            rgBreps_per_shell = rgBreps_Joined[:]
        
        
        attr = rdBrep_In.Attributes.Duplicate()
        if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
        
        
        gBreps_1Shell = map(lambda x: sc.doc.Objects.AddBrep(x, attr), rgBreps_per_shell)
        map(lambda x: x.Dispose(), rgBreps_per_shell)
        
        return gBreps_1Shell
    
    
    def addBrepOfSubsetOfFaces_RemoveAt():
        
        rgBrep1 = rgBrep_In.Duplicate()
        
        # Create list of non-extracted faces.
        # This is faster than looping through all faces while testing each.
        idx_rgFaces_ToRemove = list(
                set(range(rgBrep1.Faces.Count)) - set(idxFaces))
        idx_rgFaces_ToRemove.sort(reverse=True)
        
        stopwatch = Stopwatch()
        stopwatch.Start()
        [rgBrep1.Faces.RemoveAt(idx) for idx in idx_rgFaces_ToRemove]
        stopwatch.Stop()
        if bDebug: print "{} seconds for Faces.RemoveAt".format(stopwatch.Elapsed.TotalSeconds)
        
        # Separate any brep shells of modified brep and act based on shell quantity.
        stopwatch.Restart()
        rgBreps_per_shell = rg.Brep.CreateBooleanUnion(
            [rgBrep1], tolerance=0.0, manifoldOnly=False)
        stopwatch.Stop()
        if bDebug:
            print "{} seconds for CreateBooleanUnion".format(
                stopwatch.Elapsed.TotalSeconds)
        if rgBreps_per_shell is None:
            if bDebug: print "Error in attempting to separate brep shells.  No objects have been modified."
            return
        
        rgBrep1.Dispose()
        
        attr = rdBrep_In.Attributes.Duplicate()
        if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
        
        stopwatch.Restart()
        gBreps_1Shell = map(lambda x: sc.doc.Objects.AddBrep(x, attr), rgBreps_per_shell)
        stopwatch.Stop()
        if bDebug: print "{:.1f} seconds for AddBrep".format(stopwatch.Elapsed.TotalSeconds)
        map(lambda x: x.Dispose(), rgBreps_per_shell)
        
        return gBreps_1Shell
    
    
    nFaces = rgBrep_In.Faces.Count
    
    # If brep has only 1 face, return the brep's GUID.
    if nFaces == 1:
        return [rdBrep_In.Id]
    
    idxFaces = list(set(idxFaces))
    
    if not bAddOnlyMonofaces:
        
        stopwatch = Stopwatch()
        
        # Create brep(s) of extracted faces using method chosen by whether number of
        # extracted faces is less or more than the remaining number of faces.
        if len(idxFaces) < nFaces // 2 :
            stopwatch.Restart()
            gBreps1_Extracted = addBrepOfSubsetOfFaces_JoinBreps()
            stopwatch.Stop()
            if bDebug:
                print "{} seconds for addBrepOfSubsetOfFaces_JoinBreps".format(
                        stopwatch.Elapsed.TotalSeconds)
            if gBreps1_Extracted is not None:
                return gBreps1_Extracted
            if bDebug: print "addBrepOfSubsetOfFaces_JoinBreps returned None."
        # Since addBrepOfSubsetOfFaces_JoinBreps failed, will try RemoveAt instead.
        
        # Either the number of faces to add > half the total number of faces in the brep or addBrepOfSubsetOfFaces_JoinBreps had returned None.
        stopwatch.Restart()
        gBreps1_Extracted = addBrepOfSubsetOfFaces_RemoveAt()
        stopwatch.Stop()
        if bDebug:
            print "{} seconds for addBrepOfSubsetOfFaces_RemoveAt".format(
                    stopwatch.Elapsed.TotalSeconds)
        if gBreps1_Extracted is not None:
            return gBreps1_Extracted
        if bDebug: print "addBrepOfSubsetOfFaces_RemoveAt returned None."
    
    # Add only monoface breps.
    
    attr = rdBrep_In.Attributes.Duplicate()
    if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
    
    gBreps1_Extracted = []
    
    for idx in idxFaces:
        rgFace = rgBrep_In.Faces[idx]
        
        # Duplicate face to its own brep.
        rgBrep1 = rgFace.DuplicateFace(duplicateMeshes=True)
        if not rgBrep1.IsValid:
            gBrep1 = None
        else:
            gBrep1 = sc.doc.Objects.AddBrep(rgBrep1, attr)

        if gBrep1 is None:
            s = "Brep face {} from {} could not be added to document.".format(
                    idx, rdBrep_In.Id)
            print s
            rc = rs.MessageBox(
                s + "\nContinue extracting faces, skipping this one?",
                buttons=4,
                title="xBrepObject.addFromSubsetOfFaces")
            if rc is not None and rc == 6:
                continue
            #if not bDebug: rs.DeleteObjects(gBreps1_Extracted)
            return
        
        gBreps1_Extracted.append(gBrep1)
        rgBrep1.Dispose()
    
    return gBreps1_Extracted


def replaceGeometry(rhBrep_toReplace, rgBreps_NewGeom):
    """
    Replaces BrepObject with new geometry after separating any 'shells'.
    Parameters:
        brep_toReplace: BrepObject, ObjRef, or GUID.
        rgBreps_NewGeom: Single or list of BrepObjects or Breps in any form.
    Returns:
        list of Brep GUIDs of remaining GUIDs.
        None if input doesn't qualify.
    """
    
    #if rgBrep1.Faces.Count == 0: return
    
    rdBrep_In = coerceBrepObject(rhBrep_toReplace)

    try: rgBreps_NewGeom = list(set(rgBreps_NewGeom))
    except: rgBreps_NewGeom = [rgBreps_NewGeom]

    rgBreps_NewGeom = [coerceBrep(o) for o in rgBreps_NewGeom]
    
    for rgB in rgBreps_NewGeom:
        if not rgB.IsValid:
            print "Brep geometry was not replaced because" \
                  " some of the new geometry is invalid."
            return

    # Separate any shells of modified brep and act based on shell quantity.
    rgBreps_per_shell = []
    for brepN in rgBreps_NewGeom:
        rgBs_Unioned = rg.Brep.CreateBooleanUnion(
            [brepN], tolerance=0.0, manifoldOnly=False)
        if rgBs_Unioned is None:
            rgBreps_per_shell.append(None)
            break # rgBreps_per_shell will be disregarded later.
        for rgB_Unioned in rgBs_Unioned:
            rgBreps_per_shell.append(rgB_Unioned)
            if rgB_Unioned is None:
                break # rgBreps_per_shell will be disregarded later.

    if None in rgBreps_per_shell:
        print "Error in attempting to separate 'shells'." \
            "  Attempting to replace brep as is."
        rgBreps_ReadyToReplace = [brep.Duplicate() for brep in rgBreps_NewGeom]
    else:
        rgBreps_ReadyToReplace = [brep.Duplicate() for brep in rgBreps_per_shell]

    for brep in rgBreps_per_shell: brep.Dispose()

    if len(rgBreps_ReadyToReplace) == 1:
        # Replace original brep with its modified geometry.
        # Starting GUID will be maintained.
        if sc.doc.Objects.Replace(rdBrep_In.Id, rgBreps_ReadyToReplace[0]):
            return [rdBrep_In.Id]
        else:
            print "Replace brep {} failed.  Check this.".format(rdBrep_In.Id)
    else:
        # Multiple breps will be added.
        # Original will be deleted, not replaced,
        # so that none will inherit the original GUID.
        attr = rdBrep_In.Attributes
        gBreps_Out = [
            sc.doc.Objects.AddBrep(brep, attr) for brep in rgBreps_ReadyToReplace]
        if not any(g == Guid.Empty for g in gBreps_Out):
            if not sc.doc.Objects.Delete(rdBrep_In, True):
                print "Original brep cannot be deleted.  Check this"
        return gBreps_Out


def removeFaces(rhBrep, idxs_rgFaces):
    """
    rhBrep: GUID or BrepObject.
    idxs_rgFaces: Iterable of face indices.
    
    Returns:
        Iterable of GUIDs of breps with faces removed.
        None on fail.
    """

    try: idxs_rgFaces = sorted(list(set(idxs_rgFaces)), reverse=True)
    except: idxs_rgFaces = [idxs_rgFaces]

    rdBrep_In = coerceBrepObject(rhBrep)
    if rdBrep_In is None: return
    rgBrep_In = rdBrep_In.Geometry
    if not rgBrep_In.IsValid and rgBrep_In.Faces.Count > 1: return
    
    if rgBrep_In.Faces.Count == len(idxs_rgFaces):
        # Since all faces of original brep are to be removed,
        # delete brep instead of removing faces.
        bSuccess = sc.doc.Objects.Delete(rdBrep_In, True)
        rgBrep_In.Dispose()
        return [] if bSuccess else None
    
    rgBrep_WIP = rgBrep_In.Duplicate()
    
    # Remove faces from brep geometry.
    map(rgBrep_WIP.Faces.RemoveAt, idxs_rgFaces)
    if (rgBrep_In.Faces.Count - rgBrep_WIP.Faces.Count) != len(idxs_rgFaces):
        return None
    return replaceGeometry(rdBrep_In, rgBrep_WIP) # GUIDs of modified (remaining) Breps.


def replaceFaces(rhBrep, idxs_rgFaces, rgBreps_NewGeom, bExtract=False, fTolerance_Join=None, bEcho=False, bDebug=False):
    """Replaces faces using Brep.JoinBreps .
    
    Parameters:
        rhBrep: Guid or BrepObject
        idxs_rgFaces: Indices of faces that are to be replaced.
            ints, iterable of ints, or None (to indicate that all faces are to be replaced.)
        rgBreps_NewGeom: Brep or surface geometry.
        bEcho: Boolean
        bDebug: Boolean
    Returns on success (bExtract==True):
        list(GUIDs of Breps),
        list(GUIDs of remaining Breps)
    Returns on success (bExtract==False):
        list(GUIDs of BrepObjects)
    Returns on fail:
        None
    
    TODO: Allow BrepFace as alternative input for rgBreps_NewGeom.
    """

    rdBrep0 = coerceBrepObject(rhBrep)
    if rdBrep0 is None: return
    rgBrep0 = rdBrep0.BrepGeometry
    
    def coerceList(obj):
        try: obj = list(set(obj))
        except: obj = [obj]
        return obj
    if idxs_rgFaces is None:
        idxs_rgFaces = range(rgBrep0.Faces.Count)
        idxs_rgFaces.reverse()
    else:
        idxs_rgFaces = sorted(coerceList(idxs_rgFaces), reverse=True)
    rgBreps_NewGeom = coerceList(rgBreps_NewGeom)

    bSolidAtStart = rgBrep0.IsSolid

    # Convert any surfaces in list to breps.
    #rgBreps_New1F_WIP = [coerceBrep(o) for o in rgBreps_NewGeom] # Not using due to Dispose() problem.
    breps_NewGeom_WIP = []
    for i in xrange(len(rgBreps_NewGeom)):
        if isinstance(rgBreps_NewGeom[i], rg.Surface):
            srf = rgBreps_NewGeom[i]
            breps_NewGeom_WIP.append(srf.ToBrep())
        else:
            breps_NewGeom_WIP.append(rgBreps_NewGeom[i].Duplicate())
    
    if rgBrep0.Faces.Count == 1:
        # Monoface brep.
        def replaceGeomOfMonofaceBrep():
            gBreps_Out = replaceGeometry(rdBrep0, breps_NewGeom_WIP) # Nore or list of GUIDs.
            rgBrep0.Dispose()
            for b in breps_NewGeom_WIP: b.Dispose()
            return gBreps_Out
        if bExtract:
            return replaceGeomOfMonofaceBrep(), []
        else:
            return replaceGeomOfMonofaceBrep()

    if not bExtract:
        fTol_Join_toUse = max(
            fTolerance_Join,
            1.05 * max([edge.Tolerance for edge in rgBrep0.Edges]),
            2.1 * sc.doc.ModelAbsoluteTolerance,
        ) # None (fTolerance_Join) is acceptable input for max().


    if len(idxs_rgFaces) < rgBrep0.Faces.Count:
        # Not all faces of a polyface brep are to be replaced.
        
        # Remove faces from copy of original brep geometry.
        rgBrep0_WithFacesRemoved = rgBrep0.DuplicateBrep()
        rgBrep0.Dispose()
        for idx_rgFace_ToRemove in idxs_rgFaces:
            rgBrep0_WithFacesRemoved.Faces.RemoveAt(idx_rgFace_ToRemove)

        b, sLog = rgBrep0_WithFacesRemoved.IsValidWithLog()
        if not b:
            if bEcho: print sLog
            rgBrep0_WithFacesRemoved.Dispose()
            for b in breps_NewGeom_WIP: b.Dispose()
            return

        if bExtract:
            def replaceSome_Extract():
                gBreps_Extracted = []
                for brep in breps_NewGeom_WIP:
                    gBrep1_1Face = sc.doc.Objects.AddBrep(brep, rdBrep0.Attributes)
                    brep.Dispose()
                    if gBrep1_1Face != Guid.Empty:
                        gBreps_Extracted.append(gBrep1_1Face)
                if len(gBreps_Extracted) == len(breps_NewGeom_WIP):
                    gBreps_Modified = replaceGeometry(rdBrep0, rgBrep0_WithFacesRemoved)
                    rgBrep0_WithFacesRemoved.Dispose()
                else:
                    if bEcho:
                        s  = "{} new faces could not be added.".format(
                            len(breps_NewGeom_WIP) - len(gBreps_Extracted))
                        s += "  No faces were removed from original brep."
                        print s
                return gBreps_Extracted, gBreps_Modified
            return replaceSome_Extract()
        else:
            def replaceSome_Join():
                gBreps_Out = []

                rgBreps_Joined = rg.Brep.JoinBreps(
                        [rgBrep0_WithFacesRemoved] + breps_NewGeom_WIP,
                        tolerance=fTol_Join_toUse)

                rgBrep0_WithFacesRemoved.Dispose()
                for brep in breps_NewGeom_WIP: brep.Dispose()

                len_Joined = len(rgBreps_Joined)
                if len_Joined == 0:
                    if bEcho:
                        print "rgBreps_Joined is empty.  No objects have been modified."
                    return
                else:
                    if bEcho:
                        print "JoinBreps returned {} breps.".format(len_Joined)
                        if (
                            len_Joined == 1 and
                            bSolidAtStart and
                            not rgBreps_Joined[0].IsSolid
                        ):
                            print "Warning: Brep is no longer a solid!"
                    gBreps_Out = replaceGeometry(rdBrep0, rgBreps_Joined) # None or list of GUIDs.
                    for brep in rgBreps_Joined: brep.Dispose()
                    return gBreps_Out
            return replaceSome_Join()
    else:
        # All faces of a polyface brep are to be replaced.
        rgBrep0.Dispose()

        if bExtract:
            def replaceAll_Extract():
                # All output GUIDs will be new.
                gBreps_Extracted = []
                for brep in breps_NewGeom_WIP:
                    gBrep1_1Face = sc.doc.Objects.AddBrep(brep, rdBrep0.Attributes)
                    brep.Dispose()
                    if gBrep1_1Face != Guid.Empty:
                        gBreps_Extracted.append(gBrep1_1Face)
                if len(gBreps_Extracted) == len(breps_NewGeom_WIP):
                    sc.doc.Objects.Delete(objectId=rdBrep0.Id, quiet=False)
                else:
                    if bEcho:
                        s  = "{} new faces could not be added.".format(
                            len(breps_NewGeom_WIP) - len(gBreps_Extracted))
                        s += "  Original brep was not deleted."
                        print s
                return gBreps_Extracted, []
            return replaceAll_Extract()
        else:
            def replaceAll_Join():
                gBreps_Out = []

                rgBreps_Joined = rg.Brep.JoinBreps(
                        breps_NewGeom_WIP,
                        tolerance=fTol_Join_toUse)

                for brep in breps_NewGeom_WIP: brep.Dispose()

                len_Joined = len(rgBreps_Joined)
                if len_Joined == 0:
                    if bEcho:
                        print "rgBreps_Joined is empty.  No objects have been modified."
                    return
                else:
                    if bEcho:
                        print "JoinBreps returned {} breps.".format(len_Joined)
                        if (
                            len_Joined == 1 and
                            bSolidAtStart and
                            not rgBreps_Joined[0].IsSolid
                        ):
                            print "Warning: Brep is no longer a solid!"
                    gBreps_Out = replaceGeometry(rdBrep0, rgBreps_Joined) # None or list of GUIDs.
                    for brep in rgBreps_Joined: brep.Dispose()
                    return gBreps_Out
            return replaceAll_Join()


def extractFaces(rhBrep, idxFaces, bAddOnlyMonofaces=True, bRetainLayer=True, bRetainColor=True, bCurrentLayer=None, bByLayerColor=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBrep: rd.BrepObject, GUID of rd.BrepObject, or ObjRef of rd.BrepObject
        idxFaces: list(int(Index of face) or not for all.
    Returns:
        list(GUID of extracted brep), list(GUID of remaining brep) on success.
        None on fail.
    
    TODO: Remove support for bCurrentLayer and bByLayerColor after dependent scripts are modified.
    """

    rdBrep_In = coerceBrepObject(rhBrep)
    if rdBrep_In is None: return
    rgBrep_In = rdBrep_In.BrepGeometry
    if not rgBrep_In.IsValid:
        if bEcho: print "Invalid Brep {} will not be processed.".format(rdBrep_In.Id)
        return

    if not idxFaces:
        if rgBrep_In.Faces.Count == 1:
            return [rdBrep_In.Id], []
        idxFaces = range(rgBrep_In.Faces.Count)
    else:
        try: idxFaces = list(set(idxFaces))
        except: idxFaces = [idxFaces]

    if bCurrentLayer is not None: bRetainLayer = not bCurrentLayer
    if bByLayerColor is not None: bRetainColor = not bByLayerColor
    
    if rgBrep_In.Faces.Count == 1:
        # Brep is monoface.
        rgBrep_In.Dispose()
        
        if bRetainLayer and bRetainColor: return [rdBrep_In.Id], []
        
        # The layer and/or color is to be modified.
        attr = rdBrep_In.Attributes
        if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
        rdBrep_In.CommitChanges()
        return [rdBrep_In.Id], []
    
    
    # Brep is polyface.
    
    gBreps1_Extracted = addFromSubsetOfFaces(
            rhBrep=rdBrep_In,
            idxFaces=idxFaces,
            bRetainLayer=bRetainLayer,
            bRetainColor=bRetainColor,
            bAddOnlyMonofaces=bAddOnlyMonofaces,
            bDebug=bDebug
    )
    if gBreps1_Extracted is None: return
    
    rc = removeFaces(rdBrep_In, idxFaces)
    if rc is None:
        return
    else:
        gBs_Remaining = rc
    
    if bEcho:
        print "Extracted {} faces.".format(len(gBreps1_Extracted))
    
    return gBreps1_Extracted, gBs_Remaining


def selectFaces(rhBreps, idxFaces_perBrep=None, surfaceTestFunc=None, bEcho=True):
    if Rhino.RhinoApp.ExeVersion < 6: return

    try: rhBreps = list(rhBreps)
    except: rhBreps = [rhBreps]

    if idxFaces_perBrep:
        if isinstance(idxFaces_perBrep[0], int):
            idxFaces_perBrep = [idxFaces_perBrep]


    if idxFaces_perBrep is None and isinstance(rhBreps[0], rd.ObjRef):
        rc = sortObjrefsIntoBrepGuidsAndFaceIndices(objrefs=rhBreps)
        if rc is None: return

        rhBreps, idxFaces_perBrep = rc

    if idxFaces_perBrep is None and surfaceTestFunc is None:
        print "Parameters result in selection of all faces.  None will be be selected."
        return


    iCt_Fs_Found = 0 # Applies when using surfaceTestFunc.
    iCt_Fs_Selected = 0

    for iB, rhB_In in enumerate(rhBreps):

        idxFaces = idxFaces_perBrep[iB]
    
        rdBrep0 = coerceBrepObject(rhB_In)
        if rdBrep0 is None: return
    
        rgBrep0 = rdBrep0.Geometry
        if not rgBrep0.IsValid:
            rgBrep0.Dispose()
            return

        for idxFace in idxFaces:

            if surfaceTestFunc:
                if not surfaceTestFunc(rgBrep0.Faces[idxFace]):
                    continue

            iCt_Fs_Found += 1

            compIdx = rg.ComponentIndex(
                    rg.ComponentIndexType.BrepFace,
                    index=idxFace)
            iSelected = rdBrep0.SelectSubObject(
                    componentIndex=compIdx,
                    select=True, syncHighlight=True, persistentSelect=True)
            """
            iSelected:
                0: object is not selected
                1: object is selected
                2: object is selected persistently.
            """

        for idxFace in xrange(rgBrep0.Faces.Count):
            compIdx = rg.ComponentIndex(
                    rg.ComponentIndexType.BrepFace, idxFace)
            if rdBrep0.IsSubObjectSelected(compIdx): iCt_Fs_Selected += 1
        rgBrep0.Dispose()


    if bEcho:
        if surfaceTestFunc:
            if iCt_Fs_Found == iCt_Fs_Selected:
                print "{} faces found and are selected.".format(
                    iCt_Fs_Found)
            else:
                print "{} faces found, but only {} are selected.".format(
                    iCt_Fs_Found, iCt_Fs_Selected)
        else:
            print "{} faces are selected.".format(iCt_Fs_Selected)

    return iCt_Fs_Selected


def separateShells(rhBrep, bEcho=True):
    """
    Separate any brep shells of modified brep and act based on shell quantity.
    Returns list of GUIDs of new breps or list with GUID of only existing brep.
    """
    rdBrep0 = coerceBrepObject(rhBrep)
    if rdBrep0 is None: return
    rgBrep0 = rdBrep0.Geometry
    if rgBrep0.Faces.Count == 1:
        return [rdBrep0.Id]
    elif rgBrep0.Faces.Count > 1:
        rgBreps_per_shell = rg.Brep.CreateBooleanUnion(
            [rgBrep0], tolerance=0.0, manifoldOnly=False)
        if rgBreps_per_shell is None:
            if bDebug: print "Error in attempting to separate brep shells.  No objects have been modified."
            return
    else:
        return
    
    if rgBreps_per_shell is None: return []
    elif len(rgBreps_per_shell) == 1: return [rdBrep0.Id]
    else:
        gBreps1 = []
        attr = rdBrep0.Attributes
        for rgBrep1 in rgBreps_per_shell:
            gBrep1 = sc.doc.Objects.AddBrep(rgBrep1, attr)
            if gBrep1 != Guid.Empty:
                gBreps1.append(gBrep1)
            else:
                # Delete all new BrepObjects if any AddBrep fails.
                if bEcho:
                    s  = "Brep could not be added to document."
                    s += "  Brep {} will not have its shells separated.".format(
                            rdBrep0.Id)
                    print s
                for gBrep1 in gBreps1:
                    sc.doc.Objects.Delete(gBrep1, True)
                return
        sc.doc.Objects.Delete(rdBrep0, True)
        return gBreps1


def formatDistance(fDistance):
    if fDistance is None:
        return "(No deviation was provided.)"
    if fDistance == 0.0:
        return "Exactly zero"
    if fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-3)):
        return "{:.1e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def sortObjrefsIntoBrepGuidsAndFaceIndices(objrefs):
    """
    Sorts objrefs of entire breps or faces into 2 lists.

    Returns on success:
        list(GUIDs of breps)
        list(lists(int(Indices of faces)) per brep in other list)
    """

    gBreps = []
    idx_Faces_PerBrep = []
    
    for objref in objrefs:
        rdBrep = objref.Object()
        gBrep = objref.ObjectId
        compIndex = objref.GeometryComponentIndex.Index

        if compIndex == -1:
            rgBrep = objref.Brep()
            if not rgBrep: continue
            if gBrep in gBreps:
                if rgBrep.Faces.Count == len(idx_Faces_PerBrep[gBreps.index(gBrep)]):
                    continue
                else:
                    idx_Faces_PerBrep[gBreps.index(gBrep)] = range(rgBrep.Faces.Count)
            else:
                gBreps.append(gBrep)
                idx_Faces_PerBrep.append(range(rgBrep.Faces.Count))
        else:
            # Face of polyface brep was selected.
            rdObj = objref.Object()
            if rdObj.ObjectType == rd.ObjectType.InstanceReference:
                print "Objects in block instances are not supported."
                continue
            else:
                if gBrep in gBreps:
                    idx_gBrep = gBreps.index(gBrep)
                    idx_Faces_PerBrep[idx_gBrep].append(compIndex)
                else:
                    gBreps.append(gBrep)
                    idx_Faces_PerBrep.append([compIndex])

    return gBreps, idx_Faces_PerBrep


def modifyFacesOfBrepObjects(rhBreps_In, idx_Faces_perBrep=None, surfaceFunc=None, bExtract=False, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBreps_In, # Plural so that summary of 
    """


    if not idx_Faces_perBrep and isinstance(rhBreps_In[0], rd.ObjRef):
        rc = sortObjrefsIntoBrepGuidsAndFaceIndices(objrefs=rhBreps_In)
        if rc is None: return

        rhBreps_In, idx_Faces_perBrep = rc

        #print gBreps_In
        #for idxs in idx_Faces_perBrep: print idxs


    gBs_Out_All = []
    srf_devs_All = []
    sLogs_allBreps = []

    sc.doc.Objects.UnselectAll()


    for iB, rhBrep_In in enumerate(rhBreps_In):

        if isinstance(rhBrep_In, Guid):
            gBrep_In = rhBrep_In
            rdBrep_In = sc.doc.Objects.FindId(gBrep_In) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep_In)
        elif isinstance(rhBrep_In, rd.BrepObject):
            rdBrep_In = rhBrep_In
            gBrep_In = rdBrep_In.Id
        else:
            continue

        rgBrep_In = rdBrep_In.BrepGeometry

        if not idx_Faces_perBrep:
            idx_Faces_toProcess = range(rgBrep_In.Faces.Count)
        else:
            idx_Faces_toProcess = idx_Faces_perBrep[iB]

        s = "brep {} of {}".format(iB+1, len(rhBreps_In))
        Rhino.RhinoApp.SetCommandPrompt(prompt="Processing {} ...".format(s))

        if bDebug: print "Processing brep {}.".format(gBrep_In)


        rgBs_1F_Processed_thisB = []
        idx_Faces_Success_thisB = []
        srf_devs_thisBrep = []

        sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt

        for iF, idx_Face in enumerate(idx_Faces_toProcess):
            Rhino.RhinoApp.SetCommandPrompt(
                    prompt=sCmdPrompt_In +
                    "  Face[{}] of 0-{} ...".\
                            format(iF, len(idx_Faces_toProcess)))
            if bDebug: print Rhino.RhinoApp.CommandPrompt

            rgF_B_In = rgBrep_In.Faces[idx_Face]

            res, sLog = xBrepFace.createModifiedFace(
                rgF_B_In,
                surfaceFunc,
                bDebug=bDebug,
                )

            if res is None:
                sLogs_allBreps.append(sLog)
                continue

            rgB_1F_Created, srf_dev = res
        
            if not rgB_1F_Created.IsValid:
                print "Invalid brep geometry for face (index={}) after modification.".format(idx_Face)
                print "It will not be included in processed results."
                rgB_1F_Created.Dispose()
                continue
        
            rgBs_1F_Processed_thisB.append(rgB_1F_Created)
            idx_Faces_Success_thisB.append(idx_Face)
            srf_devs_thisBrep.append(srf_dev)
            srf_devs_All.append(srf_dev)

        if not rgBs_1F_Processed_thisB:
            continue # to next brep.


        # Modify BrepObject with successful surface conversions.

        rc = replaceFaces(
            rdBrep_In,
            idx_Faces_Success_thisB,
            rgBs_1F_Processed_thisB,
            bExtract=bExtract)
        if rc is None:
            pass
        else:
            if bExtract:
                gBs_Out_thisBrep, gRemaining = rc
            else:
                gBs_Out_thisBrep = rc
            gBs_Out_All.extend(gBs_Out_thisBrep)


    if bEcho:
        for sLog in set(sLogs_allBreps):
            print "[{}] {}".format(sLogs_allBreps.count(sLog), sLog)


    if not gBs_Out_All:
        return []


        if min(srf_devs_All) != max(srf_devs_All):
            print "Maximum deviations range from {} to {}".format(
                        formatDistance(min(srf_devs_All)),
                        formatDistance(max(srf_devs_All)))
        else:
            print "Maximum deviation: {}".format(formatDistance(max(srf_devs_All)))


    return gBs_Out_All


def createBrepObjectsOfNewSurfaces(rhSrfs_In, surfaceFunc, bEcho=True, bDebug=False):
    """
    """

    gBs_Out = []
    srf_devs_All = []
    sLogs = []
    
    sc.doc.Objects.UnselectAll()

    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt

    for iS, rhSrf_In in enumerate(rhSrfs_In):

        Rhino.RhinoApp.SetCommandPrompt(
            "Processing surface {} of {} ...".format(
                iS+1,
                len(rhSrfs_In)))

        if isinstance(rhSrf_In, Guid):
            gBrep_In = rhSrf_In
            rdBrep_In = sc.doc.Objects.FindId(gBrep_In) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep_In)
            rgBrep_In = rdBrep_In.BrepGeometry
            if rgBrep_In.Faces.Count > 1:
                sLogs.append("BrepObject input has more than 1 face and was skipped.")
                continue # to next surface.
            rgFace_In = rgBrep_In.Faces[0]
        elif isinstance(rhSrf_In, rd.BrepObject):
            rdBrep_In = rhSrf_In
            gBrep_In = rdBrep_In.Id
            rgBrep_In = rdBrep_In.BrepGeometry
            if rgBrep_In.Faces.Count > 1:
                sLogs.append("BrepObject input has more than 1 face and was skipped.")
                continue # to next surface.
            rgFace_In = rgBrep_In.Faces[0]
        elif isinstance(rhSrf_In, rd.ObjRef):
            gBrep_In = rhSrf_In.ObjectId
            rdBrep_In = rhSrf_In.Object()
            rgBrep_In = rhSrf_In.Brep()
            rgFace_In = rhSrf_In.Face()
        else:
            sLogs.append("{} passed to createBrepObjects.".format(
                rhSrf_In.GetType().Name))
            continue # to next surface.

        rgSrf_In = rgFace_In.UnderlyingSurface()

        res, sLog = surfaceFunc(rgSrf_In)

        if sLog is not None:
            sLogs.append(sLog)

        if res is None:
            continue # to next surface.

        ns_Res, srf_dev = res

        if not ns_Res:
            continue # to next brep.

        srf_devs_All.append(srf_dev)

        gB_Added = sc.doc.Objects.AddSurface(ns_Res)
        ns_Res.Dispose()
        if gB_Added == Guid.Empty:
            sLogs.append("New surface could not be added to document.")
            continue

        gBs_Out.append(gB_Added)

    for sLog in set(sLogs):
        print "{} count of: {}".format(sLogs.count(sLog), sLog)

    if gBs_Out:
        print "{} surfaces added to document.".format(len(gBs_Out))
        print "Maximum deviation from original surface: {}".format(formatDistance(max(srf_devs_All)))
    else:
        print "No surfaces added to document."

    return gBs_Out


