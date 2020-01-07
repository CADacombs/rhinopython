"""
181025: Created, starting with another script.
181108: Added bAddOnlyMonofaces to addFromSubsetOfFaces.
        Added extractFaces.
181202: Added comments.
190130: extractFaces now modifies attributes (as desired) of monoface breps.
190412: Added iterable checker to repslaceFaces.
190428-30: Changed some parameter names.
190506: Added a function.
190519: Fixed some print bugs.
190524: Added a function.
190530: replaceFaces now allows replacement of all faces.
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
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Diagnostics import Stopwatch


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


def addFromSubsetOfFaces(rhBrep0, idxFaces, bAddOnlyMonofaces=True, bRetainLayer=True, bRetainColor=True, bDebug=False):
    """
    """
    
    rdBrep0 = coerceBrepObject(rhBrep0)
    if rdBrep0 is None: return
    rgBrep0 = rdBrep0.BrepGeometry
    if not rgBrep0.IsValid: return
    gBrep0 = rs.coerceguid(rhBrep0)
    
    
    def addBrepOfSubsetOfFaces_JoinBreps():
        
        # Duplicate faces to their own breps to be joined.
        rgBreps1 = [] # Faces (breps) to be duplicated.
        
        for i in idxFaces:
            rgFace = rgBrep0.Faces[i]
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
        
        
        attr = rdBrep0.Attributes.Duplicate()
        if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
        
        
        gBreps_1Shell = map(lambda x: sc.doc.Objects.AddBrep(x, attr), rgBreps_per_shell)
        map(lambda x: x.Dispose(), rgBreps_per_shell)
        
        return gBreps_1Shell
    
    
    def addBrepOfSubsetOfFaces_RemoveAt():
        
        rgBrep1 = rgBrep0.Duplicate()
        
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
        
        attr = rdBrep0.Attributes.Duplicate()
        if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
        
        stopwatch.Restart()
        gBreps_1Shell = map(lambda x: sc.doc.Objects.AddBrep(x, attr), rgBreps_per_shell)
        stopwatch.Stop()
        if bDebug: print "{:.1f} seconds for AddBrep".format(stopwatch.Elapsed.TotalSeconds)
        map(lambda x: x.Dispose(), rgBreps_per_shell)
        
        return gBreps_1Shell
    
    
    nFaces = rgBrep0.Faces.Count
    
    # If brep has only 1 face, return the brep's GUID.
    if nFaces == 1:
        return [gBrep0]
    
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
    
    attr = rdBrep0.Attributes.Duplicate()
    if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
    
    gBreps1_Extracted = []
    
    for idx in idxFaces:
        rgFace = rgBrep0.Faces[idx]
        
        # Duplicate face to its own brep.
        rgBrep1 = rgFace.DuplicateFace(duplicateMeshes=True)
        gBrep1 = sc.doc.Objects.AddBrep(rgBrep1, attr)
        if gBrep1 is None:
            print "Brep face {} from {} could not be added to document.".format(
                    idx, gBrep0)
            rc = rs.MessageBox("Continue extracting faces, skipping this one?",
                    buttons=4,
                    title="addFromSubsetOfFaces")
            if rc is not None and rc == 6:
                continue
            #if not bDebug: rs.DeleteObjects(gBreps1_Extracted)
            return
        
        gBreps1_Extracted.append(gBrep1)
        rgBrep1.Dispose()
    
    return gBreps1_Extracted


def replaceGeometry(brep_toReplace, rgBreps_NewGeom):
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
    
    rdBrep_In = coerceBrepObject(brep_toReplace)

    try: rgBreps_NewGeom = list(set(rgBreps_NewGeom))
    except: rgBreps_NewGeom = [rgBreps_NewGeom]

    rgBreps_NewGeom = [coerceBrep(o) for o in rgBreps_NewGeom]

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


def removeFaces(rhBrep0, idxs_rgFaces):
    """
    rhBrep0: GUID or BrepObject.
    idxs_rgFaces: Iterable of face indices.
    
    Returns:
        Iterable of GUIDs of breps with faces removed.
        None on fail.
    """

    try: idxs_rgFaces = sorted(list(set(idxs_rgFaces)), reverse=True)
    except: idxs_rgFaces = [idxs_rgFaces]

    rdBrep0 = coerceBrepObject(rhBrep0)
    if rdBrep0 is None: return
    rgBrep0 = rdBrep0.Geometry
    if not rgBrep0.IsValid and rgBrep0.Faces.Count > 1: return
    gBrep0 = rs.coerceguid(rhBrep0)
    
    if rgBrep0.Faces.Count == len(idxs_rgFaces):
        # Since all faces of original brep are affected, delete brep instead of removing faces.
        bSuccess = sc.doc.Objects.Delete(rdBrep0, True)
        rgBrep0.Dispose()
        return [] if bSuccess else None
    
    rgBrep_WIP = rgBrep0.Duplicate()
    
    # Remove faces from brep geometry.
    map(rgBrep_WIP.Faces.RemoveAt, idxs_rgFaces)
    if (rgBrep0.Faces.Count - rgBrep_WIP.Faces.Count) != len(idxs_rgFaces):
        return None
    return replaceGeometry(rdBrep0, rgBrep_WIP) # GUIDs of modified (remaining) Breps.


def replaceFaces(rhBrep0, idxs_rgFaces, rgBreps_NewGeom, bExtract=False, fTolerance_Join=None, bEcho=True, bDebug=False):
    """Replaces faces using Brep.JoinBreps .
    
    Parameters:
        rhBrep0: Guid or BrepObject
        idxs_rgFaces: Indices of faces that are to be replaced.
            ints, iterable of ints, or None (to indicate that all faces are to be replaced.)
        rgBreps_NewGeom: Brep or surface geometry.
        bEcho: Boolean
        bDebug: Boolean
    Returns:
        bExtract==True
            Tuple: (List of GUIDs of Breps), (List of GUIDs of remaining Breps)
        bExtract==False
            List of GUIDs of BrepObjects.
        None on fail.
    
    TODO: Allow BrepFace as alternative input for rgBreps_NewGeom.
    """

    rdBrep0 = coerceBrepObject(rhBrep0)
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
            print sLog
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
                    print "rgBreps_Joined is empty.  No objects have been modified."
                    return
                else:
                    print "JoinBreps returned {} breps.".format(len_Joined)
                    if len_Joined == 1 and bSolidAtStart and not rgBreps_Joined[0].IsSolid:
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
                    print "rgBreps_Joined is empty.  No objects have been modified."
                    return
                else:
                    print "JoinBreps returned {} breps.".format(len_Joined)
                    if len_Joined == 1 and bSolidAtStart and not rgBreps_Joined[0].IsSolid:
                        print "Warning: Brep is no longer a solid!"
                    gBreps_Out = replaceGeometry(rdBrep0, rgBreps_Joined) # None or list of GUIDs.
                    for brep in rgBreps_Joined: brep.Dispose()
                    return gBreps_Out
            return replaceAll_Join()


def extractFaces(rhBrep0, idxFaces, bAddOnlyMonofaces=True, bRetainLayer=True, bRetainColor=True, bCurrentLayer=None, bByLayerColor=None, bEcho=True, bDebug=False):
    """
    Returns gBreps1_Extracted.
    
    TODO: Remove support for bCurrentLayer and bByLayerColor after dependent scripts are modified.
    """

    try: idxFaces = list(set(idxFaces))
    except: idxFaces = [idxFaces]

    rdBrep0 = coerceBrepObject(rhBrep0)
    if rdBrep0 is None: return
    rgBrep0 = rdBrep0.Geometry
    if not rgBrep0.IsValid:
        if bEcho: print "Invalid Brep {} will not be processed.".format(rdBrep0.Id)
        return
    gBrep0 = rs.coerceguid(rhBrep0)
    
    if bCurrentLayer is not None: bRetainLayer = not bCurrentLayer
    if bByLayerColor is not None: bRetainColor = not bByLayerColor
    
    if rgBrep0.Faces.Count == 1:
        # Brep is monoface.
        rgBrep0.Dispose()
        
        if bRetainLayer and bRetainColor: return [gBrep0], []
        
        # The layer and/or color is to be modified.
        attr = rdBrep0.Attributes
        if not bRetainLayer: attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        if not bRetainColor: attr.ColorSource = rd.ObjectColorSource.ColorFromLayer
        rdBrep0.CommitChanges()
        return [gBrep0], []
    
    
    # Brep is polyface.
    
    gBreps1_Extracted = addFromSubsetOfFaces(
            rhBrep0=rdBrep0,
            idxFaces=idxFaces,
            bRetainLayer=bRetainLayer,
            bRetainColor=bRetainColor,
            bAddOnlyMonofaces=bAddOnlyMonofaces,
            bDebug=bDebug
    )
    if gBreps1_Extracted is None: return
    
    rc = removeFaces(rdBrep0, idxFaces)
    if rc is None:
        return
    else:
        gBs_Remaining = rc
    
    if bEcho:
        print "Extracted {} faces.".format(len(gBreps1_Extracted))
    
    return gBreps1_Extracted, gBs_Remaining


def selectFaces(rhBrep0, idxFaces, bEcho=True):
    if Rhino.RhinoApp.ExeVersion < 6: return
    
    rdBrep0 = coerceBrepObject(rhBrep0)
    if rdBrep0 is None: return
    
    rgBrep0 = rdBrep0.Geometry
    if not rgBrep0.IsValid:
        rgBrep0.Dispose()
        return
    
    for idx_rgFace in idxFaces:
        compIdx = rg.ComponentIndex(
                rg.ComponentIndexType.BrepFace,
                index=idx_rgFace)
        rdBrep0.SelectSubObject(
                componentIndex=compIdx,
                select=True, syncHighlight=True, persistentSelect=True)
    ct_Faces_Selected = 0
    for idxFace in xrange(rgBrep0.Faces.Count):
        compIdx = rg.ComponentIndex(
                rg.ComponentIndexType.BrepFace, idxFace)
        if rdBrep0.IsSubObjectSelected(compIdx): ct_Faces_Selected += 1
    rgBrep0.Dispose()
    if len(idxFaces) == ct_Faces_Selected:
        if bEcho:
            print "{} faces found and are selected.".format(
                    len(idxFaces))
    else:
        if bEcho:
            s = "{} faces found, but only {} are selected.".format(
                    len(idxFaces), ct_Faces_Selected)
            print s
    
    return ct_Faces_Selected


def separateShells(rhBrep0, bEcho=True):
    """
    Separate any brep shells of modified brep and act based on shell quantity.
    Returns list of GUIDs of new breps or list with GUID of only existing brep.
    """
    rdBrep0 = coerceBrepObject(rhBrep0)
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
