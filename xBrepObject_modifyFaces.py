"""
Library module.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220317: Created as a breakaway from another script.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Diagnostics import Stopwatch

import xBrepFace
import xBrepObject


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
                print("Objects in block instances are not supported.")
                continue
            else:
                if gBrep in gBreps:
                    idx_gBrep = gBreps.index(gBrep)
                    idx_Faces_PerBrep[idx_gBrep].append(compIndex)
                else:
                    gBreps.append(gBrep)
                    idx_Faces_PerBrep.append([compIndex])

    return gBreps, idx_Faces_PerBrep


def modifyFaces(rhBreps_In, idx_Faces_perBrep=None, surfaceFunc=None, bExtract=False, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBreps_In, # Plural so that summary of 
    """


    if not idx_Faces_perBrep and isinstance(rhBreps_In[0], rd.ObjRef):
        rc = sortObjrefsIntoBrepGuidsAndFaceIndices(objrefs=rhBreps_In)
        if rc is None: return

        rhBreps_In, idx_Faces_perBrep = rc

        #print(gBreps_In)
        #for idxs in idx_Faces_perBrep: print(idxs)


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

        if bDebug: print("Processing brep {}.".format(gBrep_In))


        rgBs_1F_Processed_thisB = []
        idx_Faces_Success_thisB = []
        srf_devs_thisBrep = []

        sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt

        for iF, idx_Face in enumerate(idx_Faces_toProcess):
            Rhino.RhinoApp.SetCommandPrompt(
                    prompt=sCmdPrompt_In +
                    "  Face[{}] of 0-{} ...".\
                            format(iF, len(idx_Faces_toProcess)))
            if bDebug: print(Rhino.RhinoApp.CommandPrompt)

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
                print("Invalid brep geometry for face (index={}) after modification.".format(idx_Face))
                print("It will not be included in processed results.")
                rgB_1F_Created.Dispose()
                continue
        
            rgBs_1F_Processed_thisB.append(rgB_1F_Created)
            idx_Faces_Success_thisB.append(idx_Face)
            srf_devs_thisBrep.append(srf_dev)
            srf_devs_All.append(srf_dev)

        if not rgBs_1F_Processed_thisB:
            continue # to next brep.


        # Modify BrepObject with successful surface conversions.

        rc = xBrepObject.replaceFaces(
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
            print("[{}] {}".format(sLogs_allBreps.count(sLog), sLog))


    if not gBs_Out_All:
        return []


        if min(srf_devs_All) != max(srf_devs_All):
            print("Maximum deviations range from {} to {}".format(
                    formatDistance(min(srf_devs_All)),
                    formatDistance(max(srf_devs_All)))
                    )
        else:
            print("Maximum deviation: {}".format(formatDistance(max(srf_devs_All))))


    return gBs_Out_All
