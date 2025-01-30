"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
171107: Created.
...
171209-11: Added a function. Modularized.
171220: Bug fix in final printed result.
180625: Now only reports edge counts per brep for brep count up to 10.
180726: Changed method in a function.
181031: Corrected behavior for object preselection.
190829: Replaced GetObjects from rhinoscriptsyntax with RhinoGet.GetMultipleObjects due sometimes getting an unsolved bug.
220401: Now accepts subselection of faces.
241230: Bug fix in counting naked edges of subselected faces.
250130: Now, input of Nothing will select all normal breps. Added information to printed report.
        Refactored to better support mix of face and brep selections.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc


def _getAllNormalBreps():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))


def _get_Faces_and_Breps_from_ObjRefs(objrefs):
    rgFs = []
    gBs_perF = []
    rgBs = []
    gBs_perB = []
    for objref in objrefs:
        rgF = objref.Face()
        if rgF is not None:
            rgFs.append(rgF)
            gBs_perF.append(objref.ObjectId)
        else:
            rgBs.append(objref.Brep())
            gBs_perB.append(objref.ObjectId)

    # Remove faces whose entire breps are to be analyzed.
    for i in range(len(gBs_perF)-1, -1, -1):
        if gBs_perF[i] in gBs_perB:
            del rgFs[i]
            del gBs_perF[i]

    return rgFs, gBs_perF, rgBs, gBs_perB


def main():

    res, objrefs = Rhino.Input.RhinoGet.GetMultipleObjects(
        "Select breps and/or faces <All normal breps>",
        acceptNothing=True,
        filter=(rd.ObjectType.Brep | rd.ObjectType.Surface))

    if res == Rhino.Commands.Result.Cancel: return

    if objrefs:
        rgFs, gBs_perF, rgBs, gBs_perB = _get_Faces_and_Breps_from_ObjRefs(objrefs)
    else:
        rdBs = _getAllNormalBreps()
        rgBs = [_.BrepGeometry for _ in rdBs]
        gBs_perB = [_.Id for _ in rdBs]
        rgFs = []
        gBs_perF = []

    iCt_Es_All = iCt_Es_ThisObj = 0
    edgeValences_All = []
    iCt_NEs_All = 0

    iCt_Objs_Input = len(rgFs) + len(rgBs)
    iCt_Objs_Valid = 0

    bReportEveryObj = iCt_Objs_Input < 11

    for rgO, gB in zip(rgFs + rgBs, gBs_perF + gBs_perB):
        if isinstance(rgO, rg.BrepFace):
            rgF = rgO
            rgB = rgF.Brep
        else:
            rgB = rgO
            rgF = None

        if not rgB.IsValid:
            print("Warning. Brep {} is not valid.".format(gB),
                  "Repair before trusting edge count results.")

        if rgF:
            idx_Es = rgF.AdjacentEdges()
            iCt_Es_ThisObj = len(rgF.AdjacentEdges())
        else:
            iCt_Es_ThisObj = rgB.Edges.Count

        if iCt_Es_ThisObj is None: # Don't check for 'not' or 0 in case the brep doesn't have edges.
            s  = "Edges were not obtained from {}.".format(gB)
            s += " Check validity of brep."
            print(s)
            continue

        iCt_Objs_Valid += 1
        iCt_Es_All += iCt_Es_ThisObj

        edgeValences_ThisObj = []

        if rgF:
            edgeValences_ThisObj = [rgB.Edges[idx_E].Valence for idx_E in idx_Es]
        else:
            edgeValences_ThisObj = [rgE.Valence for rgE in rgB.Edges]

        edgeValences_All.extend(edgeValences_ThisObj)

        if bReportEveryObj:
            if rgF:
                print("{} edges ({}) in face of {}".format(
                    iCt_Es_ThisObj,
                    " + ".join(["{} {}".format(edgeValences_ThisObj.count(_), _) for _ in set(edgeValences_ThisObj)]),
                    gB))
            else:
                print("{} edges ({}) in {}".format(
                    iCt_Es_ThisObj,
                    " + ".join(["{} {}".format(edgeValences_ThisObj.count(_), _) for _ in set(edgeValences_ThisObj)]),
                    gB))

    if (iCt_Objs_Input > 1) and iCt_Es_All:
        print("{} edges ({}) in {} breps and/or faces.".format(
            iCt_Es_All,
            " + ".join(["{} {}".format(edgeValences_All.count(_), _) for _ in set(edgeValences_All)]),
            iCt_Objs_Valid))


if __name__ == '__main__': main()