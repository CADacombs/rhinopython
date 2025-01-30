"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
200603: Created.
250130: Now, input of Nothing will select all normal breps. Added information to printed report.
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

    return rgFs, rgBs


def main():

    res, objrefs = Rhino.Input.RhinoGet.GetMultipleObjects(
        "Select breps and/or faces <All normal breps>",
        acceptNothing=True,
        filter=(rd.ObjectType.Brep | rd.ObjectType.Surface))

    if res == Rhino.Commands.Result.Cancel: return

    if objrefs:
        rgFs, rgBs = _get_Faces_and_Breps_from_ObjRefs(objrefs)
    else:
        rgBs = [_.BrepGeometry for _ in _getAllNormalBreps()]
        rgFs = []

    iCtTs_All = 0
    iCt_TrimTypes = []
    iCt_IsoStatuses = []
    edgeValences = []

    for rgB in rgBs:
        iCtTs_All += rgB.Trims.Count

        for rgT in rgB.Trims:
            iCt_TrimTypes.append(rgT.TrimType)
            iCt_IsoStatuses.append(rgT.IsoStatus)
            rgE = rgT.Edge
            if rgE is not None:
                edgeValences.append(rgE.Valence)

    for rgF in rgFs:
        iCtTs_All += sum([rgL.Trims.Count for rgL in rgF.Loops])

        for rgL in rgF.Loops:
            for rgT in rgL.Trims:
                iCt_TrimTypes.append(rgT.TrimType)
                iCt_IsoStatuses.append(rgT.IsoStatus)
                rgE = rgT.Edge
                if rgE is not None:
                    edgeValences.append(rgE.Valence)


    print("{} trims in {}:".format(
        iCtTs_All,
        ", ".join([
            "{} breps".format(len(rgBs)),
            "{} selected faces".format(len(rgFs)),
            ]
            )
        ))

    print("  TrimTypes: {}".format(
        ', '.join(["{} {}".format(iCt_TrimTypes.count(_), _) for _ in set(iCt_TrimTypes)])))

    print("  IsoStatuses: {}".format(
        ', '.join(["{} {}".format(iCt_IsoStatuses.count(_), _) for _ in set(iCt_IsoStatuses)])))

    print("  Edge valences: {}".format(
        ", ".join(["{} {}".format(edgeValences.count(_), _) for _ in set(edgeValences)])))


if __name__ == '__main__': main()