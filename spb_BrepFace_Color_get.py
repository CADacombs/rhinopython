"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
200617: Created.
230331: Now allows multiple face input and processing.
250811: Now, accesses BrepFaces via BrepObjects instead of ObjRefs because the latter
        was too slow and memory intensive for some models.
        Now when only 1 face is selected, object color is also reported.
"""

import Rhino
import Rhino.DocObjects as rd


def _sortBrepObjects_and_face_indices(objrefs):
    rdBs = []
    gBs = []
    iFs_perB = []

    for objref in objrefs:
        if objref.Object().ObjectType != rd.ObjectType.Brep: continue

        if objref.GeometryComponentIndex.Index == -1:
            if objref.ObjectId in gBs:
                print("Already selected. Why is it double-selected?")
                pass
            else:
                rdBs.append(objref.Object())
                gBs.append(objref.ObjectId)
                iFs_perB.append([0])
        elif (
            objref.GeometryComponentIndex.ComponentIndexType ==
            Rhino.Geometry.ComponentIndexType.BrepFace
            ):
            if objref.ObjectId in gBs:
                iB = gBs.index(objref.ObjectId)
                if iFs_perB[iB] is None:
                    print("Already selected. Why is it double-selected?")
                    pass
                else:
                    iFs_perB[iB].append(objref.GeometryComponentIndex.Index)
            else:
                rdBs.append(objref.Object())
                gBs.append(objref.ObjectId)
                iFs_perB.append([objref.GeometryComponentIndex.Index])


    iFs_perB_Sorted = []
    for iFs in iFs_perB:
        if iFs is None:
            iFs_perB_Sorted.append(None)
        else:
            iFs_perB_Sorted.append(sorted(iFs))

    return rdBs, iFs_perB_Sorted


def _formatColor(color):
    return str(color).strip('Color [').strip(']')

def main():

    go = Rhino.Input.Custom.GetObject()
    go.SetCommandPrompt("Select faces")
    go.GeometryFilter = rd.ObjectType.Surface

    res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

    if res == Rhino.Input.GetResult.Cancel:
        go.Dispose()
        return

    objrefs = go.Objects()
    go.Dispose()

    rdBreps, idxFaces_perBrep = _sortBrepObjects_and_face_indices(objrefs)

    colors = []

    # The more faces that are selected, the slower and more memory intensive
    # this routine is than the succeeding routine.
    #for objref in objrefs:
    #    colors.append(objref.Face().PerFaceColor)

    for rdBrep, iFs in zip(rdBreps, idxFaces_perBrep):
        rgBrep = rdBrep.Geometry
        #print(rdBrep, iFs[:10])
        for iF in iFs:
            colors.append(rgBrep.Faces[iF].PerFaceColor)

    if len(objrefs) == 1:
        attr = rdBreps[0].Attributes
        colorSource = attr.ColorSource
        objectcolor = attr.ObjectColor
        if objectcolor == colors[0]:
            print("PerFaceColor & BrepObject color: {}".format(_formatColor(colors[0])))
        else:
            print("PerFaceColor:     {}".format(_formatColor(colors[0])))
            print("BrepObject color: {}".format(
                colorSource if colorSource is rd.ObjectColorSource.ColorFromObject
                else _formatColor(objectcolor)))
    else:
        for color in set(colors):
            print("{} of PerFaceColor {}".format(
                colors.count(color),
                _formatColor(color)))


if __name__ == '__main__': main()