"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
181628: Created.
191025: Now accepts invalid breps.
210319: Now also reports surface count.
250130: Now, input of Nothing will select all normal breps.
        Refactored.
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


def main():

    res, objrefs = Rhino.Input.RhinoGet.GetMultipleObjects(
        "Select breps <All normal>",
        acceptNothing=True,
        filter=rd.ObjectType.Brep)

    if res == Rhino.Commands.Result.Cancel: return

    rdBs = [_.Object() for _ in objrefs] if objrefs else _getAllNormalBreps()

    iCt_Fs_All = iCt_Ss_All = iCt_Bs = 0

    iCt_Bs_Input = len(rdBs)
    
    bReportEveryBrep = False
    if iCt_Bs_Input > 10:
        s  = "More than 10 breps are selected, so only total faces of all,"
        s += " not individual brep counts, is reported."
        print(s)
    elif iCt_Bs_Input > 1 and iCt_Bs_Input <= 10:
        bReportEveryBrep = True
    
    for rdB in rdBs:
        rgB = rdB.BrepGeometry
        
        iCt_Fs = rgB.Faces.Count
        iCt_Ss = rgB.Surfaces.Count
        
        if iCt_Fs is None:
            s  = "Faces were not obtained from {}.".format(rdB.Id)
            s += "  Check validity of brep."
            print(s)
            continue
        
        iCt_Fs_All += iCt_Fs
        iCt_Ss_All += iCt_Ss
        iCt_Bs += 1
        
        if bReportEveryBrep:
            s = "{} faces in {}, ".format(
                    iCt_Fs, rdB.Id)
            if iCt_Ss != iCt_Fs:
                s += "but {} surfaces are in brep.".format(iCt_Ss)
            else:
                s += "and surface count is the same."

            print(s)
    
    if iCt_Fs_All > 0:
        s = "{} faces in {}, ".format(
            iCt_Fs_All,
            ("{} breps".format(iCt_Bs) if bool(iCt_Bs-1) else "brep"))
        if iCt_Ss_All != iCt_Fs_All:
            s += "but {} surfaces are in brep.".format(iCt_Ss_All)
        else:
            s += "and surface count is the same."

        print(s)


if __name__ == '__main__': main()