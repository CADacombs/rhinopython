"""
All memory sizes are output in MB.
Size of each block definition is only obtained and added once.
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


def _getAllNormalObjects():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    return list(sc.doc.Objects.GetObjectList(oes))


def _getMemoryEstimate(rhObj):
    if isinstance(rhObj, rd.InstanceObject):
        iME = 0
        rdInstDef = rhObj.InstanceDefinition
        for rgO in rdInstDef.GetObjects():
            uint_MEs_ThisObj = _getMemoryEstimate(rgO)
            iME += uint_MEs_ThisObj
        return iME
    else:
        return rhObj.MemoryEstimate()


def main():

    res, iCt_toFind = Rhino.Input.RhinoGet.GetInteger(
        prompt="Number of objects with greatest memory estimates to select",
        acceptNothing=True,
        outputNumber=1,
        lowerLimit=1,
        upperLimit=-Rhino.RhinoMath.UnsetIntIndex)

    if res == Rhino.Commands.Result.Cancel: return

    rdObjs = _getAllNormalObjects()

    iMB_Divisor = 1024**2

    uint_MEs_Total = 0
    rdOs_Recorded = []
    sDef_Names = []
    iMBs = []

    for rdO in rdObjs:
        if isinstance(rdO, rd.InstanceObject):
            # Size of each block definition is only obtained and added once.
            rdInstDef = rdO.InstanceDefinition
            if rdInstDef.Name in sDef_Names:
                continue
            else:
                sDef_Names.append(rdInstDef.Name)

        uint_MEs_ThisObj = _getMemoryEstimate(rdO)

        rdOs_Recorded.append(rdO)
        iMBs.append(uint_MEs_ThisObj)

        uint_MEs_Total += uint_MEs_ThisObj

    zipped_lists = zip(iMBs, rdOs_Recorded)

    sorted_zipped_lists = sorted(zipped_lists, reverse=True)

    sorted_iMBs, sorted_rdOs_Recorded = zip(*sorted_zipped_lists)

    uints_MBs_forFound = sorted_iMBs[:iCt_toFind]
    rdOs_Found = sorted_rdOs_Recorded[:iCt_toFind]

    print("{:.3f} MB in {} objects.".format(
        sum(uints_MBs_forFound)/iMB_Divisor,
        len(rdOs_Found)))

    if sc.doc.Objects.Select(objectIds=[_.Id for _ in rdOs_Found]):
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()