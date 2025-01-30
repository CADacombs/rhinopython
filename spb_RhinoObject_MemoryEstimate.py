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

    res, objrefs = Rhino.Input.RhinoGet.GetMultipleObjects(
        "Select objects <All normal>",
        acceptNothing=True,
        filter=rd.ObjectType.AnyObject)

    if res == Rhino.Commands.Result.Cancel: return

    if objrefs:
        rdObjs = [_.Object() for _ in objrefs]
    else:
        rdObjs = _getAllNormalObjects()

    iMB_Divisor = 1024**2

    uint_MEs_Total = 0

    sDef_Names = []

    for rdO in rdObjs:
        if isinstance(rdO, rd.InstanceObject):
            # Size of each block definition is only obtained and added once.
            rdInstDef = rdO.InstanceDefinition
            if rdInstDef.Name in sDef_Names:
                continue
            else:
                sDef_Names.append(rdInstDef.Name)

        uint_MEs_ThisObj = _getMemoryEstimate(rdO)

        if len(rdObjs) < 11:
            print("{:.3f} MB in {}".format(
                uint_MEs_ThisObj/iMB_Divisor,
                rdO.GetType().Name))

        uint_MEs_Total += uint_MEs_ThisObj

    if len(rdObjs) > 1:
        print("Total: {:.3f} MB in {} objects.".format(
            uint_MEs_Total/iMB_Divisor,
            len(rdObjs)))


if __name__ == '__main__': main()