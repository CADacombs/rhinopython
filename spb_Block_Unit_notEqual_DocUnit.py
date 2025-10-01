"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250923: Created.
"""

import Rhino
import rhinoscriptsyntax as rs
import scriptcontext as sc


def main():

    rdIdefs_All = sc.doc.InstanceDefinitions.GetList(ignoreDeleted=True)

    if not rdIdefs_All:
        print("No blocks in document.")
        return
    
    rdIdefs_Found = []
    rdIdefs_Linked = []
    
    for rdIdef in rdIdefs_All:
        if rdIdef.UpdateType == Rhino.DocObjects.InstanceDefinitionUpdateType.Linked:
            rdIdefs_Linked.append(rdIdef)
            continue
        if rdIdef.UnitSystem != sc.doc.ModelUnitSystem:
            rdIdefs_Found.append(rdIdef)
    
    if rdIdefs_Linked:
        print("Skipped {} linked block definitions.".format(len(rdIdefs_Linked)))
    
    if not rdIdefs_Found:
        print("All block definitions in document have the same unit system as the document.")
        return
    else:
        print("Block definitions whose unit systems are not {}, the document's unit system:".format(
            sc.doc.ModelUnitSystem))
        for idIdef in rdIdefs_Found:
            print("{}: {}".format(rdIdef.Name, rdIdef.UnitSystem))

    sChoice = rs.GetString(
        message="Set all units to {}?".format(sc.doc.ModelUnitSystem),
        defaultString="No",
        strings=['Yes', 'No']
        )
    if sChoice != 'Yes': return

    rdIdefs_SetUnit = []
    rdIdefs_SetUnitFails = []
    
    for rdIdef in rdIdefs_Found:
        rdIdef.UnitSystem = sc.doc.ModelUnitSystem
        if rdIdef.UnitSystem == sc.doc.ModelUnitSystem:
            rdIdefs_SetUnit.append(rdIdef)
        else:
            rdIdefs_SetUnitFails.append(rdIdef)
            continue
        
        
    if not rdIdefs_SetUnitFails:
        print("All {} block definitions' unit systems have been set to {}.".format(
            len(rdIdefs_SetUnit), sc.doc.ModelUnitSystem))
    
    
    #if bDebug: sEval="brep"; print("{}: {}".format(sEval, eval(sEval)))

if __name__ == '__main__': main()