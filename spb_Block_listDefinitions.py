"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180528: Created, replacing ListDefinitionNamesOfSelBlockInstances.rvb.
180606: Testing alternative output format.
180620: Added rs.BlockInstanceCount to printed summary.
        Changed name from definitionNamesOfSelectedBlockInstances.py to listDefinitionNamesOfSelectedBlockInstances.py
181024: Added check at start for whether there are any blocks in the document.
250421: Now reports more useful output for nested counts.
        Now, can also report on all of document's blocks instead of a selection.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from clr import StrongBox


def getInput():
    
    # Create options with default values before any are changed with those in sticky.
    optT_bCount = ri.Custom.OptionToggle(True, 'No', 'Yes')
    
    # Load sticky.
    stickyKeys = (
            'bCount({})'.format(__file__),
        )
    stickyValues = [
        optT_bCount.CurrentValue,
        ] # Default values.
    for i, stickyKey in enumerate(stickyKeys):
        if sc.sticky.has_key(stickyKey): stickyValues[i] = sc.sticky[stickyKey]
    (
        optT_bCount.CurrentValue,
        ) = stickyValues
    
    # Get instance references with optional input.
    
    def addOptions():
        go.AddOptionToggle('IncludeQuantity', optT_bCount)
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select block instances")
    
    go.GeometryFilter = rd.ObjectType.InstanceReference

    go.AcceptNothing(True)
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    while True:
        go.ClearCommandOptions()
        addOptions()

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return [], optT_bCount.CurrentValue

        if res == ri.GetResult.Object:
            gInsts = [objref.ObjectId for objref in go.Objects()]
            go.Dispose()

            return (
                gInsts,
                optT_bCount.CurrentValue,
                )

        stickyValues = (optT_bCount.CurrentValue,)
        for i, stickyKey in enumerate(stickyKeys):
            sc.sticky[stickyKey] = stickyValues[i]


def main():
    rdIdefs = sc.doc.InstanceDefinitions.GetList(ignoreDeleted=True)

    if not rdIdefs:
        print("No blocks in document.")
        return
    
    rc = getInput()
    if rc is None: return

    gInsts, bCount = rc

    sBlockNames = []
    iSelInstCts = []
    rdIdefs_Out = []

    if len(gInsts) == 0:
        sBlockNames = [_.Name for _ in rdIdefs]
        iSelInstCts = [0]*len(rdIdefs)
        rdIdefs_Out = rdIdefs
    else:
        for gInst in gInsts:
            sDefName = rs.BlockInstanceName(gInst)
            if sDefName not in sBlockNames:
                sBlockNames.append(sDefName)
                iSelInstCts.append(1)
                rdIdefs_Out.append(
                    sc.doc.InstanceDefinitions.Find(
                        instanceDefinitionName=sDefName,
                        ignoreDeletedInstanceDefinitions=True))
            else:
                iSelInstCts[sBlockNames.index(sDefName)] += 1

    if not bCount:
        for name in sorted(sBlockNames):
            print("{}".format(name))
        return

    try:
        rd.InstanceDefinition.UseCount
        bRhinoVerIsAtLeast7_4 = True
    except:
        bRhinoVerIsAtLeast7_4 = False

    if bRhinoVerIsAtLeast7_4:
        topLevelReferenceCount = StrongBox[int]()
        nestedReferenceCount = StrongBox[int]()

        for name, rdIdef, qty in sorted(zip(sBlockNames, rdIdefs_Out, iSelInstCts)):
            iUseCount_Total = rdIdef.UseCount(topLevelReferenceCount, nestedReferenceCount)
            s = ""
            if qty:
                s += "{} selected of ".format(qty)
            s += "{} ({} documentTotal = {} top-level + {} nested".format(
                name,
                iUseCount_Total,
                topLevelReferenceCount,
                nestedReferenceCount,
                )
            if nestedReferenceCount.Value:
                rdIrefs = rdIdef.GetReferences(wheretoLook=2)
                sDefNames = [_.InstanceDefinition.Name for _ in rdIrefs]
                s += " in {} definitions".format(len(rdIdef.GetContainers()))
            s += ")"
            print(s)
    else:
        if iSelInstCts:
            print("(Selected count) of (Name) ( (top level refs in active doc), (top level and nested refs in active doc), (refs from other inst defs) )")
            for name, rdIdef, qty in sorted(zip(sBlockNames, rdIdefs_Out, iSelInstCts)):
                print("{} of {} ({}, {}, {})".format(
                    qty,
                    name,
                    rs.BlockInstanceCount(block_name=name, where_to_look=0),
                    rs.BlockInstanceCount(block_name=name, where_to_look=1),
                    rs.BlockInstanceCount(block_name=name, where_to_look=2),
                    ))
        else:
            print("(Name) ( (top level refs in active doc), (top level and nested refs in active doc), (refs from other inst defs) )")
            for name, rdIdef in sorted(zip(sBlockNames, rdIdefs_Out)):
                print("{} ({}, {}, {})".format(
                    name,
                    rs.BlockInstanceCount(block_name=name, where_to_look=0),
                    rs.BlockInstanceCount(block_name=name, where_to_look=1),
                    rs.BlockInstanceCount(block_name=name, where_to_look=2),
                    ))


if __name__ == '__main__': main()