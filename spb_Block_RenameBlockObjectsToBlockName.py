"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
151111-12: Created.
180830: Various code refactoring.
181204: Options added whether to name Block ... objects and to name them to Block ... container blocks.
250815: Refactored and added some options.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


_sBlockBlockNote = "'Block xx' block objects will be renamed to container block of 'Block xx', not to 'Block xx'.\nRename will not occur if container doesn't exist."


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bProcessAllBlocks'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRenameOnlyEmptyOrNone'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRenameBlockBlockObjs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bVerifyAll'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    for key in keys:
        if key not in names:
            names[key] = key[1:]


    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        #if key == 'fRadius':
        #    if cls.riOpts[key].CurrentValue < 2.0*sc.doc.ModelAbsoluteTolerance:
        #        cls.riOpts[key].CurrentValue = 0.0

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Invalid key?")
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get block instances with options.
    """

    if Opts.values['bProcessAllBlocks']:
        go = ri.Custom.GetOption()
        go.SetCommandPrompt("Set options")
        go.AcceptNothing(True)
    else:
        go = ri.Custom.GetObject()
        go.SetCommandPrompt("Select block instances")
        go.GeometryFilter = Rhino.DocObjects.ObjectType.InstanceReference
        # go.AcceptNothing(False) is already set.

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        addOption('bProcessAllBlocks')
        addOption('bRenameOnlyEmptyOrNone')
        addOption('bRenameBlockBlockObjs')
        addOption('bVerifyAll')
        addOption('bEcho')
        addOption('bDebug')


        if Opts.values['bProcessAllBlocks']:
            res = go.Get()
        else:
            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return [] # Means 'all block definitions'.

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if not go.Option():
            print("What happened?")
            continue

        # An option was selected.
        #sEval = 'go.Option().Index'; print(sEval,'=',eval(sEval))
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)

                if (
                    'bProcessAllBlocks' in idxs_Opts and
                    go.Option().Index == idxs_Opts['bProcessAllBlocks']):
                    go.Dispose()
                    if Opts.values['bProcessAllBlocks']:
                        go = ri.Custom.GetOption()
                        go.SetCommandPrompt("Set options")
                    else:
                        go = ri.Custom.GetObject()
                        go.SetCommandPrompt("Select block instances")
                        go.GeometryFilter = Rhino.DocObjects.ObjectType.InstanceReference

                elif (
                    'bRenameBlockBlockObjs' in idxs_Opts and
                    go.Option().Index == idxs_Opts['bRenameBlockBlockObjs']):

                    if Opts.values['bRenameBlockBlockObjs']:
                        print(_sBlockBlockNote)



                break


def main():

    sBlocks_All = rs.BlockNames(sort=True);
    if not sBlocks_All:
        print("No blocks in document.")
        return

    objrefs_Iref = getInput()
    if objrefs_Iref is None: return
    #sEval = 'objrefs_Iref'; print(sEval,'=',eval(sEval))

    bProcessAllBlocks = Opts.values['bProcessAllBlocks']
    bRenameOnlyEmptyOrNone = Opts.values['bRenameOnlyEmptyOrNone']
    bRenameBlockBlockObjs = Opts.values['bRenameBlockBlockObjs']
    bVerifyAll = Opts.values['bVerifyAll']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False


    #Individual block instance selection to replace contents of sBlocks_All
    if bProcessAllBlocks:
        sBlocks_ToProcess = sBlocks_All
    else:
        sBlocks_ToProcess = []
        for objref_Iref in objrefs_Iref:
            rgIref = objref_Iref.Geometry()
            rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)
            sBlocks_ToProcess.append(rdIdef.Name)


    if bRenameBlockBlockObjs:
        print(_sBlockBlockNote)

    # Process block definitions.

    print("Block name (New object name) ( (Objects renamed) of (Total objects) )")
    
    for sBlock in sBlocks_ToProcess:
        if bDebug: print(sBlock,)
        if bRenameBlockBlockObjs and sBlock[:6] == "Block " and sBlock[6].isnumeric():
            sContainers = rs.BlockContainers(sBlock)
            if bDebug: print(sContainers)
            ct_Containers = len(sContainers)
            if ct_Containers == 0:
                print("No container for {}.".format(sBlock))
                continue
            elif ct_Containers == 1:
                sObjName_Out = sContainer
            else:
                sObjName_Out = None
                for sContainer in sContainers:
                    # Verify that the block def. in sContainers is the immediate container block.
                    for g in rs.BlockObjects(sContainer):
                        if (
                                rs.ObjectType(g) == rs.filter.instance
                                and
                                rs.BlockInstanceName(g) == sBlock
                        ):
                            sObjName_Out = sContainer
                            break
                    if sObjName_Out is not None: break
                else:
                    print("No immediate container found for {}!".format(
                            sBlock))
        else:
            sObjName_Out = sBlock
        
        nObjsRenamed = 0
        gObjs = rs.BlockObjects(sBlock)
        if not gObjs:
            print("{} has no objects!".format(sBlock))
            continue
        
        gAlreadyNamedAsTarget = []
        
        for gObj in gObjs:
            sObjName_In = rs.ObjectName(gObj)
            if sObjName_In == sObjName_Out:
                gAlreadyNamedAsTarget.append(gObj)
                continue
            if (
                not bRenameOnlyEmptyOrNone or
                (not sObjName_In or sObjName_In.lower() == "none")
            ):
                if bVerifyAll:
                    sAnswer = rs.GetString(
                        "Rename {} to {}?".format(sObjName_In, sObjName_Out),
                        defaultString="Y",
                        strings=("Yes", "No"))
                    if sAnswer is None or sAnswer[0].lower() != 'y':
                        continue
                rs.ObjectName(gObj, sObjName_Out)
                sObjName1 = rs.ObjectName(gObj)
                if sObjName_In <> sObjName1:
                    nObjsRenamed += 1

                print("{} -> {}".format(sObjName_In, sObjName_Out))

        s = "In block def., {}, {} of {} objects were renamed to {}.".format(
            sBlock,
            nObjsRenamed,
            len(gObjs),
            sObjName_Out,
            )
        if gAlreadyNamedAsTarget:
            s += " {} objects already have that name.".format(len(gAlreadyNamedAsTarget))
        print(s)


if __name__ == '__main__': main()