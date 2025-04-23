"""
This script helps managing aliases.
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
161209: Created.
...
180729: Reabled the use of rs.OpenFileNames since the bug fix in Rhino 6.5.
...
211103: Added export.
220601: Added EditLines.
220602: Bug fix.
220917: Added ExportDefault and OpenAliasFolder options.  Modified a printed output.
231217: Now, output of Compare2Txt is sorted.
231222: Now, 2 alias .txt can be immediately selected in the file dialog and used.
        Added option to export ListPerSubstring to file.
        editMacro no longer asks for confirmation when only 1 macro is edited.
250210,13: Fixed codecs-related problems due to behavior of Python's readline vs. Rhino's _-Options _Aliases _Export.
250422: Some alias name/macro searches are now comma delimited AND.
"""

import Rhino
import Rhino.Input as ri
import rhinoscriptsyntax as rs

import codecs
import re


CAL = Rhino.ApplicationSettings.CommandAliasList

sCmdPrefix = 'x'
sDefaultPyFolderPath = "C:\\My\\Rhino\\PythonScript"
sDefaultAliasExportFolder = "C:\\My\\Rhino\\Aliases"

sOptions = (
    'Count',
    'AddPy',
    'AddLines',
    'DelLines',
    'DelListPick',
    'EditLines',
    'EditAlias',
    'EditMacro',
    'ListPerSubstring',
    'ListPyNotFound',
    'ListDupMacros',
    'ReassignPy',
    'Compare2Txt',
    'ExportDefault',
    'ExportBrowse',
    'OpenAliasFolder',
    )


def getOption():
    
    idxOpt_Default = 0
    
    go = ri.Custom.GetOption()
    go.SetCommandPrompt('Alias')
    go.SetCommandPromptDefault(defaultValue=sOptions[idxOpt_Default])
    go.AcceptNothing(True)
    # for idxOpt of 1 - 5, not 0.
    for sDirection in sOptions:
        go.AddOption(sDirection)
    #    for s in listPropOpts: gs.AddOption(s)
    res = go.Get()
    if res == ri.GetResult.Nothing:
        idxOpt = idxOpt_Default
    elif res == ri.GetResult.Option:
        idxOpt = go.Option().Index - 1 # '- 1' because go.Option().Index is base 1.
    else:
        return
    
    go.Dispose()
    
    Rhino.RhinoApp.SetCommandPrompt("Alias")
    
    return sOptions[idxOpt]


def fileBaseName(sFileFullPath):
    # Get full file name from full path.
    sSplit = sFileFullPath.split('\\')
    if len(sSplit) == 0: return
    sScriptFileFullName = sSplit[-1]
    # Remove '.' and extension.
    sSplit = sScriptFileFullName.split('.')
    if len(sSplit) != 2: return
    return sSplit[0]


def addAlias(sScriptFileBaseName, bOverwrite=False):
    """
    """
    
    sTitle = "Add Alias"
    
    if sScriptFileBaseName[:4] == 'spb_':
        sAlias = sCmdPrefix + sScriptFileBaseName[4].capitalize() + sScriptFileBaseName[5:]
    elif sScriptFileBaseName[0] != sCmdPrefix:
        # Add sCmdPrefix to beginning of command to distinguish them from Rhino commands, etc.
        sAlias = sCmdPrefix + sScriptFileBaseName[0].capitalize() + sScriptFileBaseName[1:]
    else:
        sAlias = sScriptFileBaseName[:]
    
    sAlias = rs.StringBox("Alias to add", sAlias, sTitle)
    if sAlias is None: return False
    
    if CAL.IsAlias(sAlias):
        sMacro = CAL.GetMacro(sAlias)
        nMsgBox = rs.MessageBox("{} exists with command macro\n{}\n\nReplace?".format(
                sAlias, sMacro), 3, sTitle)
        if nMsgBox == 2: return None
        elif nMsgBox == 7:
            print("{} was not modified.".format(sAlias))
            return False
    
    sMacro = "_NoEcho ! _-RunPythonScript " + sScriptFileBaseName
    
    sMacro = rs.StringBox("Macro for alias", sMacro, sTitle)
    if sMacro is None:
        print("{} was skipped.".format(sAlias))
        return False
    
    if CAL.Add(sAlias, sMacro):
        print("{} alias created with command macro {}".format(sAlias, sMacro))
        return sAlias, sMacro
    else:
        print("Error in creating alias for {}".format(sAlias))


def addAliases():
    """
    """
    
    # Placed in try due to 'External component has thrown an exception.' from OpenFileNames.
    try:
        sFileFullPaths = rs.OpenFileNames(
                title="Add Alias",
                filter="python scripts|*.py||",
                folder=sDefaultPyFolderPath,
        )
        if len(sFileFullPaths) == 0: return
    except:
        return
    
    sAliases = []
    sMacros = []
    
    for sFileFullPath in sFileFullPaths:
        sScriptFileBaseName = fileBaseName(sFileFullPath)
        if sScriptFileBaseName is None:
            print("Error in understanding {}   Exiting...".format(sFileFullPath))
            return
        rc = addAlias(sScriptFileBaseName, bOverwrite=True)
        if rc is None:
            print("Script canceled.")
            return
        elif rc is False:
            continue
        sAlias, sMacro = rc
        if rc is None:
            continue
        else:
            sAlias, sMacro = rc
            sAliases.append(sAlias)
            sMacros.append(sMacro)
    
    if sAliases:
        return sAliases, sMacros


def addLines_PerAliasExportFormat():
    """
    """
    
    sTitle = "Add Aliases in Alias Export format"

    bSuccess, sMultiInput = Rhino.UI.Dialogs.ShowEditBox(
        title=sTitle,
        message="Enter alias data in Alias Export format",
        defaultText=None,
        multiline=True)
    if not bSuccess: return

    for sSingleLine in sMultiInput.splitlines(keepends=False):
        sAlias, sMacro = sSingleLine.split(sep=' ', maxsplit=1)
        if sAlias is None:
            print("Format for {} is incorrect.  Alias name could not be determined.".format(
                sSingleLine))
            continue
        if sMacro is None:
            print("Format for {} is incorrect.  Marco could not be determined.".format(
                sSingleLine))
            continue

        if CAL.Add(sAlias, sMacro):
            print("'{}' alias created with command macro '{}'".format(sAlias, sMacro))
        else:
            print("Error in creating alias'{}'.".format(sAlias))


def _constructAliasLines_All():

    sNames = list(CAL.GetNames())
    sMacros = [CAL.GetMacro(s) for s in sNames]
    
    sLines_Starting = []
    
    for i in xrange(len(sNames)):
        sLine = sNames[i] + ' ' + sMacros[i]
        sLines_Starting.append(sLine)

    if not sLines_Starting:
        print("No aliases!")
        return

    return sLines_Starting


def _splitLineTo_AliasName_and_Macro(sLine):
    if ' ' not in sLine:
        print("Line has no spaces.")
        return

    return sLine.split(sep=' ', maxsplit=1)


def deleteLines_PerAliasExportFormat():
    """
    """
    
    sTitle = "Delete Aliases in Alias Export format"
    
    bSuccess, sMultiInput = Rhino.UI.Dialogs.ShowEditBox(
        title=sTitle,
        message="Enter alias data in Alias Export format",
        defaultText=None,
        multiline=True)
    if not bSuccess: return


    sLines_ToDel = sMultiInput.splitlines(keepends=False)


    sLines_Starting = _constructAliasLines_All()
    if not sLines_Starting: return


    sLines_Deleted = []

    for sLine_ToDel in sLines_ToDel:
        if sLine_ToDel in sLines_Starting:
            rc = _splitLineTo_AliasName_and_Macro(sLine_ToDel)
            if not rc: continue
            sAlias, sMacro = rc
            if CAL.Delete(alias=sAlias):
                print("Deleted {}".format(sLine_ToDel))
                sLines_Deleted.append(sLine_ToDel)
        else:
            print("{} not in aliases.".format(sLine_ToDel))
    
    
    if len(sLines_Deleted) == len(sLines_ToDel):
        print("Deleted all {} alias/macro lines.".format(len(sLines_Deleted)))
    else:
        print("Deleted {} out of {} intended alias/macro lines.".format(
            len(sLines_Deleted), len(sLines_ToDel)))


def _replace_string_case_insensitive(text, old_string, new_string):
    """ From Gemini."""
    return re.sub(re.escape(old_string), new_string, text, flags=re.IGNORECASE)


def _areAllStringsInText(text, strings):
    for string in strings:
        if not string.lower() in text.lower():
            return False

    return True


def _removeCommonMacroTextFromString(s_In):
    ss_ToRemove = (
        "noecho",
        "loadscript",
        "runpythonscript",
        "runscript",
        "scripteditor",
        )

    s_Out = s_In
    for s_ToRemove in ss_ToRemove:
        s_Out = _replace_string_case_insensitive(s_Out, s_ToRemove, "")

    return s_Out


def getAliasMacroLines(sTitle, bAllOpt=False):
    """
    """

    if bAllOpt:
        sStringForAll = "__ListAll__"
    
        sSearchSubString = rs.StringBox(
                message="Enter substring to find"
                "\n\nSearch is case insensitive."
                "\nNoEcho, RunPythonScript, etc., are ignored."
                "\nUse commas to delimit strings for an AND search.",
                default_value=sStringForAll,
                title=sTitle)
    else:
        sSearchSubString = rs.StringBox(
                message="Enter substring to find",
                default_value=None,
                title=sTitle)

    if sSearchSubString is None: return
    
    
    sLines = []
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)


    if bAllOpt and (sSearchSubString == sStringForAll):
        for sName in sNames:
            sMacro = CAL.GetMacro(sName)
            sLine = sName + ' ' + sMacro
            sLines.append(sLine)
        if not sLines:
            print("No aliases.")
            return
        print("Found {} aliases.".format(len(sLines)))
        return sLines

    ssSearchSubString = sSearchSubString.split(',')

    # Only include macro lines with entered substring.
    for sName in sNames:
        sMacro = CAL.GetMacro(sName)

        sMacro_Cleaned = _removeCommonMacroTextFromString(sMacro)

        if _areAllStringsInText(sName, ssSearchSubString):
            pass
        elif _areAllStringsInText(sMacro_Cleaned, ssSearchSubString):
            pass
        else:
            continue

        sLine = sName + ' ' + sMacro
        sLines.append(sLine)

    if not sLines:
        print("Substring not found.")
        return
    print("Found {} aliases with substring {}.".format(
        len(sLines),
        sSearchSubString))

    return sLines


def editLines():
    """
    """

    sTitle = "Edit Aliase/Macro Lines"

    sLines = getAliasMacroLines(sTitle, bAllOpt=False)
    if not sLines: return

    bSuccess, sMultiLines_FromEdit = Rhino.UI.Dialogs.ShowEditBox(
        title=sTitle,
        message="Edit alias data per Alias Export format",
        defaultText='\n'.join(sLines),
        multiline=True)
    if not bSuccess: return

    sLines_FromEdit = sMultiLines_FromEdit.splitlines(keepends=False)

    sAliasNames_BeforeEdit = list(CAL.GetNames())

    sLines_All_BeforeEdit = _constructAliasLines_All()
    if not sLines_All_BeforeEdit: return


    sLines_Edited = []

    for sLine_ToAdd in sLines_FromEdit:
        if sLine_ToAdd in sLines_All_BeforeEdit:
            # No change.
            continue
        rc = _splitLineTo_AliasName_and_Macro(sLine_ToAdd)
        if not rc: continue
        sAlias, sMacro = rc

        if CAL.Add(sAlias, sMacro):
            if sAlias in sAliasNames_BeforeEdit:
                print("Modified {}".format(sLine_ToAdd))
            else:
                print("Added {}".format(sLine_ToAdd))

            sLines_Edited.append(sLine_ToAdd)
        else:
            print("{} not added.".format(sLine_ToAdd))


    print("Edited {} out of {} alias/macro lines.".format(
        len(sLines_Edited), len(sLines_FromEdit)))


def editAlias():
    """
    """
    
    sTitle = "Edit Alias Name"
    
    sSearchSubString = rs.StringBox(
            message="Enter substring to find or Enter (or Cancel) to list all",
            default_value=None,
            title=sTitle)
    
    sItems = []
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    for sName in sNames:
        sMacro = CAL.GetMacro(sName)
        sLine = sName + ' ' + sMacro
        sItems.append(sLine)   
    
    if sSearchSubString is not None:
        sItems_WIP = []
        for sItem in sItems:
            if sSearchSubString.lower() in sItem.lower():
                sItems_WIP.append(sItem)
        if not sItems_WIP:
            print("Substring not found.")
            return
        sItems = sItems_WIP
    
    sItemsSelected = rs.MultiListBox(
            items=sItems,
            message="Pick aliases to edit",
            title=sTitle,
            defaults=None)
    if sItemsSelected is None: return
    
    for sItem in sItemsSelected:
        sAliasName_Sel = sItem.split(sep=' ')[0]
        
        sAliasName_New = rs.StringBox(
                message="Edit alias",
                default_value=sAliasName_Sel,
                title=sTitle)
        if not sAliasName_New: return
        if sAliasName_New == sAliasName_Sel: continue
        
        iMbReturn = rs.MessageBox(
                message="Replace alias\n{}\nwith\n{}\n?".format(
                        sAliasName_Sel,
                        sAliasName_New),
                buttons=3, title='')
        if iMbReturn == 2:
            # Cancel
            return
        elif iMbReturn == 6:
            # Yes
            sMacro = CAL.GetMacro(alias=sAliasName_Sel)
            if sMacro:
                if CAL.Add(alias=sAliasName_New, macro=sMacro):
                    if CAL.Delete(alias=sAliasName_Sel):
                        print("Alias {} replaced with {}.".format(
                                sAliasName_Sel,
                                sAliasName_New))
                    else:
                        print("Alias {} macro WAS NOT DELETED (PHASE 1 OF 2 OF EDIT).".format(
                                sAliasName_Sel))
            else:
                print("Alias {} macro WAS NOT ADDED (PHASE 1 OF 2 OF EDIT).".format(
                        sAliasName_New))
        elif iMbReturn == 7:
            # No
            continue


def editMacro():
    """
    """
    
    sTitle = "Edit Alias Macro"
    
    sSearchSubString = rs.StringBox(
            message="Enter substring to find or Enter (or Cancel) to list all",
            default_value=None,
            title=sTitle)
    
    sItems = []
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    for sName in sNames:
        sMacro = CAL.GetMacro(sName)
        sLine = sName + ' ' + sMacro
        sItems.append(sLine)   
    
    if sSearchSubString is not None:
        sItems_WIP = []
        for sItem in sItems:
            if sSearchSubString.lower() in sItem.lower():
                sItems_WIP.append(sItem)
        if not sItems_WIP:
            print("Substring not found.")
            return
        sItems = sItems_WIP
    
    if len(sItems) > 1:
        sItemsSelected = rs.MultiListBox(
                items=sItems,
                message="Pick aliases to edit their macros",
                title=sTitle,
                defaults=None)
        if sItemsSelected is None: return
    else:
        sItemsSelected = sItems
    
    for sItem in sItemsSelected:
        sAliasName_Sel = sItem.split(sep=' ')[0]
        
        sMacro_Old = CAL.GetMacro(alias=sAliasName_Sel)
        
        sMacro_New = rs.StringBox(
                message="Edit macro",
                default_value=sMacro_Old,
                title=sTitle)
        if not sMacro_New: return
        if sMacro_New == sMacro_Old: continue
        
        if len(sItemsSelected) == 1:
            iMbReturn = 6
        else:
            iMbReturn = rs.MessageBox(
                    message="Replace macro\n{}\nwith\n{}\n?".format(
                            sMacro_Old,
                            sMacro_New),
                    buttons=3, title='')
        if iMbReturn == 2:
            # Cancel
            return
        elif iMbReturn == 6:
            # Yes
            sMacro_Prev = CAL.GetMacro(alias=sAliasName_Sel)
            if not CAL.SetMacro(alias=sAliasName_Sel, macro=sMacro_New):
                sMacro_Prev = None
            if sMacro_Prev == sMacro_Old:
                print("Alias {} macro {} replaced with {}.".format(
                        sAliasName_Sel,
                        sMacro_Prev,
                        CAL.GetMacro(alias=sAliasName_Sel)))
            else:
                print("Alias {} macro WAS NOT EDITED.".format(sAliasName_Sel))
        elif iMbReturn == 7:
            # No
            continue


def delete_Single_OLD():
    """
    """
    
    sTitle = "Delete Alias"
    
    sSearchSubString = rs.StringBox(
            message="Enter substring to find or Enter (or Cancel) to list all",
            default_value=None,
            title=sTitle)
    
    sItems = []
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    for sName in sNames:
        sMacro = CAL.GetMacro(sName)
        sLine = sName + ' ' + sMacro
        sItems.append(sLine)   
    
    if sSearchSubString is not None:
        sItems_WIP = []
        for sItem in sItems:
            if sSearchSubString.lower() in sItem.lower():
                sItems_WIP.append(sItem)
        if not sItems_WIP:
            print("Substring not found.")
            return
        sItems = sItems_WIP
    
    sItemsSelected = rs.MultiListBox(
            items=sItems,
            message="Pick aliases to delete",
            title=sTitle,
            defaults=None)
    if sItemsSelected is None: return
    
    for sItem in sItemsSelected:
        sAliasName_Sel = sItem.split(sep=' ')[0]
        iMbReturn = rs.MessageBox(
                message="Delete alias?:\n{}".format(sAliasName_Sel),
                buttons=3, title='')
        if iMbReturn == 2:
            # Cancel
            return
        elif iMbReturn == 6:
            # Yes
            if CAL.Delete(alias=sAliasName_Sel):
                print("Alias {} was deleted.".format(sAliasName_Sel))
            else:
                print("Alias {} COULD NOT BE DELETED.".format(sAliasName_Sel))
        elif iMbReturn == 7:
            # No
            continue


def deletePickedFromList(sLines, sTitle=None):
    """
    """
    
    rc = Rhino.UI.Dialogs.ShowMultiListBox(
        title=sTitle,
        message="Select macros to delete",
        items=sLines,
        defaults=None)
    
    if rc is None:
        for s in sLines:
            print(s)
        return
    
    sLines_toDelete = list(rc)
    
    sLines_Deleted = []
    
    for sLine_toDelete in sLines_toDelete:
        sName = sLine_toDelete.split(" ",1)[0]
        #print(sName
        if CAL.Delete(alias=sName):
            print("Deleted {}".format(sLine_toDelete))
            sLines_Deleted.append(sLine_toDelete)
    
    return sLines_Deleted


def exportCustomList(sExport):
    import os
    sDesktop = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop') 
    sFilePath_Out = rs.SaveFileName(
        title="Save custom list",
        filter=None,
        folder=sDesktop,
        filename=None,
        extension="txt")
    if not sFilePath_Out: return
    with open(sFilePath_Out, 'w') as f:
        f.write('\n'.join(sExport))


def listAliasesMacrosWithSubstring():
    """
    """

    sLines = getAliasMacroLines("List Aliases", bAllOpt=True)
    if not sLines: return

    for sLine in sLines:
        print(sLine)
    
    if len(sLines) < CAL.Count:
        print("{} alias/macro lines listed out of {} total.".format(len(sLines), CAL.Count))
    else:
        print("All {} alias/macro lines listed.".format(len(sLines)))
    
    if len(sLines) > 100:
        sDecision = rs.GetString("Export to text file?", defaultString="No", strings=["Yes", "No"])
        if sDecision and sDecision.lower()[0] == 'y':
            exportCustomList(sLines)


def listAliasesMacrosWithSubstringForDelete():
    """
    """
    
    sTitle = "List Aliases for Delete"
    sStringForAll = "__ListAll__"
    
    sSearchSubString = rs.StringBox(
            message="Enter substring to find or Enter (or Cancel) to list all",
            default_value=sStringForAll,
            title=sTitle)
    
    if sSearchSubString is None: return
    
    
    sLines = []
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    
    if sSearchSubString == sStringForAll:
        for sName in sNames:
            sMacro = CAL.GetMacro(sName)
            sLine = sName + ' ' + sMacro
            sLines.append(sLine)
        if not sLines:
            print("No aliases.")
            return
        print("Found {} aliases.".format(len(sLines)))
    else:
        # Only include macro lines with entered substring.
        for sName in sNames:
            sMacro = CAL.GetMacro(sName)
            sLine = sName + ' ' + sMacro
            
            if sSearchSubString.lower() in sLine.lower():
                sLines.append(sLine)
        
        if not sLines:
            print("Substring not found.")
            return
        print("Found {} aliases with substring {}.".format(
            len(sLines),
            sSearchSubString))
    
    rc = deletePickedFromList(sLines, sTitle=sTitle)
    
    if rc: return


def listPythonScriptsNotFound():
    """
    """
    
    sTitle = "List Aliases Referencing Python Scripts Not Found"
    
    import os
    
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    list_sMacroLines_with_files_not_found = []
    
    for sName in sNames:
        
        sMacro = CAL.GetMacro(sName)
        if 'RunPythonScript' in sMacro:
            sLine = sName + ' ' + sMacro
            #print(sLine
            
            sSplit = sMacro.split("RunPythonScript ",1)
            
            if len(sSplit) < 2:
                continue
            
            sAfterRunPythonScript = sSplit[1]
            
            sScriptNameOrPath = sAfterRunPythonScript.split(" ",1)[0]
            
            if '"' in sScriptNameOrPath:
                if sScriptNameOrPath[-1] != '"':
                    continue
                sScriptNameOrPath = sScriptNameOrPath[:-1] + '.py"'
                pass 
            
            elif sScriptNameOrPath[-3:] != '.py':
                sScriptNameOrPath += '.py'
            else:
                pass
            
            if not os.path.isfile(sScriptNameOrPath):
                #print(sScriptNameOrPath
                list_sMacroLines_with_files_not_found.append(sLine)
    
    if not list_sMacroLines_with_files_not_found:
        print("No python scripts with missing files.")
        return
    
    print("Found {} aliases with Python scripts not found.".format(
        len(list_sMacroLines_with_files_not_found)))
    
    deletePickedFromList(list_sMacroLines_with_files_not_found, sTitle=sTitle)


def listDuplicateMacros():
    """
    """
    
    sTitle = "List and Optionally Delete Aliases with Macro Repeats"
    
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    sMacros = [CAL.GetMacro(s) for s in sNames]
    
    sMacros_Repeats = []
    
    sMacros_WithoutRepeats = list(set(sMacros))
    
    if len(sMacros_WithoutRepeats) == len(sMacros):
        print("No duplicate macros.")
        return
    
    sMacros_WithoutRepeats.sort(key=str.lower)
    
    for i in xrange(len(sMacros_WithoutRepeats)):
        sMacro = sMacros_WithoutRepeats[i]
        if sMacros.count(sMacro) > 1:
            if sMacro not in sMacros_Repeats:
                sMacros_Repeats.append(sMacro)
    
    sLines_with_RepeatMacros = []
    
    for sMacro_Repeat in sMacros_Repeats:
        for i in xrange(len(sNames)):
            sMacro = sMacros[i]
            if sMacro == sMacro_Repeat:
                sLine = sNames[i] + ' ' + sMacro
                sLines_with_RepeatMacros.append(sLine)
    
    print("Found {} aliases with repeat macros.".format(len(sLines_with_RepeatMacros)))
    
    deletePickedFromList(sLines_with_RepeatMacros, sTitle=sTitle)


def reAssignMacro():
    """
    """
    
    sTitle = "Reassign Alias Macro to Python Script"
    
    sSearchSubString = rs.StringBox(
            message="Enter substring to find or Enter (or Cancel) to list all",
            default_value=None,
            title=sTitle)
    
    sItems = []
    sNames = list(CAL.GetNames())
    sNames.sort(key=str.lower)
    
    for sName in sNames:
        sMacro = CAL.GetMacro(sName)
        sLine = sName + ' ' + sMacro
        sItems.append(sLine)   
    
    if sSearchSubString is not None:
        sItems_WIP = []
        for sItem in sItems:
            if sSearchSubString.lower() in sItem.lower():
                sItems_WIP.append(sItem)
        if not sItems_WIP:
            print("Substring not found.")
            return
        sItems = sItems_WIP
    
    sItemsSelected = rs.MultiListBox(
            items=sItems,
            message="Pick aliases to edit their macros",
            title=sTitle,
            defaults=None)
    if sItemsSelected is None: return
    
    for sItem in sItemsSelected:
        sAliasName_Sel = sItem.split(sep=' ')[0]
        
        sMacro_Old = CAL.GetMacro(alias=sAliasName_Sel)
        
        sFileFullPath = rs.OpenFileName(
                title="Reassign Macro for Alias {}".format(sAliasName_Sel),
                filter="python scripts|*.py||",
                folder=sDefaultPyFolderPath,
        )
        if len(sFileFullPath) == 0: return
        
        sAliases = []
        sMacros = []
        
        sScriptFileBaseName = fileBaseName(sFileFullPath)
        if sScriptFileBaseName is None:
            print("Error in understanding {}   Exiting...".format(sFileFullPath))
            return
        
        sMacro_New = "_NoEcho ! _-RunPythonScript " + sScriptFileBaseName
        
        sMacro_New = rs.StringBox("Macro for alias", sMacro_New, sTitle)
        if sMacro_New is None:
            print("{} was skipped.".format(sAliasName_Sel))
            return False
        
        if sMacro_New == sMacro_Old: continue
        
        iMbReturn = rs.MessageBox(
                message="Replace macro\n{}\nwith\n{}\n?".format(
                        sMacro_Old,
                        sMacro_New),
                buttons=3, title='')
        if iMbReturn == 2:
            # Cancel
            return
        elif iMbReturn == 6:
            # Yes
            sMacro_Prev = CAL.GetMacro(alias=sAliasName_Sel)
            if not CAL.SetMacro(alias=sAliasName_Sel, macro=sMacro_New):
                sMacro_Prev = None
            if sMacro_Prev == sMacro_Old:
                print("Alias {} macro {} replaced with {}.".format(
                        sAliasName_Sel,
                        sMacro_Prev,
                        CAL.GetMacro(alias=sAliasName_Sel)))
            else:
                print("Alias {} macro WAS NOT EDITED.".format(sAliasName_Sel))
        elif iMbReturn == 7:
            # No
            continue


def _remove_bom(text):
    if text.startswith(codecs.BOM_UTF8):
        return text.decode('utf-8-sig')
    return text


def compare2AliasExports():
    """
    """
    
    sTitle = "Select 2 alias exports to compare"
    
    # Placed in try due to 'External component has thrown an exception.' from OpenFileNames.
    try:
        sFilePaths_In = rs.OpenFileNames(
            title=sTitle,
            filter="python scripts|*.txt||",
            folder=sDefaultAliasExportFolder,
        )
    except:
        return
    
    dict_sFilePaths = {}
    
    if len(sFilePaths_In) == 0: return
    dict_sFilePaths['A'] = sFilePaths_In[0]
    if len(sFilePaths_In) > 1:
        dict_sFilePaths['B'] = sFilePaths_In[1]
    else:
        # Placed in try due to 'External component has thrown an exception.' from OpenFileNames.
        try:
            sFilePaths_In = rs.OpenFileNames(
                    title=sTitle,
                    filter="python scripts|*.txt||",
                    folder=sDefaultAliasExportFolder,
            )
            if len(sFilePaths_In) == 0: return
            sFilePath_B = sFilePaths_In[0]
        except:
            return

    dict_Files = {}
    dict_Lines = {}

    for sFile in ('A', 'B'):
        with open(dict_sFilePaths[sFile]) as dict_Files[sFile]:
            dict_Lines[sFile] = dict_Files[sFile].readlines()

    if dict_Lines['A'] == dict_Lines['B']:
        print("Files contents are identical, including line order.")
        return

    for sFile in ('A', 'B'):
        sLines_X = dict_Lines[sFile]
        sLine_Cleaned = _remove_bom(sLines_X[0])
        if sLine_Cleaned != sLines_X[0]:
            sLine_First_Save = sLines_X[0]
            sLines_X[0] = sLine_Cleaned
            if sLine_Cleaned != sLines_X[0]:
                raise Exception("Could no remove byte order mark.")
            else:
                print("Removed byte order mark from the first line of {}: {} -> {}".format(
                    dict_sFilePaths[sFile],
                    sLine_First_Save[:-1],
                    sLines_X[0][:-1]),
                    )

    dict_Lines['A'] = set(dict_Lines['A'])
    dict_Lines['B'] = set(dict_Lines['B'])

    if dict_Lines['A'] == dict_Lines['B']:
        print("Ignoring line order, files contents are identical.")
        return

    print('-'*80)
    sLines_Delta = sorted(dict_Lines['A'] - dict_Lines['B'])
    if sLines_Delta:
        print("In {} only:".format(dict_sFilePaths['A']))
        print('-'*40)
        s = ""
        for sLine in sLines_Delta:
            s += sLine
        print(s)
        print('-'*80)
    
    sLines_Delta = sorted(dict_Lines['B'] - dict_Lines['A'])
    if sLines_Delta:
        print("In {} only:".format(dict_sFilePaths['B']))
        print('-'*40)
        s = ""
        for sLine in sLines_Delta:
            s += sLine
        print(s)
        print('-'*80)


def exportDefault():
    """
    Creates file names such as Aliases-V7-Main-220917.txt
    """

    from datetime import datetime
    from os import environ

    sFileName = 'Aliases-V{}-{}-{}.txt'.format(
        Rhino.RhinoApp.Version.Major,
        environ.get("USERNAME"),
        datetime.today().strftime('%y%m%d'))

    sPath = '"' + sDefaultAliasExportFolder + '\\' + sFileName + '"'

    script = "! _-Options _Aliases _Export {} _EnterEnd".format(sPath)

    Rhino.RhinoApp.RunScript(script, echo=True)


def exportBrowse():
    script = "! _-Options _Aliases _Export _Browse _EnterEnd"
    Rhino.RhinoApp.RunScript(script, echo=True)


def openAliasFolder():
    import subprocess
    subprocess.Popen('explorer "{}"'.format(sDefaultAliasExportFolder))


def main():
    """
    """
    
    rc = getOption()
    if rc is None: return
    sOption = rc
    
    if sOption == 'Count':
        print("{} aliases".format(CAL.GetNames().Count))
    elif sOption == 'AddPy':
        addAliases()
    elif sOption == 'AddLines':
        addLines_PerAliasExportFormat()
    elif sOption == 'DelLines':
        deleteLines_PerAliasExportFormat()
    elif sOption == 'DelListPick':
        listAliasesMacrosWithSubstringForDelete()
    elif sOption == 'EditLines':
        editLines()
    elif sOption == 'EditAlias':
        editAlias()
    elif sOption == 'EditMacro':
        editMacro()
    elif sOption == 'ListPerSubstring':
        listAliasesMacrosWithSubstring()
    elif sOption == 'ListPyNotFound':
        listPythonScriptsNotFound()
    elif sOption == 'ListDupMacros':
        listDuplicateMacros()
    elif sOption == 'ReassignPy':
        reAssignMacro()
    elif sOption == 'Compare2Txt':
        compare2AliasExports()
    elif sOption == 'ExportDefault':
        exportDefault()
    elif sOption == 'ExportBrowse':
        exportBrowse()
    elif sOption == 'OpenAliasFolder':
        openAliasFolder()
    else:
        raise ValueError('')


if __name__ == '__main__': main()