"""
Wrapper for _Export with options to build file name and pick a common output folder.

IGES files are exported with current document units and absolute model tolerance (distance).

_SaveTextures is set to No when SaveSmall option is Yes.
Otherwise, _SaveTextures is not modified from its current setting.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
210505: Created.
...
240229: Changed temporary folder path.
240612: Maximum Rhino version is no long hard-coded, but the current Rhino version.
250120: Moved some command line options into a separate GetOptions.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import DateTime

import os


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iFileType'; keys.append(key)
    listValues[key] = ( 
        'Rhino3dm',
        'STEP',
        'IGES',
        'Parasolid',
        'DWG',
        'DXF',
        ) # All items must be strings.
    values[key] = 0
    names[key] = 'ExportType'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'i3dmVersion'; keys.append(key)
    values[key] = Rhino.RhinoApp.ExeVersion
    riOpts[key] = ri.Custom.OptionInteger(values[key], 2, Rhino.RhinoApp.ExeVersion)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSaveSmall'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSaveNotes'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseCurrentOpts'; keys.append(key)
    values[key] = True
    names[key] = 'FormatOpts'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Set', 'UseCurrent')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bManuallyEditName'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCustomPrefix'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'sPrefix'; keys.append(key)
    values[key] = ""
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSourceFileName'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPrecedeSFNameWith_From_'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDate'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTime'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bFromRhVer'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bToFormat'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCustomSuffix'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'sSuffix'; keys.append(key)
    values[key] = "-"
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iFolder'; keys.append(key)
    if sc.doc.Path is None:
        listValues[key] = ( 
            'TempTranslations',
            'Desktop',
            ) # All items must be strings.
        values[key] = 0
    else:
        listValues[key] = ( 
            'TempTranslations',
            'Desktop',
            'DocFolder',
            ) # All items must be strings.
        values[key] = 2
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
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, value=None):

        #if key == 'fScale':
        #    if cls.riOpts[key].CurrentValue <= 1e-9:
        #        print "Invalid input for scale value."
        #        cls.riOpts[key].CurrentValue = cls.values[key]

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        else:
            cls.values[key] = value

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getObjects():
    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        prompt="Select objects to export",
        acceptNothing=False,
        filter=rd.ObjectType.AnyObject)

    if res == Rhino.Commands.Result.Cancel: return

    return objrefs


def _addExtension(sFileName):
    sFileType = Opts.listValues['iFileType'][Opts.values['iFileType']]

    if sFileType == 'Rhino3dm':
        sFileName += ".3dm"
    elif sFileType == 'STEP':
        sFileName += ".stp"
    elif sFileType == 'IGES':
        sFileName += ".igs"
    elif sFileType == 'Parasolid':
        sFileName += ".x_t"
    elif sFileType == 'DWG':
        sFileName += ".dwg"
    elif sFileType == 'DXF':
        sFileName += ".dxf"
    else:
        raise Exception("File type {} is not supported.".format(sFileType))

    return sFileName


def createFileName(bWithExt=False, **kwargs):
    """
    File name does not include the file type extension.
    
    Parameters:
        kwargs (all optional)
            iFileType
            i3dmVersion
            bManuallyEditName
            bSaveSmall
            bSaveNotes
            bPrecedeSFNameWith_From_
            bSourceFileName
            bDate
            bTime
            bFromRhVer
            bToFormat
            iFolder
            bEcho
            bDebug
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    #iFileType = getOpt('iFileType')
    bManuallyEditName = getOpt('bManuallyEditName')
    i3dmVersion = getOpt('i3dmVersion')
    sPrefix = getOpt('sPrefix') if getOpt('bCustomPrefix') else ""
    bSourceFileName = getOpt('bSourceFileName')
    bPrecedeSFNameWith_From_ = getOpt('bPrecedeSFNameWith_From_') if bSourceFileName else False
    bDate = getOpt('bDate')
    bTime = getOpt('bTime')
    bFromRhVer = getOpt('bFromRhVer')
    bToFormat = getOpt('bToFormat')
    sSuffix = getOpt('sSuffix') if getOpt('bCustomSuffix') else ""
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    sFileType = Opts.listValues['iFileType'][Opts.values['iFileType']]

    sFileName = ""

    if sPrefix:
        sFileName += sPrefix

    if bSourceFileName:
        if bPrecedeSFNameWith_From_:
            sFileName += "_From_"
        if sc.doc.Name is None:
            sFileName += "Untitled"
        else:
            sFileName += sc.doc.Name[:-4]
    elif not (
        bDate and
        bTime
        ):
            sFileName += "Untitled"


    if bDate:
        if sFileName: sFileName += "_"
        if bTime:
            sFileName += DateTime.Now.ToString("yyMMdd-HHmmss")
        else:
            sFileName += DateTime.Now.ToString("yyMMdd")

    if bFromRhVer:
        sFileName += "_From_Rh{}".format(Rhino.RhinoApp.ExeVersion)

    if sSuffix:
        sFileName += sSuffix


    if bToFormat:
        if sFileType == 'Rhino3dm':
            sFileName += "_to_Rh{}_3dm".format(i3dmVersion)
        elif sFileType == 'STEP':
            sFileName += "_to_stp".format(i3dmVersion)
        elif sFileType == 'IGES':
            sFileName += "_to_igs".format(i3dmVersion)
        elif sFileType == 'Parasolid':
            sFileName += "_to_x_t".format(i3dmVersion)
        elif sFileType == 'DWG':
            sFileName += "_to_dwg".format(i3dmVersion)
        elif sFileType == 'DXF':
            sFileName += "_to_dxf".format(i3dmVersion)
        else:
            raise ValueError("{} file type is not supported.".format(sFileType))


    if bWithExt:
        sFileName = _addExtension(sFileName)


    return sFileName


def setOptions():
    """
    """

    sFileName = createFileName(bWithExt=True)
    print("Current file name: {}".format(sFileName))

    go = ri.Custom.GetOption()

    go.SetCommandPrompt("Export options")

    go.AcceptNothing(True)

    idxs_Opts = {}

    def addOption(ric, key): idxs_Opts[key] = Opts.addOption(ric, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        if Opts.values['bCustomPrefix']:
            print('Current custom prefix: "{}"'.format(Opts.values['sPrefix']))
        if Opts.values['bCustomSuffix']:
            print('Current custom suffix: "{}"'.format(Opts.values['sSuffix']))


        addOption(go, 'iFileType')


        if Opts.listValues['iFileType'][Opts.values['iFileType']] == 'Rhino3dm':
            addOption(go, 'i3dmVersion')
            addOption(go, 'bSaveSmall')
        else:
            addOption(go, 'bUseCurrentOpts')

        key = 'AutoNamingSettings'; idxs_Opts[key] = go.AddOption(key)
        addOption(go, 'bManuallyEditName')
        addOption(go, 'iFolder')
        addOption(go, 'bEcho')
        addOption(go, 'bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return True


        if go.Option().Index == idxs_Opts['AutoNamingSettings']:

            idxs_Opts.clear()

            go_Naming = ri.Custom.GetOption()
            go_Naming.SetCommandPrompt("Auto naming settings")

            addOption(go_Naming, 'bCustomPrefix')
            addOption(go_Naming, 'bSourceFileName')
            if Opts.values['bSourceFileName']:
                addOption(go_Naming, 'bPrecedeSFNameWith_From_')
            addOption(go_Naming, 'bDate')
            if Opts.values['bDate']:
                addOption(go_Naming, 'bTime')
            addOption(go_Naming, 'bFromRhVer')
            addOption(go_Naming, 'bToFormat')
            addOption(go_Naming, 'bCustomSuffix')

            while True:
                res = go_Naming.Get()

                if res != ri.GetResult.Option:
                    break

                key = 'bCustomPrefix'
                if go_Naming.Option().Index == idxs_Opts[key]:
                    Opts.setValue(key, go_Naming.Option().CurrentListOptionIndex)
                    if Opts.values[key]:
                        res, text = Rhino.UI.Dialogs.ShowEditBox(
                            title="Export",
                            message="Enter string for file name prefix",
                            defaultText=Opts.values['sPrefix'],
                            multiline=False)
                        if not res:
                            continue
                        Opts.setValue('sPrefix', text)
                    continue

                key = 'bCustomSuffix'
                if go_Naming.Option().Index == idxs_Opts[key]:
                    Opts.setValue(key, go_Naming.Option().CurrentListOptionIndex)
                    if Opts.values[key]:
                        res, text = Rhino.UI.Dialogs.ShowEditBox(
                            title="Export",
                            message="Enter string for file name suffix",
                            defaultText=Opts.values['sSuffix'],
                            multiline=False)
                        if not res:
                            continue
                        Opts.setValue('sSuffix', text)
                    continue


                for key in idxs_Opts:
                    if go_Naming.Option().Index == idxs_Opts[key]:
                        Opts.setValue(key, go_Naming.Option().CurrentListOptionIndex)
                        break

            go_Naming.Dispose()

            sFileName = createFileName(bWithExt=True)
            print("Current file name: {}".format(sFileName))

            continue


        # Another option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


"""
        Case "dwg"
            strExportOpts = ""
        Case "dxf"
            strExportOpts = "_Scheme 2004SPB"
        Case "igs"
            strExportOpts = Chr(34) & "SolidWorks surfaces" & Chr(34)
        Case "x_t"
        ''
        Case "stp"
            strExportOpts = "_Schema=AP214AutomotiveDesign"
        Case Else
            Call Rhino.Print("Export file extension, " & strFileExt & _
            ", not supported.  Exiting...")
            Exit Sub
        End Select
"""


def exportDocObjects(rhObjs, sFileName, **kwargs):
    """
    Parameters:
        rhObjs: RhinoDocObjects ObjRefs, RhinoDocObjects ...Objects, or GUIDs
        sFileType: str("Rhino3dm", "STEP", etc.)
        kwargs (all optional)
            i3dmVersion
            bSaveSmall
            bSaveNotes
            bSourceFileName
            bDate
            bTime
            bFromRhVer
            bToFormat
            iFolder
            bEcho
            bDebug
    """

    #if sFileType not in Opts.listValues['iFileType']:
    #    print("{} format not supported.  Supported: {}".format(
    #        Opts.listValues['iFileType']))
    #    return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    i3dmVersion = getOpt('i3dmVersion')
    bSaveSmall = getOpt('bSaveSmall')
    bSaveNotes = getOpt('bSaveNotes')
    bUseCurrentOpts = getOpt('bUseCurrentOpts')
    iFolder = getOpt('iFolder')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    sOutputPath = '"'

    if iFolder == 0 or (iFolder == 2 and sc.doc.Path is None):
        sDesktopPath = "C:/TempShared/Translations"
        sOutputPath += sDesktopPath + "\\"
    #if iFolder == 0 or (iFolder == 2 and sc.doc.Path is None):
    #    sDesktopPath = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')
    #    sOutputPath += sDesktopPath + "\\"
    elif iFolder == 1:
        sDesktopPath = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')
        sOutputPath += sDesktopPath + "\\"
    elif iFolder == 2:
        sOutputPath += sc.doc.Path[:-len(sc.doc.Name)]
    else:
        raise Exception("What happened?")


    sOutputPath += sFileName

    sOutputPath += '"'


    if os.path.exists(sOutputPath.strip('"')):
        print("{} exists.\n File will not be exported.".format(sOutputPath))
        return

    if bDebug:
        print(sOutputPath)

    sc.doc.Objects.UnselectAll()

    sc.doc.Objects.Select(rhObjs)

    sCommand = "! _-Export _Pause "

    if sFileName.split('.')[-1].lower() == '3dm':
        sCommand += "_Version={} ".format(i3dmVersion)
        sCommand += "_SaveSmall={} ".format("Yes" if bSaveSmall else "No")
        if bSaveSmall:
            sCommand += "_SaveTextures=No "
        sCommand += "_GeometryOnly=No "
        sCommand += "_SaveNotes={} ".format("Yes" if bSaveNotes else "No")
        sCommand += sOutputPath
    elif sFileName.split('.')[-1].lower() in ('igs', 'iges'):
        sCommand += sOutputPath
        sCommand += " "
        sCommand += "_UnitSystem={} ".format(sc.doc.ModelUnitSystem.ToString())
        sCommand += "_Tolerance={} ".format(sc.doc.ModelAbsoluteTolerance)
        if bUseCurrentOpts:
            sCommand += "_Enter"
    else:
        sCommand += sOutputPath
        sCommand += " "
        if bUseCurrentOpts:
            sCommand += "_Enter"

    if bDebug:
        print(sCommand)
        print("Export skipped in Debug mode.")
        return

    Rhino.RhinoApp.RunScript(sCommand, echo=bEcho)


def main():

    #if sc.doc.Path is None:
    #    print("Document doesn't have a path, so the exported file will be saved to the Desktop.")

    gObjs_Pre = [o.Id for o in
        sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]

    objrefs = getObjects()
    if objrefs is None: return

    rc = setOptions()
    if not rc: return

    sFileType = Opts.listValues['iFileType'][Opts.values['iFileType']]

    sFileName = createFileName(bWithExt=False)

    if Opts.values['bManuallyEditName']:
        res, text = Rhino.UI.Dialogs.ShowEditBox(
            title="Export",
            message="Edit file name",
            defaultText=sFileName,
            multiline=False)
        if not res: return
        sFileName = text

    if not sFileName: return

    sFileName = _addExtension(sFileName)


    exportDocObjects(objrefs, sFileName)


    if gObjs_Pre:
        sc.doc.Objects.Select(objectIds=gObjs_Pre)


if __name__ == '__main__': main()