"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
230329-30: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Input as ri
import scriptcontext as sc

from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'bRemoveColor'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iR'; keys.append(key)
    values[key] = 255
    riOpts[key] = ri.Custom.OptionInteger(values[key], 0, 255)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iG'; keys.append(key)
    values[key] = 0
    riOpts[key] = ri.Custom.OptionInteger(values[key], 0, 255)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iB'; keys.append(key)
    values[key] = 0
    riOpts[key] = ri.Custom.OptionInteger(values[key], 0, 255)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseBrepColor'; keys.append(key)
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
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
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
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """

    goB = ri.Custom.GetObject()

    goB.SetCommandPrompt("Select breps to set all their faces")

    goB.GeometryFilter = rd.ObjectType.Brep

    goB.AlreadySelectedObjectSelect = True
    goB.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    goB.EnableClearObjectsOnEntry(False) # Keep objects in goB on repeats of While loop.
    goB.EnableUnselectObjectsOnExit(False) # False to keep objects selected when picking an option.

    bPreselectedObjsChecked = False

    sQuickPicks = []
    rgbs_QP = {}
    idxs_Opts_QP = {}

    key = 'Red'
    sQuickPicks.append(key)
    rgbs_QP[key] = 255,0,0
    key = 'Green'
    sQuickPicks.append(key)
    rgbs_QP[key] = 0,255,0
    key = 'Blue'
    sQuickPicks.append(key)
    rgbs_QP[key] = 0,0,255

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(goB, key)

    while True:
        goB.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('bRemoveColor')
        if not Opts.values['bRemoveColor']:
            addOption('bUseBrepColor')
            if not Opts.values['bUseBrepColor']:
                addOption('iR')
                addOption('iG')
                addOption('iB')
                for sOpt in sQuickPicks:
                    idxs_Opts_QP[sOpt] = idxs_Opts[sOpt] = goB.AddOption(sOpt)
                idxs_Opts['Dialog'] = goB.AddOption('Dialog')
                idxs_Opts['MatchFace'] = goB.AddOption('MatchFace')
        addOption('bEcho')
        addOption('bDebug')

        res = goB.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            goB.Dispose()
            return

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to goB.GetMultiple is considered.
        if not bPreselectedObjsChecked and goB.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            goB.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Object:
            objrefs = goB.Objects()
            goB.Dispose()
            return objrefs

        # An option was selected.
        if goB.Option().Index in idxs_Opts_QP.values():
            key_QP = list(idxs_Opts_QP)[list(idxs_Opts_QP.values()).index(goB.Option().Index)]
            key = 'iR'; Opts.riOpts[key].CurrentValue = rgbs_QP[key_QP][0]; Opts.setValue(key)
            key = 'iG'; Opts.riOpts[key].CurrentValue = rgbs_QP[key_QP][1]; Opts.setValue(key)
            key = 'iB'; Opts.riOpts[key].CurrentValue = rgbs_QP[key_QP][2]; Opts.setValue(key)
            continue

        if 'Dialog' in idxs_Opts and goB.Option().Index == idxs_Opts['Dialog']:
            color = Color.FromArgb(
                red=Opts.values['iR'],
                green=Opts.values['iG'],
                blue=Opts.values['iB'])
            bSuccess, color = Rhino.UI.Dialogs.ShowColorDialog(color)
            if bSuccess:
                key = 'iR'; Opts.riOpts[key].CurrentValue = color.R; Opts.setValue(key)
                key = 'iG'; Opts.riOpts[key].CurrentValue = color.G; Opts.setValue(key)
                key = 'iB'; Opts.riOpts[key].CurrentValue = color.B; Opts.setValue(key)
            continue

        if 'MatchFace' in idxs_Opts and goB.Option().Index == idxs_Opts['MatchFace']:
            goF = ri.Custom.GetObject()
            goF.SetCommandPrompt("Pick face to obtain its color")
            goF.GeometryFilter = rd.ObjectType.Surface
            goF.DeselectAllBeforePostSelect = False # False to keep goB selections.
            goF.EnableHighlight(False)
            goB.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True) # False to not change brep of selected face to match.
            res = goF.Get()

            if res == ri.GetResult.Cancel:
                goF.Dispose()
                goB.Dispose()
                return

            if res == ri.GetResult.Object:
                objref = goF.Object(0)
                goF.Dispose()

                color = objref.Face().PerFaceColor
                key = 'iR'; Opts.riOpts[key].CurrentValue = color.R; Opts.setValue(key)
                key = 'iG'; Opts.riOpts[key].CurrentValue = color.G; Opts.setValue(key)
                key = 'iB'; Opts.riOpts[key].CurrentValue = color.B; Opts.setValue(key)

            continue


        for key in idxs_Opts:
            if goB.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, goB.Option().CurrentListOptionIndex)
                break


def _sortInput(objrefs):
    gBs = []
    rdIdefs = []

    for objref in objrefs:
        rdObj = objref.Object()
        if rdObj.ObjectType == rd.ObjectType.InstanceReference:
            rdIdef = rdObj.InstanceDefinition
            rdB = objref.InstanceDefinitionPart()
        else:
            rdB = rdObj
            rdIdef = None
        gB = rdB.Id
        if gB not in gBs:
            gBs.append(gB)
            rdIdefs.append(rdIdef)

    return gBs, rdIdefs


def _replaceInDef(rdIdef, gB, rgB):

    geoms = []
    attrs = []

    for rdObj in rdIdef.GetObjects():
        if rdObj.Id == gB:
            geoms.append(rgB)
        else:
            geoms.append(rdObj.Geometry)
        attrs.append(rdObj.Attributes)

    if not sc.doc.InstanceDefinitions.ModifyGeometry(rdIdef.Index, geoms, attrs):
        print("Could not modify block instance.")


def main():

    objrefs = getInput()
    if objrefs is None: return

    gBs, rdIdefs = _sortInput(objrefs)

    if Opts.values['bRemoveColor']:
        color = Color.Empty
    else:
        color = Color.FromArgb(
            red=Opts.values['iR'],
            green=Opts.values['iG'],
            blue=Opts.values['iB'])

    bUseBrepColor = Opts.values['bUseBrepColor'] and not Opts.values['bRemoveColor']

    iCt_Bs_WithMods = 0
    iCt_Mods_AllBs = 0

    for gB, rdIdef in zip(gBs, rdIdefs):
        rdB = sc.doc.Objects.FindId(gB)
        rgB = rdB.BrepGeometry
        iCt_Mods_ThisB = 0
        if bUseBrepColor:
            if rdB.Attributes.ColorSource == rd.ObjectColorSource.ColorFromLayer:
                layer = sc.doc.Layers.FindIndex(rdB.Attributes.LayerIndex)
                color = layer.Color
            elif rdB.Attributes.ColorSource == rd.ObjectColorSource.ColorFromObject:
                color = rdB.Attributes.ObjectColor
            else:
                print("ColorSource {} not supported (yet?).".format(
                    rdB.Attributes.ColorSource),
                      "Face colors were not changed.")
                continue
        for idxF in xrange(rgB.Faces.Count):
            if rgB.Faces[idxF].PerFaceColor != color:
                rgB.Faces[idxF].PerFaceColor = color
                if rgB.Faces[idxF].PerFaceColor == color:
                    iCt_Mods_ThisB += 1
        if rdIdef is None:
            if sc.doc.Objects.Replace(objectId=gB, brep=rgB):
                iCt_Bs_WithMods += 1
                iCt_Mods_AllBs += iCt_Mods_ThisB
        else:
            _replaceInDef(rdIdef, gB, rgB)

    if iCt_Mods_AllBs == 0:
        print("Not including instances, no face colors were modified.")
    else:
        print("Not including instances, {} face colors were modified in {} breps.".format(
            iCt_Mods_AllBs, iCt_Bs_WithMods))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()