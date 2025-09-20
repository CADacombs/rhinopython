"""
This script sets a RhinoObject's Attribute's color or color source, e.g., ByLayer,
to that of another object.
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250919: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Drawing import Color


_result_GUIDs = {
    'Modified': [],
    'AlreadyAtTarget': [],
    'FailedFilter': [],
    'FailedModification': [],
    }

_sSetResult = {
    'ObjectColorSource': None,
    'ObjectColor': None,
    }


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAcceptByLayer'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAcceptByParent'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAcceptByMaterial'; keys.append(key)
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

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList
            return

        print("Invalid key?")


def getInput_toModify():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select objects to modify their color")

    #go.GeometryFilter = rd.ObjectType.Brep

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('bAcceptByLayer')
        addOption('bAcceptByParent')
        addOption('bAcceptByMaterial')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)


        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Source():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select source object")

    #go.GeometryFilter = rd.ObjectType.Brep

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('bAcceptByLayer')
        addOption('bAcceptByParent')
        addOption('bAcceptByMaterial')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return objref

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _tryFindRhinoColorName(color):
    for rh_UI_NamedColor in Rhino.UI.NamedColorList.Default:
        if rh_UI_NamedColor.Color == color:
            return rh_UI_NamedColor.Name
    return color.ToString()


def _set_Color(rdObjs_toModify, rdObj_Source):
    color_Source = rdObj_Source.Attributes.ObjectColor

    for rdO in rdObjs_toModify:

        if (
            (rdO.Attributes.ObjectColor == color_Source) and
            (rdO.Attributes.ColorSource == colorSource_Source)
        ):
            _result_GUIDs['AlreadyAtTarget'].append(rdO.Id)
            continue

        if rdO.Attributes.ObjectColor != color_Source:
            rdO.Attributes.ObjectColor = color_Source

        if rdO.Attributes.ColorSource != rd.ObjectColorSource.ColorFromObject:
            rdO.Attributes.ColorSource = rd.ObjectColorSource.ColorFromObject

        if not rdO.CommitChanges():
            _result_GUIDs['FailedModification'].append(rdO.Id)
            continue

        if (
            (rdO.Attributes.ObjectColor == color_Source) and
            (rdO.Attributes.ColorSource == rd.ObjectColorSource.ColorFromObject)
        ):
            _result_GUIDs['Modified'].append(rdO.Id)
        else:
            _result_GUIDs['FailedModification'].append(rdO.Id)
            continue

    _sSetResult['ObjectColor'] = color_Source


def _set_ColorSource(rdObjs_toModify, rdObj_Source):

    colorSource_Source = rdObj_Source.Attributes.ColorSource

    for rdO in rdObjs_toModify:
        if rdO.Attributes.ColorSource == colorSource_Source:
            _result_GUIDs['AlreadyAtTarget'].append(rdO.Id)
            continue

        rdO.Attributes.ColorSource = colorSource_Source

        if not rdO.CommitChanges():
            _result_GUIDs['FailedModification'].append(rdO.Id)
            continue

        if rdO.Attributes.ColorSource == rd.ObjectColorSource.ColorFromLayer:
            _result_GUIDs['Modified'].append(rdO.Id)
        else:
            _result_GUIDs['FailedModification'].append(rdO.Id)
            continue

    _sSetResult['ObjectColorSource'] = colorSource_Source


def _process_ByLayer_ColorSource(rdObjs_toModify, rdObj_Source):
    color_Layer = sc.doc.Layers[rdObj_Source.Attributes.LayerIndex].Color
    sColor = _tryFindRhinoColorName(color_Layer)

    rc, bLayerColor_notByLayer = ri.RhinoGet.GetBool(
        prompt="Source object is ByLayer with layer color {}. Set object color to".format(sColor),
        acceptNothing=False,
        offPrompt="ByLayer",
        onPrompt="LayerColor",
        boolValue=True)
    if rc != Rhino.Commands.Result.Success: return

    if bLayerColor_notByLayer:
        for rdO in rdObjs_toModify:
            if (
                (rdO.Attributes.ColorSource == rd.ObjectColorSource.ColorFromObject) and
                (rdO.Attributes.ObjectColor == color_Layer)
            ):
                _result_GUIDs['AlreadyAtTarget'].append(rdO.Id)
                continue

            if rdO.Attributes.ColorSource != rd.ObjectColorSource.ColorFromObject:
                rdO.Attributes.ColorSource = rd.ObjectColorSource.ColorFromObject

            if rdO.Attributes.ObjectColor != color_Layer:
                rdO.Attributes.ObjectColor = color_Layer

            if not rdO.CommitChanges():
                _result_GUIDs['FailedModification'].append(rdO.Id)
                continue

            if (
                (rdO.Attributes.ColorSource == rd.ObjectColorSource.ColorFromObject) and
                (rdO.Attributes.ObjectColor == color_Layer)
            ):
                _result_GUIDs['Modified'].append(rdO.Id)
            else:
                _result_GUIDs['FailedModification'].append(rdO.Id)
                continue

        _sSetResult['ObjectColor'] = color_Layer

    else:
        _set_ColorSource(rdObjs_toModify, rdObj_Source)


def main():

    rdObjs_toModify = list(
        sc.doc.Objects.GetSelectedObjects(
            includeLights=False,
            includeGrips=False))

    if not rdObjs_toModify:
        #objrefs_toModify = getInput_toModify()
        #if objrefs_toModify is None: return

        rc, objrefs_toModify = ri.RhinoGet.GetMultipleObjects(
            prompt="Select objects to modify their color",
            acceptNothing=False,
            filter=rd.ObjectType.AnyObject)
        if rc != Rhino.Commands.Result.Success: return

        rdObjs_toModify = [o.Object() for o in objrefs_toModify]

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    objref_Source = getInput_Source()
    if objref_Source is None: return

    #rc, objref_Source = ri.RhinoGet.GetOneObject(
    #    prompt="Select source object",
    #    acceptNothing=False,
    #    filter=rd.ObjectType.AnyObject)
    #if rc != Rhino.Commands.Result.Success: return

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    bAcceptByLayer = Opts.values['bAcceptByLayer']
    bAcceptByParent = Opts.values['bAcceptByParent']
    bAcceptByMaterial = Opts.values['bAcceptByMaterial']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    rdObj_Source = objref_Source.Object()
    colorSource_Source = rdObj_Source.Attributes.ColorSource

    if not bAcceptByLayer:
        if colorSource_Source == rd.ObjectColorSource.ColorFromLayer:
            print("Source object color is ByLayer, so ending script.")
            return

    if not bAcceptByParent:
        if colorSource_Source == rd.ObjectColorSource.ColorFromParent:
            print("Source object color is ByParent, so ending script..")
            return

    if not bAcceptByMaterial:
        if colorSource_Source == rd.ObjectColorSource.ColorFromMaterial:
            print("Source object color is ByMaterial, so ending script.")
            return

    if colorSource_Source == rd.ObjectColorSource.ColorFromLayer:
        _process_ByLayer_ColorSource(rdObjs_toModify, rdObj_Source)
    elif colorSource_Source != rd.ObjectColorSource.ColorFromObject:
        _set_ColorSource(rdObjs_toModify, rdObj_Source)
    else:
        _set_Color(rdObjs_toModify, rdObj_Source)

    if _sSetResult['ObjectColorSource']:
        print("Resultant object color source: {}".format(_sSetResult['ObjectColor']))
    if _sSetResult['ObjectColor']:
        print("Resultant object color: {}".format(
            _tryFindRhinoColorName(_sSetResult['ObjectColor'])))


    if len(_result_GUIDs['AlreadyAtTarget']) == len(rdObjs_toModify):
        print("All {} objects are already at target color.".format(len(rdObjs_toModify)))
        return

    for key in _result_GUIDs:
        if _result_GUIDs[key]:
            print("{} {}".format(len(_result_GUIDs[key]), key))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()