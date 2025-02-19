"""
Inspired by CadQuery, this script manages build history in block defintions.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums:
https://discourse.mcneel.com/
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250214-: WIP: Created.

TODO:
    Assign random color to each layer that is created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import spb_Block_UserDict_Geom


sOptions = (
    'ReadGeomFromKVPs',
    'WriteGeomToKVPs',
    'CreateLayersForBlocks',
    'ListKVPsWithGeom',
    'RemoveKeys',
    )


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bReadAll'; keys.append(key)
    values[key] = False
    names[key] = 'ReadMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'SelectKeys', 'AllKeys')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bWriteAllNormalGeomPerLayers'; keys.append(key)
    values[key] = True
    #names[key] = 'WriteMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAppend_Not_overwrite'; keys.append(key)
    values[key] = False
    names[key] = 'WriteMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Overwrite', 'Append')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteInputOnWrite'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteLayersOnWrite'; keys.append(key)
    values[key] = True
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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getOption():

    stickyKey_DefaultOpt = 'DefaultOpt({})'.format(__file__)

    if sc.sticky.has_key(stickyKey_DefaultOpt):
        idxOpt_Default = sc.sticky[stickyKey_DefaultOpt]
        if idxOpt_Default >= len(sOptions):
            sc.sticky[stickyKey_DefaultOpt] = idxOpt_Default = 0
    else:
        idxOpt_Default = 0

    go = ri.Custom.GetOption()
    go.SetCommandPrompt("Choose option")

    go.SetCommandPromptDefault(defaultValue=sOptions[idxOpt_Default])
    go.AcceptNothing(True)

    # for idxOpt of 1 - n, not 0.

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    #    for s in listPropOpts: gs.AddOption(s)
    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        for sDirection in sOptions:
            go.AddOption(sDirection)

        #addOption('bDeleteInputOnWrite')
        #addOption('bLayers')
        #if Opts.values['bLayers']:
        #    addOption('bBlockParentLayer')
        #addOption('bEcho')
        #addOption('bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            idxOpt = idxOpt_Default
            return sOptions[idxOpt]


        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break
        else:
            idxOpt = go.Option().Index - 1 # '- 1' because go.Option().Index is base 1.
            go.Dispose()
            sc.sticky[stickyKey_DefaultOpt] = idxOpt
            return sOptions[idxOpt]


def getInput_OneBlockInst(bWrite_NotRead):

    go = ri.Custom.GetObject()

    go.SetCommandPrompt('Select block instance')

    go.GeometryFilter = rd.ObjectType.InstanceReference

    def addOption(key):
        Opts.addOption(go, key)
        keys_in_order_added.append(key)

    while True:
        go.ClearCommandOptions()
        keys_in_order_added = [None] # None is just a placeholder since option indices are base 1.

        if bWrite_NotRead:
            addOption('bWriteAllNormalGeomPerLayers')
            addOption('bAppend_Not_overwrite')
            addOption('bDeleteInputOnWrite')
            if Opts.values['bDeleteInputOnWrite']:
                addOption('bDeleteLayersOnWrite')
        else:
            addOption('bReadAll')
        #addOption('bLayers')
        #if Opts.values['bLayers']:
        #    addOption('bBlockParentLayer')
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

        Opts.setValue(keys_in_order_added[go.Option().Index])

        if Opts.values['bReadAll']:
            return getInput_MultiBlockInsts(bWrite_NotRead)


def getInput_MultiBlockInsts(bWrite_NotRead):

    go = ri.Custom.GetObject()

    go.SetCommandPrompt('Select block instances')

    go.GeometryFilter = rd.ObjectType.InstanceReference

    def addOption(key):
        Opts.addOption(go, key)
        keys_in_order_added.append(key)

    while True:
        go.ClearCommandOptions()
        keys_in_order_added = [None] # None is just a placeholder since option indices are base 1.

        if bWrite_NotRead:
            addOption('bAppend_Not_overwrite')
            addOption('bDeleteInputOnWrite')
        else:
            addOption('bReadAll')
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

        Opts.setValue(keys_in_order_added[go.Option().Index])

        if Opts.values['bReadAll']:
            return getInput_OneBlockInst(bWrite_NotRead)


def create_geometry_in_block(sBlock):
    """
    Parameters:
    Returns:
    """
    sEval = "sBlock"; print(sEval,'=',eval(sEval))

    rdIdef = sc.doc.InstanceDefinitions.Find(sBlock)
    if rdIdef is None:
        raise Exception("Definition for block '{}' doesn't exist.".format(sBlock))

#     # print(type(rdIdef))

#     # crvs_perimeter = 1/0
#     keys_with_geoms = _get_keys_for_geometries(rdIdef)
#     if not keys_with_geoms:
#         print("No keys for geometry.")
#         return

#     geoms = rdIdef.UserDictionary['sketch1']


def _get_keys_for_geometries(rgIref_or_rdIdef):
    """
    Returns one of the following:
        None for no dictionary
        Empty list for no keys with values containing geometry or iterable of geometry
        Non-empty list for keys with values containing geometry or iterable of geometry
    """

    if isinstance(rgIref_or_rdIdef, rd.InstanceDefinition):
        rdIdef = rgIref_or_rdIdef
    elif isinstance(rgIref_or_rdIdef, rg.InstanceReferenceGeometry):
        rgIref = rgIref_or_rdIdef
        rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)
    else:
        import sys
        raise Exception("{} sent to {}.".format(
            rgIref_or_rdIdef.GetType().Name,
            sys._getframe().f_code.co_name))

    if rdIdef.UserDictionary.Count == 0:
        return

    sKeys_with_geom = []
    for sKey in rdIdef.UserDictionary.Keys:
        if spb_Block_UserDict_Geom._is_key_for_geometry(rdIdef, sKey):
            sKeys_with_geom.append(sKey)

    return sKeys_with_geom


def prepareLayer(sLayerName, sLayerName_Parent=None):
    """
    Returns: str of layer path.
    """

    if sLayerName_Parent:
        sLayerPath = "{}::{}".format(sLayerName_Parent, sLayerName)
    else:
        sLayerPath = sLayerName

    if not rs.IsLayer(sLayerPath):
        return rs.AddLayer(
            name=sLayerPath,
            color=None,
            visible=True,
            locked=False,
            parent=None)

    if not rs.IsLayerVisible(sLayerPath):
        rs.LayerVisible(sLayerPath)
    if rs.IsLayerLocked(sLayerPath):
        rs.LayerLocked(sLayerPath, locked=False)

    return sLayerPath


def createLayersForBlocks():
    res, objrefs_Iref = ri.RhinoGet.GetMultipleObjects(
        "Select instances <All>",
        acceptNothing=True,
        filter=rd.ObjectType.InstanceReference)
    if res != Rhino.Commands.Result.Success: return

    if objrefs_Iref is None:
        rdIdefs = list(sc.doc.InstanceDefinitions)
    else:
        rdIdefs = []
        for objref in objrefs_Iref:
            rgIref = objref.Geometry()
            rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)
            if rdIdef not in rdIdefs:
                rdIdefs.append(rdIdef)

    iDictsWithGeomValues = 0

    for rdIdef in rdIdefs:
        #prepareLayer(rdIdef.Name)

        #sKeys = _get_keys_for_geometries(rdIdef)
        sKeys = spb_Block_UserDict_Geom._get_keys_for_geometries(rdIdef)

        if sKeys is None:
            if len(rdIdefs) == 1:
                print("There is no dictionary in {}.".format(rdIdef.Name))
            continue

        if len(sKeys) == 0:
            if len(rdIdefs) == 1:
                print("There is no geometry in {}'s dictionary.".format(rdIdef.Name))
            continue 

        iDictsWithGeomValues += 1

        for sKey in sKeys:
            if not rd.ModelComponent.IsValidComponentName(sKey):
                print("{} is not a valid layer name and will be skipped".format(sKey))
                continue
            prepareLayer(sKey, rdIdef.Name)

    if iDictsWithGeomValues == 0:
        print("There are no dictionary values of geometry in any of the {} definitions.".format(len(rdIdefs)))


def getAllNormalObjectsOnLayer(sLayerPath):
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = (
        rd.ObjectType.Brep |
        rd.ObjectType.Curve |
        rd.ObjectType.Extrusion |
        rd.ObjectType.None
        )
    idxLayer = sc.doc.Layers.FindByFullPath(
        layerPath=sLayerPath,
        ignoreDeletedLayers=True)
    #sEval = "idxLayer"; print(sEval,'=',eval(sEval))
    rdObjs_Out = []
    for rdObj in sc.doc.Objects.GetObjectList(oes):
        if rdObj.Attributes.LayerIndex == idxLayer:
            rdObjs_Out.append(rdObj)
    return rdObjs_Out


def readGeomFromKVPs():

    if Opts.values['bReadAll']:
        objrefs_Iref = getInput_MultiBlockInsts(bWrite_NotRead=False)
        if objrefs_Iref is None: return
    else:
        objref_Iref = getInput_OneBlockInst(bWrite_NotRead=False)
        if objref_Iref is None: return
        objrefs_Iref = [objref_Iref]

    bReadAll = Opts.values['bReadAll']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    #res, objref_Iref = ri.RhinoGet.GetOneObject(
    #    "Select block instance",
    #    acceptNothing=False,
    #    filter=rd.ObjectType.InstanceReference)
    #if res != Rhino.Commands.Result.Success: return

    #sEval = "sc.doc.Objects.SelectedObjectsExist(objectType=rd.ObjectType.AnyObject, checkSubObjects=False)"; print(sEval,'=',eval(sEval))
    sc.doc.Objects.UnselectAll()

    for objref_Iref in objrefs_Iref:
        rgIref = objref_Iref.Geometry()
        rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)

        keys_with_geoms = _get_keys_for_geometries(rgIref)
        if not keys_with_geoms:
            print("No keys for geometry.")
            return

        if bReadAll:
            sKeys = keys_with_geoms
        elif len(keys_with_geoms) == 1:
            if bEcho: print("Only 1 key/value pair with geometry.")
            sKeys = keys_with_geoms
        else:
            sKeys = rs.MultiListBox(
                items=keys_with_geoms,
                message="Pick keys to access their geometries.",
                title="Get Geometry per Key",
                defaults=None)
            if sKeys is None:
                return

        for sKey in sKeys:
            geoms_ret = spb_Block_UserDict_Geom._get_geoms_transformed_to_instance(rgIref, sKey)
            if not geoms_ret:
                continue

            if not rd.ModelComponent.IsValidComponentName(sKey):
                print("{} is not a valid layer name and will be skipped".format(sKey))
                continue
            if '::' in sKey:
                print("{} is not a valid layer name and will be skipped".format(sKey))
                continue

            sLayerPath = prepareLayer(
                sLayerName=sKey,
                sLayerName_Parent=rdIdef.Name)
            attr = rd.ObjectAttributes()
            idxLayer = sc.doc.Layers.FindByFullPath(
                layerPath=sLayerPath,
                ignoreDeletedLayers=True)
            attr.LayerIndex = idxLayer

            nFails = 0
            for geom_ret in geoms_ret:
                gOut = sc.doc.Objects.Add(geom_ret, attributes=attr)
                if gOut == gOut.Empty:
                    print("Could not add {}.".format(geom_ret))
                    nFails += 1

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def writeGeomToKVPs():

    objref_Iref = getInput_OneBlockInst(bWrite_NotRead=True)
    if objref_Iref is None: return

    bWriteAllNormalGeomPerLayers = Opts.values['bWriteAllNormalGeomPerLayers']
    bAppend_Not_overwrite = Opts.values['bAppend_Not_overwrite']
    bDeleteInputOnWrite = Opts.values['bDeleteInputOnWrite']
    bDeleteLayersOnWrite = Opts.values['bDeleteLayersOnWrite']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    #res, objref_Iref = ri.RhinoGet.GetOneObject(
    #    "Select block instance",
    #    acceptNothing=False,
    #    filter=rd.ObjectType.InstanceReference)
    #if res != Rhino.Commands.Result.Success: return

    sc.doc.Objects.UnselectAll()

    if bWriteAllNormalGeomPerLayers:
        rgIref = objref_Iref.Geometry()
        rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)
        if not rs.IsLayer(rdIdef.Name):
            print("Layer '{}' does not exist, so no geometry will be written to its dictionary.".format(rdIdef.Name))
            return
        if rs.LayerChildCount(rdIdef.Name) == 0:
            print("Layer '{}' has no children, so no geometry will be written to its dictionary.".format(rdIdef.Name))
            return
        for sChildLayerPath in rs.LayerChildren(rdIdef.Name):
            rdObjs_OnLayer = getAllNormalObjectsOnLayer(sChildLayerPath)

            bSuccess = spb_Block_UserDict_Geom._store_in_def_of_ref_Geoms_of_rdObjs(
                rgIref,
                rdObjs_OnLayer,
                sChildLayerPath.split('::')[-1],
                bAppend=bAppend_Not_overwrite,
                bEcho=bEcho,
                bDebug=bDebug)

            if bSuccess and bDeleteInputOnWrite:
                [sc.doc.Objects.Delete(o, quiet=False) for o in rdObjs_OnLayer]
                #rs.DeleteObjects([o.ObjectId for o in objrefs_Geom_In])
                #rs.DeleteLayer(key)

    return


    res, objrefs_Geom_In = ri.RhinoGet.GetMultipleObjects(
        "Select objects to save in block",
        acceptNothing=False,
        filter=rd.ObjectType.AnyObject)
    if res != Rhino.Commands.Result.Success: return

    objrefs_Geom_In = list(objrefs_Geom_In) # from Array.

    # If the container block is in selection, remove it.
    for objref_ToPack in objrefs_Geom_In:
        if objref_ToPack.ObjectId == objref_Iref.ObjectId:
            objrefs_Geom_In.remove(objref_ToPack)
            if len(objrefs_Geom_In) == 0:
                print("No objects.")
                return
            print("Subject instance was removed from selection.")
            break

    #if 'New' in sOption:
    #    sKey = rs.StringBox(message="Enter key", default_value=None, title="New Storage Set")
    #    if sKey is None: return
    #    rgIref = objref_Iref.Geometry()
    #if sOption in ('a', 'ad'):
    #    # WIP: Get keys.
    #    rgIref = objref_Iref.Geometry()
    #    rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)
    #    for key in rdIdef.UserDictionary.Keys:
    #        pass

    rgIref = objref_Iref.Geometry()

    #keys_with_geoms = _get_keys_for_geometries(rgIref)
    keys_with_geoms = spb_Block_UserDict_Geom._get_keys_for_geometries(rgIref)

    if not keys_with_geoms:
        sKey = rs.StringBox(message="Enter key", default_value=None, title="New Storage Set")
        if sKey is None: return
    else:
        sKey = rs.ListBox(
            items=keys_with_geoms,
            message="Pick key to write to its value or <Cancel> to type a new key",
            title="Write",
            default=None)
        sEval = "sKey"; print(sEval,'=',eval(sEval))
        if sKey is None:
            sKey = rs.StringBox(message="Enter key", default_value=None, title="New Storage Set")
            if sKey is None: return


    bSuccess = _store_in_def_of_ref_Geoms_of_objrefs(
        rgIref,
        objrefs_Geom_In,
        sKey,
        bAppend=bAppend_Not_overwrite,
        bEcho=bEcho,
        bDebug=bDebug)

    if bSuccess and bDeleteInputOnWrite:
        [sc.doc.Objects.Delete(o, quiet=False) for o in objrefs_Geom_In]
        #rs.DeleteObjects([o.ObjectId for o in objrefs_Geom_In])
        #rs.DeleteLayer(key)

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def main():

    if sc.doc.InstanceDefinitions.ActiveCount == 0:
        print("No blocks in document. Create blocks/instance(s) before running this script.")
        return

    bDebug = Opts.values['bDebug']

    if bDebug:
        spb_Block_UserDict_Geom._auditAllBlockDefsForGeometry()

    rc = getOption()
    if rc is None: return
    sOption = rc

    if sOption == 'ReadGeomFromKVPs':
        readGeomFromKVPs()
        return

    if sOption == 'WriteGeomToKVPs':
        writeGeomToKVPs()
        return

    if sOption == 'CreateLayersForBlocks':
        createLayersForBlocks()
        return


    return

    res, objref = ri.RhinoGet.GetOneObject(
        prompt="Select instance",
        acceptNothing=False,
        filter=rd.ObjectType.InstanceReference)




    sBlock = 'box'

    create_geometry_in_block(sBlock)



    return

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select curves",
        acceptNothing=False,
        filter=rd.ObjectType.Curve)
    if res != Rhino.Commands.Result.Success: return

    gObj = rs.GetObject(
        "Select object",
        filter=0,
        preselect=True,
        select=True)

    gObjs = rs.GetObjects(
        "Select objects",
        filter=0,
        preselect=True,
        select=True)


    if rs.SelectedObjects():
        gObjs = rs.GetObjects("Select objects", filter=0,
                preselect=True, select=True)
    else:
        gObjs = rs.GetObjects("Select objects", filter=0)
    if not gObjs: return


    brep = rs.coercebrep(gObj)

    bDebug = True
    if bDebug: sEval="brep"; print("{}: {}".format(sEval, eval(sEval)))

if __name__ == '__main__': main()