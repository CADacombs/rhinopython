"""
This script uses UserDictionary of a block definition to
    1. Store duplicates of geometry of objects with option to delete the source DocObjects.
    2. Retrieve geometry and add them to the document with option to delete the
    source geometry in the block definitions.
    ...

    Attributes are ignored stored, only Rhino.Geometry objects.

This script only works with Rhino 7 and up.

_Undo and _Redo do not support UserDictionary operations.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums:
https://discourse.mcneel.com/
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
220207-8: Created.
221014-15: Further development.
240926: Bug fix and command option name changes.
240929-1001: Further development.
250208-10, 18, 0417: WIP: Further development.

TODO:
    Add DeleteFromAll routine?
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List


sOptions = (
    'ReadFromBlock',
    'WriteToBlock',
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

    key = 'bRemoveKeysRead'; keys.append(key)
    values[key] = False
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

    key = 'bLayers'; keys.append(key)
    values[key] = False
    names[key] = 'LayerMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Off', 'On')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bBlockParentLayer'; keys.append(key)
    values[key] = True
    names[key] = 'BlockParentLayer'
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

        #addOption('bRemoveKeysRead')
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
            addOption('bAppend_Not_overwrite')
            addOption('bDeleteInputOnWrite')
        else:
            addOption('bReadAll')
            addOption('bRemoveKeysRead')
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
            addOption('bRemoveKeysRead')
        #addOption('bLayers')
        #if Opts.values['bLayers']:
        #    addOption('bBlockParentLayer')
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


def _auditAllBlockDefsForGeometry():
    sNames_DefsWithGeom = []

    for rdIdef in sc.doc.InstanceDefinitions:
        bHasGeom = False
        for key in rdIdef.UserDictionary.Keys:
            item = rdIdef.UserDictionary[key]
            if hasattr(item, '__iter__'):
                iCt = 0
                for item_sub in item:
                    if isinstance(item_sub, rg.GeometryBase):
                        iCt += 1
                if iCt == 1:
                    print("{}'s {} has a {}.".format(
                        rdIdef.Name, key, item_sub.GetType().Name))
                else:
                    print("{}'s {} has {} Rhino.Geometry objects.".format(
                        rdIdef.Name, key, iCt))
                if iCt > 0:
                    bHasGeom = True
            elif isinstance(item, rg.GeometryBase):
                print("{}'s {} has a {}.".format(
                    rdIdef.Name, key, item.GetType().Name))
                bHasGeom = True

        if bHasGeom: sNames_DefsWithGeom.append(rdIdef.Name)

    return sNames_DefsWithGeom


def _createTransformedGeometry(geom_In, xform):

    geoms_Out = []

    if xform.SimilarityType == rg.TransformSimilarityType.NotSimilarity:
        if isinstance(geom_In, rg.ArcCurve):
            geom_Out = geom_In.ToNurbsCurve()
        elif isinstance(geom_In, (rg.PolyCurve, rg.Brep)):
            geom_Out = geom_In.Duplicate()
            geom_Out.MakeDeformable()
    else:
        geom_Out = geom_In.Duplicate()
    geom_Out.Transform(xform)
    return geom_Out

    geoms_Out.append(geom_Out)

    return List[rg.GeometryBase]([geom_In for geom_In in geoms_Out])


def _store_in_def_of_ref_Geoms(rgIref, sKey, rgObjs_In, bAppend=True, bEcho=True):
    """
    rg input, not ObjRefs.
    """

    bSuccess, xform = rg.Transform.TryGetInverse(rgIref.Xform)

    geoms_Xformed = []

    for geom_In in rgObjs_In:
        geom_Xformed = _createTransformedGeometry(geom_In, xform)
        geoms_Xformed.append(geom_Xformed)

    geoms_Xformed = List[rg.GeometryBase]([geom for geom in geoms_Xformed])

    rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)

    if bAppend:
        if not rdIdef.UserDictionary.ContainsKey(sKey):
            pass
        elif sKey not in _get_keys_for_geometries(rgIref):
            pass
        else:
            geoms_ToStore = List[rg.GeometryBase](
                rdIdef.UserDictionary[sKey])
            geoms_ToStore.AddRange(geoms_Xformed)

            # Alternative using Python list +.
            #geoms_ToStore = List[rg.GeometryBase](
            #    list(rdIdef.UserDictionary[sKey]) + list(geoms_ToStore))

            if rdIdef.UserDictionary.Set(sKey, geoms_ToStore):
                return geoms_Xformed.Count
    # Overwrite.
    if rdIdef.UserDictionary.Set(sKey, geoms_Xformed):
        return geoms_Xformed.Count


def _convert_rdObjs_to_geoms(rdObjs):
    """
    Parameters:
        rdObjs: iterable of ObjRef, RhinoObject, or GUID
    Returns:
        .NET List[rg]
    """
    if not rdObjs:
        return
    if isinstance(rdObjs[0], rd.ObjRef):
        return List[rg.GeometryBase]([_.Geometry() for _ in rdObjs])
    if isinstance(rdObjs[0], rd.RhinoObject):
        return List[rg.GeometryBase]([_.Geometry for _ in rdObjs])
    if isinstance(rdObjs[0], Guid):
        return List[rg.GeometryBase]([sc.doc.Objects.FindId(_).Geometry for _ in rdObjs])
    raise Exception("{} not accepted input.".format(rdObjs[0].GetType()))


def _store_in_def_of_ref_Geoms_of_rdObjs(rgIref, rdObjs_In, sKey, bAppend=True, bEcho=True, bDebug=False):
    """
    Returns: bool for Success
    """

    geoms_ToStore = _convert_rdObjs_to_geoms(rdObjs_In)
    if not geoms_ToStore:
        return

    nStored = _store_in_def_of_ref_Geoms(
        rgIref,
        sKey,
        geoms_ToStore,
        bAppend=bAppend,
        bEcho=bEcho)
    if not nStored:
        if bEcho: print("UserDictionary.Set failed.")
        return

    rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)

    print("{} objects packed into block reference '{}'.".format(
        nStored, rdIdef.Name))

    return bool(nStored)


def _is_key_for_geometry(rdIdef, sKey):
    bHasKey, value = rdIdef.UserDictionary.TryGetValue(sKey)

    if not bHasKey:
        raise Exception(
            "rdIdef.UserDictionary.TryGetValue(key) returned False.")

    if isinstance(value, rg.GeometryBase):
        return True

    if not (hasattr(value, '__iter__') and not isinstance(value, str)):
        return False

    if isinstance(value[0], rg.GeometryBase):
        return True

    return False


def _get_geoms_transformed_to_instance(rgIref, sKey):

    rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)

    if not _is_key_for_geometry(rdIdef, sKey):
        return

    geoms_of_def = rdIdef.UserDictionary[sKey]

    xform = rgIref.Xform

    geoms_Out = []

    if xform.SimilarityType == rg.TransformSimilarityType.NotSimilarity:
        for geom in geoms_of_def:
            if isinstance(geom, rg.ArcCurve):
                geom_Out = geom.ToNurbsCurve()
            elif isinstance(geom, (rg.PolyCurve, rg.Brep)):
                geom_Out = geom.Duplicate()
                geom.MakeDeformable()
            else:
                geom_Out = geom.Duplicate()

            geom_Out.Transform(xform)
            geoms_Out.append(geom_Out)
    else:
        for geom in geoms_of_def:
            geom_Out = geom.Duplicate()
            geom_Out.Transform(xform)
            geoms_Out.append(geom_Out)

    return geoms_Out


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
        if _is_key_for_geometry(rdIdef, sKey):
            sKeys_with_geom.append(sKey)

    return sKeys_with_geom


def readFromBlock():

    if Opts.values['bReadAll']:
        objrefs_Iref = getInput_MultiBlockInsts(bWrite_NotRead=False)
        if objrefs_Iref is None: return
    else:
        objref_Iref = getInput_OneBlockInst(bWrite_NotRead=False)
        if objref_Iref is None: return
        objrefs_Iref = [objref_Iref]

    bReadAll = Opts.values['bReadAll']
    bRemoveKeysRead = Opts.values['bRemoveKeysRead']
    bLayers = Opts.values['bLayers']
    bBlockParentLayer = Opts.values['bBlockParentLayer']
    bLayers = Opts.values['bLayers']
    bBlockParentLayer = Opts.values['bBlockParentLayer']
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
            geoms_ret = _get_geoms_transformed_to_instance(rgIref, sKey)
            if not geoms_ret:
                continue
            nFails = 0
            for geom_ret in geoms_ret:
                gOut = sc.doc.Objects.Add(geom_ret)
                if gOut == gOut.Empty:
                    print("Could not add {}.".format(geom_ret))
                    nFails += 1
            if bRemoveKeysRead:
                if nFails:
                    print("Will not delete key '{}' from block since not all geometry could be added to the document.".format(
                        sKey))
                else:
                    print(rdIdef.UserDictionary.Remove(sKey))

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def writeToBlock():

    objref_Iref = getInput_OneBlockInst(bWrite_NotRead=True)
    if objref_Iref is None: return

    bDeleteInputOnWrite = Opts.values['bDeleteInputOnWrite']
    bAppend_Not_overwrite = Opts.values['bAppend_Not_overwrite']
    bLayers = Opts.values['bLayers']
    bBlockParentLayer = Opts.values['bBlockParentLayer']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    #res, objref_Iref = ri.RhinoGet.GetOneObject(
    #    "Select block instance",
    #    acceptNothing=False,
    #    filter=rd.ObjectType.InstanceReference)
    #if res != Rhino.Commands.Result.Success: return

    sc.doc.Objects.UnselectAll()

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

    keys_with_geoms = _get_keys_for_geometries(rgIref)

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


    bSuccess = _store_in_def_of_ref_Geoms_of_rdObjs(
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


def listKVPsWithGeom():
    res, objrefs_Iref = ri.RhinoGet.GetMultipleObjects(
        "Select block instances <All definitions>",
        acceptNothing=True,
        filter=rd.ObjectType.InstanceReference)
    if res != Rhino.Commands.Result.Success: return

    if objrefs_Iref is None:
        rdIdefs = list(sc.doc.InstanceDefinitions)
        #auditAllBlockDefsForGeometry_OLD()
        #return
    else:
        rdIdefs = []

        for objref in objrefs_Iref:
            rgIref = objref.Geometry()
            rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)
            if rdIdef not in rdIdefs:
                rdIdefs.append(rdIdef)

    iDictsWithGeomValues = 0

    for rdIdef in rdIdefs:

        sKeys = _get_keys_for_geometries(rdIdef)

        if sKeys is None:
            if len(rdIdefs) == 1:
                print("There is no dictionary in {}.".format(rdIdef.Name))
            continue

        if len(sKeys) == 0:
            if len(rdIdefs) == 1:
                print("There is no geometry in {}'s dictionary.".format(rdIdef.Name))
            continue 

        iDictsWithGeomValues += 1

        if iDictsWithGeomValues == 1:
            print('-'*80)
            print('-'*80)
            s = "BlockDefName:Key:GeomCt for the following:"
            print(s, '-'*len(s), sep='\n')

        for sKey in sKeys:
            print("'{}':'{}':{}".format(rdIdef.Name, sKey, rdIdef.UserDictionary[sKey].Count))

    if iDictsWithGeomValues == 0:
        print("There are no dictionary values of geometry in any of the {} definitions.".format(len(rdIdefs)))

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def removeKeys():
    res, objref_Iref = ri.RhinoGet.GetOneObject(
        "Select block instance",
        acceptNothing=False,
        filter=rd.ObjectType.InstanceReference)
    if res != Rhino.Commands.Result.Success: return

    rgIref = objref_Iref.Geometry()
    rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)

    keys_with_geoms = _get_keys_for_geometries(rgIref)
    if not keys_with_geoms:
        print("No keys for geometry.")
        return

    sKeys = rs.MultiListBox(
        items=keys_with_geoms,
        message="Pick keys to delete.",
        title="Delete Keys for Values with Geometry",
        defaults=None)
    if sKeys is None:
        return

    ssKeys_deleted = []
    for sKey in sKeys:
        if rdIdef.UserDictionary.Remove(sKey):
            ssKeys_deleted.append(sKey)
        else:
            print("Could not delete key, {}".format(sKey))

    print("Deleted {} keys for values with geometry.".format(len(ssKeys_deleted)))

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def main():

    if sc.doc.InstanceDefinitions.ActiveCount == 0:
        print("No blocks in document. End of script.")
        return

    bDebug = Opts.values['bDebug']

    if bDebug:
        _auditAllBlockDefsForGeometry()

    rc = getOption()
    if rc is None: return
    sOption = rc

    bRemoveKeysRead = Opts.values['bRemoveKeysRead']
    bDeleteInputOnWrite = Opts.values['bDeleteInputOnWrite']
    #bLayers = Opts.values['bLayers']
    #bBlockParentLayer = Opts.values['bBlockParentLayer']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if sOption == 'ReadFromBlock':
        readFromBlock()
        return

    if sOption == 'DumpFromBlocks':
        res, objrefs_Iref = ri.RhinoGet.GetMultipleObjects(
            "Select block instances",
            acceptNothing=False,
            filter=rd.ObjectType.InstanceReference)
        if res != Rhino.Commands.Result.Success: return

        for objref_Iref in objrefs_Iref:
            rgIref = objref_Iref.Geometry()
            rdIdef = sc.doc.InstanceDefinitions.FindId(rgIref.ParentIdefId)

            keys_with_geoms = _get_keys_for_geometries(rgIref)

            for sKey in keys_with_geoms:
                geoms_ret = _get_geoms_transformed_to_instance(rgIref, sKey)
                if not geoms_ret:
                    continue
                nFails = 0
                for geom_ret in geoms_ret:
                    gOut = sc.doc.Objects.Add(geom_ret)
                    if gOut == gOut.Empty:
                        print("Could not add {}.".format(geom_ret))
                        nFails += 1

        sc.doc.Objects.UnselectAll()
        sc.doc.Views.Redraw()

        return

    if sOption == 'WriteToBlock':
        writeToBlock()
        return

    if sOption == 'ListKVPsWithGeom':
        listKVPsWithGeom()
        return

    if sOption == 'RemoveKeys':
        removeKeys()



    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


if __name__ == '__main__': main()