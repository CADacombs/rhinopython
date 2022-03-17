"""
180521-25: Created, starting with deleteCrvsInAllNonRefBlocks.py.
190924-25: Bug fixes for better handling of blocks that may contain instances
        to delete but are not on layers to delete.
        Improved feedback.
        Refactored.
"""

import Rhino
import rhinoscriptsyntax as rs
import scriptcontext as sc


def getLayersToDelete(bEcho=True):
    """
    """
    
    layers_ToDel_Picked = rs.GetLayers(
            "Select layers to delete (Reference layers will be ignored.)")
    if not layers_ToDel_Picked: return
    
    layers_ToDel_Picked_NoRefs = []
    for layer in layers_ToDel_Picked:
        if not rs.IsLayerReference(layer):
            layers_ToDel_Picked_NoRefs.append(layer)
    
    def getLayerDescendants(layer_Parent):
        for layer_Child in rs.LayerChildren(layer_Parent):
            layers_Descendant.append(layer_Child)
            getLayerDescendants(layer_Child)
    
    layers_ToDel_WIP = layers_ToDel_Picked_NoRefs[:]
    
    # Add decendant layers if any layers have children.
    for layer in layers_ToDel_Picked_NoRefs:
        if rs.LayerChildCount(layer):
            layers_Descendant = []
            getLayerDescendants(layer)
            
            for layer_Desc in layers_Descendant:
                
                # Remove any layer already in list so that
                # descendants succeed their ancestors.
                if layer_Desc in layers_ToDel_WIP:
                    layers_ToDel_WIP.remove(layer_Desc)
                
                layers_ToDel_WIP.append(layer_Desc)
    
    if bEcho and layers_ToDel_WIP:
        print "Layers to delete:"
        for layer in layers_ToDel_WIP:
            print "  " + layer
    
    # Return the layers names to a new list in reverse order;
    # Descendants precede their ancestors.
    return layers_ToDel_WIP[::-1]


def getModifiableBlocks():
    blocks_ToEval = []
    for block in rs.BlockNames(sort=True):
        if rs.IsBlockReference(block) or not rs.IsBlockEmbedded(block): continue
        blocks_ToEval.append(block)
    return blocks_ToEval


def getBlocksToDelete(layers_ToDel, bEcho=True):
    """
    """
    
    #
    # Phase 1: Acquire names of all block definitions that contain
    # no objects on layers to keep.
    
    blocks_ToDel_P1 = []
    
    for block in getModifiableBlocks():
        bHasLayerToKeep = False
        
        for gBlockObj in rs.BlockObjects(block):
            if rs.ObjectLayer(gBlockObj) not in layers_ToDel:
                bHasLayerToKeep = True
        
        if not bHasLayerToKeep:
            blocks_ToDel_P1.append(block)
    #
    
    #
    # Phase 2: Acquire names of all block definitions that contain
    # no instances of blocks to keep and
    # no non-instance objects on layers to keep.
    
    blocks_ToDel_P2 = []
    
    for block in getModifiableBlocks():
        if block in blocks_ToDel_P1: continue
        
        bHasLayerToKeep = False
        bHasInstToKeep = False
        
        for gBlockObj in rs.BlockObjects(block):
            if rs.ObjectType(gBlockObj) == rs.filter.instance:
                block_Inst = rs.BlockInstanceName(gBlockObj)
                if block_Inst not in blocks_ToDel_P1:
                    bHasInstToKeep = True
            elif rs.ObjectLayer(gBlockObj) not in layers_ToDel:
                bHasLayerToKeep = True
            
        
        if not bHasLayerToKeep and not bHasInstToKeep:
            blocks_ToDel_P2.append(block)
    
    blocks_ToDel_All = sorted(blocks_ToDel_P1 + blocks_ToDel_P2)
    
    if bEcho and blocks_ToDel_All:
        print "Blocks to delete:"
        for block in blocks_ToDel_All:
            print "  " + block
    
    return blocks_ToDel_All


def getBlocksToModify(layers_ToDel, blocks_ToDel, bEcho=True):
    """
    Acquire names of all block definitions to modify. This includes those that
    are not set to be deleted and contain
    (objects on layers to delete) or
    (instances of block definitions to delete).
    """
    
    blocks_ToMod = []
    
    for block in getModifiableBlocks():
        if block in blocks_ToDel: continue
        
        bHasLayerToDel = False
        bHasInstToDel = False
        
        for gBlockObj in rs.BlockObjects(block):
            if rs.ObjectType(gBlockObj) == rs.filter.instance:
                if rs.BlockInstanceName(gBlockObj) in blocks_ToDel:
                    bHasInstToDel = True
            
            if rs.ObjectLayer(gBlockObj) in layers_ToDel:
                bHasLayerToDel = True
        
        if bHasLayerToDel or bHasInstToDel:
            blocks_ToMod.append(block)
    
    if bEcho and blocks_ToMod:
        print "Blocks to modify (contain combination of instances and/or" \
              " other objects to delete and not to delete):"
        for block in blocks_ToMod:
            print "  " + block
    
    return blocks_ToMod


def deleteObjectsEx(object_ids, ignore_modes=False):
    # From https://discourse.mcneel.com/t/rs-deleteobject-with-object-in-locked-layer/39988/10
    rc = 0
    object_ids = rs.coerceguidlist(object_ids)
    if object_ids:
        for id in object_ids:
            rhobj = rs.coercerhinoobject(id)
            if rhobj and sc.doc.Objects.Delete(rhobj, True, ignore_modes): rc+=1
        if rc: sc.doc.Views.Redraw()
    return rc


def addBlock(object_ids, base_point, name=None, delete_input=False):
    """
    rhinoscriptsyntax.AddBlock but allows other block definition instances to
    be included in the block.
    
    Adds a new block definition to the document
    Parameters:
      object_ids ([guid, ....]) objects that will be included in the block
      base_point (point): 3D base point for the block definition
      name (str, optional): name of the block definition. If omitted a name will be
        automatically generated
      delete_input (bool): if True, the object_ids will be deleted
    Returns:
      str: name of new block definition on success
    Example:
      import rhinoscriptsyntax as rs
      objs = rs.GetObjects("Select objects to define block")
      if objs:
          point = rs.GetPoint("Block base point")
          if point:
              block = rs.AddBlock(objs, point, None, True)
              rs.InsertBlock(block, point)
    See Also:
      InsertBlock
    """
    base_point = rs.coerce3dpoint(base_point, True)
    if not name:
        name = sc.doc.InstanceDefinitions.GetUnusedInstanceDefinitionName()
    found = sc.doc.InstanceDefinitions.Find(name)
    objects = []
    for id in object_ids:
        obj = rs.coercerhinoobject(id, True)
        if obj.IsReference: return
        ot = obj.ObjectType
        if ot==Rhino.DocObjects.ObjectType.Light: return
        if ot==Rhino.DocObjects.ObjectType.Grip: return
        if ot==Rhino.DocObjects.ObjectType.Phantom: return
        if ot==Rhino.DocObjects.ObjectType.InstanceReference and found:
            uses, nesting = obj.UsesDefinition(found.Index)
        #            if uses: return    # This is preventing blocks containing instances from being redefined.
        objects.append(obj)
    if objects:
        geometry = [obj.Geometry for obj in objects]
        attrs = [obj.Attributes for obj in objects]
        rc = 0
        if found:
          rc = sc.doc.InstanceDefinitions.ModifyGeometry(found.Index, geometry, attrs)
        else:
          rc = sc.doc.InstanceDefinitions.Add(name, "", base_point, geometry, attrs)
        if rc>=0:
            if delete_input:
                for obj in objects: sc.doc.Objects.Delete(obj, True)
            sc.doc.Views.Redraw()
    return name


def modifyBlocks(blocks_ToMod, layers_ToDel, blocks_ToDel, bEcho=True):
    """
    """
    if bEcho: print "Modifying blocks ..."
    
    ctBlocksModified = 0
    
    # Lock all objects not in blocks before modifying blocks.
    # Same list will be used later in rs.UnlockObjects.
    gObjsToLock = rs.NormalObjects()
    if gObjsToLock: rs.LockObjects(gObjsToLock)
    
    for iB, block in enumerate(blocks_ToMod):
        # Delete objects on each layer to delete.
        for gBlockObj in rs.ExplodeBlockInstance(
                rs.InsertBlock(block, (0.0, 0.0, 0.0))):
            if rs.ObjectLayer(gBlockObj) in layers_ToDel:
                rs.DeleteObject(gBlockObj)
            elif rs.ObjectType(gBlockObj) == rs.filter.instance:
                if rs.BlockInstanceName(gBlockObj) in blocks_ToDel:
                    rs.DeleteObject(gBlockObj)
        
        # Redefine block.
        block_redefined = addBlock(
                object_ids=rs.NormalObjects(),
                base_point=(0.0, 0.0, 0.0),
                name=block,
                delete_input=True)
        if block_redefined is None:
            if bEcho: print "  {} was NOT modified.".format(block)
        else:
            ctBlocksModified += 1
    
    # Unlock only objects that were previously unlocked.
    if gObjsToLock: rs.UnlockObjects(gObjsToLock)
    
    if bEcho:
        if ctBlocksModified == len(blocks_ToMod):
            print "  All {} blocks were modified.".format(ctBlocksModified)
        else:
            print "  {} out of {} blocks were modified.".format(
                    ctBlocksModified, len(blocks_ToMod))
    
    return ctBlocksModified


def deleteBlocks(blocks_ToDel, bEcho=True):
    """
    """
    
    blocks_ToDel_Removing = blocks_ToDel[:]
    if bEcho: print "Deleting blocks ..."
    ctBlocksDeleted = 0
    while blocks_ToDel_Removing:
        ctBlocksDeleted_ThisRound = 0
        for block in blocks_ToDel_Removing[:]:
            if rs.DeleteBlock(block):
                blocks_ToDel_Removing.remove(block)
                ctBlocksDeleted += 1
                ctBlocksDeleted_ThisRound += 1
            else:
                # If a block was not deleted,
                # it may be nested inside another block.
                pass
        if not ctBlocksDeleted_ThisRound:
            break
    
    if bEcho:
        if ctBlocksDeleted == len(blocks_ToDel):
            print "  All {} blocks were deleted.".format(ctBlocksDeleted)
        else:
            print "  {} out of {} blocks were deleted.".format(
                    ctBlocksDeleted, len(blocks_ToDel))
    
    return ctBlocksDeleted


def deleteLayers(layers_ToDel, bEcho=True):
    
    # Change the current layer if it is slated to be deleted.
    if rs.CurrentLayer() in layers_ToDel:
        for layer in rs.LayerNames(sort=True):
            if layer not in layers_ToDel and not rs.IsLayerReference(layer):
                rs.CurrentLayer(layer)
                break
        else:
            if bEcho:
                print "No eligible layers are available to be made current."
    
    if bEcho: print "Deleting layers ..."
    ctLayersDeleted = 0
    for layer in layers_ToDel:
        if rs.DeleteLayer(layer):
            ctLayersDeleted += 1
        else:
            if bEcho: print "  {} was NOT deleted.".format(layer)
    if bEcho:
        if ctLayersDeleted == len(layers_ToDel):
            print "  All {} layers were deleted.".format(ctLayersDeleted)
        else:
            print "  {} out of {} layers were deleted.".format(
                    ctLayersDeleted, len(layers_ToDel))
    
    return ctLayersDeleted


def main():
    
    bEcho = True
    
    layers_ToDel = getLayersToDelete(bEcho)
    if not layers_ToDel: return
    
    if not getModifiableBlocks():
        if bEcho:
            print "There are no non-reference, non-linked block definitions" \
                  " in this document."
        blocks_ToDel = blocks_ToMod = None
    else:
        blocks_ToDel = getBlocksToDelete(layers_ToDel, bEcho)
        blocks_ToMod = getBlocksToModify(layers_ToDel, blocks_ToDel, bEcho)
    
    plugin = rs.GetPlugInObject("Rhino Bonus Tools")
    layerState = "Used by deleteLayers"
    plugin.SaveLayerState(layerState)
    
    rs.EnableRedraw(False)
    
    for layer in rs.LayerNames():
        rs.LayerVisible(layer, True)
        rs.LayerLocked(layer, False)
    
    for layer in layers_ToDel:
        rc = deleteObjectsEx(rs.ObjectsByLayer(layer), ignore_modes=True)
    
    if blocks_ToMod:
        rc = modifyBlocks(blocks_ToMod, layers_ToDel, blocks_ToDel, bEcho)
    
    plugin.RestoreLayerState(layerState)
    plugin.DeleteLayerState(layerState)
    
    if blocks_ToDel:
        rc = deleteBlocks(blocks_ToDel, bEcho)
    
    rc = deleteLayers(layers_ToDel, bEcho)
    
    rs.EnableRedraw()


if __name__ == '__main__': main()