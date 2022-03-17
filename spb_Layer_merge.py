"""
190925: Created.
"""

import Rhino
import rhinoscriptsyntax as rs
import scriptcontext as sc


def getLayersToDelete(bEcho=True):
    """
    """
    
    layers_ToDel_Picked = rs.GetLayers(
            "Select layers from which to transfer objects" \
            " (Reference layers will be ignored.)")
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
    
    # Return the layers names to a new list in reverse order;
    # Descendants precede their ancestors.
    return layers_ToDel_WIP[::-1]


def getTargetLayer(layers_ToDel, bEcho=True):
    """
    """
    
    layer = None
    while layer is None:
        layer = rs.GetLayer(
                "Select layer into which to merge " \
                "(Reference layers will be ignored.)")
        if not layer: return
        
        if rs.IsLayerReference(layer):
            layer = rs.GetLayer(
                    "Reference layer chosen.  Pick again." \
                    "Select layer into which to merge " \
                    "(Reference layers will be ignored.)")
            if not layer: return
        
        if layer in layers_ToDel:
            layer = rs.GetLayer(
                    "Layer to delete was chosen.  Pick again.  " \
                    "Select layer into which to merge " \
                    "(Reference layers will be ignored.)")
            if not layer: return
    
    return layer


def getModifiableBlocks():
    blocks_ToEval = []
    for block in rs.BlockNames(sort=True):
        if rs.IsBlockReference(block) or not rs.IsBlockEmbedded(block): continue
        blocks_ToEval.append(block)
    return blocks_ToEval


def setLayers(layers_ToDel, layer_Target, bEcho):
    """
    """
    
    ctObjsTrans_Pass = 0
    
    for gObj in rs.AllObjects():
        if rs.ObjectLayer(gObj) in layers_ToDel:
            if rs.ObjectLayer(gObj, layer_Target):
                ctObjsTrans_Pass += 1
            else:
                ctObjsTrans_Fail += 1
    
    for block in getModifiableBlocks():
        ctObjsTrans_Fail = 0
        for gObj in rs.BlockObjects(block):
            if rs.ObjectLayer(gObj) in layers_ToDel:
                if rs.ObjectLayer(gObj, layer_Target):
                    ctObjsTrans_Pass += 1
                else:
                    ctObjsTrans_Fail += 1
        if bEcho and ctObjsTrans_Fail:
            print "Layers for {} objects could not be reassigned in {}.".format(
                    ctObjsTrans_Fail, block)
    
    if bEcho and ctObjsTrans_Pass:
        print "{} objects' layers were set to {}.".format(
                ctObjsTrans_Pass, layer_Target)
    
    return ctObjsTrans_Pass


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
    
    ctLayersDeleted = 0
    for layer in layers_ToDel:
        if rs.DeleteLayer(layer):
            ctLayersDeleted += 1
        else:
            if bEcho: print "{} was NOT deleted.".format(layer)
    if bEcho:
        if ctLayersDeleted == len(layers_ToDel):
            print "All {} layers were deleted.".format(ctLayersDeleted)
        else:
            print "{} out of {} layers were deleted.".format(
                    ctLayersDeleted, len(layers_ToDel))
    
    return ctLayersDeleted


def main():
    
    bEcho = True
    
    layers_ToDel = getLayersToDelete(bEcho)
    if not layers_ToDel: return
    
    layer_Target = getTargetLayer(layers_ToDel, bEcho)
    if not layer_Target: return
    
    rs.EnableRedraw(False)
    
    rc = setLayers(layers_ToDel, layer_Target, bEcho)
    
    rc = deleteLayers(layers_ToDel, bEcho)
    
    rs.EnableRedraw()


if __name__ == '__main__': main()