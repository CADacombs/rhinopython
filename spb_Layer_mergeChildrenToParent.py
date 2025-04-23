"""
This script can be used for preparing models before exporting to layer-supported
formats, e.g, STEP.

For each descendant layer (children and all x-grandchildren):
  1. ByLayer colors of objects will be changed to ObjectColor
  2. Objects' layers will be changed to selected parent.
  3. Descendant layer itself will be deleted.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
220307: Created.
250422: Added support for blocks.
"""

import Rhino
import Rhino.DocObjects as rd
import scriptcontext as sc


class logs:
    iCt_rdObjs_modified = 0
    iCt_Layers_deleted = 0


def getInput_Parents():
    """
    Reference layers are not returned.
    Layers with no children are not returned.
    """

    bSuccess, layer_indices = Rhino.UI.Dialogs.ShowSelectMultipleLayersDialog(
        defaultLayerIndices=None,
        dialogTitle="Select parent layers",
        showNewLayerButton=False)

    if not bSuccess: return
    

    def isLayerValidParent(layer):
        if layer.IsReference:
            print("Reference layer skipped.")
            return False
        if not layer.GetChildren():
            print("Layer with no children skipped.")
            return False
        return True


    if len(layer_indices) == 1:
        layer = sc.doc.Layers.FindIndex(layer_indices[0])
        if not isLayerValidParent(layer): return
        layers = [sc.doc.Layers.FindIndex(layer_indices[0])]
    else:
        layers = []

        for i in range(len(layer_indices)):
            idxA = layer_indices[i]
            layerA = sc.doc.Layers.FindIndex(idxA)
            if not isLayerValidParent(layerA): continue
            for j in range(i+1, len(layer_indices)):
                idxB = layer_indices[j]
                layerB = sc.doc.Layers.FindIndex(idxB)
                if layerA.IsChildOf(layerB) or layerB.IsChildOf(layerA):
                    print("Selection is ambiguous because it contains descendants of others."
                          "Script canceled.",
                          sep='  ')
                    return
            layers.append(layerA)

        if layers: layers.append(layerB)

    return layers


def _getChildren(layer_Parent):


    def getLayerDescendants(layer_Parent):
        """
        This function is used recursively. Leave it nested to use
        the list variable within the scope of the outer function.
        """
        layers_Children = layer_Parent.GetChildren()
        if layer_Parent.GetChildren() is None: return []
        for layer_Child in layers_Children:
            layers_Descendant.append(layer_Child)
            getLayerDescendants(layer_Child)


    layers_Descendant = []

    getLayerDescendants(layer_Parent)

    layers_Descendant.reverse() # So that descendants precede their ancestors.

    return layers_Descendant


def _changeObjectAttributes_NoBlocks(layer_from, layer_to):
    bModification = False
    for rdObj_onLayer in list(sc.doc.Objects.FindByLayer(layer_from)):
        if rdObj_onLayer.Attributes.ColorSource == rd.ObjectColorSource.ColorFromLayer:
            rdObj_onLayer.Attributes.ObjectColor = layer_from.Color
            rdObj_onLayer.Attributes.ColorSource = rd.ObjectColorSource.ColorFromObject
        rdObj_onLayer.Attributes.LayerIndex = layer_to.LayerIndex
        rdObj_onLayer.CommitChanges()
        logs.iCt_rdObjs_modified += 1
        bModification = True
    return bModification


def _changeObjectAttributes_Blocks(layer_from, layer_to):
    rdIdefs = sc.doc.InstanceDefinitions.GetList(ignoreDeleted=True)
    if not rdIdefs:
        return

    rdIdefs_Modified = []

    for idef in rdIdefs:
        iCt_AttrsMod = 0
        rdObjs = idef.GetObjects()
        geoms = []
        attrs = []
        for rdObj in rdObjs:
            geoms.append(rdObj.Geometry)
            if rdObj.Attributes.LayerIndex == layer_from.LayerIndex:
                if rdObj.Attributes.ColorSource == rd.ObjectColorSource.ColorFromLayer:
                    rdObj.Attributes.ObjectColor = layer_from.Color
                    rdObj.Attributes.ColorSource = rd.ObjectColorSource.ColorFromObject
                rdObj.Attributes.LayerIndex = layer_to.LayerIndex
                iCt_AttrsMod += 1
            attrs.append(rdObj.Attributes)

        if not iCt_AttrsMod:
            continue

        if sc.doc.InstanceDefinitions.ModifyGeometry(idef.Index, geoms, attrs):
            rdIdefs_Modified.append(idef)
            logs.iCt_rdObjs_modified += iCt_AttrsMod
        else:
            print("Could not modify block instance, {}.".format(idef.Name))

    return bool(rdIdefs_Modified)


def mergeLayers(layer_Parent):

    layers_Children = _getChildren(layer_Parent)
    if not layers_Children: return

    for layer_Child in layers_Children:
        rv = _changeObjectAttributes_NoBlocks(layer_Child, layer_Parent)

        rv = _changeObjectAttributes_Blocks(layer_Child, layer_Parent)

        # In case layer is empty, try to delete regardless of above document changes.
        if sc.doc.Layers.Delete(layer_Child, quiet=True):
            logs.iCt_Layers_deleted += 1
        else:
            print("Could not delete layer '{}'.".format(layer_Child.Name))


def main():

    layers_Parents = getInput_Parents()
    if not layers_Parents: return

    sc.doc.Views.RedrawEnabled = False

    for layer_Parent in layers_Parents:
        rc = mergeLayers(layer_Parent)
        if not rc: continue

    sc.doc.Views.RedrawEnabled = True

    print(
        "Modified the attributes of {} objects.".format(logs.iCt_rdObjs_modified),
        "Deleted {} layers.".format(logs.iCt_Layers_deleted),
        sep='  ')


if __name__ == '__main__': main()