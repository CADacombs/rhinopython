"""
160528: Created.
160618: 
"""

import rhinoscriptsyntax as rs
import scriptcontext as sc
import System


def areLayerAndAllAncestorsVisible(idx_rdLayer_Child):
    
    rdLayer_Child = sc.doc.Layers[idx_rdLayer_Child]
    if not rdLayer_Child.IsVisible: return False
    
    idLayer_Parent = rdLayer_Child.ParentLayerId
    if idLayer_Parent == System.Guid.Empty: return True
    else:
        idx_rdLayer_Parent = sc.doc.Layers.Find(idLayer_Parent, True)
        return areLayerAndAllAncestorsVisible(idx_rdLayer_Parent)


def getLayerIdsFromNames(sLayerNames):
    """
    Parameter: List of layer names
    Returns: List of layer GUID's
    """
    
    # Get all layer GUID's in document.
    idx_Layers = rs.LayerIds()
    if idx_Layers is None: return
    
    idx_Layers_Target = [idx_Layer for idx_Layer in idx_Layers if rs.LayerName(idx_Layer) in sLayerNames]
    return idx_Layers_Target


if __name__ == "__main__": print getLayerIdsFromNames(["Default"])