"""
This script is an alternative to _Block.
_Block creates blocks rotationally relative to the WorldXY plane.
This script creates blocks rotationally relative to the current CPlane.

Due to a bug,
https://discourse.mcneel.com/t/bug-instanceobject-usesdefinition-returns-false-positives/211552
in rhinoscript in Rhino versions ? - 8.23 - ? blocks cannot be redefined with AddBlock.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
230406: Created.
251110: Added StringBox to get block name. Now uses only rhinoscriptsyntax.
"""

import rhinoscriptsyntax as rs


def main():
    
    gObjs = rs.GetObjects(
        message="Select objects to define block",
        filter=0,
        preselect=True)
    if not gObjs: return
    
    basePoint = rs.GetPoint(
        message="Block base point")
    if not basePoint: return
    
    blockPlane = rs.MovePlane(
        plane=rs.ViewCPlane(),
        origin=basePoint)
    
    while True:
        name = rs.StringBox(
            message="Enter name of block",
            title="Create Block per Current CPlane")
        if not name: return
        
        if not rs.IsBlock(name):
            break
        
        print('Block with name, "{}",'.format(name), "already exists."
            "\nBlocks cannot be redfined with ? - Rhino 8.23 - ? version of "
            "rhinoscriptsyntax."
            "\nTODO: Check for fix in future versions or add a workaround.")
    
    xform = rs.XformChangeBasis(
        initial_plane=rs.WorldXYPlane(),
        final_plane=blockPlane)
    
    rs.EnableRedraw(False)
    
    gTransformed_Objs = rs.TransformObjects(
        object_ids=gObjs,
        matrix=xform,
        copy=False)
    
    iCt_XformFail = gTransformed_Objs.count(None)
    
    if iCt_XformFail:
        print("{} objects could not be transformed to WCS.  Check results.".format(
        iCt_XformFail))
        return
    
    sBlock = rs.AddBlock(
        object_ids=gObjs,
        base_point=(0,0,0),
        name=name,
        delete_input=True)
    
    gBlockInst = rs.InsertBlock(
        block_name=sBlock,
        insertion_point=(0,0,0),
        )
    
    gTransformed_Inst = rs.TransformObjects(
        object_ids=[gBlockInst],
        matrix=rs.XformInverse(xform),
        copy=False)
    
    if gTransformed_Objs is None:
        print("Block instance was not transformed.")
        return
    
    print('"{}" was successfully created with {} objects.'.format(
        name,
        len(gObjs) - iCt_XformFail))
    
    rs.EnableRedraw()


if __name__ == '__main__': main()