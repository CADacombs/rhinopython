"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250923: Created.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

from clr import StrongBox


def main():

    rc, objref = Rhino.Input.RhinoGet.GetOneObject(
        prompt="Select block instance",
        acceptNothing=False,
        filter=Rhino.DocObjects.ObjectType.InstanceReference)

    if rc != Rhino.Commands.Result.Success: return

    rdInst = objref.Object()

    xform_In = rdInst.InstanceXform

    if not xform_In.IsAffine:
        print("Transformation is not affine.")
        return

    #bSuccess, xform_Inverse = xform_In.TryGetInverse()
    #if not bSuccess:
    #    print("Cannot get the inverse of the transform.")
    #    return

    if xform_In.SimilarityType == rg.TransformSimilarityType.OrientationReversing:
        print("Warning, this instance is mirrored, so the instance will be exploded in a STEP export.")

    vector3d_translation_out = StrongBox[rg.Vector3d]()
    xform_rotation_out = StrongBox[rg.Transform]()
    xform_orthogonal_out = StrongBox[rg.Transform]()
    vector3d_diagonal_out = StrongBox[rg.Vector3d]()

    if not xform_In.DecomposeAffine(
        vector3d_translation_out,
        xform_rotation_out,
        xform_orthogonal_out,
        vector3d_diagonal_out
    ):
        print("Transform.DecomposeAffine failed.")
        return

    xform_Reset1 = rg.Transform.Translation(vector3d_translation_out.Value)
    xform_Reset1 *= xform_rotation_out.Value
    vDiagonal = vector3d_diagonal_out.Value
    if any(v == 0.0 for v in vDiagonal):
        print("0.0 value in diagonal vector is invalid. Check the instance.")
        return
    xform_Reset1 *= rg.Transform.Diagonal(
        1.0 if vDiagonal.X > 0.0 else -1.0,
        1.0 if vDiagonal.Y > 0.0 else -1.0,
        1.0 if vDiagonal.Z > 0.0 else -1.0,
        )
    
    inst_Out = rg.InstanceReferenceGeometry(
        instanceDefinitionId=rdInst.InstanceDefinition.Id,
        transform=xform_Reset1)
    
    if not sc.doc.Objects.Replace(objref, inst_Out, ignoreModes=True):
        print("Could not replace block instance geometry.")
        return
    
    #sc.doc.Objects.AddInstanceObject(
    #    instanceDefinitionIndex=rdInst.InstanceDefinition.Index,
    #    instanceXform=xform_Reset1)
    sc.doc.Views.Redraw()

    return


    # Alternative.
    double_dilation_out = StrongBox[float().GetType()]()

    transformSimilarityType = xform_In.DecomposeSimilarity(
        vector3d_translation_out,
        double_dilation_out,
        xform_rotation_out,
        tolerance=Rhino.RhinoMath.ZeroTolerance)
    print("\nTransformSimilarityType: {}".format(transformSimilarityType))
    if transformSimilarityType:
        print("\nDecomposeSimilarity:")
        print("  Vector3d translation:")
        print(vector3d_translation_out.ToString().replace(",", "\n"))
        print("  dilation:")
        print(double_dilation_out.ToString())
        print("  Transform rotation:")
        print(xform_rotation_out.ToString().replace(", ", "\n"))


if __name__ == '__main__': main()