"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250923: Created.
"""

import Rhino
import Rhino.Geometry as rg

from clr import StrongBox


def main():

    rc, objref = Rhino.Input.RhinoGet.GetOneObject(
        prompt="Select block instance",
        acceptNothing=False,
        filter=Rhino.DocObjects.ObjectType.InstanceReference)

    if rc != Rhino.Commands.Result.Success: return

    rdInst = objref.Object()

    xform = rdInst.InstanceXform

    xform_linear_out = StrongBox[rg.Transform]()
    vector3d_translation_out = StrongBox[rg.Vector3d]()
    xform_rotation_out = StrongBox[rg.Transform]()
    xform_orthogonal_out = StrongBox[rg.Transform]()
    vector3d_diagonal_out = StrongBox[rg.Vector3d]()

    if not xform.IsAffine:
        print("\nTransformation is not affine.")
    else:
        if xform.DecomposeAffine(
            xform_linear_out,
            vector3d_translation_out
        ):
            print("DecomposeAffine:")
            print("  Transform linear:")
            #print(str(xform_linear_out).split(", "))
            print(xform_linear_out.ToString().replace(", ", "\n"))
            print("  Vector3d translation:")
            print(vector3d_translation_out.ToString().replace(",", "\n"))

        if xform.DecomposeAffine(
            vector3d_translation_out,
            xform_linear_out
        ):
            print("\nDecomposeAffine:")
            print("  Vector3d translation:")
            print(vector3d_translation_out.ToString().replace(",", "\n"))
            print("  Transform linear:")
            print(xform_linear_out.ToString().replace(", ", "\n"))

        if xform.DecomposeAffine(
            vector3d_translation_out,
            xform_rotation_out,
            xform_orthogonal_out,
            vector3d_diagonal_out
        ):
            print("\nDecomposeAffine:")
            print("  Vector3d translation:")
            print(vector3d_translation_out.ToString().replace(",", "\n"))
            print("  Transform rotation:")
            print(xform_rotation_out.ToString().replace(", ", "\n"))
            print("  Transform orthogonal:")
            print(xform_orthogonal_out.ToString().replace(", ", "\n"))
            print("  Vector3d diagonal:")
            print(vector3d_diagonal_out.ToString().replace(",", "\n"))

        transformRigidType = xform.DecomposeRigid(
            vector3d_translation_out,
            xform_linear_out,
            tolerance=Rhino.RhinoMath.ZeroTolerance)
        print("\nTransformRigidType: {}".format(transformRigidType))
        if transformRigidType:
            print("n\DecomposeRigid:")
            print("  Vector3d translation:")
            print(vector3d_translation_out.ToString().replace(",", "\n"))
            print("  Transform linear:")
            print(xform_linear_out.ToString().replace(", ", "\n"))

        double_dilation_out = StrongBox[float().GetType()]()

        transformSimilarityType = xform.DecomposeSimilarity(
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

    if not rg.Transform.IsLinear:
        print("\nTransformation is not linear.")
    else:
        xform_matrix_out = StrongBox[rg.Transform]()
        if xform.DecomposeSymmetric(
            xform_matrix_out,
            vector3d_diagonal_out
        ):
            print("\nDecomposeSymmetric:")
            print("  Transform matrix:")
            print(xform_matrix_out.ToString().replace(", ", "\n"))
            print("  Vector3d diagonal:")
            print(vector3d_diagonal_out.ToString().replace(",", "\n"))
        else:
            print("\nTransformation is not symmetric.")

if __name__ == '__main__': main()