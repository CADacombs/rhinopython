"""
This script adds 3 lines, and optionally, 3 planes, and 4 points, representing
the X, Y, and Z axes.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
150909: Created.
170529: Changed default of Block option to No.
241007-08: Streamlined input. Refactored. Added options.
241017: Fixed typo.
241027: Split an option into two.
250219: Split again into a library module.
"""

import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc


def createLines_and_attributes(fLineLength):

    ptX = rg.Point3d(fLineLength, 0.0, 0.0)
    ptY = rg.Point3d(0.0, fLineLength, 0.0)
    ptZ = rg.Point3d(0.0, 0.0, fLineLength)

    geoms_Out = (
        rg.LineCurve(rg.Point3d.Origin, ptX),
        rg.LineCurve(rg.Point3d.Origin, ptY),
        rg.LineCurve(rg.Point3d.Origin, ptZ),
        )

    attr_Red = rd.ObjectAttributes()
    attr_Red.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr_Red.ObjectColor = attr_Red.ObjectColor.Red

    attr_Green = attr_Red.Duplicate()
    attr_Green.ObjectColor = attr_Green.ObjectColor.Lime

    attr_Blue = attr_Red.Duplicate()
    attr_Blue.ObjectColor = attr_Blue.ObjectColor.Blue

    attrs_Out = (
        attr_Red,
        attr_Green,
        attr_Blue,
        )

    return geoms_Out, attrs_Out


def createPlaneSurfaces_and_attributes(fLineLength):

    interval = rg.Interval(-fLineLength, fLineLength)

    # Y-normal plane is defined differently than others to make it aligned to the 'Front' view.
    geoms_Out = (
        rg.PlaneSurface(rg.Plane(rg.Point3d.Origin, normal=rg.Vector3d.XAxis), xExtents=interval, yExtents=interval),
        rg.PlaneSurface(rg.Plane(rg.Point3d.Origin, rg.Vector3d.XAxis, rg.Vector3d.ZAxis), xExtents=interval, yExtents=interval),
        rg.PlaneSurface(rg.Plane(rg.Point3d.Origin, normal=rg.Vector3d.ZAxis), xExtents=interval, yExtents=interval),
        )

    attr_Red = rd.ObjectAttributes()
    attr_Red.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr_Red.ObjectColor = attr_Red.ObjectColor.Red

    attr_Green = attr_Red.Duplicate()
    attr_Green.ObjectColor = attr_Green.ObjectColor.Lime

    attr_Blue = attr_Red.Duplicate()
    attr_Blue.ObjectColor = attr_Blue.ObjectColor.Blue

    attrs_Out = (
        attr_Red,
        attr_Green,
        attr_Blue,
        )

    return geoms_Out, attrs_Out


def createPoint3d_and_attributes_OriginPt():
    attr_White = rd.ObjectAttributes()
    attr_White.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    attr_White.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr_White.ObjectColor = attr_White.ObjectColor.White

    return rg.Point(rg.Point3d.Origin), attr_White


def createPoint3ds_and_attributes_Axes(fLineLength):

    ptX = rg.Point3d(fLineLength, 0.0, 0.0)
    ptY = rg.Point3d(0.0, fLineLength, 0.0)
    ptZ = rg.Point3d(0.0, 0.0, fLineLength)

    geoms_Out = (
        rg.Point(ptX),
        rg.Point(ptY),
        rg.Point(ptZ),
        )

    attr_Red = rd.ObjectAttributes()
    attr_Red.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr_Red.ObjectColor = attr_Red.ObjectColor.Red

    attr_Green = attr_Red.Duplicate()
    attr_Green.ObjectColor = attr_Green.ObjectColor.Lime

    attr_Blue = attr_Red.Duplicate()
    attr_Blue.ObjectColor = attr_Blue.ObjectColor.Blue

    attrs_Out = (
        attr_Red,
        attr_Green,
        attr_Blue,
        )

    return geoms_Out, attrs_Out


def addLineCurves(fLineLength):
    gOuts = []
    for Lc, attr in zip(*createLines_and_attributes(fLineLength)):
        gOuts.append(sc.doc.Objects.AddCurve(Lc, attributes=attr))
    return gOuts


def addSurfaceObjects_Planes(fLineLength):
    gOuts = []
    for ps, attr in zip(*createPlaneSurfaces_and_attributes(fLineLength)):
        gOuts.append(sc.doc.Objects.AddSurface(ps, attributes=attr))
    return gOuts


def addPointObject_Origin():
    pt, attr = createPoint3d_and_attributes_OriginPt()
    gOut = sc.doc.Objects.Add(pt, attributes=attr)
    return gOut


def addPointObjects_Axes(fLineLength):
    gOuts = []
    for pt, attr in zip(*createPoint3ds_and_attributes_Axes(fLineLength)):
        gOuts.append(sc.doc.Objects.Add(pt, attributes=attr))
    return gOuts


def main():
    print("This script is only a library. Try other scripts with similar names.")


if __name__ == '__main__': main()