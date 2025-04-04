"""
This script wraps _Rotate3D with the additional option of inferring the rotation axis
from a post-selected object, e.g., linear curve, arc curve, cylinder.
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
240731-0801: Created.
240805: All options of _Rotate3D are now available, therefore this script can be used
        in place of that command.
240901: Bug fix related to behavior of GetObject.CustomGeometryFilter.
250225: Bug fix concerning the active CPlane. Translated from Python 3 to 2.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'fLineTol'
    keys.append(key)
    values[key] = max((1e-6, 0.001 * sc.doc.ModelAbsoluteTolerance))
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = "{}({})({})".format(key, __file__, sc.doc.Name)

    key = 'fArcTol'
    keys.append(key)
    values[key] = max((1e-6, 0.001 * sc.doc.ModelAbsoluteTolerance))
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = "{}({})({})".format(key, __file__, sc.doc.Name)

    key = 'fFaceTol'
    keys.append(key)
    values[key] = max((1e-6, 0.001 * sc.doc.ModelAbsoluteTolerance))
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = "{}({})({})".format(key, __file__, sc.doc.Name)

    key = 'bEcho'
    keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], "No", "Yes")
    stickyKeys[key] = "{}({})".format(key, __file__)

    key = 'bDebug'
    keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], "No", "Yes")
    stickyKeys[key] = "{}({})".format(key, __file__)

    for key in keys:
        if key not in names:
            names[key] = key[1:]

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]

    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == "b":
                idxOpt = go.AddOptionToggle(cls.names[key], cls.riOpts[key])[0]
            elif key[0] == "f":
                idxOpt = go.AddOptionDouble(cls.names[key], cls.riOpts[key])[0]
            elif key[0] == "i":
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key]
                )[0]
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key],
            )
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt

    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fArcTol':
            if cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print(
                "Why is key, {}, here?  Value was not set or sticky-saved.".format(key)
            )
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getAxisFromArc(arc):
    xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, arc.Plane)
    line_Axis = rg.Line(
            rg.Point3d.Origin,
            rg.Point3d(0.0, 0.0, arc.Radius),
        )
    if not line_Axis.IsValid:
        sEval = "line_Axis.IsValid"; print(sEval,'=',eval(sEval))
    line_Axis.Transform(xform)
    return line_Axis


def getAxisFromCrv(crv, fLineTol, fArcTol):
    if isinstance(crv, rg.LineCurve):
        return crv.Line

    if isinstance(crv, rg.ArcCurve):
        arc = crv.Arc
        return getAxisFromArc(arc)

    if crv.IsLinear(tolerance=fLineTol):
        return rg.Line(crv.PointAtStart, crv.PointAtEnd)

    bSuccess, arc = crv.TryGetArc(tolerance=fArcTol)
    if bSuccess:
        return getAxisFromArc(arc)


def isLineOrArc(crv, fLineTol, fArcTol):
    if isinstance(crv, (rg.LineCurve, rg.ArcCurve)):
        return True

    if crv.IsLinear(tolerance=fLineTol):
        return True

    return crv.IsArc(tolerance=fArcTol)


def getAxisFromFace(face, fFaceTol):
    bSuccess, cone = face.TryGetCone(tolerance=fFaceTol)
    if bSuccess:
        return rg.Line(cone.BasePoint, cone.ApexPoint)

    bSuccess, cylinder = face.TryGetCylinder(tolerance=fFaceTol)
    if bSuccess:
        return rg.Line(cylinder.Center, cylinder.Axis)

    bSuccess, torus = face.TryGetTorus(tolerance=fFaceTol)
    if bSuccess:
        return rg.Line(torus.Plane.Origin, torus.Plane.Normal)


def getInput_ObjForRotAxis():
    """
    Get object with optional input.
    """

    def customGeomFilter(rdObj, geom, compIdx):
        # print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if isinstance(geom, rg.Curve):
            return isLineOrArc(geom, Opts.values['fLineTol'], Opts.values['fArcTol'])
        elif isinstance(geom, rg.BrepFace):
            line_Axis = getAxisFromFace(
                geom,
                Opts.values['fFaceTol']
            )
            return bool(line_Axis)
        else:
            sEval = "geom"; print(sEval,'=',eval(sEval))

        return False

    idxs_Opt = {}

    def addOption(key):
        idxs_Opt[key] = Opts.addOption(go, key)

    while True:

        go = ri.Custom.GetObject()

        go.SetCommandPrompt("Select object for rotation axis")

        go.GeometryFilter = (
            rd.ObjectType.Curve
            |
            rd.ObjectType.Surface
        )

        go.SetCustomGeometryFilter(customGeomFilter)

        go.OneByOnePostSelect = True

        go.AcceptNumber(True, acceptZero=True)

        go.AlreadySelectedObjectSelect = True
        go.DeselectAllBeforePostSelect = (
            False  # So objects won't be deselected on repeats of While loop.
        )
        go.EnableClearObjectsOnEntry(False)  # Keep objects in go on repeats of While loop.
        go.EnableUnselectObjectsOnExit(False)


        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('fLineTol')
        addOption('fArcTol')
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

        # An option was selected or a number was entered.
        if res == ri.GetResult.Number:
            key = 'fArcTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getAxisFromObjRef(objref_Axis, fLineTol, fArcTol, fFaceTol):

    crv_In = objref_Axis.Curve()

    if crv_In:
        # rgC_Dup = crv_In.DuplicateCurve()
        line_Axis = getAxisFromCrv(crv_In, fLineTol, fArcTol)
        # rgC_Dup.Dispose()
        if line_Axis is None:
            raise Exception("Cannot obtain Line or Arc from reference curve.")
    else:
        face_In = objref_Axis.Surface()
        if face_In:
            line_Axis = getAxisFromFace(
                face_In,
                Opts.values['fFaceTol']
            )
            if line_Axis is None:
                raise Exception("Cannot obtain reference axis from face.")
        else:
            raise Exception(
                "Axis extraction from {} not implemented yet.".format(face_In))

    return line_Axis


def main():

    res, objrefs_ToRotate = ri.RhinoGet.GetMultipleObjects(
        "Select objects to rotate",
        acceptNothing=False,
        filter=rd.ObjectType.AnyObject)
    if res != Rhino.Commands.Result.Success: return

    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    gp = ri.Custom.GetPoint()

    gp.SetCommandPrompt("Start of rotation axis")

    gp.AcceptNumber(False, acceptZero=False)

    idxs_Opt = {}

    idxs_Opt['SurfaceNormal'] = gp.AddOption('SurfaceNormal')
    idxs_Opt['DeriveFromObject'] = gp.AddOption('DeriveFromObject')

    res = gp.Get()

    if res == ri.GetResult.Cancel:
        gp.Dispose()
        return

    if res == ri.GetResult.Point:
        pt = gp.Point()
        gp.Dispose()

        sc.doc.Objects.UnselectAll()

        for o in objrefs_ToRotate:
            sc.doc.Objects.Select(o)

        Rhino.RhinoApp.RunScript(
            script="_Rotate3D {}".format(
                pt,
            ),
            echo=bDebug,
        )

        sc.doc.Views.Redraw()

        return


    if gp.Option().Index == idxs_Opt['SurfaceNormal']:
        gp.Dispose()

        sc.doc.Objects.UnselectAll()

        for o in objrefs_ToRotate:
            sc.doc.Objects.Select(o)

        Rhino.RhinoApp.RunScript(
            script="_Rotate3D _SurfaceNormal",
            echo=bDebug,
        )

        sc.doc.Views.Redraw()

        return


    # if gp.Option().Index == idxs_Opt['DeriveFromObject']:
    #     pass

    gp.Dispose()


    # DeriveFromObject

    sc.doc.Objects.UnselectAll()

    objref_Axis = getInput_ObjForRotAxis()
    if objref_Axis is None:
        return

    fLineTol = Opts.values['fLineTol']
    fArcTol = Opts.values['fArcTol']
    fFaceTol = Opts.values['fFaceTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    line_Axis = getAxisFromObjRef(
        objref_Axis,
        fLineTol,
        fArcTol,
        fFaceTol)

    #sc.doc.Objects.AddLine(line_Axis), sc.doc.Views.Redraw(); 1/0

    sc.doc.Objects.UnselectAll()

    for o in objrefs_ToRotate:
        sc.doc.Objects.Select(o)

    Rhino.RhinoApp.RunScript(
        script="_Rotate3D w{} w{}".format(
            line_Axis.From,
            line_Axis.To,
        ),
        echo=bDebug,
    )

    sc.doc.Views.Redraw()


if __name__ == "__main__": main()