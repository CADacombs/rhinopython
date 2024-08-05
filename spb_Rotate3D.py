"""
This script wraps _Rotate3D with the axis of rotation inferred from a
post-selected object.
"""

#! python3

"""
240731-0801: Created.
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
        print(f"{line_Axis.IsValid = }")
    line_Axis.Transform(xform)
    return line_Axis


def getAxisFromCrv(crv: rg.Curve, fLineTol: float, fArcTol: float):
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


def getAxisFromFace(face: rg.BrepFace, fFaceTol: float):
    bSuccess, cone = face.TryGetCone(tolerance=fFaceTol)
    if bSuccess:
        return rg.Line(cone.BasePoint, cone.ApexPoint)

    bSuccess, cylinder = face.TryGetCylinder(tolerance=fFaceTol)
    if bSuccess:
        return rg.Line(cylinder.Center, cylinder.Axis)

    bSuccess, torus = face.TryGetTorus(tolerance=fFaceTol)
    if bSuccess:
        torus: rg.Torus
        return rg.Line(torus.Plane.Origin, torus.Plane.Normal)


def getInput():
    """
    Get object with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select object for axis of rotation")

    go.GeometryFilter = (
        rd.ObjectType.Curve
        |
        rd.ObjectType.Surface
    )

    def customGeomFilter(rdObj, geom, compIdx):
        # print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        rgC: rg.Curve

        if isinstance(geom, rg.Curve):
            rgC_Dup = geom.DuplicateCurve()
            line_Axis = getAxisFromCrv(
                rgC_Dup,
                Opts.values['fLineTol'],
                Opts.values['fArcTol'])
            rgC_Dup.Dispose()
            return bool(line_Axis)
        elif isinstance(geom, rg.BrepFace):
            line_Axis = getAxisFromFace(
                geom,
                Opts.values['fFaceTol']
            )
            return bool(line_Axis)
        else:
            print(f"{geom = }")

        return False

    go.SetCustomGeometryFilter(customGeomFilter)

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = (
        False  # So objects won't be deselected on repeats of While loop.
    )
    go.EnableClearObjectsOnEntry(False)  # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key):
        idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('fLineTol')
        addOption('fArcTol')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue

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


def getAxisFromObjRef(objref_Axis: rd.ObjRef, fLineTol: float, fArcTol: float, fFaceTol: float):
    crv_In = objref_Axis.Curve()
    if crv_In:
        rgC_Dup = crv_In.DuplicateCurve()
        line_Axis = getAxisFromCrv(rgC_Dup, fLineTol, fArcTol)
        rgC_Dup.Dispose()
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
            raise Exception("Not implemented axis extraction from {} yet.".format(srf_In))

    return line_Axis


def main():

    res, objrefs_ToRotate = ri.RhinoGet.GetMultipleObjects(
        "Select objects to rotate",
        acceptNothing=False,
        filter=rd.ObjectType.AnyObject)
    if res != Rhino.Commands.Result.Success: return


    objref_Axis = getInput()
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

    sc.doc.Objects.UnselectAll()

    for o in objrefs_ToRotate:
        sc.doc.Objects.Select(o)

    Rhino.RhinoApp.RunScript(
        script="_Rotate3D {} {}".format(
            line_Axis.From,
            line_Axis.To,
        ),
        echo=bDebug,
    )

    sc.doc.Views.Redraw()


if __name__ == "__main__": main()
