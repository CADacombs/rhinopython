#! python3

"""
240729: Created.
240730: Added _Rotate3D feature.
240731: Now reports and selects output.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = "bIncludeRotate3D"
    keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], "No", "Yes")
    stickyKeys[key] = "{}({})".format(key, __file__)

    key = "fArcTol"
    keys.append(key)
    values[key] = max((1e-6, 0.001 * sc.doc.ModelAbsoluteTolerance))
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = "{}({})({})".format(key, __file__, sc.doc.Name)

    key = "fMinorToMajorRadiusMult"
    keys.append(key)
    values[key] = 2.0
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True, limit=1.001)
    stickyKeys[key] = "{}({})".format(key, __file__)

    key = "bEcho"
    keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], "No", "Yes")
    stickyKeys[key] = "{}({})".format(key, __file__)

    key = "bDebug"
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

        if key == "fArcTol":
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


def getInput():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select arcs and/or circles")

    go.GeometryFilter = rd.ObjectType.Curve

    def customGeomFilter(rdObj, geom, compIdx):
        # print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if isinstance(geom, rg.BrepEdge):
            # DuplicateCurve gets a Curve from the BrepEdge.
            # The curve may be a subset of the EdgeCurve.
            rgC = geom.DuplicateCurve()
        elif isinstance(geom, rg.Curve):
            rgC = geom.DuplicateCurve()
        else:
            return False

        if isinstance(rgC, rg.ArcCurve):
            rgC.Dispose()
            return True

        bSuccess, arc = rgC.TryGetArc(tolerance=Opts.values["fArcTol"])

        rgC.Dispose()

        return bSuccess

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

        addOption("bIncludeRotate3D")
        addOption("fArcTol")
        addOption("fMinorToMajorRadiusMult")
        addOption("bEcho")
        addOption("bDebug")

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

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
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected or a number was entered.
        if res == ri.GetResult.Number:
            key = "fArcTol"
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createPlaneForTorus(
    arc: rg.Arc, fMinorToMajorRadiusMult: float, bEcho=True, bDebug=False
):
    """
    Parameters:
    Returns:
    """
    plane: rg.Plane = rg.Plane.WorldXY

    majorRadius = fMinorToMajorRadiusMult * arc.Radius

    plane.Rotate(
        Rhino.RhinoMath.ToRadians(90.0),
        axis=rg.Vector3d.XAxis,
        centerOfRotation=rg.Point3d.Origin,
    )
    plane.Translate(rg.Vector3d(-majorRadius, 0.0, 0.0))
    xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, arc.Plane)
    plane.Transform(xform)

    return plane


def createTorusThruArc(
    arc: rg.Arc, fMinorToMajorRadiusMult: float, bEcho=True, bDebug=False
):
    """
    Parameters:
    Returns:
    """

    majorRadius = fMinorToMajorRadiusMult * arc.Radius

    plane = createPlaneForTorus(arc, fMinorToMajorRadiusMult, bEcho, bDebug)

    torus = rg.Torus(plane, majorRadius=majorRadius, minorRadius=arc.Radius)
    if not torus.IsValid:
        print(f"{torus.IsValid = }")
        return

    return torus


def addTorusBrepObject(
    arc: rg.Arc,
    bIncludeRotate3D: bool,
    fMinorToMajorRadiusMult: float,
    bEcho=True,
    bDebug=False,
):
    """
    Parameters:
    Returns:
    """

    xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, arc.Plane)

    lineForRotRef = rg.Line(
        rg.Point3d(0.0, 0.0, -arc.Radius),
        rg.Point3d(0.0, 0.0, arc.Radius),
    )
    if not lineForRotRef.IsValid:
        print(f"{lineForRotRef.IsValid = }")
    lineForRotRef.Transform(xform)

    # sc.doc.Objects.AddLine(lineForRotRef)

    torus = createTorusThruArc(arc, fMinorToMajorRadiusMult, bEcho, bDebug)
    if torus is None:
        return

    rgB_Torus = torus.ToBrep()
    bValid, sLog = rgB_Torus.IsValidWithLog()
    if not bValid:
        print(f"Brep is not valid: {sLog}")
        return

    gB_Torus: Guid = sc.doc.Objects.AddBrep(rgB_Torus)
    if gB_Torus == gB_Torus.Empty:
        return

    majorRadius = fMinorToMajorRadiusMult * arc.Radius

    # plane = createPlaneForTorus(arc, fMinorToMajorRadiusMult, bEcho, bDebug)

    if not bIncludeRotate3D:
        sc.doc.Objects.AddLine(lineForRotRef)
        return gB_Torus

    sc.doc.Objects.UnselectAll()
    sc.doc.Objects.Select(gB_Torus)

    Rhino.RhinoApp.RunScript(
        script="_Rotate3D {} {} {}".format(
            lineForRotRef.From,
            lineForRotRef.To,
            torus.Plane.Origin,
        ),
        echo=bDebug,
    )

    return gB_Torus


def main():

    objrefs = getInput()
    if objrefs is None:
        return

    bIncludeRotate3D = Opts.values["bIncludeRotate3D"]
    fArcTol = Opts.values["fArcTol"]
    fMinorToMajorRadiusMult = Opts.values["fMinorToMajorRadiusMult"]
    bEcho = Opts.values["bEcho"]
    bDebug = Opts.values["bDebug"]

    gBs_Tori = []

    objref: rd.ObjRef
    for objref in objrefs:
        rgC_Dup = objref.Curve().DuplicateCurve()

        arc: rg.Arc
        bSuccess, arc = rgC_Dup.TryGetArc(tolerance=fArcTol)

        if not bSuccess:
            continue

        gB_Torus = addTorusBrepObject(
            arc=arc,
            bIncludeRotate3D=bIncludeRotate3D,
            fMinorToMajorRadiusMult=fMinorToMajorRadiusMult,
            bEcho=bEcho,
            bDebug=bDebug,
        )
        if gB_Torus is not None:
            gBs_Tori.append(gB_Torus)

        rgC_Dup.Dispose()

    if gBs_Tori:
        sc.doc.Objects.UnselectAll()
        for g in gBs_Tori:
            sc.doc.Objects.Select(g)

    if len(gBs_Tori) == len(
        list(sc.doc.Objects.GetSelectedObjects(
            includeLights=False, includeGrips=False))
    ):
        print(f"{len(gBs_Tori)} tori were created and are selected.")

    sc.doc.Views.Redraw()


if __name__ == "__main__": main()