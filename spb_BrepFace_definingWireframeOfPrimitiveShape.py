"""
This script will create reference geometry for BrepFaces whose shapes are can be represented
by primitives (analytic shapes).

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
170529: Created.
171010: Bug fix: Line layer for plane representation is now set to the current layer.
190207, ..., 0529: Updated an import name.
190603: Commented out seam creation of cylinders.
190807: Import-related change.
190830: Added an import to make cylinder finite when needed.
191118: Import-related update.
201212: Added point to output for plane that can be used for snapping.
210312: Line for plane is now 10 mm long for mm models.
210316: Modified and option default value.
210327: Added Opts.  Refactored.  Changed tolerances used by TryGet... methods.  Added more info to printed output.
210328: Added normal vector to printed feedback of plane.
210703: Added some options.
250326: Now prints value of vector of a cylinder's axis.
250826-27: Bug fix. Added printed center point for sphere. Modified some options.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import xBlock
import xBrepFace_tryGetPrimitiveShape
import xPrimitiveShape

import math


_MY_ZERO = 1e-6 * Rhino.RhinoMath.UnitScale(Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)

def _decimalPlacesForZero():
    return int(abs(math.log10(abs(_MY_ZERO)))) + 1

#_DECIMAL_PLACES_FOR_0 = int(abs(math.log10(abs(_MY_ZERO)))) + 1


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fShapeTol'; keys.append(key)
    values[key] = 0.001 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAddWireframe'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bIncludeArcCrvsForRoundShapes'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)


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
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fShapeTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                print("TryGet(shape) methods appear to use ZeroTolerance" \
                    " as the minimum when less is input.")
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get face with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face")

    go.GeometryFilter = (
        rd.ObjectType.Surface |
        rd.ObjectType.InstanceReference)

    go.AcceptNumber(True, True)
    go.EnableHighlight(False)
    
    rgFace = None

    idxs_Opt = {}

    while rgFace is None:
        while True:
            go.ClearCommandOptions()

            idxs_Opt.clear()

            def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

            addOption('fShapeTol')
            addOption('bAddWireframe')
            if Opts.values['bAddWireframe']:
                addOption('bIncludeArcCrvsForRoundShapes')
            addOption('bEcho')
            addOption('bDebug')

            res = go.Get()

            if res == ri.GetResult.Cancel:
                go.Dispose()
                return

            if res == ri.GetResult.Object:
                objref = go.Object(0)
                rdObj = objref.Object()
                ptPicked = objref.SelectionPoint()
        
                if go.ObjectsWerePreselected:
                    if rdObj.ObjectType == rd.ObjectType.Brep:
                        rgObj = rdObj.BrepGeometry
                    elif rdObj.ObjectType == rd.ObjectType.Extrusion:
                        rgObj = rdObj.Geometry.ToBrep()
                    else:
                        sc.doc.Objects.UnselectAll() # Necessary if instance reference was preselected.
                        sc.doc.Views.Redraw()
                        continue
                    if rgObj.Faces.Count == 1: rgFace = rgObj.Faces[0]
                if rdObj.ObjectType == rd.ObjectType.InstanceReference:
                    rgFace, ptPicked = xBlock.tryPickedFaceOfBlock(rdObj, ptPicked)
                else: rgFace = objref.Face() # This also works for ExtrusionObject.
        
                if rgFace is None:
                    sc.doc.Objects.UnselectAll() # Necessary when gb.Get() is repeated.
                else:
                    go.Dispose()

                    return (
                        rgFace,
                        Opts.values['fShapeTol'],
                        Opts.values['bAddWireframe'],
                        Opts.values['bIncludeArcCrvsForRoundShapes'],
                        Opts.values['bEcho'],
                        Opts.values['bDebug'],
                        )

            if res == ri.GetResult.Number:
                key = 'fShapeTol'
                Opts.riOpts[key].CurrentValue = go.Number()
                Opts.setValue(key)
                continue

            # An option was selected.
            for key in idxs_Opt:
                if go.Option().Index == idxs_Opt[key]:
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)
                    break


def createReport(rgShape, fTol):
    """
    """

    if isinstance(rgShape, rg.Plane):
        return ("Plane with {:.{sd}g},{:.{sd}g},{:.{sd}g} ({sd} sig. digits) normal vector found at a tolerance of {}.".format(
            rgShape.Normal.X,
            rgShape.Normal.Y,
            rgShape.Normal.Z,
            fTol,
            sd=15,
            ))

    if isinstance(rgShape, rg.Cylinder):
        return ("R{:.{dp}f} cylinder with {:.{sd}g},{:.{sd}g},{:.{sd}g} ({sd} sig. digits) axis vector found at a tolerance of {}.".format(
            rgShape.Radius,
            rgShape.Axis.X,
            rgShape.Axis.Y,
            rgShape.Axis.Z,
            fTol,
            dp=_decimalPlacesForZero(),
            sd=15,
            ))

    if isinstance(rgShape, rg.Sphere):
        return ("R{:.{dp}f} sphere with {:.{dp}f},{:.{dp}f},{:.{dp}f} center point found at a tolerance of {}.".format(
            rgShape.Radius,
            rgShape.Center.X,
            rgShape.Center.Y,
            rgShape.Center.Z,
            fTol,
            dp=_decimalPlacesForZero(),
            ))

    if isinstance(rgShape, rg.Torus):
        return ("R{0:.{2}f},r{1:.{2}f} torus found at a tolerance of {3}.".format(
            rgShape.MajorRadius,
            rgShape.MinorRadius,
            sc.doc.ModelDistanceDisplayPrecision,
            fTol))

    if isinstance(rgShape, rg.Cone):
        return ("{:.{sd}f} degree cone with {:.{sd}g},{:.{sd}g},{:.{sd}g} ({sd} sig. digits) axis vector found at a tolerance of {}.".format(
            rgShape.AngleInDegrees(),
            rgShape.Axis.X,
            rgShape.Axis.Y,
            rgShape.Axis.Z,
            fTol,
            sd=15,
            ))

    return ("What shape is this?: {}".format(
        rgShape.GetType().Name.ToString(),
        ))


def getFaceShape(rgFace0, fShapeTol=None, bDebug=False):
    """
    """
    
    rgSrf = rgFace0.UnderlyingSurface()

    fTols_ToTry = Rhino.RhinoMath.ZeroTolerance, 1e-9, 1e-8
    
    # Check whether face's UnderlyingSurface is already a primitive.
    for fTol in fTols_ToTry:
        if fShapeTol is not None and fShapeTol < fTol: break

        if rgSrf.IsPlanar(fTol):
            intervalU = rgSrf.Domain(0)
            intervalV = rgSrf.Domain(1)
            
            rgShape = rg.Plane(rgSrf.PointAt(intervalU.Mid, intervalV.Mid),
                    rgSrf.PointAt(intervalU.Mid + intervalU.Length / 2.0, intervalV.Mid),
                    rgSrf.PointAt(intervalU.Mid, intervalV.Mid + intervalV.Length / 2.0))
            return rgShape, fTol
    
    for fTol in fTols_ToTry:
        if fShapeTol is not None and fShapeTol < fTol: break

        b, rgShape = rgSrf.TryGetCylinder(fTol)
        if b:
            return rgShape, fTol
        
        b, rgShape = rgSrf.TryGetCone(fTol)
        if b:
            return rgShape, fTol
        
        b, rgShape = rgSrf.TryGetSphere(fTol)
        if b:
            return rgShape, fTol
        
        b, rgShape = rgSrf.TryGetTorus(fTol)
        if b:
            return rgShape, fTol
    
    if fShapeTol is None: return
    
    # Try to get the primitive of face
    rc = xBrepFace_tryGetPrimitiveShape.tryGetPrimitiveShape(
        rgFace0, fTolerance=fShapeTol)

    if rc[0] is None:
        print("No matching primitive found.")
        return
    
    rgShape, fTol_PrimitiveMatch, bShrunkUsed = rc[0]
    
    typeShape = rgShape.GetType()
    
    if not rgShape.IsValid:
        print("{} is not valid.".format(typeShape))
        return
    
    if (typeShape == Rhino.Geometry.Plane or
            typeShape == Rhino.Geometry.Cone or
            typeShape == Rhino.Geometry.Sphere or
            typeShape == Rhino.Geometry.Torus):
        return rgShape, fTol_PrimitiveMatch
    
    # Perform any needed repair to the cylinder to make it usable.
    if typeShape == rg.Cylinder:
        rgCyl = rgShape
        # Fix cylinder that is infinite (has 0 height).
        if not rgCyl.IsFinite:
            print("Making cylinder finite...")
            rgBbox_Face = rgFace0.GetBoundingBox(True)
            fAddLen = 1.1 * rgBbox_Face.Min.DistanceTo(rgBbox_Face.Max)
            rgCyl.Height1 = -fAddLen
            rgCyl.Height2 = fAddLen
        
        
        #    for sAttr in dir(rgCyl):
        #        if (sAttr != 'Unset' and sAttr != '__doc__' and
        #                not callable(getattr(rgCyl, sAttr))):
        #            sPrint = 'rgCyl.' + sAttr; print(sPrint + ':', eval(sPrint)
    
    #        sPrint = 'dir(rgCyl)'; print(sPrint + ':', eval(sPrint)
        #sPrint = 'rgCyl.Unset'; print(sPrint + ':', eval(sPrint)
    #        for sAttr in dir(rgCyl):
    #            if (sAttr != 'Unset' and sAttr != '__doc__' and
    #                    not callable(getattr(rgCyl, sAttr))):
    #                sPrint = 'rgShape.' + sAttr; print(sPrint + ':', eval(sPrint)
        
        elif rgCyl.TotalHeight < 2.0*sc.doc.ModelAbsoluteTolerance:
            print("Resultant height of cylinder is {}.".format(
                    rgCyl.TotalHeight))
            print("Using 100 times model's absolute tolerance.")
            rgCyl.Height1 = -100.*sc.doc.ModelAbsoluteTolerance
            rgCyl.Height2 = 100.*sc.doc.ModelAbsoluteTolerance
        else:
            # Increase length of cylinder so that surface is always larger than trim.
            rgBbox_Face = rgFace0.GetBoundingBox(False)
            rgLines = rgBbox_Face.GetEdges()
            fMaxLen = 0
            for rgLine in rgLines:
                if rgLine.Length > fMaxLen: fMaxLen = rgLine.Length
            fAddLen = 1.1 * fMaxLen
            #fAddLen = 1.1 * rgBbox_Face.Min.DistanceTo(rgBbox_Face.Max)
            rgCyl.Height1 = -fAddLen
            rgCyl.Height2 = fAddLen
        
        return rgCyl, fTol_PrimitiveMatch
        
    else:
        print("What happened?")
        return


def addWireframeOfShape(rgShape, rgFace0, bIncludeArcCrvsForRoundShapes=False, bDebug=False):
    """
    """
    
    typeShape = rgShape.GetType()
    
    if typeShape == rg.Plane:
        plane = rgShape
        #        planeSrf = rg.PlaneSurface(plane, rg.Interval(0.0, 1.0), rg.Interval(0.0, 1.0))
        #        sc.doc.Objects.AddSurface(planeSrf)
        #sc.doc.Objects.AddLine(plane.Origin, plane.Origin + plane.XAxis)
        #sc.doc.Objects.AddLine(plane.Origin, plane.Origin + plane.YAxis)
        attr = rd.ObjectAttributes()
        attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        attr.ObjectDecoration = rd.ObjectDecoration.EndArrowhead
        sc.doc.Objects.AddPoint(plane.Origin)
        units = sc.doc.ModelUnitSystem
        fLineLen = 10.0 if units == Rhino.UnitSystem.Millimeters else 1.0
        sc.doc.Objects.AddLine(plane.Origin, plane.Origin + fLineLen * plane.ZAxis, attr)
        #sc.doc.Objects.AddPoint(plane.OriginX, plane.OriginY + rg.Point3d(plane.YAxis), plane.OriginZ)
    
    elif typeShape == rg.Cylinder:
        cyl = rgShape
        #sc.doc.Objects.AddPoint(cyl.Center)
    #        if line_Axis.Length < 2.0 * sc.doc.ModelAbsoluteTolerance:
    #            line_Axis.Length = 1.0
    #            print("Line length was too short and was set to 1.0."
        if bIncludeArcCrvsForRoundShapes:
            circle_Mid = cyl.CircleAt((cyl.Height1 + cyl.Height2) / 2.0)
            sc.doc.Objects.AddCircle(circle_Mid)

        if not cyl.IsFinite:
            xPrimitiveShape.Cylinder.extendToObjectSize(cyl, rgFace0)
        circle_Ht1 = cyl.CircleAt(cyl.Height1)
        circle_Ht2 = cyl.CircleAt(cyl.Height2)
        line_Axis = rg.Line(circle_Ht1.Center, circle_Ht2.Center)
        sc.doc.Objects.AddLine(line_Axis)
#        line_Profile = rg.Line(circle_Ht1.PointAt(0.0), circle_Ht2.PointAt(0.0))
#        sc.doc.Objects.AddLine(line_Profile)
    
    elif typeShape == rg.Cone:
        cone = rgShape
        #sc.doc.Objects.AddPoint(cone.Center)
        line_Axis = rg.Line(cone.BasePoint, cone.ApexPoint)
        if line_Axis.Length < 2.0 * sc.doc.ModelAbsoluteTolerance:
            line_Axis.Length = 1.0
            print("Line length was too short and was set to 1.0.")
        sc.doc.Objects.AddLine(line_Axis)

        if bIncludeArcCrvsForRoundShapes:
            circle = rg.Circle(cone.Plane, cone.BasePoint, cone.Radius)
            sc.doc.Objects.AddCircle(circle)
            line_Profile = rg.Line(circle.PointAt(0.0), cone.ApexPoint)
            sc.doc.Objects.AddLine(line_Profile)

    elif typeShape == rg.Sphere:
        sphere = rgShape
        sc.doc.Objects.AddPoint(sphere.Center)

        if bIncludeArcCrvsForRoundShapes:
            # sphere's EquatorialPlane may be invalid.  Just create a plane for Circle.
            circle = rg.Circle(
                plane=rg.Plane(sphere.Center, normal=rg.Vector3d.ZAxis),
                radius=sphere.Radius)
            sc.doc.Objects.AddCircle(circle)

    elif typeShape == rg.Torus:
        plane_Major = rgShape.Plane
        
        #sc.doc.Objects.AddPoint(plane_Major.Origin)
        
        if bIncludeArcCrvsForRoundShapes:
            circle_Major = rg.Circle(plane_Major, plane_Major.Origin, rgShape.MajorRadius)
            sc.doc.Objects.AddCircle(circle_Major)
            plane_Minor_Origin = (rg.Point3d(
                    plane_Major.OriginX,
                    plane_Major.OriginY,
                    plane_Major.OriginZ) +
                    rgShape.MajorRadius * rg.Point3d(plane_Major.XAxis))
        
        #sc.doc.Objects.AddPoint(plane_Minor_Origin)
        
        if bIncludeArcCrvsForRoundShapes:
            plane_Minor = rg.Plane(plane_Minor_Origin,
                    plane_Major.XAxis,
                    plane_Major.ZAxis)
            circle_Minor = rg.Circle(plane_Minor, plane_Minor_Origin, rgShape.MinorRadius)
            sc.doc.Objects.AddCircle(circle_Minor)
        
        #Add line representing axis.
        line = rg.Line(plane_Major.Origin + rgShape.MinorRadius * plane_Major.ZAxis,
                plane_Major.Origin - rgShape.MinorRadius * plane_Major.ZAxis)
        sc.doc.Objects.AddLine(line)
    
    else:
        print("What happened?")
        return
    
    sc.doc.Views.Redraw()


def main():
    """
    """
    
    rc = getInput()
    if rc is None: return
    (
        rgFace0,
        fShapeTol,
        bAddWireframe,
        bIncludeArcCrvsForRoundShapes,
        bEcho,
        bDebug,
        ) = rc

    rc = getFaceShape(
        rgFace0,
        fShapeTol=fShapeTol,
        bDebug=bDebug)
    if rc is None: return
    rgShape, fTol_PrimitiveMatch = rc


    if bEcho:
        print(createReport(rgShape, fTol_PrimitiveMatch))


    if bAddWireframe:
        addWireframeOfShape(
            rgShape,
            rgFace0,
            bIncludeArcCrvsForRoundShapes=bIncludeArcCrvsForRoundShapes,
            bDebug=bDebug)


if __name__ == '__main__': main()