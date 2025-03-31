"""
This script will create reference geometry for Rev- and SumSurfaces.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250329-30: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math


MYZERO = 1e-6 * Rhino.RhinoMath.UnitScale(
    Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAddSumCrvA'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddSumCrvB'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddRevAxis'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddRevCrv'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddArcPerRevAngle'; keys.append(key)
    values[key] = True
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

        #if key == 'fShapeTol':
        #    if cls.riOpts[key].CurrentValue < 0.0:
        #        cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
        #    elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
        #        print("TryGet(shape) methods appear to use ZeroTolerance" \
        #            " as the minimum when less is input.")
        #        cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

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

    go.SetCommandPrompt("Select face (RevSrfs & SumSrfs are filtered.)")

    go.GeometryFilter = rd.ObjectType.Surface


    def customGeomFilter(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index)

        if not isinstance(geom, rg.BrepFace):
            return False

        rgS = geom.UnderlyingSurface()

        return isinstance(rgS, (rg.RevSurface, rg.SumSurface))


    go.SetCustomGeometryFilter(customGeomFilter)

    #go.AcceptNumber(True, acceptZero=True)
    go.EnableHighlight(False)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)


    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bAddSumCrvA')
        addOption('bAddSumCrvB')
        addOption('bAddRevAxis')
        addOption('bAddRevCrv')
        addOption('bAddArcPerRevAngle')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose
            return objref

        #if res == ri.GetResult.Number:
        #    key = 'fShapeTol'
        #    Opts.riOpts[key].CurrentValue = go.Number()
        #    Opts.setValue(key)
        #    continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _formatDistance(fDistance, fPrecision=None):
    if fDistance is None:
        return "(None)"
    if fDistance == Rhino.RhinoMath.UnsetValue:
        return "(Infinite)"
    if fPrecision is None:
        fPrecision = sc.doc.ModelDistanceDisplayPrecision

    if fDistance < 10.0**(-(fPrecision-1)):
        # For example, if fDistance is 1e-5 and fPrecision == 5,
        # the end of this script would display only one digit.
        # Instead, this return displays 2 digits.
        return "{:.2e}".format(fDistance)

    return "{:.{}f}".format(fDistance, fPrecision)


class defsrfgeom:
    revAngle = None
    revAxis = None
    revCrv = None
    sumCrvA = None
    sumCrvB = None


def getDefiningGeomOfSrf(rgSrf, bDebug=False):
    """
    """

    if isinstance(rgSrf, rg.RevSurface):
        defsrfgeom.revAngle = rgSrf.Angle
        defsrfgeom.revAxis = rgSrf.Axis
        defsrfgeom.revCrv = rgSrf.Curve
        return True
        #return (
        #    rgSrf.Angle,
        #    rgSrf.Axis,
        #    rgSrf.Curve,
        #    )

    if isinstance(rgSrf, rg.SumSurface):
        defsrfgeom.sumCrvA = rgSrf.IsoCurve(direction=1, constantParameter=rgSrf.Domain(0).Min)
        defsrfgeom.sumCrvB = rgSrf.IsoCurve(direction=0, constantParameter=rgSrf.Domain(1).Min)
        return True
        #return (
        #    rgSrf.IsoCurve(direction=1, constantParameter=rgSrf.Domain(0).Min),
        #    rgSrf.IsoCurve(direction=0, constantParameter=rgSrf.Domain(1).Min),
        #    )

    print("{} is not supported.".format(rgSrf.GetType().Name))

    return False


def createReport():
    """
    """

    def decimalPlacesForZero():
        return int(abs(math.log10(abs(MYZERO)))) + 1

    def myround(x):
        return round(x, decimalPlacesForZero())

    ss = []

    if defsrfgeom.revAxis:
        line = defsrfgeom.revAxis
        ss.append("RevSurface axis vector:{},{},{}, start:{},{},{}, end:{},{},{}".format(
            myround(line.ToX - line.FromX),
            myround(line.ToY - line.FromY),
            myround(line.ToZ - line.FromZ),
            myround(line.FromX),
            myround(line.FromY),
            myround(line.FromZ),
            myround(line.ToX),
            myround(line.ToY),
            myround(line.ToZ),
            ))

    if defsrfgeom.revCrv:
        ss.append("RevSurface curve type:{}".format(
            defsrfgeom.revCrv.GetType().Name))

    if defsrfgeom.revAngle:
        ss.append("RevSurface angle interval:[{},{}]radians [{},{}]degrees".format(
            myround(defsrfgeom.revAngle.T0),
            myround(defsrfgeom.revAngle.T1),
            myround(Rhino.RhinoMath.ToDegrees(defsrfgeom.revAngle.T0)),
            myround(Rhino.RhinoMath.ToDegrees(defsrfgeom.revAngle.T1)),
            ))

    if defsrfgeom.sumCrvA:
        ss.append("SumSurface curveA type:{}".format(
            defsrfgeom.sumCrvA.GetType().Name))

    if defsrfgeom.sumCrvB:
        ss.append("SumSurface curveA type:{}".format(
            defsrfgeom.sumCrvB.GetType().Name))

    return "\n".join(ss)


def addDocObjects(bAddSumCrvA=False, bAddSumCrvB=False, bAddRevAxis=False, bAddRevCrv=False, bAddArcPerRevAngle=False, bDebug=False):
    """
    """

    if bAddSumCrvA and defsrfgeom.sumCrvA:
        sc.doc.Objects.AddCurve(defsrfgeom.sumCrvA)

    if bAddSumCrvB and defsrfgeom.sumCrvB:
        sc.doc.Objects.AddCurve(defsrfgeom.sumCrvB)

    if bAddRevAxis and defsrfgeom.revAxis:
        sc.doc.Objects.AddLine(defsrfgeom.revAxis)

    if bAddRevCrv and defsrfgeom.revCrv:
        sc.doc.Objects.AddCurve(defsrfgeom.revCrv)

    return
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
        if not bAddRevCrv:
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

        if not bAddRevCrv:
            circle = rg.Circle(cone.Plane, cone.BasePoint, cone.Radius)
            sc.doc.Objects.AddCircle(circle)
            line_Profile = rg.Line(circle.PointAt(0.0), cone.ApexPoint)
            sc.doc.Objects.AddLine(line_Profile)

    elif typeShape == rg.Sphere:
        sphere = rgShape
        sc.doc.Objects.AddPoint(sphere.Center)

        if not bAddRevAxis:
            # sphere's EquatorialPlane may be invalid.  Just create a plane for Circle.
            circle = rg.Circle(
                plane=rg.Plane(sphere.Center, normal=rg.Vector3d.ZAxis),
                radius=sphere.Radius)
            sc.doc.Objects.AddCircle(circle)

    elif typeShape == rg.Torus:
        plane_Major = rgShape.Plane
        
        #sc.doc.Objects.AddPoint(plane_Major.Origin)
        
        if not bAddRevCrv:
            circle_Major = rg.Circle(plane_Major, plane_Major.Origin, rgShape.MajorRadius)
            sc.doc.Objects.AddCircle(circle_Major)
            plane_Minor_Origin = (rg.Point3d(
                    plane_Major.OriginX,
                    plane_Major.OriginY,
                    plane_Major.OriginZ) +
                    rgShape.MajorRadius * rg.Point3d(plane_Major.XAxis))
        
        #sc.doc.Objects.AddPoint(plane_Minor_Origin)
        
        if not bAddRevCrv:
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

    objref = getInput()
    if objref is None: return
    rgF_In = objref.Face()
    if rgF_In is None:
        print("Face not obtained.")
        return

    rgS = rgF_In.UnderlyingSurface()

    bAddSumCrvA = Opts.values['bAddSumCrvA']
    bAddSumCrvB = Opts.values['bAddSumCrvB']
    bAddRevAxis = Opts.values['bAddRevAxis']
    bAddRevCrv = Opts.values['bAddRevCrv']
    bAddArcPerRevAngle = Opts.values['bAddArcPerRevAngle']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    bSuccess = getDefiningGeomOfSrf(
        rgS,
        bDebug=bDebug)
    if not bSuccess: return


    if bEcho:
        sReport = createReport()
        if sReport:
            print(sReport)

    addDocObjects(
        bAddSumCrvA=bAddSumCrvA,
        bAddSumCrvB=bAddSumCrvB,
        bAddRevAxis=bAddRevAxis,
        bAddRevCrv=bAddRevCrv,
        bAddArcPerRevAngle=bAddArcPerRevAngle,
        bDebug=bDebug)


if __name__ == '__main__': main()