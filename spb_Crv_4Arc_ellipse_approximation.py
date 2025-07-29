"""
The script creates a 4-arc approximation of an ellipse. This is useful when an
ellipse or its NURBS equivalent are not wanted or allowed, as is the case of
some CAM software.

The method used is based on https://www.researchgate.net/publication/241719740_Approximating_an_ellipse_with_four_circular_arcs
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250728-29: Created.

TODO:
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fTol_IsEllipse'; keys.append(key)
    values[key] = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    names[key] = 'IsEllipseTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True,
        limit=Rhino.RhinoMath.ZeroTolerance)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    key = 'iNumPtsToCheck'; keys.append(key)
    values[key] = 31
    riOpts[key] = ri.Custom.OptionInteger(values[key], lowerLimit=1, upperLimit=255)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPolyCrvOutput_notArcs'; keys.append(key)
    values[key] = True
    names[key] = 'OutputCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Arcs', 'Poly')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddEllipticalDeg2NurbsIfInputIsNot'; keys.append(key)
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
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
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
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        #if key == 'fSearchTol':
        #    if cls.riOpts[key].CurrentValue < 0.0:
        #        cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
        #        sc.sticky[cls.stickyKeys[key]] = cls.values[key]
        #        return

        #    sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
        #    return

        if key == 'fTol_IsEllipse':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < cls.riOpts[key].InitialValue:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def getInput(bDebug=False):
    """
    Get curves or ellipse diameters.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select elliptical-shaped curves")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    #go.GeometryAttributeFilter = (
    #    ri.Custom.GeometryAttributeFilter.

    #go.AcceptNumber(True, acceptZero=False)
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == ri.GetResult.Object 

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        #addOption('fScaleFactor')
        addOption('fTol_IsEllipse')
        addOption('iNumPtsToCheck')
        addOption('bPolyCrvOutput_notArcs')
        addOption('bAddEllipticalDeg2NurbsIfInputIsNot')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fTol_IsEllipse'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
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


def _formatRadius(radius):
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
        return _formatDistance(radius)

    return "{} {} [{} millimeters]".format(
        _formatDistance(radius),
        sc.doc.GetUnitSystemName(
            modelUnits=True,
            capitalize=False,
            singular=False,
            abbreviate=False),
        _formatDistance(radius*Rhino.RhinoMath.UnitScale(
            sc.doc.ModelUnitSystem, Rhino.UnitSystem.Millimeters),
            fPrecision=4)
        )


def create4ArcApproximation(a, b, iNumPtsToCheck):

    if a > b:
        bFlippedAB = False
    else:
        a, b = b, a
        bFlippedAB = True

    nc_RefEllipse = rg.Ellipse(rg.Plane.WorldXY, radius1=a, radius2=b).ToNurbsCurve()

    C = rg.Point3d(
        x= (a-b)/2.0,
        y=-(a-b)/2.0,
        z=0.0)

    #sEval = "C"; print(sEval,'=',eval(sEval))

    incircle = rg.Circle(center=C, radius=(a-b)/2.0)

    arc_in = rg.Arc(
        circle=incircle,
        angleIntervalRadians=rg.Interval(
            t0=Rhino.RhinoMath.ToRadians(-90.0),
            t1=0.0))
    #sc.doc.Objects.AddArc(arc_in)

    arcCrv_in = rg.ArcCurve(arc_in)
    #sc.doc.Objects.AddCurve(arcCrv_in)

    A = rg.Point3d(a, 0.0, 0.0)
    B = rg.Point3d(0.0, b, 0.0)
    E = rg.Point3d(a, b, 0.0)

    radius_conArc = (A - C).Length

    #sEval = "radius_conArc"; print(sEval,'=',eval(sEval))
    #sEval = "type(radius_conArc)"; print(sEval,'=',eval(sEval))

    concentriccircle = rg.Circle(center=C, radius=radius_conArc)

    lineBE = rg.Line(B, E)
    #sc.doc.Objects.AddLine(lineBE)

    # Alternative method to determine D.
    #rvs = rg.Intersect.Intersection.LineCircle(lineBE, concentriccircle)
    #D = rvs[4]

    D = rg.Point3d(a-b, b, 0.0)

    #sEval = "D"; print(sEval,'=',eval(sEval))

    angleIntervalRadians = rg.Interval(
       math.atan2(A.Y - C.Y, A.X - C.X),
       math.atan2(D.Y - C.Y, D.X - C.X))

    #sEval = "angleIntervalRadians"; print(sEval,'=',eval(sEval))

    arc_con = rg.Arc(
        circle=concentriccircle,
        angleIntervalRadians=angleIntervalRadians)
    #sc.doc.Objects.AddArc(arc_con)

    domainLength = arc_con.AngleDomain.Length
    domainStart = arc_con.AngleDomain.T0

    Ts = []
    devs = []
    arcCrvs_AT = []
    arcCrvs_BT = []

    for i in range(1, iNumPtsToCheck):
        sc.escape_test()

        tangle = (float(i)/float(iNumPtsToCheck+1))*domainLength + domainStart
        T = arc_con.PointAt(tangle)
        #sc.doc.Objects.AddPoint(T)

        bSuccess, t = arcCrv_in.GetLocalTangentPoint(
            testPoint=T,
            seedParmameter=arcCrv_in.Domain.Mid)

        #sEval = "t"; print(sEval,'=',eval(sEval))

        line_tan = rg.Line(T, arcCrv_in.PointAt(t))
        #sc.doc.Objects.AddLine(line)

        lineOA = rg.Line(rg.Point3d.Origin, A)
        lineOB = rg.Line(rg.Point3d.Origin, B)

        bSuccess, tOA, t_tan = rg.Intersect.Intersection.LineLine(lineOA, line_tan)

        C1 = lineOA.PointAt(tOA)
        #sEval = "C1"; print(sEval,'=',eval(sEval))

        bSuccess, tOB, t_tan = rg.Intersect.Intersection.LineLine(lineOB, line_tan)

        C2 = lineOB.PointAt(tOB)
        #sEval = "C2"; print(sEval,'=',eval(sEval))

        circle_AT = rg.Circle(center=C1, radius=(A-C1).Length)
        #sc.doc.Objects.AddCircle(circle_AT)

        angle = math.atan2(T.Y - C1.Y, T.X - C1.X)

        angleIntervalRadians = rg.Interval(-angle, angle)

        arc_AT = rg.Arc(
            circle=circle_AT,
            angleIntervalRadians=angleIntervalRadians)
        #sc.doc.Objects.AddArc(arc_AT)

        arcCrv_AT = rg.ArcCurve(arc_AT)

        circle_BT = rg.Circle(center=C2, radius=(B-C2).Length)
        #sc.doc.Objects.AddCircle(circle_BT)

        angle = math.atan2(T.Y - C2.Y, T.X - C2.X)

        angleIntervalRadians = rg.Interval(angle, math.pi-angle)
        #sEval = "angleIntervalRadians"; print(sEval,'=',eval(sEval))

        arc_BT = rg.Arc(
            circle=circle_BT,
            angleIntervalRadians=angleIntervalRadians)
        #sc.doc.Objects.AddArc(arc_BT)

        arcCrv_BT = rg.ArcCurve(arc_BT)

        rv = rg.Curve.GetDistancesBetweenCurves(arcCrv_AT, nc_RefEllipse, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
        if not rv[0]:
            raise Exception("GetDistancesBetweenCurves failed.")
        devA = rv[1]

        rv = rg.Curve.GetDistancesBetweenCurves(arcCrv_BT, nc_RefEllipse, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
        if not rv[0]:
            raise Exception("GetDistancesBetweenCurves failed.")
        devB = rv[1]

        dev = max((devA, devB))

        Ts.append(T)
        devs.append(dev)
        arcCrvs_AT.append(arcCrv_AT)
        arcCrvs_BT.append(arcCrv_BT)

        # Graph deviations.
        #sc.doc.Objects.AddLine(
        #    rg.Point3d(float(i), 0.0, 0.0),
        #    rg.Point3d(float(i), dev, 0.0))

    idxWinner = devs.index(min(devs))

    #print("idxWinner: {}".format(idxWinner))
    print("Max. dev. from ellipse: {}".format(_formatDistance(devs[idxWinner])))

    T = Ts[idxWinner]
    arcCrv_AT = arcCrvs_AT[idxWinner]
    arcCrv_BT = arcCrvs_BT[idxWinner]

    if bFlippedAB:
        xform = rg.Transform.Mirror(
            mirrorPlane=rg.Plane(
                origin=rg.Point3d.Origin,
                normal=rg.Vector3d(-1.0, 1.0, 0.0)))

        arcCrv_AT.Transform(xform)
        arcCrv_BT.Transform(xform)

    xform = rg.Transform.Rotation(
        angleRadians=math.pi,
        rotationCenter=rg.Point3d.Origin)

    arcCrv_AT_Dup = rg.ArcCurve(arcCrv_AT)

    bSuccess = arcCrv_AT_Dup.Transform(xform)

    #sc.doc.Objects.AddCurve(arcCrv_AT)
    #sc.doc.Objects.AddCurve(arcCrv_AT_Dup)

    arcCrv_BT_Dup = rg.ArcCurve(arcCrv_BT)

    bSuccess = arcCrv_BT_Dup.Transform(xform)

    #sc.doc.Objects.AddCurve(arcCrv_BT)
    #sc.doc.Objects.AddCurve(arcCrv_BT_Dup)

    if bFlippedAB:
        return arcCrv_BT, arcCrv_AT, arcCrv_BT_Dup, arcCrv_AT_Dup
    else:
        return arcCrv_AT, arcCrv_BT, arcCrv_AT_Dup, arcCrv_BT_Dup


def main():

    objrefs = getInput()
    if objrefs is None: return

    #fScaleFactor = Opts.values['fScaleFactor']
    fTol_IsEllipse = Opts.values['fTol_IsEllipse']
    iNumPtsToCheck = Opts.values['iNumPtsToCheck']
    bPolyCrvOutput_notArcs = Opts.values['bPolyCrvOutput_notArcs']
    bAddEllipticalDeg2NurbsIfInputIsNot = Opts.values['bAddEllipticalDeg2NurbsIfInputIsNot']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gOuts = []

    for objref in objrefs:
        rgC_In = objref.Curve()
        if rgC_In is None:
            raise Exception("Curve not found!")
            continue

        if isinstance(rgC_In, rg.ArcCurve):
            if bEcho and len(objrefs) == 1:
                print("Skipped ArcCurve.")
                return
            continue

        if rgC_In.IsArc():
            if bEcho and len(objrefs) == 1:
                print("Skipped arc-shaped curve.")
                return
            continue

        bSuccess, ellipse = rgC_In.TryGetEllipse()
        if bSuccess:
            bInputIsAccurateEllipse = True
        else:
            bSuccess, ellipse = rgC_In.TryGetEllipse(tolerance=fTol_IsEllipse)
            if not bSuccess:
                continue
            bInputIsAccurateEllipse = False

        if isinstance(rgC_In, rg.BrepEdge):
            bIsEdge = True
            rgC_In_NotProxy = rgC_In.DuplicateCurve()
        else:
            bIsEdge = False
            rgC_In_NotProxy = rgC_In

        bIsNurbsCrv = isinstance(rgC_In_NotProxy, rg.NurbsCurve)

        #radius = arc.Radius
        #radii.append(radius)

        a = ellipse.Radius1
        b = ellipse.Radius2

        arcCrvs_Res = create4ArcApproximation(a, b, iNumPtsToCheck)
        #sEval = "arcCrvs_Res"; print(sEval,'=',eval(sEval))

        if arcCrvs_Res is None:
            continue

        plane = ellipse.Plane
        xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, plane)

        if bPolyCrvOutput_notArcs:
            polycrv = rg.Curve.JoinCurves(arcCrvs_Res)[0]
            polycrv.Transform(xform)
            gOut = sc.doc.Objects.AddCurve(polycrv)
            if gOut != gOut.Empty:
                gOuts.append(gOut)
        else:
            for c in arcCrvs_Res:
                c.Transform(xform)
                gOut = sc.doc.Objects.AddCurve(c)
                if gOut != gOut.Empty:
                    gOuts.append(gOut)

        if bAddEllipticalDeg2NurbsIfInputIsNot:
            if not (bInputIsAccurateEllipse and bIsNurbsCrv and rgC_In_NotProxy.Degree == 2):
                gOut = sc.doc.Objects.AddEllipse(ellipse)
                if gOut != gOut.Empty:
                    gOuts.append(gOut)



    if bEcho:
        #if len(objrefs) > 1:
            #if len(radii) == 0:
            #    print("None of the curves are arc-shaped.")
            #elif len(objrefs) > 10:
            #    print("Arc-shaped curves found.")
            #else:
            #    print("Curves' radii: {}".format(
            #        ", ".join([_formatRadius(radius) for radius in radii])))

        print("Added {} curves.".format(len(gOuts)))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()