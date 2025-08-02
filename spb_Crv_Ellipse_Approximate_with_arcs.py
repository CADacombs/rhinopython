"""
The script creates a 4-arc or 8-arc approximation of an ellipse. This is useful when an
ellipse or its NURBS equivalent are not wanted or allowed, as is the case of
some CAM software.

The 4-arc method is based on https://www.researchgate.net/publication/241719740_Approximating_an_ellipse_with_four_circular_arcs
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250728-0801: Created.

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

    key = 'b8Arcs_not4'; keys.append(key)
    values[key] = True
    names[key] = 'NumOfArcsInEntireEllipse'
    riOpts[key] = ri.Custom.OptionToggle(values[key], '4', '8')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    #key = 'iNumPtsToCheck'; keys.append(key)
    #values[key] = 31
    #riOpts[key] = ri.Custom.OptionInteger(values[key], lowerLimit=1, upperLimit=255)
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPolyCrvOutput_notArcs'; keys.append(key)
    values[key] = True
    names[key] = 'OutputCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Arcs', 'Poly')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = False
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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

        addOption('fTol_IsEllipse')
        addOption('b8Arcs_not4')
        #addOption('iNumPtsToCheck')
        addOption('bPolyCrvOutput_notArcs')
        addOption('bReplace')
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


def _is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False


def _joinCrvs(curves_In):
    rv = rg.Curve.JoinCurves(curves_In)
    if len(rv) != 1:
        raise Exception("JoinCurves produced {} curves.".format(len(rv)))
    return rv[0]


def _getMaxDev(curvesA, curveB):

    if _is_iterable(curvesA):
        crv_WIP = _joinCrvs(curvesA)
    else:
        crv_WIP = curvesA.DuplicateCurve()

    if not isinstance(curveB, rg.Curve):
        curveB = curveB.ToNurbsCurve()

    rv = rg.Curve.GetDistancesBetweenCurves(crv_WIP, curveB, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)

    crv_WIP.Dispose()

    if not rv[0]:
        raise Exception("GetDistancesBetweenCurves failed.")

    return rv[1]


    devs = []

    for curveA in curvesA:
        rv = rg.Curve.GetDistancesBetweenCurves(curveA, curveB, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
        if not rv[0]:
            for c in curvesA: sc.doc.Objects.AddCurve(c)
            raise Exception("GetDistancesBetweenCurves failed.")
        devs.append(rv[1])

    return max(devs)


def create4ArcApproximation(a, b, iNumPtsToCheck=31):

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
    domainEnd = arc_con.AngleDomain.T1


    def getArcsAtT(T):
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

        return arcCrv_AT, arcCrv_BT


    ##

    # Equally-spaced method
    Ts = []
    devs = []
    arcCrvs_AT = []
    arcCrvs_BT = []


    def checkAllEquallySpacedPts():

        for i in range(1, (iNumPtsToCheck+1)):
            sc.escape_test()

            tangle = (float(i)/float(iNumPtsToCheck+1))*domainLength + domainStart
            T = arc_con.PointAt(tangle)
            #sc.doc.Objects.AddPoint(T)

            arcCrv_AT, arcCrv_BT = getArcsAtT(T)
            dev = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)


            Ts.append(T)
            devs.append(dev)
            arcCrvs_AT.append(arcCrv_AT)
            arcCrvs_BT.append(arcCrv_BT)

            # Graph deviations.
            #sc.doc.Objects.AddLine(
            #    rg.Point3d(float(i), 0.0, 0.0),
            #    rg.Point3d(float(i), dev, 0.0))

        idxWinner = devs.index(min(devs))

        print("Index {} of {} pts".format(idxWinner, len(Ts)))
        print("Max. dev. from ellipse: {}".format(_formatDistance(devs[idxWinner])))

        T = Ts[idxWinner]
        arcCrv_AT = arcCrvs_AT[idxWinner]
        arcCrv_BT = arcCrvs_BT[idxWinner]

        return arcCrv_AT, arcCrv_BT


    #arcCrv_AT, arcCrv_BT = checkAllEquallySpacedPts()

    ##


    ##

    # Binary search method

    arcCrv_con = rg.ArcCurve(arc_con)
    #sc.doc.Objects.AddCurve(arcCrv_con)

    #sEval = "arc_con.AngleDomain"; print(sEval,'=',eval(sEval))
    #sEval = "arcCrv_con.Domain"; print(sEval,'=',eval(sEval))

    t_H = domainEnd
    T_H = arc_con.PointAt(domainEnd)
    arcCrv_AT, arcCrv_BT = getArcsAtT(T_H)
    dev_H = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)

    t = arcCrv_con.DivideByLength(10.0*sc.doc.ModelAbsoluteTolerance, includeEnds=False)[0]
    T_L = arcCrv_con.PointAt(t)
    t_L = arc_con.ClosestParameter(T_L)
    arcCrv_AT, arcCrv_BT = getArcsAtT(T_L)
    dev_L = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)

    #sc.doc.Objects.AddPoint(T_L)

    i = 0

    while T_H.DistanceTo(T_L) > (0.1 * sc.doc.ModelAbsoluteTolerance):
        sc.escape_test()

        i += 1

        #print('-'*20)
        #sEval = "dev_H"; print(sEval,'=',eval(sEval))

        t_MH = 0.75*t_H + 0.25*t_L
        T_MH = arc_con.PointAt(t_MH)
        arcCrv_AT, arcCrv_BT = getArcsAtT(T_MH)
        dev_MH = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)
        #sEval = "dev_MH"; print(sEval,'=',eval(sEval))

        t_MM = 0.5*t_H + 0.5*t_L
        T_MM = arc_con.PointAt(t_MM)
        arcCrv_AT, arcCrv_BT = getArcsAtT(T_MM)
        dev_MM = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)
        #sEval = "dev_MM"; print(sEval,'=',eval(sEval))

        t_ML = 0.25*t_H + 0.75*t_L
        T_ML = arc_con.PointAt(t_ML)
        arcCrv_AT, arcCrv_BT = getArcsAtT(T_ML)
        dev_ML = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)
        #sEval = "dev_ML"; print(sEval,'=',eval(sEval))

        #sEval = "dev_L"; print(sEval,'=',eval(sEval))

        bChanged = False

        if dev_MM < dev_MH < dev_H:
            T_H = T_MH
            t_H = t_MH
            dev_H = dev_MH
            bChanged = True

        if dev_MM < dev_ML < dev_L:
            T_L = T_ML
            t_L = t_ML
            dev_L = dev_ML
            bChanged = True

        if not bChanged:
            raise Exception("No change.")

    #print("Search took {} iterations.".format(i))

    #print('-'*20)
    t_MM = 0.5*t_H + 0.5*t_L
    T_MM = arc_con.PointAt(t_MM)
    arcCrv_AT, arcCrv_BT = getArcsAtT(T_MM)

    dev_MM = _getMaxDev((arcCrv_AT, arcCrv_BT), nc_RefEllipse)
    #sEval = "dev_MM"; print(sEval,'=',eval(sEval))


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
        return (arcCrv_BT, arcCrv_AT, arcCrv_BT_Dup, arcCrv_AT_Dup), dev_MM
    else:
        return (arcCrv_AT, arcCrv_BT, arcCrv_AT_Dup, arcCrv_BT_Dup), dev_MM


def create8ArcApproximation(a, b):

    if abs(a - b) < sc.doc.ModelAbsoluteTolerance:
        return

    if a > b:
        bFlippedAB = False
    else:
        a, b = b, a
        bFlippedAB = True

    nc_RefEllipse = rg.Ellipse(rg.Plane.WorldXY, radius1=a, radius2=b).ToNurbsCurve()

    rA = (b**2)/a
    rB = (a**2)/b

    #sEval = "rA"; print(sEval,'=',eval(sEval))
    #sEval = "rB"; print(sEval,'=',eval(sEval))

    rF = 0.5*(rA + rB)

    CA = rg.Point3d(a-rA, 0.0, 0.0)
    circle_A = rg.Circle(center=CA, radius=rA)
    #sc.doc.Objects.AddCircle(circle_A)

    arc_A = rg.Arc(
        circle=circle_A,
        angleIntervalRadians=rg.Interval(
            t0=0.0,
            t1=math.pi/2.0))
    #sc.doc.Objects.AddArc(arc_A)
    arcCrv_A = rg.ArcCurve(arc_A)
    #sc.doc.Objects.AddCurve(arcCrv_A)

    CB = rg.Point3d(0.0, b-rB, 0.0)
    circle_B = rg.Circle(center=CB, radius=rB)
    #sc.doc.Objects.AddCircle(circle_B)

    arc_B = rg.Arc(
        circle=circle_B,
        angleIntervalRadians=rg.Interval(
            t0=0.0,
            t1=math.pi/2))
    #sc.doc.Objects.AddArc(arc_B)
    arcCrv_B = rg.ArcCurve(arc_B)

    if arcCrv_B.PointAtStart.X > a:
        bSuccess, t = arcCrv_B.ClosestPoint(arc_A.Center)
        if not bSuccess:
            raise Exception("ClosestPoint failed.")
        arcCrv_B = arcCrv_B.Trim(t, arcCrv_B.Domain.T1)

    #sc.doc.Objects.AddCurve(arcCrv_B)

    #arc_F = rg.Curve.CreateFillet(
    #    curve0=arcCrv_A,
    #    curve1=arcCrv_B,
    #    radius=rF,
    #    t0Base=arcCrv_A.Domain.Mid,
    #    t1Base=arcCrv_B.Domain.Mid)
    #sc.doc.Objects.AddArc(arc_F)

    A = rg.Point3d(a, 0.0, 0.0)
    B = rg.Point3d(0.0, b, 0.0)


    def createFilletCurves(radius):
        for point0 in arcCrv_A.PointAtMid, arcCrv_A.PointAtStart, arcCrv_A.PointAtEnd:
            for point1 in arcCrv_B.PointAtMid, arcCrv_B.PointAtEnd, arcCrv_B.PointAtStart:
                rv = rg.Curve.CreateFilletCurves(
                    curve0=arcCrv_A,
                    point0=point0,
                    curve1=arcCrv_B,
                    point1=point1,
                    radius=radius,
                    join=False,
                    trim=True,
                    arcExtension=True,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
                    angleTolerance=0.1*sc.doc.ModelAngleToleranceRadians)
                if len(rv) == 0:
                    raise Exception("No curves.")
                for c in rv:
                    if not c.IsValid:
                        raise Exception("Curve from CreateFilletCurves is invalid.")
                for c in rv:
                    if c.AngleDegrees > 90.0:
                        #print("Curve of {} degrees found in CreateFilletCurves with radius of {} and point0 at {} and point1 at {}.".format(
                        #    c.AngleDegrees, radius, point0, point1))
                        for c in rv:
                            c.Dispose()
                        break
                else:
                    # Success.
                    return rv

        raise Exception("CreateFilletCurves failed to provide all good curves < 90 degrees.")


    #sEval = "rv"; print(sEval,'=',eval(sEval))

    #sc.doc.Objects.AddCurve(rv[2])

    r_H = (a**2 + a*b) / (2*b)
    r_L = (b**2 + a*b) / (2*a)

    crvs_H = createFilletCurves(r_H)
    dev_H = _getMaxDev(crvs_H, nc_RefEllipse)

    crvs_L = createFilletCurves(r_L)
    dev_L = _getMaxDev(crvs_L, nc_RefEllipse)

    i = 0

    while (r_H - r_L) > 0.1 * sc.doc.ModelAbsoluteTolerance:
        sc.escape_test()

        i += 1

        r_MH = 0.75*r_H + 0.25*r_L
        crvs_MH = createFilletCurves(r_MH)
        dev_MH = _getMaxDev(crvs_MH, nc_RefEllipse)

        r_MM = 0.5*r_H + 0.5*r_L
        crvs_MM = createFilletCurves(r_MM)
        dev_MM = _getMaxDev(crvs_MM, nc_RefEllipse)

        r_ML = 0.25*r_H + 0.75*r_L
        crvs_ML = createFilletCurves(r_ML)
        dev_ML = _getMaxDev(crvs_ML, nc_RefEllipse)

        bChanged = False

        if dev_MM < dev_MH < dev_H:
            r_H = r_MH
            dev_H = dev_MH
            for c in crvs_H: c.Dispose()
            crvs_H = crvs_MH
            bChanged = True
        elif dev_MM < dev_H and dev_MM < dev_MH:
            r_H = r_MH
            dev_H = dev_MH
            for c in crvs_H: c.Dispose()
            crvs_H = crvs_MH
            bChanged = True

        if dev_MM < dev_L and dev_MM < dev_ML:
            r_L = r_ML
            dev_L = dev_ML
            for c in crvs_L: c.Dispose()
            crvs_L = crvs_ML
            bChanged = True

        if not bChanged:
            if dev_ML < dev_MH:
                r_H = r_MH
                dev_H = dev_MH
                for c in crvs_H: c.Dispose()
                crvs_H = crvs_MH
                bChanged = True
            elif dev_MH < dev_ML:
                r_L = r_ML
                dev_L = dev_ML
                for c in crvs_L: c.Dispose()
                crvs_L = crvs_ML
                bChanged = True

            if not bChanged:
                sEval = "dev_H"; print(sEval,'=',eval(sEval))
                sEval = "dev_MH"; print(sEval,'=',eval(sEval))
                sEval = "dev_MM"; print(sEval,'=',eval(sEval))
                sEval = "dev_ML"; print(sEval,'=',eval(sEval))
                sEval = "dev_L"; print(sEval,'=',eval(sEval))
                
                #for c in crvs_H: sc.doc.Objects.AddCurve(c)
                
                print("Iterations: {}".format(i))
                
                print("Eccentricy of ellipse may be too large for this script.")
                return

        for c in crvs_MH: c.Dispose()
        for c in crvs_MM: c.Dispose()
        for c in crvs_ML: c.Dispose()


    print("Search took {} iterations.".format(i))

    r_MM = 0.5*r_H + 0.5*r_L

    #sc.doc.Objects.AddCurve(arcCrv_A)
    #sc.doc.Objects.AddCurve(arcCrv_B)
    #return

    rv = createFilletCurves(r_MM)

    if len(rv) != 3:
        raise Exception("len(rv) : ".format(len(rv)))

    dev_MM = _getMaxDev(rv, nc_RefEllipse)


    #for c in rv:
    #    sc.doc.Objects.AddCurve(c)

    arcCrv_A, arcCrv_B, arcCrv_F = rv

    if bFlippedAB:
        xform = rg.Transform.Mirror(
            mirrorPlane=rg.Plane(
                origin=rg.Point3d.Origin,
                normal=rg.Vector3d(-1.0, 1.0, 0.0)))

        arcCrv_A.Transform(xform)
        arcCrv_B.Transform(xform)
        arcCrv_F.Transform(xform)



    arc_A = arcCrv_A.Arc
    arc_B = arcCrv_B.Arc

    #for c in rv:
    #    sc.doc.Objects.AddCurve(c)
    #return

    arc_A.StartAngle = -arc_A.EndAngle
    arcCrv_A = rg.ArcCurve(arc_A)
    #sc.doc.Objects.AddCurve(arcCrv_A)

    arc_B.EndAngle = math.pi - arc_B.StartAngle
    arcCrv_B = rg.ArcCurve(arc_B)
    #sc.doc.Objects.AddCurve(arcCrv_B)

    arcCrv_A_Dup = rg.ArcCurve(arcCrv_A)
    arcCrv_B_Dup = rg.ArcCurve(arcCrv_B)
    arcCrv_F_Dup1 = rg.ArcCurve(arcCrv_F)

    xform = rg.Transform.Mirror(
        mirrorPlane=rg.Plane(
            origin=rg.Point3d.Origin,
            normal=rg.Vector3d(0.0, 1.0, 0.0)))

    bSuccess = arcCrv_B_Dup.Transform(xform)
    bSuccess = arcCrv_F_Dup1.Transform(xform)

    arcCrv_F_Dup2 = rg.ArcCurve(arcCrv_F)
    arcCrv_F_Dup3 = rg.ArcCurve(arcCrv_F_Dup1)

    xform = rg.Transform.Mirror(
        mirrorPlane=rg.Plane(
            origin=rg.Point3d.Origin,
            normal=rg.Vector3d(1.0, 0.0, 0.0)))

    bSuccess = arcCrv_A_Dup.Transform(xform)
    bSuccess = arcCrv_F_Dup2.Transform(xform)
    bSuccess = arcCrv_F_Dup3.Transform(xform)

    return ((
        arcCrv_A,
        arcCrv_F,
        arcCrv_B,
        arcCrv_F_Dup2,
        arcCrv_A_Dup,
        arcCrv_F_Dup3,
        arcCrv_B_Dup,
        arcCrv_F_Dup1
        ),
            dev_MM
            )


def main():

    objrefs = getInput()
    if objrefs is None: return

    fTol_IsEllipse = Opts.values['fTol_IsEllipse']
    b8Arcs_not4 = Opts.values['b8Arcs_not4']
    #iNumPtsToCheck = Opts.values['iNumPtsToCheck']
    bPolyCrvOutput_notArcs = Opts.values['bPolyCrvOutput_notArcs']
    bAddEllipticalDeg2NurbsIfInputIsNot = Opts.values['bAddEllipticalDeg2NurbsIfInputIsNot']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    fDevs = []
    gAdded = []
    gReplaced = []

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

        if b8Arcs_not4:
            rv = create8ArcApproximation(a, b)
        else:
            rv = create4ArcApproximation(a, b)

        if rv is None:
            continue

        arcCrvs_Res, fDev = rv
        fDevs.append(fDev)

        xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, ellipse.Plane)

        if bPolyCrvOutput_notArcs:
            polycrv = rg.Curve.JoinCurves(arcCrvs_Res)[0]
            polycrv.Transform(xform)
            if bReplace and not bIsEdge:
                if sc.doc.Objects.Replace(objref, curve=polycrv):
                    gReplaced.append(objref.ObjectId)
                else:
                    print("Could not replace curve, {}.".format(objref.ObjectId))
            else:
                gOut = sc.doc.Objects.AddCurve(polycrv)
                if gOut != gOut.Empty:
                    gAdded.append(gOut)
                else:
                    print("Could not add curve approximation of {}.".format(objref.ObjectId))
        else:
            bFailedToAddSomeArcs = False
            for c in arcCrvs_Res:
                c.Transform(xform)
                gOut = sc.doc.Objects.AddCurve(c)
                if gOut != gOut.Empty:
                    gAdded.append(gOut)
                else:
                    bFailedToAddSomeArcs = True
                    print("Could not add ArcCurve for approximation of {}.".format(objref.ObjectId))
                if bReplace and not bIsEdge:
                    if not bFailedToAddSomeArcs:
                        if sc.doc.Objects.Delete(objref, quiet=False):
                            gReplaced.append(objref.ObjectId)

        if bAddEllipticalDeg2NurbsIfInputIsNot:
            if not (bInputIsAccurateEllipse and bIsNurbsCrv and rgC_In_NotProxy.Degree == 2):
                gOut = sc.doc.Objects.AddEllipse(ellipse)
                if gOut != gOut.Empty:
                    gAdded.append(gOut)

    if bEcho:
        ss = []
        if gAdded:
            ss.append("Added {} curves.".format(len(gAdded)))
        if gReplaced:
            ss.append("Replaced {} curves.".format(len(gReplaced)))
        if fDevs:
            ss.append("{} maximum deviation.".format(_formatDistance(max(fDevs))))

        print(*ss)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()