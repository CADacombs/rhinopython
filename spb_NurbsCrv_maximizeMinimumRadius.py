"""
200112-14: Created from another script.
200128: Import-related update.
200422: Modified pre- and post-object selection behavior.
200504-05: Added function to adjust weight of degree 2 Bezier curves.  Refactored various code.
        Added bSetMinTargetRad, fRadius_AcceptableMin, and bAdjWeightInDeg2Bezier.
210317: Bug fix.
210725-27: Improvements in Degree 2 Bezier weight adjustment routine.
210814-16: Improvements in tangent control point routines.
210908,12: Major improvements in tangent control point routines for single spanned curves.
221122: Import-related update.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System.Diagnostics import Stopwatch

import spb_Crv_inflections

if Rhino.RhinoApp.ExeVersion < 7:
    import spb_Crv_radiusMinima


stopwatch = Stopwatch() # One instance will be used for all tests.


class Opts():
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'bSetMinTargetRad'; keys.append(key)
    values[key] = False
    names[key] = 'SetMinTargetRadius'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fRadius_AcceptableMin'; keys.append(key)
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 0.15
    else:
        values[key] = 4.0 * Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters, to=sc.doc.ModelUnitSystem)
    names[key] = 'Radius'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bLimitCrvDev'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'Dev'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAdjWeightInDeg2Bezier'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve and parameter with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    
    #go.DisablePreSelect()
    
    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    while True:
        Opts.riAddOpts['bSetMinTargetRad'](go)
        if Opts.values['bSetMinTargetRad']: Opts.riAddOpts['fRadius_AcceptableMin'](go)
        Opts.riAddOpts['bLimitCrvDev'](go)
        if Opts.values['bLimitCrvDev']: Opts.riAddOpts['fDevTol'](go)
        Opts.riAddOpts['bAdjWeightInDeg2Bezier'](go)
        Opts.riAddOpts['bReplace'](go)
        Opts.riAddOpts['bEcho'](go)
        Opts.riAddOpts['bDebug'](go)
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return
        else:
            # An option was selected or a number was entered.

            if res == ri.GetResult.Number:
                Opts.riOpts['fDevTol'].CurrentValue = go.Number()

            if Opts.riOpts['bLimitCrvDev'].CurrentValue:
                key = 'fDevTol'
                if Opts.riOpts[key].CurrentValue < 0.0:
                    Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()


def getMaximumDeviation(rgCrvA, rgCrvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            rgCrvA,
            rgCrvB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
    if rc[0]:
        return rc[1]


def formatDistance(fDistance):
    if fDistance is None: 
        return "(No deviation provided)"
    if fDistance == 0.0:
        return 0.0
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def createNcsWithAdjTanCpSpreadForMaxMinRad_OLD(nc0, bT1WorkEnd, fDevTol=None, nc_forDevComp=None, bDebug=False):
    """
    returns list of new NurbsCurves with minimum radii larger than the original.
    """

    idxCp_Pos_A = (nc0.Points.Count - 1) if bT1WorkEnd else 0
    idxCp_Tan_A = (nc0.Points.Count - 2) if bT1WorkEnd else 1

    fRadius_Min_In = getMinimumRadius(nc0)

    if fRadius_Min_In is None: return
    if bDebug: sEval='fRadius_Min_In'; print sEval+': ',eval(sEval)

    ncs_WithLargerMinRadii = []
    
    # Try longer and shorter tangent control point spreads.
    for scaleIncr in 0.1, -0.1:
        nc_A_Pre = nc0.Duplicate()
        scale = 1.0 + scaleIncr

        fRadius_MaxMin_ThisScaleDir = None

        while True:
            sc.escape_test()
            nc_WIP = nc0.Duplicate()
            if bDebug: sEval='scale'; print sEval+': ',eval(sEval),
            xform = rg.Transform.Scale(
                    anchor=nc_WIP.Points[idxCp_Pos_A].Location,
                    scaleFactor=scale)
            pt_Target = nc_WIP.Points[idxCp_Tan_A].Location
            pt_Target.Transform(xform)
            nc_WIP.Points[idxCp_Tan_A] = pt_Target
            fRadius_Min_WIP = getMinimumRadius(nc_WIP)
            if bDebug: sEval='fRadius_Min_WIP'; print sEval+': ',eval(sEval),
            if fRadius_Min_WIP <= fRadius_Min_In:
                if bDebug: sEval='fRadius_Min_WIP <= fRadius_Min_In'; print sEval+': ',eval(sEval)
                # Minimum radius is not increasing,
                # so stop and capture the previous curve.
                break

            if fRadius_MaxMin_ThisScaleDir is not None and fRadius_Min_WIP <= fRadius_MaxMin_ThisScaleDir:
                if bDebug: sEval='fRadius_Min_WIP <= fRadius_MaxMin_ThisScaleDir'; print sEval+': ',eval(sEval)
                # Minimum radius is now decreasing from best found in this scale direction,
                # so stop and capture the previous curve.
                break

            if fDevTol is not None and nc_forDevComp is not None:
                dev = getMaximumDeviation(nc_forDevComp, nc_WIP)
                if bDebug:
                    sEval='dev'; print sEval+': ',eval(sEval),
                    sEval='dev > fDevTol'; print sEval+': ',eval(sEval),
                if dev > fDevTol:
                    # Deviation is out of tolerance,
                    # so capture the previous curve.
                    break
            elif bDebug: print

            nc_A_Pre.Dispose()
            nc_A_Pre = nc_WIP
            scale += scaleIncr
            fRadius_MaxMin_ThisScaleDir = fRadius_Min_WIP

        nc_WIP.Dispose()

        if scale != 1.0 + scaleIncr:
            # Scale is of nc_A_Pre.
            ncs_WithLargerMinRadii.append(nc_A_Pre)
        else:
            # Fail had occurred at first scale in this direction.
            nc_A_Pre.Dispose()
    
    return ncs_WithLargerMinRadii


def getMinimumRadius(nc, epsilon=1e-6, bDebug=False):
    """
    Returns:
        float(maximum minimum radius)
        None for no minimum or for unacceptable curve shapes.
    """


    if Rhino.RhinoApp.ExeVersion < 7:
        return spb_Crv_radiusMinima.getMinimumRadius(nc_In)


    pts = nc.MaxCurvaturePoints() # "An array of points if successful, null if not successful or on error." from https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_MaxCurvaturePoints.htm
    if pts is None or len(pts) == 0:
        bIsArc, arc = nc.TryGetArc(tolerance=epsilon)
        if bIsArc:
            return arc.Radius
        #sc.doc.Objects.AddCurve(nc)
        raise ValueError("Why?")
    #        for pt in pts:
    #            sc.doc.Objects.AddPoint(pt)

    #if len(pts) > 1 and nc.Degree == 2 and nc.SpanCount==1:
    #    if bDebug:
    #        print "More than 1 minimum point for Degree 2 Bezier.  Curve may be near linear."
    #        #sc.doc.Objects.AddCurve(nc)
    #        [sc.doc.Objects.AddPoint(pt) for pt in pts]
    #        sc.doc.Views.Redraw(); 1/0
    #    return None

    ts = [nc.ClosestPoint(pt)[1] for pt in pts]

    min_rad = None
    t_min_rad = None

    for t in ts:
        curvature = nc.CurvatureAt(t).Length
        if curvature <= Rhino.RhinoMath.ZeroTolerance:
            raise Exception("curvature <= Rhino.RhinoMath.ZeroTolerance")
            return None

        rad = 1.0 / nc.CurvatureAt(t).Length
        
        if min_rad is None or rad < min_rad:
            min_rad = rad
            t_min_rad = t

    return min_rad


def getMinimaRadii_V7(nc, epsilon=1e-8, bDebug=False):
    """
    Returns:
        list(float(maximum minimum radius))
        None for unacceptable curve shapes.
    """
    pts = nc.MaxCurvaturePoints() # "An array of points if successful, null if not successful or on error." from https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_MaxCurvaturePoints.htm
    if pts is None or len(pts) == 0:
        bIsArc, arc = nc.TryGetArc(tolerance=epsilon)
        if bIsArc:
            return arc.Radius
        #sc.doc.Objects.AddCurve(nc)
        raise ValueError("Why?")
    #        for pt in pts:
    #            sc.doc.Objects.AddPoint(pt)

    if len(pts) > 1 and nc.Degree == 2 and nc.SpanCount==1:
        if bDebug:
            print "More than 1 minimum point for Degree 2 Bezier.  Curve may be near linear."
            #sc.doc.Objects.AddCurve(nc)
            [sc.doc.Objects.AddPoint(pt) for pt in pts]
            sc.doc.Views.Redraw(); 1/0
        return None, None

    ts = [nc.ClosestPoint(pt)[1] for pt in pts]

    rads = []

    for t in ts:
        curvature = nc.CurvatureAt(t).Length
        if curvature <= Rhino.RhinoMath.ZeroTolerance:
            raise Exception("curvature <= Rhino.RhinoMath.ZeroTolerance")

        rad = 1.0 / nc.CurvatureAt(t).Length

        rads.append(rad)

    return rads


def adjustWeight(nc_In, fRadius_AcceptableMin=None, fDevTol=None, nc_forDevComp=None, epsilon=1e-8, bDebug=False):
    """
    Returns:
        Success: (rg.NurbsCurve(New), float(minRadius), float(deviation)), None
        Fail: None, sLog
    """

    if nc_In.Degree != 2 or nc_In.Points.Count > 3:
        return None, "Not a degree 2 Bezier."
    
    if nc_In.IsArc(tolerance=epsilon):
        return (
            None,
            "NurbsCurve is already arc-shaped (within {})," \
                " and thus has a maximum minimum radius.".format(epsilon))

    idxCp_Pos_A = 0
    idxCp_Tan_A = 1

    fRadius_Min_In = getMinimumRadius(nc_In)

    if fRadius_Min_In is None: return None, "Minimum radius could not be obtained."
    if bDebug: sEval='fRadius_Min_In'; print sEval+': ',eval(sEval)
    
    if fRadius_AcceptableMin and fRadius_Min_In >= fRadius_AcceptableMin:
        return None, "Curve's minimum radius is already >= {}.".format(
            fRadius_AcceptableMin)

    ncs_Winner_perScaleDir = [None, None]
    fRadii_MaxMin_perScaleDir = [None, None]
    fDev_perScaleDir = [None, None]
    dev = None


    def setWeight(nc, weight):
        """
        This modifies input NurbsCurve.
        Do not use NurbsCurvePointList.SetWeight instead since it relocates the point in R3.
        https://discourse.mcneel.com/t/nurbs-definition/13114
        """
        nc.Points[1] = rg.ControlPoint(nc.Points[1].Location, weight)
        if nc.Points[1].Weight != weight:
            raise ValueError("Weight was not set!")


    def doesCurveHaveMaximumMinimumRadius(nc):

        if getMinimumRadius(nc) is None:
            return False

        nc_WIP_L = nc.Duplicate()
        nc_WIP_M = nc.Duplicate()
        nc_WIP_R = nc.Duplicate()

        weight_toTry_L = None # This line keeps the variables in a particular order in the Rhino Python Debugger.
        weight_toTry_M = nc.Points[1].Weight

        weight_delta = 1e-6 #Rhino.RhinoMath.ZeroTolerance
        
        while True:
            sc.escape_test()
            weight_toTry_L = weight_toTry_M - weight_delta
            weight_toTry_R = weight_toTry_M + weight_delta

            # 
            nc_WIP_L.Points[1] = rg.ControlPoint(nc_WIP_L.Points[1].Location, weight=weight_toTry_L)
            nc_WIP_M.Points[1] = rg.ControlPoint(nc_WIP_M.Points[1].Location, weight=weight_toTry_M)
            nc_WIP_R.Points[1] = rg.ControlPoint(nc_WIP_R.Points[1].Location, weight=weight_toTry_R)

            fRad_Min_WIP_L = getMinimumRadius(nc_WIP_L)
            fRad_Min_WIP_M = getMinimumRadius(nc_WIP_M)
            fRad_Min_WIP_R = getMinimumRadius(nc_WIP_R)

            if bDebug:
                print "Weight delta: {}".format(weight_delta)
                print "W:", weight_toTry_L, weight_toTry_M, weight_toTry_R
                print "R:", fRad_Min_WIP_L, fRad_Min_WIP_M, fRad_Min_WIP_R

            if (
                (fRad_Min_WIP_L is not None and fRad_Min_WIP_L < fRad_Min_WIP_M) and
                (fRad_Min_WIP_R is not None and fRad_Min_WIP_R < fRad_Min_WIP_M)
            ):
                return True

            if (
                fRad_Min_WIP_L is not None and fRad_Min_WIP_M is not None and
                not abs(fRad_Min_WIP_L - fRad_Min_WIP_M) <= epsilon
            ):
                return False

            if (fRad_Min_WIP_R is not None and fRad_Min_WIP_M is not None and
                not abs(fRad_Min_WIP_R - fRad_Min_WIP_M) <= epsilon
            ):
                return False

            weight_delta *= 10.0


    if doesCurveHaveMaximumMinimumRadius(nc_In):
        return None, "Weight is already set for maximum minimum radius."


    def doEitherEndHaveRadius(nc, radius):
        if areAbsEqual(radius, 1.0 / nc.CurvatureAt(nc.Domain.T0).Length):
            return True
        if areAbsEqual(radius, 1.0 / nc.CurvatureAt(nc.Domain.T1).Length):
            return True
        return False


    def find_left_edge_weight(nc, starting_weight):
        """
        wL will be first weight resulting in 2 minima radii.
        """

        nc_WIP = nc.DuplicateCurve()

        wL = 10.0
        while True:
            sc.escape_test()
            if wL <= 0.1:
                wL = 0.1
                break
            setWeight(nc_WIP, wL)
            if len(getMinimaRadii_V7(nc_WIP)) > 1:
                break
            wL /= 2.0

        wR = starting_weight

        while True:
            sc.escape_test()
            wM = 0.5*wL + 0.5*wR
            if areRelEqual(wM, wL):
                setWeight(nc_WIP, wL)
                rL = getMinimumRadius(nc_WIP)
                if bDebug: print wL, rL
                #if rL is not None:
                #    raise Exception("rL should be None!")

                setWeight(nc_WIP, wM)
                rM = getMinimumRadius(nc_WIP)
                if bDebug: print wM, rM
                if rM is not None:
                    return wM
                wR = wM
                while True:
                    sc.escape_test()
                    wR += epsilon
                    setWeight(nc_WIP, wR)
                    rR = getMinimumRadius(nc_WIP)
                    if bDebug: print wR, rR
                    if rR is not None:
                        nc_WIP.Dispose()
                        return wR

            setWeight(nc_WIP, wM)
            if getMinimumRadius(nc_WIP) is None:
                wL = wM
            else:
                wR = wM


    def binarySearch(nc):

        # Determine direction from current weight to search.

        nc_WIP_L = nc.Duplicate()
        nc_WIP_M = nc.Duplicate()
        nc_WIP_R = nc.Duplicate()

        good_weights = []
        good_rads = []

        rad = getMinimumRadius(nc)

        if rad is not None:
            good_weights.append(nc.Points[1].Weight)
            good_rads.append(rad)

        for p in range(-3,7):
            weight_toTry_M = 10.0**p
            nc_WIP_M.Points[1] = rg.ControlPoint(nc_WIP_M.Points[1].Location, weight=weight_toTry_M)
            rad = getMinimumRadius(nc_WIP_M)
            if rad is not None and rad >= sc.doc.ModelAbsoluteTolerance:
                good_weights.append(weight_toTry_M)
                good_rads.append(rad)

        fRad_Min_WIP_M = max(good_rads)
        weight_toTry_M = good_weights[good_rads.index(max(good_rads))]

        weight_toTry_L = find_left_edge_weight(nc, starting_weight=weight_toTry_M)

        weight_toTry_R = 10.0


        setWeight(nc_WIP_L, weight_toTry_L)
        fRad_Min_WIP_L = getMinimumRadius(nc_WIP_L)

        setWeight(nc_WIP_R, weight_toTry_R)
        fRad_Min_WIP_R = getMinimumRadius(nc_WIP_R)

        if fRad_Min_WIP_L > fRad_Min_WIP_M:
            weight_toTry_M = weight_toTry_L + epsilon
            setWeight(nc_WIP_M, weight_toTry_M)
            fRad_Min_WIP_M = getMinimumRadius(nc_WIP_M)
            dev = getMaximumDeviation(nc_forDevComp, nc_WIP_L)
            if fRad_Min_WIP_L > fRad_Min_WIP_M:
                return (nc_WIP_L, fRad_Min_WIP_L, dev), None
            weight_toTry_M = 0.5*weight_toTry_L + 0.5*weight_toTry_R

        if fRad_Min_WIP_R > fRad_Min_WIP_M:
            raise ValueError("R radius should not be greater than M.")


        i = 0

        while True:
            sc.escape_test()
            
            # Do not use SetWeight.  https://discourse.mcneel.com/t/nurbs-definition/13114
            nc_WIP_L.Points[1] = rg.ControlPoint(nc_WIP_L.Points[1].Location, weight=weight_toTry_L)
            nc_WIP_M.Points[1] = rg.ControlPoint(nc_WIP_M.Points[1].Location, weight=weight_toTry_M)
            nc_WIP_R.Points[1] = rg.ControlPoint(nc_WIP_R.Points[1].Location, weight=weight_toTry_R)


            fRad_Min_WIP_L = getMinimumRadius(nc_WIP_L)
            fRad_Min_WIP_M = getMinimumRadius(nc_WIP_M)
            fRad_Min_WIP_R = getMinimumRadius(nc_WIP_R)

            if bDebug:
                print "i{} W:".format(i), weight_toTry_L, weight_toTry_M, weight_toTry_R
                print "i{} R:".format(i), fRad_Min_WIP_L, fRad_Min_WIP_M, fRad_Min_WIP_R

            if areRelEqual(weight_toTry_L, weight_toTry_R):
                dev = getMaximumDeviation(nc_forDevComp, nc_WIP_M)
                nc_WIP_L.Dispose()
                nc_WIP_R.Dispose()
                return (nc_WIP_M, fRad_Min_WIP_M, dev), None

            if fRad_Min_WIP_L is not None:
                if areAbsEqual(fRad_Min_WIP_L, fRad_Min_WIP_R):
                    if areAbsEqual(fRad_Min_WIP_M, fRad_Min_WIP_R):
                        dev = getMaximumDeviation(nc_forDevComp, nc_WIP_M)
                        nc_WIP_L.Dispose()
                        nc_WIP_R.Dispose()
                        return (nc_WIP_M, fRad_Min_WIP_M, dev), None
            
            if fRad_Min_WIP_M > fRad_Min_WIP_L and fRad_Min_WIP_M > fRad_Min_WIP_R:
                if fRad_Min_WIP_L > fRad_Min_WIP_R:
                    weight_toTry_R = 0.75*weight_toTry_R + 0.25*weight_toTry_M
                elif fRad_Min_WIP_L < fRad_Min_WIP_R:
                    weight_toTry_L = 0.75*weight_toTry_L + 0.25*weight_toTry_M
                else:
                    raise ValueError("asdf")
                #weight_toTry_L = 0.75*weight_toTry_L + 0.25*weight_toTry_M
                #weight_toTry_R = 0.75*weight_toTry_R + 0.25*weight_toTry_M
            elif fRad_Min_WIP_L > fRad_Min_WIP_R and fRad_Min_WIP_M > fRad_Min_WIP_R:
                weight_toTry_R = weight_toTry_M
            elif fRad_Min_WIP_R > fRad_Min_WIP_L and fRad_Min_WIP_M > fRad_Min_WIP_L:
                weight_toTry_L = weight_toTry_M
            else:
                print "Is this possible?:"
                sEval='fRadius_Min_WIP_L'; print '  '+sEval+': ',eval(sEval)
                sEval='fRadius_Min_WIP_C'; print '  '+sEval+': ',eval(sEval)
                sEval='fRadius_Min_WIP_R'; print '  '+sEval+': ',eval(sEval)
            weight_toTry_M = 0.5*weight_toTry_L + 0.5*weight_toTry_R
            
            i += 1


    rc = binarySearch(nc_In)

    if rc[0] is None:
        return rc

    if fRadius_AcceptableMin is None and fDevTol is None:
        return rc

    nc_Ret, fRadius_Min_Res, fDev_Res = rc[0]

    if fRadius_AcceptableMin is None and fDevTol is not None:
        if fDev_Res <= fDevTol:
            return rc

        # Will try code below.

    elif fRadius_AcceptableMin is not None and fDevTol is None:
        if fRadius_Min_Res >= fRadius_AcceptableMin:
            return rc
        else:
            return None, "Minimum radius cannot be achieved."

    elif fRadius_AcceptableMin is not None and fDevTol is not None:
        sLogs = []
        if fRadius_Min_Res < fRadius_AcceptableMin:
            return None, "Minimum radius cannot be achieved."
        if fDev_Res <= fDevTol:
            return rc

        # Will try code below.


    # TODO: Rewrite the following code by testing weights between
    # starting weight of nc_In and weight of nc_Ret.

    # Try weights less than and greater than starting weight.
    for iDir, weightDelta in enumerate((-0.01, 0.01)):
        weight_toTry = nc_In.Points[1].Weight + weightDelta

        while True:
            sc.escape_test()
            nc_WIP = nc_In.Duplicate()

            if bDebug:
                sEval='weight_toTry'; print sEval+': ',eval(sEval)

            nc_WIP.Points[1] = rg.ControlPoint(nc_WIP.Points[1].Location, weight=weight_toTry)

            if bDebug: stopwatch.Restart()
            fRadius_Min_WIP = getMinimumRadius(nc_WIP)
            if bDebug:
                stopwatch.Stop()
                timeElapsed = stopwatch.Elapsed.TotalSeconds
                s  = "{:.2f} seconds for ".format(timeElapsed)
                s += "getMinimumRadius"
                print s
                sEval='fRadius_Min_WIP'; print '  '+sEval+': ',eval(sEval),

            if fRadius_Min_WIP <= fRadius_Min_In:
                if bDebug:
                    sEval='fRadius_Min_WIP <= fRadius_Min_In'; print sEval+': ',eval(sEval)
                    print "  Minimum radius is not increasing," \
                          " so break."
                break

            if fDevTol is not None and nc_forDevComp is not None:
                if bDebug: stopwatch.Restart()
                dev = getMaximumDeviation(nc_forDevComp, nc_WIP)
                if bDebug:
                    stopwatch.Stop()
                    timeElapsed = stopwatch.Elapsed.TotalSeconds
                    s  = "{:.2f} seconds for ".format(timeElapsed)
                    s += "getMaximumDeviation"
                    print s
                    sEval='dev'; print sEval+': ',eval(sEval),
                    sEval='dev > fDevTol'; print sEval+': ',eval(sEval)
                if dev > fDevTol:
                    if bDebug:
                        print "  Deviation is out of tolerance, so break."
                    break

            if fRadius_AcceptableMin and fRadius_Min_WIP >= fRadius_AcceptableMin:
                print "Curve's minimum radius, {} is >= {}.".format(
                    fRadius_Min_WIP, fRadius_AcceptableMin)
                return (nc_WIP, fRadius_Min_WIP, fDevTol), None
        
            if bDebug:
                print "Curve minimum radius is larger than that of the starting curve."

            if ncs_Winner_perScaleDir[iDir] is None:
                ncs_Winner_perScaleDir[iDir] = nc_WIP
                fRadii_MaxMin_perScaleDir[iDir] = fRadius_Min_WIP
                fDev_perScaleDir[iDir] = dev
                
                weight_toTry += weightDelta
                continue

            if fRadius_Min_WIP > fRadii_MaxMin_perScaleDir[iDir]:
                if bDebug:
                    sEval='fRadius_Min_WIP > fRadii_MaxMin_perScaleDir[iDir]'; print sEval+': ',eval(sEval)
                    print "Minimum radius is increasing."
                ncs_Winner_perScaleDir[iDir].Dispose()
                
                ncs_Winner_perScaleDir[iDir] = nc_WIP
                fRadii_MaxMin_perScaleDir[iDir] = fRadius_Min_WIP
                fDev_perScaleDir[iDir] = dev
                
                weight_toTry += weightDelta
                continue


            if bDebug:
                print "Curve minimum radius is no longer increasing, so break"

            nc_WIP.Dispose()

            break

    if not (ncs_Winner_perScaleDir[0] or ncs_Winner_perScaleDir[1]):
        return None, "No curves with increased minimum radius found at given parameters."

    if ncs_Winner_perScaleDir[0] and ncs_Winner_perScaleDir[1]:
        iDir_Winner = fRadii_MaxMin_perScaleDir.index(max(fRadii_MaxMin_perScaleDir))
    else:
        iDir_Winner = 0 if ncs_Winner_perScaleDir[0] else 1

    return (
        (
            ncs_Winner_perScaleDir[iDir_Winner],
            fRadii_MaxMin_perScaleDir[iDir_Winner],
            fDev_perScaleDir[iDir_Winner]
        ),
        None)


def adjustTanCpSpread_OneEndOnly(nc_In, bT1WorkEnd, fRadius_AcceptableMin=None, fDevTol=None, nc_forDevComp=None, bDebug=False):
    """
    Returns:
        Success: rg.NurbsCurve(New), float(minRadius), float(deviation)
        Fail: None, str(Feedback)
    
    This is not the latest function.
    """

    idxCp_Pos_A = (nc_In.Points.Count - 1) if bT1WorkEnd else 0
    idxCp_Tan_A = (nc_In.Points.Count - 2) if bT1WorkEnd else 1

    fRadius_Min_In = getMinimumRadius(nc_In)

    if fRadius_Min_In is None:
        return None, "Minimum radius could not be obtained."
    if bDebug: sEval='fRadius_Min_In'; print sEval+': ',eval(sEval)
    
    if fRadius_AcceptableMin and fRadius_Min_In >= fRadius_AcceptableMin:
        return None, "Curve's minimum radius is already >= {}.".format(
            fRadius_AcceptableMin)

    ncs_Winner_perScaleDir = [None, None]
    fRadii_MaxMin_perScaleDir = [None, None]
    fDev_perScaleDir = [None, None]
    dev = None
    
    vTan = nc_In.Points[idxCp_Tan_A].Location - nc_In.Points[idxCp_Pos_A].Location
    vTan.Unitize()

    radius_epsilon = (1e-3) * sc.doc.ModelAbsoluteTolerance

    # Try longer and shorter tangent control point spreads.
    for iDir, scaleDelta in enumerate((-0.01, 0.01)):
        #scale = 1.0 + scaleDelta

        i = 1

        while True:
            sc.escape_test()

            nc_WIP = nc_In.Duplicate()

            #if bDebug: sEval='scale'; print sEval+': ',eval(sEval),

            #xform = rg.Transform.Scale(
            #        anchor=nc_WIP.Points[idxCp_Pos_A].Location,
            #        scaleFactor=scale)
            #pt_Target = nc_WIP.Points[idxCp_Tan_A].Location
            #pt_Target.Transform(xform)

            vTrans = float(i) * sc.doc.ModelAbsoluteTolerance * vTan
            if iDir == 0:
                vTrans = rg.Vector3d.Negate(vTrans)

            xform = rg.Transform.Translation(vTrans)
            pt_Target = nc_In.Points[idxCp_Tan_A].Location
            pt_Target.Transform(xform)

            nc_WIP.Points[idxCp_Tan_A] = pt_Target
            if bDebug: stopwatch.Restart()
            fRadius_Min_WIP = getMinimumRadius(nc_WIP)
            if bDebug:
                stopwatch.Stop()
                timeElapsed = stopwatch.Elapsed.TotalSeconds
                s  = "{:.2f} seconds for ".format(timeElapsed)
                s += "getMinimumRadius"
                print s
                sEval='fRadius_Min_WIP'; print '  '+sEval+': ',eval(sEval),

            if (fRadius_Min_WIP - radius_epsilon) <= fRadius_Min_In:
                #sc.doc.Objects.AddCurve(nc_WIP)#;1/0
                if bDebug:
                    sEval='fRadius_Min_WIP <= fRadius_Min_In'; print sEval+': ',eval(sEval)
                    print "  Minimum radius is not increasing," \
                          " so break."
                break

            if fDevTol is not None and nc_forDevComp is not None:
                if bDebug: stopwatch.Restart()
                dev = getMaximumDeviation(nc_forDevComp, nc_WIP)
                if bDebug:
                    stopwatch.Stop()
                    timeElapsed = stopwatch.Elapsed.TotalSeconds
                    s  = "{:.2f} seconds for ".format(timeElapsed)
                    s += "getMaximumDeviation"
                    print s
                    sEval='dev'; print sEval+': ',eval(sEval),
                    sEval='dev > fDevTol'; print sEval+': ',eval(sEval)
                if dev > fDevTol:
                    if bDebug:
                        print "  Deviation is out of tolerance, so break."
                    break

            if fRadius_AcceptableMin:
                if (fRadius_Min_WIP + radius_epsilon) >= fRadius_AcceptableMin:
                    print "Curve's minimum radius, {} is >= {}.".format(
                        fRadius_Min_WIP, fRadius_AcceptableMin)
                    return (nc_WIP, fRadius_Min_WIP, fDevTo), None

            if bDebug:
                print "Curve minimum radius is larger than that of the starting curve."

            if ncs_Winner_perScaleDir[iDir] is None:
                ncs_Winner_perScaleDir[iDir] = nc_WIP
                fRadii_MaxMin_perScaleDir[iDir] = fRadius_Min_WIP
                fDev_perScaleDir[iDir] = dev

                #scale += scaleDelta
                i += 1
                continue

            if (fRadius_Min_WIP - radius_epsilon) > fRadii_MaxMin_perScaleDir[iDir]:
                if bDebug:
                    sEval='fRadius_Min_WIP > fRadii_MaxMin_perScaleDir[iDir]'; print sEval+': ',eval(sEval)
                    print "Minimum radius is increasing."
                ncs_Winner_perScaleDir[iDir].Dispose()
                
                ncs_Winner_perScaleDir[iDir] = nc_WIP
                fRadii_MaxMin_perScaleDir[iDir] = fRadius_Min_WIP
                fDev_perScaleDir[iDir] = dev

                #scale += scaleDelta
                i += 1
                continue


            if bDebug:
                print "Curve minimum radius is no longer increasing, so break"

            nc_WIP.Dispose()

            break

    if not (ncs_Winner_perScaleDir[0] or ncs_Winner_perScaleDir[1]):
        return None, "No winner."

    if ncs_Winner_perScaleDir[0] and ncs_Winner_perScaleDir[1]:
        iDir_Winner = fRadii_MaxMin_perScaleDir.index(max(fRadii_MaxMin_perScaleDir))
    else:
        iDir_Winner = 0 if ncs_Winner_perScaleDir[0] else 1

    return (
        (
            ncs_Winner_perScaleDir[iDir_Winner],
            fRadii_MaxMin_perScaleDir[iDir_Winner],
            fDev_perScaleDir[iDir_Winner]
        ),
        None)


def adjustTanCpSpread_BothEndsSimultaneously(nc_In, fRadius_AcceptableMin=None, fDevTol=None, nc_forDevComp=None, bDebug=False):
    """
    Returns:
        Success: rg.NurbsCurve(New), float(minRadius), float(deviation)
        Fail: None, str(Feedback)
    
    210908: New.
    """

    idx_G0_T0 = 0
    idx_G1_T0 = 1
    idx_G0_T1 = nc_In.Points.Count - 1
    idx_G1_T1 = nc_In.Points.Count - 2


    fRadius_Min_In = getMinimumRadius(nc_In)

    if fRadius_Min_In is None:
        return None, "Minimum radius could not be obtained."
    if bDebug: sEval='fRadius_Min_In'; print sEval+': ',eval(sEval)
    
    if fRadius_AcceptableMin and fRadius_Min_In >= fRadius_AcceptableMin:
        return None, "Curve's minimum radius is already >= {}.".format(
            fRadius_AcceptableMin)

    ncs_Winner_perScaleDir = []
    fRadii_MaxMin_perScaleDir = []
    fDev_perScaleDir = []
    dev = None

    vTan_T0 = nc_In.Points[idx_G1_T0].Location - nc_In.Points[idx_G0_T0].Location
    vTan_T0.Unitize()

    vTan_T1 = nc_In.Points[idx_G0_T1].Location - nc_In.Points[idx_G1_T1].Location
    vTan_T1.Unitize()

    radius_epsilon = (1e-3) * sc.doc.ModelAbsoluteTolerance

    # Try longer and shorter tangent control point spreads.
    for iDir_T0, iDir_T1 in zip((0,0,1,1),(0,1,0,1)):
        sc.escape_test()

        ncs_Winner_thisScaleDir = None
        fRadii_MaxMin_thisScaleDir = None
        fDev_thisScaleDir = None

        i = 1

        while True:
            sc.escape_test()

            if bDebug: sEval='i'; print sEval+': ',eval(sEval)

            nc_WIP = nc_In.Duplicate()

            vTrans_T0 = float(i) * sc.doc.ModelAbsoluteTolerance * vTan_T0
            vTrans_T1 = float(i) * sc.doc.ModelAbsoluteTolerance * vTan_T1

            if iDir_T0 == 0:
                vTrans_T0 = rg.Vector3d.Negate(vTrans_T0)
            if iDir_T1 == 0:
                vTrans_T1 = rg.Vector3d.Negate(vTrans_T1)

            xform_G1_T0 = rg.Transform.Translation(vTrans_T0)
            pt_Target = nc_In.Points[idx_G1_T0].Location
            pt_Target.Transform(xform_G1_T0)

            nc_WIP.Points[idx_G1_T0] = pt_Target

            xform_G1_T1 = rg.Transform.Translation(vTrans_T1)
            pt_Target = nc_In.Points[idx_G1_T1].Location
            pt_Target.Transform(xform_G1_T1)

            nc_WIP.Points[idx_G1_T1] = pt_Target


            fRadius_Min_WIP = getMinimumRadius(nc_WIP)
            if bDebug:
                sEval='fRadius_Min_WIP'; print '  '+sEval+': ',eval(sEval),

            if (fRadius_Min_WIP - radius_epsilon) <= fRadius_Min_In:
                #sc.doc.Objects.AddCurve(nc_WIP)#;1/0
                if bDebug:
                    sEval='fRadius_Min_WIP <= fRadius_Min_In'; print sEval+': ',eval(sEval)
                    print "  Minimum radius is not increasing," \
                          " so break."
                break

            if fDevTol is not None and nc_forDevComp is not None:
                dev = getMaximumDeviation(nc_forDevComp, nc_WIP)
                if bDebug:
                    sEval='dev'; print sEval+': ',eval(sEval),
                    sEval='dev > fDevTol'; print sEval+': ',eval(sEval)
                if dev > fDevTol:
                    if bDebug:
                        print "  Deviation is out of tolerance, so break."
                    break

            if fRadius_AcceptableMin:
                if (fRadius_Min_WIP + radius_epsilon) >= fRadius_AcceptableMin:
                    print "Curve's minimum radius, {} is >= {}.".format(
                        fRadius_Min_WIP, fRadius_AcceptableMin)
                    return (nc_WIP, fRadius_Min_WIP, fDevTo), None

            if bDebug:
                print "Curve minimum radius is larger than that of the starting curve."

            if ncs_Winner_thisScaleDir is None:
                ncs_Winner_thisScaleDir = nc_WIP
                fRadii_MaxMin_thisScaleDir = fRadius_Min_WIP
                fDev_thisScaleDir = dev

                i += 1
                continue

            if (fRadius_Min_WIP - radius_epsilon) > fRadii_MaxMin_thisScaleDir:
                if bDebug:
                    sEval='fRadius_Min_WIP > fRadii_MaxMin_thisScaleDir'; print sEval+': ',eval(sEval)
                    print "Minimum radius is increasing."
                ncs_Winner_thisScaleDir.Dispose()
                
                ncs_Winner_thisScaleDir = nc_WIP
                fRadii_MaxMin_thisScaleDir = fRadius_Min_WIP
                fDev_thisScaleDir = dev

                i += 1
                continue


            if bDebug:
                print "Curve minimum radius is no longer increasing, so break"

            nc_WIP.Dispose()

            break


    if ncs_Winner_thisScaleDir is not None:
        ncs_Winner_perScaleDir.append(ncs_Winner_thisScaleDir)
        fRadii_MaxMin_perScaleDir.append(fRadii_MaxMin_thisScaleDir)
        fDev_perScaleDir.append(fDev_thisScaleDir)


    if len(ncs_Winner_perScaleDir) == 0:
        return None, "No winner."


    if len(ncs_Winner_perScaleDir) == 1:
        return (
            (
                ncs_Winner_perScaleDir[0],
                fRadii_MaxMin_perScaleDir[0],
                fDev_perScaleDir[0]
            ),
            None)


    idx_Winner = fRadii_MaxMin_perScaleDir.index(max(fRadii_MaxMin_perScaleDir))

    return (
        (
            ncs_Winner_perScaleDir[idx_Winner],
            fRadii_MaxMin_perScaleDir[idx_Winner],
            fDev_perScaleDir[idx_Winner]
        ),
        None)


def areAbsEqual(a, b, decimal_tolerance=None):
    if decimal_tolerance is None: decimal_tolerance = 1e-6
    return abs(a-b) <= decimal_tolerance


def areRelEqual(a, b, decimal_tolerance=None):
    if decimal_tolerance is None: decimal_tolerance = 1e-6
    return abs(a-b)/max(a,b) <= decimal_tolerance


def adjustTanCpSpread_IsolatedEnd(nc_In, bT1WorkEnd, fRadius_AcceptableMin=None, fDevTol=None, nc_forDevComp=None, epsilon=1e-8, bDebug=False):
    """
    Returns:
        Success: tuple(rg.NurbsCurve(New), float(minRadius), float(deviation)), None
        Fail: None, str(Feedback)
    
    Recreated using epsilon and binary search and trimming curve to study only the working end.
    """

    if Rhino.RhinoApp.ExeVersion < 7:
        raise Exception("Unsupported routine for Rhino < 7.")

    fRads_Min_In_Whole = getMinimaRadii_V7(nc_In, bDebug=bDebug)

    if fRads_Min_In_Whole is None:
        return None, "Minimum radii could not be obtained."
    if bDebug: sEval='fRads_Min_In_Whole'; print sEval+': ',eval(sEval)
    
    if fRadius_AcceptableMin and fRads_Min_In_Whole >= fRadius_AcceptableMin:
        return None, "Curve's minimum radius is already >= {}.".format(
            fRadius_AcceptableMin)

    if nc_In.SpanCount == 1:
        print "Single-spanned curve."
        t_Trim_T0 = t_Trim_T1 = None
    else:
        if bT1WorkEnd:
            m = nc_In.Knots.KnotMultiplicity(nc_In.Knots.Count - nc_In.Degree - 1)
            if m > 1:
                # Last segment.
                t_Trim_T0 = nc_In.Knots[nc_In.Knots.Count - nc_In.Degree - 1]
                t_Trim_T1 = nc_In.Domain.T1
            else:
                # m == 1
                if nc_In.SpanCount == 2:
                    t_Trim_T0 = t_Trim_T1 = None
                else:
                    # Last 2 segments.
                    t_Trim_T0 = nc_In.Knots[nc_In.Knots.Count - nc_In.Degree - 2]
                    t_Trim_T1 = nc_In.Domain.T1
        else:
            m = nc_In.Knots.KnotMultiplicity(nc_In.Degree)
            if m > 1:
                # First segment.
                t_Trim_T0 = nc_In.Domain.T0
                t_Trim_T1 = nc_In.Knots[nc_In.Degree]
            else:
                # m == 1
                if nc_In.SpanCount == 2:
                    t_Trim_T0 = t_Trim_T1 = None
                else:
                    # First 2 segments.
                    t_Trim_T0 = nc_In.Domain.T0
                    t_Trim_T1 = nc_In.Knots[nc_In.Degree + 1]

    if bT1WorkEnd:
        idxCp_G0 = nc_In.Points.Count - 1
        idxCp_G1 = nc_In.Points.Count - 2
        idxCp_G2 = nc_In.Points.Count - 3
        vPosToTan = -nc_In.TangentAtEnd
    else:
        idxCp_G0 = 0
        idxCp_G1 = 1
        idxCp_G2 = 2
        vPosToTan = nc_In.TangentAtStart


    ncs_Winner_perScaleDir = [None, None]
    fRadii_MaxMin_perScaleDir = [None, None]
    fDev_perScaleDir = [None, None]
    dev = None

    epsilon_pt_loc = 0.1 * sc.doc.ModelAbsoluteTolerance

    # New code.  TODO: Describe this.
    def createEvalCurve(nc):
        if t_Trim_T0 is None:
            return nc.DuplicateCurve()
        return nc.Trim(t0=t_Trim_T0, t1=t_Trim_T1)


    nc_In_Eval = createEvalCurve(nc_In)
    #sc.doc.Objects.AddCurve(nc_In_Eval); #sc.doc.Views.Redraw(); 1/0

    fRads_Mins_In_Eval = getMinimaRadii_V7(nc_In_Eval, bDebug=bDebug)

    if nc_In.SpanCount == 1 and len(fRads_Mins_In_Eval) == 2:
        #print abs(fRads_Mins_In_Eval[0] - fRads_Mins_In_Eval[1])
        if areAbsEqual(fRads_Mins_In_Eval[0], fRads_Mins_In_Eval[1]):
            return None, "Single-spanned curve already contains 2 equal minima radius." \
                "  TODO: Simultaneously modify 2nd control points on both ends."


    def largerRadiusDirection(nc):
        """
        Returns: -1 for shorter, +1 for longer, 0 for no change.
        """

        nc_0_Eval = createEvalCurve(nc)

        fRads_0_Eval = getMinimaRadii_V7(nc_0_Eval, bDebug=bDebug)

        nc_Longer = nc.DuplicateCurve()
        nc_Longer.Points[idxCp_G1] = rg.Point3d(
            nc.Points[idxCp_G1].Location + epsilon_pt_loc*vPosToTan)
        #sc.doc.Objects.AddCurve(nc_Longer)

        nc_Longer_Eval = createEvalCurve(nc_Longer)
        #sc.doc.Objects.AddCurve(nc_Longer_Eval); #sc.doc.Views.Redraw(); 1/0

        fRads_Longer_Eval = getMinimaRadii_V7(nc_Longer_Eval, bDebug=bDebug)

        nc_Shorter = nc.DuplicateCurve()
        nc_Shorter.Points[idxCp_G1] = rg.Point3d(
            nc.Points[idxCp_G1].Location - epsilon_pt_loc*vPosToTan)
        #sc.doc.Objects.AddCurve(nc_Shorter); sc.doc.Views.Redraw(); 1/0

        nc_Shorter_Eval = createEvalCurve(nc_Shorter)
        #sc.doc.Objects.AddCurve(nc_Shorter_Eval); sc.doc.Views.Redraw(); 1/0

        fRads_Shorter_Eval = getMinimaRadii_V7(nc_Shorter_Eval, bDebug=bDebug)

        #sEval='fRads_Shorter_Eval'; print sEval+': ',eval(sEval)
        #sEval='fRads_0_Eval'; print sEval+': ',eval(sEval)
        #sEval='fRads_Longer_Eval'; print sEval+': ',eval(sEval)

        fRad_Min_0_Eval = min(fRads_0_Eval)
        fRad_Min_Shorter_Eval = min(fRads_Shorter_Eval)
        fRad_Min_Longer_Eval = min(fRads_Longer_Eval)

        if (
            areAbsEqual(fRad_Min_0_Eval, fRad_Min_Longer_Eval) and
            areAbsEqual(fRad_Min_0_Eval, fRad_Min_Shorter_Eval)
        ):
            return 0
        elif fRad_Min_Longer_Eval > fRad_Min_0_Eval and fRad_Min_Shorter_Eval > fRad_Min_0_Eval:
            raise Exception("How can both modifications result in larger minimum radius?")
        else:
            if fRad_Min_Longer_Eval > fRad_Min_0_Eval:
                return 1
            elif fRad_Min_Shorter_Eval > fRad_Min_0_Eval:
                return -1
            else:
                return 0
                sEval='fRad_Min_Shorter_Eval'; print sEval+': ',eval(sEval)
                sEval='fRad_Min_0_Eval'; print sEval+': ',eval(sEval)
                sEval='fRad_Min_Longer_Eval'; print sEval+': ',eval(sEval)
                raise Exception("What happened?")

    iDir = largerRadiusDirection(nc_In)


    if iDir == 0:
        return None, "Minimum radius is already at its maximum."

    dist_G0_to_G1 = nc_In.Points[idxCp_G0].Location.DistanceTo(nc_In.Points[idxCp_G1].Location)
    dist_G0_to_G2 = nc_In.Points[idxCp_G0].Location.DistanceTo(nc_In.Points[idxCp_G2].Location)

    if iDir == 1:
        dist_L = dist_G0_to_G1
        # TODO: Is there something better than 0.99?
        dist_H = 0.99 * dist_G0_to_G2
    elif iDir == -1:
        # TODO: Is there something better than 0.01?
        dist_L = 0.01 * dist_G0_to_G2
        dist_H = dist_G0_to_G1



    while not areAbsEqual(dist_L, dist_H):
        if sc.escape_test(throw_exception=False):
            sc.doc.Objects.AddCurve(nc_WIP_M); sc.doc.Views.Redraw()
            raise Exception("Escape key pressed.  Current WIP curve was added.")

        nc_WIP_L = nc_In.DuplicateCurve()
        nc_WIP_L.Points[idxCp_G1] = rg.Point3d(nc_In.Points[idxCp_G0].Location + dist_L*vPosToTan)
        fRad_Min_L = getMinimumRadius(createEvalCurve(nc_WIP_L))

        nc_WIP_H = nc_In.DuplicateCurve()
        nc_WIP_H.Points[idxCp_G1] = rg.Point3d(nc_In.Points[idxCp_G0].Location + dist_H*vPosToTan)
        fRad_Min_H = getMinimumRadius(createEvalCurve(nc_WIP_H))

        dist_M = 0.5*(dist_L + dist_H)

        nc_WIP_M = nc_In.DuplicateCurve()
        nc_WIP_M.Points[idxCp_G1] = rg.Point3d(nc_In.Points[idxCp_G0].Location + dist_M*vPosToTan)
        fRad_Min_M = getMinimumRadius(createEvalCurve(nc_WIP_M))

        if bDebug:
            print '-'*20
            sEval='dist_L'; print sEval+': ',eval(sEval)
            sEval='  fRad_Min_L'; print sEval+': ',eval(sEval)
            sEval='dist_M'; print sEval+': ',eval(sEval)
            sEval='  fRad_Min_M'; print sEval+': ',eval(sEval)
            sEval='dist_H'; print sEval+': ',eval(sEval)
            sEval='  fRad_Min_H'; print sEval+': ',eval(sEval)

        if areAbsEqual(fRad_Min_L, fRad_Min_H):
            if bDebug:
                sEval='dist_L'; print sEval+': ',eval(sEval)
                sEval='dist_H'; print sEval+': ',eval(sEval)

            return (nc_WIP_M, fRad_Min_M, None), None


        if fRad_Min_M > fRad_Min_H > fRad_Min_L:
            if bDebug: print "M > H > L"
            dist_L = dist_M

        elif fRad_Min_M > fRad_Min_L > fRad_Min_H:
            if bDebug: print "M > L > H"
            dist_L = dist_M

        elif fRad_Min_M > fRad_Min_L and fRad_Min_M > fRad_Min_H:
            if bDebug: print "Mid has largest radius."

            iDir = largerRadiusDirection(nc_WIP_M)

            if iDir == 0:
                raise Exception("No change in radius.")
            if iDir == 1:
                dist_L = dist_M
            elif iDir == -1:
                dist_H = dist_M

        elif fRad_Min_H > fRad_Min_M > fRad_Min_L:
            if bDebug: print "H > M > L"
            dist_L = dist_M

        elif fRad_Min_L > fRad_Min_M > fRad_Min_H:
            if bDebug: print "L > M > H"
            dist_H = dist_M


    return (nc_WIP_M, fRad_Min_M, None), None



    if not (ncs_Winner_perScaleDir[0] or ncs_Winner_perScaleDir[1]):
        return None, "No winner."

    if ncs_Winner_perScaleDir[0] and ncs_Winner_perScaleDir[1]:
        iDir_Winner = fRadii_MaxMin_perScaleDir.index(max(fRadii_MaxMin_perScaleDir))
    else:
        iDir_Winner = 0 if ncs_Winner_perScaleDir[0] else 1

    return (
        (
            ncs_Winner_perScaleDir[iDir_Winner],
            fRadii_MaxMin_perScaleDir[iDir_Winner],
            fDev_perScaleDir[iDir_Winner]
        ),
        None)


def createNcsWithAdjTanCpSpreadForMaxMinRad_ANY_WIP(nc0, bT1WorkEnd, fDevTol=None, nc_forDevComp=None, bDebug=False):
    """
    returns list of new NurbsCurves with ANY minimum radii larger than the original.
    """

    idxCp_Pos_A = (nc0.Points.Count - 1) if bT1WorkEnd else 0
    idxCp_Tan_A = (nc0.Points.Count - 2) if bT1WorkEnd else 1

    fRadii_Min_0 = getMinimumRadii(nc0)
    if fRadii_Min_0 is None: return
    if bDebug: sEval='fRadii_Min_0'; print sEval+': ',eval(sEval)

    ncs_WithLargerMinRadii = []
    
    # Try longer and shorter tangent control point spreads.
    for scaleIncr in 0.1, -0.05:
        nc_A_Pre = nc0.Duplicate()
        scale = 1.0 + scaleIncr

        fRadii_Minima_ThisScaleDir = None

        while True:
            sc.escape_test()
            nc_WIP = nc0.Duplicate()
            if bDebug: sEval='scale'; print sEval+': ',eval(sEval),
            xform = rg.Transform.Scale(
                    anchor=nc_WIP.Points[idxCp_Pos_A].Location,
                    scaleFactor=scale)
            pt_Target = nc_WIP.Points[idxCp_Tan_A].Location
            pt_Target.Transform(xform)
            nc_WIP.Points[idxCp_Tan_A] = pt_Target
            if bDebug: stopwatch.Restart()
            fRadius_Min_WIP = getMinimumRadius(nc_WIP)
            fRadii_Min_WIP = getMinimumRadii(nc_WIP)
            if bDebug: sEval='fRadii_Min_WIP'; print sEval+': ',eval(sEval),

            def hasAnyMinRadiusDecreased(fRadii_WIP, fRadii_CompareWith):
                for iR in range(len(fRadii_WIP)):
                    if fRadii_WIP[iR] <= fRadii_CompareWith[iR]:
                        return True
                    if bDebug: sEval='fRadii_WIP[iR] <= fRadii_CompareWith[iR]'; print sEval+': ',eval(sEval)
                return False

            if hasAnyMinRadiusDecreased(fRadii_Min_WIP, fRadii_Min_0):
                # Minimum radius is not increasing,
                # so stop and capture the previous curve.
                break

            if fRadii_Minima_ThisScaleDir is not None:
                if hasAnyMinRadiusDecreased(fRadii_Min_WIP, fRadii_MaxMin_ThisScaleDir):
                    # Minimum radius is now decreasing from best found in this scale direction,
                    # so stop and capture the previous curve.
                    break

            if fDevTol is not None and nc_forDevComp is not None:
                dev = getMaximumDeviation(nc_forDevComp, nc_WIP)
                if bDebug:
                    sEval='dev'; print sEval+': ',eval(sEval),
                    sEval='dev > fDevTol'; print sEval+': ',eval(sEval),
                if dev > fDevTol:
                    # Deviation is out of tolerance,
                    # so capture the previous curve.
                    break
            elif bDebug: print

            nc_A_Pre.Dispose()
            nc_A_Pre = nc_WIP
            scale += scaleIncr
            fRadii_MaxMin_ThisScaleDir = fRadii_Min_WIP

        nc_WIP.Dispose()

        if scale != 1.0 + scaleIncr:
            ncs_WithLargerMinRadii.append(nc_A_Pre)
        else:
            nc_A_Pre.Dispose()
    
    return ncs_WithLargerMinRadii


def createNurbsCurve(nc_In, fRadius_AcceptableMin=None, fDevTol=None, bAdjWeightInDeg2Bezier=True, bDebug=False):
    """
    Returns:
        Success: tuple(rgNurbsCurve, fMinRadius, fDeviation), None
        Fail: None, sLog
    """
    
    if not isinstance(nc_In, rg.NurbsCurve):
        return None, "Not a NurbsCurve."



    def adjustCpLocs_WithWhileLoop():
        nc_Result = None
        radius_Winner = None
        dev_Winner = None
    
        nc_WIP = nc_In.DuplicateCurve()
    
        while True:
            sc.escape_test()

            if bDebug: stopwatch.Restart()
            
            fRadius_Min_StartOfWhile = getMinimumRadius(nc_WIP)


            if bDebug:
                stopwatch.Stop()
                timeElapsed = stopwatch.Elapsed.TotalSeconds
                s  = "{:.2f} seconds for ".format(timeElapsed)
                s += "getMinimumRadius"
                print s
            
            if fRadius_AcceptableMin and fRadius_Min_StartOfWhile >= fRadius_AcceptableMin:
                return None, "Curve's minimum radius is already >= {}.".format(
                    fRadius_AcceptableMin)
    
            for bT1WorkEnd in False, True:
                rc = adjustTanCpSpread_OneEndOnly(
                    nc_WIP,
                    bT1WorkEnd,
                    fRadius_AcceptableMin,
                    fDevTol,
                    nc_In,
                    bDebug=bDebug,
                    )

                if rc[0] is not None:
                    if bDebug: print "Some improvement was made."
                    nc_Result, radius_Winner, dev_Winner = rc[0]
                    nc_WIP.Dispose()
                    nc_WIP = nc_Result

            if bDebug: stopwatch.Restart()
            
            rc = adjustTanCpSpread_BothEndsSimultaneously(
                nc_WIP,
                fRadius_AcceptableMin,
                fDevTol,
                nc_In,
                bDebug=bDebug,
                )

            if rc[0] is not None:
                if bDebug: print "Some improvement was made."
                nc_Result, radius_Winner, dev_Winner = rc[0]
                nc_WIP.Dispose()
                nc_WIP = nc_Result

            if not nc_Result:
                nc_WIP.Dispose()
                return None, "Curve already has maximized minimum radius."
            
            fRadius_Min_EndOfWhile = getMinimumRadius(nc_WIP)
            
            
            if bDebug:
                stopwatch.Stop()
                timeElapsed = stopwatch.Elapsed.TotalSeconds
                s  = "{:.2f} seconds for ".format(timeElapsed)
                s += "getMinimumRadius"
                print s
    
            if bDebug: print fRadius_Min_StartOfWhile, fRadius_Min_EndOfWhile
    
            if abs(fRadius_Min_StartOfWhile - fRadius_Min_EndOfWhile) <= 0.1*sc.doc.ModelAbsoluteTolerance:
                break # out of while loop.

        return (nc_Result, radius_Winner, dev_Winner), None


    def adjustCpLocs_NoWhileLoop():
        nc_Result = None
        radius_Winner = None
        dev_Winner = None
    
        nc_WIP = nc_In.DuplicateCurve()

        if bDebug: stopwatch.Restart()

        fRadius_Min_StartOfWhile = getMinimumRadius(nc_WIP)

        #print getMinimaRadii_V7(nc_WIP); 1/0


        if bDebug:
            stopwatch.Stop()
            timeElapsed = stopwatch.Elapsed.TotalSeconds
            s  = "{:.2f} seconds for ".format(timeElapsed)
            s += "getMinimumRadius"
            print s
            
        if fRadius_AcceptableMin and fRadius_Min_StartOfWhile >= fRadius_AcceptableMin:
            return None, "Curve's minimum radius is already >= {}.".format(
                fRadius_AcceptableMin)
    
        for bT1WorkEnd in False, True:
            #rc = adjustTanCpSpread(
            #    nc_WIP,
            #    bT1WorkEnd,
            #    fRadius_AcceptableMin,
            #    fDevTol,
            #    nc0,
            #    bDebug=bDebug,
            #    )

            rc = adjustTanCpSpread_IsolatedEnd(
                nc_WIP,
                bT1WorkEnd,
                fRadius_AcceptableMin,
                fDevTol,
                nc_In,
                bDebug=bDebug,
                )

            if rc[0] is not None:
                if bDebug: print "Some improvement was made."
                nc_Result, radius_Winner, dev_Winner = rc[0]
                #if bT1WorkEnd:
                #    sc.doc.Objects.AddCurve(nc_Result); sc.doc.Views.Redraw(); 1/0
                nc_WIP.Dispose()
                nc_WIP = nc_Result
    
        if bDebug: stopwatch.Restart()
            
        if not nc_Result:
            nc_WIP.Dispose()
            return None, "Curve already has minimum radius."
            
        fRadius_Min_EndOfWhile = getMinimumRadius(nc_WIP)
            
            
        if bDebug:
            stopwatch.Stop()
            timeElapsed = stopwatch.Elapsed.TotalSeconds
            s  = "{:.2f} seconds for ".format(timeElapsed)
            s += "getMinimumRadius"
            print s
    
        if bDebug: print fRadius_Min_StartOfWhile, fRadius_Min_EndOfWhile
    
        return (nc_Result, radius_Winner, dev_Winner), None


    nc_In = nc_In.ToNurbsCurve()

    if nc_In.Degree == 2 and nc_In.SpanCount == 1:
        if bAdjWeightInDeg2Bezier:
            rc = adjustWeight(
                nc_In,
                fRadius_AcceptableMin=fRadius_AcceptableMin,
                fDevTol=fDevTol,
                nc_forDevComp=nc_In,
                bDebug=bDebug)
        else:
            nc_In.Dispose()
            return None, "Curve doesn't have enough control points to be modified."

    else:
        fRads_Mins_In = getMinimaRadii_V7(nc_In, bDebug=bDebug)

        if nc_In.SpanCount == 1: # and len(fRads_Mins_In) == 2:
            #print abs(fRads_Mins_In_Eval[0] - fRads_Mins_In_Eval[1])
            rc = adjustCpLocs_WithWhileLoop()

            #if areAbsEqual(fRads_Mins_In[0], fRads_Mins_In[1]):
                #return None, "Single-spanned curve" \
                #    " already contains 2 equal minima radii." \
                #    "  TODO: Simultaneously modify 2nd control points on both ends."


        elif nc_In.SpanCount == 2 and nc_In.Knots.KnotMultiplicity(nc_In.Degree) == 1:
            rc = adjustCpLocs_WithWhileLoop()
            #if areAbsEqual(fRads_Mins_In[0], fRads_Mins_In[1]):
                #return None, "Double-spanned curve with simple knot between" \
                #    " already contains 2 equal minima radii." \
                #    "  TODO: Simultaneously modify 2nd control points on both ends."

        else:
            rc = adjustCpLocs_NoWhileLoop()

    if not rc[0]:
        nc_In.Dispose()
        return None, rc[1]

    (nc_Result, radius_Winner, dev_Winner), sLog = rc

    if nc_Result.EpsilonEquals(nc_In, 1e-8):
        nc_In.Dispose()
        nc_Result.Dispose()
        return None, "New curve is within 1e-8 of old and will not be modified."

    nc_In.Dispose()

    return (nc_Result, radius_Winner, dev_Winner), sLog


def processCurveObjects(objrefs, fRadius_AcceptableMin=None, fDevTol=None, bAdjWeightInDeg2Bezier=True, bReplace=True, bEcho=True, bDebug=False):

    gs1 = []
    sLogs = []

    for objref in objrefs:

        gCrv0_A = objref.ObjectId
        nc_In = objref.Curve()

        rc = createNurbsCurve(
                nc_In=nc_In,
                fRadius_AcceptableMin=fRadius_AcceptableMin,
                fDevTol=fDevTol,
                bAdjWeightInDeg2Bezier=bAdjWeightInDeg2Bezier,
                bDebug=bDebug,
        )
        if rc[0] is None:
            sLogs.append(rc[1])
            continue


        def processResults(gCrv0_ToMod, rgCrv0_ToMod, winningResults):

            rgNurbsCrv1, fRadius_Min_New, fDev = winningResults

            g1 = None

            if bReplace:
                if sc.doc.Objects.Replace(gCrv0_ToMod, rgNurbsCrv1):
                    if bEcho: print "Curve was replaced."
                    g1 = gCrv0_ToMod
                else:
                    g1 = None
            else:
                g1 = sc.doc.Objects.AddCurve(rgNurbsCrv1)
                if g1 == g1.Empty:
                    g1 = None 
                else:
                    if bEcho: print "Curve was added."

            sc.doc.Views.Redraw()

            s = ""

            rc = spb_Crv_inflections.getInflectionParameters(rgCrv0_ToMod)
            i_Infl_Ct_Pre = len(rc) if rc else None
            rc = spb_Crv_inflections.getInflectionParameters(rgNurbsCrv1)
            i_Infl_Ct_Post = len(rc) if rc else None
            s = "InflectionCount:{}->{}".format(i_Infl_Ct_Pre, i_Infl_Ct_Post)

            rc = getMinimumRadius(rgCrv0_ToMod)
            fRadius_Min_Original = rc if rc is not None else None
            sRadius_Min_Original = '{:.{}f}'.format(fRadius_Min_Original, sc.doc.ModelDistanceDisplayPrecision)
            sRadius_Min_New = '{:.{}f}'.format(fRadius_Min_New, sc.doc.ModelDistanceDisplayPrecision)
            if s: s += "  "
            s += "  MinimumRadius:{}->{}".format(sRadius_Min_Original, sRadius_Min_New)

            if fDev is None:
                fDev = getMaximumDeviation(rgCrv0_ToMod, rgNurbsCrv1)
            s += "  Deviation:{}".format(formatDistance(fDev))
            print s
            
            return g1


        winningResults = rc[0]

        g1 = processResults(gCrv0_A, nc_In, winningResults)
        if g1: gs1.append(g1)
        winningResults[0].Dispose()


        nc_In.Dispose()

    return gs1, sLogs


def main():
    
    rc = getInput()
    if rc is None: return
    (
        objrefs,
        bSetMinTargetRad,
        fRadius_AcceptableMin,
        bLimitCrvDev,
        fDevTol,
        bAdjWeightInDeg2Bezier,
        bReplace,
        bEcho,
        bDebug,
    ) = rc

    if not bDebug:
        sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gCs_Modified, sLogs = processCurveObjects(
        objrefs,
        fRadius_AcceptableMin if bSetMinTargetRad else None,
        fDevTol=fDevTol if bLimitCrvDev else None,
        bAdjWeightInDeg2Bezier=bAdjWeightInDeg2Bezier,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
    )

    for sLog in set(sLogs):
        print "[{}] {}".format(sLogs.count(sLog), sLog)

    for gC in gCs_Modified:
        sc.doc.Objects.Select(gC)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()