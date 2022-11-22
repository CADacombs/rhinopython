"""
190703-04: Created. 
190708: Refactored.
190908: Added bAddDot and iDotHeight.
...
200121-28: Created a more accurate routine that finds parameters.
        Added bIncludeArcs and function that finds and returns mid points of arc-shaped segments.
200307: Changed a tolerance value from 1e-12 to 1e-9 for determining curvature in a linear section.
200319: Stopped using ZeroTolerance in many places since that value is different in V7 than it is in V5 & V6.
200420: Further improved parameters/points found by refining tolerances used.
200803: Fixed variable name typo.
210315: Now processes ArcCurves in PolyCurves.

The intention of the script is to provide minimum radius (maximum curvature) data
within a practical (design for manufacturability) tolerance,
not to provide the most accurate results.

TODO:
    Test a routine with 1e-12 in V7.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
from Rhino import RhinoMath as rm # This style allows IntelliSense, etc. in VS editor.
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System.Collections.Generic import List
from System.Drawing import Color


class Opts:

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


    key = 'fMinRadToReport'; keys.append(key)
    values[key] = 0.0
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fMaxRadToReport'; keys.append(key)
    values[key] = 0.0
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fRadEqTol'; keys.append(key)
    values[key] = 10.0*sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bIncludeArcs'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddPt'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHeight'; keys.append(key)
    values[key] = 11
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddOnlyMinimumPt'; keys.append(key)
    values[key] = True
    names[key] = "At"
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='AllMinima', onValue='Minimum')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddOnlyMinimumPtOfAllCrvs'; keys.append(key)
    values[key] = True
    names[key] = "Of"
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='EachCrv', onValue='AllCrvs')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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
    Get edges and options values.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = rd.ObjectType.Curve # Curve is also used for brep edges.

    go.AcceptNumber(True, True)

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    bPreselectedObjsChecked = False

    print "0 for {} or {} will disable that limit.".format(
        Opts.names['fMinRadToReport'], Opts.names['fMaxRadToReport'])

    idxs_Opts = {}

    while True:
        key = 'fMinRadToReport'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'fMaxRadToReport'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'fRadEqTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bIncludeArcs'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bAddPt'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bAddDot'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'iDotHeight'; idxs_Opts[key] = (Opts.riAddOpts[key](go) if Opts.values['bAddDot']
                                              else None)
        key = 'bAddOnlyMinimumPt'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bAddOnlyMinimumPtOfAllCrvs'; idxs_Opts[key] = (Opts.riAddOpts[key](go)
                                                              if Opts.values['bAddOnlyMinimumPt']
                                                              else None)
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return
        else:
            # An option was selected or a number was entered.
            key = 'fMinRadToReport'
            if res == ri.GetResult.Number:
                Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 0.0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
            key = 'fMaxRadToReport'
            if Opts.riOpts[key].CurrentValue < 0.0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
            key = 'fRadEqTol'
            if Opts.riOpts[key].CurrentValue < 2.0**(-53):
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()


def tryGetArcDataFromWholeCrv(crv, tolerance=None):

    if tolerance is None:
        tolerance = 0.001 * sc.doc.ModelAbsoluteTolerance

    if isinstance(crv, rg.ArcCurve):
        arc = crv.Arc
        ts = crv.DivideByCount(segmentCount=2, includeEnds=False)
        return crv.PointAt(ts[0]), arc.Radius

    bSuccess, arc = crv.TryGetArc(tolerance=tolerance)
    if bSuccess:
        #sc.doc.Objects.AddArc(arc); sc.doc.Views.Redraw(); 1/0
        ts = crv.DivideByCount(segmentCount=2, includeEnds=False)
        return crv.PointAt(ts[0]), arc.Radius
    else:
        pass
        #sc.doc.Objects.AddCurve(crv)
        #sc.doc.Views.Redraw()


def getMinimumRadiiData_ArcSegments(rgCrv0, bDebug=False):
    """
    Returns:
        On success: list(Point3ds of midpoint of arc intervals), list(floats of arc radii)
        On fail: None
    """


    if not isinstance(rgCrv0, rg.Curve):
        raise ValueError("{} not an accepted input for getArcData.".format(rgCrv0.GetType().Name))

    sType = rgCrv0.GetType().Name

    if sType == 'BrepEdge':
        crvTemp = rgCrv0.EdgeCurve
        sType = crvTemp.GetType().Name
    else:
        crvTemp = rgCrv0.DuplicateCurve()

    if sType in ('LineCurve', 'PolylineCurve'):
        crvTemp.Dispose()
        return

    if crvTemp.IsLinear(1.0 / 2.0**32):
        crvTemp.Dispose()
        return


    def tryGetArcDataFromNurbsSpans(nc):
        pOut = []
        rOut = []
        for iSpan in range(nc.SpanCount):
            seg = nc.Trim(domain=nc.SpanDomain(iSpan))
            rc = tryGetArcDataFromWholeCrv(seg)
            if rc:
                pOut.append(rc[0])
                rOut.append(rc[1])
            seg.Dispose()
        return pOut, rOut



    rc = tryGetArcDataFromWholeCrv(crvTemp)
    if rc:
        crvTemp.Dispose()
        return [rc[0]], [rc[1]]


    pts_Out = []
    radii_Out = []

    if sType == 'NurbsCurve':
        rc = tryGetArcDataFromNurbsSpans(crvTemp)
        if rc[0]:
            pts_Out.extend(rc[0])
            radii_Out.extend(rc[1])
    else:
        #'PolyCurve'
        crvTemp.RemoveNesting()
        for seg in crvTemp.DuplicateSegments():
            rc = tryGetArcDataFromWholeCrv(seg)
            if rc:
                pts_Out.append(rc[0])
                radii_Out.append(rc[1])
            elif isinstance(seg, rg.NurbsCurve):
                rc = tryGetArcDataFromNurbsSpans(seg)
                if rc[0]:
                    pts_Out.extend(rc[0])
                    radii_Out.extend(rc[1])
            seg.Dispose()

    crvTemp.Dispose()

    return pts_Out, radii_Out


def getMinimumRadiiData_NurbsCrv(rgCrv0, bDebug=False):
    """
    Analyzes 6 divisions between each Greville parameter.

    Returns:
        On success: list(float(ts_Maxima)), list(float(radii))
        On fail: None, str(Feedback)
    """


    def isAcceptableNurbsCurve(rgCrv0):
        sType = rgCrv0.GetType().Name
        if sType in ("ArcCurve", "LineCurve", "PolylineCurve"): return False
        if rgCrv0.IsArc(1.0 / 2.0**32): return False
        if rgCrv0.IsLinear(1.0 / 2.0**32): return False
        return True
    if not isAcceptableNurbsCurve(rgCrv0): return


    def f2s(f):
        return "{:.17e}".format(f)


    cross = rg.Vector3d.CrossProduct


    def curvatureAt(t, bRightSide):
        side = (
            rg.CurveEvaluationSide.Above if bRightSide
            else rg.CurveEvaluationSide.Below)
        ds = nc.DerivativeAt(t, derivativeCount=2, side=side)
        return cross(cross(ds[1], ds[2]), ds[1]) / ds[1].Length**4.0


    # rm.ZeroTolerance ==
    #   V5 & V6: 1e-12
    #   V7: 2.3283064365386962890625e-10 (1.0 / 2.0**32)
    eps_ON = 1.0 / 2.0**32 # OpenNURBS ZeroTolerance.
    eps_param = 1e-9

    # Machine epsilon per C, C++, and Python ( https://en.wikipedia.org/wiki/Machine_epsilon )
    eps_Mach = 2.0**(-52) #2.2204460492503130808472633361816e-16

    # eps_toUse was determined through trial and error.
    # Based on results of future executions of this script,
    # it may be determined that this value will need to be increased.
    eps_toUse = 2.0**(-42) # 2.27373675443232059478759765625e-13


    def areEpsilonEqual(a, b, epsilon):
        if abs(a) <= epsilon and abs(b) <= epsilon: return True
        delta = abs(a-b)
        #if 0.55147058823529 < a < 0.551470588235295:
        #    print f2s(delta)
        if delta <= epsilon: return True
        # A relative comparison is used for other values.
        fRelComp = delta/max(abs(a), abs(b))
        #print f2s(fRelComp)
        return fRelComp <= epsilon


    def findParameter(tL_Start, tR_Start):

        tL = tL_Start
        tM = None # For debugging.
        tR = tR_Start
        
        tM = 0.5*tL + 0.5*tR

        j_forDebug = 18

        j = 0

        while True:
            if sc.escape_test(False):
                print "Break at iT of {} and j of {}.".format(iT, j)
                sEval = 'tL_Start'; print sEval+':',f2s(eval(sEval))
                sEval = 'tL'; print sEval+':',f2s(eval(sEval))
                sEval = 'tR'; print sEval+':',f2s(eval(sEval))
                sEval = 'tR_Start'; print sEval+':',f2s(eval(sEval))
                return

            if tR > 0.5294:
                pass

            if iT == iT_forDebug and j == j_forDebug:
                pass

            Kl = curvatureAt(tL, True)
            if not Kl.IsValid:
                print "Curvature vector of left test point is invalid," \
                    " so no results will be returned."
                return

            Km = curvatureAt(tM, False) # False or True shouldn't make a difference.
            if not Km.IsValid:
                print "Curvature vector of left test point is invalid," \
                    " so no results will be returned."
                return

            Kr = curvatureAt(tR, False)
            if not Kr.IsValid:
                print "Curvature vector of right test point is invalid," \
                    " so no results will be returned."
                return

            kL = Kl.Length
            kM = Km.Length
            kR = Kr.Length

            rad_kL = 1.0/(kL) if kL > eps_param else None
            rad_kM = 1.0/(kM) if kM > eps_param else None
            rad_kR = 1.0/(kR) if kR > eps_param else None

            if rad_kL is None and rad_kR is None:
                raise ValueError("Both rad_kL and rad_kR are None!")

            if rad_kL is None or rad_kR is None:
                #sEval = 'tL_Start'; print sEval+':',eval(sEval)
                #sc.doc.Objects.AddPoint(nc.PointAt(tL_Start))
                #sEval = 'tR_Start'; print sEval+':',eval(sEval)
                #sc.doc.Objects.AddPoint(nc.PointAt(tR_Start))
                if rad_kL is None:
                    #sEval = 'tL'; print sEval+':',eval(sEval)
                    #sc.doc.Objects.AddPoint(nc.PointAt(tL))
                    tL = 0.75*tL + 0.25*tR
                elif rad_kR is None:
                    #sEval = 'tR'; print sEval+':',eval(sEval)
                    #sc.doc.Objects.AddPoint(nc.PointAt(tR))
                    tR = 0.25*tL + 0.75*tR
                #sc.doc.Views.Redraw()
                continue

            #if iT == iT_forDebug and j >= j_forDebug:
                #print '='*20
                #tX = 0.24137988639995100 # For debug.
                #print "{} {} {} {} {}".format(f2s(tL), f2s(rad_kL)                          , f2s(kL), f2s(nc.CurvatureAt(tL).Length), f2s(nc.CurvatureAt(tL)*nc.CurvatureAt(tL)))
                #print "{} {} {} {} {}".format(f2s(tX), f2s(1.0/curvatureAt(tX, False).Length), f2s(0.0), f2s(nc.CurvatureAt(tX).Length), f2s(nc.CurvatureAt(tX)*nc.CurvatureAt(tX)))
                #print "{} {} {} {} {}".format(f2s(tR), f2s(rad_kR)                          , f2s(kR), f2s(nc.CurvatureAt(tR).Length), f2s(nc.CurvatureAt(tR)*nc.CurvatureAt(tR)))
                #if tL > tX or tR < tX:
                #    pass

            # Example found where tL is showing to be larger than tR is where
            # difference is about 3.79696274421804e-14 larger than kR.


            if areEpsilonEqual(tR, tL, eps_toUse):
                #if tL > 0.5294:
                #    print f2s(tR)
                #    print f2s(tL)
                return tL, tR, kL, kR


            # Not accurate enough for uniform (or with only internal monoknot?) curves.
            #if abs(rad_kL - rad_kR) < (0.1 * sc.doc.ModelAbsoluteTolerance):
            #    return tL, tR, kL, kR

            #if areEpsilonEqual(rad_kL, rad_kR, eps_toUse):
            #    return tL, tR, kL, kR


            if areEpsilonEqual(Kl*Kl, Kr*Kr, eps_toUse):
                return tL, tR, kL, kR


            # Example found where kL is showing to be larger than kR is where
            # kL is about 2.63618080391892e-15 larger than kR.

            if areEpsilonEqual(kL, kR, eps_toUse):
                return tL, tR, kL, kR

            if kL > kM > kR:
                tR = tM
            elif kR > kM > kL:
                tL = tM
            elif kM > kL > kR:
                tR = tM
            elif kM > kR > kL:
                tL = tM
            elif kL > kR and kM == kR:
                tR = tM
            elif kR > kL and kM == kL:
                tL = tM
            else:
                # Possibly approaching a curvature minimum / radius maximum.
                return tL, tR, kL, kR
                print '='*20
                print "Fine results of curvature, etc., cannot be trusted."
                print "tL:{} kL:{}".format(f2s(tL), f2s(kL))
                print "{} {}".format(f2s(tM), f2s(kM))
                print "tR:{} kR{}".format(f2s(tR), f2s(kR))
                print "delta t:{}".format(f2s(tR-tL))
                print "delta k:{}".format(f2s(abs(kR-kL)))
                print "delta rad:{}".format(f2s(abs(rad_kR-rad_kL)))
                if kR > kL:
                    return tR, tR, kR, kR
                else:
                    return tL, tL, kL, kL

            tM = 0.5*tL + 0.5*tR

            # For debug.
            if tL > 0.241379886:
                pass

            j += 1


    nc = rgCrv0.ToNurbsCurve()
    domain_notNormal = nc.Domain

    # Work with parameters in domain [0.0, 1.0]
    # to avoid certain floating point accuracy errors.
    nc.Domain = rg.Interval(0.0,1.0)


    def getParametersToCheck(nc):
        if nc.IsPeriodic:
            ts_Grev = [t for t in nc.GrevilleParameters() if 0.0 <= t <= 1.0]
        else:
            ts_Grev = nc.GrevilleParameters()
    
    
        # Determine fractional divisions between Grevilles.
        total_divs_btwGrevs = 6
        precision_of_divs = 0.125 # Use an exact floating point value, e.g. 0.125, not 0.1.
        divs_btwGrevs = []
        for i in range(1, total_divs_btwGrevs):
            divs_btwGrevs.append(
                round(
                ((1.0 / precision_of_divs) * float(i) / float(total_divs_btwGrevs)), 0) *
                precision_of_divs)
    
    
        ts_toUse = []
        for iT in range(len(ts_Grev)-1):
            ts_toUse.append(ts_Grev[iT])
            for div in divs_btwGrevs:
                ts_toAdd = (1.0-div)*ts_Grev[iT] + div*ts_Grev[iT+1]
                #sc.doc.Objects.AddPoint(rg.NurbsCurve.PointAt(nc, t=ts_toAdd))
                ts_toUse.append(ts_toAdd)
        ts_toUse.append(ts_Grev[-1])
    
    
        # Add non-full polyknots parameters since radius can be different on each side.
        ts_NonFullPolyKnot = []
        iKnot = 0
        t = 0.0
        while t < 1.0:
            t = nc.Knots[iKnot]
            if t < 0.0: continue
            multi = nc.Knots.KnotMultiplicity(iKnot)
            if 1 < multi < nc.Degree:
                ts_NonFullPolyKnot.append(t)
            iKnot += multi
        ts_toUse += ts_NonFullPolyKnot
        ts_toUse = sorted(ts_toUse)
    
        return ts_toUse


    ts_toUse = getParametersToCheck(nc)

    #for t in ts_toUse:
        #sc.doc.Objects.AddPoint(rgCrv0.PointAt(domain_notNormal.ParameterAt(t)))
    #sc.doc.Views.Redraw(); return


    ts_Pass_Norm = []
    radii_Out = []
    
    bLastFoundWasOnRight = False
    bFoundAtStartOfPeriodic = False
    
    # Furthest left parameter is not continuous at T0
    # and at a non-continuous break.
    bOkToAddOnLeft = True

    iT_forDebug = 4

    for iT in range(len(ts_toUse)-1):
        
        # For bDebug:
        if iT == iT_forDebug:
            pass
        
        #print "{}:".format(iT),
        tL_Start = ts_toUse[iT]
        tL = tL_Start
        tR = None # This is here for its placement in debugging variable list.
        tR_Start = ts_toUse[iT+1]
        tR = tR_Start

        if tL_Start > 0.2413:
            pass

        Kl = curvatureAt(tL, True)
        if not Kl.IsValid:
            print "Curvature vector of left test point is invalid," \
                " so no results will be returned."
            return
        
        Kr = curvatureAt(tR, False)
        if not Kr.IsValid:
            print "Curvature vector of right test point is invalid," \
                " so no results will be returned."
            return
        
        kL_Start = None
        kL = Kl.Length
        kR = Kr.Length

        kL_Start = kL
        kR_Start = kR
        kL_Start_NextCrv = None
        
        # Test for linear interval.
        if abs(kL_Start) < eps_param or abs(kR_Start) < eps_param:
            bOkToAddOnLeft = True
            continue

        fRadL_Start = 1.0/kL_Start
        fRadR_Start = 1.0/kR_Start

        # Testing for circular arc interval.
        if areEpsilonEqual(fRadL_Start, fRadR_Start, 0.001*sc.doc.ModelAbsoluteTolerance):
            bOkToAddOnLeft = True
            continue

        kL_Start_NextCrv = curvatureAt(tR_Start, True).Length
        
        #sEval = 'kL_Start_NextCrv'; print sEval+':',eval(sEval)
        
        if kL_Start_NextCrv < eps_toUse:
            # Start of linear interval.
            # Greatest curvature is already kL_Start.
            pass
        else:
            # For bDebug:
            if iT == iT_forDebug:
                pass

            fRadL_Start_NextCrv = 1.0 / kL_Start_NextCrv # This variable is used for viewing in the debugger.

            rc = findParameter(tL_Start, tR_Start)
            if not rc: return

            tL, tR, kL, kR = rc


        #print abs(tL-tL_Start)
        #print abs(tL-tL_Start)/max(abs(tL),abs(tL_Start))


        #if tL_Start > 0.5294:
        #    print f2s(tL)
        #    print f2s(tL_Start)
        #    pass

        # For bDebug:
        if iT == iT_forDebug:
            pass

        if areEpsilonEqual(tL, tL_Start, eps_ON):
            #print f2s(abs(tL - tL_Start))
            if bOkToAddOnLeft:
                if bDebug: print "{}: Left start".format(iT)
                ts_Pass_Norm.append(tL_Start)
                rad = 1.0/kL_Start
                radii_Out.append(rad)
            bLastFoundWasOnRight = False
            
            if iT == 0 and nc.IsPeriodic:
                bFoundAtStartOfPeriodic = True
            elif bFoundAtStartOfPeriodic and iT == len(ts_toUse)-2:
                ts_Pass_Norm = ts_Pass_Norm[1:]
                radii_Out = radii_Out[1:]
                print "T0 removed from passing parameters" \
                    " since a smaller radius was found" \
                    " in the last division of the curve."
        elif areEpsilonEqual(tR, tR_Start, eps_ON):
            #print f2s(abs(tR - tR_Start))
            if iT == len(ts_toUse)-2:
                if nc.IsPeriodic:
                    print "Passing parameter already found" \
                        " at beginning of periodic domain" \
                        " and will not be added to other parameters."
                elif nc.IsClosed and areEpsilonEqual(1.0/kR_Start, radii_Out[0], eps_toUse):
                    print "Passing parameter already found" \
                        " at beginning of closed domain" \
                        " at the same radius minima" \
                        " and will not be added to other parameters."
                else:
                    if bDebug: print "{}: Right start at end of entire curve.".format(iT)
                    ts_Pass_Norm.append(tR_Start)
                    rad = 1.0/kR_Start
                    radii_Out.append(rad)
            else:
                if not areEpsilonEqual(kR_Start, kL_Start_NextCrv, eps_toUse):
                    if bDebug:
                        print "{}: Right start at kink.".format(iT)
                        #print f2s(kR_Start)
                        #print f2s(kL_Start_NextCrv)
                        #print '*'*4
                        #print f2s(abs(kR_Start - kL_Start_NextCrv))
                        #print '*'*4
                    ts_Pass_Norm.append(tR_Start)
                    rad = 1.0/kR_Start
                    radii_Out.append(rad)
            bLastFoundWasOnRight = True
        else:
            if bDebug: print "{}: Not at a starting left or right parameter.".format(iT)
            if abs(tL-tR) <= eps_toUse:
                ts_Pass_Norm.append(tR)
                rad = 1.0/kR
                radii_Out.append(rad)
            elif eps_toUse < abs(tL-tR) < eps_ON:
                print "eps_toUse < abs(tL-tR) < eps_ON"
                ts_Pass_Norm.append(0.5*tL + 0.5*tR)
                rad = 1.0/(0.5*kL + 0.5*kR)
                radii_Out.append(rad)
            else:
                length = nc.GetLength(subdomain=rg.Interval(tL, tR))
                if bDebug:
                    s  = "Normalized parameter space between L and R: {}".format(
                        f2s(abs(tL-tR)))
                    s += "  Length: {}".format(f2s(length))
                    print s
                if length >= sc.doc.ModelAbsoluteTolerance:
                    pass
                    #sc.doc.Objects.AddPoint(nc.PointAt(tL))
                    #sc.doc.Objects.AddPoint(nc.PointAt(tR))
                    #sc.doc.Views.Redraw()
                    
                    #print "Curve length between points having same curvature is {}.".format(length)
                ts_Pass_Norm.append(0.5*tL + 0.5*tR)
                rad = 1.0/(0.5*kL + 0.5*kR)
                radii_Out.append(rad)
            bLastFoundWasOnRight = False

            if bFoundAtStartOfPeriodic and iT == len(ts_toUse)-2:
                ts_Pass_Norm = ts_Pass_Norm[1:]
                radii_Out = radii_Out[1:]
                print "T0 removed from passing parameters" \
                    " since a smaller radius was found" \
                    " in the last division of the curve."

        bOkToAddOnLeft = (
            bLastFoundWasOnRight or
            not areEpsilonEqual(kR_Start, kL_Start_NextCrv, eps_toUse))

    ts_Out = [domain_notNormal.ParameterAt(t) for t in ts_Pass_Norm]

    nc.Dispose()

    return ts_Out, radii_Out


def getMinimumRadiiPoints(curve, bIncludeArcs=True, bDebug=False):

    rc = tryGetArcDataFromWholeCrv(curve)
    if rc:
        return [rc[0]], [rc[1]]

    rcNurbs= getMinimumRadiiData_NurbsCrv(curve, bDebug)

    rc_forArcs = getMinimumRadiiData_ArcSegments(curve, bDebug) if bIncludeArcs else None

    if not ((rcNurbs and rcNurbs[0]) or (rc_forArcs and rc_forArcs[0])): return


    pts_Out = []
    radii_Out = []

    if rcNurbs and rcNurbs[0]:
        ts, radii_Out = rcNurbs
        for t in ts:
            pts_Out.append(curve.PointAt(t))

    if rc_forArcs and rc_forArcs[0]:
        pts_Out.extend(rc_forArcs[0])
        radii_Out.extend(rc_forArcs[1])

    return pts_Out, radii_Out


def getMinimumRadii(rgCrv0, bIncludeArcs=True, bDebug=False):
    """
    """


    sType_rgCrv0 = rgCrv0.GetType().Name

    if sType_rgCrv0 == "ArcCurve":
        return [rgCrv0.Radius] if bIncludeArcs else None
    elif sType_rgCrv0 == "LineCurve":
        return [Rhino.RhinoMath.UnsetValue]
    elif sType_rgCrv0 == "PolylineCurve":
        return [Rhino.RhinoMath.UnsetValue]
    
    rc = getMinimumRadiiData_NurbsCrv(rgCrv0, bDebug)
    if not rc or not rc[0]: return

    ts, radii = rc

    return radii


def getMinimumRadius(rgCrv0, bIncludeArcs=True, bDebug=False):
    """
    """

    rc = getMinimumRadii(rgCrv0, bDebug)

    if rc is None: return

    return min(rc)


def getFormattedDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def processCurves(curvesAndEdges0, **kwargs):
    """
    """
    
    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fMinRadToReport = getOpt('fMinRadToReport')
    fMaxRadToReport = getOpt('fMaxRadToReport')
    fRadEqTol = getOpt('fRadEqTol')
    bIncludeArcs = getOpt('bIncludeArcs')
    bAddPt = getOpt('bAddPt')
    bAddDot = getOpt('bAddDot')
    iDotHeight = getOpt('iDotHeight')
    bAddOnlyMinimumPt = getOpt('bAddOnlyMinimumPt')
    bAddOnlyMinimumPtOfAllCrvs = getOpt('bAddOnlyMinimumPtOfAllCrvs')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    gCrvs0 = []
    for curveOrEdge0 in curvesAndEdges0:
        gCrv0 = rs.coerceguid(curveOrEdge0)
        if gCrv0:
            gCrvs0.append(gCrv0)
    
    rgCrvs0_Found_All = []
    pts_Pass_perCrv = []
    radii_Pass_perCrv = []
    
    sFails = []
    
    radii_Accepted_AllCrvs = []
    
    len_curvesAndEdges0 = len(curvesAndEdges0)
    
    idxs_AtTenths = [int(round(0.1*i*len_curvesAndEdges0,0)) for i in range(10)]
    
    for iC, curveOrEdge0 in enumerate(curvesAndEdges0):
        if iC in idxs_AtTenths:
            Rhino.RhinoApp.SetCommandPrompt(
                    "Processing curve {} ...".format(
                    "" if len_curvesAndEdges0 == 1 else "{} of {} ".format(
                        iC+1, len_curvesAndEdges0)))
        
        rgCrv0 = rs.coercecurve(curveOrEdge0) # Will return various rg.Curves, including rg.BrepEdge.
        if rgCrv0 is None: return None, "Geometry for {} not found!".format(rgCrv0)
        
        if isinstance(rgCrv0, rg.BrepEdge):
            bDeleteInput = False
        else:
            gCrv0 = rs.coerceguid(curveOrEdge0)
            rdCrv0 = rs.coercerhinoobject(curveOrEdge0)
        
        if isinstance(curveOrEdge0, rd.ObjRef):
            rd.ObjRef.GeometryComponentIndex
        
        sType_rgCrv0 = rgCrv0.GetType().Name
        
        
        
        rc = getMinimumRadiiPoints(rgCrv0, bIncludeArcs, bDebug)
        if rc is None:
            sFails.append("getMinMaxRadiusPoints returned None.")
            rgCrv0.Dispose()
            continue

        pts_Minima, radii_Maxima = rc

        if not pts_Minima:
            sFails.append("getMinMaxRadiusPoints returned no points.")
            rgCrv0.Dispose()
            continue
        

        if fMinRadToReport or fMaxRadToReport:
            pts_Pass_ThisCrv = []
            radii_Pass_ThisCrv = []
            for iT in range(len(pts_Minima)):
                pt = pts_Minima[iT]
                radius = radii_Maxima[iT]
                if fMinRadToReport and radius >= (fMinRadToReport - fRadEqTol):
                    pts_Pass_ThisCrv.append(pt)
                    radii_Pass_ThisCrv.append(radius)
                elif fMaxRadToReport and radius <= (fMaxRadToReport + fRadEqTol):
                    pts_Pass_ThisCrv.append(pt)
                    radii_Pass_ThisCrv.append(radius)
            if radii_Pass_ThisCrv:
                radii_Accepted_AllCrvs.extend(radii_Pass_ThisCrv)
                rgCrvs0_Found_All.append(rgCrv0)
                pts_Pass_perCrv.append(pts_Pass_ThisCrv)
                radii_Pass_perCrv.append(radii_Pass_ThisCrv)
            else:
                sFail  = "No radii found"
                if fMinRadToReport:
                    sFail += " above {}".format(fMinRadToReport)
                if fMaxRadToReport:
                    if fMinRadToReport:
                        sFail += ", and"
                    sFail += " below {}".format(fMaxRadToReport)
                sFail += "."
                sFails.append(sFail)
                rgCrv0.Dispose()
                continue
        else:
            rgCrvs0_Found_All.append(rgCrv0)
            pts_Pass_perCrv.append(pts_Minima)
            radii_Pass_perCrv.append(radii_Maxima)
    

    if not pts_Pass_perCrv:
        s  = "No minimum radii found"
        if fMinRadToReport:
            s += " above {}".format(fMinRadToReport)
            if fMaxRadToReport: s += ","
        if fMaxRadToReport:
            s += " below {}".format(fMaxRadToReport)
        s += "."
        print s
        return

    if bEcho:
        if len_curvesAndEdges0 == 1:
            s =  "{} radius minima found.".format(len(radii_Maxima))
            if fMinRadToReport or fMaxRadToReport:
                s += "  {} are".format(
                        sum([len(rads) for rads in radii_Pass_perCrv]))
                if fMinRadToReport:
                    s += " above {}".format(fMinRadToReport)
                    if fMaxRadToReport: s += ","
                if fMaxRadToReport:
                    s += " below {}".format(fMaxRadToReport)
                s += "."
            s += "  Minimum is {}.".format(
                    getFormattedDistance(min(radii_Maxima)))
            if sFails: s += '\n' + sFails[0]
            print s
        elif len_curvesAndEdges0 > 1:
            s = "Out of {} curves selected:".format(len_curvesAndEdges0)
            for sFail in set(sFails):
                s += "\n[{}] {}".format(sFails.count(sFail), sFail)
            if fMinRadToReport:
                s += "  {0} radius minima [{2:.{4}f},{3:.{4}f}] are less than {1}.".format(
                        len(radii_Accepted_AllCrvs), fMinRadToReport,
                        min(radii_Accepted_AllCrvs), max(radii_Accepted_AllCrvs),
                        sc.doc.ModelDistanceDisplayPrecision+1)
            else:
                s += "\n{} radius minima found.".format(
                        sum([len(_) for _ in radii_Pass_perCrv]))
            print s

    if not (bAddDot or bAddPt): return

    if bAddDot:
        attrib_Red = rd.ObjectAttributes()
        attrib_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
        attrib_Red.ObjectColor = Color.Red

    if radii_Pass_perCrv:
        pts_toAdd = []
        dots_toAdd = []
        
        if bAddOnlyMinimumPt:
            if bAddOnlyMinimumPtOfAllCrvs:
                radius_Mininum_AllCrvs = min([__ for _ in radii_Pass_perCrv for __ in _])
                sMinRad = "R{0:.{1}f}".format(
                        radius_Mininum_AllCrvs,
                        sc.doc.ModelDistanceDisplayPrecision)
                
                # Include all minima within tolerance to minimum.
                for iC in range(len(rgCrvs0_Found_All)):
                    rgCrv0 = rgCrvs0_Found_All[iC]
                    pts_Minima = pts_Pass_perCrv[iC]
                    radii_Minima_PerCrv = radii_Pass_perCrv[iC]

                    for iR in range(len(radii_Minima_PerCrv)):
                        radius = radii_Minima_PerCrv[iR]
                        if abs(radius - radius_Mininum_AllCrvs) <= fRadEqTol:
                            pt = pts_Minima[iR]
                            if pt not in pts_toAdd:
                                pts_toAdd.append(pt)
                            if bAddDot:
                                rgDot = rg.TextDot(sMinRad, pt)
                                rgDot.FontHeight = iDotHeight
                                sc.doc.Objects.AddTextDot(rgDot, attrib_Red)
                s  = "{} points added".format(len(pts_toAdd))
                s += " {} added".format(
                    'points' if (bAddPt and not bAddDot) else 'dots')
                s += " at {}, minimum radius of all curves.".format(sMinRad)
            else:
                # Add a minimum of each curve.
                for iC in range(len(rgCrvs0_Found_All)):
                    rgCrv0 = rgCrvs0_Found_All[iC]
                    pts_Minima = pts_Pass_perCrv[iC]
                    radii_Minima_PerCrv = radii_Pass_perCrv[iC]

                    for iR in range(len(radii_Minima_PerCrv)):
                        radius = radii_Minima_PerCrv[iR]
                        if abs(radius - min(radii_Minima_PerCrv)) <= fRadEqTol:
                            pt = pts_Minima[iR]
                            if bAddPt:
                                if pt not in pts_toAdd:
                                    pts_toAdd.append(pt)
                            if bAddDot:
                                sMinRad = "R{0:.{1}f}".format(
                                        radius,
                                        sc.doc.ModelDistanceDisplayPrecision)
                                rgDot = rg.TextDot(sMinRad, pt)
                                rgDot.FontHeight = iDotHeight
                                dots_toAdd.append(rgDot)
                s  = "{} points added at each curve's minimum radius.".format(
                        len(pts_toAdd),
                        'points' if (bAddPt and not bAddDot) else 'dots')
        else:
            for iC in range(len(rgCrvs0_Found_All)):
                rgCrv0 = rgCrvs0_Found_All[iC]
                pts_Minima = pts_Pass_perCrv[iC]
                radii_Minima_PerCrv = radii_Pass_perCrv[iC]

                for iR in range(len(radii_Minima_PerCrv)):
                    radius = radii_Minima_PerCrv[iR]
                    pt = pts_Minima[iR]
                    if bAddPt:
                        # TODO: Test this in V7 since its ZeroTolerance is different.
                        if not pts_toAdd or not rg.Point3d.EpsilonEquals(pt, pts_toAdd[0], epsilon=1.0e-12):
                            pts_toAdd.append(pt)
                    if bAddDot:
                        sMinRad = "R{0:.{1}f}".format(
                                radius,
                                sc.doc.ModelDistanceDisplayPrecision)
                        rgDot = rg.TextDot(sMinRad, pt)
                        rgDot.FontHeight = iDotHeight
                        dots_toAdd.append(rgDot)

    if bAddDot:
        for dot in dots_toAdd:
            sc.doc.Objects.AddTextDot(dot, attrib_Red)
        if bAddOnlyMinimumPt:
            if bAddOnlyMinimumPtOfAllCrvs:
                s  = "{} dot".format(len(pts_toAdd))
                s += " added at minimum radius of all curves: "
                s += sMinRad
                s += "."
            else:
                s  = "{} total dots added at minimum radius of each curve.".format(
                        len(dots_toAdd))
        else:
            s  = "{} total dots added at all radius minima.".format(
                    len(dots_toAdd))

    if bAddPt:
        sc.doc.Objects.AddPoints(List[rg.Point3d](pts_toAdd)) # Python list has to be converted to .NET list for AddPoints.
        if bAddOnlyMinimumPt:
            if bAddOnlyMinimumPtOfAllCrvs:
                s  = "{} points".format(len(pts_toAdd))
                s += " added at {}, minimum radius of all curves.".format(sMinRad)
            else:
                s  = "{} total points added at minimum radius of each curve.".format(
                        len(pts_toAdd))
        else:
            s  = "{} total point(s) added at all radius minima.".format(
                    len(pts_toAdd))

    if bEcho: print s


def main():
    
    rc = getInput()
    if rc is None: return

    objrefs = rc[0]
    
    if Opts.values['bDebug']:
        pass
    else:
        sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll()
    
    rc = processCurves(curvesAndEdges0=objrefs)
    if not rc: return
    
    rgCrvs0_Found_All = rc
    
    if rgCrvs0_Found_All:
        [sc.doc.Objects.Select(objectId=__) for _ in rgCrvs0_Found_All for __ in _]
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()