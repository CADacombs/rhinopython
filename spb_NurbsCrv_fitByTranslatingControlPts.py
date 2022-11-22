"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190920, 23-24: Created.
191021: Changed some variable names.
220823, 0908-10, 1007: Further development.

Notes:
Greville points aren't directly translated because any curve tangency is not maintained,
and results are squiggly.

No longer applicable for this script: It doesn't seem that Greville points can be adjusted within 1e-12 accuracy, so use a minimum of 1e-11.

TODO:
    What to do about DistanceBetweenCurves bug where maximum distance is measured at
    perpendiculars which are not near shortest distances?


    Create alternative control point translation: Move Greville point on a temp curve,
    then adjust just the relative control point on the WIP curve.

    Upon one test, using a scale weight of 1.5 in createNcWithTransCp_ByGrevClosestPt
    was better than 1.0 or 2.0 or using createNcWithTransCp_CpFromGrevTrans.
    More testing is needed.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


debugHelp = {'gCrv' : None}


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fResolution'; keys.append(key)
    values[key] = 1e-6
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPreserveEndTans'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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

        if key == 'fResolution':
            if cls.riOpts[key].CurrentValue <= 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(objref_ToMod=None):
    """
    Get curve with optional input.
    """

    go = ri.Custom.GetObject()

    go.GeometryFilter = rd.ObjectType.Curve

    if objref_ToMod is None:
        go.SetCommandPrompt("Select curve to modify")
    else:
        rdObj_ToMod = objref_ToMod.Object()
        edge = objref_ToMod.Edge()
        if edge is None:
            wire = objref_ToMod.Curve()
        else:
            wire = None
        if Opts.values['bDebug']: print(edge, wire)

        go.SetCommandPrompt("Select reference curve")
        def customGeometryFilter(rdObj, rgObj, compIdx):
            if rdObj.Id != rdObj_ToMod.Id:
                return True
            if wire and not isinstance(rgObj, rg.BrepEdge):
                return False
            if edge and isinstance(rgObj, rg.BrepEdge):
                if edge.EdgeIndex == rgObj.EdgeIndex:
                    return False
            return True
        go.SetCustomGeometryFilter(customGeometryFilter)


    #go.SubObjectSelect = True

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('fResolution')
        addOption('bPreserveEndTans')
        addOption('bReplace')
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

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def formatDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def _createCrvWithTranslatedCp(nc_In, idx, vect):
    nc_Out = nc_In.Duplicate()
    nc_Out.Points.SetPoint(
        idx,
        nc_In.Points[idx].Location + vect,
        nc_In.Points[idx].Weight)
    return nc_Out


def _getMaxDev(rgCrvA, rgCrvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            rgCrvA,
            rgCrvB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
    if rc[0]:
        return rc[1]


def _func_calcMaxDev_wrapper(nc_WIP, idxCP, vTrans_Unit, crv_Ref):
    def _calcMaxDev_wrapped(dist):
        nc_New = _createCrvWithTranslatedCp(nc_WIP, idxCP, dist*vTrans_Unit)
        dev = _getMaxDev(nc_New, crv_Ref)
        nc_New.Dispose()
        return dev
    return _calcMaxDev_wrapped


def _translationVectorForG1CPLocationAlongTangent(nc_In, bT1WorkEnd, nc_forDevComp, fResolution, bDebug=False):
    """
    returns:
        Success: New NurbsCurve with deviation from nc_forDevComp less than nc_In
        Fail: None
    """


    def determineDir(dist_Start, fResolution):
        """
        Positive vector is from end CP[0 or -1] to next to end CP [1 or -2].
        
        Returns:
            dir: float(1.0 for along vector, -1.0 for against vector)
            dist: float(Sign is direction relative to vector.)
            dev: float(Returned to avoid having to recalculate this later.)
        """

        dev_Incr = calcMaxDev(dist_Start + fResolution)
        dev_Decr = calcMaxDev(dist_Start - fResolution)

        if dev_Start < dev_Incr and dev_Start < dev_Decr:
            sEval = "dev_Start"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "dev_Incr"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "dev_Decr"; print("{}: {}".format(sEval, eval(sEval)))
            return
            raise Exception("dev_Start < dev_Incr and dev_Start < dev_Decr")

        if dev_Incr == dev_Decr:
            print("Increase and decrease have the same result.  Is this the minimum?")
            return

        if dev_Incr < dev_Decr:
            #if bDebug: print("Spread should increase.")
            return 1.0, dist_Start+fResolution, dev_Incr
        else:
            #if bDebug: print("Spread should decrease.")
            return -1.0, dist_Start-fResolution, dev_Decr


    def getBinarySearchExtremes(fDir, dist_Min, dev_Min):

        #if bDebug:
        #    print('-'*30)
        #    print('getBinarySearchExtremes:')

        dist_Incr = fDir

        dist_Max = dist_Incr
        dev_Max = calcMaxDev(dist_Max)

        dist_BeyondMax = 2.0 * dist_Incr

        while True:
            sc.escape_test()

            dev_BeyondMax = calcMaxDev(dist_BeyondMax)

            if bDebug:
                sEval = "dist_Min"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_Max"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_BeyondMax"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Min"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Max"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_BeyondMax"; print("{}: {}".format(sEval, eval(sEval)))

            if dev_Max < dev_Min and dev_Max < dev_BeyondMax:
                return (
                    dist_Min,
                    dev_Min,
                    dist_BeyondMax,
                    dev_BeyondMax,
                    )

            if dev_BeyondMax > dev_Min and dev_BeyondMax > dev_Max:
                if dev_Max < dev_Min:
                    return (
                        dist_Min,
                        dev_Min,
                        dist_BeyondMax,
                        dev_BeyondMax,
                        )
                else:
                    return (
                        dist_Min,
                        dev_Min,
                        dist_Max,
                        dev_Max,
                        )

            if dev_BeyondMax < dev_Min:
                dist_Min = dist_Max
                dev_Min = dev_Max
                dist_Max = dist_BeyondMax
                dev_Max = dev_BeyondMax
            else:
                sEval = "dist_Min"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_Max"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_BeyondMax"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Min"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Max"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_BeyondMax"; print("{}: {}".format(sEval, eval(sEval)))

                raise Exception("Should dev_BeyondMax be as shown above?")

            dist_BeyondMax += dist_Incr


    def binarySearch(dist_Min, dev_Min, dist_Max, dev_Max):
        
        if bDebug:
            print('-'*30)
            print('binarySearch:')

        while True:
            sc.escape_test()

            dist_Mid = 0.5*dist_Min + 0.5*dist_Max

            #if abs(dist_Mid-dist_Min) <= 0.0001:#Rhino.RhinoMath.ZeroTolerance:
            #    return dist_Mid

            nc_WIP = _createCrvWithTranslatedCp(nc_In, idxCp_G1, dist_Mid*v_01)
            #sc.doc.Objects.AddCurve(nc_WIP)
            dev_Mid = _getMaxDev(nc_forDevComp, nc_WIP)
            #sEval = "dist_Mid"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "dev_Mid"; print("{}: {}".format(sEval, eval(sEval)))

            if bDebug:
                print('-'*20)
                sEval = "dist_Min"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_Mid"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_Max"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Min"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Mid"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev_Max"; print("{}: {}".format(sEval, eval(sEval)))


            if dev_Mid > dev_Min and dev_Mid > dev_Max:
                if bDebug:
                    print("Mid is greater than both extremes.")
                    sc.doc.Objects.AddCurve(nc_forDevComp)
                    sc.doc.Objects.AddCurve(nc_WIP)
                1/0
                return

            nc_WIP.Dispose()

            if (
                (abs(dist_Mid-dist_Min) > Rhino.RhinoMath.ZeroTolerance)
                and
                (abs(dev_Mid-dev_Min) <= Rhino.RhinoMath.ZeroTolerance)
            ):
                if bDebug:
                    print("Span of minimal deviations may mean that adjusting"
                    "tangency for deviation minimization may be moot at this time.")
                return

            if dev_Mid < dev_Min and dev_Mid < dev_Max:

                rc = determineDir(dist_Mid, fResolution)
                if rc is None:
                    return

                if rc[0] == -1.0:
                    if dist_Min < 0.0:
                        dist_Min = rc[1]
                        dev_Min = rc[2]
                    else:
                        dist_Max = rc[1]
                        dev_Max = rc[2]
                elif rc[0] == 1.0:
                    if dist_Min > 0.0:
                        dist_Min = rc[1]
                        dev_Min = rc[2]
                    else:
                        dist_Max = rc[1]
                        dev_Max = rc[2]
                else:
                    raise Exception("Exception!")

                #if bDebug: print(
                #    "Minimum deviation may be at either side of mid."
                #    "  dist_Mid will be returned.")
                #return dist_Mid

            elif dev_Mid < dev_Max < dev_Min:
                dist_Min = dist_Mid
                dev_Min = dev_Mid
            elif dev_Mid < dev_Min < dev_Max:
                dist_Max = dist_Mid
                dev_Max = dev_Mid
            elif dev_Mid < dev_Min:
                dist_Min = dist_Mid
                dev_Min = dev_Mid
            elif dev_Mid < dev_Max:
                dist_Max = dist_Mid
                dev_Max = dev_Mid
            else:
                1/0

            delta_Min_Max_Dist = abs(dist_Min-dist_Max)

            if delta_Min_Max_Dist <= fResolution:
                if bDebug:
                    sEval = "delta_Min_Max_Dist"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "dist_Mid"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "dev_Mid"; print("{}: {}".format(sEval, eval(sEval)))
                return dist_Mid


            delta_Min_Max_Dev = abs(dev_Min-dev_Max)

            if delta_Min_Max_Dist <= fResolution:
                if bDebug:
                    sEval = "delta_Min_Max_Dev"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "dist_Mid"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "dev_Mid"; print("{}: {}".format(sEval, eval(sEval)))
                return dist_Mid



    dev_Start = _getMaxDev(nc_forDevComp, nc_In)
    if bDebug: sEval='dev_Start'; print("{}: {}".format(sEval, eval(sEval)))

    idxCp_G0 = (nc_In.Points.Count - 1) if bT1WorkEnd else 0
    idxCp_G1 = (nc_In.Points.Count - 2) if bT1WorkEnd else 1

    v_01 = nc_In.Points[idxCp_G1].Location - nc_In.Points[idxCp_G0].Location
    v_01.Unitize()

    calcMaxDev = _func_calcMaxDev_wrapper(nc_In, idxCp_G1, v_01, nc_forDevComp)

    nc_Out = None

    rc = determineDir(0.0, fResolution)
    if rc is None:
        return

    fDir, dist_Min, dev_Min = rc

    if bDebug:
        print('-'*10)
        print('-'*10)
        print("From determineMinDist:")
        sEval = "dist_Min"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dev_Min"; print("{}: {}".format(sEval, eval(sEval)))
        print('-'*10)
        print('-'*10)

    if dist_Min == 0:
        return rg.Vector3d.Zero


    rc = getBinarySearchExtremes(fDir, dist_Min, dev_Min)
    (
        dist_Min,
        dev_Min,
        dist_Max,
        dev_Max,
        ) = rc

    if bDebug:
        print('-'*10)
        print('-'*10)
        print("From getBinarySearchExtremes:")
        sEval = "dist_Min"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dist_Max"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dev_Min"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dev_Max"; print("{}: {}".format(sEval, eval(sEval)))
        print('-'*10)
        print('-'*10)


    dist = binarySearch(dist_Min, dev_Min, dist_Max, dev_Max)
    if dist is None: return

    if bDebug:
        print('-'*10)
        print('-'*10)
        print("From binarySearch:")
        sEval = "dist"; print("{}: {}".format(sEval, eval(sEval)))
        print('-'*10)
        print('-'*10)

    return v_01 * dist


def knotMultiplicityPattern(nc):
    mp = []
    iK = 0
    while iK < nc.Knots.Count:
        m = nc.Knots.KnotMultiplicity(iK)
        mp.append(m)
        iK += m
    return mp


def fit_Adjust_G1_CP_spreads(nc_In, nc_forDevComp, fResolution, iCallCt, bDebug=False):
    """
    returns:
        Success: New NurbsCurve with deviation from nc_forDevComp less than nc_In
        Fail: None

    When adjusting curve ends that do not affect opposite ends, minimum knot multiplicity patterns are:
    31113, 323, 51115, 525, etc.
    """

    if nc_In.Points.Count == 3: return


    kmp = knotMultiplicityPattern(nc_In)
    if bDebug: sEval = "kmp"; print("{}: {}".format(sEval ,eval(sEval)))

    kmp_Inr = kmp[1:-1]
    if len(kmp_Inr) == 0:
        # No interior knots.
        bIndependentEndSpans = False
    elif all(m==1 for m in kmp_Inr):
        # All interior knots are simple.
        bIndependentEndSpans = len(kmp_Inr) >= 3
    else:
        # Some non-simple knots.
        bIndependentEndSpans = True

    if bDebug: sEval = "bIndependentEndSpans"; print("{}: {}".format(sEval ,eval(sEval)))

    # WIP: Looking for mv's sweet spots per iteration to this function.
    #mv = 1.0 - 0.5**iCallCt # Multiplier for vector.
    #if iCallCt < 3:
    #    mv = 0.5
    #else:
    #    mv = 0.5
    mv = 0.5

    if bDebug: sEval = "mv"; print("{}: {}".format(sEval ,eval(sEval)))

    if bDebug: print('-'*80,"\nCalculate vector near T0\n",'-'*40)

    v_G0G1_NearT0 = _translationVectorForG1CPLocationAlongTangent(
        nc_In,
        bT1WorkEnd=False,
        nc_forDevComp=nc_forDevComp.Duplicate().Trim(nc_forDevComp.Domain.T0, nc_forDevComp.Domain.Mid),
        fResolution=fResolution,
        bDebug=bDebug)

    if bDebug: print('-'*80,"\nCalculate vector near T1\n",'-'*40)

    v_G0G1_NearT1 = _translationVectorForG1CPLocationAlongTangent(
        nc_In,
        bT1WorkEnd=True,
        nc_forDevComp=nc_forDevComp.Duplicate().Trim(nc_forDevComp.Domain.Mid, nc_forDevComp.Domain.T1),
        fResolution=fResolution,
        bDebug=bDebug)

    #print(type(v_G0G1_NearT0))
    #return

    if v_G0G1_NearT0 is None and v_G0G1_NearT1 is None:
        return

    nc_Out = nc_In.Duplicate()

    if v_G0G1_NearT0 is not None:
        idxG1 = 1
        pt_New = nc_In.Points[idxG1].Location + mv*v_G0G1_NearT0
        nc_Out.Points.SetPoint(idxG1, pt_New, nc_In.Points[idxG1].Weight)

    if v_G0G1_NearT1 is not None:
        idxG1 = nc_In.Points.Count - 2
        pt_New = nc_In.Points[idxG1].Location + mv*v_G0G1_NearT1
        nc_Out.Points.SetPoint(idxG1, pt_New, nc_In.Points[idxG1].Weight)

    return nc_Out


def fitCurve(nc_toDeform, nc_forDevComp, **kwargs):
    """
    returns:
        Success: New NurbsCurve with deviation less than nc0
        Fail: None
    """

    if not nc_toDeform.IsValid: return

    if nc_toDeform.Points.Count < 3: return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fResolution = getOpt('fResolution')
    bPreserveEndTans = getOpt('bPreserveEndTans')
    bDebug = getOpt('bDebug')


    if bPreserveEndTans and nc_toDeform.Points.Count == 3: return


    fDev0 = _getMaxDev(nc_forDevComp, nc_toDeform)
    if bDebug: sEval='fDev0'; print(sEval+': ',eval(sEval))
    if fDev0 is None: return
    #if dev <= fDevTol:
    #    if bDebug: print("Deviation is already within tolerance.")
    #    # Deviation is out of tolerance,
    #    # so capture the previous curve.
    #    break

    nc_WIP = nc_toDeform.Duplicate() # Keep copy of original because nc_WIP will be modified per iteration.

    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt


    def calculateNewCPLocation_GrevilleToClosestPt(nc_toDeform, iCp, nc_Target, bDebug=False):
        pt_g0 = nc_toDeform.GrevillePoint(iCp)
        b, t = nc_Target.ClosestPoint(pt_g0)
        pt_Target = nc_Target.PointAt(t)
        dist = pt_Target.DistanceTo(pt_g0)
        if bDebug: sEval='dist'; print(sEval+': ',eval(sEval))
        if dist <= Rhino.RhinoMath.ZeroTolerance:
            return # No change.
        pt_cp_1 = (
            nc_toDeform.Points[iCp].Location +
            (pt_Target - pt_g0)
            )
        return pt_cp_1


    def translateCP_GrevilleToClosestPt(nc_toDeform, iCp, nc_Target, bDebug=False):
        """
        returns:
            Success: New NurbsCurve with deviation less than 1e-6.
            Fail: None
        """

        #Rhino.RhinoApp.CommandPrompt = (
        #        sCmdPrompt0 +
        #        "Searching for adjusted control points with less curve deviation ...")

        nc_toReturn = nc_toDeform.Duplicate()

        fWeight = 1.0#5

        pt_g0 = nc_toReturn.GrevillePoint(iCp)
        b, t = nc_Target.ClosestPoint(pt_g0)
        pt_Target = nc_Target.PointAt(t)
        dist = pt_Target.DistanceTo(pt_g0)
        if bDebug: sEval='dist'; print(sEval+': ',eval(sEval))
        if dist <= Rhino.RhinoMath.ZeroTolerance:
            return # No change.
        pt_cp_1 = (
            nc_toReturn.Points[iCp].Location +
            fWeight * (pt_Target - pt_g0)
            )
        nc_toReturn.Points.SetPoint(iCp, point=pt_cp_1)
        return nc_toReturn


        for i in xrange(100):
            sc.escape_test()
            pt_g0 = nc_toReturn.GrevillePoint(iCp)
            b, t = nc_Target.ClosestPoint(pt_g0)
            pt_Target = nc_Target.PointAt(t)
            dist = pt_Target.DistanceTo(pt_g0)
            if bDebug: sEval='dist'; print(sEval+': ',eval(sEval))
            if dist <= Rhino.RhinoMath.ZeroTolerance:
                if i == 0: return # No change.
                #print("{} iterations to adjust CP (index {}).".format(i, iCp))
                return nc_toReturn
            else:
                pt_cp_1 = (
                    nc_toReturn.Points[iCp].Location +
                    fWeight * (pt_Target - pt_g0)
                    )
                nc_toReturn.Points.SetPoint(iCp, point=pt_cp_1)
        else:
            print('#'*80)
            print("100 iterations in createNcWithTransCp_ByGrevClosestPt.  Check this!")
            print('#'*80)

        return nc_toReturn


    bWIP_was_adjusted = True # True so that while loop runs at least once.

    iFullAdjLoop = 1

    vs = [None]*nc_WIP.Points.Count # Translation vectors per CP.

    if bPreserveEndTans:
        # Unlike other vectors, these will not change during deformation of curves.
        vs[1] = nc_WIP.Points[1].Location - nc_WIP.Points[0].Location
        vs[1].Unitize()
        i = nc_WIP.Points.Count - 2
        vs[i] = nc_WIP.Points[i].Location - nc_WIP.Points[i+1].Location
        vs[i].Unitize()


    while iFullAdjLoop < 100:
        sc.escape_test()
        if bDebug: sEval='iFullAdjLoop'; print(sEval+': ',eval(sEval))

        bWIP_was_adjusted = False

        if bPreserveEndTans:
            # Move end tangent control points along tangency line.
            nc_AdjTan = fit_Adjust_G1_CP_spreads(
                    nc_In=nc_WIP,
                    nc_forDevComp=nc_forDevComp,
                    fResolution=fResolution,
                    iCallCt=iFullAdjLoop,
                    bDebug=bDebug)
            if nc_AdjTan is None:
                if bDebug: print("No G1 CP adjustment.")
            else:
                nc_WIP.Dispose()
                nc_WIP = nc_AdjTan
                bWIP_was_adjusted = True

                if debugHelp['gCrv'] is None:
                    debugHelp['gCrv'] = sc.doc.Objects.AddCurve(nc_WIP)
                    rdCrv = sc.doc.Objects.FindId(debugHelp['gCrv'])
                    rdCrv.GripsOn = True
                else:
                    sc.doc.Objects.Replace(objectId=debugHelp['gCrv'], curve=nc_WIP)
                
                sc.doc.Views.Redraw()

                res, s = ri.RhinoGet.GetString(
                    "Enter to continue / Esc to cancel",
                    acceptNothing=True,
                    outputString="")
                if res == Rhino.Commands.Result.Cancel: return

                #rc = Rhino.UI.Dialogs.ShowMessage(
                #    "Continue?",
                #    "Translate Control Points",
                #    buttons=Rhino.UI.ShowMessageButton.OKCancel,
                #    icon=Rhino.UI.ShowMessageIcon.None)
                #if rc != Rhino.UI.ShowMessageResult.OK:
                #    return

        # Translate controlPoints.
        if bPreserveEndTans:
            start = 2
            stop = nc_WIP.Points.Count-2
        else:
            start = 1
            stop = nc_WIP.Points.Count-1

        pts_TranslateTo = []

        for iCp in range(start, stop):
            pt = calculateNewCPLocation_GrevilleToClosestPt(
                nc_WIP, iCp, nc_forDevComp, bDebug=False)
            pts_TranslateTo.append(pt)

        if bDebug:
            [print(pt) for pt in pts_TranslateTo]

        if all(_ is None for _ in pts_TranslateTo):
            if not bWIP_was_adjusted:
                break

        for j, iCp in enumerate(range(start, stop)):
            point = pts_TranslateTo[j]
            if point is None: continue
            nc_WIP.Points.SetPoint(iCp, point=point)

        #sc.doc.Objects.AddCurve(nc_WIP); sc.doc.Views.Redraw()

        iFullAdjLoop += 1



        continue










        #for iCp in range(start, stop):
        #    if bDebug: sEval='iCp'; print(sEval+': ',eval(sEval),)

        #    nc_Trans = translateCP_GrevilleToClosestPt(
        #        nc_toDeform=nc_WIP,
        #        iCp=iCp,
        #        nc_Target=rgNurbsCrv_forDevComp,
        #        bDebug=bDebug)

        #    if bDebug:
        #        sEval='nc_Trans'; print(sEval+': ',eval(sEval))
        #        sc.doc.Objects.AddCurve(nc_Trans); sc.doc.Views.Redraw()
        #        #1/0

        #    if nc_Trans is not None:
        #        # Check deviation.
        #        fDev = getMaximumDeviation(rgNurbsCrv_forDevComp, nc_Trans)
        #        if bDebug: sEval='fDev'; print(sEval+': ',eval(sEval))
        #        if fDev <= (fDev0 + 1e-6):
        #            nc_WIP.Dispose()
        #            nc_WIP = nc_Trans
        #            bWIP_was_adjusted = True
        #            if bDebug:
        #                print("Reduced deviation by adjusting control point {}.".format(iCp))
        #        else:
        #            nc_Trans.Dispose()
        #            if bDebug:
        #                sEval='fDev > fDev0'; print(sEval+': ',eval(sEval))
        #                print("Morphing increased the deviation.")


        #iFullAdjLoop += 1
        
        #sc.doc.Objects.AddCurve(nc_WIP); sc.doc.Views.Redraw()#; return
    else:
        print("For loop completed.")

    print("After {} iterations.".format(iFullAdjLoop))

    if nc_WIP.EpsilonEquals(nc_toDeform, epsilon=fResolution):
        nc_WIP.Dispose()
        return

    return nc_WIP


def main():
    """
    """

    objref_ToMod = getInput()
    if objref_ToMod is None: return

    sc.doc.Objects.UnselectAll()

    objref_Ref = getInput(objref_ToMod)
    if objref_Ref is None: return

    fResolution = Opts.values['fResolution']
    bPreserveEndTans = Opts.values['bPreserveEndTans']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    edge_A = objref_ToMod.Edge()
    crv_A = objref_ToMod.Curve()
    crv_R = objref_Ref.Curve()


    if bDebug:
        Rhino.RhinoApp.ClearCommandHistoryWindow()


    Rhino.RhinoApp.CommandPrompt = "Working ..."

    nc_Res = fitCurve(
        crv_A,
        crv_R,
        fResolution=fResolution,
        bPreserveEndTans=bPreserveEndTans,
        bDebug=bDebug)


    if nc_Res is None:
        print("Curve was not created.")
        return

    if edge_A or not bReplace:
        # Only add curve.
        gOut = sc.doc.Objects.AddCurve(nc_Res)
        if gOut == gOut.Empty:
            print("Curve could not be added.")
        return

    if not sc.doc.Objects.Replace(objref_ToMod, curve=nc_Res):
        print("Curve could not be replaced.")

    return


    rc = getInput_OLD()
    if rc is None: return
    objrefs = rc[0]

    if Opts.values['bDebug']:
        pass

    processCurveObjects_OLD(
            rhCrvs0=objrefs,
    )


if __name__ == '__main__':
    #main()

    try:
        main()
    except ZeroDivisionError as e:
        print(e)
    except Exception as e:
        import traceback
        print('Traceback: {}'.format(traceback.format_exc()))
