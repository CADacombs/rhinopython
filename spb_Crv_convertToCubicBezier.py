"""
End conditions will automatically be tangency-matched to input curve.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221005: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select reference curve")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    def customGeometryFilter(rdObj, rgObj, compIdx):
        if isinstance(rgObj, rg.BrepEdge):
            crv = rgObj.DuplicateCurve()
        else:
            crv = rgObj
        if not isinstance(crv, rg.NurbsCurve):
            return True
        nc = crv
        if nc.Degree == 3 and nc.SpanCount == 1 and not nc.IsRational:
            #print("Curve is already a non-rational, cubic Bezier.")
            return False
        return True
    go.SetCustomGeometryFilter(customGeometryFilter)


    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

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


def _formatDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    #if fDistance < 0.001:
    #    return "{:.2e}".format(fDistance)
    #else:
    return "{:.{}f}".format(fDistance, max(6, sc.doc.ModelDistanceDisplayPrecision))


def _getMaxDev(rgCrvA, rgCrvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            rgCrvA,
            rgCrvB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
    if rc[0]:
        return rc[1]


def createCurve(crv_In, **kwargs):
    """
    returns:
        Success: New NurbsCurve with deviation less than nc0
        Fail: None
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if isinstance(crv_In, rg.LineCurve):
        nc_Out = crv_In.Rebuild(pointCount=4, degree=3, preserveTangents=False)
        if nc_Out:
            print("Input is a LineCurve.")
            return nc_Out
        return
    if crv_In.IsLinear(Rhino.RhinoMath.ZeroTolerance):
        nc_Out = rg.LineCurve(crv_In.PointAtStart, crv_In.PointAtEnd).Rebuild(pointCount=4, degree=3, preserveTangents=False)
        if nc_Out:
            print("Input is a linear {}.".format(crv_In.GetType().Name))
            return nc_Out
        return

    if isinstance(crv_In, rg.NurbsCurve) and crv_In.Degree==2 and crv_In.SpanCount==1 and not crv_In.IsRational:
        nc_Out = crv_In.DuplicateCurve()
        if rg.NurbsCurve.IncreaseDegree(nc_Out, 3):
            print("Input is a degree-2 NurbsCurve.")
            return nc_Out
        return

    if isinstance(crv_In, rg.BrepEdge):
        print("Input is a BrepEdge of a {}.".format(crv_In.DuplicateCurve().GetType().Name))
    else:
        print("Input is a {}.".format(crv_In.GetType().Name))



    arc = None

    if isinstance(crv_In, rg.ArcCurve):
        arc = rg.ArcCurve.Arc
    rc = rg.Curve.TryGetArc(crv_In, tolerance=Rhino.RhinoMath.ZeroTolerance)
    if rc[0]:
        arc = rc[1]

    if arc is not None:
        
        def createPoints_matchArcAtMidPoint(arc):
        
            if (arc.Angle - Rhino.RhinoMath.ZeroTolerance) > math.pi:
                return
        
        
            # Per https://pomax.github.io/bezierinfo/#circles_cubic
            k = 4.0 * math.tan(arc.Angle/4.0) / 3.0
        
            return Rhino.Collections.Point3dList(
                arc.StartPoint,
                arc.StartPoint + k * arc.TangentAt(arc.AngleDomain.T0) * arc.Radius,
                arc.EndPoint   - k * arc.TangentAt(arc.AngleDomain.T1) * arc.Radius,
                arc.EndPoint)

        point3dList = createPoints_matchArcAtMidPoint(arc)

        return rg.Curve.CreateControlPointCurve(point3dList, degree=3)



    fBulgeUnitDist = crv_In.PointAtStart.DistanceTo(crv_In.PointAtEnd)

    pts = [crv_In.PointAtStart, None, None, crv_In.PointAtEnd]

    vA = crv_In.TangentAtStart # Already length of 1.0
    vB = crv_In.TangentAtEnd # Already length of 1.0


    def findCrv_Symmetrical():

        fmin = 0.0
        fmax = 1.0

        min_devs = None

        fDivs = res_WIP = 0.1

        ncs_Res = []

        while True:
            sc.escape_test()

            for nc in ncs_Res: nc.Dispose()

            if bDebug:
                sEval = "fmin"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "fmax"; print("{}: {}".format(sEval, eval(sEval)))


            ncs_Res = []
            mABs = []
            devs = []

            for iAB in range(int(1.0/fDivs) + 1):
                mAB = fmin * (1 - iAB * fDivs) + fmax * (iAB * fDivs)
                if mAB == 0.0: continue
                pts[1] = pts[0] + mAB * fBulgeUnitDist * vA
                pts[2] = pts[3] - mAB * fBulgeUnitDist * vB

                nc_WIP = rg.NurbsCurve.CreateControlPointCurve(pts, degree=3)

                #sc.doc.Objects.AddCurve(nc_WIP)

                dev = _getMaxDev(nc_WIP, crv_In)
                if not dev:
                    if bDebug:
                        print("Deviation could not be determined for mA, mB:{}".format(mAB))
                    #sc.doc.Objects.AddCurve(nc_WIP); 1/0
                    continue

                ncs_Res.append(nc_WIP)
                mABs.append(mAB)
                devs.append(dev)

            if not ncs_Res:
                print("No curves were generated.")
                fmin *= 10.0
                fmax *= 10.0
                continue


            min_devs = min(devs)

            s_min_devs = _formatDistance(min_devs)

            if bDebug:
                print("Minimum deviation: {}".format(s_min_devs))

            if devs.count(min_devs) > 1:
                print("More than one curve with deviation {}.".format(s_min_devs))

            idx_Winner = devs.index(min_devs)

            if bDebug:
                sEval = "mABs[idx_Winner]"; print("{}: {}".format(sEval, eval(sEval)))

            if res_WIP < 1.1e-6:
                return ncs_Res[idx_Winner], min_devs

            fmin = mABs[idx_Winner] - res_WIP
            fmax = mABs[idx_Winner] + res_WIP

            res_WIP *= 0.1

            fDivs = 0.05


        return ncs_Res[idx_Winner], min_devs


    def findCrv_NonSymmetrical():

        fmin_mA = fmin_mB = 0.0
        fmax_mA = fmax_mB = 1.0

        min_devs_Prev = min_devs = None

        fDivs = res_WIP = 0.1

        ncs_Res = []

        while True:
            sc.escape_test()

            for nc in ncs_Res: nc.Dispose()

            if bDebug:
                sEval = "fmin_mA"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "fmax_mA"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "fmin_mB"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "fmax_mB"; print("{}: {}".format(sEval, eval(sEval)))


            ncs_Res = []
            mAs = []
            mBs = []
            devs = []

            for iA in range(int(1.0/fDivs) + 1):
                mA = fmin_mA * (1 - iA * fDivs) + fmax_mA * (iA * fDivs)
                if mA == 0.0: continue
                pts[1] = pts[0] + mA * fBulgeUnitDist * vA
                for iB in range(int(1.0/fDivs) + 1):
                    mB = fmin_mB * (1 - iB * fDivs) + fmax_mB * (iB * fDivs)
                    if mB == 0.0: continue
                    pts[2] = pts[3] - mB * fBulgeUnitDist * vB

                    nc_WIP = rg.NurbsCurve.CreateControlPointCurve(pts, degree=3)

                    #sc.doc.Objects.AddCurve(nc_WIP)

                    dev = _getMaxDev(nc_WIP, crv_In)
                    if not dev:
                        if bDebug:
                            print("Deviation could not be determined for mA:{} , mB:{}".format(mA, mB))
                        #sc.doc.Objects.AddCurve(nc_WIP); 1/0
                        continue

                    ncs_Res.append(nc_WIP)
                    mAs.append(mA)
                    mBs.append(mB)
                    devs.append(dev)

            if not ncs_Res:
                print("No curves were generated.")
                fmin_mA *= 10.0
                fmax_mA *= 10.0
                fmin_mB *= 10.0
                fmax_mB *= 10.0
                continue


            min_devs = min(devs)

            s_min_devs = _formatDistance(min_devs)

            if bDebug:
                print("Minimum deviation: {}".format(s_min_devs))

            if devs.count(min_devs) > 1:
                print("More than one curve with deviation {}.".format(s_min_devs))

            idx_Winner = devs.index(min_devs)

            if bDebug:
                sEval = "mAs[idx_Winner]"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "mBs[idx_Winner]"; print("{}: {}".format(sEval, eval(sEval)))

            if min_devs_Prev is not None and abs(min_devs_Prev - min_devs) < 1e-6: #0.1*sc.doc.ModelAbsoluteTolerance:
                return ncs_Res[idx_Winner], min_devs

            min_devs_Prev = min_devs

            if res_WIP < 0.00101 and abs(mAs[idx_Winner] - mBs[idx_Winner]) <= 0.001:
                if bDebug: print("Go with symmetrical solution.")
                return

            fmin_mA = mAs[idx_Winner] - res_WIP
            fmax_mA = mAs[idx_Winner] + res_WIP
            fmin_mB = mBs[idx_Winner] - res_WIP
            fmax_mB = mBs[idx_Winner] + res_WIP

            res_WIP *= 0.1

            fDivs = 0.05

        return ncs_Res[idx_Winner], min_devs


    if bDebug: print("Symmetrical search:")
    rc = findCrv_Symmetrical()
    (nc_Sym, dev_Sym) = rc if rc else (None, None)

    if bDebug: print("Non-symmetrical search:")
    rc = findCrv_NonSymmetrical()
    (nc_NonSym, dev_NonSym) = rc if rc else (None, None)

    if dev_Sym is None and dev_NonSym is None:
        return

    if dev_Sym is None:
        if bEcho:
            print("Non-symmetrical solution with deviation of {} from input.".format(_formatDistance(dev_NonSym)))
        return nc_NonSym

    if dev_NonSym is None:
        if bEcho:
            print("Symmetrical solution with deviation of {} from input.".format(_formatDistance(dev_Sym)))
        return nc_Sym

    if dev_Sym <= dev_NonSym:
        if bEcho:
            print("Symmetrical solution with deviation of {} from input.".format(_formatDistance(dev_Sym)))
        return nc_Sym

    if bEcho:
        print("Non-symmetrical solution with deviation of {} from input.".format(_formatDistance(dev_NonSym)))
    return nc_NonSym


def main():
    """
    """

    objref_In = getInput()
    if objref_In is None: return

    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    crv_In = objref_In.Curve()


    #if bDebug:
    #    Rhino.RhinoApp.ClearCommandHistoryWindow()


    Rhino.RhinoApp.CommandPrompt = "Working ..."

    nc_Res = createCurve(
        crv_In,
        bEcho=bEcho,
        bDebug=bDebug)


    if nc_Res is None:
        print("Curve was not created.")
        return


    if objref_In.Edge() or not bReplace:
        gOut = sc.doc.Objects.AddCurve(nc_Res)
        if gOut != gOut.Empty:
            print("Curve was added.")
            sc.doc.Views.Redraw()
            return
        else:
            print("Curve could not be added.")
            return


    if sc.doc.Objects.Replace(objref_In, curve=nc_Res):
        print("Curve was replaced.")
        sc.doc.Views.Redraw()
    else:
        print("Curve could not be replaced.")


if __name__ == '__main__': main()