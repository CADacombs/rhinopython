"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190920, 23-24: Created.
191021: Changed some variable names.
220823, 0908-10: Further development.
221007: Branched script off of another to focus on specific type of control point translation.

Notes:
Greville points aren't directly translated because any curve tangency is not maintained,
and results are squiggly.

No longer applicable for this script: It doesn't seem that Greville points can be adjusted within 1e-12 accuracy, so use a minimum of 1e-11.

TODO:
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


    key = 'iPreserveEndCont'; keys.append(key)
    values[key] = 1
    names[key] = 'PreserveEndCContinuity'
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fResolution'; keys.append(key)
    values[key] = 1e-6
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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


def knotMultiplicityPattern(nc):
    mp = []
    iK = 0
    while iK < nc.Knots.Count:
        m = nc.Knots.KnotMultiplicity(iK)
        mp.append(m)
        iK += m
    return mp


def getInput(objref_ToMod=None):
    """
    Get curve with optional input.
    """

    go = ri.Custom.GetObject()

    go.GeometryFilter = rd.ObjectType.Curve

    if objref_ToMod is None:
        go.SetCommandPrompt("Select NURBS curve to modify")
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

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('iPreserveEndCont')
        addOption('fResolution')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)

            if objref_ToMod is not None:
                go.Dispose()
                return objref


            # Verify that curve to modify is acceptable.

            nc = go.Object(0).Curve().DuplicateCurve() # In case Curve is a BrepEdge.
            if not isinstance(nc, rg.NurbsCurve):
                print("Not a NURBS curve.")
                nc.Dispose()
                sc.doc.Objects.UnselectAll()
                go.ClearObjects()
                sc.doc.Views.Redraw()
                continue

            iK = nc.Degree
            if any(m > 1 for m in knotMultiplicityPattern(nc)[1:-1]):
                print(
                    "NURBS curve has some interior knots with multiplicity > 1."
                    "This is not supported.")
                nc.Dispose()
                sc.doc.Objects.UnselectAll()
                go.ClearObjects()
                sc.doc.Views.Redraw()
                continue

            iPreserveEndCont = Opts.values['iPreserveEndCont']
            if iPreserveEndCont >= 0:
                if nc.Points.Count < (2*(iPreserveEndCont + 1) + 1):
                    print(
                        "NURBS curve has only {} control points.".format(nc.Points.Count),
                        "Need at least {} for parametric continuity preservation.".format(2*(iPreserveEndCont + 1) + 1))
                    nc.Dispose()
                    sc.doc.Objects.UnselectAll()
                    go.ClearObjects()
                    sc.doc.Views.Redraw()
                    continue

            # All good.

            go.Dispose()
            return objref

        if res == ri.GetResult.Number:
            key = 'iPreserveEndCont'
            Opts.riOpts[key].CurrentValue = int(go.Number())
            Opts.setValue(key)
            continue

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


def _OLD_CODE():

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

    pass


def fitCurve(nc_toDeform, nc_forDevComp, **kwargs):
    """
    returns:
        Success: New NurbsCurve with deviation less than nc0
        Fail: None
    """

    if not nc_toDeform.IsValid: return

    if nc_toDeform.Points.Count < 3: return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iPreserveEndCont = getOpt('iPreserveEndCont')
    fResolution = getOpt('fResolution')
    bDebug = getOpt('bDebug')


    if iPreserveEndCont >= 0:
        if nc_toDeform.Points.Count < (2*(iPreserveEndCont + 1) + 1):
            return


    def getMaxDev(rgCrvA, rgCrvB):
        rc = rg.Curve.GetDistancesBetweenCurves(
                rgCrvA,
                rgCrvB,
                tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
        if rc[0]:
            return rc[1]


    fDev_Start = getMaxDev(nc_forDevComp, nc_toDeform)
    if bDebug: sEval='fDev0'; print(sEval+': ',eval(sEval))
    if fDev_Start is None: return
    #if dev <= fDevTol:
    #    if bDebug: print("Deviation is already within tolerance.")
    #    # Deviation is out of tolerance,
    #    # so capture the previous curve.
    #    break

    nc_WIP = nc_toDeform.Duplicate() # Keep copy of original because nc_WIP will be modified per iteration.

    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt


    def calculateGrevilleToClosestPtVector(nc_toDeform, iCp, nc_Target, bDebug=False):
        pt_Gr_In = nc_toDeform.GrevillePoint(iCp)#; sc.doc.Objects.AddPoint(pt_Gr_In)
        b, t = nc_Target.ClosestPoint(pt_Gr_In)
        pt_Target = nc_Target.PointAt(t)#; sc.doc.Objects.AddPoint(pt_Target)
        dist = pt_Target.DistanceTo(pt_Gr_In)
        if bDebug: sEval='dist'; print(sEval+': ',eval(sEval))
        if dist <= max(1e-6, 0.1*sc.doc.ModelAbsoluteTolerance):
            return # No change.
        return pt_Target - pt_Gr_In


    iFullAdjLoop = 1


    if iPreserveEndCont < 0:
        start_idxCP = 0
        stop_idxCP = nc_WIP.Points.Count
    else:
        start_idxCP = iPreserveEndCont + 1
        stop_idxCP = nc_WIP.Points.Count - iPreserveEndCont -1


    ncs_Res = []
    devs = []


    while iFullAdjLoop < 1000:
        sc.escape_test()
        if bDebug: sEval='iFullAdjLoop'; print(sEval+': ',eval(sEval))

        # Calculate control point translations.
        vts_TranslateToward = []

        for iCp in range(start_idxCP, stop_idxCP):
            vt = calculateGrevilleToClosestPtVector(
                nc_WIP, iCp, nc_forDevComp, bDebug=False)
            vts_TranslateToward.append(vt)

        if bDebug:
            [print(vt) for vt in vts_TranslateToward]

        if all(_ is None for _ in vts_TranslateToward):
            print("No more control points need to be translated.")
            break

        iCt_toTrans = (len(vts_TranslateToward) - vts_TranslateToward.count(None))
        #print(iCt_toTrans)

        # TODO: Iterate various weight calculations and later choose curve with least deviation.
        # Why does 8.0 sometimes work well (so far)?  Goal is achieved faster than lower values,
        # but 16.0 is too large.
        weight = 2.00 / iCt_toTrans
        #print(weight)

        # Apply control point translations.
        for j, iCp in enumerate(range(start_idxCP, stop_idxCP)):
            vt = vts_TranslateToward[j]
            if vt is None:
                continue
            vt *= weight
            pt_Start = nc_WIP.Points[iCp].Location
            #sc.doc.Objects.AddLine(rg.Line(pt_Start, span=vt))
            pt_End = pt_Start + vt
            nc_WIP.Points.SetPoint(iCp, point=pt_End)

        dev = getMaxDev(nc_forDevComp, nc_WIP)

        #sc.doc.Views.Redraw()
        #return

        #sc.doc.Objects.AddCurve(nc_WIP); sc.doc.Views.Redraw()

        ncs_Res.append(nc_WIP.Duplicate())
        devs.append(dev)

        iFullAdjLoop += 1

    else:
        print("For loop completed.")

    print("After {} iterations.".format(iFullAdjLoop))

    print("Latest deviation: {}".format(formatDistance(devs[-1])))
    min_dev = min(devs)
    print("Minimum deviation: {}".format(formatDistance(min_dev)))

    count_min_dev = devs.count(min_dev)
    if count_min_dev > 1:
        print("Minimum deviation for {} curves.".format(count_min_dev))

    idx_Winner = devs.index(min_dev)

    nc_Out = ncs_Res[idx_Winner].Duplicate()

    for nc in ncs_Res: nc.Dispose()
    nc_WIP.Dispose()

    if nc_Out.EpsilonEquals(nc_toDeform, epsilon=fResolution):
        nc_Out.Dispose()
        return

    return nc_Out


def main():
    """
    """

    objref_ToMod = getInput()
    if objref_ToMod is None: return

    sc.doc.Objects.UnselectAll()

    objref_Ref = getInput(objref_ToMod)
    if objref_Ref is None: return

    iPreserveEndCont = Opts.values['iPreserveEndCont']
    fResolution = Opts.values['fResolution']
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
        iPreserveEndCont=iPreserveEndCont,
        fResolution=fResolution,
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
    main()

#    try:
#        main()
#    except ZeroDivisionError as e:
#        print(e)
#    except Exception as e:
#        import traceback
#        print('Traceback: {}'.format(traceback.format_exc()))
