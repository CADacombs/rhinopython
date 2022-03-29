"""
This script is an alternative to _RemoveMultiKnots on wire curves.

Like _RemoveMultiKnots, no multiplicity is reduced to 0.
Unlike _RemoveMultiKnots, it calculates various curves and chooses the one,
within deviation tolerance, if enabled, and with the least number of knots.


NurbsCurveKnotList.RemoveMultipleKnots' minimumMultiplicity and maximumMultiplicity are,
respectively,
the minimum and maximum allowed (not to remove) knot multiplicities.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
210603: Created starting from a split from another script.
220327: Modified main routine for finding curve to return.  Modified options.  Refactored.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bReduceIfRemoveFail'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iPreserveEnd'; keys.append(key)
    values[key] = 1
    listValues[key] = Enum.GetNames(rg.NurbsCurve.NurbsCurveEndConditionType)[1:]
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitDev'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bReplace'; keys.append(key)
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

        if key == 'fDevTol':
            if cls.riOpts[key].CurrentValue == 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = Rhino.RhinoMath.ZeroTolerance
            elif cls.riOpts[key].CurrentValue <= 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
        else:
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

    # Get curves with optional input.

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select wire curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bReduceIfRemoveFail')
        addOption('iPreserveEnd')
        addOption('bLimitDev')
        if Opts.values['bLimitDev']:
            addOption('fDevTol')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

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

        if res == ri.GetResult.Number:
            if Opts.values['bLimitDev']:
                key = 'fDevTol'
                Opts.riOpts[key].CurrentValue = go.Number()

            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def formatDistance(fDistance, iPrecision=None):
    if iPrecision is None: iPrecision = sc.doc.ModelDistanceDisplayPrecision

    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, iPrecision)


def knotMultiplicityList(knots):
    """Returns a list."""
    i = 0
    iMulties = []
    fKnotTs_Unique = []
    while True:
        knot = knots[i]
        fKnotTs_Unique.append(knot)
        iMulti = knots.KnotMultiplicity(index=i)
        iMulties.append(iMulti)
        #print("{} at {:.4f}".format(iMulti, knot),
        i += iMulti
        if i >= knots.Count:
            break
    return iMulties


def filterTargetInput(rgCrv):
    """
    Returns on non-target curve: None, str
    Returns on target curve: rg.NurbsCurve, None
    """

    if isinstance(rgCrv, rg.NurbsCurve):
        nc_In = rgCrv
    elif isinstance(rgCrv, rg.PolyCurve):
        nc_In = rgCrv.ToNurbsCurve()
    else:
        return None, "{} skipped.".format(rgCrv.GetType().Name)

    knots = nc_In.Knots
    degree = nc_In.Degree

    # Test whether any knots can be removed by RemoveMultipleKnots.
    if nc_In.IsPeriodic:
        if knots.Count == 2*degree + degree -1:
            return None, "Curve has no interior knots that can be removed by NurbsCurveKnotList.RemoveMultipleKnots."
    elif knots.Count == 2 * degree:
        return None, "Curve has no interior knots that can be removed by NurbsCurveKnotList.RemoveMultipleKnots."

    iKnot = degree
    iKnot_Stop = knots.Count - degree

    while iKnot < iKnot_Stop:
        iMult = knots.KnotMultiplicity(iKnot)

        if iMult > 1:
            return nc_In, None

        iKnot += 1 # or iMult.

    return None, "Curve has no interior knots that can be removed by NurbsCurveKnotList.RemoveMultipleKnots."


def getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        return None

    return rc[1]


def processCurve(rgCrv, **kwargs):
    """
    Parameters:
        rgCrv
        fDevTol:
            None for no deviation limit,
            Positive value for limit,
            Negative value for removing as many knots as possible within absolute value.
        bDebug
        
    Initial tests show that removing all knots produces the same result as an
    iterative approach, e.g., remove knots starting at full multiplicity.
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bReduceIfRemoveFail = getOpt('bReduceIfRemoveFail')
    iPreserveEnd = getOpt('iPreserveEnd')
    fDevTol = getOpt('fDevTol') if getOpt('bLimitDev') else None
    bDebug = getOpt('bDebug')


    rc = filterTargetInput(rgCrv)
    if rc[0] is None:
        return rc

    nc_In = rc[0]

    degree = nc_In.Degree

    #iMinMultiplicity = 1 # Remove all polyknots.


    def preserveEndCondition(nc):
        if iPreserveEnd == 1:
            bSuccess = nc.SetEndCondition(
                bSetEnd=False,
                continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
                point=nc.PointAtStart,
                tangent=nc_In.TangentAtStart)
            if not bSuccess:
                print("SetEndCondition failed.")
            bSuccess = nc.SetEndCondition(
                bSetEnd=True,
                continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
                point=nc.PointAtEnd,
                tangent=nc_In.TangentAtEnd)
            if not bSuccess:
                print("SetEndCondition failed.")
        elif iPreserveEnd == 2:
            bSuccess = nc.SetEndCondition(
                bSetEnd=False,
                continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                point=nc.PointAtStart,
                tangent=nc_In.TangentAtStart,
                curvature=nc_In.CurvatureAt(nc_In.Domain.T0))
            if not bSuccess:
                print("SetEndCondition failed.")
            bSuccess = nc.SetEndCondition(
                bSetEnd=True,
                continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                point=nc.PointAtEnd,
                tangent=nc_In.TangentAtEnd,
                curvature=nc_In.CurvatureAt(nc_In.Domain.T1))
            if not bSuccess:
                print("SetEndCondition failed.")


    # Not setting tolerance parameter to other than Rhino.RhinoMath.UnsetValue
    # due to bugs in RC 6.x - 7.6 (+):
    #   In pre-7.6: Float value is ignored.
    #   In 7.6: Return value does not reflect actual number of knots removed.
    #     https://discourse.mcneel.com/t/nurbscurveknotlist-removemultipleknots-return-value-bug/125314
    # For V6 compatibility, deviation results will be checked instead.


    # RemoveMultipleKnots parameters:
    #   minimumMultiplicity is the lowest value not to change.
    #   maximumMultiplicity is the highest value not to change.

    # Setting minimumMultiplicity to a value less than 1 is the same as setting it to 1.



    if fDevTol is None:
        nc_WIP = nc_In.DuplicateCurve()

        iCt_KnotsRemoved = nc_WIP.Knots.RemoveMultipleKnots(
            minimumMultiplicity=1,
            maximumMultiplicity=degree+1,
            tolerance=Rhino.RhinoMath.UnsetValue)

        preserveEndCondition(nc_WIP)

        fDev = getDistancesBetweenCurves(nc_In, nc_WIP)

        if fDev is None:
            return (nc_WIP, None), "GetDistancesBetweenCurves failed."

        return (nc_WIP, fDev), None


    if not bReduceIfRemoveFail:
        nc_WIP = nc_In.DuplicateCurve()

        iCt_KnotsRemoved = nc_WIP.Knots.RemoveMultipleKnots(
            minimumMultiplicity=1,
            maximumMultiplicity=degree+1,
            tolerance=Rhino.RhinoMath.UnsetValue)

        preserveEndCondition(nc_WIP)

        fDev = getDistancesBetweenCurves(nc_In, nc_WIP)

        if fDev is None:
            nc_WIP.Dispose()
            return None, "GetDistancesBetweenCurves failed."

        if fDev > fDevTol:
            nc_WIP.Dispose()
            return (None, fDev), None

        return (nc_WIP, fDev), None



    ncs_Pass = []
    fDevs = []
    iCts_Knots = []


    for minimumMultiplicity in range(1, degree):
        for maximumMultiplicity in range(degree+1, minimumMultiplicity+1, -1):

            nc_WIP = nc_In.DuplicateCurve()

            if bDebug:
                sEval='minimumMultiplicity'; print("{}: {}".format(sEval, eval(sEval)))
                sEval='maximumMultiplicity'; print("{}: {}".format(sEval, eval(sEval)))

            iCt_KnotsRemoved = nc_WIP.Knots.RemoveMultipleKnots(
                minimumMultiplicity,
                maximumMultiplicity,
                tolerance=Rhino.RhinoMath.UnsetValue)


            if bDebug: sEval='iCt_KnotsRemoved'; print("{}: {}".format(sEval, eval(sEval)))

            if iCt_KnotsRemoved == 0:
                nc_WIP.Dispose()
                if bDebug: print("No knots were removed.")
                continue

            preserveEndCondition(nc_WIP)

            fDev = getDistancesBetweenCurves(nc_In, nc_WIP)

            if bDebug: sEval='fDev'; print("{}: {}".format(sEval, eval(sEval)))

            if fDev is None:
                if bDebug: print("GetDistancesBetweenCurves failed.")
                nc_WIP.Dispose()
                continue


            if fDev > fDevTol:
                nc_WIP.Dispose()
                continue

            ncs_Pass.append(nc_WIP)
            fDevs.append(fDev)
            iCts_Knots.append(nc_WIP.Knots.Count)

            if bDebug:
                print(
                    "  Dev:{}".format(formatDistance(fDev)),
                    "  KnotCt:{}".format(nc_WIP.Knots.Count),
                    "  Multies:({})".format(
                        ",".join(str(i) for i in knotMultiplicityList(nc_WIP.Knots)))
                    )

    if bDebug: print(iCts_Knots)

    if not ncs_Pass:
        return None, "Curve could not be modified within tolerance."

    idx_Winner = iCts_Knots.index(min(iCts_Knots))

    for i, nc in enumerate(ncs_Pass):
        if i == idx_Winner:
            nc_Out = nc
            continue # for loop.
        nc.Dispose()


    sOut = None


    if nc_In.IsPeriodic and not nc_Out.IsPeriodic:
        sOut = "Periodic curve was converted to non-periodic." \
              "  This is a limitation of NurbsCurveKnotList.RemoveMultipleKnots."


    return (nc_Out, fDevs[idx_Winner]), sOut


def processCurveObjects(objrefs_In, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bReduceIfRemoveFail = getOpt('bReduceIfRemoveFail')
    iPreserveEnd = getOpt('iPreserveEnd')
    fDevTol = getOpt('fDevTol') if getOpt('bLimitDev') else None
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    gCrvs_Out = []
    fDevs = []
    sLogs = []

    for objref_In in objrefs_In:
        rdCrv = objref_In.Object()
        rgCrv_In = rdCrv.CurveGeometry
        nc_ToMod = rgCrv_In.ToNurbsCurve()

        #if not isinstance(rgCrv_In, rg.NurbsCurve):
        #    print("{} skipped.".format(rgCrv_In.GetType().Name)
        #    continue

        rc = processCurve(
            rgCrv=nc_ToMod,
            bReduceIfRemoveFail=bReduceIfRemoveFail,
            iPreserveEnd=iPreserveEnd,
            fDevTol=fDevTol,
            bDebug=bDebug,
            )

        if rc[0] is None:
            sLogs.append(rc[1])
            continue

        nc_Ret, fDev = rc[0]

        fDevs.append(fDev)

        if bReplace:
            if sc.doc.Objects.Replace(objectId=rdCrv.Id, curve=nc_Ret):
                gCrvs_Out.append(rdCrv.Id)
                sLogs.append("Curve was replaced.")
            else:
                sLogs.append("Curve could not be replaced.")
                continue
        else:
            gC_Out = sc.doc.Objects.AddCurve(nc_Ret)
            if gC_Out == gC_Out.Empty:
                sLogs.append("CurveObject could not be added.")
                continue
            gCrvs_Out.append(gC_Out)
            sLogs.append("Curve was added.")

    if bEcho:
        if len(gCrvs_Out) == 1:
            print(sLogs[0])
            s = "KnotMultiplicities:({})->({})".format(
                ",".join(str(i) for i in knotMultiplicityList(nc_ToMod.Knots)),
                ",".join(str(i) for i in knotMultiplicityList(nc_Ret.Knots)))
            s += "  Deviation:{}".format(formatDistance(fDev))
            print(s)
        elif len(sLogs) == 1:
            print(sLogs[0])
        else:
            for sLog in set(sLogs):
                print("Ct:{} {}".format(sLogs.count(sLog), sLog))
            if fDevs:
                print("Devs:[{},{}]".format(
                    formatDistance(min(fDevs)),
                    formatDistance(max(fDevs))))

    return gCrvs_Out


def main():

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script is supported only in Rhino V6 and above.")
        return

    objrefs_In = getInput()
    if objrefs_In is None: return

    bReduceIfRemoveFail = Opts.values['bReduceIfRemoveFail']
    iPreserveEnd = Opts.values['iPreserveEnd']
    fDevTol = Opts.values['fDevTol'] if Opts.values['bLimitDev'] else None
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObjects(
        objrefs_In=objrefs_In,
        bReduceIfRemoveFail=bReduceIfRemoveFail,
        iPreserveEnd=iPreserveEnd,
        fDevTol=fDevTol,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if gC_Res is None:
        return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()