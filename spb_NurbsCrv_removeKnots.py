"""
This script removes curves knots in a few interative ways,
keeping results with least deviation from input curves.

Observation:
    Removing multiple knots at one time at a multiplicity sometimes results in
    less deviation than removing one knot at a time, a la _RemoveKnot.
    This means that _RemoveKnot is not optimal.


Why this script doesn't work on Periodic curves:
    RemoveKnots converts periodic curves to non-periodic because it modifies
    the overlapping control point locations independently.


Why NurbsCurveKnotList.RemoveMultipleKnots isn't used:
#1 reason: RemoveMultipleKnots converts all knots as possible to reduce multiplicity to 1.
This may not be necessary for, say, a degree-5 curve that only needs to contain
multiplicities to guarantee G2 (maximum multiplicity of 3 in this case).
Still, RemoveMultipleKnots doesn't seem to respect tolerance parameter.
This may be problematic only when degree != 3.
Regardless, RemoveMultipleKnots seems to work the same as removing
knots starting at the end of the domain working toward the start.  No advantage was
witnessed using this method.
This bug was corrected: https://discourse.mcneel.com/t/nurbscurveknotlist-removemultipleknots-return-value-bug/125314
but tolerance value has been found to not be honored.


TODO:
    Investigate how periodic curves can be supported.  (See notes above.)  Try RemoveKnotAt.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
170819-23: Created.
...
210530-0604: Script revamped.
210712: Implemented RemoveMultipleKnots.  Modified an option default value.
220326: Now loops through continuities to remove from (to removes all internal knots) through target min. continuity.
        Disabled use of RemoveKnotAt since its results seem to just duplicate RemoveKnots (singular removal).
        Refactored.
220327: Added function to remove all knots in one go.
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


    key = 'iMinTargetGContinuity'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(values[key])
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
    #fModelToInches = Rhino.RhinoMath.UnitScale(
    #     sc.doc.ModelUnitSystem, Rhino.UnitSystem.Inches)
    #fInchesToModel = Rhino.RhinoMath.UnitScale(
    #     Rhino.UnitSystem.Inches, sc.doc.ModelUnitSystem)
    #values[key] = max(
    #    0.1 * sc.doc.ModelAbsoluteTolerance,
    #    (sc.doc.ModelAbsoluteTolerance * fModelToInches - 0.0001) * fInchesToModel
    #    )
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bRemoveAsManyAsPossibleWithinTol'; keys.append(key)
    values[key] = True
    names[key] = 'RemoveKnots'
    riOpts[key] = ri.Custom.OptionToggle(
        values[key], 'AllOrNothing', 'AsManyAsPossibleWithinTol')
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
        elif key == 'iMinTargetGContinuity':
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

    print("MinTargetGContinuity <= 0 will use the curve degree as the minimum.")

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iMinTargetGContinuity')
        addOption('iPreserveEnd')
        addOption('bLimitDev')
        if Opts.values['bLimitDev']:
            addOption('fDevTol')
            addOption('bRemoveAsManyAsPossibleWithinTol')
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
            else:
                key = 'iMinTargetGContinuity'
                Opts.riOpts[key].CurrentValue = int(go.Number())

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


def filterTargetInput(rgCrv, iMinTargetGContinuity):
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

    # Return if no knots to remove.
    if nc_In.IsPeriodic:
        if knots.Count == 2*degree + degree -1:
            return None, "Minimum-knot periodic curve skipped."
        return (
            None,
            "Periodic curves not supported yet."
            "  Workarounds: Split curve then rerun this script or use _RemoveKnot."
            )

    if knots.Count == 2 * degree:
        return None, "Bezier curve skipped."


    if iMinTargetGContinuity is None:
        return nc_In, None

    # Further test whether any knots need to be removed.
    knotCtMin = degree - iMinTargetGContinuity

    #if nc_In.IsPeriodic:
    #    iKnot = 0
    #    iKnot_Stop = knots.Count
    #else:
    iKnot = degree
    iKnot_Stop = knots.Count - degree

    while iKnot < iKnot_Stop:
        iMult = knots.KnotMultiplicity(iKnot)

        if iMult > knotCtMin:
            return nc_In, None

        iKnot += iMult

    return None, "Curve has no interior knots that need to be removed."


def getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        return None

    return rc[1]


def removeKnotAt_StartToEnd(nc_In, iMinMultiplicity):
    nc_Out = nc_In.DuplicateCurve()

    knots = nc_Out.Knots

    iK = nc_In.Degree

    bAnyRemoved = False

    while iK < knots.Count - nc_In.Degree: #iK >= nc_In.Degree:

        #print(iK, knots.KnotMultiplicity(iK))

        m = knots.KnotMultiplicity(iK)

        if m > iMinMultiplicity:
            bRemoved = knots.RemoveKnotAt(knots[iK])
            if bRemoved:
                bAnyRemoved = True
                #print("RemoveKnotAt {}".format(iK_toRemove)
            else:
                print("RemoveKnotAt failed.")
                nc_Out.Dispose()
                return

            continue

        iK += 1

    if bAnyRemoved:
        return nc_Out

    nc_Out.Dispose()


def removeKnotAt_EndToStart(nc_In, iMinMultiplicity):
    nc_Out = nc_In.DuplicateCurve()

    knots = nc_Out.Knots

    iK = knots.Count - nc_In.Degree - 1

    bAnyRemoved = False

    while iK >= nc_In.Degree:

        #print(iK, knots.KnotMultiplicity(iK)

        m = knots.KnotMultiplicity(iK)

        if m > iMinMultiplicity:
            iK_ToRemove_Start = iK
            iK_ToRemove_Stop = iK - m + iMinMultiplicity
            #sEval = 'iK_ToRemove_Stop'; print(sEval+':',eval(sEval)
            for iK_toRemove in range(iK_ToRemove_Start, iK_ToRemove_Stop, -1):
                bRemoved = knots.RemoveKnotAt(knots[iK_toRemove])
                if bRemoved:
                    bAnyRemoved = True
                    #print("RemoveKnotAt {}".format(iK_toRemove)
                else:
                    print("RemoveKnotAt fail.")
                    nc_Out.Dispose()
                    return

            iK = iK_ToRemove_Stop
            continue

        iK -= 1

    if bAnyRemoved:
        return nc_Out

    nc_Out.Dispose()


def removeKnots_All_possible(nc_In):
    """
    """

    nc_Out = nc_In.DuplicateCurve()
    degree = nc_Out.Degree

    if nc_In.IsPeriodic:
        index0 = degree + 1
        index1 = nc_Out.Knots.Count - nc_Out.Degree - 1
    else:
        index0 = degree
        index1 = nc_Out.Knots.Count - nc_Out.Degree

    nc_Out.Knots.RemoveKnots(index0, index1)

    return nc_Out


def removeKnots_Plural_StartToEnd(nc_In, iMinTargetGContinuity, fDev_Tol=None):
    """
    Plural means multiple knots may be removed per call to RemoveKnots.
    """

    nc_Out = nc_In.DuplicateCurve()
    degree = nc_Out.Degree
    iK = degree
    iMinMultiplicity = max((0, (degree - iMinTargetGContinuity)))

    while iK < nc_Out.Knots.Count - degree:
        # nc_Out.Knots.Count - degree is the stop index.
        sc.escape_test()

        m = nc_Out.Knots.KnotMultiplicity(iK)

        if m <= iMinMultiplicity:
            iK += m
            continue

        index0 = iK # Index of first knot to remove.
        index1 = iK + m - iMinMultiplicity # Stop index (last index + 1) to remove.
        #sEval = 'm'; print(sEval+':',eval(sEval))
        #sEval = 'index0'; print(sEval+':',eval(sEval))
        #sEval = 'index1'; print(sEval+':',eval(sEval))

        if fDev_Tol is None:
            bRemoved = nc_Out.Knots.RemoveKnots(index0, index1)

            if not bRemoved:
                print("RemoveKnots fail.")
                nc_Out.Dispose()
                return

            iK += m - (index1-index0)
            continue # while loop.


        # Remove within fDev_Tol.

        nc_Out_Save = nc_Out.DuplicateCurve()

        while True:
            sc.escape_test()

            bRemoved = nc_Out.Knots.RemoveKnots(index0, index1)

            if not bRemoved:
                print("RemoveKnots fail.")
                nc_Out_Save.Dispose()
                nc_Out.Dispose()
                return

            fDev = getDistancesBetweenCurves(nc_In, nc_Out)

            if fDev <= fDev_Tol:
                break # out of while loop.

            nc_Out.Dispose()
            nc_Out = nc_Out_Save.DuplicateCurve()
            index1 -= 1
            if index1 == index0:
                break # out of while loop.

        m = nc_Out.Knots.KnotMultiplicity(iK)
        iK += m

        nc_Out_Save.Dispose()


    if nc_Out.Knots.Count < nc_In.Knots.Count:
        #sc.doc.Objects.AddCurve(nc_Out); sc.doc.Views.Redraw(); 1/0
        return nc_Out

    nc_Out.Dispose()


def removeKnots_Plural_EndToStart(nc_In, iMinTargetGContinuity, fDev_Tol=None):
    """
    Plural means multiple knots may be removed per call to RemoveKnots.
    """

    nc_Out = nc_In.DuplicateCurve()
    degree = nc_Out.Degree
    iK = nc_Out.Knots.Count - degree - 1
    iMinMultiplicity = max((0, (degree - iMinTargetGContinuity)))

    while iK >= degree:
        # degree is the last index to process.
        sc.escape_test()

        m = nc_Out.Knots.KnotMultiplicity(iK)

        if m <= iMinMultiplicity:
            iK -= m
            continue

        index0 = iK - m + iMinMultiplicity + 1 # Index of first knot to remove.
        index1 = iK + 1 # Stop index, not the last to remove.
        #sEval = 'm'; print(sEval+':',eval(sEval))
        #sEval = 'index0'; print(sEval+':',eval(sEval)
        #sEval = 'index1'; print(sEval+':',eval(sEval)

        if fDev_Tol is None:
            bRemoved = nc_Out.Knots.RemoveKnots(index0, index1)

            if not bRemoved:
                print("RemoveKnots fail.")
                nc_Out.Dispose()
                return

            iK -= m
            continue # while loop.


        # Remove within fDev_Tol.

        nc_Out_Save = nc_Out.DuplicateCurve()

        while True:
            sc.escape_test()

            bRemoved = nc_Out.Knots.RemoveKnots(index0, index1)

            if not bRemoved:
                print("RemoveKnots fail.")
                nc_Out_Save.Dispose()
                nc_Out.Dispose()
                return

            fDev = getDistancesBetweenCurves(nc_In, nc_Out)

            if fDev <= fDev_Tol:
                break # out of while loop.

            nc_Out.Dispose()
            nc_Out = nc_Out_Save.DuplicateCurve()
            index0 += 1
            if index0 == index1:
                break # out of while loop.

        m = nc_Out.Knots.KnotMultiplicity(iK)
        iK -= m

        nc_Out_Save.Dispose()


    if nc_Out.Knots.Count < nc_In.Knots.Count:
        #sc.doc.Objects.AddCurve(nc_Out); sc.doc.Views.Redraw(); 1/0
        return nc_Out

    nc_Out.Dispose()


def removeKnots_Single_StartToEnd(nc_In, iMinTargetGContinuity, fDev_Tol=None):
    """
    Single means only one knot is removed for each call to RemoveKnots.
    """

    nc_Out = nc_In.DuplicateCurve()
    degree = nc_Out.Degree
    iK = nc_In.Degree
    iMinMultiplicity = max((0, (degree - iMinTargetGContinuity)))

    while iK < nc_Out.Knots.Count - degree:
        # nc_Out.Knots.Count - degree is the stop index.
        sc.escape_test()

        m = nc_Out.Knots.KnotMultiplicity(iK)

        if m <= iMinMultiplicity:
            iK += m
            continue

        #sEval = 'iK'; print(sEval+':',eval(sEval))
        #sEval = 'm'; print(sEval+':',eval(sEval))

        if fDev_Tol is None:
            bRemoved = nc_Out.Knots.RemoveKnots(iK, iK+1)

            if not bRemoved:
                print("RemoveKnots fail.")
                nc_Out.Dispose()
                return

            continue # while loop.


        nc_Out_Save = nc_Out.DuplicateCurve()

        bRemoved = nc_Out.Knots.RemoveKnots(iK, iK+1)

        if not bRemoved:
            print("RemoveKnots fail.")
            nc_Out_Save.Dispose()
            nc_Out.Dispose()
            return

        fDev = getDistancesBetweenCurves(nc_In, nc_Out)

        #sEval = 'fDev'; print(sEval+':',eval(sEval))

        if fDev <= fDev_Tol:
            iK += 1
            nc_Out_Save.Dispose()
            continue # while loop.

        nc_Out.Dispose()
        nc_Out = nc_Out_Save.DuplicateCurve()

        iK += 1

    if nc_Out.Knots.Count < nc_In.Knots.Count:
        #sc.doc.Objects.AddCurve(nc_Out); sc.doc.Views.Redraw(); 1/0
        return nc_Out

    nc_Out.Dispose()


def removeKnots_Single_EndToStart(nc_In, iMinTargetGContinuity, fDev_Tol=None):
    """
    Single means only one knot is removed for each call to RemoveKnots.
    """

    nc_Out = nc_In.DuplicateCurve()
    degree = nc_Out.Degree
    iK = nc_Out.Knots.Count - nc_Out.Degree - 1
    iMinMultiplicity = degree - iMinTargetGContinuity

    while iK >= nc_Out.Degree:
        # degree is the last index to process.
        sc.escape_test()

        m = nc_Out.Knots.KnotMultiplicity(iK)

        if m <= iMinMultiplicity:
            iK -= m
            continue

        #sEval = 'iK'; print(sEval+':',eval(sEval))
        #sEval = 'm'; print(sEval+':',eval(sEval))

        if fDev_Tol is None:
            bRemoved = nc_Out.Knots.RemoveKnots(iK, iK+1)

            if not bRemoved:
                print("RemoveKnots fail.")
                nc_Out.Dispose()
                return

            iK -= 1
            continue # while loop.


        nc_Out_Save = nc_Out.DuplicateCurve()

        bRemoved = nc_Out.Knots.RemoveKnots(iK, iK+1)

        if not bRemoved:
            print("RemoveKnots fail.")
            nc_Out_Save.Dispose()
            nc_Out.Dispose()
            return

        fDev = getDistancesBetweenCurves(nc_In, nc_Out)

        #sEval = 'fDev'; print(sEval+':',eval(sEval))

        if fDev <= fDev_Tol:
            iK -= 1
            nc_Out_Save.Dispose()
            continue # while loop.

        nc_Out.Dispose()
        nc_Out = nc_Out_Save.DuplicateCurve()

        iK -= 1

    if nc_Out.Knots.Count < nc_In.Knots.Count:
        #sc.doc.Objects.AddCurve(nc_Out); sc.doc.Views.Redraw(); 1/0
        return nc_Out

    nc_Out.Dispose()


def processCurve(rgCrv, iMinTargetGContinuity=2, iPreserveEnd=True, fDev_Tol=None, bRemoveAsManyAsPossibleWithinTol=False, bEcho=True, bDebug=False):
    """
    Parameters:
        rgCrv
        iMinTargetGContinuity
        iPreserveEnd,
        fDevTol:
            None for no deviation limit,
            Positive value for limit,
            Negative value for removing as many knots as possible within absolute value.
        bDebug
    """

    rc = filterTargetInput(rgCrv, iMinTargetGContinuity if fDev_Tol is None else None)
    if rc[0] is None:
        return rc

    nc_In = rc[0]

    degree = nc_In.Degree

    if iMinTargetGContinuity > 0:
        iMinTargetG_Actual = iMinTargetGContinuity
        iMinMult_ActualTarget = max((0, (degree - iMinTargetGContinuity)))
    else:
        iMinTargetG_Actual = degree
        iMinMult_ActualTarget = 0


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


    def overTargetKnotMultiplicityCount(nc):
        ct = 0
        iK = nc.Degree
        while iK < (nc.Knots.Count - nc.Degree):
            m = nc.Knots.KnotMultiplicity(iK)
            if m > iMinMult_ActualTarget:
                ct += 1
            iK += m
        return ct


    def processResult(nc_Ret, sFunc):
        if nc_Ret is None: return

        preserveEndCondition(nc_Ret)

        fDev = getDistancesBetweenCurves(nc_In, nc_Ret)
        if fDev_Tol is not None and fDev > fDev_Tol:
            if bDebug: print("Result with deviation {} skipped.".format(fDev))
            nc_Ret.Dispose()
            return

        ncs_Ret.append(nc_Ret)
        fDevs.append(fDev)
        iCts_Knots.append(nc_Ret.Knots)
        iCts_OverMs.append(overTargetKnotMultiplicityCount(nc_Ret))

        if bDebug:
            print("{}".format(sFunc),
                  "  Dev:{}".format(formatDistance(fDev)),
                  "  KnotCt:{}".format(nc_Ret.Knots.Count),
                  "  OverMultCt:{}".format(iCts_OverMs[-1]),
                  "  Multies:({})".format(
                      ",".join(str(i) for i in knotMultiplicityList(nc_Ret.Knots)))
                  )


    def printKnotMultiesChange(nc_In, nc_Out):
        s = "KnotMultiplicities:({})->({})".format(
            ",".join(str(i) for i in knotMultiplicityList(nc_In.Knots)),
            ",".join(str(i) for i in knotMultiplicityList(nc_Out.Knots)))
        print(s)



    if fDev_Tol:
        ncs_Ret = []
        fDevs = []
        iCts_Knots = []
        iCts_OverMs = []

        rc = removeKnots_All_possible(nc_In)
        processResult(rc, 'removeKnots_All_possible')
        if len(ncs_Ret):
            print("Solution acquired by removing all possible knots.")
            return (ncs_Ret[0], fDevs[0]), None



    if fDev_Tol is None:
        iMinTargetG_WIP = iMinTargetG_Actual
    else:
        iMinTargetG_WIP = degree



    while iMinTargetG_WIP >= iMinTargetG_Actual:

        ncs_Ret = []
        fDevs = []
        iCts_Knots = []
        iCts_OverMs = []


        if 1:
            rc = removeKnots_Plural_StartToEnd(nc_In, iMinTargetG_WIP, fDev_Tol=None)
            processResult(rc, 'removeKnots_Plural_StartToEnd')


        if 1:
            rc = removeKnots_Plural_EndToStart(nc_In, iMinTargetG_WIP, fDev_Tol=None)
            processResult(rc, 'removeKnots_Plural_EndToStart')


        if 1:
            rc = removeKnots_Single_StartToEnd(nc_In, iMinTargetG_WIP, fDev_Tol=None)
            processResult(rc, 'removeKnots_Single_StartToEnd')


        if 1:
            rc = removeKnots_Single_EndToStart(nc_In, iMinTargetG_WIP, fDev_Tol=None)
            processResult(rc, 'removeKnots_Single_StartToEnd')


        if len(ncs_Ret) == 0:
            if bDebug: print("No knots were removed by any routine.")
        else:
            fDev_Min = min(fDevs)

            idx_Dev_Min = fDevs.index(min(fDevs))
            if bDebug: print("WinningIndex:{}".format(idx_Dev_Min))

            for i, nc_Ret in enumerate(ncs_Ret):
                if i == idx_Dev_Min:
                    nc_Out = nc_Ret
                    continue # for loop.

                nc_Ret.Dispose()


            if fDev_Tol is None:
                if bDebug:
                    printKnotMultiesChange(nc_In, nc_Out)
                return (nc_Out, fDev_Min), None

            if fDev_Min <= fDev_Tol:
                if bDebug:
                    printKnotMultiesChange(nc_In, nc_Out)
                return (nc_Out, fDev_Min), None



        try: nc_Out.Dispose()
        except: pass


        if not bRemoveAsManyAsPossibleWithinTol:
            iMinTargetG_WIP -= 1
            continue # while loop.



        if bDebug:
            print("Removing all target multiplicities failed."
                  "  Now will try removing as many knots as possible within tolerance.")

        ncs_Ret = []
        fDevs = []
        iCts_Knots = []
        iCts_OverMs = []


        if 1:
            rc = removeKnots_Plural_StartToEnd(nc_In, iMinTargetG_WIP, fDev_Tol=fDev_Tol)
            processResult(rc, 'removeKnots_Plural_StartToEnd')


        if 1:
            rc = removeKnots_Plural_EndToStart(nc_In, iMinTargetG_WIP, fDev_Tol=fDev_Tol)
            processResult(rc, 'removeKnots_Plural_EndToStart')


        if 1:
            rc = removeKnots_Single_StartToEnd(nc_In, iMinTargetG_WIP, fDev_Tol=fDev_Tol)
            processResult(rc, 'removeKnots_Single_StartToEnd')


        if 1:
            rc = removeKnots_Single_EndToStart(nc_In, iMinTargetG_WIP, fDev_Tol=fDev_Tol)
            processResult(rc, 'removeKnots_Single_StartToEnd')


        if len(ncs_Ret) == 0:
            if bDebug: print("No knots were removed by any routine.")
        else:
            iCt_OverM_Min = min(iCts_OverMs)
            idxs_OverM_Min = [i for i in range(len(iCts_OverMs)) if iCts_OverMs[i] == iCt_OverM_Min]

            if len(idxs_OverM_Min) == 1 and iCt_OverM_Min == 0:
                # Single, true winner.

                idx_Winner = idxs_OverM_Min[0]

                if bDebug: print("WinningIndex:{}".format(idx_Winner))

                for i, nc_Ret in enumerate(ncs_Ret):
                    if i == idx_Winner:
                        nc_Out = nc_Ret
                        continue # for loop.
                    nc_Ret.Dispose()

                if bEcho:
                    printKnotMultiesChange(nc_In, nc_Out)

                return (nc_Out, fDevs[idx_Winner]), None

            elif iMinTargetG_WIP == iMinTargetG_Actual:
                # Last iteration, so take best.

                idx_Winner = idxs_OverM_Min[0]
                fDev_Min_In_Min = fDevs[idx_Winner]
                for i in range(1, len(idxs_OverM_Min)):
                    idx = idxs_OverM_Min[i]
                    if fDevs[idxs_OverM_Min[i]] < fDev_Min_In_Min:
                        idx_Winner = idx
                        fDev_Min_In_Min = fDevs[idx_Winner]

                if bDebug: print("WinningIndex:{}".format(idx_Winner))

                for i, nc_Ret in enumerate(ncs_Ret):
                    if i == idx_Winner:
                        nc_Out = nc_Ret
                        continue # for loop.
                    nc_Ret.Dispose()

                if bEcho:
                    printKnotMultiesChange(nc_In, nc_Out)

                return (nc_Out, fDevs[idx_Winner]), None

            else:
                # Try again.
                for nc in ncs_Ret: nc.Dispose()


        iMinTargetG_WIP -= 1


    return None, "Curve could not be modified within tolerance."


def processCurveObjects(objrefs_In, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iMinTargetGContinuity = getOpt('iMinTargetGContinuity')
    iPreserveEnd = getOpt('iPreserveEnd')
    fDevTol = getOpt('fDevTol') if getOpt('bLimitDev') else None
    bRemoveAsManyAsPossibleWithinTol = getOpt('bRemoveAsManyAsPossibleWithinTol')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    gCrvs_Out = []
    fDevs = []
    sLogs = []

    for objref_In in objrefs_In:
        rdCrv = objref_In.Object()
        rgCrv_In = rdCrv.CurveGeometry

        rc = processCurve(
            rgCrv=rgCrv_In,
            iMinTargetGContinuity=iMinTargetGContinuity,
            iPreserveEnd=iPreserveEnd,
            fDev_Tol=fDevTol,
            bRemoveAsManyAsPossibleWithinTol=bRemoveAsManyAsPossibleWithinTol,
            bEcho=bEcho,
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


    iMinTargetGContinuity = Opts.values['iMinTargetGContinuity']
    iPreserveEnd = Opts.values['iPreserveEnd']
    fDevTol = Opts.values['fDevTol'] if Opts.values['bLimitDev'] else None
    bRemoveAsManyAsPossibleWithinTol = Opts.values['bRemoveAsManyAsPossibleWithinTol']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObjects(
        objrefs_In=objrefs_In,
        iMinTargetGContinuity=iMinTargetGContinuity,
        iPreserveEnd=iPreserveEnd,
        fDevTol=fDevTol,
        bRemoveAsManyAsPossibleWithinTol=bRemoveAsManyAsPossibleWithinTol,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if gC_Res is None:
        return

    if gC_Res:
        sc.doc.Objects.UnselectAll()
        sc.doc.Objects.Select(objectIds=gC_Res)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()