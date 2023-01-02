"""
This script is an alternative to _MakeUniform for curves in that it includes
options to
    1. Reject results outside of a distance tolerance,
    2. Preserve end curvatures,
    3. Filter some curve types,
    4. Add instead of only modify input.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
200507: Created.
221122: Import-related update.
230102: Added option to preserve curvature of each ends.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bLimitCrvDev'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPreserveEndCurvatures'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bModifyArcCrvs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bModifyPolyCrvs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bModify_NotAdd'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Modify')
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

        if key == 'fDevTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve and parameter with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    #go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    idxs_Opt = {}

    bPreselectedObjsChecked = False

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bLimitCrvDev')
        if Opts.values['bLimitCrvDev']:
            addOption('fDevTol')
            go.AcceptNumber(True, acceptZero=True)
        else:
            go.AcceptNumber(False, acceptZero=True)
        addOption('bPreserveEndCurvatures')
        addOption('bModifyArcCrvs')
        addOption('bModifyPolyCrvs')
        addOption('bModify_NotAdd')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fDevTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def formatDistance(fDistance):
    if fDistance is None: 
        return "(No deviation provided)"
    if fDistance == 0.0:
        return 0.0
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def areEpsilonEqual(a, b, epsilon):
    # This is a relative comparison.
    delta = abs(a-b)
    #print f2s(delta),
    fRelComp = delta / max(abs(a), abs(b))
    #print f2s(fRelComp)
    return fRelComp < epsilon


def isUniform(nc):
    if nc.Points.Count == nc.Degree + 1:
        return True

    start = 0 if nc.IsPeriodic else nc.Degree - 1
    end = nc.Knots.Count - (0 if nc.IsPeriodic else nc.Degree - 1) - 1

    #print start, end

    span0 = nc.Knots[start+1] - nc.Knots[start]

    #print f2s(span0)

    for i in range(start+1, end):
        if nc.Knots.KnotMultiplicity(i) > 1:
            return False
        #print f2s(nc.Knots[i+1] - nc.Knots[i])
        if not areEpsilonEqual(
            span0, nc.Knots[i+1] - nc.Knots[i],
            epsilon=1e-9):
                return False

    return True


def matchEndCurvatures(nc_ToMod, nc_Ref):
    """
    nc_ToMod is modified.
    """
    bSuccess = nc_ToMod.SetEndCondition(
        bSetEnd=False,
        continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
        point=nc_ToMod.PointAtStart,
        tangent=nc_Ref.TangentAtStart,
        curvature=nc_Ref.CurvatureAt(nc_Ref.Domain.T0))
    if not bSuccess:
        print("SetEndCondition failed.")
    bSuccess = nc_ToMod.SetEndCondition(
        bSetEnd=True,
        continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
        point=nc_ToMod.PointAtEnd,
        tangent=nc_Ref.TangentAtEnd,
        curvature=nc_Ref.CurvatureAt(nc_Ref.Domain.T1))
    if not bSuccess:
        print("SetEndCondition failed.")

    return True


def getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None

    return rc[1]


def createUniformNurbsCurve(nc_In, fDevTol=None, bPreserveEndCurvatures=False, bDebug=False):
    """
    Parameters:
        fDevTol: None for no limit.
        bPreserveEndCurvatures
        bDebug

    Returns:
        Success: rgNurbsCurve, fMinRadius, fDeviation, None
        Fail: None, None, None, sLog
    """
    
    if not isinstance(nc_In, rg.NurbsCurve):
        return None, "Not a NURBS curve."


    if isUniform(nc_In):
        return None, 'Curve is already uniform.'

    nc_Out = nc_In.ToNurbsCurve()


    if nc_In.IsPeriodic:

        nc_Out.Knots.CreatePeriodicKnots(knotSpacing=1.0) # Modifies Knots.

        #fK = float(-nc_In.Degree)

        #for iK in range(nc_In.Knots.Count):
        #    fK += 1.0

        #    nc_Out.Knots[iK] = fK
    else:

        nc_Out.Knots.CreateUniformKnots(knotSpacing=1.0) # Modifies Knots.

        #fK = 0.0

        #for iK in range(nc_In.Knots.Count):
        #    if nc_In.Degree <= iK <= (nc_In.Knots.Count - nc_In.Degree):
        #        fK += 1.0

        #    nc_Out.Knots[iK] = fK

    if bPreserveEndCurvatures:
        if not matchEndCurvatures(nc_Out, nc_In):
            raise Exception("Preserve end curvatures failed.")

    if fDevTol is not None:
        dev = getDistancesBetweenCurves(nc_In, nc_Out)
        if dev > fDevTol:
            nc_Out.Dispose()
            return None, "Result is not within allowed deviation."

    return nc_Out, None


def processCurveObject(rhCrv_In, **kwargs):
    """
    Parameters:
        fDevTol
        bPreserveEndCurvatures
        bModifyArcCrvs
        bModifyPolyCrvs
        bModify_NotAdd
        bEcho
        bDebug
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bPreserveEndCurvatures = getOpt('bPreserveEndCurvatures')
    bModifyArcCrvs = getOpt('bModifyArcCrvs')
    bModifyPolyCrvs = getOpt('bModifyPolyCrvs')
    bModify_NotAdd = getOpt('bModify_NotAdd')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getCurve(rhObj):
        if isinstance(rhObj, rd.CurveObject):
            return rhObj, rhObj.CurveGeometry

        if isinstance(rhObj, rg.Curve):
            return None, rhObj

        if isinstance(rhObj, rg.GeometryBase):
            rdObj = None
            rgObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
            rgObj = rhObj.Geometry()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
            rgObj = rdObj.Geometry
        else:
            return

        if isinstance(rgObj, rg.Curve):
            return rdObj, rgObj


    rdCrv_In, rgCrv_In = getCurve(rhCrv_In)

    if isinstance(rgCrv_In, rg.NurbsCurve):
        nc_Start = rgCrv_In.DuplicateCurve()
    elif bModifyArcCrvs and isinstance(rgCrv_In, rg.ArcCurve):
        nc_Start = rgCrv_In.ToNurbsCurve()
    elif bModifyPolyCrvs and isinstance(rgCrv_In, rg.PolyCurve):
        nc_Start = rgCrv_In.ToNurbsCurve()
    else:
        return None, "{} is not an accepted input.".format(rgCrv_In.GetType().Name)

    if isUniform(nc_Start):
        nc_Start.Dispose()
        return None, 'Curve is already uniform.'

    gCrv_In = rdCrv_In.Id

    rc = createUniformNurbsCurve(
        nc_In=nc_Start,
        fDevTol=fDevTol,
        bPreserveEndCurvatures=bPreserveEndCurvatures,
        bDebug=bDebug,
        )
    if rc[0] is None:
        return rc

    nc_Res, sLog = rc

    if sLog:
        raise Exception("sLog should be None but is {}".format(sLog))

    nc_Res = rc[0]

    if bModify_NotAdd:
        if sc.doc.Objects.Replace(gCrv_In, nc_Res):
            return gCrv_In, "Curve was replaced."
        else:
            return None, "Curve could not be replaced."
    else:
        g1 = sc.doc.Objects.AddCurve(nc_Res)
        if g1 != Guid.Empty:
            return gCrv_In, "Curve was added."
        else:
            return None, "Curve could not be added."

    nc_Start.Dispose()


    return gs1, sLogs


def main():

    objrefs = getInput()
    if objrefs is None: return

    fDevTol = Opts.values['fDevTol'] if Opts.values['bLimitCrvDev'] else None
    bPreserveEndCurvatures = Opts.values['bPreserveEndCurvatures']
    bModifyArcCrvs = Opts.values['bModifyArcCrvs']
    bModifyPolyCrvs = Opts.values['bModifyPolyCrvs']
    bModify_NotAdd = Opts.values['bModify_NotAdd']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if not bDebug:
        sc.doc.Views.RedrawEnabled = False

    gCs_Result = []
    sLogs = []

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    for objref in objrefs:

        gC_Result, sLog = processCurveObject(
            objref,
            fDevTol=Opts.values['fDevTol'] if Opts.values['bLimitCrvDev'] else None,
            bPreserveEndCurvatures=bPreserveEndCurvatures,
            bModifyArcCrvs=bModifyArcCrvs,
            bModifyPolyCrvs=bModifyPolyCrvs,
            bModify_NotAdd=bModify_NotAdd,
            bEcho=bEcho,
            bDebug=bDebug
            )

        if gC_Result:
            gCs_Result.append(gC_Result)

        if sLog:
            sLogs.append(sLog)

    if bEcho:
        for sLog in set(sLogs):
            print("[{}] {}".format(sLogs.count(sLog), sLog))


    for gC in gCs_Result:
        sc.doc.Objects.Select(gC)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()