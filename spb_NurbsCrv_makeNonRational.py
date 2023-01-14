"""
This script will set all weights of a rational NURBS curve to 1.0 and
IsRational (RhinoCommon) to False.  Option is included for allowed curve deviation.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
230114: Created.
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


    key = 'bLimitDistDev'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bModifyArcCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bModifyPolyCrvs'; keys.append(key)
    values[key] = False
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

        if key == 'fDistTol':
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

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    bPreselectedObjsChecked = False

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bLimitDistDev')
        if Opts.values['bLimitDistDev']:
            addOption('fDistTol')
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

        if res == ri.GetResult.Number and Opts.values['bLimitDistDev']:
            key = 'fDistTol'
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


def hasRationalWeights(nc):
    for cp in nc.Points:
        weight = cp.Weight
        if weight == 1.0: continue
        return True
    return False


def getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
        crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)
    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None
    return rc[1]


def createNonRationalNurbsCurve(nc_In, fDistTol=None, bDebug=False):
    """
    Parameters:
        fDistTol: None for no limit.
        bDebug

    Returns:
        Success: rgNurbsCurve, fMinRadius, fDeviation, None
        Fail: None, None, None, sLog
    """
    
    if not isinstance(nc_In, rg.NurbsCurve):
        return None, "Not a NURBS curve."


    if not nc_In.IsRational:
        if hasRationalWeights(nc_In):
            return None, "Curve is flagged as being non-rational, but some of its weights are not 1.0!"
        else:
            return None, "Curve is flagged as being non-rational, and all its weights are 1.0."

    nc_Out = rg.NurbsCurve(degree=nc_In.Degree, pointCount=nc_In.Points.Count)

    for iK, k in enumerate(nc_In.Knots):
        nc_Out.Knots[iK] = k

    for iCp, cp in enumerate(nc_In.Points):
        nc_Out.Points[iCp] = nc_In.Points[iCp].Location


    if fDistTol is not None:
        dev = getDistancesBetweenCurves(nc_In, nc_Out)
        if dev > fDistTol:
            nc_Out.Dispose()
            return None, "Result is not within allowed distance deviation of {}".format(fDistTol)

    return nc_Out, None


def processCurveObject(rhCrv_In, **kwargs):
    """
    Parameters:
        fDistTol
        bModifyArcCrvs
        bModifyPolyCrvs
        bModify_NotAdd
        bEcho
        bDebug
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDistTol = getOpt('fDistTol')
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

    if not nc_Start.IsRational:
        if hasRationalWeights(nc_Start):
            nc_Start.Dispose()
            return None, "Curve is flagged as being non-rational, but some of its weights are not 1.0!"
        else:
            nc_Start.Dispose()
            return None, "Curve is flagged as being non-rational, and all its weights are 1.0."


    gCrv_In = rdCrv_In.Id

    rc = createNonRationalNurbsCurve(
        nc_In=nc_Start,
        fDistTol=fDistTol,
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

    fDistTol = Opts.values['fDistTol'] if Opts.values['bLimitDistDev'] else None
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
            fDistTol=fDistTol,
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