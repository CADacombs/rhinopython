"""
This script uses the UI's _RebuildCrvNonUniform since an equivalent method is
not available in RhinoCommon as least up through 8.2.

_RebuildCrvNonUniform only returns degree-3 NURBS curves.

Reasons why this script was created
1.  _RebuildCrvNonUniform complies with the MaxPointCount option at the expense 
    of the RequestedTolerance option.  This script outputs curves only within the 
    requested tolerance / curve deviation.
2.  _RebuildCrvNonUniform works best on one curve at a time since it applies the 
    same number of control points to all input curves when multiple curves are 
    selected.  This script applies _RebuildCrvNonUniform to one curve at a time.
3.  Deviation reported by _RebuildCrvNonUniform has been found to be erroneous
    as compared to _CrvDeviation.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
170807-08: Created.
...
240105: Refactored. Now, skips uniform curve input. Changed an option.
        Added span count to printed report.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iMaxPointCount'; keys.append(key)
    values[key] = 103
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=5)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iPolyCrv_Opp'; keys.append(key)
    listValues[key] = 'No', 'InWhole' #, 'PerCrv' # All items must be strings.
    values[key] = 1
    names[key] = 'PolyCrv'
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
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Why is key, {}, here?  Value was not set or sticky-saved.".format(key))
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve wires with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=True)
    #go.EnableHighlight(False)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('fDevTol')
        addOption('iMaxPointCount')
        addOption('iPolyCrv_Opp')
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
        return "(None)"
    elif fDistance == 0.0:
        return "exactly 0".format(fDistance)
    elif fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-2)):
        return "{:.1e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def processNurbsCurve(nc, fDevTol, iMaxPointCount, bDebug=False):
    """
    Returns on success: rg.NurbsCurve, float(deviation)
    Returns on fail: None, str(log)
    """
    
    nc_In = nc
    
    if not isinstance(nc_In, rg.NurbsCurve):
        return None, "{} not processed.".format(nc_In.GetType().Name)

    if nc_In.SpanCount == 1:
        return None, "Bezier curve skipped.".format(nc_In.GetType().Name)

    if nc_In.Knots.KnotStyle in (rg.KnotStyle.Uniform, rg.KnotStyle.QuasiUniform):
        return None, "Uniform NURBS curve skipped.".format(nc_In.GetType().Name)

    fTol_WIP = fDevTol

    fDev_Min = None

    while fTol_WIP >= 1e-6:
        sc.escape_test()

        if bDebug:
            print('-'*80)
            sEval = "fTol_WIP"; print("{}: {}".format(sEval, eval(sEval)))

        gC_Out = sc.doc.Objects.AddCurve(nc_In)
        if gC_Out == gC_Out.Empty:
            return None, "Could not add starting curve to document."

        sc.doc.Objects.UnselectAll()
    
        if not sc.doc.Objects.Select(objectId=gC_Out, select=True):
            raise ValueError("Could not select WIP curve.")
    
        sCmd  = "_-RebuildCrvNonUniform "
        sCmd += "_RequestedTolerance={} ".format(fTol_WIP)
        sCmd += "_MaxPointCount={} ".format(iMaxPointCount)
        sCmd += "_Quarters=No "
        sCmd += "_DeleteInput=Yes "
        sCmd += "_EnterEnd"

        sc.doc.UndoRecordingEnabled = False

        Rhino.RhinoApp.RunScript(script=sCmd, echo=bDebug)

        rdC_Out = sc.doc.Objects.FindId(gC_Out)
        nc_Out = rdC_Out.Geometry

        sc.doc.Objects.Delete(objectId=gC_Out, quiet=False)

        sc.doc.UndoRecordingEnabled = True

        bEqual = nc_Out.EpsilonEquals(other=nc_In, epsilon=Rhino.RhinoMath.ZeroTolerance)

        if bEqual:
            return None, "Input and output are identical."

        rc = rg.Curve.GetDistancesBetweenCurves(
            nc_In, nc_Out, 0.1*min((fDevTol, sc.doc.ModelAbsoluteTolerance)))

        if bDebug:
            sPrint = 'nc_Out.Points.Count'; print(sPrint + ':', eval(sPrint))
            sPrint = 'nc_Out.Degree'; print(sPrint + ':', eval(sPrint))
    
        if not rc[0]:
            return None, "Curve deviation could not be determined."

        fDev = rc[1]

        if bDebug:
            sPrint = 'fDev'; print(sPrint + ':', eval(sPrint))

        if fDev <= fDevTol:
            return nc_Out, fDev

        if fDev_Min is None:
            fDev_Min = fDev
        else:
            if fDev > fDev_Min:
                if bDebug: print("Deviation started to increase.")
                break

        fTol_WIP /= 2.0


    return None, "Curve could not be rebuilt within tolerance."



def processCurveObjects(objrefs_In, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    iMaxPointCount = getOpt('iMaxPointCount')
    iPolyCrv_Opp = getOpt('iPolyCrv_Opp')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    gCs_Out = []
    sLogs = []


    def isKnotVectorUniform(knots):
        return (
            (knots.KnotStyle == rg.KnotStyle.Uniform) or
            (knots.KnotStyle == rg.KnotStyle.QuasiUniform) or
            (
                (knots.KnotStyle == rg.KnotStyle.PiecewiseBezier) and
                knots.Count == knots.KnotMultiplicity(0) * 2)
            )


    for objref_In in objrefs_In:
        rdObj_In = objref_In.Object()

        rgC_In = objref_In.Curve()
        if rgC_In is None:
            sLogs.append("{} skipped.".format(rdObj_In.ObjectType.Name))
            continue # to next object.

        if isinstance(rgC_In, rg.PolyCurve):
            if iPolyCrv_Opp == 0:
                sLogs.append("PolyCurve skipped.")
                continue # to next curve.
            elif iPolyCrv_Opp == 1:
                nc_In = rgC_In.ToNurbsCurve()
            else:
                # TODO:
                sLogs.append("Processing subcurves of PolyCurves is not yet supported.")
                continue # to next curve.
        else:
            nc_In = rgC_In.ToNurbsCurve() # Regardless if rgC_In is already a NurbsCurve since it needs to be duplicated anyway.

        rc = processNurbsCurve(
            nc=nc_In,
            fDevTol=fDevTol,
            iMaxPointCount=iMaxPointCount,
            bDebug=bDebug)

        if rc[0] is None:
            sLogs.append(rc[1])
            continue

        nc_Ret, fDev = rc

        if bReplace and rdObj_In.ObjectType == rd.ObjectType.Curve:
            if not sc.doc.Objects.Replace(rdObj_In.Id, nc_Ret):
                raise ValueError("Curve could not be replaced!")
            gCs_Out.append(rdObj_In.Id)
        else:
            # This includes objrefs of Edges.
            gC_Out = sc.doc.Objects.AddCurve(nc_Ret)
            if gC_Out == gC_Out.Empty:
                raise ValueError("Curve could not be added!")
            gCs_Out.append(gC_Out)
        
        if len(objrefs_In) == 1:
            s = "Prop:I->O"
            s += "  {}:{}->{}".format("PtCt", nc_In.Points.Count, nc_Ret.Points.Count)
            s += "  {}:{}->{}".format("SpanCt", nc_In.SpanCount, nc_Ret.SpanCount)
            if Rhino.RhinoApp.ExeVersion >= 7:
                s += "  {}:{}->{}".format("IsUniform",
                        str(isKnotVectorUniform(nc_In.Knots))[0],
                        str(isKnotVectorUniform(nc_Ret.Knots))[0])
            s += "  {}:{}->{}".format("IsPeriodic",
                    str(nc_In.IsPeriodic)[0],
                    str(nc_Ret.IsPeriodic)[0])
            s += "\nDeviation: {}".format(formatDistance(fDev))

            print(s)

        nc_In.Dispose()
        nc_Ret.Dispose()

    if sLogs and len(objrefs_In) == 1:
        print(sLogs[0])
    else:
        for sLog in set(sLogs):
            print("[{}] {}".format(sLogs.count(sLog), sLog))
    
    if len(objrefs_In) > 1:
        print("{} curves were successfully processed.".format(len(gCs_Out)))

    return gCs_Out


def main():

    objrefs_In = getInput()
    if objrefs_In is None: return

    fDevTol = Opts.values['fDevTol']
    iMaxPointCount = Opts.values['iMaxPointCount']
    iPolyCrv_Opp = Opts.values['iPolyCrv_Opp']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug:
        sc.doc.Views.RedrawEnabled = False

    processCurveObjects(
        objrefs_In,
        fDevTol=fDevTol,
        iMaxPointCount=iMaxPointCount,
        iPolyCrv_Opp=iPolyCrv_Opp,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()