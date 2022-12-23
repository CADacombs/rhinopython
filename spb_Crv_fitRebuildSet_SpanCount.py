"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221221-22: Created, starting with another script.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fDevTol'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPreserveEndG1'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDegree'; keys.append(key)
    values[key] = 3
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=2)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iPowerOf2ForMaxKnotSpanCt'; keys.append(key)
    values[key] = 6
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

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
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-9:
                cls.values[key] = cls.riOpts[key].CurrentValue = 1e-9
            else:
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
    Get curves with optional input
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = rd.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fDevTol')
        addOption('bPreserveEndG1')
        addOption('iDegree')
        addOption('iPowerOf2ForMaxKnotSpanCt')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=2, maximumNumber=0)

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

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getMaximumDeviation(rgCrvA, rgCrvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            rgCrvA,
            rgCrvB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
    if rc[0]:
        return rc[1]


def rebuildFit(rgCrv_In, spanCount, degree, preserveTangents, tolerance):
    """
    Returns on success:
        rg.NurbsCurve,
        fDeviation
    Returns on fail:
        None
    """

    pointCount = degree + spanCount

    nc_Res = rgCrv_In.Rebuild(pointCount, degree, preserveTangents)

    if nc_Res is None: return

    dev = getMaximumDeviation(nc_Res, rgCrv_In)
    if dev is None or dev > tolerance:
        nc_Res.Dispose()
        return

    return nc_Res, dev


def processCurves(rgCrvs_In, **kwargs):
    """
    rgCrvs_In = rg.Curve
    fDevTol
    iDegree
    bPreserv
    iPowerOf2ForMaxKnotSpanCt
    bEcho
    bDebug
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bPreserveEndG1 = getOpt('bPreserveEndG1')
    iDegree = getOpt('iDegree')
    iPowerOf2ForMaxKnotSpanCt = getOpt('iPowerOf2ForMaxKnotSpanCt')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    ncs_Out = []
    devs_Out = []


    iPossibleKnotCts = [2**p for p in range(iPowerOf2ForMaxKnotSpanCt+1)]


    for i, nc0 in enumerate(rgCrvs_In):

        if isinstance(nc0, rg.NurbsCurve):
            if nc0.SpanCount == 1:
                ncs_Out.append(nc0)
                devs_Out.append(0.0)
                continue # to next curve.

        # Try single knot span count.
        rc = rebuildFit(nc0, 1, iDegree, bPreserveEndG1, fDevTol)
        if rc is not None:
            nc_Res, dev = rc
            ncs_Out.append(nc_Res)
            devs_Out.append(dev)
            continue # to next curve.

        # Try maximum knot span counts.  Failure results in failure for all.
        rc = rebuildFit(nc0, iPossibleKnotCts[-1], iDegree, bPreserveEndG1, fDevTol)
        if rc is None:
            print("Need more than {} knot spans for curve[{}].".format(
                iPossibleKnotCts[-1], i))
            return

        nc_Res, dev = rc


        # Binary search
        nc_Hi = nc_Res
        dev_Hi = dev

        # Po2 == Power of 2

        iPo2_SpanCt_Lo = 0
        iPo2_SpanCt_Hi = iPowerOf2ForMaxKnotSpanCt

        while True:
            sc.escape_test()

            iPo2_SpanCt_Md = (iPo2_SpanCt_Hi + iPo2_SpanCt_Lo) // 2

            if iPo2_SpanCt_Md in (iPo2_SpanCt_Lo, iPo2_SpanCt_Hi):
                ncs_Out.append(nc_Hi)
                devs_Out.append(dev_Hi)
                break # out of while loop.

            rc = rebuildFit(nc0, 2**iPo2_SpanCt_Md, iDegree, bPreserveEndG1, fDevTol)
            if rc is None:
                iPo2_SpanCt_Lo = iPo2_SpanCt_Md
                continue

            nc_Res, dev = rc

            nc_Hi.Dispose()
            nc_Hi = nc_Res
            dev_Hi = dev
            iPo2_SpanCt_Hi = iPo2_SpanCt_Md

    return ncs_Out, devs_Out


def processCurveObjects(curvesAndEdges0, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bPreserveEndG1 = getOpt('bPreserveEndG1')
    iDegree = getOpt('iDegree')
    iPowerOf2ForMaxKnotSpanCt = getOpt('iPowerOf2ForMaxKnotSpanCt')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    list_bIsWire = []
    gObjs_In = []
    rgCrvs_In = []


    for curveOrEdge0 in curvesAndEdges0:
        rdObj_In = rs.coercerhinoobject(curveOrEdge0)
        if isinstance(rdObj_In, rd.CurveObject):
            list_bIsWire.append(True)
        elif isinstance(rdObj_In, rd.BrepObject):
            list_bIsWire.append(False)
        else:
            raise ValueError("Invalid input: {}".format(rdObj_In))
        
        gObjs_In.append(rdObj_In.Id)
        
        rgCrvs_In.append(rs.coercecurve(curveOrEdge0))


    rc = processCurves(
        rgCrvs_In=rgCrvs_In,
        fDevTol=fDevTol,
        bPreserveEndG1=bPreserveEndG1,
        iDegree=iDegree,
        iPowerOf2ForMaxKnotSpanCt=iPowerOf2ForMaxKnotSpanCt,
        bEcho=bEcho,
        bDebug=bDebug,
        )


    if rc is None: return

    ncs_Res, devs_Res = rc


    ct_Replaced = 0
    ct_Added = 0

    for gObj_In, nc2, bIsWire in zip(gObjs_In, ncs_Res, list_bIsWire):
        if bReplace and bIsWire:
            if sc.doc.Objects.Replace(gObj_In, nc2):
                ct_Replaced += 1
        else:
            gCrv1 = sc.doc.Objects.AddCurve(nc2)
            if gCrv1 != gCrv1.Empty:
                ct_Added += 1

    s = ""

    if bReplace and any(list_bIsWire):
        if ct_Replaced == 0:
            s += "No wires were replaced.".format(ct_Replaced)
        elif ct_Replaced == len(gObjs_In):
            s += "All {} wires were replaced.".format(ct_Replaced)
        else:
            s +=  "Only {} out of {} wires were replaced.".format(ct_Replaced, len(gObjs_In))

    if ct_Added == len(gObjs_In):
        s += "All {} curves were added.".format(ct_Added)
    elif not all(list_bIsWire):
        s += "Only {} out of {} curves were added.".format(ct_Added, len(gObjs_In))

    s += "  Degree: {}".format(nc2.Degree)
    s += "  CpCt: {}".format(nc2.Points.Count)
    s += "  Maximum deviation: {0:.{1}f}".format(
            max(devs_Res), sc.doc.ModelDistanceDisplayPrecision)
    
    print(s)


def main():

    objrefs = getInput()
    if objrefs is None: return


    fDevTol = Opts.values['fDevTol']
    bPreserveEndG1 = Opts.values['bPreserveEndG1']
    iDegree = Opts.values['iDegree']
    iPowerOf2ForMaxKnotSpanCt = Opts.values['iPowerOf2ForMaxKnotSpanCt']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    processCurveObjects(
        curvesAndEdges0=objrefs,
        fDevTol=fDevTol,
        bPreserveEndG1=bPreserveEndG1,
        iDegree=iDegree,
        iPowerOf2ForMaxKnotSpanCt=iPowerOf2ForMaxKnotSpanCt,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()