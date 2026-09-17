#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
Send any questions, comments, or script development service needs to
@spb on the McNeel Forums ( https://discourse.mcneel.com/ ).
"""

"""
190417: This script started as a split from another.
191019-200121: Import-related updates.
200210: Bug fix.  Import-related update.
210316: Bug fix for BrepEdge support.
220328: Import-related update.
221221-22: Import-related update.  Removed an option.  Refactored.
260915-16: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import spb_RebuildCrvUniform


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
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=1)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iMaxCpCt'; keys.append(key)
    values[key] = 40
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=2)
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

    s  = "For MaxCpCt, \"0\" will use the curve's current control point count."
    s += "  Other values are the absolute maximum to allow."
    print(s)


    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fDevTol')
        addOption('bPreserveEndG1')
        addOption('iDegree')
        addOption('iMaxCpCt')
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


def processCurves(rgCrvs_In, **kwargs):#fDevTol=None, iDegree=None, bPreserveEndG1=None, iMaxCpCt=None, bReplace=None, bEcho=None, bDebug=None):
    """
    rgCrvs_In = rg.Curve
    fDevTol
    iDegree
    bPreserv
    iMaxCpCt
    bEcho
    bDebug
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bPreserveEndG1 = getOpt('bPreserveEndG1')
    iDegree = getOpt('iDegree')
    iMaxCpCt = getOpt('iMaxCpCt')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    ncs1 = []
    devs1 = []

    for i, nc0 in enumerate(rgCrvs_In):
        nc1, dev1, sLog = spb_RebuildCrvUniform.findRebuild(
                rgCurve0=nc0,
                fDevTol=fDevTol,
                iDegree=iDegree,
                iPreserveEndG=1 if bPreserveEndG1 else 0,
                bFurtherTranslateCps=False,
                iMinCpCt=None,
                iMaxCpCt=iMaxCpCt,
                bDebug=bDebug,
        )
        if nc1 is None:
            if bEcho:
                print("Dev:{}  Log:{} for Curve[{}]".format(
                    dev1, sLog, i))
            return
        
        ncs1.append(nc1)
        devs1.append(dev1)
    
    iMax_cpCt_ncs1 = max([nc1.Points.Count for nc1 in ncs1])
    if iMax_cpCt_ncs1 > iMaxCpCt:
        s  = "Required control point count to create all curves within tolerance:"
        s += " {}".format(iMax_cpCt_ncs1)
        for c in ncs1: c.Dispose()
        return

    if min([nc1.Points.Count for nc1 in ncs1]) == iMax_cpCt_ncs1:
        # NurbsCurve are already at same control point count.
        ncs2 = [nc1.DuplicateCurve() for nc1 in ncs1]
        devs2 = devs1[:]
    else:
        pointCount = iMax_cpCt_ncs1
        while pointCount <= iMaxCpCt:
            sc.escape_test()
            ncs2 = []
            devs2 = []
            for nc0 in rgCrvs_In:
                nc2 = nc0.Rebuild(
                        pointCount=pointCount,
                        degree=iDegree,
                        preserveTangents=bPreserveEndG1)

                if nc2 is None:
                    raise ValueError("NurbsCurve.Rebuild returned None.")
            
                ncs2.append(nc2)

                rc = rg.Curve.GetDistancesBetweenCurves(nc0, nc2, 0.1*fDevTol)
                if rc[0]:
                    dev2 = rc[1]
                    if dev2 > fDevTol:
                        for c in ncs2: c.Dispose()
                        pointCount += 1
                        break # for loop back to while to try to Rebuild all curves to new pointCount.
                    devs2.append(dev2)
            else:
                # All curves were rebuilt within fDevTol.
                break

    return ncs2, devs2


def processCurveObjects(curvesAndEdges0, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bPreserveEndG1 = getOpt('bPreserveEndG1')
    iDegree = getOpt('iDegree')
    iMaxCpCt = getOpt('iMaxCpCt')
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
        iMaxCpCt=iMaxCpCt,
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
    iMaxCpCt = Opts.values['iMaxCpCt']
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
        iMaxCpCt=iMaxCpCt,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()