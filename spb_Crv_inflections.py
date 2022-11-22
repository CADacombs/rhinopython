"""
190703: Created.
190707: Now, inflection points found within ModelAbsoluteTolerance of ends of curves are skipped.
190708: Curves are now processed only as NurbsCurves.
191010: Now, preselected curves are immediately processed.
191203: Modified number input in UI and modified printed output.
191208: Modified output of getInflectionParameters for ArcCurves and LineCurves from None to [].
191219: Bug fix that caused some curves to be deleted.  Modified some return value options of processCurve.
200114: Replaced main routine to find maximums with code from McNeel.  (Search for 'github' for link.)
        Same/similar code is supposed to be available in V7's RhinoCommon.
        Refactored.
200419: Bug fix.

TODO:
Add code for and enable bIncludeEndsOfLinearSpans?  Add a linear tolerance?
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System.Collections.Generic import List
from System.Drawing import Color
from System import Guid
from System import Math


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'bIncludeEndsOfLinearSpans'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddPts'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSplitCrv'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDeleteInput'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """Get edges with optional input."""
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select curves")
    go.GeometryFilter = rd.ObjectType.Curve # Curve is also used for brep edges.
    
    go.AcceptNumber(True, True)
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    idxs_Opts = {}

    while True:
        #key = 'bIncludeEndsOfLinearSpans'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bAddPts'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bSplitCrv'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        if Opts.values['bSplitCrv']:
            key = 'bDeleteInput'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInflectionParameters(curve):
    """
    Returns:
        list(float(Parameters)) on success.
        None on failure.
    This is a translation and modification of code at
    https://github.com/mcneel/rhino-developer-samples/blob/6/rhinocommon/cs/SampleCsCommands/SampleCsExtractInflectionPoints.cs
    """



    def findInflection(nurb, t0, t1, k0, k1, epsilon):
        """
        Returns:
            Success: float(Parameter)
            Failure: None
        """

        if not nurb: return

        t_Out = None

        while True:
            sc.escape_test()

            t = (t0 + t1) * 0.5
            k = nurb.CurvatureAt(t)
            if not k.IsValid:
                break

            if k.IsTiny() or t1 - t0 < epsilon:
                t_Out = t
                break

            if k * k0 < 0.0:
                t1 = t
                k1 = k
                continue

            if k * k1 < 0.0:
                t0 = t
                k0 = k
                continue

            break

        return t_Out



    on_zero_tolerance = 2.3283064365386962890625e-10 # Why is this used instead of 

    nurb = curve.ToNurbsCurve()
    if not nurb: return

    count = nurb.Points.Count * 4
    mul = 1.0 / count
    epsilon = (on_zero_tolerance if nurb.Domain.Length > 1.0
               else nurb.Domain.Length * on_zero_tolerance)

    t0 = 0.0
    k0 = rg.Vector3d.Unset
    start_set = False
    ts = []

    for i in xrange(count+1):
        t1 = nurb.Domain.ParameterAt(i * mul)
        k1 = nurb.CurvatureAt(t1)
        if k1.IsValid:
            if k1.IsTiny(): continue

            if not start_set:
                t0 = t1
                k0 = k1
                start_set = True
                continue

            if k0 * k1 < 0.0:
                t = findInflection(nurb, t0, t1, k0, k1, epsilon)
                if t is not None: ts.Add(t)

            k0 = k1
            t0 = t1

    return ts


def getInflectionPoints(curve):

    ts = getInflectionParameters(curve)
    if not ts: return ts

    return [curve.PointAt(t) for t in ts]


def processCurves(curvesAndEdges0, **kwargs):
    """
    """



    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bIncludeEndsOfLinearSpans = getOpt('bIncludeEndsOfLinearSpans')
    bAddPts = getOpt('bAddPts')
    bSplitCrv = getOpt('bSplitCrv')
    bDeleteInput = getOpt('bDeleteInput')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')



    def splitCurve():
        """
        Returns:
            Success: list(GUIDs of curves from split), None
            Fail: None, str(Failure)
        """
        gCrvs_Split = []
        rgNurbsCrvs2_Split = rgCrv0.Split(ts_Inflection)
        rgCrv0.Dispose()
        if rgNurbsCrvs2_Split is None:
            return None, "Curve could not be split."
        if len(rgNurbsCrvs2_Split) == 1:
            return None, "Curve was split into just one curve.  Inflection point must be too close to an endpoint."
        
        attr0 = rdCrv0.Attributes if bDeleteInput else None
        bFail = False
        for rgCrv_Split in rgNurbsCrvs2_Split:
            gCrv_Split = sc.doc.Objects.AddCurve(rgCrv_Split, attr0)
            rgCrv_Split.Dispose()
            if gCrv_Split != Guid.Empty:
                gCrvs_Split.append(gCrv_Split)
            else:
                bFail = True
        if not bFail:
            sc.doc.Objects.Delete(gCrv0, quiet=False)
            if bEcho: print "Curve was successfully split into {}.".format(len(rgNurbsCrvs2_Split))
        else:
            if bDeleteInput:
                s  = "Curve was not deleted because"
                s += "not all split curves could be added."
                print s
            else:
                print "Not all split curves could be added to the document."
        for rgCrv_Split in rgNurbsCrvs2_Split: rgCrv_Split.Dispose()
        return gCrvs_Split, None


    gCrvs0 = []
    for curveOrEdge0 in curvesAndEdges0:
        gCrv0 = rs.coerceguid(curveOrEdge0)
        if gCrv0:
            gCrvs0.append(gCrv0)
    
    gNcs_Split_perCrv0 = []
    pts_Inflection_All = []
    gPts_Inflection_Nested = []
    
    sFails = []
    
    idxs_AtTenths = [int(round(0.1*i*len(curvesAndEdges0),0)) for i in range(10)]

    bEcho_inLoop = bEcho if len(curvesAndEdges0) == 1 else False
    
    for iC, curveOrEdge0 in enumerate(curvesAndEdges0):
        if iC in idxs_AtTenths:
            Rhino.RhinoApp.SetCommandPrompt(
                    "Processing curve {} ...".format(
                    "" if len(curvesAndEdges0) == 1 else "{} of {} ".format(
                        iC+1, len(curvesAndEdges0))))

        rgCrv0 = rs.coercecurve(curveOrEdge0) # Will return various rg.Curves, including rg.BrepEdge.
        if rgCrv0 is None:
            sFails.append("Geometry not found for {}.".format(curveOrEdge0))
            continue

        if isinstance(rgCrv0, rg.BrepEdge):
            bDeleteInput = False
        else:
            gCrv0 = rs.coerceguid(curveOrEdge0)
            rdCrv0 = rs.coercerhinoobject(curveOrEdge0)
    
        ts_Inflection = getInflectionParameters(rgCrv0)

        if ts_Inflection is None:
            s  = "getInflectionParameters returned None"
            sFails.append(s)
            s += " for {}.".format(curveOrEdge0)
            s += "  Bad input Curve?"
            if bEcho_inLoop: print s
            continue
    
        if bEcho_inLoop:
            print "{} inflection(s) found.".format(len(ts_Inflection))

        if not ts_Inflection: continue

        pts_Inflection = [rgCrv0.PointAt(t) for t in ts_Inflection]

        pts_Inflection_All.extend(pts_Inflection)

        if bAddPts:
            gPts = []
            for pt in pts_Inflection:
                gPt = sc.doc.Objects.AddPoint(pt)
                if gPt != Guid.Empty: gPts.append(gPt)
            gPts_Inflection_Nested.append(gPts)

        if bSplitCrv:
            gNcs_Split, sFail = splitCurve()
            if gNcs_Split: gNcs_Split_perCrv0.append(gNcs_Split)
            if sFail: sFails.append(sFail)
            gNcs_Split_perCrv0.append(gNcs_Split)

        rgCrv0.Dispose()


    if bEcho and len(curvesAndEdges0) > 1:
        s = "Out of {} curves selected:".format(len(curvesAndEdges0))
        for sFail in set(sFails):
            s += "\n[{}] {}".format(sFails.count(sFail), sFail)
        if pts_Inflection_All:
            s += "\n{} inflection point(s) found.".format(pts_Inflection_All)
        if gPts_Inflection_Nested:
            s += "\n{} points added.".format(
                len(gNcs_Split_perCrv0), sum([len(gs) for gs in gNcs_Split_perCrv0]))
        if gNcs_Split_perCrv0:
            s += "\n{} curves were replaced with {}.".format(
                len(gNcs_Split_perCrv0), sum([len(gs) for gs in gNcs_Split_perCrv0]))
        if not gNcs_Split_perCrv0:
            s += "\nNo curves were added or replaced."
        print s
    
    return gNcs_Split_perCrv0


def main():
    
    rc = getInput()
    if rc is None: return
    objrefs = rc[0]
    
    if Opts.values['bDebug']:
        pass
    else:
        sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll()
    
    rc = processCurves(curvesAndEdges0=objrefs)
    if not rc: return
    
    gNcs_Split_perCrv0 = rc
    
    if gNcs_Split_perCrv0:
        for gs in gNcs_Split_perCrv0:
            for g in gs:
                sc.doc.Objects.Select(objectId=g)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()