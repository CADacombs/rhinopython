"""
This script will convert a curve to a
Piecewise Uniform Rational Basis Spline (URBS) curve.

Periodic curves' domains will (most likely?) be modified due to moving their seams
to a parameter to be split.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

#! python 2

"""
240819: Created starting with the 'Bezier' script.
240828: Bug fix to correctly handle closed and periodic curves.
        Fixed crash when no output.
        Added support for BrepEdge input. New curves are added, the edge curves are not replaced within the Brep.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

sc.escape_test(throw_exception=False, reset=True) # Added due to a bug in RhinoCode (at least up through 8.10).


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bExplode'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bNurbs_NotPoly'; keys.append(key)
    values[key] = True
    names[key] = 'Curve'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Poly', 'Nurbs')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteInput'; keys.append(key)
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

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

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = rd.ObjectType.Curve
    # go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bExplode')
        if not Opts.values['bExplode']:
            addOption('bNurbs_NotPoly')
        addOption('bDeleteInput')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def formatDistance(fDistance, iPrecision=15):
    if fDistance is None: return "(No deviation provided)"
    
    if fDistance < 0.01:
        return "{:.{}e}".format(fDistance, iPrecision)
    
    if fDistance < 0.1:
        return "{:.{}g}".format(fDistance, iPrecision+1)
    
    return "{:.{}g}".format(fDistance, iPrecision)


def createCrvs(rgCrv_In, bExplode, bNurbs_NotPoly, bDebug=False):
    """
    Returns on success: rg.NurbsCurve, None
    Returns on fail: None, str(Log)
    """

    if isinstance(rgCrv_In, rg.PolyCurve):
        rgC_WIP = rgCrv_In.CleanUp()
        if rgC_WIP is None:
            rgC_WIP = rgCrv_In.DuplicateCurve()
    else:
        rgC_WIP = rgCrv_In.DuplicateCurve()

    if isinstance(rgC_WIP, (rg.ArcCurve, rg.LineCurve, rg.PolylineCurve)):
        rgC_WIP.Dispose()
        return [], "{} was skipped.".format(rgC_WIP.GetType().Name)

    if rgC_WIP.SpanCount == 1:
        rgC_WIP.Dispose()
        return [], "{} with single span skipped.".format(rgC_WIP.GetType().Name)


    # Only PolyCurves and multi-span NurbsCurves should remain.

    # Check whether PolyCurve has consecutive NurbsCurves.
    if isinstance(rgC_WIP, rg.NurbsCurve):
        nc_WIP = rgC_WIP
    elif isinstance(rgC_WIP, rg.PolyCurve):
        idx_LastNurbs = None
        for i in range(rgC_WIP.SegmentCount):
            seg = rgC_WIP.SegmentCurve(i)
            if isinstance(seg, rg.NurbsCurve):
                if seg.Knots.KnotStyle not in (
                    rg.KnotStyle.QuasiUniform,
                    rg.KnotStyle.Uniform
                    ):
                    break # out of for loop.
        else:
            rgC_WIP.Dispose()
            return [], "No segments in PolyCurve are non-uniform NURBS."
        nc_WIP = rgC_WIP.ToNurbsCurve()
        rgC_WIP.Dispose()
    else:
        raise Exception("{}!".format(rgC_WIP.GetType().Name))

    # nc_WIP is now a non-uniform NurbsCurve.

    def parametersBetweenUniformSpans(nc):
        if nc.IsClosed and not nc.IsPeriodic:
            ts = [nc.Knots[0]]
        else:
            ts = []
        
        if nc.IsPeriodic:
            iK = nc.Degree - 1
            k = nc.Knots[iK]
            m = nc.Knots.KnotMultiplicity(iK)
            if m > 1:
                ts.append(nc.Knots[iK])

            iK = nc.Degree # nc.Knots.KnotMultiplicity(0) # First interior knot, not starting end.
            k = nc.Knots[iK]
            fSpanLength_Ref = nc.Knots[iK] - nc.Knots[iK-1]
        else:
            iK = nc.Degree # nc.Knots.KnotMultiplicity(0) # First interior knot, not starting end.
            fSpanLength_Ref = None
        while iK < nc.Knots.Count:
            sc.escape_test(throw_exception=True, reset=True)

            k = nc.Knots[iK]
            m = nc.Knots.KnotMultiplicity(iK)
            fSpanLength_Now = nc.Knots[iK] - nc.Knots[iK-1]

            if fSpanLength_Ref is None:
                fSpanLength_Ref = fSpanLength_Now
            elif abs(fSpanLength_Now - fSpanLength_Ref) > Rhino.RhinoMath.ZeroTolerance:
                ts.append(nc.Knots[iK-1])

            if nc.Knots[iK] == nc.Domain.T1:
                break

            if m > 1:
                ts.append(nc.Knots[iK])
                fSpanLength_Ref = None
            else:
                fSpanLength_Ref = fSpanLength_Now

            iK += m

        if len(ts) == 0:
            return []

        # if nc.IsPeriodic:
        #     if nc.Domain.T0 not in ts:
        #         nc.ChangeClosedCurveSeam(ts[0])

        return ts


    ts = parametersBetweenUniformSpans(nc_WIP)
    # sEval = "ts"; print("{}: {}".format(sEval, eval(sEval)))
    if len(ts) == 0:
        return [], "No non-uniform spans found."

    # TODO: Implement moving the seam to a split parameter.
    #       Through V8.10, using the following code generates a knot vector - curve domain mismatch.
    if nc_WIP.IsPeriodic:
        if nc_WIP.Domain.T0 not in ts:
            print(nc_WIP.ChangeClosedCurveSeam(ts[0]))
            # Generate the list of parameters again since the domain has been changed.
            ts = parametersBetweenUniformSpans(nc_WIP)
            ts.append(nc_WIP.Domain.T1) # Since T1 (T0) was confirmed to be a parameter for splitting.

    # if nc_WIP.IsPeriodic:
    #     if nc_WIP.Domain.T0 not in ts:
    #         print(nc_WIP.ChangeClosedCurveSeam(ts[0]))
    #         sc.doc.Objects.AddCurve(nc_WIP); sc.doc.Views.Redraw(); 1/0


    ncs_URBS = nc_WIP.Split(ts)

    if bExplode:
        return ncs_URBS, None

    pcs = rg.Curve.JoinCurves(ncs_URBS, joinTolerance=1e-6, preserveDirection=True)

    if len(pcs) != 1:
        raise Exception("{} curves returned from Curve.JoinCurves.".format(len(pcs)))

    pc_Res = pcs[0]

    if not bNurbs_NotPoly:
        return [pc_Res], None

    nc_Out = pc_Res.ToNurbsCurve()
    pc_Res.Dispose()

    return [nc_Out], None


def processCurveObject(rhCrv_In, bExplode, bNurbs_NotPoly, bDeleteInput, bEcho=True, bDebug=False):
    """
    """

    if isinstance(rhCrv_In, rd.ObjRef):
        rgCrv_In = rhCrv_In.Curve()
        if rgCrv_In is None:
            return [], "Curve could not be obtained from ObjRef."
        bDeleteInput = False
    elif isinstance(rhCrv_In, rd.CurveObject):
        rgCrv_In = rhCrv_In.CurveGeometry
        if rgCrv_In is None:
            return [], "Curve could not be obtained from CurveObject."
    else:
        raise Exception("{} was passed to processCurveObject".format(rgCrv_In.GetType().Name))
    # rdCrv_In = rs.coercerhinoobject(rhCrv_In) # IsDocumentControlled.
    # rgCrv_In = rdCrv_In.Geometry # IsDocumentControlled.

    cs_Res, sLog = createCrvs(
        rgCrv_In,
        bExplode=bExplode,
        bNurbs_NotPoly=bNurbs_NotPoly,
        bDebug=bDebug,
        )
    if len(cs_Res) == 0:
        return [], [] if sLog is None else [sLog]

    gOuts = []
    sOuts = []

    bAllAddSuccess = True

    if bDeleteInput:
        if len(cs_Res) == 1:
            if not sc.doc.Objects.Replace(objectId=rdCrv_In.Id, curve=cs_Res[0]):
                sOut = "Replace failed."
                if bEcho: print(sOut)
                return None, sOut
            sOuts.append("Curve was replaced")
        else:
            for c in cs_Res:
                gOut = sc.doc.Objects.AddCurve(c)
                if gOut != gOut.Empty:
                    gOuts.append(gOut)
                    sOuts.append("Curve was added.")
                else:
                    bAllAddSuccess = False
                    sOuts.append("AddCurve failed.")
            if bAllAddSuccess:
                sc.doc.Objects.Delete(rdCrv_In)
                # sOut = "Curve was replaced."
    else:
        # bDeleteInput == False
        for c in cs_Res:
            gOut = sc.doc.Objects.AddCurve(c)
            if gOut != gOut.Empty:
                gOuts.append(gOut)
                sOuts.append("Curve was added.")
            else:
                bAllAddSuccess = False
                sOuts.append("AddCurve failed.")

    if not bAllAddSuccess:
        return gOuts, sOuts

    # if sc.doc.Objects.Select(gOut):
    #     sOut += " and is selected."
    # else:
    #     sOut += " but could not be selected."

    # if bEcho: print(sOut)

    return gOuts, sOuts


def main():

    objrefs = getInput()
    if objrefs is None: return

    bExplode = Opts.values['bExplode']
    bNurbs_NotPoly = Opts.values['bNurbs_NotPoly']
    bDeleteInput = Opts.values['bDeleteInput']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gCrvs_Res = []
    sLogs = []

    for objref in objrefs:
        gCrvs_Ret, sLogs_Ret = processCurveObject(
            objref,
            bExplode=bExplode,
            bNurbs_NotPoly=bNurbs_NotPoly,
            bDeleteInput=bDeleteInput,
            bEcho=bEcho if len(objrefs) == 1 else False,
            bDebug=bDebug,
            )
        gCrvs_Res.extend(gCrvs_Ret)
        sLogs.extend(sLogs_Ret)

    sc.doc.Views.RedrawEnabled = True

    if bEcho and sLogs:
        if len(sLogs) == 1:
            print(sLogs[0])
        else:
            for sLog in set(sLogs):
                print("[{}] {}".format(sLogs.count(sLog), sLog))


if __name__ == '__main__': main()