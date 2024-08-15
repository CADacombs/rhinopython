"""
This script will convert a curve in a Piecewise-Bezier NurbsCurve and set all
span domains so that consecutive spans have parametrically-matched domains as so:

Starting with the second span/segment (M), with R as the previous span/segment,
this script will match the domain of M so that 
D(M) = D(R) * (L(M) / L(R))
where
D is the end span domain length
L is the distance between the end control point and next control point

The resultant NurbsCurve will allow knot removal with less curve deviation than
a non-domain matched version.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

#! python 2

"""
240107: Created.
240116: Bug fix.
240813: Raised tolerance threshold before parameter matching.
240814: Now can output PolyCurves or individual Bezier (single-spanned) NurbsCurves.

TODO: Add option for resultant curve to be a PolyCurve instead of a NurbsCurve.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bExplode'; keys.append(key)
    values[key] = True
    # names[key] = 'ExplodeToBeziers'
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

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

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
        return None, "{} was skipped.".format(rgC_WIP.GetType().Name)

    if rgC_WIP.SpanCount == 1:
        rgC_WIP.Dispose()
        return None, "{} with single span skipped.".format(rgC_WIP.GetType().Name)


    # Only PolyCurves and multi-span NurbsCurves should remain.

    # Check whether PolyCurve has consecutive NurbsCurves.
    if isinstance(rgC_WIP, rg.NurbsCurve):
        nc_WIP = rgC_WIP
    elif isinstance(rgC_WIP, rg.PolyCurve):
        idx_LastNurbs = None
        for i in range(rgC_WIP.SegmentCount):
            seg = rgC_WIP.SegmentCurve(i)
            if isinstance(seg, rg.NurbsCurve):
                if seg.SpanCount > 1:
                    break # out of for loop.
                else:
                    if idx_LastNurbs is None:
                        idx_LastNurbs = i
                    else:
                        if idx_LastNurbs == (i - 1):
                            # Consecutive NurbsCurves found.
                            break # out of for loop.
        else:
            rgC_WIP.Dispose()
            return None, "No segments in PolyCurve are multi-span NURBS."
        nc_WIP = rgC_WIP.ToNurbsCurve()
        rgC_WIP.Dispose()
    else:
        raise Exception("{}!".format(rgC_WIP.GetType().Name))

    # nc_WIP is now a multi-span NurbsCurve.

    ncs_Beziers = [bc.ToNurbsCurve() for bc in rg.BezierCurve.CreateBeziers(nc_WIP)]


    for i in range(1, len(ncs_Beziers)):
        nc_R = ncs_Beziers[i-1] # Reference.
        nc_M = ncs_Beziers[i] # To modify.

        deg_R = nc_R.Degree
        deg_M = nc_M.Degree

        #sEval = "deg_R"; print(sEval+':',eval(sEval))
        #sEval = "deg_M"; print(sEval+':',eval(sEval))

        domain_R = nc_R.Domain
        idxCp_Pos_R = nc_R.Points.Count - 1
        idxCp_Tan_R = nc_R.Points.Count - 2
        pt_Pos_R = nc_R.Points[idxCp_Pos_R].Location
        pt_Tan_R = nc_R.Points[idxCp_Tan_R].Location
        fDist_CPs_R = pt_Pos_R.DistanceTo(pt_Tan_R)

        domain_M_Pre = nc_M.Domain
        idxCp_Pos_M = 0
        idxCp_Tan_M = 1
        pt_Pos_M = nc_M.Points[idxCp_Pos_M].Location
        pt_Tan_M = nc_M.Points[idxCp_Tan_M].Location
        fDist_CPs_M = pt_Pos_M.DistanceTo(pt_Tan_M)


        # D(A) = D(R) * (L(A) / L(R))

        m = (
            (domain_R.Length / domain_M_Pre.Length) * 
            (fDist_CPs_M / fDist_CPs_R)
            )

        if bDebug:
            sEval = "m"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))
            sEval = "m - 1.0"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))

        if abs(m - 1.0) == 0.0:
            if bDebug: print("Domains already exactly match.")
            continue
        if abs(m - 1.0) <= 2**-53:
            if bDebug:
                print("Domains already match within {} (machine epsilon/2.0).".format(
                    formatDistance(2**-53)))
            continue
        if abs(m - 1.0) <= 2**-52:
            if bDebug:
                print("Domains already match within {} (machine epsilon).".format(
                    formatDistance(2**-52)))
            continue
        if abs(m - 1.0) <= 2**-51:
            if bDebug:
                print("Domains already match within {}.".format(
                    formatDistance(2**-51)))
            continue
        if abs(m - 1.0) <= 2**-50:
            if bDebug:
                print("Domains already match within {}.".format(
                    formatDistance(2**-50)))
            continue
        #if abs(m - 1.0) <= Rhino.RhinoMath.ZeroTolerance:
        #    if bDebug:
        #        print("Domains already match within {}.".format(
        #            formatDistance(abs(m - 1.0))))
        #    continue

        if bDebug:
            print("Domain length multiplier to apply: {}".format(
                formatDistance(m)))

        length_Domain_M_Out = domain_M_Pre.Length * m
        nc_M.Domain = rg.Interval(0.0, length_Domain_M_Out)
        domain_M_Post = nc_M.Domain

    if bExplode:
        return ncs_Beziers, None

    pcs = rg.Curve.JoinCurves(ncs_Beziers, joinTolerance=1e-6, preserveDirection=True)

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

    rdCrv_In = rs.coercerhinoobject(rhCrv_In) # IsDocumentControlled.
    rgCrv_In = rdCrv_In.Geometry # IsDocumentControlled.

    cs_Res, sLog = createCrvs(
        rgCrv_In,
        bExplode=bExplode,
        bNurbs_NotPoly=bNurbs_NotPoly,
        bDebug=bDebug,
        )
    if not cs_Res:
        return None, sLog

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
        rc = processCurveObject(
            objref,
            bExplode=bExplode,
            bNurbs_NotPoly=bNurbs_NotPoly,
            bDeleteInput=bDeleteInput,
            bEcho=bEcho if len(objrefs) == 1 else False,
            bDebug=bDebug,
            )
        if rc is None: continue
        gCrvs_Res.extend(rc[0])
        sLogs.extend(rc[1])

    sc.doc.Views.RedrawEnabled = True

    if bEcho and sLogs:
        if len(sLogs) == 1:
            print(sLogs[0])
        else:
            for sLog in set(sLogs):
                print("[{}] {}".format(sLogs.count(sLog), sLog))


if __name__ == '__main__': main()