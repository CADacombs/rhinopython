"""
This script is an alternative to _Connect.

Unlike _Connect,
1. It never converts curves to another type, e.g., from NURBS to lines or arcs.
2. It can solve more NURBS-to-NURBS connects.
3. If the intersection does not exist, the curves are still extended to their closest points.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250106-07: Created.

Typically, do not use
Curve Extend(CurveEnd side, CurveExtensionStyle style, Point3d endPoint)
because it creates a NurbsCurve extension to a Polycurve for non-NurbsCurve input.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


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
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves near their ends")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    go.OneByOnePostSelect = True
    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

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


def formatDistance(fDistance, iPrecision=6):
    if fDistance is None:
        return "(No deviation provided)"

    if fDistance < 10**-iPrecision:
        return "{:.{}e}".format(fDistance, 0)

    if fDistance < 0.1:
        return "{:.{}g}".format(fDistance, iPrecision)

    return "{:.{}g}".format(fDistance, iPrecision)


def createWorkingSegment(rgC_In, curveEnd):
    seg = rgC_In.Duplicate()
    if isinstance(seg, rg.PolyCurve):
        seg.RemoveNesting()
        idx_Seg = seg.SegmentCount-1 if curveEnd==rg.CurveEnd.End else 0
        seg = seg.SegmentCurve(idx_Seg) # Shouldn't be a PolyCurve.
        if isinstance(seg, rg.PolyCurve):
            raise Exception("Shouldn't be a PolyCurve.")
        #sEval = "rgC_WIP.GetType().Name"; print(sEval,'=', eval(sEval))
    if isinstance(seg, rg.NurbsCurve):
        if seg.SpanCount > 1:
            idx_Span = seg.SpanCount-1 if curveEnd==rg.CurveEnd.End else 0
            domain = seg.SpanDomain(idx_Span)
            seg = seg.Trim(domain)
    return seg


def distanceBetweenEnds(rgCrv_A, rgCrv_B, curveEnd_A, curveEnd_B):
    ptA = rgCrv_A.PointAtEnd if curveEnd_A == rg.CurveEnd.End else rgCrv_A.PointAtStart
    ptB = rgCrv_B.PointAtEnd if curveEnd_B == rg.CurveEnd.End else rgCrv_B.PointAtStart
    return ptA.DistanceTo(ptB)


def processCurves(rgCrv_In_A, rgCrv_In_B, curveEnd_A, curveEnd_B, bEcho=True, bDebug=False):
    """
    Returns tuple of string, string: (strShortContinuityDescription, strLongContinuityDescription)
    """

    # Set up to use loops.
    rgCs_In = rgCrv_In_A, rgCrv_In_B
    curveEnds = curveEnd_A, curveEnd_B
    #sEval = "curveEnds"; print(sEval,'=', eval(sEval))
    rgCs_WIP = [None, None]
    planes = [None, None]

    for iC, rgC_In in enumerate(rgCs_In):
        curveEnd = curveEnds[iC]
        rgCs_WIP[iC] = createWorkingSegment(rgC_In, curveEnds[iC])
        #sc.doc.Objects.AddCurve(rgCs_WIP[iC])

        bSuccess, planes[iC] = rgCs_WIP[iC].PerpendicularFrameAt(
            rgCs_WIP[iC].Domain.T1 if curveEnd==rg.CurveEnd.End else rgCs_WIP[iC].Domain.T0)
        if not bSuccess: raise Exception("Failed to obtain frame from curve [{}].".format(iC))

    fDistBetweenEnds_In = planes[0].Origin.DistanceTo(planes[1].Origin)
    sEval = "fDistBetweenEnds_In"; print(sEval,'=', eval(sEval))

    rgCs_Extension = [None, None]

    for iC, rgC_WIP in enumerate(rgCs_WIP):
        curveEnd = curveEnds[iC]
        plane_Opp = planes[int(not bool(iC))]
        ps = rg.PlaneSurface(
            plane_Opp,
            xExtents=rg.Interval(-1.0*fDistBetweenEnds_In, 1.0*fDistBetweenEnds_In),
            yExtents=rg.Interval(-1.0*fDistBetweenEnds_In, 1.0*fDistBetweenEnds_In))
        #sc.doc.Objects.AddSurface(ps)

        rgC_WIP_Extnsn = rgC_WIP.Extend(
            side=curveEnd,
            style=rg.CurveExtensionStyle.Smooth,
            geometry=(ps,))
        #sc.doc.Objects.AddCurve(rgC_WIP_Extnsn)
        #sEval = "rgC_WIP_Extnsn"; print(sEval,'=', eval(sEval))
        rgCs_Extension[iC] = rgC_WIP_Extnsn

    pts = [None, None]
    ts = [None, None]

    bSuccess, pts[0], pts[1] = rgCs_Extension[0].ClosestPoints(rgCs_Extension[1])
    if not bSuccess:
        print("Closest points not found.")
        return

    for iC, rgC_WIP_Extnsn in enumerate(rgCs_Extension):
        bSuccess, t = rgC_WIP_Extnsn.ClosestPoint(pts[iC])
        if bSuccess:
            ts[iC] = t

    for _ in rgCs_Extension: _.Dispose()

    rgCs_Out = [None, None]

    for iC, rgC_In in enumerate(rgCs_In):
        #sEval = "rgC_In.Domain"; print(sEval,'=', eval(sEval))
        domain_Ext = rg.Interval(rgC_In.Domain)
        #sEval = "ts[iC]"; print(sEval,'=', eval(sEval))
        domain_Ext.Grow(ts[iC])
        #sEval = "domain_Ext"; print(sEval,'=', eval(sEval))
        rgC_Final_Extnsn = rgC_In.Extend(domain=domain_Ext)
        #rgC_Final_Extnsn = rgC_In.Extend(
        #    side=curveEnds[iC],
        #    style=rg.CurveExtensionStyle.Smooth,
        #    endPoint=pts[iC])
        rgCs_Out[iC] = rgC_Final_Extnsn


    #sEval = "rvs"; print(sEval,'=', eval(sEval))

    return rgCs_Out


def processCurveObjects(rhCrv_In_A, rhCrv_In_B, curveEnd_A, curveEnd_B, bReplace=True, bEcho=True, bDebug=False):
    """
    """

    rdC_In_A = rs.coercerhinoobject(rhCrv_In_A)
    rdC_In_B = rs.coercerhinoobject(rhCrv_In_B)
    rgC_In_A = rdC_In_A.Geometry
    rgC_In_B = rdC_In_B.Geometry

    rgCs_Res = processCurves(
        rgC_In_A,
        rgC_In_B,
        curveEnd_A=curveEnd_A,
        curveEnd_B=curveEnd_B,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rgCs_Res is None:
        return

    if not any(rgCs_Res):
        print("Neither curve was extended.")
        return


    rdCs_In = (rdC_In_A, rdC_In_B)


    gOuts = [None, None]

    for iC, rgC_Out in enumerate(rgCs_Res):
        if not rgC_Out:
            continue

        if bReplace:
            gOut = rdCs_In[iC].Id
            bReplaced = sc.doc.Objects.Replace(gOut, rgC_Out)
            if bReplace:
                gOuts[iC] = gOut
            else:
                if bEcho: print("Curve could not be replaced.")
        else:
            gOut = sc.doc.Objects.AddCurve(rgC_Out)
            if gOut == Guid.Empty:
                if bEcho: print("Curve could not be added.")
            else:
                gOuts[iC] = gOut


    if bEcho:
        print("Distance between ends: {} -> {}".format(
            formatDistance(
                distanceBetweenEnds(rgC_In_A, rgC_In_B, curveEnd_A, curveEnd_B)
                ),
            formatDistance(
                distanceBetweenEnds(
                    rgCs_Res[0] if rgCs_Res[0] else rgC_In_A,
                    rgCs_Res[1] if rgCs_Res[1] else rgC_In_B,
                    curveEnd_A,
                    curveEnd_B)
                )
            )
              )

    #sEval = "gOuts"; print(sEval,'=', eval(sEval))

    #sc.doc.Objects.UnselectAll()
    #if sc.doc.Objects.Select(rdC_In_A.Id):
    #    s += " and is selected."
    #else:
    #    s += " but could not be selected."

    #print(s)

    return gOuts


def curveEnd_closestToPick(objref_Crv):
    rgCrv = objref_Crv.Curve()

    bSuccess, t = rgCrv.ClosestPoint(objref_Crv.SelectionPoint())
    if not bSuccess:
        return

    return rg.CurveEnd.End if t >= rgCrv.Domain.Mid else rg.CurveEnd.Start


def main():

    objrefs = getInput()
    if objrefs is None: return

    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    curveEnd_A = curveEnd_closestToPick(objrefs[0])
    curveEnd_B = curveEnd_closestToPick(objrefs[1])

    processCurveObjects(
        rhCrv_In_A=objrefs[0],
        rhCrv_In_B=objrefs[1],
        curveEnd_A=curveEnd_A,
        curveEnd_B=curveEnd_B,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()