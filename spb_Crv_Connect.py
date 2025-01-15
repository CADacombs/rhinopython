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
250112-13: Added routine to support curves that already overlap. Refactored.
250114: Added routine to lengthen each curve when length to PlaneSurface fails for both curves.

TODO:
    Add coincident tolerance, e.g., max((0.1*sc.ModelDistanceTolerance, 1e-4))

Typically, do not use
Curve Extend(CurveEnd side, CurveExtensionStyle style, Point3d endPoint)
because it creates a NurbsCurve extension to a Polycurve for non-NurbsCurve input.
This is probably what the _ToPoint option of _Extend uses.
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


def _formatDistance(fDistance, iPrecision=6):
    if fDistance is None:
        return "(No deviation provided)"

    if fDistance < 10**-iPrecision:
        return "{:.{}e}".format(fDistance, 0)

    if fDistance < 0.1:
        return "{:.{}g}".format(fDistance, iPrecision)

    return "{:.{}g}".format(fDistance, iPrecision)


def _createWorkingHalf(rgC_In, curveEnd):

    t_MidLength = rgC_In.DivideByCount(segmentCount=2, includeEnds=False)[0]

    return rgC_In.Trim(
        rgC_In.Domain.T0 if curveEnd==rg.CurveEnd.Start else t_MidLength,
        t_MidLength if curveEnd==rg.CurveEnd.Start else rgC_In.Domain.T1)


def _createWorkingSegment(rgC_In, curveEnd):
    if isinstance(rgC_In, (rg.ArcCurve, rg.LineCurve)):
        return rgC_In.Duplicate()

    if isinstance(rgC_In, rg.PolyCurve):
        pc_WIP = rgC_In.Duplicate()
        pc_WIP.RemoveNesting()
        idx_Seg = pc_WIP.SegmentCount-1 if curveEnd==rg.CurveEnd.End else 0
        seg = pc_WIP.SegmentCurve(idx_Seg) # Shouldn't be a PolyCurve.
        if isinstance(seg, rg.PolyCurve):
            raise Exception("Shouldn't be a PolyCurve.")
        #sEval = "rgC_WIP.GetType().Name"; print(sEval,'=', eval(sEval))
        if isinstance(seg, (rg.ArcCurve, rg.LineCurve)):
            return seg

    if isinstance(rgC_In, rg.PolylineCurve):
        return rgC_In.DuplicateSegments()[-1]

    if isinstance(rgC_In, rg.NurbsCurve):
        if seg.SpanCount > 1:
            idx_Span = seg.SpanCount-1 if curveEnd==rg.CurveEnd.End else 0
            domain = seg.SpanDomain(idx_Span)
            seg = seg.Trim(domain)

    return seg


def _arePtsWithinTolerance(pts, tolerance=1e-6):
    fDistBetweenEnd_and_ClosestPt = pts[0].DistanceTo(pts[1])
    return fDistBetweenEnd_and_ClosestPt <= tolerance


def _pt_atWorkingEnd(rgCrv, curveEnd):
    if curveEnd == rg.CurveEnd.Start:
        return rgCrv.PointAtStart
    if curveEnd == rg.CurveEnd.End:
        return rgCrv.PointAtEnd
    raise Exception("Wrong argument(s) passed to _pt_atWorkingEnd?")


def _vector_atWorkingEnd(rgCrv, curveEnd):
    if curveEnd == rg.CurveEnd.Start:
        return rgCrv.TangentAtStart
    if curveEnd == rg.CurveEnd.End:
        return rgCrv.TangentAtEnd
    raise Exception("Wrong argument(s) passed to _pt_atWorkingEnd?")


def _distanceBetweenEnds(rgCrvs, curveEnds):
    ptA = _pt_atWorkingEnd(rgCrvs[0], curveEnds[0])
    ptB = _pt_atWorkingEnd(rgCrvs[1], curveEnds[1])
    return ptA.DistanceTo(ptB)


def _isPtAtEndOfCrv(rgCrv, curveEnd, pt, tolerance=1e-6):
    pt_WorkingEnd = _pt_atWorkingEnd(rgCrv, curveEnd)
    return _arePtsWithinTolerance((pt, pt_WorkingEnd), tolerance=tolerance)


def _currentClosestPtParameters(rgCrvs, curveEnds, tolerance=1e-6):
    pts = [None, None]

    bSuccess, pts[0], pts[1] = rgCrvs[0].ClosestPoints(rgCrvs[1])
    if not bSuccess:
        raise Exception("Closest points not found for current curves. Why?")

    fDistBetweenClosestPts_Segs_In = pts[0].DistanceTo(pts[1])
    if fDistBetweenClosestPts_Segs_In > tolerance:
        return

    rgCs_Out = [None, None]

    ts = [None, None]

    for i in 0,1:

        if _isPtAtEndOfCrv(rgCrvs[i], curveEnds[i], pts[i], tolerance=tolerance):
            continue

        bSuccess, t = rgCrvs[i].ClosestPoint(pts[i])
        if not bSuccess:
            continue

        ts[i] = t

    return ts


def _trimCurveToParameter(rgCrv_In, curveEnd, t):
    if curveEnd == rg.CurveEnd.Start:
        return rgCrv_In.Trim(t, rgCrv_In.Domain.T1)

    if curveEnd == rg.CurveEnd.End:
        return rgCrv_In.Trim(rgCrv_In.Domain.T0, t)


def _dist_pt_to_closestPt_on_Crv(pt, crv):
    bSuccess, t = crv.ClosestPoint(pt)
    if not bSuccess:
        raise Exception("Closest point not found!")
    pt_Closest = rg.Curve.PointAt(crv, t)
    #sc.doc.Objects.AddPoint(pt_Closest)
    return pt.DistanceTo(pt_Closest)


def _extendByLength(rgCrv_toExtend, curveEnd_toExtend, rgCrv_Ref):
    # Closest points between end of curve and opposite curve.
    pt_End = _pt_atWorkingEnd(rgCrv_toExtend, curveEnd_toExtend)
    #sc.doc.Objects.AddPoint(pt_End)

    fDistFromCrvAEnd_toCrvB = _dist_pt_to_closestPt_on_Crv(pt_End, rgCrv_Ref)
    sEval = "fDistFromCrvAEnd_toCrvB"; print(sEval,'=', eval(sEval))
    #sEval = "curveEnd_toExtend"; print(sEval,'=', eval(sEval))


    rgC_Extended = rg.Curve.Extend(
        rgCrv_toExtend,
        side=curveEnd_toExtend,
        length=2.0*fDistFromCrvAEnd_toCrvB,
        style=rg.CurveExtensionStyle.Smooth)

    return rgC_Extended



    #line = rg.Line(
    #    start=pt_End,
    #    span=fDistFromCrvAEnd_toCrvB*_vector_atWorkingEnd(rgCrv_toExtend, curveEnd_toExtend))
    #lc = rg.LineCurve(line)
    #sc.doc.Objects.AddCurve(lc)
    #sc.doc.Views.Redraw(); 1/0

    # Extend by line.
    #rgC_LineExtended = rg.Curve.Extend(
    #    rgCrv_toExtend,
    #    side=curveEnd_toExtend,
    #    length=2.0*fDistFromCrvAEnd_toCrvB,
    #    style=rg.CurveExtensionStyle.Line)
    #sc.doc.Objects.AddCurve(rgC_LineExtended)
    #sc.doc.Views.Redraw(); 1/0

    # Closest points with line and opposite curve.



    # Use the distance between closest point on line to lengthen its origin curve. (1.5*distance)
    # Return the curve.
    #pass

    pass # In Visual Studio, a non-comment below comments allows code folding.


def processCurves(rgCrvs_In, curveEnds, tolerance=1e-6, bDebug=False):
    """
    Parameters:
        rgCrvs_In: list, tuple, etc., of 2 rg.Curves
        curveEnds: list, tuple, etc., of 2 rg.CurveEnds
        tolerance: float
        bDebug: bool

    Returns: tuple of 2 items:
        list(rg.Curves) or None
        str(Feedback on fail, etc.) or None
    """

    # First, test whether trim or extend is even needed.
    fDistBetwnEnds = _distanceBetweenEnds(rgCrvs_In, curveEnds)
    if fDistBetwnEnds <= tolerance:
        return None, "Distance between ends is already {}. Curve were not modified.".format(_formatDistance(fDistBetwnEnds))


    # Set up to use in loops.
    rgCs_In = rgCrvs_In
    #sEval = "curveEnds"; print(sEval,'=', eval(sEval))

    rgCs_WIP = [None, None]

    for i in 0,1:
        #rgCs_WIP[i] = _createWorkingSegment(rgCs_In[i], curveEnds[i])
        rgCs_WIP[i] = _createWorkingHalf(rgCs_In[i], curveEnds[i])
        #sc.doc.Objects.AddCurve(rgCs_WIP[i])


    # Before extending, try closest points of current segments.
    rv = _currentClosestPtParameters(rgCs_WIP, curveEnds)
    if rv:
        ts = rv
        rgCs_Out = [None, None]
        for i in 0,1:
            if ts[i] is None:
                continue
            rgCs_Out[i] = _trimCurveToParameter(rgCs_In[i], curveEnds[i], ts[i])

        return rgCs_Out, None


    # Extend then check for closest point.

    planes = [None, None]

    for i in 0,1:
        bSuccess, planes[i] = rgCs_WIP[i].PerpendicularFrameAt(
            rgCs_WIP[i].Domain.T1 if curveEnds[i]==rg.CurveEnd.End else rgCs_WIP[i].Domain.T0)
        if not bSuccess: raise Exception("Failed to obtain frame from curve [{}].".format(i))

    fDistBetweenEnds_In = planes[0].Origin.DistanceTo(planes[1].Origin)
    sEval = "fDistBetweenEnds_In"; print(sEval,'=', eval(sEval))

    rgCs_PostExtend = [None, None]
    # PostExtend just means after either or both have been extended.
    # If not both, a duplicate of the working half of the original is used instead.

    for i in 0,1:
        curveEnd = curveEnds[i]
        plane_Opp = planes[int(not bool(i))]
        ps = rg.PlaneSurface(
            plane_Opp,
            xExtents=rg.Interval(-1.0*fDistBetweenEnds_In, 1.0*fDistBetweenEnds_In),
            yExtents=rg.Interval(-1.0*fDistBetweenEnds_In, 1.0*fDistBetweenEnds_In))
        #sc.doc.Objects.AddSurface(ps)

        rgC_WIP_Extnsn = rgCs_WIP[i].Extend(
            side=curveEnd,
            style=rg.CurveExtensionStyle.Smooth,
            geometry=(ps,))
        #sc.doc.Objects.AddCurve(rgC_WIP_Extnsn)
        sEval = "rgC_WIP_Extnsn"; print(sEval,'=', eval(sEval))
        rgCs_PostExtend[i] = rgC_WIP_Extnsn


    if not any(rgCs_PostExtend):
        print("Neither curve was extended. WIP: Creating new routine for this condition, e.g., measure from end to opposite curve and Curve.Extend by a multiple of that distance.")

    if not all(rgCs_PostExtend):
        for i in 0,1:
            if not rgCs_PostExtend[i]:
                rgCs_PostExtend[i] = _extendByLength(rgCs_WIP[i], curveEnds[i], rgCs_WIP[(i+1)%2])

    pts = [None, None]

    bSuccess, pts[0], pts[1] = rgCs_PostExtend[0].ClosestPoints(rgCs_PostExtend[1])
    if not bSuccess:
        return None, "Closest points not found."

    ts = [None, None]

    for i in 0,1:
        bSuccess, t = rgCs_PostExtend[i].ClosestPoint(pts[i])
        if bSuccess:
            ts[i] = t

    for _ in rgCs_PostExtend: _.Dispose()

    rgCs_Out = [None, None]

    for i in 0,1:
        if rgCs_In[i].Domain.IncludesParameter(ts[i]):
            rgCs_Out[i] = _trimCurveToParameter(rgCs_In[i], curveEnds[i], ts[i])
        else:
            #sEval = "rgCs_In[i].Domain"; print(sEval,'=', eval(sEval))
            domain_Ext = rg.Interval(rgCs_In[i].Domain)
            #sEval = "ts[iC]"; print(sEval,'=', eval(sEval))
            domain_Ext.Grow(ts[i])
            #sEval = "domain_Ext"; print(sEval,'=', eval(sEval))
            rgC_Final_Extnsn = rgCs_In[i].Extend(domain=domain_Ext)
            #rgC_Final_Extnsn = rgCs_In[i].Extend(
            #    side=curveEnds[iC],
            #    style=rg.CurveExtensionStyle.Smooth,
            #    endPoint=pts[iC])
            rgCs_Out[i] = rgC_Final_Extnsn


    #sEval = "rvs"; print(sEval,'=', eval(sEval))

    return rgCs_Out, None


def processCurveObjects(rhCrv_In_A, rhCrv_In_B, curveEnd_A, curveEnd_B, bReplace=True, bEcho=True, bDebug=False):
    """
    """

    rdC_In_A = rs.coercerhinoobject(rhCrv_In_A)
    rdC_In_B = rs.coercerhinoobject(rhCrv_In_B)
    rgC_In_A = rdC_In_A.Geometry
    rgC_In_B = rdC_In_B.Geometry

    rv = processCurves(
        rgCrvs_In=(rgC_In_A, rgC_In_B),
        curveEnds=(curveEnd_A, curveEnd_B),
        tolerance=1e-6,
        bDebug=bDebug,
        )
    if rv is None:
        return

    rgCs_Res, sLog = rv

    if not any(rgCs_Res):
        print("Neither curve was modified.")
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
            _formatDistance(
                _distanceBetweenEnds((rgC_In_A, rgC_In_B), (curveEnd_A, curveEnd_B))
                ),
            _formatDistance(
                _distanceBetweenEnds(
                    (rgCs_Res[0] if rgCs_Res[0] else rgC_In_A,
                    rgCs_Res[1] if rgCs_Res[1] else rgC_In_B),
                    (curveEnd_A, curveEnd_B)))
            ))

    #sEval = "gOuts"; print(sEval,'=', eval(sEval))

    #sc.doc.Objects.UnselectAll()
    #if sc.doc.Objects.Select(rdC_In_A.Id):
    #    s += " and is selected."
    #else:
    #    s += " but could not be selected."

    #print(s)

    return gOuts


def _curveEnd_closestToPick(objref_Crv):
    rgCrv = objref_Crv.Curve()

    bSuccess, t = rgCrv.ClosestPoint(objref_Crv.SelectionPoint())
    if not bSuccess:
        return

    t_MidLength = rgCrv.DivideByCount(segmentCount=2, includeEnds=False)[0]

    return rg.CurveEnd.End if t >= t_MidLength else rg.CurveEnd.Start


def main():

    objrefs = getInput()
    if objrefs is None: return

    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    curveEnd_A = _curveEnd_closestToPick(objrefs[0])
    curveEnd_B = _curveEnd_closestToPick(objrefs[1])

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