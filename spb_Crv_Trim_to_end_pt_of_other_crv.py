"""
This script trims one curve at the end point of a reference curve closest to
the picked point on the latter.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250323-24: Created.
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


def getInput_Ref():
    """
    Get curve to use its end for trimming.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve near its end point to use")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return objref

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_to_trim(objref_Ref):
    """
    Get curve to trim.
    """

    edge_Ref = objref_Ref.Edge()
    idxE_Ref = None if edge_Ref is None else edge_Ref.EdgeIndex
    gRef = objref_Ref.ObjectId

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve to trim")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    def customGeomFilter(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if rdObj.Id == gRef:
            if isinstance(geom, rg.BrepEdge):
                # GUID may be the same as reference curve, so only edge index need be different.
                if geom.EdgeIndex == idxE_Ref:
                    #print("Same edge index {}".format(idxE_Ref))
                    return False
            else:
                return False

        return True

    go.SetCustomGeometryFilter(customGeomFilter)

    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return objref

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


def trimCurve(rgCrv_TrimMe, t_on_trim_side, t_closest_to_trimming_pt, bDebug=False):
    """
    Parameters:
        rgCrv_TrimMe
        t_on_trim_side
        bDebug: bool

    Returns: tuple of 2 items:
        list(rg.Curves) or None
        str(Feedback on fail, etc.) or None
    """

    if t_on_trim_side > t_closest_to_trimming_pt:
        return rgCrv_TrimMe.Trim(
            t0=rgCrv_TrimMe.Domain.T0,
            t1=t_closest_to_trimming_pt)

    return rgCrv_TrimMe.Trim(
        t0=t_closest_to_trimming_pt,
        t1=rgCrv_TrimMe.Domain.T1)


def _pt_at_curve_end_closest_to_pick(objref_Crv):
    rgCrv = objref_Crv.Curve()

    bSuccess, t = rgCrv.ClosestPoint(objref_Crv.SelectionPoint())
    if not bSuccess:
        return

    t_MidLength = rgCrv.DivideByCount(segmentCount=2, includeEnds=False)[0]

    return rgCrv.PointAtEnd if t >= t_MidLength else rgCrv.PointAtStart


def _closest_pt(rgCurve, testPoint):
    bSuccess, t = rgCurve.ClosestPoint(
        testPoint=testPoint)
    if not bSuccess: raise Exception("ClosestPoint failed.")
    return t


def processObjRefs(objref_TrimMe, objref_Ref, bReplace=True, bEcho=True, bDebug=False):
    """
    """

    edge_TrimMe = objref_TrimMe.Edge()

    pt_on_trim_side = objref_TrimMe.SelectionPoint()
    rgC_TrimMe, t_on_trim_side = objref_TrimMe.CurveParameter()
    pt_Trimming = _pt_at_curve_end_closest_to_pick(objref_Ref)
    t_closest_to_trimming_pt = _closest_pt(
        rgC_TrimMe,
        testPoint=pt_Trimming)
    fDist_between_pts = rgC_TrimMe.PointAt(t_closest_to_trimming_pt).DistanceTo(pt_Trimming)

    #maximumDistance = sc.doc.ModelAbsoluteTolerance if tolerance is None else tolerance
    fMaxDistFromCrv = sc.doc.ModelAbsoluteTolerance

    print("End point is {} from curve to trim.".format(_formatDistance(fDist_between_pts)))

    if fDist_between_pts > fMaxDistFromCrv:
        print("Curve not trimmed.")
        return


    rdC_TrimMe = objref_TrimMe.Object()
    #rgC_TrimMe = rdC_TrimMe.Geometry

    rgC_Out = trimCurve(
        rgCrv_TrimMe=rgC_TrimMe,
        t_on_trim_side=t_on_trim_side,
        t_closest_to_trimming_pt=t_closest_to_trimming_pt,
        bDebug=bDebug,
        )
    if rgC_Out is None:
        return


    if bReplace and edge_TrimMe is None:
        gOut = rdC_TrimMe.Id
        bReplaced = sc.doc.Objects.Replace(gOut, rgC_Out)
        if bReplaced:
            if bEcho: print("Curve was replaced.")
        else:
            gOut = None
            if bEcho: print("Curve could not be replaced.")
    else:
        gOut = sc.doc.Objects.AddCurve(rgC_Out)
        if gOut == Guid.Empty:
            gOut = None
            if bEcho: print("Curve could not be added.")
        else:
            if bEcho: print("Curve was added.")

    return gOut


def _curveEnd_closestToPick(objref_Crv):
    rgCrv = objref_Crv.Curve()

    bSuccess, t = rgCrv.ClosestPoint(objref_Crv.SelectionPoint())
    if not bSuccess:
        return

    t_MidLength = rgCrv.DivideByCount(segmentCount=2, includeEnds=False)[0]

    return rg.CurveEnd.End if t >= t_MidLength else rg.CurveEnd.Start


def main():

    objref_Ref = getInput_Ref()
    if objref_Ref is None: return

    objref_TrimMe = getInput_to_trim(objref_Ref)
    if objref_TrimMe is None: return

    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    processObjRefs(
        objref_TrimMe=objref_TrimMe,
        objref_Ref=objref_Ref,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()