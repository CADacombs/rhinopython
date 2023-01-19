"""
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220905-06, 08-09: Created.
"""

import Rhino
import Rhino.DocObjects as rd
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


    key = 'iDegree'; keys.append(key)
    values[key] = 3
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iPointCt'; keys.append(key)
    values[key] = 4
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPreserveEndTangentsForNonPolylines'; keys.append(key)
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

        if key == 'iDegree':
            if cls.riOpts[key].CurrentValue == cls.values[key]:
                # No change.
                return
            if cls.riOpts[key].CurrentValue <= 0:
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            cls.values['iPointCt'] = cls.riOpts['iPointCt'].CurrentValue = cls.values['iDegree'] + 1
            sc.sticky[cls.stickyKeys['iPointCt']] = cls.values['iPointCt']
            return

        if key == 'iPointCt':
            if cls.riOpts[key].CurrentValue < (cls.values['iDegree'] + 1):
                cls.values[key] = cls.riOpts[key].CurrentValue = (cls.values['iDegree'] + 1)
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
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    #    def customGeometryFilter(rdObj, rgObj, compIdx):
    #        return not rgObj.IsClosed
    #
    #    go.SetCustomGeometryFilter(customGeometryFilter)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('iDegree')
        addOption('iPointCt')
        addOption('bPreserveEndTangentsForNonPolylines')
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

        if res == ri.GetResult.Number:
            key = 'iDegree'
            Opts.riOpts[key].CurrentValue = int(abs(go.Number()))
            Opts.setValue(key)
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def approximateCurveWithPolyline(crv, target_pt_ct=6):

    min = 0.1 * sc.doc.ModelAbsoluteTolerance
    plc_WIP = crv.ToPolyline(
        tolerance=min,
        angleTolerance=0.0,
        minimumLength=0.0,
        maximumLength=0.0)
    ct_Min = plc_WIP.PointCount
    if ct_Min == target_pt_ct:
        return plc_WIP.ToPolyline()

    plc_WIP.Dispose()

    max = 10.0 * min
    while True:
        sc.escape_test()
        plc_WIP = crv.ToPolyline(
            tolerance=max,
            angleTolerance=0.0,
            minimumLength=0.0,
            maximumLength=0.0)
        ct_Max = plc_WIP.PointCount
        if ct_Max == target_pt_ct:
            return plc_WIP.ToPolyline()
        elif ct_Max < target_pt_ct:
            break
        max *= 10.0

    plc_WIP.Dispose()

    while True:
        sc.escape_test()
        
        mid = 0.5 * (min + max)
        
        #if abs(mid-min) <= Rhino.RhinoMath.ZeroTolerance:
        #    return
        
        plc_WIP = crv.ToPolyline(
            tolerance=mid,
            angleTolerance=0.0,
            minimumLength=0.0,
            maximumLength=0.0)
        
        ct_Mid = plc_WIP.PointCount
        if ct_Mid == target_pt_ct:
            return plc_WIP.ToPolyline()
        if ct_Mid > target_pt_ct:
            min = mid
            ct_Min = ct_Mid
        elif ct_Mid < target_pt_ct:
            max = mid
            ct_Max = ct_Mid
        else:
            raise Exception("What?")
        
        plc_WIP.Dispose()
        
        if abs(max-min) <= Rhino.RhinoMath.ZeroTolerance:
            plc_Max = crv.ToPolyline(
                tolerance=max,
                angleTolerance=0.0,
                minimumLength=0.0,
                maximumLength=0.0)
            return plc_Max.ToPolyline()


def samplePointsOnArcCurve(ac, pt_ct):
    ts = rg.ArcCurve.DivideByCount(
        ac, segmentCount=pt_ct-1, includeEnds=True)
    return Rhino.Collections.Point3dList(
        [rg.ArcCurve.PointAt(ac, t) for t in ts])


def simplifyPolyline(pl, target_pt_ct=6):
    ct_In = pl.Count
    if ct_In <= target_pt_ct:
        return
    
    min = 0.1 * sc.doc.ModelAbsoluteTolerance
    pl_WIP = pl.Duplicate()
    ct_Min = ct_In - pl_WIP.ReduceSegments(tolerance=min)
    if ct_Min == target_pt_ct:
        return pl_WIP
    
    max = 10.0 * min
    while True:
        sc.escape_test()
        pl_WIP = pl.Duplicate()
        ct_Max = ct_In - pl_WIP.ReduceSegments(tolerance=max)
        if ct_Max == target_pt_ct:
            return pl_WIP
        elif ct_Max < target_pt_ct:
            break
        max *= 10.0

    while True:
        sc.escape_test()
        
        mid = 0.5 * (min + max)
        
        #if abs(mid-min) <= Rhino.RhinoMath.ZeroTolerance:
        #    return
        
        pl_WIP = pl.Duplicate()
        ct_WIP = ct_In - pl_WIP.ReduceSegments(tolerance=mid)
        
        if ct_WIP == target_pt_ct:
            return pl_WIP
        if ct_WIP > target_pt_ct:
            min = mid
        elif ct_WIP < target_pt_ct:
            max = mid
        else:
            raise Exception("What?")
        
        if abs(max-min) <= Rhino.RhinoMath.ZeroTolerance:
            return


def getStartingCpLocations(rgCrv_In, pt_ct):

    if pt_ct == 2:
        return rg.Polyline([rgCrv_In.PointAtStart, rgCrv_In.PointAtEnd])

    if isinstance(rgCrv_In, rg.PolylineCurve):
        pl = rgCrv_In.ToPolyline()
        return simplifyPolyline(pl, target_pt_ct=pt_ct)
    elif isinstance(rgCrv_In, rg.NurbsCurve):
        return approximateCurveWithPolyline(rgCrv_In, target_pt_ct=pt_ct)
    elif isinstance(rgCrv_In, rg.ArcCurve):
        return samplePointsOnArcCurve(rgCrv_In, pt_ct)


def setEndConditions(rgCrv_ToMod, rgCrv_Target, bTanEnds):
    if bTanEnds:
        bSuccess_Start = rgCrv_ToMod.SetEndCondition(
            bSetEnd=False,
            continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
            point=rgCrv_Target.PointAtStart,
            tangent=rgCrv_Target.TangentAtStart)
        bSuccess_End = rgCrv_ToMod.SetEndCondition(
            bSetEnd=True,
            continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
            point=rgCrv_Target.PointAtEnd,
            tangent=rgCrv_Target.TangentAtEnd)
    else:
        bSuccess_Start = rgCrv_ToMod.SetEndCondition(
            bSetEnd=False,
            continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Position,
            point=rgCrv_Target.PointAtStart,
            tangent=rg.Vector3d.Unset)
        bSuccess_End = rgCrv_ToMod.SetEndCondition(
            bSetEnd=True,
            continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Position,
            point=rgCrv_Target.PointAtEnd,
            tangent=rg.Vector3d.Unset)

    return bSuccess_Start or bSuccess_End


def createCurve(rgCrv_In, degree, pt_ct, bTanEnds=True):
    """
    Parameters:
        rgCrv
        degree
        pt_ct
        bPreserveTans: Only for non-polylines
    Returns:
    """

    points = getStartingCpLocations(rgCrv_In, pt_ct)
    if not points: return

    if points.Count == pt_ct:
        nc = rg.NurbsCurve.Create(periodic=False, degree=degree, points=points)
    else:
        nc = rg.NurbsCurve.Create(
            periodic=False,
            degree=degree-(pt_ct-points.Count),
            points=points)
        nc.IncreaseDegree(degree)

    if not bTanEnds:
        return nc

    if points.Count == 2:
        # Nothing else can be done with this curve.
        return nc

    if degree == 2 and points.Count == 3:
        print("TODO: Create parabolic degree 2.")

    if not setEndConditions(nc, rgCrv_In, bTanEnds):
        print("Failed setting end condition (continuity).")
        return

    return nc


def main():

    objrefs_In = getInput()
    if objrefs_In is None: return

    iDegree = Opts.values['iDegree']
    iPointCt = Opts.values['iPointCt']
    bPreserveEndTangentsForNonPolylines = Opts.values['bPreserveEndTangentsForNonPolylines']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    for objref_In in objrefs_In:
        rgCrv = objref_In.Curve()
        rc_Res = createCurve(
            rgCrv,
            degree=iDegree,
            pt_ct=iPointCt,
            bTanEnds=bPreserveEndTangentsForNonPolylines)
        if rc_Res:
            sc.doc.Objects.AddCurve(rc_Res)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()