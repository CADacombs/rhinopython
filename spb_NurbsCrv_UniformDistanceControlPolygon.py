"""
This script is an alternative to _Smooth and _Fair but for a different purpose
or when those commands produce less desirable results.
It sets control point locations of NURBS curves to uniform distances
while preserving end clamp continuities per input.
Curve deviation is not limited as it is with _Fair.
If curve (to which the points are to be blended) curls toward itself, as in the
case of some S-curves, the resultant distances will likely be less uniform.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
230101-02: Created.

TODO: Produce more uniform spacing when points are blended to S-curves.
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


    key = 'bCurvedNotLinear'; keys.append(key)
    values[key] = True
    names[key] = 'Shape'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='Linear', onValue='Curved')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEqualEndClamps'; keys.append(key)
    values[key] = True
    names[key] = 'EqualEndClampCurvatureOrder'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClampDerivativeOrder'; keys.append(key)
    values[key] = 1
    names[key] = 'Order'
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClampStartDerivativeOrder'; keys.append(key)
    values[key] = 1
    names[key] = 'StartOrder'
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClampEndDerivativeOrder'; keys.append(key)
    values[key] = 1
    names[key] = 'EndOrder'
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bModify_NotAdd'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='Add', onValue='Modify')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'iClampDerivativeOrder':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = 0

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

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve


    def customGeometryFilter(rdObj, rgObj, compIdx):
        return isinstance(rgObj, rg.NurbsCurve)

    go.SetCustomGeometryFilter(customGeometryFilter)


    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    #go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False

    print("For ClampOrder(s), value represents maximum index"
          " (from relative end) of control points whose locations are not modified."
          "\nFor example, 2 means that first 3 controls point locations are preserved.")

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bCurvedNotLinear')
        addOption('bEqualEndClamps')
        if Opts.values['bEqualEndClamps']:
            addOption('iClampDerivativeOrder')
        else:
            addOption('iClampStartDerivativeOrder')
            addOption('iClampEndDerivativeOrder')
        addOption('bModify_NotAdd')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

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
            if Opts.values['bEqualEndClamps']:
                key = 'iClampDerivativeOrder'
            else:
                key = 'iClampStartDerivativeOrder'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getPreselectedCurves():
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    if gObjs_Preselected:
        gCrvs_Preselected = []
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
            if rdRhinoObject.Id in gObjs_Preselected:
                if rdRhinoObject.ObjectType == rd.ObjectType.Curve:
                    gCrvs_Preselected.append(rdRhinoObject.Id)
        if gCrvs_Preselected:
            if Opts.values['bEcho']:
                s  = "({} curves".format(len(gCrvs_Preselected))
                s += " were preselected.)"
                print(s)
            return tuple(gCrvs_Preselected)


def _createBlendCurve(nc, iClampStartDerivativeOrder, iClampEndDerivativeOrder):
    if iClampStartDerivativeOrder == 0:
        lcA = rg.LineCurve(
            rg.Point3d.Origin,
            nc.Points[iClampStartDerivativeOrder].Location)
        continuity0 = rg.BlendContinuity.Position
    else:
        lcA = rg.LineCurve(
            nc.Points[iClampStartDerivativeOrder-1].Location,
            nc.Points[iClampStartDerivativeOrder].Location)
        continuity0 = rg.BlendContinuity.Tangency

    if iClampEndDerivativeOrder == 0:
        lcB = rg.LineCurve(
            nc.Points[nc.Points.Count-1-iClampEndDerivativeOrder].Location,
            rg.Point3d.Origin)
        continuity1 = rg.BlendContinuity.Position
    else:
        lcB = rg.LineCurve(
            nc.Points[nc.Points.Count-1-iClampEndDerivativeOrder].Location,
            nc.Points[nc.Points.Count-1-iClampEndDerivativeOrder+1].Location)
        continuity1 = rg.BlendContinuity.Tangency

    return rg.Curve.CreateBlendCurve(
        curve0=lcA,
        t0=lcA.Domain.T1,
        reverse0=False,
        continuity0=continuity0,
        curve1=lcB,
        t1=lcB.Domain.T0,
        reverse1=True,
        continuity1=continuity1
        )


def _divideCurveEquidistantPerCount(rgCurve, point_ct, epsilon=1e-6, recurse=False, bDebug=False):
    """
    Parameters:
        rgCurve: rg.Curve
        point_ct: int(Target count of output points)
        epsilon: float(Accuracy of binary search.  Due to accuracy limits of 
            DivideEquidistant, values smaller than 1e-6 may be pointless.)
        recurse: bool(Whether this function call was called from itself
            If so, (2*point_ct)-1 points were passed for point_ct and half distance is returned
            Note that .)
    Returns:
        list(rg.Point3d)
    
    How it works:
        Finds distances on each side of target, then binary searches the distance.
        If last point doesn't match end of curve, then function calls itself at
        a higher point count.  Points returned are not exactly uniformly spaced.
    """


    def findEnds():
        
        dist_L = None
        dist_M = 1.0
        dist_H = None
        
        pts = rgCurve.DivideEquidistant(dist_M)
        iCt_M = len(pts)
        
        if len(pts) < iCt_Target:
            dist_H = dist_M
            multiplier = 0.5
        else:
            dist_L = dist_M
            multiplier = 2.0
        
        
        while True:
            sc.escape_test()
            
            dist_M *= multiplier
            
            pts = rgCurve.DivideEquidistant(dist_M)
            iCt_M = len(pts)
            
            if dist_L is None:
                if len(pts) >= iCt_Target:
                    return dist_M, dist_H
            else:
                if len(pts) < iCt_Target:
                    return dist_L, dist_M


    iCt_Target = point_ct

    dist_L = None
    dist_M = 1.0
    dist_H = None

    dist_L, dist_H = findEnds()


    iW = 0

    while True:
        sc.escape_test()
        
        if iW == 1000:
            for pt in pts:
                sc.doc.Objects.AddPoint(pt)
            raise Exception("Solution not found after 1000 iterations.  Points have been added.")

        dist_M = (dist_L + dist_H)/2.0
        
        if (dist_M - dist_L) <= epsilon or (dist_H - dist_M) <= epsilon:
            pts = rgCurve.DivideEquidistant(dist_L) # Because result from dist_M may be too long.

            dist_LastPtToCrvEnd = pts[-1].DistanceTo(rgCurve.PointAtEnd)
            if dist_LastPtToCrvEnd > epsilon:
                return _divideCurveEquidistantPerCount(
                    rgCurve,
                    2*iCt_Target - 1,
                    epsilon,
                    recurse=True)[::2]
            return pts
        
        pts = rgCurve.DivideEquidistant(dist_M) # Start point is included.
        iCt_M = len(pts)
        
        if len(pts) < iCt_Target:
            dist_H = dist_M
        else:
            dist_L = dist_M
        
        iW += 1


def modifyControlPolygon_Curved(nc, iClampStartDerivativeOrder, iClampEndDerivativeOrder, bDebug=False):
    """
    """

    if bDebug:
        sEval = "nc.Points.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "nc.Degree"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "iClampStartDerivativeOrder"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "iClampEndDerivativeOrder"; print("{}: {}".format(sEval, eval(sEval)))


    if nc.Points.Count <= (iClampStartDerivativeOrder + iClampEndDerivativeOrder + 2):
        return False, "Clamping continuity is too high for curve's point count."


    nc_Blend = _createBlendCurve(nc, iClampStartDerivativeOrder, iClampEndDerivativeOrder)

    if bDebug:
        sc.doc.Objects.AddCurve(nc_Blend)

    point_ct = nc.Points.Count - iClampEndDerivativeOrder - iClampStartDerivativeOrder

    pts_OnBlend = _divideCurveEquidistantPerCount(nc_Blend, point_ct)

    if len(pts_OnBlend) != point_ct:
        raise Exception("pts_OnBlend != point_ct")


    bModified = False

    for iPt, pt in enumerate(pts_OnBlend[1:-1]):
        if bDebug: sEval = "iPt"; print("{}: {}".format(sEval, eval(sEval)))

        if bDebug:
            sc.doc.Objects.AddPoint(pt)

        idxCp = iClampStartDerivativeOrder + 1 + iPt

        pt_Old = nc.Points[idxCp].Location
        dist = pt.DistanceTo(pt_Old)

        if dist <= Rhino.RhinoMath.ZeroTolerance:
            continue

        nc.Points.SetPoint(
            index=idxCp,
            point=pt,
            weight=nc.Points.GetWeight(idxCp))

        bModified = True

    if not bModified:
        return False, "All points are already at target locations."


    return True, None


def modifyControlPolygon_Linear(nc, iClampStartDerivativeOrder, iClampEndDerivativeOrder, bDebug=False):
    """
    """

    if bDebug:
        sEval = "nc.Points.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "nc.Degree"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "iClampStartDerivativeOrder"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "iClampEndDerivativeOrder"; print("{}: {}".format(sEval, eval(sEval)))


    if nc.Points.Count <= (iClampStartDerivativeOrder + iClampEndDerivativeOrder + 2):
        return False, "Clamping continuity is too high for curve's point count."

    pt_Start = nc.Points[iClampStartDerivativeOrder].Location
    pt_End = nc.Points[nc.Points.Count - iClampEndDerivativeOrder - 1].Location

    point_ct = nc.Points.Count - iClampEndDerivativeOrder - iClampStartDerivativeOrder

    pts_Uniform = [pt_Start]
    for iPt in range(1, point_ct-1):
        perunum = float(iPt)/float(point_ct-1)
        pt = (1.0 - perunum)*pt_Start + perunum*pt_End
        pts_Uniform.append(pt)
    pts_Uniform.append(pt_End)


    bModified = False

    for iPt, pt in enumerate(pts_Uniform[1:-1]):
        if bDebug: sEval = "iPt"; print("{}: {}".format(sEval, eval(sEval)))

        if bDebug:
            sc.doc.Objects.AddPoint(pt)

        idxCp = iClampStartDerivativeOrder + 1 + iPt

        pt_Old = nc.Points[idxCp].Location
        dist = pt.DistanceTo(pt_Old)

        if dist <= Rhino.RhinoMath.ZeroTolerance:
            continue

        nc.Points.SetPoint(
            index=idxCp,
            point=pt,
            weight=nc.Points.GetWeight(idxCp))

        bModified = True

    if not bModified:
        return False, "All points are already at target locations."


    return True, None


def main():

    gCs_Preselected = _getPreselectedCurves()

    objrefs = getInput()
    if objrefs is None: return

    bCurvedNotLinear = Opts.values['bCurvedNotLinear']
    if Opts.values['bEqualEndClamps']:
        iClampStartDerivativeOrder = Opts.values['iClampDerivativeOrder']
        iClampEndDerivativeOrder = Opts.values['iClampDerivativeOrder']
    else:
        iClampStartDerivativeOrder = Opts.values['iClampStartDerivativeOrder']
        iClampEndDerivativeOrder = Opts.values['iClampEndDerivativeOrder']
    bModify_NotAdd = Opts.values['bModify_NotAdd']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    gCs_Modified = []
    gCs_Added = []
    sLogs = []
    

    bUseCurvedFunction = (
        bCurvedNotLinear and
        iClampStartDerivativeOrder != 0 and
        iClampEndDerivativeOrder !=0)

    for objref in objrefs:
        rdC_In = objref.Object()

        nc_In = rdC_In.Geometry

        if nc_In.IsLinear():
            print("Linear curve skipped.")
            continue

        if bModify_NotAdd:
            nc = nc_In

            if bUseCurvedFunction:
                bModifiedPts, sLog = modifyControlPolygon_Curved(
                    nc=nc,
                    iClampStartDerivativeOrder=iClampStartDerivativeOrder,
                    iClampEndDerivativeOrder=iClampEndDerivativeOrder,
                    bDebug=bDebug)
            else:
                bModifiedPts, sLog = modifyControlPolygon_Linear(
                    nc=nc,
                    iClampStartDerivativeOrder=iClampStartDerivativeOrder,
                    iClampEndDerivativeOrder=iClampEndDerivativeOrder,
                    bDebug=bDebug)

            if sLog is not None:
                sLogs.append(sLog)

            if not bModifiedPts:
                continue

            if rdC_In.CommitChanges():
                gCs_Modified.append(rdC_In.Id)

        else:
            nc = nc_In.DuplicateCurve()

            if bUseCurvedFunction:
                bModifiedPts, sLog = modifyControlPolygon_Curved(
                    nc=nc,
                    iClampStartDerivativeOrder=iClampStartDerivativeOrder,
                    iClampEndDerivativeOrder=iClampEndDerivativeOrder,
                    bDebug=bDebug)
            else:
                bModifiedPts, sLog = modifyControlPolygon_Linear(
                    nc=nc,
                    iClampStartDerivativeOrder=iClampStartDerivativeOrder,
                    iClampEndDerivativeOrder=iClampEndDerivativeOrder,
                    bDebug=bDebug)

            if sLog is not None:
                sLogs.append(sLog)

            if not bModifiedPts:
                continue

            gC_Added = sc.doc.Objects.AddCurve(nc)
            if gC_Added != gC_Added.Empty:
                gCs_Added.append(gC_Added)


    for sLog in set(sLogs):
        print("[{}]: {}".format(sLogs.count(sLog), sLog))

    if bModify_NotAdd:
        print("{} curves modified.".format(len(gCs_Modified)))
    else:
        print("{} curves modified.".format(len(gCs_Added)))


    if gCs_Preselected:
        [sc.doc.Objects.Select(objectId=_) for _ in gCs_Preselected]
    else:
        if gCs_Modified:
            [sc.doc.Objects.Select(objectId=_) for _ in gCs_Modified]
        if gCs_Added:
            [sc.doc.Objects.Select(objectId=_) for _ in gCs_Added]

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()