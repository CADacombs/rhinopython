"""
This script is an alternative to _CrvDeviation.

'0' for LwrLimit, UprLimit, or MaxMinDistToRegard will disable the option.

MaxMinDistToRegard is the lowest value, over which the distances will be ignored.


Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

#! python 2

"""
170927: Created.
...
171204: Now polycurves and polylines are "exploded" into segments.  If they are not,
        Curve.GetDistancesBetweenCurves sometimes fails or produces erroneous results.
        Curves from curve set A are not selectable for set B even though they previously were removed anyway.
...
240903-10:  Added optional 2 curve input routine. Removed some functions. Refactored.
            Add display conduit similar to _CrvDeviation for when Mode=Abs.
            Trialing a UX different than _CrvDeviation for when Mode=Abs,
            that gives the user the option to skip leaving marks.

TODO:
    Clean code in spb_GDBCs_1Way that eliminates false positives at curve ends.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System.Drawing import Color


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bInputSets'; keys.append(key)
    values[key] = False
    names[key] = 'Input2'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Crvs', 'Sets')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExplodeCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fLocAlongCrvTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    # Using ModelUnitSystem in case sc.doc.Name is None.
    stickyKeys[key] = '{}({})({})({})'.format(key, __file__, sc.doc.Name, sc.doc.ModelUnitSystem)

    key = 'bOnlyPerp'; keys.append(key)
    values[key] = True
    names[key] = 'ClosestPtType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Any', 'OnlyPerp')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDist_max_to_regard'; keys.append(key)
    names[key] = 'MaxDistToRegard'
    values[key] = 200.0 * sc.doc.ModelAbsoluteTolerance
    #if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
    #    values[key] = 0.125
    #elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
    #    values[key] = 3.0
    #else:
    #    values[key] = 1000.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bLimitMode'; keys.append(key)
    values[key] = False
    names[key] = 'Mode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Abs', 'Limit')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fUprLimit'; keys.append(key)
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 0.026
    elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
        values[key] = 3.6
    else:
        values[key] = 1000.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fLwrLimit'; keys.append(key)
    values[key] = (
        0.014 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches
        else (
            2.4 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters
            else 10.0*sc.doc.ModelAbsoluteTolerance)
        )
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMarkMax'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMarkMin'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bVerifyAddMarks'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddLine'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotDecPlaces'; keys.append(key)
    values[key] = sc.doc.ModelDistanceDisplayPrecision - 1
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iDotFontHt'; keys.append(key)
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
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

        if key == 'fLocAlongCrvTol':
            if cls.riOpts[key].CurrentValue <= 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue <= max((1e-6, 0.001*sc.doc.ModelAbsoluteTolerance)):
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

        if key == 'fUprLimit':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]

            if (cls.values[key] > 0.0) and (cls.values[key] > cls.values['fDist_max_to_regard']):
                cls.values['fDist_max_to_regard'] = cls.riOpts['fDist_max_to_regard'].CurrentValue = cls.values[key]
                sc.sticky[cls.stickyKeys['fDist_max_to_regard']] = cls.values['fDist_max_to_regard']
            return

        if key == 'fLwrLimit':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
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


def getPreselectedCurves():
    gObjs_Preselected = []
    for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False):
        gObjs_Preselected.append(rdObj.Id)
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
        if len(gCrvs_Preselected) == 2:
            if Opts.values['bEcho']:
                s  = "({} curves".format(len(gCrvs_Preselected))
                s += " were preselected and will thus be the selection set.)"
                print(s)
            return tuple(gCrvs_Preselected)


def _addCommonOptions(go):

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    addOption('bInputSets')
    addOption('bExplodeCrvs')
    addOption('fLocAlongCrvTol')
    addOption('bOnlyPerp')
    addOption('fDist_max_to_regard')
    addOption('bLimitMode')
    if Opts.values['bLimitMode']:
        addOption('fUprLimit')
        addOption('fLwrLimit')
    else:
        addOption('bMarkMax')
        addOption('bMarkMin')
        if Opts.values['bMarkMax'] or Opts.values['bMarkMin']:
            addOption('bVerifyAddMarks')
    if (Opts.values['bLimitMode'] or Opts.values['bMarkMax'] or Opts.values['bMarkMin']):
        if not Opts.values['bAddLine'] and not Opts.values['bAddDot']:
            Opts.riOpts['bAddLine'].CurrentValue = True
            Opts.setValue('bAddLine')
            Opts.riOpts['bAddDot'].CurrentValue = True
            Opts.setValue('bAddDot')
        addOption('bAddLine')
        addOption('bAddDot')
        if Opts.values['bAddDot']:
            addOption('iDotDecPlaces')
            addOption('iDotFontHt')
    addOption('bEcho')
    addOption('bDebug')

    return idxs_Opts


def getInput_2Crvs():
    """
    Get 2 curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves")

    go.GeometryFilter = rd.ObjectType.Curve


    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.GroupSelect = True
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False

    idxs_Opts = {}
    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)


    while True:
        if Opts.values['bInputSets']:
            go.Dispose()
            sc.doc.Objects.UnselectAll()
            sc.doc.Views.Redraw()
            return getInput_2Sets()


        go.ClearCommandOptions()
        idxs_Opts.clear()
        idxs_Opts.update(_addCommonOptions(go))

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        # Use bPreselectedObjsChecked so that only selected objects before the
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
            return [objrefs[0].Curve()], [objrefs[1].Curve()]

        # An option was selected or a number was entered.

        if res == ri.GetResult.Number:
            key = 'fLocAlongCrvTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_2Sets():
    """
    Get 2 sets of curves with optional input.
    """

    rgCrvs_lists = [[],[]]

    go = ri.Custom.GetObject()
    sCommandPromptAdd = 'A', 'B'

    go.GeometryFilter = rd.ObjectType.Curve


    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.GroupSelect = True
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False
    
    idxs_Opts = {}
    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    # Get input for each curve set using a loop.
    for iCrvSet in 0,1:

        go.SetCommandPrompt("Select curve set {}".format(sCommandPromptAdd[iCrvSet]))

        while True:
            if not Opts.values['bInputSets']:
                go.Dispose()
                sc.doc.Objects.UnselectAll()
                sc.doc.Views.Redraw()
                return getInput_2Crvs()


            go.ClearCommandOptions()
            idxs_Opts.clear()
            idxs_Opts.update(_addCommonOptions(go))

            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
            
            # Use bPreselectedObjsChecked so that only selected objects before the
            # first call to go.GetMultiple is considered.
            if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
                bPreselectedObjsChecked = True
                go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
                continue

            if res == ri.GetResult.Cancel:
                go.Dispose()
                return

            if res == ri.GetResult.Object:
                break

            # An option was selected or a number was entered.

            if res == ri.GetResult.Number:
                key = 'fLocAlongCrvTol'
                Opts.riOpts[key].CurrentValue = go.Number()
                Opts.setValue(key)
                continue

            for key in idxs_Opts:
                if go.Option().Index == idxs_Opts[key]:
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)
                    break

        sc.doc.Objects.UnselectAll()
        sc.doc.Views.Redraw()

        if go.ObjectCount == 0: return
        
        rgCrvs_lists[iCrvSet] = [o.Geometry() for o in go.Objects()]
        
        if iCrvSet == 0:
            gCrvsA = [o.ObjectId for o in go.Objects()]
        
        # Custom geometry filter to only allow selection of curves not in first set.
        def curvesNotIn1stSetGeomFilter(rdObj, geom, compIdx):
            # TODO: Fix this.
            
            # Wires.
            if (
                compIdx.Index == -1 and
                not isinstance(geom, rg.BrepEdge) and
                rdObj.Id not in gCrvsA
            ):
                return True
            if isinstance(geom, rg.BrepEdge):
                return True
            return not rdObj.Id in gCrvsA
        go.SetCustomGeometryFilter(curvesNotIn1stSetGeomFilter)    
    
    go.Dispose()
    
    # Remove first set of curves from second.
    rgCrvs_lists[1] = rgCrvs_lists[1][len(rgCrvs_lists[0]):]
    
    if len(rgCrvs_lists[1]) == 0: return # Second set of objects were not selected.
    
    return tuple(
        ([rgCrvs_lists[0]]) +
        ([rgCrvs_lists[1]]) +
        [Opts.values[key] for key in Opts.keys])


def isMaxClosestDistBtwn2CrvsWithinTol(rgCrv_A, rgCrv_B, tolerance):
    """
    Alternative to Curve.GetDistancesBetweenCurves for better results when
    curves contain loops, etc.

    Returns:
        False (If not within tolerance parameter)
        float(Largest deviation found)
    """


    def isOutsideOfTolerance(rgC_Cat, rgC_Dog, ts_Cat):

        for iT_Cat in xrange(len(ts_Cat)):

            t_Cat = ts_Cat[iT_Cat]

            pt_Cat = rgC_Cat.PointAt(t_Cat)

            bSuccess, t_Dog = rgC_Dog.ClosestPoint(pt_Cat)

            if not bSuccess:
                raise ValueError("Closest point could not be calculated.")

            pt_Dog = rgC_Dog.PointAt(t_Dog)

            dist = pt_Cat.DistanceTo(pt_Dog)

            if dist > tolerance:
                return True

            fDevs.append(dist)

        return False

    fDivLength = 10.0*sc.doc.ModelAbsoluteTolerance


    fDevs = []

    # First, check span ends.
    ts_A = []
    for iSpan in range(rgCrv_A.SpanCount):
        spanDomain = rgCrv_A.SpanDomain(iSpan)
        ts_A.append(spanDomain.T0)
    if spanDomain.T1 not in ts_A:
        ts_A.append(spanDomain.T1)

    if isOutsideOfTolerance(rgCrv_A, rgCrv_B, ts_A):
        return False

    ts_B = []
    for iSpan in range(rgCrv_B.SpanCount):
        spanDomain = rgCrv_B.SpanDomain(iSpan)
        ts_B.append(spanDomain.T0)
    if spanDomain.T1 not in ts_B:
        ts_B.append(spanDomain.T1)

    if isOutsideOfTolerance(rgCrv_B, rgCrv_A, ts_B):
        return False


    for M in 1000.0, 10.0:
        fDivLength = M * sc.doc.ModelAbsoluteTolerance

        ts_A = []

        rc = rgCrv_A.DivideByLength(
            segmentLength=fDivLength,
            includeEnds=True)
        if rc:
            ts_A = rc
            if isOutsideOfTolerance(rgCrv_A, rgCrv_B, ts_A):
                return False

        ts_B = []

        rc = rgCrv_B.DivideByLength(
            segmentLength=fDivLength,
            includeEnds=True)
        if rc:
            ts_B = rc
            if isOutsideOfTolerance(rgCrv_B, rgCrv_A, ts_B):
                return False

    return max(fDevs)


def spb_GDBCs_1Way(curve_TestPts, curve_ClosestPt, fLocAlongCrvTol, bOnlyPerp=True, bDebug=False):
    """
    Alternative to Curve.GetDistancesBetweenCurves for more accuracy.

    Parameters:
        curve_TestPts: rg.Curve that will be divided to obtain the testPoints for ClosestPoint.
        curve_ClosestPt: rg.Curve that is the object of the ClosestPoint call.
        segmentLength: float : Division length of curve_TestPts to obtain some testPoints for ClosestPoint.
        bDebug: bool

    Returns:
        The same as Curve.GetDistancesBetweenCurves:
            bool: success
            float: maxDistance
            float: maxDistanceParameterA
            float: maxDistanceParameterB
            float: minDistance
            float: minDistanceParameterA
            float: minDistanceParameterB
    """

    segmentLength = 1000.0 * fLocAlongCrvTol

    if bDebug:
        sEval = "curve_TestPts.Domain.T0"; print(sEval, '=', eval(sEval))
        sEval = "curve_ClosestPt.Domain.T1"; print(sEval, '=', eval(sEval))
        sEval = "segmentLength"; print(sEval, '=', eval(sEval))


    def generate_list_of_curve_parameters_for_ClosestPoint(curve, segmentLength):
        if bDebug: print("generate_list_of_curve_parameters_for_ClosestPoint")

        ts_Out = []

        rc = curve.DivideByLength(
            segmentLength=segmentLength,
            includeEnds=True)
        if rc:
            ts_Out.extend(rc)

        # For open curve, add the ends of the curve.
        if not curve.IsClosed:
            # DivideByLength doesn't add the T1 segment
            # even when includeEnds == True.
            # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
            # shows the parameter labeled as 'includeStart'.
            if curve.Domain.T1 not in ts_Out:
                ts_Out.append(curve.Domain.T1)

        # Add parameters of all knots at full multiplicity.
        nc_Temp = curve.ToNurbsCurve()
        iK = 0
        while iK < nc_Temp.Knots.Count:
            sc.escape_test()
            m = nc_Temp.Knots.KnotMultiplicity(iK)
            if m == nc_Temp.Degree:
                k = nc_Temp.Knots[iK]
                if k not in ts_Out:
                    ts_Out.append(k)
            iK += m
        nc_Temp.Dispose()

        if len(ts_Out) != len(set(ts_Out)):
            sEval = "len(ts_Out)"; print(sEval, '=', eval(sEval))
            sEval = "len(set(ts_Out))"; print(sEval, '=', eval(sEval))
            raise ValueError("Duplicate parameters?  Check getClosestDistsBtwn2Crvs.")

        ts_Out.sort()

        if bDebug:
            sEval = "ts_Out[:10]"; print(sEval, '=', eval(sEval))
            sEval = "ts_Out[-10:]"; print(sEval, '=', eval(sEval))

        return ts_Out


    ts_A_FullCrv = generate_list_of_curve_parameters_for_ClosestPoint(curve_TestPts, segmentLength)

    if bDebug:
        sEval = "len(ts_A_FullCrv)"; print(sEval, '=', eval(sEval))


    def calc_parameters_and_distances(ts_A_In, curveA, curveB, bOnlyPerp):
        if bDebug: print("calc_parameters_and_distances")

        ts_B = []
        dists_per_ts_A_Out = []

        rads_90degs = Rhino.RhinoMath.ToRadians(90.0)

        for i_t_A, t_A in enumerate(ts_A_In):

            pt_A = curveA.PointAt(t_A)

            bSuccess, t_B = curveB.ClosestPoint(pt_A)

            if not bSuccess:
                raise ValueError("Closest point could not be calculated.")

            pt_B = curveB.PointAt(t_B)

            v_tan_B = curveB.TangentAt(t_B)
            v_dist = pt_A - pt_B

            if not bOnlyPerp:
                dists_per_ts_A_Out.append(pt_A.DistanceTo(pt_B))
                ts_B.append(t_B)
                continue

            # bOnlyPerp == True

            if v_dist.IsTiny():
                dists_per_ts_A_Out.append(pt_A.DistanceTo(pt_B))
                ts_B.append(t_B)
                continue

            angle_between = rg.Vector3d.VectorAngle(v_tan_B, v_dist)
            angle_from90 = abs(rads_90degs - angle_between)

            #if bDebug:
            #    sEval = "pt_A"; print(sEval, '=', eval(sEval))
            #    sEval = "pt_B"; print(sEval, '=', eval(sEval))
            #    sEval = "v_tan_B"; print(sEval, '=', eval(sEval))
            #    sEval = "v_dist"; print(sEval, '=', eval(sEval))
            #    sEval = "Rhino.RhinoMath.ToDegrees(angle_between)"; print(sEval, '=', eval(sEval))
            #    sEval = "Rhino.RhinoMath.ToDegrees(angle_from90)"; print(sEval, '=', eval(sEval))
            #    sEval = "v_dist.IsTiny()"; print(sEval, '=', eval(sEval))
                #if not v_dist.IsTiny():
                #sc.doc.Objects.AddPoint(pt_A)
                #sc.doc.Objects.AddPoint(pt_B)

            if angle_from90 > sc.doc.ModelAngleToleranceRadians:
                dists_per_ts_A_Out.append(None)
            else:
                dists_per_ts_A_Out.append(pt_A.DistanceTo(pt_B))

            ts_B.append(t_B)

        if bDebug:
            sEval = "len(ts_B)"; print(sEval, '=', eval(sEval))
            sEval = "len(dists_per_ts_A_Out)"; print(sEval, '=', eval(sEval))

        return ts_B, dists_per_ts_A_Out


    ts_B, dists_per_ts_A = calc_parameters_and_distances(
        ts_A_FullCrv,
        curve_TestPts,
        curve_ClosestPt,
        bOnlyPerp=bOnlyPerp)
    if bDebug:
        sEval = "len(dists_per_ts_A)"; print(sEval, '=', eval(sEval))
        sEval = "dists_per_ts_A[:10]"; print(sEval, '=', eval(sEval))
        sEval = "dists_per_ts_A[-10:]"; print(sEval, '=', eval(sEval))

    if all(d is None for d in dists_per_ts_A):
        return False, [], [], [], [], [], []

    # Get max distance.
    dist_Max = max(dists_per_ts_A)
    idx_MaxDist = dists_per_ts_A.index(dist_Max)
    t_A_MaxDist = ts_A_FullCrv[idx_MaxDist]
    t_B_MaxDist = ts_B[idx_MaxDist]

    if bDebug:
        sEval = "dist_Max"; print(sEval, '=', eval(sEval))
        sEval = "idx_MaxDist"; print(sEval, '=', eval(sEval))
        sEval = "t_A_MaxDist"; print(sEval, '=', eval(sEval))
        sEval = "t_B_MaxDist"; print(sEval, '=', eval(sEval))

    # Get min distance.
    dist_Min = min(d for d in dists_per_ts_A if d is not None)
    idx_MinDist = dists_per_ts_A.index(dist_Min)
    t_A_MinDist = ts_A_FullCrv[idx_MinDist]
    t_B_MinDist = ts_B[idx_MinDist]

    if bDebug:
        sEval = "dist_Min"; print(sEval, '=', eval(sEval))
        sEval = "idx_MinDist"; print(sEval, '=', eval(sEval))
        sEval = "t_A_MinDist"; print(sEval, '=', eval(sEval))
        sEval = "t_B_MinDist"; print(sEval, '=', eval(sEval))


    if bDebug:
        print("Iterate in smaller group of division points about the current winner",
              "to find a more accurate winner.")


    def findMoreAccurateWinner(ts_A_In, curveA, idx_Winner_In, segmentLength_In, fLocAlongCrvTol, bFindMax_NotMin):
        if bDebug: print("findMoreAccurateWinner")

        ts_A_WIP = ts_A_In[:]
        cA_WIP = curveA.Duplicate()
        idx_Winner_WIP = idx_Winner_In
        segmentLength = segmentLength_In

        t_A_Winner = ts_A_WIP[idx_Winner_In]

        while True:
            sc.escape_test()

            segmentLength *= 0.1

            if segmentLength < (fLocAlongCrvTol - 1e-6):
                break

            if bDebug: sEval = "segmentLength"; print(sEval, '=', eval(sEval))


            if idx_Winner_WIP > 0:
                t0 = ts_A_WIP[idx_Winner_WIP - 1]
            else:
                t0 = ts_A_WIP[0]

            if idx_Winner_WIP < (len(ts_A_WIP) - 1):
                t1 = ts_A_WIP[idx_Winner_WIP + 1]
            else:
                t1 = ts_A_WIP[len(ts_A_WIP) - 1]

            if bDebug:
                sEval = "t0"; print(sEval, '=', eval(sEval))
                sEval = "t1"; print(sEval, '=', eval(sEval))

            cA_WIP = cA_WIP.Trim(rg.Interval(t0, t1))
            #sc.doc.Objects.AddCurve(cA_WIP); sc.doc.Views.Redraw(); 1/0
            if bDebug: sEval = "cA_WIP.GetLength()"; print(sEval, '=', eval(sEval))

            ts_A_WIP = generate_list_of_curve_parameters_for_ClosestPoint(cA_WIP, segmentLength)
            if bDebug:
                sEval = "len(ts_A_WIP)"; print(sEval, '=', eval(sEval))
                sEval = "ts_A_WIP[:10]"; print(sEval, '=', eval(sEval))
                sEval = "ts_A_WIP[-10:]"; print(sEval, '=', eval(sEval))




            ts_B_WIP, dists_per_ts_A_WIP = calc_parameters_and_distances(
                ts_A_WIP,
                cA_WIP,
                curve_ClosestPt,
                bOnlyPerp=bOnlyPerp)

            # Get winning distance.
            if bFindMax_NotMin:
                dist_Winner = max(dists_per_ts_A_WIP)
            else:
                dist_Winner = min(d for d in dists_per_ts_A_WIP if d is not None)
            idx_Winner_WIP = dists_per_ts_A_WIP.index(dist_Winner)
            t_A_Winner = ts_A_WIP[idx_Winner_WIP]
            t_B_Winner = ts_B_WIP[idx_Winner_WIP]

            if bDebug:
                sEval = "dist_Winner"; print(sEval, '=', eval(sEval))
                sEval = "idx_Winner_WIP"; print(sEval, '=', eval(sEval))
                sEval = "t_A_Winner"; print(sEval, '=', eval(sEval))
                sEval = "t_B_Winner"; print(sEval, '=', eval(sEval))

        return dist_Winner, t_A_Winner, t_B_Winner


    rc = findMoreAccurateWinner(
        ts_A_In=ts_A_FullCrv,
        curveA=curve_TestPts,
        idx_Winner_In=idx_MaxDist,
        segmentLength_In=segmentLength,
        fLocAlongCrvTol=fLocAlongCrvTol,
        bFindMax_NotMin=True)

    dist_Max, t_A_MaxDist, t_B_MaxDist = rc


    if dist_Min > 0.0:
        rc = findMoreAccurateWinner(
            ts_A_In=ts_A_FullCrv,
            curveA=curve_TestPts,
            idx_Winner_In=idx_MinDist,
            segmentLength_In=segmentLength,
            fLocAlongCrvTol=fLocAlongCrvTol,
            bFindMax_NotMin=False)

        dist_Min, t_A_MinDist, t_B_MinDist = rc


    if (
        dist_Max is not None and
        dist_Max > 1e-6 and
        ((t_A_MaxDist - curve_TestPts.Domain.T0) <= 1e-6)
        or
        ((t_A_MaxDist - curve_TestPts.Domain.T1) <= 1e-6)
        ):
        v_dist = curve_TestPts.PointAt(t_A_MaxDist) - curve_ClosestPt.PointAt(t_B_MaxDist)

        if v_dist.IsTiny():
            pass
        else:
            rads_90degs = Rhino.RhinoMath.ToRadians(90.0)

            v_tan_A = curve_TestPts.TangentAt(t_A_MaxDist)

            angle_between = rg.Vector3d.VectorAngle(v_tan_A, v_dist)
            angle_from90 = abs(rads_90degs - angle_between)

            if angle_from90 > Rhino.RhinoMath.ToRadians(22.5):
                dist_Max = None
                t_A_MaxDist = None
                t_B_MaxDist = None


    if (
        dist_Min is not None and
        dist_Min > 1e-6 and
        ((t_A_MinDist - curve_TestPts.Domain.T0) <= 1e-6) or
        ((t_A_MinDist - curve_TestPts.Domain.T1) <= 1e-6)
        ):
        v_dist = curve_TestPts.PointAt(t_A_MinDist) - curve_ClosestPt.PointAt(t_B_MinDist)

        if v_dist.IsTiny():
            pass
        else:
            rads_90degs = Rhino.RhinoMath.ToRadians(90.0)

            v_tan_A = curve_TestPts.TangentAt(t_A_MinDist)

            angle_between = rg.Vector3d.VectorAngle(v_tan_A, v_dist)
            angle_from90 = abs(rads_90degs - angle_between)

            if angle_from90 > Rhino.RhinoMath.ToRadians(22.5):
                dist_Min = None
                t_A_MinDist = None
                t_B_MinDist = None



    return (
        True,
        dist_Max, 
        t_A_MaxDist,
        t_B_MaxDist,
        dist_Min,
        t_A_MinDist,
        t_B_MinDist,
        )


def spb_GDBCs_BothWays(curveA, curveB, fLocAlongCrvTol=None, bOnlyPerp=True, bDebug=False):
    """
    Alternative to Curve.GetDistancesBetweenCurves for more accurate results when
    curves contain loops, etc.
    Returns:
        The same as Curve.GetDistancesBetweenCurves:
            bool: success
            float: maxDistance
            float: maxDistanceParameterA
            float: maxDistanceParameterB
            float: minDistance
            float: minDistanceParameterA
            float: minDistanceParameterB
    """


    if fLocAlongCrvTol is None:
        fLocAlongCrvTol = 100.0*sc.doc.ModelAbsoluteTolerance
    if bDebug: sEval = "fLocAlongCrvTol"; print(sEval, '=', eval(sEval))

    Rhino.RhinoApp.Wait()
    rc = spb_GDBCs_1Way(
        curve_TestPts=curveA,
        curve_ClosestPt=curveB,
        fLocAlongCrvTol=fLocAlongCrvTol,
        bOnlyPerp=bOnlyPerp,
        bDebug=bDebug)
    (
        bSuccess_onB,
        dist_Max_ClosestPt_on_B,
        tA_MaxDist_ClosestPt_on_B,
        tB_MaxDist_ClosestPt_on_B,
        dist_Min_ClosestPt_on_B,
        tA_MinDist_ClosestPt_on_B,
        tB_MinDist_ClosestPt_on_B,
       ) = rc
    if bDebug: sEval = "rc"; print(sEval, '=', eval(sEval))

    Rhino.RhinoApp.Wait()

    # Notice that curveA and curveB are reversed.
    rc = spb_GDBCs_1Way(
        curve_TestPts=curveB,
        curve_ClosestPt=curveA,
        fLocAlongCrvTol=fLocAlongCrvTol,
        bOnlyPerp=bOnlyPerp,
        bDebug=bDebug)
    (
        bSuccess_onA,
        dist_Max_ClosestPt_on_A,
        tB_MaxDist_ClosestPt_on_A,
        tA_MaxDist_ClosestPt_on_A,
        dist_Min_ClosestPt_on_A,
        tB_MinDist_ClosestPt_on_A,
        tA_MinDist_ClosestPt_on_A,
       ) = rc
    if bDebug: sEval = "rc"; print(sEval, '=', eval(sEval))


    if not bSuccess_onA and not bSuccess_onB:
        return

    if not bSuccess_onA and bSuccess_onB:
        dist_Max = dist_Max_ClosestPt_on_B
        tA_Max_Out = tA_MaxDist_ClosestPt_on_B
        tB_Max_Out = tB_MaxDist_ClosestPt_on_B
        dist_Min = dist_Min_ClosestPt_on_B
        tA_Min_Out = tA_MinDist_ClosestPt_on_B
        tB_Min_Out = tB_MinDist_ClosestPt_on_B
        #sc.doc.Objects.AddCurve(curveA)
        #sc.doc.Objects.AddCurve(curveB)
        #sc.doc.Views.Redraw()
        #raise Exception("not bSuccess_onA and bSuccess_onB")

    elif bSuccess_onA and not bSuccess_onB:
        dist_Max = dist_Max_ClosestPt_on_A
        tA_Max_Out = tA_MaxDist_ClosestPt_on_A
        tB_Max_Out = tB_MaxDist_ClosestPt_on_A
        dist_Min = dist_Min_ClosestPt_on_A
        tA_Min_Out = tA_MinDist_ClosestPt_on_A
        tB_Min_Out = tB_MinDist_ClosestPt_on_A
        #sc.doc.Objects.AddCurve(curveA)
        #sc.doc.Objects.AddCurve(curveB)
        #sc.doc.Views.Redraw()
        #raise Exception("bSuccess_onA and not bSuccess_onB")

    else:
        if dist_Max_ClosestPt_on_B > dist_Max_ClosestPt_on_A:
            dist_Max = dist_Max_ClosestPt_on_B
            tA_Max_Out = tA_MaxDist_ClosestPt_on_B
            tB_Max_Out = tB_MaxDist_ClosestPt_on_B
            #ptA_Max = curveA.PointAt(tA_MaxDist_onB_from_pt_on_A)
            #ptB_Max = curveB.PointAt(tB_MaxDist_onB_from_pt_on_A)
        else:
            dist_Max = dist_Max_ClosestPt_on_A
            tA_Max_Out = tA_MaxDist_ClosestPt_on_A
            tB_Max_Out = tB_MaxDist_ClosestPt_on_A
            #ptA_Max = curveA.PointAt(tA_MaxDist_onA_from_pt_on_B)
            #ptB_Max = curveB.PointAt(tB_MaxDist_onA_from_pt_on_B)

        if dist_Min_ClosestPt_on_B < dist_Min_ClosestPt_on_A:
            dist_Min = dist_Min_ClosestPt_on_B
            tA_Min_Out = tA_MinDist_ClosestPt_on_B
            tB_Min_Out = tB_MinDist_ClosestPt_on_B
            #ptA_Min = curveA.PointAt(tA_MinDist_onB_from_pt_on_A)
            #ptB_Min = curveB.PointAt(tB_MinDist_onB_from_pt_on_A)
        else:
            dist_Min = dist_Min_ClosestPt_on_A
            tA_Min_Out = tA_MinDist_ClosestPt_on_A
            tB_Min_Out = tB_MinDist_ClosestPt_on_A
            #ptA_Min = curveA.PointAt(tA_MinDist_onA_from_pt_on_B)
            #ptB_Min = curveB.PointAt(tB_MinDist_onA_from_pt_on_B)


    # Check for intersection to replace any minimum distance values.
    if dist_Min > 1e-6:
        intersections = rg.Intersect.Intersection.CurveCurve(
            curveA,
            curveB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
            overlapTolerance=0.0)
        if intersections:
            dist_Min = 0.0
            tA_Min_Out = intersections[0].ParameterA
            tB_Min_Out = intersections[0].ParameterB
            #ptA_Min = intersections[0].PointA
            #ptB_Min = intersections[0].PointB


    return (
        dist_Max,
        tA_Max_Out,
        tB_Max_Out,
        dist_Min,
        tA_Min_Out,
        tB_Min_Out
        )


def getDevsBtwn2Sets(rgCrvs_SetA, rgCrvs_SetB, **kwargs):
    """
    Returns:
        fDists_Max,
        ptsA_Max,
        ptsB_Max,
        fDists_Min,
        ptsA_Min,
        ptsB_Min
        Unlike Curve.GetDistancesBetweenCurves, doesn't return a boolean value for success.
        In the case of fail, empty lists are the values of each variable.
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bExplodeCrvs = getOpt('bExplodeCrvs')
    fLocAlongCrvTol = getOpt('fLocAlongCrvTol')
    bOnlyPerp = getOpt('bOnlyPerp')
    fDist_max_to_regard = getOpt('fDist_max_to_regard')
    bDebug = getOpt('bDebug')


    try: rgCrvs_SetA = list(rgCrvs_SetA)
    except: rgCrvs_SetA = [rgCrvs_SetA]

    try: rgCrvs_SetB = list(rgCrvs_SetB)
    except: rgCrvs_SetB = [rgCrvs_SetB]

    fDists_Max = []
    ptsA_Max = []
    ptsB_Max = []

    fDists_Min = []
    ptsA_Min = []
    ptsB_Min = []


    # Explode parent curves into segments.


    if not bExplodeCrvs:
        rgCs_A = rgCrvs_SetA
        rgCs_B = rgCrvs_SetB
    else:
        rgCs_A = []
        for c in rgCrvs_SetA:
            rgCs = c.DuplicateSegments()
            if rgCs:
                rgCs_A.extend(rgCs)
            else:
                 # DuplicateSegments returns an empty array when
                 # curve is not composed of multiple segments.
                rgCs_A.append(c)

        rgCs_B = []
        for c in rgCrvs_SetB:
            rgCs = c.DuplicateSegments()
            if rgCs:
                rgCs_B.extend(rgCs)
                 # DuplicateSegments returns an empty array when
                 # curve is not composed of multiple segments.
            else:
                rgCs_B.append(c)


    # Due to a bug in V5 (SR12) bug
    # (fixed in V6 per https://mcneel.myjetbrains.com/youtrack/issue/RH-41877),
    # bad parameter values are returned from GetDistancesBetweenCurves for ArcCurves,
    # so convert any found to a NurbsCurve.
    if Rhino.RhinoApp.ExeVersion == 5:
        for i, cA in enumerate(rgCs_A):
            if isinstance(cA, rg.ArcCurve):
                rgCs_A[i] = cA.ToNurbsCurve()
        for i, cB in enumerate(rgCs_B):
            if isinstance(cB, rg.ArcCurve):
                rgCs_B[i] = cB.ToNurbsCurve()


    for i_cA, cA in enumerate(rgCs_A):
        Rhino.RhinoApp.Wait()
        Rhino.RhinoApp.SetCommandPrompt(
            "Working ... Curve A {} of {} on each of {} Curve Bs".format(
                i_cA, len(rgCs_A), len(rgCs_B)))

        for i_cB, cB in enumerate(rgCs_B):
            sc.escape_test()

            rc = spb_GDBCs_BothWays(
                cA,
                cB,
                fLocAlongCrvTol=fLocAlongCrvTol,
                bOnlyPerp=bOnlyPerp,
                bDebug=bDebug)

            if rc is None:
                if bDebug: print("No valid deviation found.")
                continue

            (
                fDist_Max,
                tA_Max,
                tB_Max,
                fDist_Min,
                tA_Min,
                tB_Min
               ) = rc

            if fDist_Min is not None:
                if (
                    fDist_max_to_regard and
                    fDist_Min > fDist_max_to_regard
                    ):
                        if bDebug: print("Skipped distance {}.".format(fDist_Min))
                else:
                    ptA_Min = None if tA_Min is None else cA.PointAt(tA_Min)
                    ptB_Min = None if tB_Min is None else cB.PointAt(tB_Min)
    
                    fDists_Min.append(fDist_Min)
                    ptsA_Min.append(ptA_Min)
                    ptsB_Min.append(ptB_Min)

            if fDist_Max is not None:
                if (
                    fDist_max_to_regard and
                    fDist_Max > fDist_max_to_regard
                    ):
                        if bDebug: print("Skipped distance {}.".format(fDist_Max))
                else:
                    ptA_Max = None if tA_Max is None else cA.PointAt(tA_Max)
                    ptB_Max = None if tB_Max is None else cB.PointAt(tB_Max)
    
                    fDists_Max.append(fDist_Max)
                    ptsA_Max.append(ptA_Max)
                    ptsB_Max.append(ptB_Max)


    return (
        fDists_Max,
        ptsA_Max,
        ptsB_Max,
        fDists_Min,
        ptsA_Min,
        ptsB_Min
        )


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.lines = []
        self.line_colors = []
        self.dots = []
        self.dot_colors = []
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        for line in self.lines:
            bbox = line.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

        for dot in self.dots:
            bbox = dot.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PreDrawObjects(self, drawEventArgs):

        for i, line in enumerate(self.lines):
            drawEventArgs.Display.DrawLine(
                line=line,
                color=self.line_colors[i],
                thickness=self.crv_thk)

        for i, dot in enumerate(self.dots):
            drawEventArgs.Display.DrawDot(
                worldPosition=dot.Point,
                text=dot.Text,
                dotColor=self.dot_colors[i],
                textColor=Color.White if self.dot_colors[i] == Color.Red else Color.Black)

            # TODO: Why doesn't this work?
            #drawEventArgs.Display.DrawDot(
            #    dot=dot,
            #    fillColor=self.dot_colors[i],
            #    dotColor=self.dot_colors[i],
            #    textColor=Color.White,
            #    borderColor=self.dot_colors[i])


def _addTextDot(text, pt, iDotFontHt, attrib):
    rgDot = rg.TextDot(text, pt)
    rgDot.FontHeight = iDotFontHt
    sc.doc.Objects.AddTextDot(rgDot, attrib)


def _addPreview(conduit, fDists_Max, ptsA_Max, ptsB_Max, fDists_Min, ptsA_Min, ptsB_Min, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAddLine = getOpt('bAddLine')
    bAddDot = getOpt('bAddDot')
    iDotDecPlaces = getOpt('iDotDecPlaces')
    iDotFontHt = getOpt('iDotFontHt')


    if not (bAddLine or bAddDot):
        return

    # Maximum deviation(s).
    for i, fDist_Max in enumerate(fDists_Max):
        ptA_Max = ptsA_Max[i]
        ptB_Max = ptsB_Max[i]
        if None in (fDist_Max, ptA_Max, ptB_Max):
            raise Exception("None in {}, {} ,{}".format(fDist_Max, ptA_Max, ptB_Max))
        if bAddLine:
            conduit.lines.append(rg.Line(ptA_Max, ptB_Max))
            conduit.line_colors.append(Color.Red)
        if bAddDot:
            conduit.dots.append(
                rg.TextDot(
                    text="{0:.{1}f}".format(fDist_Max, iDotDecPlaces),
                    location=(ptA_Max + ptB_Max) / 2.0))
            conduit.dot_colors.append(Color.Red)

    # Minimum deviation(s).
    for i, fDist_Min in enumerate(fDists_Min):
        ptA_Min = ptsA_Min[i]
        ptB_Min = ptsB_Min[i]
        if None in (fDist_Min, ptA_Min, ptB_Min):
            raise Exception("None in {}, {} ,{}".format(fDist_Max, ptA_Min, ptB_Min))
        if bAddLine and fDist_Min:
            conduit.lines.append(rg.Line(ptA_Min, ptB_Min))
            conduit.line_colors.append(Color.Lime)
        if bAddDot:
            conduit.dots.append(
                rg.TextDot(
                    text="{0:.{1}f}".format(fDist_Min, iDotDecPlaces),
                    location=(ptA_Min + ptB_Min) / 2.0))
            conduit.dot_colors.append(Color.Lime)


def _addAnnotation(fDists_Max, ptsA_Max, ptsB_Max, fDists_Min, ptsA_Min, ptsB_Min, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAddLine = getOpt('bAddLine')
    bAddDot = getOpt('bAddDot')
    iDotDecPlaces = getOpt('iDotDecPlaces')
    iDotFontHt = getOpt('iDotFontHt')


    if not (bAddLine or bAddDot):
        return

    attr_Red = rd.ObjectAttributes()
    attr_Red.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    #attr.ObjectDecoration = rd.ObjectDecoration.BothArrowhead
    attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr_Red.ObjectColor = attr_Red.ObjectColor.Red

    attr_Green = attr_Red.Duplicate()
    attr_Green.ObjectColor = attr_Green.ObjectColor.Lime

    # Maximum deviation(s).
    if fDists_Max:
        for i, fDist_Max in enumerate(fDists_Max):
            ptA_Max = ptsA_Max[i]
            ptB_Max = ptsB_Max[i]
            if None in (fDist_Max, ptA_Max, ptB_Max):
                raise Exception("None in {}, {} ,{}".format(fDist_Max, ptA_Max, ptB_Max))
            #if (fDist_max is not None
            #        and ptA_max is not None and ptB_max is not None
            #):
            if bAddLine:
                sc.doc.Objects.AddLine(ptA_Max, ptB_Max, attr_Red)
            if bAddDot:
                _addTextDot(
                    "{0:.{1}f}".format(fDist_Max, iDotDecPlaces),
                    (ptA_Max + ptB_Max) / 2.0,
                    iDotFontHt,
                    attr_Red)

    # Minimum deviation(s).
    if fDists_Min:
        for i, fDist_Min in enumerate(fDists_Min):
            ptA_Min = ptsA_Min[i]
            ptB_Min = ptsB_Min[i]
            if None in (fDist_Min, ptA_Min, ptB_Min):
                raise Exception("None in {}, {} ,{}".format(fDist_Min, ptA_Min, ptB_Min))
            #if (fDist_min is not None
            #        and ptA_min is not None and ptB_min is not None):
            if bAddLine and fDist_Min:
                sc.doc.Objects.AddLine(ptA_Min, ptB_Min, attr_Green)
            if bAddDot:
                _addTextDot(
                    "{0:.{1}f}".format(fDist_Min, iDotDecPlaces),
                    (ptA_Min + ptB_Min) / 2.0,
                    iDotFontHt,
                    attr_Green)


def main():

    sk_conduit = 'conduit({})'.format(__file__) # StickyKey
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
        sc.doc.Views.Redraw()


    rgCrvsA = []
    rgCrvsB = []

    rc = getPreselectedCurves()
    if rc:
        rgCrvsA = [sc.doc.Objects.FindId(rc[0]).Geometry]
        rgCrvsB = [sc.doc.Objects.FindId(rc[1]).Geometry]
    else:
        bInputSets = Opts.values['bInputSets']
        if bInputSets:
            rc = getInput_2Sets()
        else:
            rc = getInput_2Crvs()
        if rc is None: return
        #print(rc)
        rgCrvsA = rc[0]
        rgCrvsB = rc[1]

    bInputSets = Opts.values['bInputSets']
    bExplodeCrvs = Opts.values['bExplodeCrvs']
    fLocAlongCrvTol = Opts.values['fLocAlongCrvTol']
    bOnlyPerp = Opts.values['bOnlyPerp']
    fDist_max_to_regard = Opts.values['fDist_max_to_regard']
    bLimitMode = Opts.values['bLimitMode']
    fUprLimit = Opts.values['fUprLimit']
    fLwrLimit = Opts.values['fLwrLimit']
    bMarkMax = Opts.values['bMarkMax']
    bMarkMin = Opts.values['bMarkMin']
    bVerifyAddMarks = Opts.values['bVerifyAddMarks']
    bAddLine = Opts.values['bAddLine']
    bAddDot = Opts.values['bAddDot']
    iDotDecPlaces = Opts.values['iDotDecPlaces']
    iDotFontHt = Opts.values['iDotFontHt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if bDebug:
        sEval = "rgCrvsA"; print(sEval, '=', eval(sEval))
        sEval = "len(rgCrvsA)"; print(sEval, '=', eval(sEval))
        sEval = "rgCrvsB"; print(sEval, '=', eval(sEval))
        sEval = "len(rgCrvsB)"; print(sEval, '=', eval(sEval))


    rc = getDevsBtwn2Sets(
        rgCrvsA,
        rgCrvsB,
        bExplodeCrvs=bExplodeCrvs,
        fLocAlongCrvTol=fLocAlongCrvTol,
        bOnlyPerp=bOnlyPerp,
        fDist_max_to_regard=fDist_max_to_regard,
        bDebug=bDebug,
        )

    (
        fDists_Max,
        ptsA_Max,
        ptsB_Max,
        fDists_Min,
        ptsA_Min,
        ptsB_Min
        ) = rc



    if (
        all([d is None for d in fDists_Max])
        and
        all([d is None for d in fDists_Min])
        ):
        print("No results.")
        return


    fDist_Max_All = max(fDists_Max) if fDists_Max else None
    fDist_Min_All = min(fDists_Min) if fDists_Min else None

    if bDebug:
        sEval = "fDist_Max_All"; print(sEval, '=', eval(sEval))
        sEval = "fDist_Min_All"; print(sEval, '=', eval(sEval))


    if bEcho:
        s = ""
        if fDist_Min_All is not None:
            s += "Minimum deviation = {0:.{1}g}".format(
                    fDist_Min_All, 6)
        if bLimitMode and fDist_Min_All >= fLwrLimit:
            s += "  None are below {}.".format(fLwrLimit)
        print(s)

        s = ""
        if fDist_Max_All is not None:
            s += "Maximum deviation = {0:.{1}g}".format(
                    fDist_Max_All, 6)
        if bLimitMode and fDist_Max_All <= fUprLimit:
            s += "  None are above {}.".format(fUprLimit)

        print(s)



    if not bLimitMode and bVerifyAddMarks and (bMarkMax or bMarkMin):
        if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
            conduit = None
            #    conduit = DrawConduit()
            #    sc.sticky[sk_conduit] = conduit
            #if not ((sk_conduit in sc.sticky) and sc.sticky[sk_conduit]):
            #    conduit = DrawConduit()
            #    sc.sticky[sk_conduit] = conduit
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit


        if bMarkMax and bMarkMin:
            _addPreview(
                conduit,
                fDists_Max, ptsA_Max, ptsB_Max,
                fDists_Min, ptsA_Min, ptsB_Min,
                bAddLine=bAddLine,
                bAddDot=bAddDot,
                iDotDecPlaces=iDotDecPlaces,
                iDotFontHt=iDotFontHt,
                )
        elif bMarkMax:
            _addPreview(
                conduit,
                fDists_Max, ptsA_Max, ptsB_Max,
                [], [], [],
                bAddLine=bAddLine,
                bAddDot=bAddDot,
                iDotDecPlaces=iDotDecPlaces,
                iDotFontHt=iDotFontHt,
                )
        elif bMarkMin:
            _addPreview(
                conduit,
                [], [], [],
                fDists_Min, ptsA_Min, ptsB_Min,
                bAddLine=bAddLine,
                bAddDot=bAddDot,
                iDotDecPlaces=iDotDecPlaces,
                iDotFontHt=iDotFontHt,
                )
        else:
            raise Exception("What happened?")

        conduit.Enabled = True
        sc.doc.Views.Redraw()

        go = ri.Custom.GetOption()
        go.SetCommandPrompt("Press Enter to skip marks or")
        idx_Opt = go.AddOption("KeepMarks")

        res = go.Get()

        conduit.Enabled = False

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if not go.Option():
            return



    if bLimitMode:
        if (fDist_Min_All >= fLwrLimit) and (fDist_Max_All <= fUprLimit):
            return


        def filter_results_per_limits(
            fDists_Max_In,
            ptsA_Max_In,
            ptsB_Max_In,
            fDists_Min_In,
            ptsA_Min_In,
            ptsB_Min_In
            ):

            fDists_Max_Out = []
            ptsA_Max_Out = []
            ptsB_Max_Out = []
            fDists_Min_Out = []
            ptsA_Min_Out = []
            ptsB_Min_Out = []

            for i, fDist_Min in enumerate(fDists_Min_In):
                fDist_Max = fDists_Max_In[i]
                ptA_Max = ptsA_Max_In[i]
                ptB_Max = ptsB_Max_In[i]
                ptA_Min = ptsA_Min_In[i]
                ptB_Min = ptsB_Min_In[i]

                if fUprLimit and (fDist_Max > fUprLimit):
                    fDists_Max_Out.append(fDist_Max)
                    ptsA_Max_Out.append(ptA_Max)
                    ptsB_Max_Out.append(ptB_Max)

                if fLwrLimit and (fDist_Min < fLwrLimit):
                    fDists_Min_Out.append(fDist_Min)
                    ptsA_Min_Out.append(ptA_Min)
                    ptsB_Min_Out.append(ptB_Min)

            return (
                fDists_Max_Out,
                ptsA_Max_Out,
                ptsB_Max_Out,
                fDists_Min_Out,
                ptsA_Min_Out,
                ptsB_Min_Out
                )


        rc = filter_results_per_limits(*rc)

        (
            fDists_Max,
            ptsA_Max,
            ptsB_Max,
            fDists_Min,
            ptsA_Min,
            ptsB_Min
            ) = rc


        if not fDists_Max and not fDists_Min:
            if bEcho:
                print("None of the curves deviate outside of limits.")
            return



    if bLimitMode or (bMarkMax and bMarkMin):
        _addAnnotation(
            fDists_Max, ptsA_Max, ptsB_Max,
            fDists_Min, ptsA_Min, ptsB_Min,
            bAddLine=bAddLine,
            bAddDot=bAddDot,
            iDotDecPlaces=iDotDecPlaces,
            iDotFontHt=iDotFontHt,
            )
    elif bMarkMax:
        _addAnnotation(
            fDists_Max, ptsA_Max, ptsB_Max,
            [], [], [],
            bAddLine=bAddLine,
            bAddDot=bAddDot,
            iDotDecPlaces=iDotDecPlaces,
            iDotFontHt=iDotFontHt,
            )
    elif bMarkMin:
        _addAnnotation(
            [], [], [],
            fDists_Min, ptsA_Min, ptsB_Min,
            bAddLine=bAddLine,
            bAddDot=bAddDot,
            iDotDecPlaces=iDotDecPlaces,
            iDotFontHt=iDotFontHt,
            )
    else:
        return

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()