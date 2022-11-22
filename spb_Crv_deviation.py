"""
170927: Created.
171005: Bug fix: Line layer is now set to the current layer.
171009: Added workaround for bug in V5's GetDistancesBetweenCurves (see comments in getDistancesBetweenCurveSets).
171024: Added sc.escape_test() in main loop in getDistancesBetweenCurveSets.
171104: Now has option to add none, one, or all mininums found under a specified tolerance.
        Enabled GroupSelect.
171112: Now has option to add none, one, or all maximums found over a specified tolerance.
        Bug fix: Accumulation of minimums didn't work when fDist_Tol_Min==0.
171116: Moved maximum and minimum tolerance check from getDistancesBetweenCurveSets to main.
171204: Now polycurves and polylines are "exploded" into segments.  If they are not,
        Curve.GetDistancesBetweenCurves sometimes fails or produces erroneous results.
        Curves from curve set A are not selectable for set B even though they previously were removed anyway.
180202: Bug fix in addAnnotation.
180228: Default tolerance is now unique for an inch unit model.
180427: No longer adds arrowheads to lines.
181102: Changed name from crvDeviation to curveDeviation.  Changed option default.
181120: Fixed bug in sc.sticky routine.
190829: Now, ignores maximum results where the distances are not perpendicular to either curve.
        Testing having fDist_maxToAccept disabled.
190831: Renabled fDist_maxToAccept.
200216: Added comments to a function.
200409-10: Refactored.  Added some functions.  WIP: Added bMustBePerpAtOpenEnds.
200507: Modified minimum default division length.

TODO: Create function to check deviations at equal distance divisions from one curve
    for when GetDistances fails on looped curves.
    Finish implementing bMustBePerpAtOpenEnds.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts():

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


    key = 'bMustBePerpAtOpenEnds'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDist_Tol_Max'; keys.append(key)
    values[key] = (
        0.026 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches
        else (
            3.6 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters
            else 1000.0*sc.doc.ModelAbsoluteTolerance)
        )
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fDist_Tol_Min'; keys.append(key)
    values[key] = (
        0.014 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches
        else (
            2.4 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters
            else 10.0*sc.doc.ModelAbsoluteTolerance)
        )
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fDist_maxToAccept'; keys.append(key)
    values[key] = 10.0 * values['fDist_Tol_Max']
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAddLine'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotDecPlaces'; keys.append(key)
    values[key] = sc.doc.ModelDistanceDisplayPrecision - 1
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=0)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iDotFontHt'; keys.append(key)
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iAddMaxOpt'; keys.append(key)
    values[key] = 2
    riAddOpts[key] = addOptionList(key, names, ('None', 'One', 'All'), values)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iAddMinOpt'; keys.append(key)
    values[key] = 2
    riAddOpts[key] = addOptionList(key, names,  ('None', 'One', 'All'), values)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
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
                print s
            return tuple(gCrvs_Preselected)


def getInput():
    """
    Get curves with optional input.
    """
    
    # For iAddMaxOpt and iAddMaxOpt, 0 == None, 1 == 1, 2 == All.
    
    rgCrvs_lists = [[],[]]
    dictOptListIdxs = {'idxAddMaxOpt': None, 'idxAddMinOpt': None}
    
    go = ri.Custom.GetObject()
    sCommandPromptAdd = 'A', 'B'
    
    go.GeometryFilter = rd.ObjectType.Curve
    
    
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.GroupSelect = True
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    go.AcceptNumber(True, True)

    bPreselectedObjsChecked = False
    
    print "'0' for MinTol will accept all minimum curve deviations." \
            "  Likewise for MaxTol."
    
    idxs_Opts = {}

    # Get input for each curve set using a loop.
    for iCrvSet in 0,1:

        go.SetCommandPrompt("Select curve set {}".format(sCommandPromptAdd[iCrvSet]))

        while True:
            for key in Opts.keys: idxs_Opts[key] = None
            key = 'bMustBePerpAtOpenEnds'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            key = 'fDist_Tol_Max'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            key = 'fDist_Tol_Min'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            key = 'fDist_maxToAccept'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            key = 'bAddLine'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            key = 'bAddDot'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            if Opts.values['bAddDot']:
                key = 'iDotDecPlaces'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
                key = 'iDotFontHt'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            if Opts.values['bAddLine'] or Opts.values['bAddDot']:
                key = 'iAddMaxOpt'; idxs_Opts[key] = Opts.riAddOpts[key](go)
                key = 'iAddMinOpt'; idxs_Opts[key] = Opts.riAddOpts[key](go)
            key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
            key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]

            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
            
            # Use bPreselectedObjsChecked so that only objects before the
            # first call to go.GetMultiple is considered.
            if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
                bPreselectedObjsChecked = True
                go.EnablePreSelect(False, True)
                continue
            
            if res == ri.GetResult.Cancel:
                go.Dispose()
                return
            
            if res == ri.GetResult.Object:
                break

            # An option was selected or a number was entered.

            if res == ri.GetResult.Number:
                Opts.riOpts['fDist_Tol_Max'].CurrentValue = abs(go.Number())
            elif res == ri.GetResult.Option:
                if go.Option().Index == idxs_Opts['iAddMaxOpt']:
                    Opts.values['iAddMaxOpt'] = (
                        go.Option().CurrentListOptionIndex)
                elif go.Option().Index == idxs_Opts['iAddMinOpt']:
                    Opts.values['iAddMinOpt'] = (
                        go.Option().CurrentListOptionIndex)

            if Opts.riOpts['fDist_Tol_Max'].CurrentValue == 0.0:
                Opts.riOpts['fDist_Tol_Max'].CurrentValue = (
                        Opts.riOpts['fDist_Tol_Max'].InitialValue)
            
            if (
                Opts.riOpts['fDist_maxToAccept'].CurrentValue <
                Opts.riOpts['fDist_Tol_Max'].CurrentValue
            ):
                Opts.riOpts['fDist_maxToAccept'].CurrentValue = (
                    10.0 * Opts.riOpts['fDist_Tol_Max'].CurrentValue)
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue


            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()

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


def getClosestDistsBtwn2Crvs(rgCrv_A, rgCrv_B, tolerance=None):
    """
    Alternative to Curve.GetDistancesBetweenCurves for better results when
    curves contain loops, etc.
    Returns:
    """


    def func_toRepeat(rgC_Cat, rgC_Dog):
        ts_Cat = []

        rc = rgC_Cat.DivideByLength(
            segmentLength=fDivLength,
            includeEnds=True)
        if rc: ts_Cat.extend(rc)
        if not rgC_Cat.IsClosed:
            # DivideByLength doesn't add the T1 segment
            # even though includeEnds == True.
            # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
            # shows the parameter labeled as 'includeStart'.
            ts_Cat.append(rgC_Cat.Domain.T1)

        # Add parmeters of internal span ends.
        if rgC_Cat.SpanCount > 1:
            for iSpan in range(1, rgC_Cat.SpanCount):
                spanDomain = rgC_Cat.SpanDomain(iSpan)
                if spanDomain.T0 not in ts_Cat: ts_Cat.append(spanDomain.T0)

        if len(ts_Cat) != len(set(ts_Cat)):
            raise ValueError("Duplicate parameters?  Check getClosestDistsBtwn2Crvs.")

        dists_per_tCat = []
        ts_Cat_per_dist_per_tCat = []
        ts_Dog_per_dist_per_tCat = []

        for iT_Cat in xrange(len(ts_Cat)):

            t_Cat = ts_Cat[iT_Cat]

            pt_Cat = rgC_Cat.PointAt(t_Cat)

            bSuccess, t_Dog = rgC_Dog.ClosestPoint(pt_Cat)

            if not bSuccess:
                raise ValueError("Closest point could not be calculated.")

            ts_Cat_per_dist_per_tCat.append(t_Cat)

            pt_Dog = rgC_Dog.PointAt(t_Dog)
            dists_per_tCat.append(pt_Cat.DistanceTo(pt_Dog))
            ts_Dog_per_dist_per_tCat.append(t_Dog)

        # Get max distance.
        dist_Max = max(dists_per_tCat)
        idx_MaxDist = dists_per_tCat.index(dist_Max)
        t_Cat_MaxDist = ts_Cat_per_dist_per_tCat[idx_MaxDist]
        t_Dog_MaxDist = ts_Dog_per_dist_per_tCat[idx_MaxDist]
        #sc.doc.Objects.AddLine(rg.Line(rgC_Cat.PointAt(t_Cat_MaxDist), rgC_Dog.PointAt(t_Dog_MaxDist)))

        # Get min distance.
        dist_Min = min(dists_per_tCat)
        idx_MinDist = dists_per_tCat.index(dist_Min)
        t_Cat_MinDist = ts_Cat_per_dist_per_tCat[idx_MinDist]
        t_Dog_MinDist = ts_Dog_per_dist_per_tCat[idx_MinDist]

        return (
            dist_Max, 
            t_Cat_MaxDist,
            t_Dog_MaxDist,
            dist_Min,
            t_Cat_MinDist,
            t_Dog_MinDist,
            )


    fDivLength = 10.0*sc.doc.ModelAbsoluteTolerance if tolerance is None else tolerance


    rc = func_toRepeat(rgCrv_A, rgCrv_B)
    (
        dist_Max_BperA,
        tA_MaxDist_BperA,
        tB_MaxDist_BperA,
        dist_Min_BperA,
        tA_MinDist_BperA,
        tB_MinDist_BperA,
       ) = rc

    rc = func_toRepeat(rgCrv_B, rgCrv_A)
    (
        dist_Max_AperB,
        tB_MaxDist_AperB,
        tA_MaxDist_AperB,
        dist_Min_AperB,
        tB_MinDist_AperB,
        tA_MinDist_AperB,
       ) = rc


    if dist_Max_BperA > dist_Max_AperB:
        dist_Max = dist_Max_BperA
        tA_Max = tA_MaxDist_BperA
        tB_Max = tB_MaxDist_BperA
    else:
        dist_Max = dist_Max_AperB
        tA_Max = tA_MaxDist_AperB
        tB_Max = tB_MaxDist_AperB

    if dist_Min_BperA < dist_Min_AperB:
        dist_Min = dist_Min_BperA
        tA_Min = tA_MinDist_BperA
        tB_Min = tB_MinDist_BperA
    else:
        dist_Min = dist_Min_AperB
        tA_Min = tA_MinDist_AperB
        tB_Min = tB_MinDist_AperB


    # Check for intersection.
    intersections = rg.Intersect.Intersection.CurveCurve(
        rgCrv_A,
        rgCrv_B,
        tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
        overlapTolerance=0.0)
    if intersections:
        dist_Min = 0.0
        tA_Min = intersections[0].ParameterA
        tB_Min = intersections[0].ParameterB


    #sc.doc.Objects.AddPoint(rgCrv_A.PointAt(tA_Max))
    #sc.doc.Objects.AddPoint(rgCrv_B.PointAt(tB_Max))
    #sc.doc.Objects.AddPoint(rgCrv_A.PointAt(tA_Min))
    #if dist_Min > 0.0:
    #    sc.doc.Objects.AddPoint(rgCrv_B.PointAt(tB_Min))
    #sc.doc.Views.Redraw()

    return True, dist_Max, tA_Max, tB_Max, dist_Min, tA_Min, tB_Min


def getPerpDistsBtwn2Crvs(rgCrv_A, rgCrv_B, tolerance=None, **kwargs):
    """
    Alternative to Curve.GetDistancesBetweenCurves for better results when
    curves contain loops, etc.
    Returns:
    """


    def func_toRepeat(rgC_Cat, rgC_Dog):
        ts_Cat = []

        rc = rgC_Cat.DivideByLength(
            segmentLength=fDivLength,
            includeEnds=True)
        if rc: ts_Cat.extend(rc)
        if not rgC_Cat.IsClosed:
            # DivideByLength doesn't add the T1 segment
            # even though includeEnds == True.
            # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
            # shows the parameter labeled as 'includeStart'.
            ts_Cat.append(rgC_Cat.Domain.T1)

        # Add parmeters of internal span ends.
        if rgC_Cat.SpanCount > 1:
            for iSpan in range(1, rgC_Cat.SpanCount):
                spanDomain = rgC_Cat.SpanDomain(iSpan)
                if spanDomain.T0 not in ts_Cat: ts_Cat.append(spanDomain.T0)

        if len(ts_Cat) != len(set(ts_Cat)):
            raise ValueError("Duplicate parameters?  Check getPerpDistsBtwn2Crvs.")

        dists_per_tCat = []
        ts_Cat_per_dist_per_tCat = []
        ts_Dog_per_dist_per_tCat = []

        for iT_Cat in xrange(len(ts_Cat)):

            t_Cat = ts_Cat[iT_Cat]

            bSuccess, plane_A = rgC_Cat.PerpendicularFrameAt(t=t_Cat)
            if not bSuccess:
                raise ValueError("Perpendicular plane could not be calculated.")

            intersections = rg.Intersect.Intersection.CurvePlane(
                rgC_Dog,
                plane_A,
                tolerance=0.1*sc.doc.ModelAbsoluteTolerance)

            if not intersections:
                # No intersection, so don't bother proceeding with this t_A.
                #dists_per_tCat.append(None)
                #ts_Cat_per_dist_per_tCat.append(None)
                #ts_Dog_per_dist_per_tCat.append(None)
                continue # to next t_A.

            pt_Cat = rgC_Cat.PointAt(t_Cat)

            ts_Cat_per_dist_per_tCat.append(t_Cat)

            if len(intersections) == 1:
                pt_Dog = intersections[0].PointA
                dists_per_tCat.append(pt_Cat.DistanceTo(pt_Dog))
                ts_Dog_per_dist_per_tCat.append(intersections[0].ParameterA)
            else:
                ts_Dog = []
                dists = []
                for intersection in intersections:
                    pt_Dog = intersection.PointA
                    dists.append(pt_Cat.DistanceTo(pt_Dog))
                    ts_Dog.append(intersection.ParameterA)
                    
                # Using min() so get closest intersection of looped curves.
                dists_per_tCat.append(min(dists))
                ts_Dog_per_dist_per_tCat.append(ts_Dog[dists.index(min(dists))])

        # Get max distance.
        dist_Max = max(dists_per_tCat)
        idx_MaxDist = dists_per_tCat.index(dist_Max)
        t_Cat_MaxDist = ts_Cat_per_dist_per_tCat[idx_MaxDist]
        t_Dog_MaxDist = ts_Dog_per_dist_per_tCat[idx_MaxDist]
        #sc.doc.Objects.AddLine(rg.Line(rgC_Cat.PointAt(t_Cat_MaxDist), rgC_Dog.PointAt(t_Dog_MaxDist)))

        # Get min distance.
        dist_Min = min(dists_per_tCat)
        idx_MinDist = dists_per_tCat.index(dist_Min)
        t_Cat_MinDist = ts_Cat_per_dist_per_tCat[idx_MinDist]
        t_Dog_MinDist = ts_Dog_per_dist_per_tCat[idx_MinDist]

        return (
            dist_Max, 
            t_Cat_MaxDist,
            t_Dog_MaxDist,
            dist_Min,
            t_Cat_MinDist,
            t_Dog_MinDist,
            )


    fDivLength = 10.0*sc.doc.ModelAbsoluteTolerance if tolerance is None else tolerance


    rc = func_toRepeat(rgCrv_A, rgCrv_B)
    (
        dist_Max_BperA,
        tA_MaxDist_BperA,
        tB_MaxDist_BperA,
        dist_Min_BperA,
        tA_MinDist_BperA,
        tB_MinDist_BperA,
       ) = rc

    rc = func_toRepeat(rgCrv_B, rgCrv_A)
    (
        dist_Max_AperB,
        tB_MaxDist_AperB,
        tA_MaxDist_AperB,
        dist_Min_AperB,
        tB_MinDist_AperB,
        tA_MinDist_AperB,
       ) = rc


    if dist_Max_BperA > dist_Max_AperB:
        dist_Max = dist_Max_BperA
        tA_Max = tA_MaxDist_BperA
        tB_Max = tB_MaxDist_BperA
    else:
        dist_Max = dist_Max_AperB
        tA_Max = tA_MaxDist_AperB
        tB_Max = tB_MaxDist_AperB

    if dist_Min_BperA < dist_Min_AperB:
        dist_Min = dist_Min_BperA
        tA_Min = tA_MinDist_BperA
        tB_Min = tB_MinDist_BperA
    else:
        dist_Min = dist_Min_AperB
        tA_Min = tA_MinDist_AperB
        tB_Min = tB_MinDist_AperB


    # Check for intersection.
    intersections = rg.Intersect.Intersection.CurveCurve(
        rgCrv_A,
        rgCrv_B,
        tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
        overlapTolerance=0.0)
    if intersections:
        dist_Min = 0.0
        tA_Min = intersections[0].ParameterA
        tB_Min = intersections[0].ParameterB


    sc.doc.Objects.AddPoint(rgCrv_A.PointAt(tA_Max))
    sc.doc.Objects.AddPoint(rgCrv_B.PointAt(tB_Max))
    sc.doc.Objects.AddPoint(rgCrv_A.PointAt(tA_Min))
    if dist_Min > 0.0:
        sc.doc.Objects.AddPoint(rgCrv_B.PointAt(tB_Min))
    sc.doc.Views.Redraw()

    return True, dist_Max, tA_Max, tB_Max, dist_Min, tA_Min, tB_Min


def getMaxPerpDistBtwn2Sets(rgCrvs_SetA, rgCrvs_SetB, fCalculationTol=None, **kwargs):
    """
    Returns:
        fDists_max, ptsA_max, ptsB_max, fDists_min, ptsA_min, ptsB_min
        Unlike Curve.GetDistancesBetweenCurves, doesn't return a boolean value for success.
        In the case of fail, empty lists are the values of each variable.
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDist_maxToAccept = getOpt('fDist_maxToAccept')
    bDebug = getOpt('bDebug')


    try: rgCrvs_SetA = list(rgCrvs_SetA)
    except: rgCrvs_SetA = [rgCrvs_SetA]

    try: rgCrvs_SetB = list(rgCrvs_SetB)
    except: rgCrvs_SetB = [rgCrvs_SetB]

    if fCalculationTol is None:
        fCalculationTol = sc.doc.ModelAbsoluteTolerance

    fDists_max = []; ptsA_max = []; ptsB_max = []
    fDists_min = []; ptsA_min = []; ptsB_min = []


    # Due to a bug in V5 (SR12) bug
    # (fixed in V6 per https://mcneel.myjetbrains.com/youtrack/issue/RH-41877),
    # bad parameter values are returned from GetDistancesBetweenCurves for ArcCurves,
    # so convert any found to a NurbsCurve.
    if Rhino.RhinoApp.ExeVersion == 5:
        for i, rgC_A in enumerate(rgCrvs_SetA):
            if isinstance(rgC_A, rg.ArcCurve):
                rgCrvs_SetA[i] = rgC_A.ToNurbsCurve()
        for i, rgC_B in enumerate(rgCrvs_SetB):
            if isinstance(rgC_B, rg.ArcCurve):
                rgCrvs_SetB[i] = rgC_B.ToNurbsCurve()


    fDivisionLength = 10.0*sc.doc.ModelAbsoluteTolerance

    dists_Max_perA = []
    pts_B_MaxDist_perA = []

    for rgC_A in rgCrvs_SetA:

        dists_per_tA = []
        pts_B_Dists_perA = []

        ts_Crv_A = []

        rc = rgC_A.DivideByLength(
            segmentLength=fDivisionLength,
            includeEnds=True)
        if rc: ts_Crv_A.extend(rc)
        if not rgC_A.IsClosed:
            # DivideByLength doesn't add the T1 segment
            # even though includeEnds == True.
            # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
            # shows the parameter labeled as 'includeStart'.
            ts_Crv_A.append(rgC_A.Domain.T1)

        # Add parmeters of segment ends.


        planes_A = []

        for t in ts_Crv_A:
            pt = rgC_A.PointAt(t)
            #sc.doc.Objects.AddPoint(pt)

            bSuccess, plane = rgC_A.PerpendicularFrameAt(t=t)
            if not bSuccess:
                raise ValueError("Perpendicular plane could not be calculated.")
            planes_A.append(plane)


        for rgC_B in rgCrvs_SetB:
            sc.escape_test()

            for iT_A in xrange(len(ts_Crv_A)):

                t_A = ts_Crv_A[iT_A]
                plane_A = planes_A[iT_A]
                
                intersections = rg.Intersect.Intersection.CurvePlane(
                    rgC_B,
                    plane_A,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if not intersections:
                    # No intersection, so don't bother proceeding with this nc_B.
                    dists_per_tA.append(None)
                    pts_B_Dists_perA.append(None)
                    continue # to next rgC_B.

                pt_A = rgC_A.PointAt(t_A)

                if len(intersections) == 1:
                    pt_B = intersections[0].PointA
                    dists_per_tA.append(pt_A.DistanceTo(pt_B))
                    pts_B_Dists_perA.append(pt_B)
                else:
                    pts_B = []
                    dists = []
                    for intersection in intersections:
                        pt_B = intersection.PointA # PointA is not in reference to rgCrvs_A.
                        dists.append(pt_A.DistanceTo(pt_B))
                        pts_B.append(pt_B)
                    
                    # Using min() so get closest intersection of looped curves.
                    dists_per_tA.append(min(dists))
                    pts_B_Dists_perA.append(pts_B[dists_per_tA.index(min(dists))])

        # Get max of all ts_A of this rgC_A for all rgCrvs_SetB.
        dists_Max_perA.append(max(dists_per_tA))
        print max(dists_per_tA)
        print dists_per_tA.index(max(dists_per_tA))
        print pts_B_Dists_perA
        pts_B_MaxDist_perA.append(pts_B_Dists_perA[dists_per_tA.index(max(dists_per_tA))])

    for pt_B in pts_B_MaxDist_perA:
        sc.doc.Objects.AddPoint(pt_B)
    sc.doc.Views.Redraw()




    1/0

    return fDists_max, ptsA_max, ptsB_max, fDists_min, ptsA_min, ptsB_min


def getDevsBtwn2Sets(rgCrvs_SetA, rgCrvs_SetB, fCalculationTol=None, **kwargs):
    """
    Returns:
        fDists_max, ptsA_max, ptsB_max, fDists_min, ptsA_min, ptsB_min
        Unlike Curve.GetDistancesBetweenCurves, doesn't return a boolean value for success.
        In the case of fail, empty lists are the values of each variable.
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDist_maxToAccept = getOpt('fDist_maxToAccept')
    bDebug = getOpt('bDebug')


    try: rgCrvs_SetA = list(rgCrvs_SetA)
    except: rgCrvs_SetA = [rgCrvs_SetA]

    try: rgCrvs_SetB = list(rgCrvs_SetB)
    except: rgCrvs_SetB = [rgCrvs_SetB]

    if fCalculationTol is None:
        fCalculationTol = sc.doc.ModelAbsoluteTolerance

    fDists_max = []; ptsA_max = []; ptsB_max = []
    fDists_min = []; ptsA_min = []; ptsB_min = []


    # Explode parent curves into segments.

    rgCs_Segs_A = []
    rgCs_Segs_B = []

    for rgC_SetA in rgCrvs_SetA:
        rgCs = rgC_SetA.DuplicateSegments()
        if rgCs:
            rgCs_Segs_A.extend(rgCs)
        else:
             # DuplicateSegments returns an empty array when
             # curve is not composed of multiple segments.
            rgCs_Segs_A.append(rgC_SetA)

    for rgC_SetB in rgCrvs_SetB:
        rgCs = rgC_SetB.DuplicateSegments()
        if rgCs: rgCs_Segs_B.extend(rgCs)
        else: rgCs_Segs_B.append(rgC_SetB)


    # Due to a bug in V5 (SR12) bug
    # (fixed in V6 per https://mcneel.myjetbrains.com/youtrack/issue/RH-41877),
    # bad parameter values are returned from GetDistancesBetweenCurves for ArcCurves,
    # so convert any found to a NurbsCurve.
    if Rhino.RhinoApp.ExeVersion == 5:
        for i, cA in enumerate(rgCs_Segs_A):
            if isinstance(cA, rg.ArcCurve):
                rgCs_Segs_A[i] = cA.ToNurbsCurve()
        for i, cB in enumerate(rgCs_Segs_B):
            if isinstance(cB, rg.ArcCurve):
                rgCs_Segs_B[i] = cB.ToNurbsCurve()


    for cA in rgCs_Segs_A:
        for cB in rgCs_Segs_B:
            sc.escape_test()
            rc = rg.Curve.GetDistancesBetweenCurves(
                    cA, cB, tolerance=fCalculationTol)
            
            if not rc[0]:
                pass
                #print "GetDistancesBetweenCurves didn't calculate."
            else:
                (
                    fDist_max, tA_Max, tB_Max,
                    fDist_min, tA_Min, tB_Min
                ) = rc[1:]
                
                #sc.doc.Objects.AddCurve(cA)
                #sc.doc.Objects.AddCurve(cB)
                #sc.doc.Objects.AddPoint(cA.PointAt(tA_Max))
                #sc.doc.Objects.AddPoint(cB.PointAt(tB_Max))
                

                # Eliminate any maximum results of distances not perpendicular to curves.
                # This often happens at end points of pairs of staggered curves.
                ptA_FromGDBC = cA.PointAt(tA_Max)
                ptB_FromGDBC = cB.PointAt(tB_Max)
                #line = rg.Line(ptA_FromGDBC, ptB_FromGDBC)
                #sc.doc.Objects.AddLine(line)

                if fDist_maxToAccept is not None and fDist_max > fDist_maxToAccept:
                    if bDebug: print "Skipped distance {}.".format(fDist_max)
                else:
                    # Eliminate any maximum results of distances not mutually
                    # the closest point between the curves.
                    # This often happens at end points of pairs of staggered curves.
                    ptA_FromGDBC = cA.PointAt(tA_Max)
                    ptB_FromGDBC = cB.PointAt(tB_Max)
                    rc = cA.ClosestPoint(ptB_FromGDBC)
                    if rc[0]:
                        ptOnA_ClosestTo_ptB = cA.PointAt(rc[1])
                    rc = cB.ClosestPoint(ptA_FromGDBC)
                    if rc[0]:
                        ptOnB_ClosestTo_ptA = cB.PointAt(rc[1])
                    dist_ptA_to_B = ptA_FromGDBC.DistanceTo(ptOnB_ClosestTo_ptA)
                    dist_ptB_to_A = ptB_FromGDBC.DistanceTo(ptOnA_ClosestTo_ptB)
                    if (
                            abs(dist_ptA_to_B - fDist_max) <= 0.1*sc.doc.ModelAbsoluteTolerance
                            and
                            abs(dist_ptB_to_A - fDist_max) <= 0.1*sc.doc.ModelAbsoluteTolerance
                    ):
                        fDists_max.append(fDist_max)
                        ptsA_max.append(ptA_FromGDBC)
                        ptsB_max.append(ptB_FromGDBC)
                
                fDists_min.append(fDist_min)
                ptsA_min.append(cA.PointAt(tA_Min))
                ptsB_min.append(cB.PointAt(tB_Min))

    return fDists_max, ptsA_max, ptsB_max, fDists_min, ptsA_min, ptsB_min


def addTextDot(text, pt, iDotFontHt=14):
    rgDot = rg.TextDot(text, pt)
    rgDot.FontHeight = iDotFontHt
    sc.doc.Objects.AddTextDot(rgDot)


def addAnnotation(fDists_max, ptsA_max, ptsB_max, fDists_min, ptsA_min, ptsB_min, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAddLine = getOpt('bAddLine')
    bAddDot = getOpt('bAddDot')
    iDotDecPlaces = getOpt('iDotDecPlaces')
    iDotFontHt = getOpt('iDotFontHt')
    iAddMaxOpt = getOpt('iAddMaxOpt')
    iAddMinOpt = getOpt('iAddMinOpt')
    fDist_Tol_Max = getOpt('fDist_Tol_Max')
    fDist_Tol_Min = getOpt('fDist_Tol_Min')


    fDist_max_all = max(fDists_max) if fDists_max else None
    fDist_min_all = min(fDists_min) if fDists_min else None
    
    if not (    (bAddLine or bAddDot) and
            (iAddMaxOpt or iAddMinOpt) or
            (fDist_max_all or fDist_min_all)):
        return
    
    if bAddLine:
        attr = rd.ObjectAttributes()
        attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex # Why is this necessary?
        #attr.ObjectDecoration = rd.ObjectDecoration.BothArrowhead
    
    if bAddLine or bAddDot:
        
        # Maximum deviation(s).
        if iAddMaxOpt > 0 and fDists_max:
            if iAddMaxOpt == 1:
                fDist_max_all = max(fDists_max)
                idxDist_max_all = fDists_max.index(fDist_max_all)
                ptA_max_all = ptsA_max[idxDist_max_all]
                ptB_max_all = ptsB_max[idxDist_max_all]
                if bAddLine:
                    sc.doc.Objects.AddLine(ptA_max_all, ptB_max_all, attr)
                if bAddDot:
                    addTextDot('{0:.{1}f}'.format(fDist_max_all, iDotDecPlaces),
                            (ptA_max_all + ptB_max_all) / 2.0, iDotFontHt)
            elif iAddMaxOpt == 2:
                for i, fDist_max in enumerate(fDists_max):
                    if fDist_max > fDist_Tol_Max:
                        ptA_max = ptsA_max[i]
                        ptB_max = ptsB_max[i]
                        if (fDist_max is not None
                                and ptA_max is not None and ptB_max is not None):
                            if bAddLine:
                                sc.doc.Objects.AddLine(ptA_max, ptB_max, attr)
                            if bAddDot:
                                addTextDot('{0:.{1}f}'.format(fDist_max, iDotDecPlaces),
                                        (ptA_max + ptB_max) / 2.0, iDotFontHt)
        
        # Minimum deviation(s).
        if iAddMinOpt > 0 and fDists_min:
            if iAddMinOpt == 1:
                fDist_min_all = min(fDists_min)
                idxDist_min_all = fDists_min.index(fDist_min_all)
                ptA_min_all = ptsA_min[idxDist_min_all]
                ptB_min_all = ptsB_min[idxDist_min_all]
                if bAddLine:
                    sc.doc.Objects.AddLine(ptA_min_all, ptB_min_all, attr)
                if bAddDot:
                    
                    addTextDot('{0:.{1}f}'.format(fDist_min_all, iDotDecPlaces),
                            (ptA_min_all + ptB_min_all) / 2.0, iDotFontHt)
            elif iAddMinOpt == 2:
                for i, fDist_min in enumerate(fDists_min):
                    if fDist_Tol_Min == 0 or fDist_min < fDist_Tol_Min:
                        ptA_min = ptsA_min[i]
                        ptB_min = ptsB_min[i]
                        if (fDist_min is not None
                                and ptA_min is not None and ptB_min is not None):
                            if bAddLine:
                                sc.doc.Objects.AddLine(ptA_min, ptB_min, attr)
                            if bAddDot:
                                addTextDot('{0:.{1}f}'.format(fDist_min, iDotDecPlaces),
                                        (ptA_min + ptB_min) / 2.0, iDotFontHt)


def summaryLog(fDists_max, fDists_min):
    
    fDist_max_all = max(fDists_max) if fDists_max else None
    fDist_min_all = min(fDists_min) if fDists_min else None
    
    if not fDists_max and not fDists_min:
        s = "Curves do not overlap."
    else:
        s = ""
        if fDist_min_all:
            s += "Minimum deviation = {0:.{1}f}".format(
                    fDist_min_all, sc.doc.ModelDistanceDisplayPrecision)
        if fDist_max_all:
            s += "\nMaximum deviation = {0:.{1}f}".format(
                    fDist_max_all, sc.doc.ModelDistanceDisplayPrecision)
    return s


def main(bDebug=False):
    
    rc = getPreselectedCurves()
    if rc:
        rgCrvsA = [sc.doc.Objects.FindId(rc[0]).Geometry]
        rgCrvsB = [sc.doc.Objects.FindId(rc[1]).Geometry]
    else:
        rc = getInput()
        if rc is None: return
        print rc

        rgCrvsA = rc[0]
        rgCrvsB = rc[1]

    rc = getDevsBtwn2Sets(
            rgCrvsA[0],
            rgCrvsB[0],
            )

    fDists_max, ptsA_max, ptsB_max, fDists_min, ptsA_min, ptsB_min = rc
    
    addAnnotation(
        fDists_max, ptsA_max, ptsB_max,
        fDists_min, ptsA_min, ptsB_min,
        )
    
    print summaryLog(fDists_max, fDists_min)
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main(bDebug=bool(0))