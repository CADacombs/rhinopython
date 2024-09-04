"""
This script is an alternative to _CrvDeviation.

'0' for MinTol will accept all minimum curve deviations. Likewise for MaxTol.
MaxDistToAccept is the lowest value, over which the distances will be ignored.


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
190829: Now, ignores maximum results where the distances are not perpendicular to either curve.
        Testing having fDist_maxToAccept disabled.
190831: Renabled fDist_maxToAccept.
200216: Added comments to a function.
200409-10: Refactored.  Added some functions.  WIP: Added bMustBePerpAtOpenEnds.
200507: Modified minimum default division length.
230927: Modified some printed output.
240903-04: Removed some functions. Refactored.


TODO:
    Create function to check deviations at equal distance divisions from one curve
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
    listValues = {}
    stickyKeys = {}


    key = 'bCrvSets'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMustBePerpAtOpenEnds'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDivLength'; keys.append(key)
    values[key] = 100.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fDist_maxToAccept'; keys.append(key)
    names[key] = 'MaxDistToAccept'
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 1.0
    elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
        values[key] = 25
    else:
        values[key] = 1000.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iAnnotations'; keys.append(key)
    listValues[key] = 'None', 'MaxOnly', 'MinOnly', 'Both'
    values[key] = 0
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

    key = 'iAddMaxOpt'; keys.append(key)
    listValues[key] = 'None', 'One', 'All'
    values[key] = 2
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iAddMinOpt'; keys.append(key)
    listValues[key] = 'None', 'One', 'All'
    values[key] = 2
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDist_Tol_Max'; keys.append(key)
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 0.026
    elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
        values[key] = 3.6
    else:
        values[key] = 1000.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'MaxTol'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fDist_Tol_Min'; keys.append(key)
    values[key] = (
        0.014 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches
        else (
            2.4 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters
            else 10.0*sc.doc.ModelAbsoluteTolerance)
        )
    names[key] = 'MinTol'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

        if key == 'fDivLength':
            if cls.riOpts[key].CurrentValue <= 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue <= sc.doc.ModelAbsoluteTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

        if key == 'fDist_Tol_Max':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            cls.values[key] = cls.riOpts[key].CurrentValue
            if (cls.values[key] > 0.0) and (cls.values[key] > cls.values['fDist_maxToAccept']):
                cls.values['fDist_maxToAccept'] = cls.riOpts['fDist_maxToAccept'].CurrentValue = cls.values[key]
                sc.sticky[cls.stickyKeys['fDist_maxToAccept']] = cls.values['fDist_maxToAccept']
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fDist_Tol_Min':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
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


def getInput_2Crvs():
    """
    Get 2 curves with optional input.
    """

    # For iAddMaxOpt and iAddMaxOpt, 0 == None, 1 == 1, 2 == All.

    dictOptListIdxs = {'idxAddMaxOpt': None, 'idxAddMinOpt': None}
    
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
        if Opts.values['bCrvSets']:
            go.Dispose()
            sc.doc.Objects.UnselectAll()
            sc.doc.Views.Redraw()
            return getInput_2Sets()


        go.ClearCommandOptions()
        idxs_Opts.clear()

        addOption('bCrvSets')
        #addOption('bMustBePerpAtOpenEnds')
        addOption('fDivLength')
        addOption('fDist_maxToAccept')
        addOption('iAnnotations')
        if Opts.values['iAnnotations'] > 0:
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
            if Opts.values['iAnnotations'] in (1,3):
                addOption('iAddMaxOpt')
                addOption('fDist_Tol_Max')
            if Opts.values['iAnnotations'] in (2,3):
                addOption('iAddMinOpt')
                addOption('fDist_Tol_Min')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        # Use bPreselectedObjsChecked so that only selected objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
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
            key = 'fDivLength'
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
    
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False
    
    idxs_Opts = {}
    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    # Get input for each curve set using a loop.
    for iCrvSet in 0,1:

        go.SetCommandPrompt("Select curve set {}".format(sCommandPromptAdd[iCrvSet]))

        while True:
            if not Opts.values['bCrvSets']:
                go.Dispose()
                sc.doc.Objects.UnselectAll()
                sc.doc.Views.Redraw()
                return getInput_2Crvs


            go.ClearCommandOptions()
            idxs_Opts.clear()

            addOption('bCrvSets')
            #addOption('bMustBePerpAtOpenEnds')
            addOption('fDist_Tol_Max')
            addOption('fDist_Tol_Min')
            addOption('fDist_maxToAccept')
            addOption('bAddLine')
            addOption('bAddDot')
            if Opts.values['bAddDot']:
                addOption('iDotDecPlaces')
                addOption('iDotFontHt')
            if Opts.values['bAddLine'] or Opts.values['bAddDot']:
                addOption('iAddMaxOpt')
                addOption('iAddMinOpt')
            addOption('bEcho')
            addOption('bDebug')

            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
            
            # Use bPreselectedObjsChecked so that only selected objects before the
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
                if ('iAddMaxOpt' in idxs_Opts) and (go.Option().Index == idxs_Opts['iAddMaxOpt']):
                    Opts.values['iAddMaxOpt'] = (
                        go.Option().CurrentListOptionIndex)
                elif ('iAddMinOpt' in idxs_Opts) and (go.Option().Index == idxs_Opts['iAddMinOpt']):
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
            #Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue


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


def GDBCs_perClosestPoint_1Dir(curveA, curveB, fDivLength_Starting, bDebug=False):
    """
    Alternative to Curve.GetDistancesBetweenCurves using Curve.ClosestPoint.

    curveA is divided per segmentLength to obtain some testPoints for ClosestPoint.
    End points of segments of curveA are also added.

    curveB is the object of ClosestPoint.

    Parameters:
        curveA: rg.Curve
        curveB: rg.Curve
        segmentLength: float : Division length of curveA to obtain some testPoints for ClosestPoint.
        bDebug: bool
    """

    segmentLength = fDivLength_Starting

    def generate_parameters_for_ClosestPoint(curve, segmentLength):

        ts_Out = []

        rc = curve.DivideByLength(
            segmentLength=segmentLength,
            includeEnds=True)
        if rc:
            ts_Out.extend(rc)

        # Add parmeters of start of each span.
        if curve.SpanCount > 1:
            for iSpan in range(1, curve.SpanCount):
                spanDomain = curve.SpanDomain(iSpan)
                if spanDomain.T0 not in ts_Out: ts_Out.append(spanDomain.T0)

        # For open curveA, add the end of the curve.
        if not curve.IsClosed:
            # DivideByLength doesn't add the T1 segment
            # even though includeEnds == True.
            # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
            # shows the parameter labeled as 'includeStart'.
            if curve.Domain.T1 not in ts_Out:
                ts_Out.append(curve.Domain.T1)

        if len(ts_Out) != len(set(ts_Out)):
            sEval = "len(ts_Out)"; print(sEval, '=', eval(sEval))
            sEval = "len(set(ts_Out))"; print(sEval, '=', eval(sEval))
            raise ValueError("Duplicate parameters?  Check getClosestDistsBtwn2Crvs.")

        return ts_Out


    ts_A = generate_parameters_for_ClosestPoint(curveA, segmentLength)

    if bDebug: sEval = "len(ts_A)"; print(sEval, '=', eval(sEval))


    def calc_parameters_and_distances(ts_A, curveA, curveB):

        ts_B = []
        dists_per_ts_A = []

        for i_t_A in xrange(len(ts_A)):

            t_A = ts_A[i_t_A]

            pt_A = curveA.PointAt(t_A)

            bSuccess, t_B = curveB.ClosestPoint(pt_A)

            if not bSuccess:
                raise ValueError("Closest point could not be calculated.")

            pt_Dog = curveB.PointAt(t_B)
            dists_per_ts_A.append(pt_A.DistanceTo(pt_Dog))
            ts_B.append(t_B)

        return ts_B, dists_per_ts_A


    ts_B, dists_per_ts_A = calc_parameters_and_distances(ts_A, curveA, curveB)

    # Get max distance.
    dist_Max = max(dists_per_ts_A)
    idx_MaxDist = dists_per_ts_A.index(dist_Max)
    t_A_MaxDist = ts_A[idx_MaxDist]
    t_B_MaxDist = ts_B[idx_MaxDist]
    #sc.doc.Objects.AddLine(rg.Line(curveA.PointAt(t_Cat_MaxDist), curveB.PointAt(t_Dog_MaxDist)))

    # Get min distance.
    dist_Min = min(dists_per_ts_A)
    idx_MinDist = dists_per_ts_A.index(dist_Min)
    t_A_MinDist = ts_A[idx_MinDist]
    t_B_MinDist = ts_B[idx_MinDist]


    # Iterate lower increments about current maximum to find more accurate maximum.

    ts_A_WIP = ts_A[:]
    cA_WIP = curveA.Duplicate()
    idx_MaxDist_WIP = idx_MaxDist
    if bDebug: sEval = "dist_Max"; print(sEval, '=', eval(sEval))

    while True:
        sc.escape_test()

        segmentLength *= 0.1

        if segmentLength < (0.1 * sc.doc.ModelAbsoluteTolerance - 1e-6):
            break

        if bDebug: sEval = "segmentLength"; print(sEval, '=', eval(sEval))

        t0 = ts_A_WIP[idx_MaxDist_WIP - 1]
        t1 = ts_A_WIP[idx_MaxDist_WIP + 1]

        cA_WIP = cA_WIP.Trim(rg.Interval(t0, t1))
        #sc.doc.Objects.AddCurve(cA_WIP); sc.doc.Views.Redraw(); 1/0

        ts_A_WIP = generate_parameters_for_ClosestPoint(cA_WIP, segmentLength)
        if bDebug: sEval = "len(ts_A_WIP)"; print(sEval, '=', eval(sEval))

        ts_B_WIP, dists_per_ts_A_WIP = calc_parameters_and_distances(ts_A_WIP, cA_WIP, curveB)

        # Get max distance.
        dist_Max = max(dists_per_ts_A_WIP)
        idx_MaxDist_WIP = dists_per_ts_A_WIP.index(dist_Max)
        t_A_MaxDist = ts_A_WIP[idx_MaxDist_WIP]
        t_B_MaxDist = ts_B_WIP[idx_MaxDist_WIP]

        if bDebug: sEval = "dist_Max"; print(sEval, '=', eval(sEval))



    return (
        True,
        dist_Max, 
        t_A_MaxDist,
        t_B_MaxDist,
        dist_Min,
        t_A_MinDist,
        t_B_MinDist,
        )


def getDevs_between_2Crvs_perClosestPoint(rgCrv_A, rgCrv_B, fDivLength=None, bDebug=False):
    """
    Alternative to Curve.GetDistancesBetweenCurves for better results when
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


    if fDivLength is None:
        fDivLength = 100.0*sc.doc.ModelAbsoluteTolerance
    if bDebug: sEval = "fDivLength"; print(sEval, '=', eval(sEval))


    rc = GDBCs_perClosestPoint_1Dir(rgCrv_A, rgCrv_B, fDivLength, bDebug=bDebug)
    (
        bSuccess,
        dist_Max_BperA,
        tA_MaxDist_BperA,
        tB_MaxDist_BperA,
        dist_Min_BperA,
        tA_MinDist_BperA,
        tB_MinDist_BperA,
       ) = rc
    if bDebug: sEval = "rc"; print(sEval, '=', eval(sEval))


    rc = GDBCs_perClosestPoint_1Dir(rgCrv_B, rgCrv_A, fDivLength, bDebug=bDebug)
    (
        bSuccess,
        dist_Max_AperB,
        tB_MaxDist_AperB,
        tA_MaxDist_AperB,
        dist_Min_AperB,
        tB_MinDist_AperB,
        tA_MinDist_AperB,
       ) = rc
    if bDebug: sEval = "rc"; print(sEval, '=', eval(sEval))


    if dist_Max_BperA > dist_Max_AperB:
        dist_Max = dist_Max_BperA
        ptA_Max = rgCrv_A.PointAt(tA_MaxDist_BperA)
        ptB_Max = rgCrv_B.PointAt(tB_MaxDist_BperA)
    else:
        dist_Max = dist_Max_AperB
        ptA_Max = rgCrv_A.PointAt(tA_MaxDist_AperB)
        ptB_Max = rgCrv_B.PointAt(tB_MaxDist_AperB)

    if dist_Min_BperA < dist_Min_AperB:
        dist_Min = dist_Min_BperA
        ptA_Min = rgCrv_A.PointAt(tA_MinDist_BperA)
        ptB_Min = rgCrv_B.PointAt(tB_MinDist_BperA)
    else:
        dist_Min = dist_Min_AperB
        ptA_Min = rgCrv_A.PointAt(tA_MinDist_AperB)
        ptB_Min = rgCrv_B.PointAt(tB_MinDist_AperB)


    # Check for intersection to replace any minimum distance values.
    if dist_Min > 2**-(52):
        intersections = rg.Intersect.Intersection.CurveCurve(
            rgCrv_A,
            rgCrv_B,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
            overlapTolerance=0.0)
        if intersections:
            dist_Min = 0.0
            ptA_Min = intersections[0].PointA
            ptB_Min = intersections[0].PointB


    return (
        dist_Max,
        ptA_Max,
        ptB_Max,
        dist_Min,
        ptA_Min,
        ptB_Min
        )


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

    fDists_max = []
    ptsA_max = []
    ptsB_max = []

    fDists_min = []
    ptsA_min = []
    ptsB_min = []


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
                #print("GetDistancesBetweenCurves didn't calculate."
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
                    if bDebug: print("Skipped distance {}.".format(fDist_max))
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
                                and ptA_max is not None and ptB_max is not None
                        ):
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


def _summaryLog(fDists_max, fDists_min, fDist_Tol_Max, fDist_Tol_Min):
    
    fDist_max_all = max(fDists_max) if fDists_max else None
    fDist_min_all = min(fDists_min) if fDists_min else None
    
    if not fDists_max and not fDists_min:
        s = "Curves do not overlap."
    else:
        s = ""
        if fDist_min_all:
            s += "Minimum deviation = {0:.{1}f}".format(
                    fDist_min_all, sc.doc.ModelDistanceDisplayPrecision)
            if fDist_min_all > fDist_Tol_Min:
                s += "  None are below {}.".format(fDist_Tol_Min)
        if fDist_max_all:
            s += "\nMaximum deviation = {0:.{1}f}".format(
                    fDist_max_all, sc.doc.ModelDistanceDisplayPrecision)
            if fDist_max_all < fDist_Tol_Max:
                s += "  None are above {}.".format(fDist_Tol_Max)
    return s


def main():

    rc = getPreselectedCurves()
    if rc:
        rgCrvsA = [sc.doc.Objects.FindId(rc[0]).Geometry]
        rgCrvsB = [sc.doc.Objects.FindId(rc[1]).Geometry]
    else:
        bCrvSets = Opts.values['bCrvSets']
        if bCrvSets:
            rc = getInput_2Sets()
        else:
            rc = getInput_2Crvs()
        if rc is None: return
        #print(rc)

        rgCrvsA = rc[0]
        rgCrvsB = rc[1]


    bCrvSets = Opts.values['bCrvSets']
    #bMustBePerpAtOpenEnds = Opts.values['bMustBePerpAtOpenEnds']
    fDivLength = Opts.values['fDivLength']
    fDist_Tol_Max = Opts.values['fDist_Tol_Max']
    fDist_Tol_Min = Opts.values['fDist_Tol_Min']
    fDist_maxToAccept = Opts.values['fDist_maxToAccept']
    iAnnotations = Opts.values['iAnnotations']
    bAddLine = Opts.values['bAddLine']
    bAddDot = Opts.values['bAddDot']
    iDotDecPlaces = Opts.values['iDotDecPlaces']
    iDotFontHt = Opts.values['iDotFontHt']
    iAddMaxOpt = Opts.values['iAddMaxOpt']
    iAddMinOpt = Opts.values['iAddMinOpt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    if len(rgCrvsA) == 1 and len(rgCrvsB) == 1:
        rgCrv_A = rgCrvsA[0]
        rgCrv_B = rgCrvsB[0]
        rc = getDevs_between_2Crvs_perClosestPoint(
            rgCrv_A,
            rgCrv_B,
            fDivLength=fDivLength,
            bDebug=bDebug)
        print(rc)
        if rc is None:
            print("Deviation failed.")
            return

        fDists_max = [rc[0]]
        ptsA_max = [rc[1]]
        ptsB_max = [rc[2]]
        fDists_min = [rc[3]]
        ptsA_min = [rc[4]]
        ptsB_min = [rc[5]]

        print("Minimum deviation = {}".format(rc[3]))
        print("Maximum deviation = {}".format(rc[0]))

        #return

    else:
        rc = getDevsBtwn2Sets(
                rgCrvsA[0],
                rgCrvsB[0]
                )

        (
            fDists_max,
            ptsA_max,
            ptsB_max,
            fDists_min,
            ptsA_min,
            ptsB_min
           ) = rc

        print(_summaryLog(fDists_max, fDists_min, fDist_Tol_Max, fDist_Tol_Min))


    if iAnnotations:
        addAnnotation(
            fDists_max, ptsA_max, ptsB_max,
            fDists_min, ptsA_min, ptsB_min,
            )
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()