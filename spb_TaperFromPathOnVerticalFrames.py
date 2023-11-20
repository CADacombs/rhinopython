"""
This script is an alternative to _ExrudeCrvTapered.

Although it doesn't have all _ExrudeCrvTapered's options, e.g., Corners,
it has a few additional features:
    It can create loose extrudes.
    It can create variable tapers, where the start and end angles are
        at the ends of each connected path profile.
    It can align the ends of the tapered-to curves with the curves to taper-from curves
        to aid in constructing contiguous breps with G1 continuity.

Limitations:
    There is no equivalent option to _ExtrudeCrvTapered's Corners.
    Periodic curves are not supported.
    Variable taper is disabled for closed, non-periodic, curves.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190604: Created from a split from another script.  Added bAddTaperEndCrvs.
190624: Now, command waits for Enter or options when curves are preselected.
...
220828: Added an option.  Refactored.
231117-20: Added bAlignEndDirs.  Removed bRebuildPath.  Added bAddTaperStartCrvs.
        Replaced an import with some of its code, then simplified it.
        Input of curves can no longer be modified during script execution after their first selection.
        Instead, left clicks will cycle through distance and angle signs.
        Replaced adding DocObjects for pending geometry with DrawConduit until results are accepted.

TODO: 
    Support periodic curves.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math

from System import Enum


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bNumResForDist_NotAngle'; keys.append(key)
    values[key] = True
    names[key] = 'NumberEntry'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'DraftAngle', 'Dist')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistance'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bProjDist'; keys.append(key)
    values[key] = True
    names[key] = 'DistType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'True', 'Projected')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTaper_Start_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'TaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bVariableTaper'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTaper_End_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'EndTaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTaperChangePerCrvParam'; keys.append(key)
    values[key] = False
    names[key] = 'TaperChangePerCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Length', 'Param')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCPlane'; keys.append(key)
    values[key] = True
    names[key] = 'DirPerZAxisOf'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitPathsAtG2PlusKnots'; keys.append(key)
    values[key] = False
    names[key] = 'OnlyBezierOut'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAlignEndDirs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtGrevilles'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtEqualDivisions'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDivisionCt'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(
            initialValue=values[key],
            setLowerLimit=True,
            limit=2)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddTaperStartCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddTaperedLines'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddTaperEndCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddBrep'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iBrepMethod'; keys.append(key)
    listValues[key] = 'Loft2Crvs', 'LoftSectionLines', 'Sweep2A', 'Sweep2B', 'Network' # All items must be strings.
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftType'; keys.append(key)
    listValues[key] = Enum.GetNames(rg.LoftType)
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fBrepTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=1e-6)
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
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
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

        if key == 'fDistance':
            if -1e-6 < cls.riOpts[key].CurrentValue < 1e-6:
                print("fDistance input is too small.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fBrepTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Why is key, {}, here?  Value was not set or sticky-saved.".format(key))
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def _addCommonOptions(go):
    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    addOption('bNumResForDist_NotAngle')
    addOption('fDistance')
    if not Opts.values['bVariableTaper']:
        addOption('bProjDist')
    addOption('fTaper_Start_Deg')
    addOption('bVariableTaper')
    if Opts.values['bVariableTaper']:
        addOption('fTaper_End_Deg')
        idxs_Opt['SwapAngles'] = go.AddOption('SwapAngles')
        addOption('bTaperChangePerCrvParam')
    idxs_Opt['FlipDir'] = go.AddOption('FlipDir')
    idxs_Opt['FlipAngle'] = go.AddOption('FlipAngle')
    addOption('bCPlane')
    addOption('bSplitPathsAtG2PlusKnots')
    addOption('bAlignEndDirs')
    addOption('bAtGrevilles')
    addOption('bAtKnots')
    addOption('bAtEqualDivisions')
    if Opts.values['bAtEqualDivisions']:
        addOption('iDivisionCt')
    addOption('bAddTaperStartCrvs')
    addOption('bAddTaperedLines')
    addOption('bAddTaperEndCrvs')
    addOption('bAddBrep')
    if Opts.values['bAddBrep']:
        addOption('iBrepMethod')
        if Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'LoftSectionLines':
            addOption('iLoftType')
        addOption('fBrepTol')
    addOption('bEcho')
    addOption('bDebug')

    return idxs_Opt


def _getInput_Crvs():
    """
    Get objects with optional input.
    
    Returns
        None to cancel.
        ObjRefs[]
    """


    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select curves to extrude")
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=False)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()
        idxs_Opt = _addCommonOptions(go)

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return 

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fDistance' if Opts.values['bNumResForDist_NotAngle'] else 'fTaper_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
        elif Opts.values['bVariableTaper'] and go.OptionIndex() == idxs_Opt['SwapAngles']:
            Opts.riOpts['fTaper_Start_Deg'].CurrentValue, Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                    Opts.riOpts['fTaper_End_Deg'].CurrentValue, Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
            Opts.setValue('fTaper_Start_Deg')
            Opts.setValue('fTaper_End_Deg')
        elif go.OptionIndex() == idxs_Opt['FlipAngle']:
            for key in 'fTaper_Start_Deg', 'fTaper_End_Deg':
                Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
                Opts.setValue(key)
        elif go.OptionIndex() == idxs_Opt['FlipDir']:
            for key in 'fDistance', 'fTaper_Start_Deg', 'fTaper_End_Deg':
                Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
                Opts.setValue(key)
        elif Opts.values['bAddBrep'] and go.OptionIndex() == idxs_Opt['iBrepMethod']:
            Opts.setValue('iBrepMethod', go.Option().CurrentListOptionIndex)
            #Opts.values['iBrepMethod'] = go.Option().CurrentListOptionIndex
        elif Opts.values['bAddBrep'] and Opts.values['iBrepMethod'] == 1 and go.OptionIndex() == idxs_Opt['iLoftType']:
            #Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
            Opts.setValue('iLoftType', go.Option().CurrentListOptionIndex)
        else:
            for key in idxs_Opt:
                if go.Option().Index == idxs_Opt[key]:
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)
                    break


def _getInput_Opts():
    """
    Get objects with optional input.
    
    Returns
        None to cancel.
        bool(Generate new geometry), bool(Accept results):
    """


    gp = ri.Custom.GetPoint()

    bDefaultsForClickDirFlipOnly = (
        (not Opts.values['bVariableTaper'] and
            Opts.values['fTaper_Start_Deg'] in (0.0, 90.0, 180.0, 270.0))
        or
        (Opts.values['bVariableTaper'] and
            Opts.values['fTaper_Start_Deg'] in (0.0, 90.0, 180.0, 270.0) and
            Opts.values['fTaper_End_Deg'] == Opts.values['fTaper_Start_Deg'])
        )

    if bDefaultsForClickDirFlipOnly:
        gp.SetCommandPrompt("Left click to flip direction")
    else:
        gp.SetCommandPrompt("Left click to cycle 4 angle and direction combos")


    gp.SetCommandPromptDefault("Accept results")
    gp.AcceptNothing(True)
    gp.AcceptNumber(True, acceptZero=False)

    idxs_Opt = _addCommonOptions(gp)


    res = gp.Get()

    if res == ri.GetResult.Cancel:
        gp.Dispose()
        return 

    if res == ri.GetResult.Nothing:
        gp.Dispose()
        return False, True

    #sEval = "bDefaultsForClickDirFlipOnly"; print("{}: {}".format(sEval, eval(sEval)))
    #sEval = "Opts.values['bVariableTaper']"; print("{}: {}".format(sEval, eval(sEval)))
    #sEval = "Opts.values['fDistance']"; print("{}: {}".format(sEval, eval(sEval)))
    #sEval = "Opts.values['fTaper_Start_Deg']"; print("{}: {}".format(sEval, eval(sEval)))
    #sEval = "Opts.values['fTaper_End_Deg']"; print("{}: {}".format(sEval, eval(sEval)))


    if res == ri.GetResult.Point:
        if bDefaultsForClickDirFlipOnly:
            Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
            Opts.setValue('fDistance')
        else:
            if (
                Opts.riOpts['fDistance'].CurrentValue > 0.0 and
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue >= 0.0
            ):
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
                Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
                Opts.setValue('fTaper_Start_Deg')
                Opts.setValue('fTaper_End_Deg')
            elif (
                Opts.riOpts['fDistance'].CurrentValue < 0.0 and
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue >= 0.0
            ):
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
                Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
                Opts.setValue('fTaper_Start_Deg')
                Opts.setValue('fTaper_End_Deg')
            elif (
                Opts.riOpts['fDistance'].CurrentValue < 0.0 and
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue <= 0.0
            ):
                Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
                Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
                Opts.setValue('fDistance')
                Opts.setValue('fTaper_Start_Deg')
                Opts.setValue('fTaper_End_Deg')
            elif (
                Opts.riOpts['fDistance'].CurrentValue > 0.0 and
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue <= 0.0
            ):
                Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
                Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
                Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
                Opts.setValue('fDistance')
                Opts.setValue('fTaper_Start_Deg')
                Opts.setValue('fTaper_End_Deg')
            else:
                raise Exception("What happened?")
        gp.Dispose()
        return True, False

    if res == ri.GetResult.Number:
        key = 'fDistance' if Opts.values['bNumResForDist_NotAngle'] else 'fTaper_Start_Deg'
        Opts.riOpts[key].CurrentValue = gp.Number()
        Opts.setValue(key)
    elif Opts.values['bVariableTaper'] and gp.OptionIndex() == idxs_Opt['SwapAngles']:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue, Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                Opts.riOpts['fTaper_End_Deg'].CurrentValue, Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
        Opts.setValue('fTaper_Start_Deg')
        Opts.setValue('fTaper_End_Deg')
    elif gp.OptionIndex() == idxs_Opt['FlipAngle']:
        for key in 'fTaper_Start_Deg', 'fTaper_End_Deg':
            Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
            Opts.setValue(key)
    elif gp.OptionIndex() == idxs_Opt['FlipDir']:
        for key in 'fDistance', 'fTaper_Start_Deg', 'fTaper_End_Deg':
            Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
            Opts.setValue(key)
    elif Opts.values['bAddBrep'] and gp.OptionIndex() == idxs_Opt['iBrepMethod']:
        Opts.setValue('iBrepMethod', gp.Option().CurrentListOptionIndex)
        #Opts.values['iBrepMethod'] = go.Option().CurrentListOptionIndex
    elif Opts.values['bAddBrep'] and Opts.values['iBrepMethod'] == 1 and gp.OptionIndex() == idxs_Opt['iLoftType']:
        #Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
        Opts.setValue('iLoftType', gp.Option().CurrentListOptionIndex)
    else:
        for key in idxs_Opt:
            if gp.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, gp.Option().CurrentListOptionIndex)
                break

    return True, False


def continuityVectorsAt(nc, t, side=rg.CurveEvaluationSide.Default):
    """
    Returns: Tuple of these 4 items:
        3D point as a vector,
        Unit tangent vector,
        Curvature vector,
        Vector for comparing G3 continuity, not the G3 vector itself
        
        None for any of aforementioned vectors on fail.
    """
    
    if not isinstance(nc, rg.NurbsCurve): return
    
    vs = nc.DerivativeAt(
            t,
            derivativeCount=3,
            side=side)
    
    # Not using rg.Curve.TangentAt since it doesn't take into account CurveEvaluationSide.
    vTangency = vs[1]/vs[1].Length

    cross = rg.Vector3d.CrossProduct

    # For R3
    vCurvature = cross(cross(vs[1], vs[2]), vs[1]) / vs[1].Length**4

    return vs[0], vTangency, vCurvature


def getG2Discontinuities(rgCrv_In):
    """
    """

    fAngleTol_Rad = sc.doc.ModelAngleToleranceRadians
    fRadiusTol = sc.doc.ModelAbsoluteTolerance

    nc = rgCrv_In.ToNurbsCurve()

    t0 = nc.Domain.T0
    t1 = nc.Domain.T1

    ts_discontinuities = []

    bG3_discontinuousFound = False

    # This will also skip simple knot overlaps of Periodic curves.
    if nc.IsClosed and not nc.IsPeriodic:
        iK = 0
    else:
        iK = nc.Degree
    iK_Stop = nc.Knots.Count - nc.Degree

    while iK < iK_Stop:
        sc.escape_test()

        m = nc.Knots.KnotMultiplicity(iK)

        if m <= nc.Degree - 3:
            # Continuity at this knot is at least G3.
            iK += m
            continue

        if iK == 0:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                ) = continuityVectorsAt(
                    nc,
                    nc.Knots[nc.Knots.Count - 1],
                    rg.CurveEvaluationSide.Below)
        else:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                ) = continuityVectorsAt(
                    nc,
                    nc.Knots[iK],
                    rg.CurveEvaluationSide.Below)


        (
            vG0_Above,
            vG1_Above,
            vG2_Above,
            ) = continuityVectorsAt(
                nc,
                nc.Knots[iK],
                rg.CurveEvaluationSide.Above)


        # G1.

        if m > nc.Degree - 1:

            iParallel = vG1_Below.IsParallelTo(
                other=vG1_Above, angleTolerance=fAngleTol_Rad)

            if not iParallel == 1:
                ts_discontinuities.append(nc.Knots[iK])
                iK += m
                continue


        # G2.

        if m > nc.Degree - 2:

            if vG2_Below.IsTiny() and vG2_Above.IsTiny():
                # Linear.
                pass
            else:

                iParallel = vG2_Below.IsParallelTo(
                    other=vG2_Above, angleTolerance=fAngleTol_Rad)

                fCrvDelta_Degrees = Rhino.RhinoMath.ToDegrees(
                    rg.Vector3d.VectorAngle(vG2_Below, vG2_Above))

                if iParallel != 1:
                    ts_discontinuities.append(nc.Knots[iK])
                    iK += m
                    continue
                else:
                    kappa_Below = vG2_Below.Length
                    kappa_Above = vG2_Above.Length
                    ratio_of_curvature = (abs(kappa_Below-kappa_Above) /
                                        max(kappa_Below, kappa_Above))
                    delta_radius = abs(1.0/kappa_Below-1.0/kappa_Above)
                    if delta_radius > fRadiusTol:
                        ts_discontinuities.append(nc.Knots[iK])
                        iK += m
                        continue

        iK += m

    nc.Dispose()

    return ts_discontinuities


def tallyCurveTypes(crvs):
    sCrvTypes = []
    for crv in crvs:
        sCrvTypes.append(crv.GetType().Name)
    sOuts = []
    for sCrvType in sorted(set(sCrvTypes)):
        sOuts.append("{} of {}".format(sCrvTypes.count(sCrvType), sCrvType))
    return sOuts


def _prepareCurves(rgCrvs_In, bSplitPathsAtG2PlusKnots=False, bMakeDeg_1_Deformable=True, bMakeDeg_2_Deformable=True):
    """
    Prepare curves, including splitting per options.

    Returns: list of (profile) lists of new curves
    """

    # JoinCurves groups curves for further processing
    # and also aligns direction within each PolyCurve.
    if len(rgCrvs_In) == 1:
        rgCrvs_Joined = [rgCrvs_In[0].DuplicateCurve()]
    else:
        rgCrvs_Joined = rg.Curve.JoinCurves(
            rgCrvs_In,
            joinTolerance=sc.doc.ModelAbsoluteTolerance)
        if not rgCrvs_Joined: return

    rgCs_Out = [] # Connected paths.

    for rgCrv_Joined in rgCrvs_Joined:
        segs = rgCrv_Joined.DuplicateSegments()

        rgCrvs_ExplodedPoly = []

        for seg in segs:
            nc = seg.ToNurbsCurve()
            seg.Dispose()

            if nc.Degree != 2:
                rgCrvs_ExplodedPoly.append(nc)
                continue

            if nc.SpanCount == 1:
                rgCrvs_ExplodedPoly.append(nc)
                continue

            bcs = rg.BezierCurve.CreateBeziers(nc)
            for bc in bcs:
                rgCrvs_ExplodedPoly.append(bc.ToNurbsCurve())
        
        #print("  ".join(tallyCurveTypes(rgCrvs_ExplodedPoly)))

        rgCrv_Joined.Dispose()

        # Split each curve at each G2 discontinuity.
        rgCs_SplitAtNonG2_AllPolycrvs = []
        for rgCrv_SplitPoly in rgCrvs_ExplodedPoly:
            ts = getG2Discontinuities(rgCrv_SplitPoly)
            rgCs_SplitAtNonG2_ThisPolycrv = rgCrv_SplitPoly.Split(
                [rgCrv_SplitPoly.Domain.T0] + ts + [rgCrv_SplitPoly.Domain.T1])
            rgCs_SplitAtNonG2_AllPolycrvs.extend(rgCs_SplitAtNonG2_ThisPolycrv)

        if not bSplitPathsAtG2PlusKnots:
            rgCs_AfterG2PlusSplits = rgCs_SplitAtNonG2_AllPolycrvs[:] # Referencing same Geometry, so do not Dispose yet.
        else:
            rgCs_AfterG2PlusSplits = []
            for c in rgCs_SplitAtNonG2_AllPolycrvs:
                ts_SpanBoundaries = [c.Domain.T0]
                for iSpan in xrange(c.SpanCount):
                    ts_SpanBoundaries.append(c.SpanDomain(iSpan).T1)
                rgCrvs_SplitAtKnots = c.Split(ts_SpanBoundaries)
                rgCs_AfterG2PlusSplits.extend(rgCrvs_SplitAtKnots)
                c.Dispose()

        rgCs_Out_1Profile = []

        if not (bMakeDeg_1_Deformable or bMakeDeg_2_Deformable):
            rgCs_Out_1Profile.extend(rgCs_AfterG2PlusSplits)
        else:
            for c in rgCs_AfterG2PlusSplits:
                if bMakeDeg_1_Deformable and c.Points.Count == 2:
                    if not c.IncreaseDegree(3):
                        print("Could not raise degree of {}.".format(c.GetType().Name))
                elif bMakeDeg_2_Deformable and c.Points.Count == 3:
                    if not c.IncreaseDegree(3):
                        print("Could not raise degree of {}.".format(c.GetType().Name))
                rgCs_Out_1Profile.append(c)

        if len(rgCs_Out_1Profile) == 1:
            # Besides unnecessarily using JoinCurves,
            # JoinCurves will also convert a degree-1 NurbsCurve into a PolylineCurve.
            if rgCs_Out_1Profile[0].IsPeriodic:
                print("Periodic path profile skipped.")
                continue
            rgCs_Out.append(rgCs_Out_1Profile[0])
            continue

        rgCs_Joined_Out = rg.Curve.JoinCurves(
            rgCs_Out_1Profile,
            joinTolerance=sc.doc.ModelAbsoluteTolerance)
        if len(rgCs_Joined_Out) != 1:
            raise Exception("JoinCurves should have returned 1 curves, but returned {} instead.".format(len(rgCs_Joined_Out)))

        rgCs_Out.append(rgCs_Joined_Out[0])

    return rgCs_Out


def _getGrevilleParametersWithinDomain(nc):
    ts = list(nc.GrevilleParameters())
    #print(nc.Domain)
    #print(ts)
    if nc.IsClosed:
        while ts[0] < nc.Domain.T0:
            sc.escape_test()
            del ts[0]
        while ts[-1] >= nc.Domain.T1:
            sc.escape_test()
            del ts[-1]
    #print(ts)
    return ts


def _createArrayedGeometry(rgCrv_Path, rgLc_ToArray, plane_Proj, fTaper_Start_Deg, fTaper_End_Deg, bTaperChangePerCrvParam, bAtGrevilles, bAtKnots, bAtEqualDivisions, iDivisionCt=None, bDebug=False):
    """
    """

    xform_Proj = rg.Transform.PlanarProjection(plane_Proj)

    rgLcs_Arrayed = []

    nc_PathProfile = rgCrv_Path.ToNurbsCurve()
    if nc_PathProfile is None:
        print("NurbsCurve could not be calculated from curve.")
        return
    #sc.doc.Objects.AddCurve(nc2_Path)
        
    ts = []

    if bAtGrevilles:
        ts.extend(_getGrevilleParametersWithinDomain(nc_PathProfile))

    if bAtEqualDivisions:
        rc = nc_PathProfile.DivideByCount(
                segmentCount=iDivisionCt,
                includeEnds=True)
        if rc: ts.extend(rc)
        if nc_PathProfile.IsClosed:
            ts.append(nc_PathProfile.Domain.T1)
    
    if bAtKnots:
        rc = set(nc_PathProfile.Knots)
        if rc: ts.extend(rc)
        
    if ts is None:
        print("No parameters were obtained.")
        return

    ts = sorted(set(ts)) # Remove duplicates and sort.

    # Remove overlaps for closed (including periodic) curves.
    if nc_PathProfile.IsClosed:
        if bDebug:
            print(nc_PathProfile.Domain)
            print(ts)
        ts_WIP = []
        for t in ts:
            if t >= nc_PathProfile.Domain.T0 and t < nc_PathProfile.Domain.T1:
                ts_WIP.append(t)
        ts = ts_WIP
        if bDebug: print(ts)

    if fTaper_Start_Deg == fTaper_End_Deg:
        angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
    elif not bTaperChangePerCrvParam:
        length_Full = rgCrv_Path.GetLength()

    rgCrv_Path_Flattened = nc_PathProfile.Duplicate()
    rgCrv_Path_Flattened.Transform(xform_Proj)

    #sc.doc.Objects.AddCurve(rgCrv_Path_Flattened); return
    rgLcs_Arrayed = []
    for iT, t in enumerate(ts):
        bSuccess, frame = rgCrv_Path_Flattened.PerpendicularFrameAt(t=t)
        if not bSuccess:
            print("Perpendicular frame could not be calculated.")
            continue
                
        angle_StraightenFrame_Rad = rg.Vector3d.VectorAngle(frame.YAxis, plane_Proj.ZAxis, frame)

        if fTaper_Start_Deg != fTaper_End_Deg:
            if t == rgCrv_Path.Domain.T0:
                angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
            elif t == rgCrv_Path.Domain.T1:
                angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_End_Deg)
            #elif len(ts) >= 4 and iT == 1:
            #    # Same angle so that the brep starts at a tangent.
            #    angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
            #elif len(ts) >= 4 and iT == len(ts) - 2:
            #    # Same angle so that the brep end at a tangent.
            #    angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_End_Deg)

            else:
                if bTaperChangePerCrvParam:
                    t_Normalized = rgCrv_Path.Domain.NormalizedParameterAt(t)
                    angle_Rad = Rhino.RhinoMath.ToRadians(
                            fTaper_Start_Deg * (1.0 - t_Normalized) +
                            fTaper_End_Deg * t_Normalized)
                else:
                    length_to_t = rgCrv_Path.GetLength(
                            subdomain=rg.Interval(rgCrv_Path.Domain.T0, t))
                    angle_Rad = Rhino.RhinoMath.ToRadians(
                            fTaper_Start_Deg * (1.0 - length_to_t/length_Full) +
                            fTaper_End_Deg * length_to_t/length_Full)


            if bDebug: print(angle_Rad)
        frame.Rotate(angle=angle_StraightenFrame_Rad+angle_Rad, axis=frame.ZAxis)

        # Debug frame orientation.
        #                sc.doc.Objects.AddPoint(frame.PointAt(0.0,0.0,0.0))
        #                attr_Red = rd.ObjectAttributes()
        #                attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
        #                attr_Red.ObjectColor = Color.Red
        #                sc.doc.Objects.AddPoint(frame.PointAt(1.0,0.0,0.0), attr_Red)
        #                attr_Green = rd.ObjectAttributes()
        #                attr_Green.ColorSource = rd.ObjectColorSource.ColorFromObject
        #                attr_Green.ObjectColor = Color.Lime
        #                sc.doc.Objects.AddPoint(frame.PointAt(0.0,1.0,0.0), attr_Green)
                
        xform1 = rg.Transform.PlaneToPlane(plane_Proj, frame)
        xform2 = rg.Transform.Translation(
                nc_PathProfile.PointAt(t)-rgCrv_Path_Flattened.PointAt(t))
        xform3 = xform2 * xform1

        rgLc_Arrayed_WIP = rgLc_ToArray.Duplicate()
        rgLc_Arrayed_WIP.Transform(xform3)
        rgLcs_Arrayed.append(rgLc_Arrayed_WIP)

    nc_PathProfile.Dispose()
    rgCrv_Path_Flattened.Dispose()

    return rgLcs_Arrayed, ts


def _matchCrvEndDirs(nc_ToMod, nc_Ref):
    """
    nc_ToMod is modified.
    """
    
    if nc_ToMod.Points.Count < 4:
        return False
    
    bSuccess = nc_ToMod.SetEndCondition(
        bSetEnd=False,
        continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
        point=nc_ToMod.PointAtStart,
        tangent=nc_Ref.TangentAtStart)
    if not bSuccess:
        print("SetEndCondition failed.")
        return False
    bSuccess = nc_ToMod.SetEndCondition(
        bSetEnd=True,
        continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
        point=nc_ToMod.PointAtEnd,
        tangent=nc_Ref.TangentAtEnd)
    if not bSuccess:
        print("SetEndCondition failed.")
        return False

    return True


def createBrep(iBrepMethod, iLoftType, fBrepTol, rgCrv_Path, rgNurbsCrv_TaperEnd_1Seg, rgLineCrvs_Arrayed):
    """
    """

    rgNurbsCrv1_PathSeg = rgCrv_Path.ToNurbsCurve()

    rgBs_Out = []
    
    if Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Loft2Crvs':
        rgBs_Out = rg.Brep.CreateFromLoft(
                curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg],
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=rg.LoftType.Straight,
                closed=False)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'LoftSectionLines':
        for L in rgLineCrvs_Arrayed:
            print(L.PointAtStart, L.PointAtEnd)
        print(rgLineCrvs_Arrayed[0].PointAtEnd.EpsilonEquals(rgLineCrvs_Arrayed[-1].PointAtEnd, epsilon=1e-12))
        rgBs_Out = rg.Brep.CreateFromLoft(
                curves=rgLineCrvs_Arrayed,
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=Enum.ToObject(rg.LoftType, iLoftType),
                closed=rgNurbsCrv1_PathSeg.IsClosed)
        print(rgBs_Out)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Sweep2A':
        rgBs_Out = rg.Brep.CreateFromSweep(
                rail1=rgNurbsCrv1_PathSeg,
                rail2=rgNurbsCrv_TaperEnd_1Seg,
                shapes=rgLineCrvs_Arrayed,
                closed=rgNurbsCrv1_PathSeg.IsClosed,
                tolerance=fBrepTol)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Sweep2B':
        rgSweep2 = rg.SweepTwoRail()
        #rgSweep2.AngleToleranceRadians
        rgSweep2.ClosedSweep = rgNurbsCrv1_PathSeg.IsClosed
        rgSweep2.MaintainHeight = False
        rgSweep2.SweepTolerance = fBrepTol
        rgBs_Out = rgSweep2.PerformSweep(
                rail1=rgNurbsCrv1_PathSeg,
                rail2=rgNurbsCrv_TaperEnd_1Seg,
                crossSections=rgLineCrvs_Arrayed)
    elif Opts.listValues['iBrepMethod'][Opts.values['iBrepMethod']] == 'Network':
        rgNurbsSrf, iError = rg.NurbsSurface.CreateNetworkSurface(
                curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg]+rgLineCrvs_Arrayed,
                continuity=1,
                edgeTolerance=fBrepTol,
                interiorTolerance=fBrepTol,
                angleTolerance=0.1*sc.doc.ModelAngleToleranceDegrees)
        if iError:
            print("CreateNetworkSurface error code: {}".format(iError))
        else:
            rgBs_Out = [rgNurbsSrf.ToBrep()]
            rgNurbsSrf.Dispose()

    rgNurbsCrv1_PathSeg.Dispose()

    return rgBs_Out


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.breps = []
        self.crvs = []
        self.lines = []
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        for brep in self.breps:
            bbox = brep.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

        for crv in self.crvs:
            bbox = crv.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

        for line in self.lines:
            bbox = line.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PreDrawObjects(self, drawEventArgs):

        breps = self.breps
        color = sc.doc.Layers.CurrentLayer.Color

        for brep in breps:
            displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
            da = displayMode.DisplayAttributes
            if da.ShadingEnabled:
                drawEventArgs.Display.DrawBrepShaded(
                    brep=brep,
                    material=Rhino.Display.DisplayMaterial(diffuse=color))
            drawEventArgs.Display.DrawBrepWires(
                brep=brep,
                color=color,
                wireDensity=1)

        for crv in self.crvs:
            drawEventArgs.Display.DrawCurve(
                curve=crv,
                color=color,
                thickness=self.crv_thk)

        if self.lines:
            drawEventArgs.Display.DrawLines(
                lines=self.lines,
                color=color,
                thickness=self.crv_thk)

    def clearGeometry(self):
        for rgOs in self.breps, self.crvs, self.lines:
            if rgOs:
                for o in rgOs: o.Dispose()
            list
            rgOs.Clear() # .NET method.
        #print("All conduit's geometry are disposed and cleared.")


def main():

    sk_conduit = 'conduit({})'.format(__file__)

    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.clearGeometry()
        conduit.Enabled = False
        sc.doc.Views.Redraw()
        if Opts.values['bDebug']:
            conduit = None
            conduit = DrawConduit()
            sc.sticky[sk_conduit] = conduit
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit


    objrefs_Paths = _getInput_Crvs()
    if objrefs_Paths is None: return


    rgCs_Paths_In = []
    for objref_Path in objrefs_Paths:
        c = objref_Path.Curve()
        rgCs_Paths_In.append(c)

    bGenerateNew = True
    bAcceptResults = False

    fDistance = Opts.values['fDistance']
    bProjDist = Opts.values['bProjDist']
    fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
    bVariableTaper = Opts.values['bVariableTaper']
    fTaper_End_Deg = Opts.values['fTaper_End_Deg']
    bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
    bAlignEndDirs = Opts.values['bAlignEndDirs']
    bCPlane = Opts.values['bCPlane']
    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtKnots = Opts.values['bAtKnots']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bSplitPathsAtG2PlusKnots = Opts.values['bSplitPathsAtG2PlusKnots']
    bAddTaperStartCrvs = Opts.values['bAddTaperStartCrvs']
    bAddTaperedLines = Opts.values['bAddTaperedLines']
    bAddTaperEndCrvs = Opts.values['bAddTaperEndCrvs']
    bAddBrep = Opts.values['bAddBrep']
    iBrepMethod = Opts.values['iBrepMethod']
    iLoftType = Opts.values['iLoftType']
    fBrepTol = Opts.values['fBrepTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    # Only setting these variables for code that checks them.
    rgLine_ToArray = None
    rgCs_Paths_Prepared = []


    while True:

        while not (bAtGrevilles or bAtKnots or bAtEqualDivisions):
            print("No path point sampling is enabled.")
            sc.doc.Views.Redraw()

            rc = _getInput_Opts()
            if rc is None:
                if rgLine_ToArray: rgLine_ToArray.Dispose()
                for c in rgCs_Paths_Prepared: c.Dispose()
                return

            bGenerateNew, bAcceptResults = rc

            if bAcceptResults:
                if rgLine_ToArray: rgLine_ToArray.Dispose()
                for c in rgCs_Paths_Prepared: c.Dispose()
                return

            fDistance = Opts.values['fDistance']
            bProjDist = Opts.values['bProjDist']
            fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
            bVariableTaper = Opts.values['bVariableTaper']
            fTaper_End_Deg = Opts.values['fTaper_End_Deg']
            bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
            bAlignEndDirs = Opts.values['bAlignEndDirs']
            bCPlane = Opts.values['bCPlane']
            bAtGrevilles = Opts.values['bAtGrevilles']
            bAtKnots = Opts.values['bAtKnots']
            bAtEqualDivisions = Opts.values['bAtEqualDivisions']
            iDivisionCt = Opts.values['iDivisionCt']
            bSplitPathsAtG2PlusKnots = Opts.values['bSplitPathsAtG2PlusKnots']
            bAddTaperStartCrvs = Opts.values['bAddTaperStartCrvs']
            bAddTaperedLines = Opts.values['bAddTaperedLines']
            bAddTaperEndCrvs = Opts.values['bAddTaperEndCrvs']
            bAddBrep = Opts.values['bAddBrep']
            iBrepMethod = Opts.values['iBrepMethod']
            iLoftType = Opts.values['iLoftType']
            fBrepTol = Opts.values['fBrepTol']
            bEcho = Opts.values['bEcho']
            bDebug = Opts.values['bDebug']


        if bAlignEndDirs:
            bMakeDeg_1_Deformable = fTaper_End_Deg and (fTaper_End_Deg != fTaper_Start_Deg)
            bMakeDeg_2_Deformable = True
        else:
            bMakeDeg_1_Deformable = False
            bMakeDeg_2_Deformable = False

        Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")

        rgCs_Paths_Prepared = _prepareCurves(
            rgCrvs_In=rgCs_Paths_In,
            bSplitPathsAtG2PlusKnots=bSplitPathsAtG2PlusKnots,
            bMakeDeg_1_Deformable=bMakeDeg_1_Deformable,
            bMakeDeg_2_Deformable=bMakeDeg_2_Deformable)
        if not rgCs_Paths_Prepared:
            print("No valid path curves were provided.")
            conduit = None
            return


        if bProjDist and not bVariableTaper:
            rgLine_ToArray = rg.Line(
                rg.Point3d(0.0, 0.0, 0.0),
                rg.Point3d(
                    0.0,
                    fDistance/math.cos(math.radians(fTaper_Start_Deg)),
                    0.0))
        else:
            rgLine_ToArray = rg.Line(
                rg.Point3d(0.0, 0.0, 0.0),
                rg.Point3d(0.0, fDistance, 0.0))

        if bCPlane:
            view_Active = sc.doc.Views.ActiveView
            plane_Proj = view_Active.ActiveViewport.ConstructionPlane()
        
            xform1 = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, plane_Proj)
            rgLine_ToArray.Transform(xform1)
        else:
            plane_Proj = rg.Plane.WorldXY

        rgLineCrv_ToArray = rg.LineCurve(rgLine_ToArray)

        rgCs_TaperStart_All = []
        rgLcs_Arrayed_All = []
        rgCs_TaperEnd_All = []
        rgBs_Res_All = []

        Rhino.RhinoApp.SetCommandPrompt("Creating geometry ...")

        for i_Path_Profile, rgC_Path_Prepared in enumerate(rgCs_Paths_Prepared):

            rgBs_1Profile = []

            if bVariableTaper:
                bDisableVariableTaper = rgC_Path_Prepared.IsClosed and fTaper_End_Deg != fTaper_Start_Deg
                print("Variable taper disabled for closed path profile.")
            else:
                bDisableVariableTaper = False

            # These arrayed lines are per options.
            # Since they are of entire path profiles, they won't contain duplicates.
            rc = _createArrayedGeometry(
                rgCrv_Path=rgC_Path_Prepared,
                rgLc_ToArray=rgLineCrv_ToArray,
                plane_Proj=plane_Proj,
                fTaper_Start_Deg=fTaper_Start_Deg,
                fTaper_End_Deg=fTaper_End_Deg if (bVariableTaper and not bDisableVariableTaper) else fTaper_Start_Deg,
                bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                bAtGrevilles=bAtGrevilles,
                bAtKnots=bAtKnots,
                bAtEqualDivisions=bAtEqualDivisions,
                iDivisionCt=iDivisionCt,
                bDebug=bDebug)
            if rc is None:
                raise Exception("Nothing returned from createArrayedGeometry.")

            rgLcs_Arrayed_PerPathProfile, ts_Lcs_Arrayed_PerProfile = rc

            if bAddTaperedLines:
                rgLcs_Arrayed_All.extend(rgLcs_Arrayed_PerPathProfile)
            else:
                # bAddTaperedLines == False.
                if bAtKnots or bAtEqualDivisions:
                    if bAtKnots and iBrepMethod == 0:
                        s  = "AtKnots"
                        s += " only affects added arrayed lines"
                        s += " when BrepMethod == {},".format(Opts.listValues['iBrepMethod'][iBrepMethod])
                        s += " but AddArrayed option is disabled."
                        print(s)
                    if bAtEqualDivisions and iBrepMethod == 0:
                        s  = "AtEqualDivisions"
                        s += " only affects added arrayed lines"
                        s += " when BrepMethod == {},".format(Opts.listValues['iBrepMethod'][iBrepMethod])
                        s += " but AddArrayed option is disabled."
                        print(s)

            if bAddTaperEndCrvs or bAddBrep:
                # Taper end segments require both taper start segments and tapered lines.
                # Breps require taper start segments, tapered lines, and/or taper end segments,
                # depending on the surface-creation method used.


                segs = rgC_Path_Prepared.DuplicateSegments() # All are NurbsCurves.
                #print(rgC_Path_Prepared.Domain)
                #for seg in segs:
                #    print(seg.Domain)
                #1/0

                def getArrayedLineCrvsAtGrevillesOnly(rgC_Path_1Seg, rgLcs_In, ts_In):
                    rgLcs_Out = []
                    ts = _getGrevilleParametersWithinDomain(rgC_Path_1Seg)
                    for i, t in enumerate(ts):
                        if t in ts_In:
                            rgLcs_Out.append(rgLcs_In[ts_In.index(t)])
                        else:
                            raise Exception("Parameter, {}, is not in list: {}".format(t, ts_In))
                    return rgLcs_Out

                for i_Path_1Seg, rgCrv1_Path_1Seg in enumerate(segs):

                    rgCs_TaperStart_All.append(rgCrv1_Path_1Seg)

                    # These arrayed lines are only at Greville points.

                    rgLineCrvs_Arrayed_1PathSeg_GrevsOnly = getArrayedLineCrvsAtGrevillesOnly(
                        rgCrv1_Path_1Seg, 
                        rgLcs_Arrayed_PerPathProfile,
                        ts_Lcs_Arrayed_PerProfile)

                    #rgLineCrvs_Arrayed_1PathSeg_GrevsOnly = _createArrayedGeometry(
                    #    rgCrv_Path=rgCrv1_Path_1Seg,
                    #    rgLc_ToArray=rgLineCrv_ToArray,
                    #    plane_Proj=plane_Proj,
                    #    fTaper_Start_Deg=fTaper_Start_Deg,
                    #    fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                    #    bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                    #    bAtGrevilles=True,
                    #    bAtKnots=False,
                    #    bAtEqualDivisions=False,
                    #    iDivisionCt=0,
                    #    bDebug=bDebug)
                    #if rgLineCrvs_Arrayed_1PathSeg_GrevsOnly is None:
                    #    raise Exception("Nothing returned from createArrayedGeometry.")

                    # Create taper end curve.
                    pts_AtEndsOf_Lcs_Arrayed = []
                    for rgLc in rgLineCrvs_Arrayed_1PathSeg_GrevsOnly:
                        pts_AtEndsOf_Lcs_Arrayed.append(rgLc.PointAtEnd)
                    nc_OppOfPath = rgCrv1_Path_1Seg.DuplicateCurve()
                    nc_OppOfPath.SetGrevillePoints(pts_AtEndsOf_Lcs_Arrayed)

                    if bAlignEndDirs:
                        _matchCrvEndDirs(nc_OppOfPath, rgCrv1_Path_1Seg)


                    rgCs_TaperEnd_All.append(nc_OppOfPath)

                    if bAddBrep:
                        rc = createBrep(
                            iBrepMethod=iBrepMethod,
                            iLoftType=iLoftType,
                            fBrepTol=fBrepTol,
                            rgCrv_Path=rgCrv1_Path_1Seg,
                            rgNurbsCrv_TaperEnd_1Seg=nc_OppOfPath,
                            rgLineCrvs_Arrayed=rgLineCrvs_Arrayed_1PathSeg_GrevsOnly)
                        if rc is None:
                            print("Cannot create brep(s).  Check input.")
                        else:
                            # Save to later join, adding only the joined breps to the conduit.
                            rgBs_1Profile.extend(rc)

            elif bAddTaperStartCrvs:
                # And bAddTaperEndCrvs and bAddBrep are both False.
                rgCs_TaperStart_All = rgC_Path_Prepared.DuplicateSegments()

            if bAddBrep and rgBs_1Profile:
                rgBs_Joined = rg.Brep.JoinBreps(
                    rgBs_1Profile,
                    tolerance=2.0*fBrepTol)
                for b in rgBs_1Profile: b.Dispose()

                if rgBs_Joined:
                    rgBs_Res_All.extend(rgBs_Joined)


        if bAddTaperStartCrvs:
            conduit.crvs.extend(rgCs_TaperStart_All)

        if bAddTaperedLines:
            conduit.crvs.extend(rgLcs_Arrayed_All)

        if bAddTaperEndCrvs:
            conduit.crvs.extend(rgCs_TaperEnd_All)

        if bAddBrep:
            conduit.breps.extend(rgBs_Res_All)

        conduit.Enabled = True

        sc.doc.Views.Redraw()


        rc = _getInput_Opts()
        if rc is None:
            conduit.clearGeometry()
            conduit.Enabled = False
            sc.doc.Views.Redraw()
            return

        bGenerateNew, bAcceptResults = rc

        if bAcceptResults:
            break

        if not bGenerateNew:
            continue

        conduit.clearGeometry()
        conduit.Enabled = False
        sc.doc.Views.Redraw()

        fDistance = Opts.values['fDistance']
        bProjDist = Opts.values['bProjDist']
        fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
        bVariableTaper = Opts.values['bVariableTaper']
        fTaper_End_Deg = Opts.values['fTaper_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bCPlane = Opts.values['bCPlane']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPathsAtG2PlusKnots = Opts.values['bSplitPathsAtG2PlusKnots']
        bAddTaperStartCrvs = Opts.values['bAddTaperStartCrvs']
        bAddTaperedLines = Opts.values['bAddTaperedLines']
        bAddTaperEndCrvs = Opts.values['bAddTaperEndCrvs']
        bAddBrep = Opts.values['bAddBrep']
        iBrepMethod = Opts.values['iBrepMethod']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


    if bAddBrep:
        gBs_Out = []
        iCt_Faces_Added = 0
        for rgB in conduit.breps:
            gB_Out = sc.doc.Objects.AddBrep(rgB)
            if gB_Out != gB_Out.Empty:
                gBs_Out.append(gB_Out)
                if bEcho:
                    rdB_Out = sc.doc.Objects.FindId(gB_Out)
                    iCt_Faces_Added += rdB_Out.BrepGeometry.Faces.Count
        if bEcho:
            print("{} breps with {} faces created.".format(
                len(gBs_Out), iCt_Faces_Added))

    if bAddTaperedLines:
        gLineCrvs1_Arrayed = []
        for rgLc in rgLcs_Arrayed_All:
            gL = sc.doc.Objects.AddCurve(rgLc)
            if gL != gL.Empty:
                gLineCrvs1_Arrayed.append(gL)


    if bAddTaperStartCrvs:
        gCrvs_TaperStart_All = []
        for nc in rgCs_TaperStart_All:
            gCrv_TaperStart_1Seg = sc.doc.Objects.AddCurve(nc)
            if gCrv_TaperStart_1Seg != gCrv_TaperStart_1Seg.Empty:
                gCrvs_TaperStart_All.append(gCrv_TaperStart_1Seg)

    if bAddTaperEndCrvs:
        gCrvs_TaperEnd_All = []
        for nc in rgCs_TaperEnd_All:
            gCrv_TaperEnd_1Seg = sc.doc.Objects.AddCurve(nc)
            if gCrv_TaperEnd_1Seg != gCrv_TaperEnd_1Seg.Empty:
                gCrvs_TaperEnd_All.append(gCrv_TaperEnd_1Seg)

    conduit.clearGeometry()

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()