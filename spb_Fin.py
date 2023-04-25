"""
This script is an alternative to _Fin.

Routine SPB
    Calculates translated Greville points.
    Has options to specify an angle from normal, constant or variable.

Routine RC
    Uses RhinoCommon's OffsetNormalToSurface method and is more similar to _Fin.

Both methods creates a loft between the curve on surface and the offset curve.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190625-27: Created, starting with another script.
...
230419-20: Added use of OffsetNormalToSurface.  Refactored.
230423-25: Various major changes.  Added AlignEndDirs option.  Now uses DrawConduit for preview.

TODO:
    ? Replace _curveWithSpansCompletelyOnTheFace with a routine that instead trims curve to the face.
    Check accuracy of offset curve and add knots until within tolerance.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bUseFaceOfSelNakedEdge'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRC_NotSPB'; keys.append(key)
    values[key] = False
    names[key] = 'Routine'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'SPB', 'RC')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAlignEndDirs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExplodePolyCrv'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitAtPolyKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSimplifyCrv'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fAngle_Start_Deg'; keys.append(key)
    values[key] = 0.0
    names[key] = 'AngleFromNormal'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bVariableAngle'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fAngle_End_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'EndAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAngleChangePerCrvParam_NotLength'; keys.append(key)
    values[key] = False
    names[key] = 'AngleChangePerCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Length', 'Param')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistance'; keys.append(key)
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 1.0
    else:
        values[key] = 10.0 * Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTargetTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAddSrf'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftType'; keys.append(key)
    listValues[key] = Enum.GetNames(rg.LoftType) # All items must be strings.
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddFinEndCrv'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddLines'; keys.append(key)
    values[key] = False
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
    values[key] = 10
    riOpts[key] = ri.Custom.OptionInteger(
            initialValue=values[key], setLowerLimit=True, limit=1)
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fTargetTol':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

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


def _addCommonOptions(go):
    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    addOption('bRC_NotSPB')
    addOption('bAlignEndDirs')
    addOption('bExplodePolyCrv')
    addOption('bSplitAtPolyKnots')
    addOption('bSimplifyCrv')
    if not Opts.values['bRC_NotSPB']:
        addOption('fAngle_Start_Deg')
        addOption('bVariableAngle')
        if Opts.values['bVariableAngle']:
            addOption('fAngle_End_Deg')
            addOption('bAngleChangePerCrvParam_NotLength')
    addOption('fDistance')
    addOption('fTargetTol')
    addOption('bAddSrf')
    addOption('bAddFinEndCrv')
    addOption('bAddLines')
    if Opts.values['bAddLines']:
        addOption('bAtGrevilles')
        addOption('bAtKnots')
        addOption('bAtEqualDivisions')
        if Opts.values['bAtEqualDivisions']:
            addOption('iDivisionCt')
    addOption('bEcho')
    addOption('bDebug')

    return idxs_Opt


def _getInput_Curve():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve on face")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    #bPreselectedObjsChecked = False

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()
        addOption('bUseFaceOfSelNakedEdge')
        idxs_Opt.update(_addCommonOptions(go))

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref_CrvOnFace = go.Object(0)
            go.Dispose()

            return objref_CrvOnFace

        bDefaultsForDir_NotAngle = (
            Opts.values['bRC_NotSPB']
            or
            (not Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 180.0))
            or
            (not Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 180.0) and
                Opts.values['fAngle_End_Deg'] == Opts.values['fAngle_Start_Deg'])
            )

        if res == ri.GetResult.Number:
            if bDefaultsForDir_NotAngle:
                key = 'fDistance'
            else:
                key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        if 'iLoftType' in idxs_Opt:
            Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_Face():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select base face")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    go.AcceptNumber(True, acceptZero=True)
    
    idxs_Opt = {}

    bPreselectedObjsChecked = False

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()
        idxs_Opt.update(_addCommonOptions(go))

        res = go.Get()

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
            objref_Face = go.Object(0)

            go.Dispose()
    
            return objref_Face


        bDefaultsForDir_NotAngle = (
            Opts.values['bRC_NotSPB']
            or
            (not Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 180.0))
            or
            (not Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 180.0) and
                Opts.values['fAngle_End_Deg'] == Opts.values['fAngle_Start_Deg'])
            )


        if res == ri.GetResult.Number:
            if bDefaultsForDir_NotAngle:
                key = 'fDistance'
            else:
                key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        if Opts.values['bAddSrf'] and not Opts.values['bLoose'] and Opts.values['bUseSectionLines']:
            if (
                Opts.values['bLoft_NotSweep'] and
                go.OptionIndex() == idxs_Opt['iLoftType']
            ):
                Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
            elif Opts.riOpts['fTargetTol'].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                Opts.riOpts['fTargetTol'].CurrentValue = Opts.riOpts['fTargetTol'].InitialValue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_Click():
    """
    Click to toggle angle and/or direction with optional input.

    Returns:
        True: To recalculate and reloop
        False: To not recalculate and break out of loop with current output.
        None: To not recalculate and return without output.
    """

    go = ri.Custom.GetPoint()

    bDefaultsForDir_NotAngle = (
        Opts.values['bRC_NotSPB']
        or
        (not Opts.values['bVariableAngle'] and
            Opts.values['fAngle_Start_Deg'] in (0.0, 180.0))
        or
        (not Opts.values['bVariableAngle'] and
            Opts.values['fAngle_Start_Deg'] in (0.0, 180.0) and
            Opts.values['fAngle_End_Deg'] == Opts.values['fAngle_Start_Deg'])
        )

    if bDefaultsForDir_NotAngle:
        go.SetCommandPrompt("Left click to flip direction")
    else:
        go.SetCommandPrompt("Left click to flip angle")

    go.AcceptNumber(True, acceptZero=True)
    go.AcceptNothing(True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    key = 'FlipDir'; idxs_Opt[key] = go.AddOption(key)
    if not bDefaultsForDir_NotAngle and Opts.values['bVariableAngle']:
        key = 'FlipAngle'; idxs_Opt[key] = go.AddOption(key)
        key = 'SwapAngles'; idxs_Opt[key] = go.AddOption(key)

    idxs_Opt.update(_addCommonOptions(go))

    res = go.Get()

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Nothing:
        go.Dispose()
        return False

    if res == ri.GetResult.Point:
        if bDefaultsForDir_NotAngle:
            Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
            Opts.setValue('fDistance')
            #if Opts.riOpts['fAngle_Start_Deg'].CurrentValue == 0.0:
            #    Opts.riOpts['fAngle_Start_Deg'].CurrentValue = 180.0
            #elif Opts.riOpts['fAngle_Start_Deg'].CurrentValue < 0.0:
            #    Opts.riOpts['fAngle_Start_Deg'].CurrentValue = (
            #            -180.0 + Opts.riOpts['fAngle_Start_Deg'].CurrentValue) % 180.0
            #else:
            #    Opts.riOpts['fAngle_Start_Deg'].CurrentValue = (
            #            180.0 - Opts.riOpts['fAngle_Start_Deg'].CurrentValue) % 180.0
            
            #if Opts.riOpts['fAngle_End_Deg'].CurrentValue == 0.0:
            #    Opts.riOpts['fAngle_End_Deg'].CurrentValue = 180.0
            #elif Opts.riOpts['fAngle_End_Deg'].CurrentValue < 0.0:
            #    Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
            #            -180.0 + Opts.riOpts['fAngle_End_Deg'].CurrentValue) % 180.0
            #else:
            #    Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
            #            180.0 - Opts.riOpts['fAngle_End_Deg'].CurrentValue) % 180.0
        else:
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
            Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
            Opts.setValue('fAngle_Start_Deg')
            Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return True

    if res == ri.GetResult.Number:
        if bDefaultsForDir_NotAngle:
            key = 'fDistance'
        else:
            key = 'fAngle_Start_Deg'
        Opts.riOpts[key].CurrentValue = go.Number()
        Opts.setValue(key)
        go.Dispose()
        return True

    # An option was selected.

    if Opts.values['bVariableAngle'] and go.OptionIndex() == idxs_Opt['SwapAngles']:
        Opts.riOpts['fAngle_Start_Deg'].CurrentValue, Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
                Opts.riOpts['fAngle_End_Deg'].CurrentValue, Opts.riOpts['fAngle_Start_Deg'].CurrentValue)
        Opts.setValue('fAngle_Start_Deg')
        Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return True

    if not bDefaultsForDir_NotAngle and 'FlipAngle' in idxs_Opt and go.OptionIndex() == idxs_Opt['FlipAngle']:
        #Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
        #Opts.setValue('fDistance')
        Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
        Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
        Opts.setValue('fAngle_Start_Deg')
        Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return True

    if go.OptionIndex() == idxs_Opt['FlipDir']:
        Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
        Opts.setValue('fDistance')
    #    if Opts.riOpts['fAngle_Start_Deg'].CurrentValue == 0.0:
    #        Opts.riOpts['fAngle_Start_Deg'].CurrentValue = 180.0
    #    elif Opts.riOpts['fAngle_Start_Deg'].CurrentValue < 0.0:
    #        Opts.riOpts['fAngle_Start_Deg'].CurrentValue = (
    #                -180.0 + Opts.riOpts['fAngle_Start_Deg'].CurrentValue) % 180.0
    #    else:
    #        Opts.riOpts['fAngle_Start_Deg'].CurrentValue = (
    #                180.0 - Opts.riOpts['fAngle_Start_Deg'].CurrentValue) % 180.0
            
    #    if Opts.riOpts['fAngle_End_Deg'].CurrentValue == 0.0:
    #        Opts.riOpts['fAngle_End_Deg'].CurrentValue = 180.0
    #    elif Opts.riOpts['fAngle_End_Deg'].CurrentValue < 0.0:
    #        Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
    #                -180.0 + Opts.riOpts['fAngle_End_Deg'].CurrentValue) % 180.0
    #    else:
    #        Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
    #                180.0 - Opts.riOpts['fAngle_End_Deg'].CurrentValue) % 180.0

    #    Opts.setValue('fAngle_Start_Deg')
    #    Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return True

    if 'iLoftType' in idxs_Opt:
        Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex

    for key in idxs_Opt:
        if go.Option().Index == idxs_Opt[key]:
            Opts.setValue(key, go.Option().CurrentListOptionIndex)
            break

    go.Dispose()
    return True


def _createArrayedLines(rgCrv, rgSrf, fDistance, fAngle_Start_Deg, fAngle_End_Deg, rgCrv_Full, bAngleChangePerCrvParam_NotLength, bAtGrevilles, bAtKnots, bAtEqualDivisions, iDivisionCt):
    """
    rgSurface0: rg.Surface, including rg.BrepFace
    rgCrv_Full: For variable angle calculation
    """

    rgCrv_In = rgCrv
    rgSrf_In = rgSrf

    rgNurbsCrv_rgCrv0 = rgCrv_In.ToNurbsCurve()

    rgLines_Arrayed = []

    ts = []
    if bAtGrevilles:
        ts.extend(rgNurbsCrv_rgCrv0.GrevilleParameters())
    if bAtKnots:
        ts.append(rgNurbsCrv_rgCrv0.Domain.T0)
        for iSpan in xrange(rgNurbsCrv_rgCrv0.SpanCount):
            ts.append(rgNurbsCrv_rgCrv0.SpanDomain(iSpan).T1)
        ts.append(rgNurbsCrv_rgCrv0.Domain.T1)
    if bAtEqualDivisions:
        if rgNurbsCrv_rgCrv0.IsClosed: iDivisionCt -= 1
        rc = rgNurbsCrv_rgCrv0.DivideByCount(
                segmentCount=iDivisionCt,
                includeEnds=True)
        if rc:
            ts.extend(rc)

    ts = sorted(set(ts))
    # TODO: Remove near duplicate parameters.

    angle_Const_Rad = angle_Var_Rad = None
    if fAngle_End_Deg is None or fAngle_End_Deg == fAngle_Start_Deg:
        angle_Const_Rad = Rhino.RhinoMath.ToRadians(fAngle_Start_Deg)
    elif not bAngleChangePerCrvParam_NotLength:
        length_Full = rgCrv_Full.GetLength()

    for t in ts:
        pt_Start = rgNurbsCrv_rgCrv0.PointAt(t)
        bSuccess, u, v = rgSrf_In.ClosestPoint(pt_Start)
        if not bSuccess: continue
        vect_Normal = rgSrf_In.NormalAt(u, v)
        bSuccess, frame = rgNurbsCrv_rgCrv0.PerpendicularFrameAt(t)
        if not bSuccess: continue
        
        pt_NoAngle = pt_Start + vect_Normal * fDistance

        pt_End = rg.Point3d(pt_NoAngle)

        if angle_Const_Rad is None:
            if bAngleChangePerCrvParam_NotLength:
                t_Normalized = rgCrv_Full.Domain.NormalizedParameterAt(t)
                angle_Var_Rad = Rhino.RhinoMath.ToRadians(
                    fAngle_Start_Deg * (1.0 - t_Normalized) +
                    fAngle_End_Deg * t_Normalized)
            else:
                length_to_t = rgCrv_Full.GetLength(
                    subdomain=rg.Interval(rgCrv_Full.Domain.T0, t))
                angle_Var_Rad = Rhino.RhinoMath.ToRadians(
                    fAngle_Start_Deg * (1.0 - length_to_t/length_Full) +
                    fAngle_End_Deg * length_to_t/length_Full)

        xform_Rotation = rg.Transform.Rotation(
                angleRadians=angle_Const_Rad if angle_Var_Rad is None else angle_Var_Rad,
                rotationAxis=frame.ZAxis,
                rotationCenter=frame.Origin)

        pt_End.Transform(xform_Rotation)

        rgLines_Arrayed.append(rg.Line(pt_Start, pt_End))

    if rgLines_Arrayed[0].EpsilonEquals(rgLines_Arrayed[-1], epsilon=1.0 / (2**32)):
        rgLines_Arrayed.pop()

    return rgLines_Arrayed


def _crvWithSpansCompletelyOnFace(rgCrv, rgFace, t_Crv_Pick, fTargetTol, bDebug=False):
    """
    Only process spans of the curve whose spans start and ends are on the Face.
    """
    iSpans_OnFace = []
    for iSpan in xrange(rgCrv.SpanCount):
        t = rgCrv.SpanDomain(iSpan).T0
        pt_OnCrv = rgCrv.PointAt(t)
        bSuccess, u, v = rgFace.ClosestPoint(pt_OnCrv)
        if not bSuccess: continue
        pt_OnFace = rgFace.PointAt(u,v)
        dist = pt_OnFace.DistanceTo(pt_OnCrv)
        if dist > fTargetTol:
            if bDebug:
                sEval='dist'; print(sEval+': ',eval(sEval))
                print("PointAtStart not on underlying surface.")
            continue
        
        t = rgCrv.SpanDomain(iSpan).T1
        pt_OnCrv = rgCrv.PointAt(t)
        bSuccess, u, v = rgFace.ClosestPoint(pt_OnCrv)
        if not bSuccess: continue
        pt_OnFace = rgFace.PointAt(u,v)
        dist = pt_OnFace.DistanceTo(pt_OnCrv)
        if dist > fTargetTol:
            if bDebug:
                #sc.doc.Objects.AddPoint(pt_OnFace)
                sEval='dist'; print(sEval+': ',eval(sEval))
                print("PointAtEnd not on underlying surface.")
            continue
        
        iSpans_OnFace.append(iSpan)
    
    if not iSpans_OnFace:
        print("None of the spans of the curve are completely on the face.")
        rgFace.Brep.Dispose()
        return

    if len(iSpans_OnFace) == rgCrv.SpanCount:
        return rgCrv.Duplicate()

    if len(iSpans_OnFace) == 1:
        return rgCrv.Trim(
            rgCrv.SpanDomain(iSpans_OnFace[0]).T0,
            rgCrv.SpanDomain(iSpans_OnFace[0]).T1)

    # Create nested lists of contiguous spans so that curves that go off and
    # on the face can be correctly processed.
    iSpans_Contiguous_nests = [[iSpans_OnFace[0]]]
    for iSpan in iSpans_OnFace[1:]:
        if iSpan == iSpans_Contiguous_nests[-1][-1] + 1:
            iSpans_Contiguous_nests[-1].append(iSpan)
        else:
            iSpans_Contiguous_nests.append([iSpan])
        
    if rgCrv.IsClosed and len(iSpans_Contiguous_nests) > 1:
        if (
                iSpans_Contiguous_nests[0][0] == 0 and
                iSpans_Contiguous_nests[-1][-1] == rgCrv.SpanCount-1
        ):
            iSpans_Contiguous_nests[0] = iSpans_Contiguous_nests[-1] + iSpans_Contiguous_nests[0]
            iSpans_Contiguous_nests.pop()
    
    if len(iSpans_Contiguous_nests) == 1:
        rgC_Out = rgCrv.Trim(
                rgCrv.SpanDomain(iSpans_Contiguous_nests[0][0]).T0,
                rgCrv.SpanDomain(iSpans_Contiguous_nests[0][-1]).T1)
    elif len(iSpans_Contiguous_nests) > 1:
        for iSpan_NestIndex, iSpans_Contiguous in enumerate(iSpans_Contiguous_nests):
            for iSpan in iSpans_Contiguous:
                if rgCrv.SpanDomain(iSpan).T0 <= t_Crv_Pick <= rgCrv.SpanDomain(iSpan).T1:
                    rgC_Out = rgCrv.Trim(
                            rgCrv.SpanDomain(iSpans_Contiguous_nests[iSpan_NestIndex][0]).T0,
                            rgCrv.SpanDomain(iSpans_Contiguous_nests[iSpan_NestIndex][-1]).T1)
                else:
                    print("Curve was not picked within the face.")
                    rgCrv.Dispose()
                    return
        
    if bDebug:
        sc.doc.Objects.AddCurve(rgC_Out)

    return rgC_Out


def _getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None

    return rc[1]


def _matchCrvEndDirs(nc_ToMod, nc_Ref):
    """
    nc_ToMod is modified.
    """
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


def _rebuildCrv(rgCrv_In, tol):

    for exponent in xrange(9):
        sc.escape_test()

        spanCount = 2**exponent

        pointCount = spanCount + 3

        nc_WIP = rgCrv_In.Rebuild(
            pointCount,
            degree=3,
            preserveTangents=True)

        nc_WIP.Domain = rgCrv_In.Domain

        dev = _getDistancesBetweenCurves(rgCrv_In, nc_WIP)
        if dev > tol:
            nc_WIP.Dispose()
            continue

        return nc_WIP


def _simplifyCrvForLoft(rgCrv_In, fTargetTol, bDebug=False):
    """
    Output a degree-3 curve with only simple internal knots.
    Only output degree-3 for similar limation of _Loft and RC's Loft.
    """

    #if isinstance(rgCrv_In, rg.PolylineCurve):
    #    raise Exception("PolylineCurve is not supported.")
    if isinstance(rgCrv_In, (rg.LineCurve, rg.ArcCurve)):
        return rgCrv_In.DuplicateCurve()

    nc_WIP = rgCrv_In.ToNurbsCurve()

    if nc_WIP.Degree == 3 and nc_WIP.SpanCount == 1:
        if bDebug: print("Rail curve has only 1 span.")
        return nc_WIP

    def knotMultiplicityList(knots):
        """Returns a list."""
        i = 0
        iMulties = []
        fKnotTs_Unique = []
        while True:
            knot = knots[i]
            fKnotTs_Unique.append(knot)
            iMulti = knots.KnotMultiplicity(index=i)
            iMulties.append(iMulti)
            #print("{} at {:.4f}".format(iMulti, knot),
            i += iMulti
            if i >= knots.Count:
                break
        return iMulties


    if nc_WIP.Degree < 3:
        nc_WIP.IncreaseDegree(3)
        if nc_WIP.SpanCount == 1:
            return nc_WIP
    elif nc_WIP.Degree == 3:
        if (nc_WIP.Knots.KnotStyle in (
            rg.KnotStyle.QuasiUniform,
            rg.KnotStyle.Uniform)
        ):
            return nc_WIP

        # Non-uniform but only simple internal knots is disabled for now.
        #ms = knotMultiplicityList(nc_WIP.Knots)
        #if not(nc_WIP.IsClosed and nc_WIP.IsPeriodic):
        #    ms = ms[1:-1]
        #if all([m == 1 for m in ms]):
        #    return nc_WIP

    fTol_Rail = 0.1 * fTargetTol # Smaller tolerance in case starting curve is not exactly on surface.

    if nc_WIP.Degree == 3:

        if nc_WIP.IsPeriodic:
            nc_WIP.Knots.CreatePeriodicKnots(knotSpacing=1.0) # Modifies existing Knots.
        else:
            nc_WIP.Knots.CreateUniformKnots(knotSpacing=1.0) # Modifies existing Knots.

        nc_WIP.Domain = rgCrv_In.Domain

        dev = _getDistancesBetweenCurves(rgCrv_In, nc_WIP)
        if dev <= fTol_Rail:
            return nc_WIP
        if bDebug: print("'MakeUniform' routine result is not within {}.".format(fTol_Rail))

    nc_WIP.Dispose()

    return _rebuildCrv(rgCrv_In, fTol_Rail)


def _splitNurbsCrvsAtPolyKnots(ncs):
    ncs_Out = []
    for nc in ncs:
        if not isinstance(nc, rg.NurbsCurve):
            ncs_Out.append(nc.DuplicateCurve())
            continue

        ts_polyknots = []

        if nc.IsPeriodic:
            iKs = range(nc.Knots.Count)
        elif nc.IsClosed:
            iKs = range(nc.Knots.Count - nc.Degree)
        else:
            iKs = range(nc.Degree, nc.Knots.Count - nc.Degree)

        for iK in iKs:
            if nc.Knots.KnotMultiplicity(iK) > 1:
                ts_polyknots.append(nc.Knots[iK])

        if not ts_polyknots:
            ncs_Out.append(nc.DuplicateCurve())
            continue
        rc = nc.Split(ts_polyknots)
        if not rc:
            print("Check input.")
        else:
            ncs_Out.extend(rc)
    return ncs_Out


def _prepareCrvToFin(rgCrv_In, bSimplifyCrv, bExplodePolyCrv, bSplitAtPolyKnots, bMakeDeformable, fTargetTol, bDebug=False):
    if bExplodePolyCrv and isinstance(rgCrv_In, rg.PolyCurve):
        ncs_WIP = [_.ToNurbsCurve() for _ in rgCrv_In.Explode()]
    else:
        ncs_WIP = [rgCrv_In.ToNurbsCurve()]

    if bSimplifyCrv:
        for i in range(len(ncs_WIP)):
            rc = _simplifyCrvForLoft(ncs_WIP[i], fTargetTol, bDebug)
            if rc:
                ncs_WIP[i].Dispose()
                ncs_WIP[i] = rc
    elif bSplitAtPolyKnots:
        rc = _splitNurbsCrvsAtPolyKnots(ncs_WIP)
        for _ in ncs_WIP: _.Dispose()
        ncs_WIP = rc

    if bMakeDeformable:
        for i in range(len(ncs_WIP)):
            if ncs_WIP[i].Degree < 3 and ncs_WIP[i].Points.Count < 4:
                ncs_WIP[i].IncreaseDegree(3)

    return ncs_WIP


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

        color = sc.doc.Layers.CurrentLayer.Color

        for brep in self.breps:

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


def _createGeometryInteractively():
    """
    """

    objref_CrvToFin = _getInput_Curve()
    if objref_CrvToFin is None: return


    bUseFaceOfSelNakedEdge = Opts.values['bUseFaceOfSelNakedEdge']


    rgEdge = objref_CrvToFin.Edge()

    if rgEdge and bUseFaceOfSelNakedEdge and rgEdge.Valence == rg.EdgeAdjacency.Naked:
        idxF = objref_CrvToFin.Edge().AdjacentFaces()[0]
        rgF_In = rgEdge.Brep.Faces[idxF]
    else:
        sc.doc.Objects.UnselectAll()

        objref_Face = _getInput_Face()
        if objref_Face is None: return

        sc.doc.Objects.UnselectAll()


        rgF_In = objref_Face.Face()


    bRC_NotSPB = Opts.values['bRC_NotSPB']
    bAlignEndDirs = Opts.values['bAlignEndDirs']
    bExplodePolyCrv = Opts.values['bExplodePolyCrv']
    bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
    bSimplifyCrv = Opts.values['bSimplifyCrv']
    fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
    fAngle_End_Deg = Opts.values['fAngle_End_Deg'] if Opts.values['bVariableAngle'] else None
    bAngleChangePerCrvParam_NotLength = Opts.values['bAngleChangePerCrvParam_NotLength']
    fDistance = Opts.values['fDistance']
    fTargetTol = Opts.values['fTargetTol']
    bAddSrf = Opts.values['bAddSrf']
    iLoftType = Opts.values['iLoftType']
    bAddFinEndCrv = Opts.values['bAddFinEndCrv']
    bAddLines = Opts.values['bAddLines']
    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtKnots = Opts.values['bAtKnots']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rgC_In, t_Crv0_Pick = objref_CrvToFin.CurveParameter()


    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()


    rgC_In_TrimmedToFace = _crvWithSpansCompletelyOnFace(
        rgC_In, rgF_In, t_Crv0_Pick, fTargetTol, bDebug)
    if rgC_In_TrimmedToFace is None: return


    if rgC_In_TrimmedToFace.IsClosed and fAngle_End_Deg:
        fAngle_End_Deg = Opts.values['fAngle_End_Deg'] = sc.sticky[Opts.stickyKeys['fAngle_End_Deg']] = None
        bVariableAngle = Opts.values['bVariableAngle'] = sc.sticky[Opts.stickyKeys['bVariableAngle']] = False


    sk_conduit = 'conduit({})'.format(__file__) # StickyKey
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit



    while True:
        sc.escape_test()

        bMakeDeformable = bAlignEndDirs or (fAngle_End_Deg and fAngle_End_Deg != fAngle_Start_Deg)


        ncs_ToFin = _prepareCrvToFin(
            rgC_In_TrimmedToFace,
            bSimplifyCrv,
            bExplodePolyCrv,
            bSplitAtPolyKnots,
            bMakeDeformable,
            fTargetTol,
            bDebug)

        rgBs_FromLoft = []
        rgBs_JoinedLofts = []
        ncs_FinEnd = []
        rgLs_ForOut = []


        for nc_FinStart in ncs_ToFin:

            if bRC_NotSPB:
                # Most similar to _Fin _Direction=Normal.
                nc_FinEnd = nc_FinStart.OffsetNormalToSurface(
                    surface=rgF_In, height=fDistance)
            else:
                rc = _createArrayedLines(
                    rgCrv=nc_FinStart,
                    rgSrf=rgF_In,
                    fDistance=fDistance,
                    fAngle_Start_Deg=fAngle_Start_Deg,
                    fAngle_End_Deg=fAngle_End_Deg,
                    rgCrv_Full=rgC_In_TrimmedToFace,
                    bAngleChangePerCrvParam_NotLength=bAngleChangePerCrvParam_NotLength,
                    bAtGrevilles=True,
                    bAtKnots=False,
                    bAtEqualDivisions=False,
                    iDivisionCt=iDivisionCt)
                if not rc: continue
                rgLs_GrevillesOnly = rc

                pts_Arrayed = [line.To for line in rgLs_GrevillesOnly]
                #for pt in pts_Arrayed: sc.doc.Objects.AddPoint(pt)
                #sc.doc.Views.Redraw()

                # Create curve at other end of fin.
                print(len(nc_FinStart.GrevillePoints(False)))
                #1/0
                nc_FinEnd = nc_FinStart.Duplicate()
                nc_FinEnd.SetGrevillePoints(pts_Arrayed)

            if (
                not nc_FinEnd.IsPeriodic and
                bAlignEndDirs and
                not _matchCrvEndDirs(nc_FinEnd, nc_FinStart)
            ):
                if bEcho: print("Alignment of end tangent failed.")
                nc_WIP.Dispose()
                continue

            ncs_FinEnd.append(nc_FinEnd)

            if bAddSrf:
                rc = rg.Brep.CreateFromLoft(
                    curves=[nc_FinStart, nc_FinEnd],
                    start=rg.Point3d.Unset,
                    end=rg.Point3d.Unset,
                    loftType=rg.LoftType.Straight,
                    closed=False)
                rgBs_FromLoft.extend(rc)


            # Lines out are independent of offset method.
            if bAddLines and (bAtGrevilles or bAtKnots or bAtEqualDivisions):
                if bAtGrevilles and not(bAtKnots or bAtEqualDivisions):
                    rc = rgLs_GrevillesOnly
                else:
                    rc = _createArrayedLines(
                        rgCrv=nc_FinStart,
                        rgSrf=rgF_In,
                        fDistance=fDistance,
                        fAngle_Start_Deg=fAngle_Start_Deg,
                        fAngle_End_Deg=fAngle_End_Deg,
                        rgCrv_Full=rgC_In_TrimmedToFace,
                        bAngleChangePerCrvParam_NotLength=bAngleChangePerCrvParam_NotLength,
                        bAtGrevilles=bAtGrevilles,
                        bAtKnots=bAtKnots,
                        bAtEqualDivisions=bAtEqualDivisions,
                        iDivisionCt=iDivisionCt)
                if rc:
                    rgLs_ForOut.extend(rc)


            nc_FinStart.Dispose()



        if rgBs_JoinedLofts:
            for rgB in rgBs_JoinedLofts: rgB.Dispose()
            rgBs_JoinedLofts = None
        gBs_Out = []
        if rgBs_FromLoft:
            rgBs_JoinedLofts = rg.Brep.JoinBreps(
                rgBs_FromLoft,
                tolerance=4.0*fTargetTol)
            for _ in rgBs_FromLoft: _.Dispose()


        conduit.breps = rgBs_JoinedLofts
        conduit.crvs = ncs_FinEnd if bAddFinEndCrv else []
        conduit.lines = rgLs_ForOut

        conduit.Enabled = True

        sc.doc.Views.Redraw()

        if bEcho:
            sOut = []
            if len(rgBs_JoinedLofts) > 1: sOut.append(("{} brep(s)".format(len(rgBs_JoinedLofts))))
            if len(ncs_FinEnd) > 1: sOut.append("{} fin end curves".format(len(ncs_FinEnd)))
            if rgLs_ForOut: sOut.append("{} lines".format(len(rgLs_ForOut)))
            if sOut:
                print("Calculated {}.".format(", ".join(sOut)))


        rc = _getInput_Click()

        conduit.Enabled = False

        if rc is None:
            for _ in rgBs_JoinedLofts: _.Dispose()
            for _ in ncs_FinEnd: _.Dispose()
            return

        if not rc:
            return (
                rgBs_JoinedLofts,
                ncs_FinEnd if bAddFinEndCrv else [],
                rgLs_ForOut,
                bEcho)


        for _ in rgBs_JoinedLofts: _.Dispose()
        for _ in ncs_FinEnd: _.Dispose()



        bRC_NotSPB = Opts.values['bRC_NotSPB']
        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bExplodePolyCrv = Opts.values['bExplodePolyCrv']
        bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
        bSimplifyCrv = Opts.values['bSimplifyCrv']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        fAngle_End_Deg = Opts.values['fAngle_End_Deg'] if Opts.values['bVariableAngle'] else None
        bAngleChangePerCrvParam_NotLength = Opts.values['bAngleChangePerCrvParam_NotLength']
        fDistance = Opts.values['fDistance']
        fTargetTol = Opts.values['fTargetTol']
        bAddSrf = Opts.values['bAddSrf']
        iLoftType = Opts.values['iLoftType']
        bAddFinEndCrv = Opts.values['bAddFinEndCrv']
        bAddLines = Opts.values['bAddLines']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


def main():

    rc = _createGeometryInteractively()
    if rc is None: return

    rgBs, rgCs, rgLs, bEcho = rc

    gBs = []
    for rgB in rgBs:
        gB = sc.doc.Objects.AddBrep(rgB)
        rgB.Dispose()
        if gB != gB.Empty:
            gBs.append(gB)

    gCs = []
    for rgC in rgCs:
        gC = sc.doc.Objects.AddCurve(rgC)
        if gC != gC.Empty:
            gCs.append(gC)

    gLs = [] # LineCurves
    for line in rgLs:
        gL = sc.doc.Objects.AddLine(line)
        if gL != gL.Empty:
            gLs.append(gL)

    if bEcho:
        sOut = []
        if gBs:
            iCt_Fs = sum(sc.doc.Objects.FindId(g).BrepGeometry.Faces.Count for g in gBs)
            sOut.append(("{} brep(s) with {} faces".format(
                len(gBs), iCt_Fs)))
        if gCs: sOut.append("{} fin end curves".format(len(gCs)))
        if gLs: sOut.append("{} lines".format(len(gLs)))
        print("Added {}".format(", ".join(sOut)))


if __name__ == '__main__': main()