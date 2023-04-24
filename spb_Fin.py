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
220316: Bug fixes in option settings and Brep creation.  Refactored.
230419-20: Added use of OffsetNormalToSurface.  Refactored.
230423-24: Various major changes.  Added AlignEndDirs option.

TODO: Check accuracy of offset curve and add knots until within tolerance.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum
from System import Guid


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

    key = 'bTaperChangePerCrvParam'; keys.append(key)
    values[key] = False
    names[key] = 'TaperChangePerCrv'
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

    key = 'fTol'; keys.append(key)
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
    values[key] = 2
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

        if key == 'fTol':
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
            addOption('bTaperChangePerCrvParam')
    addOption('fDistance')
    if not Opts.values['bAlignEndDirs']:
        addOption('fTol')
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


def getInput_Curve():
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


def getInput_Face():
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
            elif Opts.riOpts['fTol'].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                Opts.riOpts['fTol'].CurrentValue = Opts.riOpts['fTol'].InitialValue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Click():
    """
    Click to toggle angle and/or direction with optional input.

    Returns:
        True, False, or None.
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
        return True

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
        return False

    if res == ri.GetResult.Number:
        if bDefaultsForDir_NotAngle:
            key = 'fDistance'
        else:
            key = 'fAngle_Start_Deg'
        Opts.riOpts[key].CurrentValue = go.Number()
        Opts.setValue(key)
        go.Dispose()
        return False

    # An option was selected.

    if Opts.values['bVariableAngle'] and go.OptionIndex() == idxs_Opt['SwapAngles']:
        Opts.riOpts['fAngle_Start_Deg'].CurrentValue, Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
                Opts.riOpts['fAngle_End_Deg'].CurrentValue, Opts.riOpts['fAngle_Start_Deg'].CurrentValue)
        Opts.setValue('fAngle_Start_Deg')
        Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return False

    if not bDefaultsForDir_NotAngle and 'FlipAngle' in idxs_Opt and go.OptionIndex() == idxs_Opt['FlipAngle']:
        #Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
        #Opts.setValue('fDistance')
        Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
        Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
        Opts.setValue('fAngle_Start_Deg')
        Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return False

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
        return False

    if 'iLoftType' in idxs_Opt:
        Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex

    for key in idxs_Opt:
        if go.Option().Index == idxs_Opt[key]:
            Opts.setValue(key, go.Option().CurrentListOptionIndex)
            break

    go.Dispose()
    return False


def _createArrayedLines(rgCurve0, rgSurface0, fDistance, fAngle_Start_Deg, bVariableAngle, fAngle_End_Deg, rgCrv0_FullForVariableTaperCalc, bTaperChangePerCrvParam, bAtGrevilles, bAtKnots, bAtEqualDivisions, iDivisionCt):
    """
    rgSurface0: rg.Surface, including rg.BrepFace
    """

    rgNurbsCrv_rgCrv0 = rgCurve0.ToNurbsCurve()

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

    if not bVariableAngle:
        angle_Rad = Rhino.RhinoMath.ToRadians(fAngle_Start_Deg)
    elif not bTaperChangePerCrvParam:
        length_Full = rgCrv0_FullForVariableTaperCalc.GetLength()

    for t in ts:
        pt_Start = rgNurbsCrv_rgCrv0.PointAt(t)
        bSuccess, u, v = rgSurface0.ClosestPoint(pt_Start)
        if not bSuccess: continue
        vect_Normal = rgSurface0.NormalAt(u, v)
        bSuccess, frame = rgNurbsCrv_rgCrv0.PerpendicularFrameAt(t)
        if not bSuccess: continue
        
        pt_NoAngle = pt_Start + vect_Normal * fDistance

        pt_End = rg.Point3d(pt_NoAngle)

        if bVariableAngle:
            if bTaperChangePerCrvParam:
                t_Normalized = rgCrv0_FullForVariableTaperCalc.Domain.NormalizedParameterAt(t)
                angle_Rad = Rhino.RhinoMath.ToRadians(
                        fAngle_Start_Deg * (1.0 - t_Normalized) +
                        fAngle_End_Deg * t_Normalized)
            else:
                length_to_t = rgCrv0_FullForVariableTaperCalc.GetLength(
                        subdomain=rg.Interval(rgCrv0_FullForVariableTaperCalc.Domain.T0, t))
                angle_Rad = Rhino.RhinoMath.ToRadians(
                        fAngle_Start_Deg * (1.0 - length_to_t/length_Full) +
                        fAngle_End_Deg * length_to_t/length_Full)

        xform_Rotation = rg.Transform.Rotation(
                angleRadians=angle_Rad,
                rotationAxis=frame.ZAxis,
                rotationCenter=frame.Origin)

        pt_End.Transform(xform_Rotation)

        rgLines_Arrayed.append(rg.Line(pt_Start, pt_End))

    if rgLines_Arrayed[0].EpsilonEquals(rgLines_Arrayed[-1], epsilon=1.0 / (2**32)):
        rgLines_Arrayed.pop()

    return rgLines_Arrayed


def _curveWithSpansCompletelyOnTheFace(rgCrv, rgFace, t_Crv_Pick, fTol, bDebug=False):
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
        if dist > fTol:
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
        if dist > fTol:
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
        rgCrv0_ToFin_Full_TrimmedToFace = rgCrv.Trim(
                rgCrv.SpanDomain(iSpans_Contiguous_nests[0][0]).T0,
                rgCrv.SpanDomain(iSpans_Contiguous_nests[0][-1]).T1)
    elif len(iSpans_Contiguous_nests) > 1:
        for iSpan_NestIndex, iSpans_Contiguous in enumerate(iSpans_Contiguous_nests):
            for iSpan in iSpans_Contiguous:
                if rgCrv.SpanDomain(iSpan).T0 <= t_Crv_Pick <= rgCrv.SpanDomain(iSpan).T1:
                    rgCrv0_ToFin_Full_TrimmedToFace = rgCrv.Trim(
                            rgCrv.SpanDomain(iSpans_Contiguous_nests[iSpan_NestIndex][0]).T0,
                            rgCrv.SpanDomain(iSpans_Contiguous_nests[iSpan_NestIndex][-1]).T1)
                else:
                    print("Curve was not picked within the face.")
                    rgCrv.Dispose()
                    return
        
    if bDebug:
        sc.doc.Objects.AddCurve(rgCrv0_ToFin_Full_TrimmedToFace)

    return rgCrv0_ToFin_Full_TrimmedToFace


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

        if not _matchCrvEndDirs(nc_WIP, rgCrv_In):
            nc_WIP.Dispose()
            continue

        dev = _getDistancesBetweenCurves(rgCrv_In, nc_WIP)
        if dev > tol:
            nc_WIP.Dispose()
            continue

        return nc_WIP


def _simplifyCrvForLoft(rgCrv_In, fTol, bDebug=False):
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

    fTol_Rail = 0.1 * fTol # Smaller tolerance in case starting curve is not exactly on surface.

    if nc_WIP.Degree == 3:

        if nc_WIP.IsPeriodic:
            nc_WIP.Knots.CreatePeriodicKnots(knotSpacing=1.0) # Modifies existing Knots.
        else:
            nc_WIP.Knots.CreateUniformKnots(knotSpacing=1.0) # Modifies existing Knots.

        dev = _getDistancesBetweenCurves(rgCrv_In, nc_WIP)
        if dev <= fTol_Rail:
            return nc_WIP
        if bDebug: print("'MakeUniform' routine result is not within {}.".format(fTol_Rail))

    nc_WIP.Dispose()

    return _rebuildCrv(rgCrv_In, fTol_Rail)


def _splitNurbsCrvsAtPolyKnots(ncs):
    ncs_Out = []
    for nc in ncs:
        ts_polyknots = []

        if nc.IsPeriodic:
            iKs = range(nc.Knots.Count)
        else:
            iKs = range(nc.Degree, nc.Knots.Count - nc.Degree)

        for iK in iKs:
            if nc.Knots.KnotMultiplicity(iK) > 1:
                ts_polyknots.append(nc.Knots[iK])

        if not ts_polyknots:
            ncs_Out.append(nc)
            continue
        rc = nc.Split(ts_polyknots)
        if not rc:
            print("Check input.")
        else:
            ncs_Out.extend(rc)
    return ncs_Out


def _prepareCrvToFin(rgCrv_In, bSimplifyCrv, bExplodePolyCrv, bSplitAtPolyKnots, fTol, bDebug=False):
    if bExplodePolyCrv and isinstance(rgCrv_In, rg.PolyCurve):
        rgCrvs_Out = rgCrv_In.Explode()
    else:
        rgCrvs_Out = [rgCrv_In.ToNurbsCurve()]

    if bSimplifyCrv:
        for i in range(len(rgCrvs_Out)):
            rc = _simplifyCrvForLoft(rgCrvs_Out[i], fTol, bDebug)
            if rc:
                rgCrvs_Out[i] = rc

    elif bSplitAtPolyKnots:
        return _splitNurbsCrvsAtPolyKnots(rgCrvs_Out)

    return rgCrvs_Out


def main():

    objref_CrvToFin = getInput_Curve()
    if objref_CrvToFin is None: return


    bAlignEndDirs = Opts.values['bAlignEndDirs']
    bUseFaceOfSelNakedEdge = Opts.values['bUseFaceOfSelNakedEdge']
    bExplodePolyCrv = Opts.values['bExplodePolyCrv']
    bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
    bSimplifyCrv = Opts.values['bSimplifyCrv']
    bRC_NotSPB = Opts.values['bRC_NotSPB']
    fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
    bVariableAngle = Opts.values['bVariableAngle']
    fAngle_End_Deg = Opts.values['fAngle_End_Deg']
    bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
    fDistance = Opts.values['fDistance']
    fTol = sc.doc.ModelAbsoluteTolerance if Opts.values['bAlignEndDirs'] else Opts.values['fTol']
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


    rgEdge = objref_CrvToFin.Edge()

    if rgEdge and bUseFaceOfSelNakedEdge and rgEdge.Valence == rg.EdgeAdjacency.Naked:
        idxF = objref_CrvToFin.Edge().AdjacentFaces()[0]
        rgF_In = rgEdge.Brep.Faces[idxF]
    else:
        sc.doc.Objects.UnselectAll()

        objref_Face = getInput_Face()
        if objref_Face is None: return


        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bExplodePolyCrv = Opts.values['bExplodePolyCrv']
        bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
        bSimplifyCrv = Opts.values['bSimplifyCrv']
        bRC_NotSPB = Opts.values['bRC_NotSPB']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        bVariableAngle = Opts.values['bVariableAngle']
        fAngle_End_Deg = Opts.values['fAngle_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        fDistance = Opts.values['fDistance']
        fTol = sc.doc.ModelAbsoluteTolerance if Opts.values['bAlignEndDirs'] else Opts.values['fTol']
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


        rgF_In = objref_Face.Face()


    rgC_In, t_Crv0_Pick = objref_CrvToFin.CurveParameter()


    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()


    rgCrv0_ToFin_Full_TrimmedToFace = _curveWithSpansCompletelyOnTheFace(
        rgC_In, rgF_In, t_Crv0_Pick, fTol, bDebug)
    if rgCrv0_ToFin_Full_TrimmedToFace is None: return


    rgC_In.Dispose()


    while True:
        sc.escape_test()

        rgCrvs_ToFin = _prepareCrvToFin(
            rgCrv0_ToFin_Full_TrimmedToFace,
            bSimplifyCrv,
            bExplodePolyCrv,
            bSplitAtPolyKnots,
            fTol,
            bDebug)


        rgCs_FinEnd = []
        rgBs_FromLoft = []
        rgLs_ForOut = []


        for rgCrv_ToFin in rgCrvs_ToFin:
            nc_OnFace = rgCrv_ToFin.ToNurbsCurve()


            if bRC_NotSPB:
                # Most similar to _Fin _Direction=Normal.
                rgNC_FinEnd = rgCrv_ToFin.OffsetNormalToSurface(
                    surface=rgF_In, height=fDistance)
            else:
                rc = _createArrayedLines(
                    rgCurve0=nc_OnFace,
                    rgSurface0=rgF_In,
                    fDistance=fDistance,
                    fAngle_Start_Deg=fAngle_Start_Deg,
                    bVariableAngle=bVariableAngle,
                    fAngle_End_Deg=fAngle_End_Deg,
                    rgCrv0_FullForVariableTaperCalc=rgCrv0_ToFin_Full_TrimmedToFace,
                    bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                    bAtGrevilles=True,
                    bAtKnots=False,
                    bAtEqualDivisions=False,
                    iDivisionCt=iDivisionCt)
                if not rc: continue
                rgLs_GrevillesOnly = rc

                pts_Arrayed = [line.To for line in rgLs_GrevillesOnly]

                # Create curve at other end of fin.
                rgNC_FinEnd = nc_OnFace.Duplicate()
                rgNC_FinEnd.SetGrevillePoints(pts_Arrayed)

            if bAlignEndDirs and not _matchCrvEndDirs(rgNC_FinEnd, rgCrv_ToFin):
                if bEcho: print("Alignment of end tangent failed.")
                nc_WIP.Dispose()
                continue

            if bAddFinEndCrv:
                rgCs_FinEnd.append(rgNC_FinEnd)

            if bAddSrf:
                rc = rg.Brep.CreateFromLoft(
                    curves=[nc_OnFace, rgNC_FinEnd],
                    start=rg.Point3d.Unset,
                    end=rg.Point3d.Unset,
                    loftType=rg.LoftType.Straight,
                    closed=False)
                rgBs_FromLoft.extend(rc)
            rgNC_FinEnd.Dispose()


            nc_OnFace.Dispose()


            # Lines out are independent of offset method.
            if bAddLines and (bAtGrevilles or bAtKnots or bAtEqualDivisions):
                if bAtGrevilles and not(bAtKnots or bAtEqualDivisions):
                    rc = rgLs_GrevillesOnly
                else:
                    rc = _createArrayedLines(
                        rgCurve0 = rgCrv_ToFin,
                        rgSurface0=rgF_In,
                        fDistance=fDistance,
                        fAngle_Start_Deg=fAngle_Start_Deg,
                        bVariableAngle=bVariableAngle,
                        fAngle_End_Deg=fAngle_End_Deg,
                        rgCrv0_FullForVariableTaperCalc=rgCrv0_ToFin_Full_TrimmedToFace,
                        bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                        bAtGrevilles=bAtGrevilles,
                        bAtKnots=bAtKnots,
                        bAtEqualDivisions=bAtEqualDivisions,
                        iDivisionCt=iDivisionCt)
                if rc:
                    rgLs_ForOut.extend(rc)


        sOut = []


        gBs_Out = []
        if rgBs_FromLoft:
            rgBs_JoinedLofts = rg.Brep.JoinBreps(
                rgBs_FromLoft,
                tolerance=4.0*fTol)
            for _ in rgBs_FromLoft: _.Dispose()
            
            if rgBs_JoinedLofts:
                gBs_Out = []
                for rgBrep1_Joined in rgBs_JoinedLofts:
                    gB = sc.doc.Objects.AddBrep(rgBrep1_Joined)
                    rgBrep1_Joined.Dispose()
                    if gB != Guid.Empty:
                        gBs_Out.append(gB)
                if bEcho:
                    sOut.append(("{} brep(s) with {} face(s)".format(
                        len(gBs_Out), len(rgBs_FromLoft))))

        gCs_FinEnd = []
        if bAddFinEndCrv:
            rgCs_FinEnd.append(rgNC_FinEnd)
            for rgC in rgCs_FinEnd:
                gCrv_FinEnd = sc.doc.Objects.AddCurve(rgNC_FinEnd)
                if gCrv_FinEnd != Guid.Empty:
                    gCs_FinEnd.append(gCrv_FinEnd)
            if bEcho:
                sOut.append("{} fin end curves".format(len(gCs_FinEnd)))

        gLs_Out = []
        if rgLs_ForOut:
            for line in rgLs_ForOut:
                gL = sc.doc.Objects.AddLine(line)
                if gL != Guid.Empty:
                    gLs_Out.append(gL)
            if bEcho:
                sOut.append("{} lines".format(len(gLs_Out)))


        sc.doc.Views.Redraw()

        if sOut:
            print("Output: {}".format(", ".join(sOut)))


        rc = getInput_Click()
        if rc is None:
            for gL in gLs_Out:
                sc.doc.Objects.Delete(gL, True)
            for gCrv_FinEnd in gCs_FinEnd:
                sc.doc.Objects.Delete(gCrv_FinEnd, True)
            for gB in gBs_Out:
                sc.doc.Objects.Delete(gB, True)
            break

        if rc: break

        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bExplodePolyCrv = Opts.values['bExplodePolyCrv']
        bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
        bSimplifyCrv = Opts.values['bSimplifyCrv']
        bRC_NotSPB = Opts.values['bRC_NotSPB']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        bVariableAngle = Opts.values['bVariableAngle']
        fAngle_End_Deg = Opts.values['fAngle_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        fDistance = Opts.values['fDistance']
        fTol = sc.doc.ModelAbsoluteTolerance if Opts.values['bAlignEndDirs'] else Opts.values['fTol']
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



        for gL in gLs_Out:
            sc.doc.Objects.Delete(gL, True)
        for gCrv_FinEnd in gCs_FinEnd:
            sc.doc.Objects.Delete(gCrv_FinEnd, True)
        for gB in gBs_Out:
            sc.doc.Objects.Delete(gB, True)


    rgF_In.Brep.Dispose()
    rgCrv0_ToFin_Full_TrimmedToFace.Dispose()


if __name__ == '__main__': main()