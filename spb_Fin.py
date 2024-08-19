"""
This script is an alternative to _Fin.
Unlike _Fin, it can:
  * Modify the input curve, including simplifying it based on an input tolerance.
  * Apply any angle from the base surface, not just to be normal (0 deg.) or tangent (90 deg.).
  * Align the ends of the surfaces so that transitions to consecutives surfaces are G1.
  * Optionally output section curves and optionally the curve based on the input curve.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190625-27: Created, starting with another script.
...
230419-20: Added use of OffsetNormalToSurface.  Refactored.
230423-25: Various major changes.  Added AlignEndDirs option.  Now uses DrawConduit for preview.
230925: Changed the click behavior when angle is not variable and not 0, 90, 180, 270.
240816,18: Modified simplification of input curve routine.
        Removed use of OffsetNormalToSurface since it provides the same solution as SPB's routine.
        Added BothDirs option.

TODO:
    Review and modify or delete bSplitAtPolyKnots.
    Add support for bAlignEndDirs when the Output is SectionCrvs.
    ?: Replace curveWithSpansCompletelyOnFace with a routine that instead trims curve to the face.
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

    key = 'bLoose'; keys.append(key)
    values[key] = True
    names[key] = 'Output'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'SectionCrvs', 'LooseSrfs')
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

    key = 'fSimplifyCrvTol'; keys.append(key)
    values[key] = 0.1*sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

    key = 'bDeg3InLoftDir_Not1'; keys.append(key)
    values[key] = False
    names[key] = 'DegInLoftDir'
    riOpts[key] = ri.Custom.OptionToggle(values[key], '1', '3')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bBothDirs'; keys.append(key)
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

    key = 'bAddCrv'; keys.append(key)
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

        if key == 'fSimplifyCrvTol':
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

    addOption('bLoose')
    if Opts.values['bLoose']:
        addOption('bAlignEndDirs')
    addOption('bExplodePolyCrv')
    #addOption('bSplitAtPolyKnots')
    addOption('bSimplifyCrv')
    if Opts.values['bSimplifyCrv']:
        addOption('fSimplifyCrvTol')
    Opts.names['fAngle_Start_Deg'] = 'AngleFromNormal'
    addOption('bVariableAngle')
    if Opts.values['bVariableAngle']:
        Opts.names['fAngle_Start_Deg'] = 'StartAngle'
    addOption('fAngle_Start_Deg')
    if Opts.values['bVariableAngle']:
        addOption('fAngle_End_Deg')
    if Opts.values['bLoose']:
        addOption('bAngleChangePerCrvParam_NotLength')
    addOption('fDistance')
    if Opts.values['bLoose']:
        addOption('bDeg3InLoftDir_Not1')
        addOption('bBothDirs')
    if not Opts.values['bLoose']:
        addOption('bAddCrv')
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

        bNumberIsForDistDir = (
            not Opts.values['bLoose']
            or
            (not Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 90.0, 180.0, 270.0))
            or
            (Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 90.0, 180.0, 270.0) and
                Opts.values['fAngle_End_Deg'] == Opts.values['fAngle_Start_Deg'])
            )

        if res == ri.GetResult.Number:
            if bNumberIsForDistDir:
                key = 'fDistance'
            else:
                key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

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


        bNumberIsForDistDir = (
            not Opts.values['bLoose']
            or
            (not Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 90.0, 180.0, 270.0))
            or
            (Opts.values['bVariableAngle'] and
                Opts.values['fAngle_Start_Deg'] in (0.0, 90.0, 180.0, 270.0) and
                Opts.values['fAngle_End_Deg'] == Opts.values['fAngle_Start_Deg'])
            )


        if res == ri.GetResult.Number:
            if bNumberIsForDistDir:
                key = 'fDistance'
            else:
                key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

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

    bDefaultsForClickDirFlipOnly = (
        not Opts.values['bLoose']
        or
        (not Opts.values['bVariableAngle'] and
            Opts.values['fAngle_Start_Deg'] in (0.0, 90.0, 180.0, 270.0))
        or
        (Opts.values['bVariableAngle'] and
            Opts.values['fAngle_Start_Deg'] in (0.0, 90.0, 180.0, 270.0) and
            Opts.values['fAngle_End_Deg'] == Opts.values['fAngle_Start_Deg'])
        )

    if bDefaultsForClickDirFlipOnly:
        go.SetCommandPrompt("Left click to flip direction")
    else:
        go.SetCommandPrompt("Left click to cycle the 4 angle and direction flips")

    go.SetCommandPromptDefault("Accept result")

    go.AcceptNumber(True, acceptZero=True)
    go.AcceptNothing(True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    key = 'FlipDir'; idxs_Opt[key] = go.AddOption(key)
    if not bDefaultsForClickDirFlipOnly:
        key = 'FlipAngle'; idxs_Opt[key] = go.AddOption(key)
        if Opts.values['bVariableAngle']:
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
        if Opts.values['bDebug']:
            sEval = "bDefaultsForClickDirFlipOnly"; print("{}: {}".format(sEval, eval(sEval)))

        if bDefaultsForClickDirFlipOnly:
            Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
            Opts.setValue('fDistance')
            go.Dispose()
            return True

        if Opts.values['bDebug']:
            sEval = "Opts.values['fDistance']"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "Opts.values['fAngle_Start_Deg']"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "Opts.values['fAngle_End_Deg']"; print("{}: {}".format(sEval, eval(sEval)))

        if (
            Opts.riOpts['fDistance'].CurrentValue > 0.0 and
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue > 0.0
        ):
            #print("+D +A")
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
            Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
            Opts.setValue('fAngle_Start_Deg')
            Opts.setValue('fAngle_End_Deg')
        elif (
            Opts.riOpts['fDistance'].CurrentValue < 0.0 and
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue > 0.0
        ):
            #print("-D +A Good")
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
            Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
            Opts.setValue('fAngle_Start_Deg')
            Opts.setValue('fAngle_End_Deg')
        elif (
            Opts.riOpts['fDistance'].CurrentValue < 0.0 and
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue <= 0.0
        ):
            #print("-D -A Good")
            Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
            Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
            Opts.setValue('fDistance')
            Opts.setValue('fAngle_Start_Deg')
            Opts.setValue('fAngle_End_Deg')
        elif (
            Opts.riOpts['fDistance'].CurrentValue > 0.0 and
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue <= 0.0
        ):
            #print("+D -A")
            Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
            Opts.riOpts['fAngle_End_Deg'].CurrentValue = -Opts.riOpts['fAngle_End_Deg'].CurrentValue
            Opts.setValue('fDistance')
            Opts.setValue('fAngle_Start_Deg')
            Opts.setValue('fAngle_End_Deg')
        else:
            raise Exception("What happened?")
        go.Dispose()
        return True

    if res == ri.GetResult.Number:
        if bDefaultsForClickDirFlipOnly:
            key = 'fDistance'
        else:
            key = 'fAngle_Start_Deg'
        Opts.riOpts[key].CurrentValue = go.Number()
        Opts.setValue(key)
        go.Dispose()
        return True

    # An option was selected.

    if Opts.values['bVariableAngle'] and 'SwapAngles' in idxs_Opt and go.OptionIndex() == idxs_Opt['SwapAngles']:
        Opts.riOpts['fAngle_Start_Deg'].CurrentValue, Opts.riOpts['fAngle_End_Deg'].CurrentValue = (
                Opts.riOpts['fAngle_End_Deg'].CurrentValue, Opts.riOpts['fAngle_Start_Deg'].CurrentValue)
        Opts.setValue('fAngle_Start_Deg')
        Opts.setValue('fAngle_End_Deg')
        go.Dispose()
        return True

    if not bDefaultsForClickDirFlipOnly and 'FlipAngle' in idxs_Opt and go.OptionIndex() == idxs_Opt['FlipAngle']:
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


    #angle_Const_Rad = angle_Var_Rad = None
    if fAngle_End_Deg is None or fAngle_End_Deg == fAngle_Start_Deg:
        angles_per_param = [Rhino.RhinoMath.ToRadians(fAngle_Start_Deg)] * len(ts)
        #angle_Const_Rad = Rhino.RhinoMath.ToRadians(fAngle_Start_Deg)
    else:
        if not bAngleChangePerCrvParam_NotLength:
            length_Full = rgCrv_Full.GetLength()
        angles_per_param = []
        for t in ts:
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
            angles_per_param.append(angle_Var_Rad)

    #sEval='angles_per_param[:2]'; print(sEval+':',eval(sEval))
    #sEval='angles_per_param[-2:]'; print(sEval+':',eval(sEval))
    #sEval='angles_per_param'; print(sEval+':',eval(sEval))

    for iT, t in enumerate(ts):
        pt_Start = rgNurbsCrv_rgCrv0.PointAt(t)
        bSuccess, u, v = rgSrf_In.ClosestPoint(pt_Start)
        if not bSuccess: continue
        vect_Normal = rgSrf_In.NormalAt(u, v)
        bSuccess, frame = rgNurbsCrv_rgCrv0.PerpendicularFrameAt(t)
        if not bSuccess: continue
        
        pt_NoAngle = pt_Start + vect_Normal * fDistance

        pt_End = rg.Point3d(pt_NoAngle)

        #if angle_Const_Rad is None:
        #    if bAngleChangePerCrvParam_NotLength:
        #        t_Normalized = rgCrv_Full.Domain.NormalizedParameterAt(t)
        #        angle_Var_Rad = Rhino.RhinoMath.ToRadians(
        #            fAngle_Start_Deg * (1.0 - t_Normalized) +
        #            fAngle_End_Deg * t_Normalized)
        #    else:
        #        length_to_t = rgCrv_Full.GetLength(
        #            subdomain=rg.Interval(rgCrv_Full.Domain.T0, t))
        #        angle_Var_Rad = Rhino.RhinoMath.ToRadians(
        #            fAngle_Start_Deg * (1.0 - length_to_t/length_Full) +
        #            fAngle_End_Deg * length_to_t/length_Full)

        #sEval='angle_Const_Rad'; print(sEval+':',eval(sEval))
        #sEval='angle_Var_Rad'; print(sEval+':',eval(sEval))


        #xform_Rotation = rg.Transform.Rotation(
        #        angleRadians=angle_Const_Rad if angle_Var_Rad is None else angle_Var_Rad,
        #        rotationAxis=frame.ZAxis,
        #        rotationCenter=frame.Origin)

        xform_Rotation = rg.Transform.Rotation(
                angleRadians=angles_per_param[iT],
                rotationAxis=frame.ZAxis,
                rotationCenter=frame.Origin)



        pt_End.Transform(xform_Rotation)

        rgLines_Arrayed.append(rg.Line(pt_Start, pt_End))

    if rgLines_Arrayed[0].EpsilonEquals(rgLines_Arrayed[-1], epsilon=1.0 / (2**32)):
        rgLines_Arrayed.pop()

    return rgLines_Arrayed


def _curveWithSpansCompletelyOnFace(rgCrv, rgFace, t_Crv_Pick, fTol, bDebug=False):
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


def _tangentAngleDifference(cA, cB, bDebug=False):
    tanDiffAtStart = rg.Vector3d.VectorAngle(cA.TangentAtStart, cB.TangentAtStart)
    tanDiffAtEnd = rg.Vector3d.VectorAngle(cA.TangentAtEnd, cB.TangentAtEnd)
    if bDebug:
        sEval='Rhino.RhinoMath.ToDegrees(tanDiffAtStart)'; print(sEval+':',eval(sEval))
        sEval='Rhino.RhinoMath.ToDegrees(tanDiffAtEnd)'; print(sEval+':',eval(sEval))


def _rebuild_to_Bezier(nc_In, iDegs, fTol, bDebug=False):
    """
    More strict fTol used on degree 2 since that Bezier doesn't allow
    tangency matching of both ends when they are not already aligned.

    Returns on success: rg.NurbsCurve, float(deviation)
    Returns on fail: None
    """
    if bDebug: print("rebuild_Bezier with fTol={}:".format(fTol))

    for iDeg in iDegs:
        pointCount = iDeg + 1

        rebuilt = nc_In.Rebuild(
            pointCount=pointCount,
            degree=iDeg,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            nc_In, rebuilt, tolerance=0.1*fTol)[:2]

        if bSuccess:
            if iDeg==2:
                if fDistMax < 1e-6:
                    if bDebug:
                        print("Rebuilt within {}  Deg:{}  PtCt:{}".format(
                            fDistMax, iDeg, pointCount))
                        _tangentAngleDifference(nc_In, rebuilt, bDebug)
                    rebuilt.Domain = nc_In.Domain # Needed for variable angles.
                    return rebuilt, fDistMax
            elif fDistMax <= fTol:
                if bDebug:
                    print("Rebuilt within {}  Deg:{}  PtCt:{}".format(
                        fDistMax, iDeg, pointCount))
                    _tangentAngleDifference(nc_In, rebuilt, bDebug)
                    rebuilt.Domain = nc_In.Domain # Needed for variable angles.
                return rebuilt, fDistMax

        if bDebug:
            print("Rebuild requires {}  Deg:{}  PtCt:{}".format(
                fDistMax, iDeg, pointCount))
        rebuilt.Dispose()


def _rebuild_to_MultiSpan(nc_In, iDegs, fTol, bDebug=False):
    """
    Will not bother rebuilding to degree-2 since multiple spans are divided by
    G1-likely knots.
    """
    if bDebug: print("rebuild_MultiSpan with fTol={}:".format(fTol))

    iCt_MaxSpans = nc_In.SpanCount

    #iCt_MaxCp = int(round(nc.GetLength() / (100.0 * sc.doc.ModelAbsoluteTolerance)))
    #if bDebug: sEval='iCt_MaxCp'; print(sEval+':',eval(sEval))


    # Rebuild at maximum control point count.
    #pointCount = iCt_MaxCp

    iDegs_ForRebuildSearch = []
    rebuilts_LastSuccess = [] # per iDeg.
    fDevs = [] # per iDeg.

    for iDeg in iDegs:

        if iDeg < 3: continue

        pointCount = iDeg + iCt_MaxSpans

        rebuilt = nc_In.Rebuild(
            pointCount=pointCount,
            degree=iDeg,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            nc_In, rebuilt, tolerance=0.1*fTol)[:2]

        if bSuccess and fDistMax <= fTol:
            iDegs_ForRebuildSearch.append(iDeg)
            rebuilts_LastSuccess.append(rebuilt)
            fDevs.append(fDistMax)
            if bDebug:
                print("Rebuilt within {}  Deg:{}  PtCt:{}".format(
                    fTol, iDeg, pointCount))
        else:
            rebuilt.Dispose()


    if not iDegs_ForRebuildSearch:
        if bDebug:
            print("Not rebuilt within {} at max. span ct. of {}, so quitting rebuilding.".format(
                fTol, iCt_MaxSpans))
        return


    if bDebug: print("Binary search.")

    for i, iDeg in enumerate(iDegs_ForRebuildSearch):

        iCt_MaxCp = iDeg + iCt_MaxSpans

        iCts_Cps_Tried = [iDeg + 1, iCt_MaxCp]

        iCt_Cp_Try = (iCt_MaxCp + iDeg + 1) // 2

        while iCt_Cp_Try not in iCts_Cps_Tried:
            sc.escape_test()

            rebuilt = nc_In.Rebuild(
                pointCount=iCt_Cp_Try,
                degree=iDeg,
                preserveTangents=True)

            bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
                nc_In, rebuilt, tolerance=0.1*fTol)[:2]

            if bDebug:
                print("Degree:{}  CPtCt:{}  fDistMax:{}  WithinTol:{}".format(
                    iDeg,
                    iCt_Cp_Try,
                    fDistMax if bSuccess else bSuccess,
                    fDistMax <= fTol if bSuccess else bSuccess
                    ))


            iCts_Cps_Tried.append(iCt_Cp_Try)
            iCts_Cps_Tried.sort()

            if bSuccess and fDistMax <= fTol:
                rebuilts_LastSuccess[i].Dispose()
                rebuilts_LastSuccess[i] = rebuilt
                fDevs[i] = fDistMax
                # Bisect left.
                iCt_Cp_Try = (
                    (iCt_Cp_Try +
                        iCts_Cps_Tried[iCts_Cps_Tried.index(iCt_Cp_Try)-1]) // 2)
            else:
                rebuilt.Dispose()
                # Bisect right.
                iCt_Cp_Try = (
                    (iCt_Cp_Try +
                        iCts_Cps_Tried[iCts_Cps_Tried.index(iCt_Cp_Try)+1]) // 2)

    if len(rebuilts_LastSuccess) == 1:
        rebuilt_Winner = rebuilts_LastSuccess[0]
        fDev_Winner = fDevs[0]
    else:
        # TODO: Change the winner to that with minimum deviation instead of least CPs?
        iCts_Pts = [nc_In.Points.Count for nc_In in rebuilts_LastSuccess]
        idx_Winner = iCts_Pts.index(min(iCts_Pts))
        rebuilt_Winner = rebuilts_LastSuccess[idx_Winner]
        fDev_Winner = fDevs[idx_Winner]
    rebuilt_Winner.Domain = nc_In.Domain # Needed for variable angles.
    _tangentAngleDifference(nc_In, rebuilt_Winner, bDebug)
    return rebuilt_Winner, fDev_Winner


def _simplifyCrv(rgCrv_In, fTol, bEcho=True, bDebug=False):
    """
    Output a degree-3 curve with only simple internal knots.
    Only output degree-3 for similar limation of _Loft and RC's Loft.
    """

    #if isinstance(rgCrv_In, rg.PolylineCurve):
    #    raise Exception("PolylineCurve is not supported.")
    if isinstance(rgCrv_In, (rg.LineCurve, rg.ArcCurve)):
        return rgCrv_In.DuplicateCurve()

    nc_WIP = rgCrv_In.ToNurbsCurve()

    if nc_WIP.SpanCount == 1:
        if bDebug: print("Curve has only 1 span.")
        return nc_WIP


    rc = _rebuild_to_Bezier(
        rgCrv_In,
        iDegs=(2,3,5),
        fTol=fTol,
        bDebug=bDebug)
    if rc:
        return rc[0]

    rc = _rebuild_to_MultiSpan(
        rgCrv_In,
        iDegs=(3,5),
        fTol=fTol,
        bDebug=bDebug)
    if rc:
        return rc[0]

    #def knotMultiplicityList(knots):
    #    """Returns a list."""
    #    i = 0
    #    iMulties = []
    #    fKnotTs_Unique = []
    #    while True:
    #        knot = knots[i]
    #        fKnotTs_Unique.append(knot)
    #        iMulti = knots.KnotMultiplicity(index=i)
    #        iMulties.append(iMulti)
    #        #print("{} at {:.4f}".format(iMulti, knot),
    #        i += iMulti
    #        if i >= knots.Count:
    #            break
    #    return iMulties


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

    if nc_WIP.IsPeriodic:
        nc_WIP.Knots.CreatePeriodicKnots(knotSpacing=1.0) # Modifies existing Knots.
    else:
        nc_WIP.Knots.CreateUniformKnots(knotSpacing=1.0) # Modifies existing Knots.

    nc_WIP.Domain = rgCrv_In.Domain

    dev = _getDistancesBetweenCurves(rgCrv_In, nc_WIP)
    if dev <= fTol:
        return nc_WIP
    if bDebug: print("'MakeUniform' routine result is {} (not within {}).".format(
        dev, fTol))

    nc_WIP.Dispose()

    if bEcho:
        print("Curve could not be simplified within {} deviation.".format(fTol))


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


def updateList(list_, func, **kwargs):
    """
    Only replace items in list when the result from function is not None.
    Otherwise, copy old item into new list.
    When replaced, the old item is Disposed.

    Returns: bool indicating whether the list has been updated.
    """

    bUpdated = False

    for i, item in enumerate(list_):

        rc = func(item, **kwargs)
        if rc is not None:
            list_[i].Dispose()
            list_[i] = rc
            bUpdated = True

    return bUpdated


def _prepareCrvToFin(rgCrv_In, bSimplifyCrv, bExplodePolyCrv, bSplitAtPolyKnots, bMakeDeformable, fTol, bEcho=True, bDebug=False):

    if bExplodePolyCrv and isinstance(rgCrv_In, rg.PolyCurve):
        ncs_WIP = [_.ToNurbsCurve() for _ in rgCrv_In.Explode()]
    else:
        ncs_WIP = [rgCrv_In.ToNurbsCurve()]

    if bSimplifyCrv:
        bSimplified = updateList(ncs_WIP, _simplifyCrv, fTol=fTol, bEcho=bEcho, bDebug=bDebug)
    #elif bSplitAtPolyKnots:
    #    for i in range(len(ncs_WIP)):
    #        rc = _splitNurbsCrvAtPolyKnots(ncs_WIP[i], fTol, bDebug)
    #        if rc:
    #            ncs_WIP[i].Dispose()
    #            ncs_WIP[i] = rc
    #    rc = _splitNurbsCrvAtPolyKnots(ncs_WIP)
    #    for _ in ncs_WIP: _.Dispose()
    #    ncs_WIP = rc

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

    sk_conduit = 'conduit({})'.format(__file__) # StickyKey
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
        sc.doc.Views.Redraw()


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


    bLoose = Opts.values['bLoose']
    bAlignEndDirs = Opts.values['bAlignEndDirs']
    bExplodePolyCrv = Opts.values['bExplodePolyCrv']
    bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
    bSimplifyCrv = Opts.values['bSimplifyCrv']
    fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
    fAngle_End_Deg = Opts.values['fAngle_End_Deg'] if Opts.values['bVariableAngle'] else None
    bAngleChangePerCrvParam_NotLength = Opts.values['bAngleChangePerCrvParam_NotLength']
    fDistance = Opts.values['fDistance']
    fSimplifyCrvTol = Opts.values['fSimplifyCrvTol']
    bDeg3InLoftDir_Not1 = Opts.values['bDeg3InLoftDir_Not1']
    bBothDirs = Opts.values['bBothDirs']
    bAddCrv = Opts.values['bAddCrv']
    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtKnots = Opts.values['bAtKnots']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    #sEval='fDistance'; print(sEval+':',eval(sEval))
    #sEval='fAngle_Start_Deg'; print(sEval+':',eval(sEval))
    #sEval='fAngle_End_Deg'; print(sEval+':',eval(sEval))

    rgC_In, t_Crv0_Pick = objref_CrvToFin.CurveParameter()


    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()


    rgC_In_TrimmedToFace = _curveWithSpansCompletelyOnFace(
        rgC_In, rgF_In, t_Crv0_Pick, sc.doc.ModelAbsoluteTolerance, bDebug)
    if rgC_In_TrimmedToFace is None: return


    if rgC_In_TrimmedToFace.IsClosed and fAngle_End_Deg:
        fAngle_End_Deg = Opts.values['fAngle_End_Deg'] = sc.sticky[Opts.stickyKeys['fAngle_End_Deg']] = None
        bVariableAngle = Opts.values['bVariableAngle'] = sc.sticky[Opts.stickyKeys['bVariableAngle']] = False


    if not ((sk_conduit in sc.sticky) and sc.sticky[sk_conduit]):
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit



    while True:
        sc.escape_test()

        rgBs_FromLoft = []
        rgBs_JoinedLofts = []
        ncs_FinEnd = []
        rgLs_ForOut = []

        if not bLoose and not (bAddCrv or bAtGrevilles or bAtKnots or bAtEqualDivisions):
            print("No output has been enabled. Pick some.")
        else:
            bMakeDeformable = (
                (bLoose and bAlignEndDirs)
                or
                (fAngle_End_Deg
                 and
                 (fAngle_End_Deg != fAngle_Start_Deg))
                )

            ncs_FinStart = _prepareCrvToFin(
                rgC_In_TrimmedToFace,
                bSimplifyCrv,
                bExplodePolyCrv,
                bSplitAtPolyKnots,
                bMakeDeformable,
                fTol=fSimplifyCrvTol,
                bEcho=bEcho,
                bDebug=bDebug)

            #for nc in ncs_FinStart:
            #    sc.doc.Objects.AddCurve(nc)


            if bLoose:
                for nc_FinStart in ncs_FinStart:
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
                    #print(len(nc_FinStart.GrevillePoints(False)))
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

                    breps_Loft = rg.Brep.CreateFromLoft(
                        curves=[nc_FinStart, nc_FinEnd],
                        start=rg.Point3d.Unset,
                        end=rg.Point3d.Unset,
                        loftType=rg.LoftType.Normal if bDeg3InLoftDir_Not1 else rg.LoftType.Straight,
                        closed=False)
                    def extendOtherDir(breps):
                        for i, brep in enumerate(breps):
                            # Changing U is loft direction.
                            ns = brep.Surfaces[0]
                            interval = ns.Domain(0)
                            interval.Reverse()
                            if not ns.Extend(direction=0, interval=interval):
                                continue
                            ns = brep.Surfaces[0]
                            breps[i] = ns.ToBrep()
                    if breps_Loft and bBothDirs:
                        breps_Loft = list(breps_Loft)
                        extendOtherDir(breps_Loft)
                    rgBs_FromLoft.extend(breps_Loft)


                    nc_FinStart.Dispose()

                if rgBs_JoinedLofts:
                    for rgB in rgBs_JoinedLofts: rgB.Dispose()
                    rgBs_JoinedLofts = None
                gBs_Out = []
                if rgBs_FromLoft:
                    rgBs_JoinedLofts = rg.Brep.JoinBreps(
                        rgBs_FromLoft,
                        tolerance=2.0*sc.doc.ModelAbsoluteTolerance)
                    for _ in rgBs_FromLoft: _.Dispose()


                conduit.breps = rgBs_JoinedLofts
                conduit.crvs = []
                conduit.lines = []

            else:
                # Not Loose.
                if bAtGrevilles or bAtKnots or bAtEqualDivisions:
                    for nc_FinStart in ncs_FinStart:
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

                conduit.breps = []
                conduit.crvs = ncs_FinStart if bAddCrv else []
                conduit.lines = rgLs_ForOut

            conduit.Enabled = True

            sc.doc.Views.Redraw()

            if bEcho:
                sOut = []
                if len(rgBs_JoinedLofts) > 1: sOut.append(("{} brep(s)".format(len(rgBs_JoinedLofts))))
                if len(ncs_FinStart) > 1: sOut.append("{} fin start curves".format(len(ncs_FinEnd)))
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
                ncs_FinStart if bAddCrv else [],
                rgLs_ForOut,
                bEcho)


        for _ in rgBs_JoinedLofts: _.Dispose()
        for _ in ncs_FinEnd: _.Dispose()



        bLoose = Opts.values['bLoose']
        bAlignEndDirs = Opts.values['bAlignEndDirs']
        bExplodePolyCrv = Opts.values['bExplodePolyCrv']
        bSplitAtPolyKnots = Opts.values['bSplitAtPolyKnots']
        bSimplifyCrv = Opts.values['bSimplifyCrv']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        fAngle_End_Deg = Opts.values['fAngle_End_Deg'] if Opts.values['bVariableAngle'] else None
        bAngleChangePerCrvParam_NotLength = Opts.values['bAngleChangePerCrvParam_NotLength']
        fDistance = Opts.values['fDistance']
        fSimplifyCrvTol = Opts.values['fSimplifyCrvTol']
        bDeg3InLoftDir_Not1 = Opts.values['bDeg3InLoftDir_Not1']
        bBothDirs = Opts.values['bBothDirs']
        bAddCrv = Opts.values['bAddCrv']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        #sEval='fDistance'; print(sEval+':',eval(sEval))
        #sEval='fAngle_Start_Deg'; print(sEval+':',eval(sEval))
        #sEval='fAngle_End_Deg'; print(sEval+':',eval(sEval))


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