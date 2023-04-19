"""
This script is an alternative to _Fin.

Loose=Yes creates a second curve based on translated Greville points of the first,
and lofts between the 2 curves.

NormalDir=Yes, Loose=No, and UseSectionLines=No creates a second curve using
rg.Curve.OffsetNormalToSurface, and lofts between the 2 curves.  This is most
similar to _Loft.

There are also options to specify an angle from normal, constant or variable.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190625-27: Created, starting with another script.
...
220316: Bug fixes in option settings and Brep creation.  Refactored.
230419: Added use of OffsetNormalToSurface.  Refactored.
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

    key = 'bNormalDir'; keys.append(key)
    values[key] = True
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

    key = 'bLoose'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseSectionLines'; keys.append(key)
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
            initialValue=values[key], setLowerLimit=True, limit=1)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitPolyCrvToSegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitNurbsCrvsAtPolyKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddLines'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddFinEndCrv'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddBrep'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLoft_NotSweep'; keys.append(key)
    values[key] = False
    names[key] = 'BrepMethod'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Sweep1', 'Loft')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftType'; keys.append(key)
    listValues[key] = Enum.GetNames(rg.LoftType) # All items must be strings.
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fBrepTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
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

        #if key == 'fScale':
        #    if cls.riOpts[key].CurrentValue <= 1e-9:
        #        print("Invalid input for scale value.")
        #        cls.riOpts[key].CurrentValue = cls.values[key]

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

    addOption('bNormalDir')
    if not Opts.values['bNormalDir']:
        addOption('fAngle_Start_Deg')
        addOption('bVariableAngle')
        if Opts.values['bVariableAngle']:
            addOption('fAngle_End_Deg')
            addOption('bTaperChangePerCrvParam')
    addOption('fDistance')
    addOption('bSplitPolyCrvToSegs')
    addOption('bSplitNurbsCrvsAtPolyKnots')
    addOption('bLoose')
    if Opts.values['bLoose']:
        addOption('bAddFinEndCrv')
        addOption('bAddBrep')
    else:
        if Opts.values['bNormalDir']:
            addOption('bUseSectionLines')
        if not Opts.values['bNormalDir'] or Opts.values['bUseSectionLines']:
            addOption('bAtGrevilles')
            addOption('bAtKnots')
            addOption('bAtEqualDivisions')
            if Opts.values['bAtEqualDivisions']:
                addOption('iDivisionCt')
            addOption('bAddLines')
            addOption('bAddBrep')
            if Opts.values['bAddBrep']:
                addOption('bLoft_NotSweep')
                if Opts.values['bLoft_NotSweep']:
                    addOption('iLoftType')
                addOption('fBrepTol')
        if Opts.values['bNormalDir'] and not Opts.values['bUseSectionLines']:
            # Most similar to _Fin.
            addOption('bAddFinEndCrv')
            addOption('bAddBrep')
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

        if res == ri.GetResult.Number:
            key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        if Opts.values['bAddBrep'] and not Opts.values['bLoose'] and Opts.values['bUseSectionLines']:
            if (
                Opts.values['bLoft_NotSweep'] and
                go.OptionIndex() == idxs_Opt['iLoftType']
            ):
                Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
            elif Opts.riOpts['fBrepTol'].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                Opts.riOpts['fBrepTol'].CurrentValue = Opts.riOpts['fBrepTol'].InitialValue

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

        if res == ri.GetResult.Number:
            key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        if Opts.values['bAddBrep'] and not Opts.values['bLoose'] and Opts.values['bUseSectionLines']:
            if (
                Opts.values['bLoft_NotSweep'] and
                go.OptionIndex() == idxs_Opt['iLoftType']
            ):
                Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
            elif Opts.riOpts['fBrepTol'].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                Opts.riOpts['fBrepTol'].CurrentValue = Opts.riOpts['fBrepTol'].InitialValue

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

    if Opts.values['bNormalDir']:
        bFlippingDirInsteadOfAngle = True
    else:
        bFlippingDirInsteadOfAngle = (
            not Opts.values['bVariableAngle'] and
            (   Opts.values['fAngle_Start_Deg'] == 0.0 or
                Opts.values['fAngle_Start_Deg'] == 180.0)
            )

    if bFlippingDirInsteadOfAngle:
        go.SetCommandPrompt("Left click to flip direction")
    else:
        go.SetCommandPrompt("Left click to flip angle")

    go.AcceptNumber(True, acceptZero=True)
    go.AcceptNothing(True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    key = 'FlipDir'; idxs_Opt[key] = go.AddOption(key)
    if not Opts.values['bNormalDir'] and Opts.values['bVariableAngle']:
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
        if bFlippingDirInsteadOfAngle:
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

    if not bFlippingDirInsteadOfAngle and 'FlipAngle' in idxs_Opt and go.OptionIndex() == idxs_Opt['FlipAngle']:
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

    if Opts.values['bAddBrep'] and not Opts.values['bLoose'] and Opts.values['bUseSectionLines']:
        if (
            Opts.values['bLoft_NotSweep'] and
            go.OptionIndex() == idxs_Opt['iLoftType']
        ):
            Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
        elif Opts.riOpts['fBrepTol'].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
            Opts.riOpts['fBrepTol'].CurrentValue = Opts.riOpts['fBrepTol'].InitialValue

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


def createSweep1(fBrepTol, nc_OnFace, rgLines_Arrayed):
    """
    """

    rgLineCrvs = [rg.LineCurve(_) for _ in rgLines_Arrayed]


    rgBreps_CFS = rg.Brep.CreateFromSweep(
        rail=nc_OnFace,
        shapes=rgLineCrvs,
        closed=nc_OnFace.IsClosed,
        tolerance=fBrepTol)

    if not rgBreps_CFS:
        print("CreateFromSweep failed.")


    rgSweep1 = rg.SweepOneRail()
    #print(Rhino.RhinoMath.ToDegrees(rgSweep1.AngleToleranceRadians))
    #rgSweep1.AngleToleranceRadians
    #print(rgSweep1.GlobalShapeBlending) # Default is False.
    #rgSweep1.GlobalShapeBlending = True
    #print(rgSweep1.IsFreeform) # Default is True.
    #print(rgSweep1.IsRoadlike)
    #print(rgSweep1.IsRoadlikeFront)
    #print(rgSweep1.IsRoadlikeTop)
    #print(rgSweep1.IsRoadlineRight)
    rgSweep1.ClosedSweep = nc_OnFace.IsClosed
    rgSweep1.SweepTolerance = fBrepTol

    rgBreps_SOR = rgSweep1.PerformSweep(
        rail=nc_OnFace,
        crossSections=rgLineCrvs)

    if not rgBreps_SOR:
        print("SweepOneRail failed.")
    else:
        return rgBreps_SOR

    for _ in rgLineCrvs: _.Dispose()

    if rgBreps_CFS and rgBreps_SOR:
        iCt_pts_CFS = 0
        for brep in rgBreps_CFS:
            for ns in brep.Surfaces:
                iCt_pts_CFS += ns.Points.CountU * ns.Points.CountV

        iCt_pts_SOR = 0
        for brep in rgBreps_SOR:
            for ns in brep.Surfaces:
                iCt_pts_SOR += ns.Points.CountU * ns.Points.CountV

        #print(iCt_pts_CFS, iCt_pts_SOR)


        if iCt_pts_SOR < iCt_pts_CFS:
            return rgBreps_SOR
        return rgBreps_CFS

    elif rgBreps_CFS:
        return rgBreps_CFS

    return rgBreps_SOR


def createLoft(iLoftType, fBrepTol, nc_OnFace, rgLines_Arrayed):
    """
    """

    rgLineCrvs = [rg.LineCurve(_) for _ in rgLines_Arrayed]
    if nc_OnFace.IsPeriodic:
        closed = True
    elif nc_OnFace.IsClosed:
        if rgLineCrvs[0].Line.From.DistanceTo(rgLineCrvs[-1].Line.From) <= (1.0/2**32):
            # 2 lines are at start/end of curve and are not colinear.
            closed = False
        else:
            closed = True
    else:
        # Open curves.
        closed = False
    rgBreps1 = rg.Brep.CreateFromLoft(
        curves=rgLineCrvs,
        start=rg.Point3d.Unset,
        end=rg.Point3d.Unset,
        loftType=Enum.ToObject(rg.LoftType, iLoftType),
        closed=closed)
    for _ in rgLineCrvs: _.Dispose()

    return rgBreps1


def splitNurbsCrvsAtPolyKnots(ncs):
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


def _curveWithSpansCompletelyOnTheFace(rgCrv, rgFace, t_Crv_Pick, bDebug=False):
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
        if dist > sc.doc.ModelAbsoluteTolerance:
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
        if dist > sc.doc.ModelAbsoluteTolerance:
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


def main():
    
    objref_CrvToFin = getInput_Curve()
    if objref_CrvToFin is None: return


    bUseFaceOfSelNakedEdge = Opts.values['bUseFaceOfSelNakedEdge']
    bNormalDir = Opts.values['bNormalDir']
    fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
    bVariableAngle = Opts.values['bVariableAngle']
    fAngle_End_Deg = Opts.values['fAngle_End_Deg']
    bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
    fDistance = Opts.values['fDistance']
    bLoose = Opts.values['bLoose']
    bUseSectionLines = Opts.values['bUseSectionLines']
    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtKnots = Opts.values['bAtKnots']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
    bSplitNurbsCrvsAtPolyKnots = Opts.values['bSplitNurbsCrvsAtPolyKnots']
    bAddLines = Opts.values['bAddLines']
    bAddFinEndCrv = Opts.values['bAddFinEndCrv']
    bAddBrep = Opts.values['bAddBrep']
    bLoft_NotSweep = Opts.values['bLoft_NotSweep']
    iLoftType = Opts.values['iLoftType']
    fBrepTol = Opts.values['fBrepTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rgEdge = objref_CrvToFin.Edge()

    if (
        rgEdge is not None and
        bUseFaceOfSelNakedEdge and
        rgEdge.Valence == rg.EdgeAdjacency.Naked
    ):
        idxF = objref_CrvToFin.Edge().AdjacentFaces()[0]
        rgF_In = rgEdge.Brep.Faces[idxF]
    else:
        sc.doc.Objects.UnselectAll()

        objref_Face = getInput_Face()
        if objref_Face is None: return


        bNormalDir = Opts.values['bNormalDir']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        bVariableAngle = Opts.values['bVariableAngle']
        fAngle_End_Deg = Opts.values['fAngle_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        fDistance = Opts.values['fDistance']
        bLoose = Opts.values['bLoose']
        bUseSectionLines = Opts.values['bUseSectionLines']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
        bSplitNurbsCrvsAtPolyKnots = Opts.values['bSplitNurbsCrvsAtPolyKnots']
        bAddLines = Opts.values['bAddLines']
        bAddFinEndCrv = Opts.values['bAddFinEndCrv']
        bAddBrep = Opts.values['bAddBrep']
        bLoft_NotSweep = Opts.values['bLoft_NotSweep']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


        rgF_In = objref_Face.Face()


    #if bAddLines and not bLoose and bUseSectionLines:
    #    if bAtKnots:
    #        s  = "AtKnots"
    #        s += " only affects added arrayed lines,"
    #        s += " but AddArrayed option is disabled."
    #        print(s)
    #    if bAtEqualDivisions:
    #        s  = "AtEqualDivisions"
    #        s += " only affects added arrayed lines,"
    #        s += " but AddArrayed option is disabled."
    #        print(s)

    rgC_In, t_Crv0_Pick = objref_CrvToFin.CurveParameter()


    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()


    rgCrv0_ToFin_Full_TrimmedToFace = _curveWithSpansCompletelyOnTheFace(
        rgC_In, rgF_In, t_Crv0_Pick, bDebug)
    if rgCrv0_ToFin_Full_TrimmedToFace is None: return


    rgC_In.Dispose()


    while True:
        sc.escape_test()

        if bNormalDir:
            fAngle_Start_Deg = 0.0
        elif not bLoose:
            bUseSectionLines = True

        def createCrvsToFin():
            if bSplitPolyCrvToSegs and isinstance(rgCrv0_ToFin_Full_TrimmedToFace, rg.PolyCurve):
                rgCrvs_ToFin_Segs = rgCrv0_ToFin_Full_TrimmedToFace.Explode()
            else:
                rgCrvs_ToFin_Segs = [rgCrv0_ToFin_Full_TrimmedToFace.ToNurbsCurve()]

            if bSplitNurbsCrvsAtPolyKnots:
                return splitNurbsCrvsAtPolyKnots(rgCrvs_ToFin_Segs)

            return rgCrvs_ToFin_Segs
        rgCrvs_ToFin = createCrvsToFin()


        gLines_Arrayed = []
        gCrvs_FinEnd = []
        rgBreps1 = []


        for rgCrv_ToFin in rgCrvs_ToFin:
            nc_OnFace = rgCrv_ToFin.ToNurbsCurve()

            if bLoose:
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
                rgLines_Arrayed_GrevillesOnly = rc

                pts_Arrayed = [line.To for line in rgLines_Arrayed_GrevillesOnly]

                # Create curve at other end of fin.
                rgNC_FinEnd = nc_OnFace.Duplicate()
                rgNC_FinEnd.SetGrevillePoints(pts_Arrayed)

                if bAddFinEndCrv:
                    gCrv_FinEnd = sc.doc.Objects.AddCurve(rgNC_FinEnd)
                    if gCrv_FinEnd != Guid.Empty:
                        gCrvs_FinEnd.append(gCrv_FinEnd)

                if bAddBrep:
                    rc = rg.Brep.CreateFromLoft(
                        curves=[nc_OnFace, rgNC_FinEnd],
                        start=rg.Point3d.Unset,
                        end=rg.Point3d.Unset,
                        loftType=rg.LoftType.Straight,
                        closed=False)
                    rgBreps1.extend(rc)
                rgNC_FinEnd.Dispose()
            elif bUseSectionLines:
                if bAtGrevilles or bAtKnots or bAtEqualDivisions:
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
                    if not rc: continue
                    rgLines_Arrayed = rc

                    if bAddLines:
                        for line in rgLines_Arrayed:
                            gLine_Arrayed = sc.doc.Objects.AddLine(line)
                            if gLine_Arrayed != Guid.Empty:
                                gLines_Arrayed.append(gLine_Arrayed)

                    else:
                        gLines_Arrayed = []

                    if (
                        (bAtGrevilles and not bAtKnots and not bAtEqualDivisions) and
                        (bAddFinEndCrv or (bAddBrep and bLoft_NotSweep))
                    ):
                        nc_OnFace = rgCrv_ToFin.ToNurbsCurve()

                    if bAddBrep and rgLines_Arrayed:
                        if bLoft_NotSweep:
                            rc = createLoft(
                                iLoftType=iLoftType,
                                fBrepTol=fBrepTol,
                                nc_OnFace=nc_OnFace,
                                rgLines_Arrayed=rgLines_Arrayed)
                        else:
                            rc = createSweep1(
                                fBrepTol=fBrepTol,
                                nc_OnFace=nc_OnFace,
                                rgLines_Arrayed=rgLines_Arrayed)
                        if rc is None:
                            print("Cannot create brep(s).")
                        else:
                            rgBreps1.extend(rc)
            else:
                # Most similar to _Fin.
                rgNC_FinEnd = rgCrv_ToFin.OffsetNormalToSurface(
                    surface=rgF_In, height=fDistance)

                if bAddFinEndCrv:
                    gCrv_FinEnd = sc.doc.Objects.AddCurve(rgNC_FinEnd)
                    if gCrv_FinEnd != Guid.Empty:
                        gCrvs_FinEnd.append(gCrv_FinEnd)

                if bAddBrep:
                    rc = rg.Brep.CreateFromLoft(
                        curves=[nc_OnFace, rgNC_FinEnd],
                        start=rg.Point3d.Unset,
                        end=rg.Point3d.Unset,
                        loftType=rg.LoftType.Straight,
                        closed=False)
                    rgBreps1.extend(rc)
                rgNC_FinEnd.Dispose()


            nc_OnFace.Dispose()

        if rgBreps1:
            rgBreps1_Joined = rg.Brep.JoinBreps(
                    rgBreps1,
                    tolerance=8.0*sc.doc.ModelAbsoluteTolerance)
            for _ in rgBreps1: _.Dispose()
            
            if rgBreps1_Joined:
                gBreps1 = []
                for rgBrep1_Joined in rgBreps1_Joined:
                    gBrep1 = sc.doc.Objects.AddBrep(rgBrep1_Joined)
                    rgBrep1_Joined.Dispose()
                    if gBrep1 != Guid.Empty:
                        gBreps1.append(gBrep1)
                if bEcho:
                    print("{} brep(s) with {} face(s) created.".format(
                    len(gBreps1), len(rgBreps1)))
        else:
            gBreps1 = []
        
        sc.doc.Views.Redraw()
        

        rc = getInput_Click()
        if rc is None:
            for gLine_Arrayed in gLines_Arrayed:
                sc.doc.Objects.Delete(gLine_Arrayed, True)
            for gCrv_FinEnd in gCrvs_FinEnd:
                sc.doc.Objects.Delete(gCrv_FinEnd, True)
            for gBrep1 in gBreps1:
                sc.doc.Objects.Delete(gBrep1, True)
            break

        if rc: break

        bUseFaceOfSelNakedEdge = Opts.values['bUseFaceOfSelNakedEdge']
        bNormalDir = Opts.values['bNormalDir']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        bVariableAngle = Opts.values['bVariableAngle']
        fAngle_End_Deg = Opts.values['fAngle_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        fDistance = Opts.values['fDistance']
        bLoose = Opts.values['bLoose']
        bUseSectionLines = Opts.values['bUseSectionLines']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
        bSplitNurbsCrvsAtPolyKnots = Opts.values['bSplitNurbsCrvsAtPolyKnots']
        bAddLines = Opts.values['bAddLines']
        bAddFinEndCrv = Opts.values['bAddFinEndCrv']
        bAddBrep = Opts.values['bAddBrep']
        bLoft_NotSweep = Opts.values['bLoft_NotSweep']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']



        for gLine_Arrayed in gLines_Arrayed:
            sc.doc.Objects.Delete(gLine_Arrayed, True)
        for gCrv_FinEnd in gCrvs_FinEnd:
            sc.doc.Objects.Delete(gCrv_FinEnd, True)
        for gBrep1 in gBreps1:
            sc.doc.Objects.Delete(gBrep1, True)


    rgF_In.Brep.Dispose()
    rgCrv0_ToFin_Full_TrimmedToFace.Dispose()


if __name__ == '__main__': main()