"""
This script is an alternative to _Project.
Segments of polycurves can be projected individually, avoiding merges through joints.
Projection of curves onto planar surfaces of polyface breps is more similar to loose projection.
For increased accuracy, all curves projected to planar faces are projected loose
to the TryGetPlane plane, then split at the intersections with the monoface brep's edges.
With AttemptRebuild enabled, curve will be projected at half the tolerance,
and a rebuild (to uniform) of the projected curve will be attempted at half tolerance.
If not successful, the projection at full tolerance will be used.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
191124-25: Created.
...
220821: Modified an option default value.
230720: Now passes fTol_Proj to a method instead of ModelAbsoluteTolerance.
231106: Improved selection routine for faces.
240712-15: Refactored and modified behavior of simplification routines.
241107, 241225: Fixed bugs created during 240712-15 refactoring.
241226-28: Fixed bugs, added another tolerance option, and refactored.

TODO:
    Reviewing tolerance values passed to:
        rg.Intersect.Intersection.CurveCurve
        rg.Curve.ProjectToBrep
        removeShortCrvsInList
        convertLinearNurbsToLinesInList
        rg.Intersect.Intersection.ProjectPointsToBreps

    Determine solution for when a loose projected curve doesn't lie within
    the projection boundary of the brep.
    Possible solutions:
        1. Keep brep-missing Greville points at their non-projected locations.
        2. Trim the curve before projection.
    
    Add support for projections of curves in overlapping directions, e.g., projecting a circle parallel to its plane.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bProjectCrvSegs'; keys.append(key)
    values[key] = True
    names[key] = "Project"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'WholeCrv', 'CrvSegs')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_Proj'; keys.append(key)
    values[key] = 0.25 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'ProjTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    key = 'bLoose'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPostProcess'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_MinLength'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    names[key] = 'MinLenTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    key = 'bOnlyLinesAndCubicNurbs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryGetArcs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAcceptRational'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryRebuildOthersUniform'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bJoinPerInputCrv'; keys.append(key)
    values[key] = True
    #names[key] = 'Join'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteInput'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iOutputLayer'; keys.append(key)
    values[key] = 1
    listValues[key] = (
        'Input',
        'Current',
        'TargetObject',
        )
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDirection'; keys.append(key)
    values[key] = 2
    listValues[key] = (
        'CPlaneX',
        'CPlaneY',
        'CPlaneZ',
        'WorldX',
        'WorldY',
        'WorldZ',
        'View',
        'Custom',
        )
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'vectCustom'; keys.append(key)
    values[key] = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


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

        if key == 'fTol_Proj':
            # Rhino.RhinoMath.ZeroTolerance == 2**-32 in RC >= V7. It is 1.0e-12 in previous versions.
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            if cls.riOpts[key].CurrentValue <= 2**-32:
                cls.values[key] = cls.riOpts[key].CurrentValue = 2**-32
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return


        if key == 'fTol_MinLength':
            # Rhino.RhinoMath.ZeroTolerance == 2**-32 in RC >= V7. It is 1.0e-12 in previous versions.
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            if cls.riOpts[key].CurrentValue <= 2**-32:
                cls.values[key] = cls.riOpts[key].CurrentValue = 2**-32
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(rdObjs_toHighlight, sPrompt, rdGeomFilter):
    """
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt(prompt=sPrompt)

    go.GeometryFilter = rdGeomFilter

    #go.AlreadySelectedObjectSelect = True # Default is False
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    bPreselectedObjsChecked = False
    go.EnablePreSelect(True, ignoreUnacceptablePreselectedObjects=True)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bProjectCrvSegs')
        addOption('fTol_Proj')
        addOption('bLoose')
        if not Opts.values['bLoose']:
            addOption('bPostProcess')
            if Opts.values['bPostProcess']:
                addOption('fTol_MinLength')
                addOption('bOnlyLinesAndCubicNurbs')
                if not Opts.values['bOnlyLinesAndCubicNurbs']:
                    addOption('bTryGetArcs')
                    addOption('bAcceptRational')
                    addOption('bTryRebuildOthersUniform')
        addOption('bJoinPerInputCrv')
        addOption('bDeleteInput')
        addOption('iOutputLayer')
        addOption('iDirection')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked:
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            bPreselectedObjsChecked = True
            if go.ObjectsWerePreselected:
                continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected or a number was entered.
        if res == ri.GetResult.Number:
            key = 'fTol_Proj'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['iOutputLayer']:
            Opts.values['iOutputLayer'] = (
                    go.Option().CurrentListOptionIndex)
        elif go.Option().Index == idxs_Opt['iDirection']:
            Opts.values['iDirection'] = (
                    go.Option().CurrentListOptionIndex)

            if Opts.listValues['iDirection'][go.Option().CurrentListOptionIndex] == 'Custom':
                rc = rs.GetLine(
                    mode=1, point=None,
                    message1='Projection direction',
                    message3='Second direction point',
                    )
                if not rc:
                    Opts.values['iDirection'] = 0
                else:
                    Opts.values['vectCustom'] = rg.Vector3d(rc[1] - rc[0])
                    Opts.values['vectCustom'].Unitize()

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.Brep):
        return rhObj
    elif isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        geom = rhObj.Geometry()
        if rhObj.GeometryComponentIndex.ComponentIndexType == rg.ComponentIndexType.BrepFace:
            geom = geom.DuplicateFace(duplicateMeshes=False)
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    if isinstance(geom, rg.Brep):
        return geom


def formatDistance(fDistance):
    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def projectCurve_Loose(rgCrv_In, rgB, vectDir, fTol_Proj=None, bDebug=False):
    """
    Parameters:
        rgC0
        rgB
        vectDir
        fTol_Proj: float = Tolerance of projected control points, not the maximum
                           deviation of curves to surfaces.
        bDebug: bool
    """

    #TODOif fTol_Proj is None: fTol_Proj = 1e-6#0.1 * sc.doc.ModelAbsoluteTolerance


    rgCs_Out = []

    # Duplicate curve and translate Greville points to projected locations.
    nc_toProj = rgCrv_In.ToNurbsCurve()

    grPts1 = []
    for gr in nc_toProj.GrevillePoints(all=False):
        pts_Proj = rg.Intersect.Intersection.ProjectPointsToBreps(
            breps=[rgB],
            points=[gr],
            direction=vectDir,
            tolerance=fTol_Proj)
        if pts_Proj:
            grPts1.append(pts_Proj[0])
        else:
            grPts1.append(gr)
            
        #sc.doc.Objects.AddPoint(pts_Proj[0])
    #sc.doc.Views.Redraw(); 1/0

    nc_toProj.SetGrevillePoints(grPts1)

    rgCs_Out.append(nc_toProj)


    return rgCs_Out

    # TODO: Move the following elsewhere.



    # This will clean up any short segments that had passed the previous
    # checks, but have been deformed to short segments after JoinCurves.
    # _Project doesn't address this problem.
    removeShortSegmentsInEachCrv_inList(rgCs_Proj_Joined_ThisC0, fTol_MinLength)

    return rgCs_Proj_Joined_ThisC0


def removeShortCrvsInList(rgCrvs, tolerance):
    bFound = False
    i = len(rgCrvs) - 1
    while i >= 0:
        sc.escape_test()

        rgC = rgCrvs[i]

        length = rgC.GetLength()

        if length < tolerance:
            bFound = True
            del rgCrvs[i]
            rgC.Dispose()
            #print("{}-long curve created at {} tolerance ignored.".format(
            #    formatDistance(length),
            #    formatDistance(tolerance)))

        i -= 1

    return bFound


def explodePolyCrv(pc):
    pc.RemoveNesting()
    rv = pc.Explode()
    if not rv:
        raise Exception("{} resulted from rg.PolyCurve.Explode.".format(rv))
    return rv


def explodePolyCrvsInList(rgCrvs):
    bFound = False
    i = len(rgCrvs) - 1
    while i >= 0:
        sc.escape_test()

        rgC = rgCrvs[i]

        if isinstance(rgC, rg.PolyCurve):
            rv = explodePolyCrv(rgC)
            if not rv:
                1/0
            else:
                bFound = True
                rgCrvs[i:i+1] = rv
                rgC.Dispose()

        i -= 1

    return bFound


def convertCrvsForPeriodicInList(rgCrvs):
    bFound = False
    for i, rgC in enumerate(rgCrvs):
        if rgC.IsArc():
            print(rgC.GetType().Name)
            1/0
            bFound = True
        elif isinstance(rgC, rg.NurbsCurve):
            continue
        else:
            print(rgC.GetType().Name)
            1/0
            rgCrvs[i] = rgC.ToNurbsCurve()
            rgC.Dispose()
            bFound = True

    return bFound


def convertLinearNurbsToLinesInList(rgCrvs, tolerance=1e-6):
    bFound = False
    for i, c in enumerate(rgCrvs):
        if isinstance(c, rg.NurbsCurve):
            if not c.IsClosed:
                if c.IsLinear(tolerance):
                    rgCrvs[i] = rg.LineCurve(
                        c.PointAtStart,
                        c.PointAtEnd)
                    #sEval = "c.IsDocumentControlled"; print(sEval,'=',eval(sEval))
                    c.Dispose()
                    bFound = True
    return bFound


def removeShortSegmentsInEachCrv_inList(rgCrvs, tolerance):
    # rgCrvs: list or array
    bFound = False
    for i, c in enumerate(rgCrvs):
        if c.RemoveShortSegments(tolerance):
            bFound = True
    return bFound


def projectCurve_Greville_pts_to_plane(rgC_In, rgPlane, vectDir):
    rgC_Working = rgC_In.ToNurbsCurve()

    grPts1 = []

    for gr0 in rgC_Working.GrevillePoints(all=False):

        rgLine = rg.Line(gr0, vectDir)
        bSuccess, tLine = Rhino.Geometry.Intersect.Intersection.LinePlane(
            rgLine, rgPlane);
        if not bSuccess:
            s  = "Point missed the plane when projecting loose to a plane."
            s += "  Check results for accuracy."
            print(s)
            sEval = "rgC_Working.IsDocumentControlled"; print(sEval,'=',eval(sEval))
            rgC_Working.Dispose()
            return

        gr1 = rgLine.PointAt(tLine);

        grPts1.append(gr1)

    rgC_Working.SetGrevillePoints(grPts1)
    
    return rgC_Working


def trimCurves_to_brep_edges(rgC_fromProjectLoose, rgB, fTol_Proj):
    """
    Returns list(rg.Curve)
    """

    ts_atSplits = []

    for rgEdge in rgB.Edges:
        crvinters = rg.Intersect.Intersection.CurveCurve(
            rgC_fromProjectLoose,
            rgEdge,
            tolerance=fTol_Proj,
            overlapTolerance=0.0)
        if crvinters.Count == 0: continue # to next curve.
        for crvinter in crvinters:
            ts_atSplits.append(crvinter.ParameterA)

    if not ts_atSplits: 
        # Presuming projected curve lies completely on the face.
        return [rgC_fromProjectLoose]

    segs2_PostSplit = rgC_fromProjectLoose.Split(ts_atSplits)

    # Determine which of the split project-to-plane curves should be kept.

    rgCs_Out = []

    for rgC in segs2_PostSplit:
        midDomainPt = rgC.PointAt(rgC.Domain.Mid)
        closestPt = rgB.ClosestPoint(midDomainPt)

        dist = midDomainPt.DistanceTo(closestPt)
        #print(dist)

        if dist <= 1e-6:
            rgCs_Out.append(rgC)

    return rgCs_Out



    # Old method. Need curves projected to brep of face.

    pts_toDetermineTrims = []
    for rgC in rgCs_Proj_This_face_or_plane:
        pts_toDetermineTrims.append(rgC.PointAt(rgC.Domain.Mid))
        rgC.Dispose()

    for rgC in segs2_PostSplit:
        for pt in pts_toDetermineTrims:
            bSuccess, t = rgC.ClosestPoint(pt)
            if not bSuccess: continue # To next point.
            if rgC.PointAt(t).DistanceTo(pt) <= 0.1 * tolerance:
                # Tight tolerance above is to avoid grabbing segments
                # adjacent to correct results that also happen to be short.
                rgCs_Out.append(rgC)
                break # To next curve from split.


    #for rgC in rgCs_Proj_ThisSeg:
    #    sc.doc.Objects.AddCurve(rgC)
    #sc.doc.Views.Redraw(); 1/0

    removeShortSegmentsInEachCrv_inList(rgCs_Out, tolerance=fTol_MinLength)

    return rgCs_Out


def projectCurve_to_plane_and_trim_to_face(rgC_In, rgPlane, rgB, vectDir, fTol_Proj):
    """
    Returns list(rg.Curve)
    """

    rgCs_Out = []

    rgC_fromProjectLoose = projectCurve_Greville_pts_to_plane(rgC_In, rgPlane, vectDir)
    if rgC_fromProjectLoose is None:
        return

    return(trimCurves_to_brep_edges(
        rgC_fromProjectLoose,
        rgB=rgB,
        fTol_Proj=fTol_Proj,
        ))


def cleanProjectedCrvs_inList(rgCs_toMod, fTol_MinLength, bDebug=False):
    """
    rgCs_Out: list, not array
    Returns list(list(rg.Curve))
    """

    #rgCs_Out = list(rgCs_Out) # Need to convert from Array to allow list slicing.

    bModified = False

    lenList = len(rgCs_toMod)

    if removeShortCrvsInList(rgCs_toMod, tolerance=1e-6):
        bModified = True
        lenList = len(rgCs_toMod)
    if removeShortSegmentsInEachCrv_inList(rgCs_toMod, tolerance=1e-6):
        bModified = True
        if bDebug: print("Short segments removed before exploding polycurves.")
    if explodePolyCrvsInList(rgCs_toMod):
        bModified = True
        lenList = len(rgCs_toMod)
    if removeShortCrvsInList(rgCs_toMod, tolerance=fTol_MinLength):
        bModified = True
        lenList = len(rgCs_toMod)
    if convertLinearNurbsToLinesInList(rgCs_toMod, 1e-6):
        bModified = True
        lenList = len(rgCs_toMod)
    if removeShortSegmentsInEachCrv_inList(rgCs_toMod, tolerance=fTol_MinLength):
        bModified = True
        if bDebug: print("Short segments removed at end of routine.")

    return bModified


def shouldCrvBeProjectedToPlaneInsteadOfFace(crv):
    """
    It's been found that while LineCurves, PolylineCurves, and 2-point NurbsCurves
    are projected accurately to planar faces,
    ArcCurves and NurbCurves with more than 2 points are often not.
    """

    #print(crv.GetType().Name)

    if isinstance(crv, rg.ArcCurve):
        return True
    if isinstance(crv, rg.NurbsCurve):
        if crv.Points.Count > 2:
            return True

    return False


def projectCrv_to_1faceBreps_and_planes(rgCrv_In, rgBs_1Face, rgPlanes, vectDir, fTol_Proj):
    """
    Returns list(list(rg.Curve))
    """

    rgCs_Out = []

    rgPlanes_ = rgPlanes # For debugging.

    for rgB_1F, rgPlane in zip(rgBs_1Face, rgPlanes_):

        rgCs_Proj_to_1FB = rg.Curve.ProjectToBrep(
            curve=rgCrv_In,
            brep=rgB_1F,
            direction=vectDir,
            tolerance=fTol_Proj)

        if not rgCs_Proj_to_1FB:
            continue

        if not rgPlane:
            rgCs_Out.extend(rgCs_Proj_to_1FB)
            #sc.doc.Objects.AddBrep(rgB_1F)
            continue

        # Since rgCs_Proj_to_1FB exists, project to the current rgPlane.

        if shouldCrvBeProjectedToPlaneInsteadOfFace(rgCrv_In):
            rgCs_Proj_to_plane = projectCurve_to_plane_and_trim_to_face(
                rgCrv_In,
                rgPlane,
                rgB_1F,
                vectDir,
                fTol_Proj=fTol_Proj)
            if rgCs_Proj_to_plane:
                rgCs_Out.extend(rgCs_Proj_to_plane)
                for _ in rgCs_Proj_to_1FB: _.Dispose()
                continue

        rgCs_Out.extend(rgCs_Proj_to_1FB)

    return rgCs_Out


def isUniformNurbsCurve(nc):
    if not isinstance(nc, rg.NurbsCurve): return


    # Bezier?
    if nc.Points.Count == nc.Degree + 1:
        return True


    # Any internal polyknots?
    start = 0 if nc.IsPeriodic else nc.Degree
    end = nc.Knots.Count - (0 if nc.IsPeriodic else nc.Degree) - 1

    for i in range(start, end):
        if nc.Knots.KnotMultiplicity(i) > 1:
            return False


    # Any non-uniform knot spans?
    start = 0 if nc.IsPeriodic else nc.Degree - 1
    end = nc.Knots.Count - (0 if nc.IsPeriodic else nc.Degree - 1) - 1

    span0 = nc.Knots[start+1] - nc.Knots[start]

    def areEpsilonEqual(a, b, epsilon):
        # This is a relative comparison.
        delta = abs(a - b)
        fRelComp = delta / max(abs(a), abs(b))
        return fRelComp < epsilon


    for i in range(start+1, end):
        if nc.Knots.KnotMultiplicity(i) > 1:
            return False
        if not areEpsilonEqual(
            span0, nc.Knots[i+1] - nc.Knots[i],
            epsilon=2**-32):
                return False

    return True


def are_all_crvs_lines_or_uniformNonrationalCubic(rgCrvs):
    for c in rgCrvs:
        if isinstance(c, rg.LineCurve):
            continue
        if not isinstance(c, (rg.NurbsCurve)):
            # Is Arc, Polyline?, PolyCurve?
            return False
        if c.IsRational:
            return False
        if not c.Degree == 3:
            return False
        if not isUniformNurbsCurve(c):
            return False
    return True


def rebuildUniformNonrationalCubic(rgCrv, tolerance, bDebug=False):

    # Try rebuilding as a Bezier.
    rebuilt = rgCrv.Rebuild(
        pointCount=4,
        degree=3,
        preserveTangents=True)

    bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
        rgCrv, rebuilt, tolerance=0.1*tolerance)[:2]

    if bSuccess and fDistMax <= tolerance:
        return rebuilt

    rebuilt.Dispose()


    iCt_MaxCp = round(rgCrv.GetLength() / (100.0 * sc.doc.ModelAbsoluteTolerance))
    if bDebug: sEval='iCt_MaxCp';  print("{}: {}".format(sEval, eval(sEval)))


    # Try rebuilding at maximum allowed control point count.

    rebuilt = rgCrv.Rebuild(
        pointCount=iCt_MaxCp,
        degree=3,
        preserveTangents=True)

    bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
        rgCrv, rebuilt, tolerance=0.1*tolerance)[:2]

    if not bSuccess or fDistMax > tolerance:
        # Fail at maximum control point count, so quit searching.
        rebuilt.Dispose()
        return


    # Binary search to find least control point count.
    rebuilt_LastSuccess = rebuilt
    iCts_Cps_Tried = [4, iCt_MaxCp]

    iCt_Cp_Try = (iCt_MaxCp + 4) // 2

    if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

    while iCt_Cp_Try not in iCts_Cps_Tried:
        sc.escape_test()


        rebuilt = rgCrv.Rebuild(
            pointCount=iCt_Cp_Try,
            degree=3,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            rgCrv, rebuilt, tolerance=0.1*tolerance)[:2]


        iCts_Cps_Tried.append(iCt_Cp_Try)
        iCts_Cps_Tried.sort()

        if bSuccess and fDistMax <= tolerance:
            rebuilt_LastSuccess.Dispose()
            rebuilt_LastSuccess = rebuilt
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

        if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

    return rebuilt_LastSuccess


def reportCrvTypes(rgCrvs):
    sTypes = [c.GetType().Name for c in rgCrvs]
    for sType in set(sTypes):
        print("[{}] {}".format(sTypes.count(sType), sType))


def tryConvertNurbsToArcs_inList(rgCrvs_toMod, tol):
    """
    rgCrvs_toMod is modified.
    Returns: None
    """
    fTol_MaxRad = 1e3
        
    for i, c in enumerate(rgCrvs_toMod):
        if isinstance(c, rg.NurbsCurve):
            if c.IsClosed:
                b, circle = c.TryGetCircle(tol)
                if b:
                    if circle.Radius <= fTol_MaxRad:
                        rgCrvs_toMod[i] = rg.ArcCurve(circle)
                        c.Dispose()
            else:
                b, arc = c.TryGetArc(tol)
                if b:
                    if arc.Radius <= fTol_MaxRad:
                        rgCrvs_toMod[i] = rg.ArcCurve(arc)
                        c.Dispose()


def rebuild(nc, tolerance, bDebug=False):

    # Try rebuilding as a Bezier.
    for iDeg in 2, 3, 5:
        rebuilt = nc.Rebuild(
            pointCount=iDeg + 1,
            degree=iDeg,
            preserveTangents=True)

        bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
            nc, rebuilt, tolerance=0.1*tolerance)[:2]

        if bSuccess and fDistMax <= tolerance:
            return rebuilt

        rebuilt.Dispose()


    iCt_MaxCp = round(nc.GetLength() / (100.0 * sc.doc.ModelAbsoluteTolerance))


    # Try rebuilding at Degree 5 and maximum allowed control point count.
    rebuilt = nc.Rebuild(
        pointCount=iCt_MaxCp,
        degree=5,
        preserveTangents=True)

    bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
        nc, rebuilt, tolerance=0.1*tolerance)[:2]

    if not bSuccess or fDistMax > tolerance:
        # Fail at maximum control point count, so quit searching.
        rebuilt.Dispose()
        return


    # Binary search.
    rebuilt_LastSuccess = rebuilt
    iCts_Cps_Tried = [4, iCt_MaxCp]

    iCt_Cp_Try = (iCt_MaxCp + 4) // 2

    if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

    while iCt_Cp_Try not in iCts_Cps_Tried:
        sc.escape_test()


        for iDeg in 3, 5:

            if iCt_Cp_Try < iDeg+1:
                pass

            rebuilt = nc.Rebuild(
                pointCount=iCt_Cp_Try,
                degree=iDeg,
                preserveTangents=True)

            bSuccess, fDistMax = rg.Curve.GetDistancesBetweenCurves(
                nc, rebuilt, tolerance=0.1*tolerance)[:2]

            if bSuccess and fDistMax <= tolerance:
                break


        iCts_Cps_Tried.append(iCt_Cp_Try)
        iCts_Cps_Tried.sort()

        if bSuccess and fDistMax <= tolerance:
            rebuilt_LastSuccess.Dispose()
            rebuilt_LastSuccess = rebuilt
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

        if bDebug: sEval='iCt_Cp_Try';  print("{}: {}".format(sEval, eval(sEval)))

    return rebuilt_LastSuccess


def convertPolysToNurbs(rgCrvs_toMod):
    for i, c in enumerate(rgCrvs_toMod):
        if isinstance(c, rg.PolyCurve):
            rgCrvs_toMod[i] = c.ToNurbsCurve()
            c.Dispose()


def projectCrv_Try_to_output_only_lines_and_cubicNurbs(fTol_Proj_Total, projectCrvToGeom_wrapped, bDebug=False):

    projs_FullTol = projectCrvToGeom_wrapped(tolerance=fTol_Proj_Total)
    if not projs_FullTol:
        return

    if bDebug:
        print('-'*20)
        reportCrvTypes(projs_FullTol)

    if are_all_crvs_lines_or_uniformNonrationalCubic(projs_FullTol):
        return projs_FullTol

    # Project to a different tolerance and rebuild.

    for fTol_Proj_WIP in (0.5*fTol_Proj_Total, 0.1*fTol_Proj_Total):

        projs_atTrialTol = projectCrvToGeom_wrapped(tolerance=fTol_Proj_WIP)

        if are_all_crvs_lines_or_uniformNonrationalCubic(projs_atTrialTol):
            for proj in projs_FullTol: proj.Dispose()
            return projs_atTrialTol

        for i, proj in enumerate(projs_atTrialTol):
            if are_all_crvs_lines_or_uniformNonrationalCubic([proj]):
                continue
            rebuilt = rebuildUniformNonrationalCubic(
                proj,
                tolerance=fTol_Proj_Total-fTol_Proj_WIP, 
                bDebug=bDebug)

            if rebuilt:
                projs_atTrialTol[i] = rebuilt
                proj.Dispose()
                continue
            else:
                for c in projs_atTrialTol: c.Dispose()
                break # out of loop of projs_atTrialTol to next tolerance for rebuild.

        else:
            for proj in projs_FullTol: proj.Dispose()
            return projs_atTrialTol
    else:
        if bDebug:
            print("All rebuild tolerances failed for segment.")
        return projs_FullTol


def projectCrv_Try_to_output_per_arguments(fTol_Proj, projectCrvToGeom_wrapped, fTol_MinLength, bTryGetArcs, bAcceptRational, bTryRebuildOthersUniform, bDebug=False):
    """
    If necessary, reproject through various tolerances.
    """

    projs_FullTol = projectCrvToGeom_wrapped(tolerance=fTol_Proj)
    if not projs_FullTol:
        return

    if bDebug:
        print('-'*20)
        reportCrvTypes(projs_FullTol)


    if bTryGetArcs:
        tryConvertNurbsToArcs_inList(projs_FullTol, 1e-6)

    if cleanProjectedCrvs_inList(projs_FullTol, fTol_MinLength, bDebug=bDebug):
        if bDebug: print("Curve was cleaned.")

    if bTryGetArcs:
        tryConvertNurbsToArcs_inList(projs_FullTol, 1e-6)


    def areAllCurvesSimplified(rgCrvs, bTryRebuildOthersUniform):
        for c in rgCrvs:
            if isinstance(c, rg.LineCurve):
                pass
            elif isinstance(c, rg.ArcCurve):
                if not bTryGetArcs:
                    return False
            elif isinstance(c, rg.NurbsCurve):
                if c.IsRational:
                    if not bAcceptRational:
                        return False
                if bTryRebuildOthersUniform and not isUniformNurbsCurve(c):
                    return False
            else:
                raise ValueError("{} in areAllCurvesSimplified.".format(c.GetType().Name))
    
        return True


    if areAllCurvesSimplified(projs_FullTol, bTryRebuildOthersUniform):
        return projs_FullTol

    # Project to a different tolerance and rebuild.

    for fTol_Proj_WIP in (0.5*fTol_Proj, 0.1*fTol_Proj):

        projs_atTrialTol = projectCrvToGeom_wrapped(tolerance=fTol_Proj_WIP)

        if bTryGetArcs:
            tryConvertNurbsToArcs_inList(projs_atTrialTol, fTol_Proj-1e-6)

        if cleanProjectedCrvs_inList(projs_atTrialTol, fTol_MinLength, bDebug=bDebug):
            if bDebug: print("Curve was cleaned.")

        if bTryGetArcs:
            tryConvertNurbsToArcs_inList(projs_atTrialTol, 1e-6)

        if areAllCurvesSimplified(projs_atTrialTol, bTryRebuildOthersUniform):
            for proj in projs_FullTol: proj.Dispose()
            return projs_atTrialTol

        potentials = []
        for proj in projs_atTrialTol:
            if isinstance(proj, (rg.LineCurve, rg.ArcCurve)):
                potentials.append(proj)
                continue
            
            if not isinstance(proj, rg.NurbsCurve):
                raise ValueError("{}, not a NurbsCurve!".format(proj.GetType().Name))

            if proj.IsRational and rg.Curve.IsEllipse(proj, tolerance=fTol_Proj-fTol_Proj_WIP):
                potentials.append(proj)
            elif not bTryRebuildOthersUniform and proj.IsRational and isUniformNurbsCurve(proj):
                potentials.append(proj)
            elif isUniformNurbsCurve(proj):
                potentials.append(proj)
            else:
                # Rebuild to complement of projection tolerance.
                rebuilt = rebuild(
                    proj,
                    tolerance=fTol_Proj-fTol_Proj_WIP,
                    bDebug=bDebug)
                if rebuilt:
                    potentials.append(rebuilt)
                    continue
                else:
                    for rgC in potentials: rgC.Dispose()
                    break # out of for loop to next tolerance for rebuild.

        else:
            # Successful simplification.
            for proj in projs_FullTol: proj.Dispose()
            return potentials
    else:
        if bDebug: print("All rebuild tolerances failed for segment.")
        return projs_FullTol


def getJoints(rgCs_In, joinTolerance):
    """
    """

    joined = rg.Curve.JoinCurves(rgCs_In, joinTolerance=joinTolerance)

    joints = []
    for rgC in joined:
        if not isinstance(rgC, rg.PolyCurve): continue

        if rgC.IsClosed:
            joints.append(rgC.PointAtStart)

        for iSeg in range(1, rgC.SegmentCount):
            seg = rgC.SegmentCurve(iSeg)
            joints.append(seg.PointAtStart)
            seg.Dispose()

        rgC.Dispose()

    return joints


def joinCurves(rgCs_In, tolerance, bDebug=False):
    """
    Since JoinCurves, as well as _Join, can skew the knot vector,
    join the curves by only moving control points as needed.
    """

    joints = getJoints(rgCs_In, tolerance)
    if bDebug:
        for _ in joints: sc.doc.Objects.AddPoint(_)
    sc.doc.Views.Redraw()

    rgCrvs_toJoin = [rgC.DuplicateCurve() for rgC in rgCs_In]

    iJ = 0

    # Before JoinCurves, adjust ends of NurbsCurves.

    for joint in joints:
        sc.escape_test()

        iCt_EndsProcessed = 0

        for rgCrv in rgCrvs_toJoin:
                
            # Only NurbsCurves should have their control points translated when needed.
            if not isinstance(rgCrv, rg.NurbsCurve):
                #print("Curve is a {}.".format(rgCrv.GetType().Name)
                #sc.doc.Objects.AddCurve(rgCrv)
                #sc.doc.Views.Redraw(); 1/0
                continue # to next curve to include in JoinCurves.

            if rgCrv.PointAtStart.DistanceTo(joint) < tolerance:
                rgCrv.Points.SetPoint(
                    index=0,
                    point=joint,
                    weight=rgCrv.Points[0].Weight)
                #sc.doc.Objects.AddCurve(rgCrv); sc.doc.Views.Redraw(); #1/0
                iCt_EndsProcessed += 1
            elif rgCrv.PointAtEnd.DistanceTo(joint) < tolerance:
                rgCrv.Points.SetPoint(
                    index=rgCrv.Points.Count-1,
                    point=joint,
                    weight=rgCrv.Points[0].Weight)
                #sc.doc.Objects.AddCurve(rgCrv); sc.doc.Views.Redraw(); #1/0
                iCt_EndsProcessed += 1

            if iCt_EndsProcessed == 2:
                break # for loop to next joint.

    return list(rg.Curve.JoinCurves(rgCrvs_toJoin, joinTolerance=1e-6))


def projectCurve_NotLoose(rgC_In, rgBs_1Face, rgPlanes, vectDir, **kwargs):
    """
    Parameters:
        fTol_Proj
        bPostProcess
        fTol_MinLength
        bOnlyLinesAndCubicNurbs
        bTryGetArcs
        bAcceptRational
        bTryRebuildOthersUniform


    Projecting to breps of individual faces results in simpler curves for any
    planar (not only PlaneSurface) faces.
    For more accurate results, the curves's control points should be projected to
    the Plane themselves.
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTol_Proj = getOpt('fTol_Proj')
    bPostProcess = getOpt('bPostProcess')
    fTol_MinLength = getOpt('fTol_MinLength')
    bOnlyLinesAndCubicNurbs = getOpt('bOnlyLinesAndCubicNurbs')
    bTryGetArcs = getOpt('bTryGetArcs')
    bAcceptRational = getOpt('bAcceptRational')
    bTryRebuildOthersUniform = getOpt('bTryRebuildOthersUniform')
    bDebug = getOpt('bDebug')


    sc.doc.Objects.UnselectAll() # For debug.


    #TODO
    #fTol_Proj_Total = 0.1*sc.doc.ModelAbsoluteTolerance if fTol_Proj is None else fTol_Proj
    fTol_Proj_Total = fTol_Proj


    if not bPostProcess:
        return projectCrv_to_1faceBreps_and_planes(
            rgC_In,
            rgBs_1Face,
            rgPlanes,
            vectDir,
            fTol_Proj=fTol_Proj)


    # Wrapping the function since only the tolerance argument will change.
    def projectCrvToGeom_wrapped(tolerance):
        return projectCrv_to_1faceBreps_and_planes(
            rgC_In,
            rgBs_1Face,
            rgPlanes,
            vectDir,
            fTol_Proj=fTol_Proj,
            )


    if bOnlyLinesAndCubicNurbs:
        return projectCrv_Try_to_output_only_lines_and_cubicNurbs(
            fTol_Proj_Total=fTol_Proj,
            projectCrvToGeom_wrapped=projectCrvToGeom_wrapped,
            bDebug=bDebug)


    return projectCrv_Try_to_output_per_arguments(
        fTol_Proj=fTol_Proj,
        projectCrvToGeom_wrapped=projectCrvToGeom_wrapped,
        fTol_MinLength=fTol_MinLength,
        bTryGetArcs=bTryGetArcs,
        bAcceptRational=bAcceptRational,
        bTryRebuildOthersUniform=bTryRebuildOthersUniform,
        bDebug=bDebug)


def duplicateSegments(rgC_In):
    rgCs_Segs = rgC_In.DuplicateSegments()
    if not rgCs_Segs:
        1/0
        rgCs_Segs = [rgC_In.DuplicateCurve()]
    return rgCs_Segs


def processDocObjects(rhObjs_toProj, rhBreps, vect, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bProjectCrvSegs = getOpt('bProjectCrvSegs')
    fTol_Proj = getOpt('fTol_Proj')
    bLoose = getOpt('bLoose')
    bPostProcess = getOpt('bPostProcess')
    fTol_MinLength = getOpt('fTol_MinLength')
    bOnlyLinesAndCubicNurbs = getOpt('bOnlyLinesAndCubicNurbs')
    bTryGetArcs = getOpt('bTryGetArcs')
    bAcceptRational = getOpt('bAcceptRational')
    bTryRebuildOthersUniform = getOpt('bTryRebuildOthersUniform')
    bJoinPerInputCrv = getOpt('bJoinPerInputCrv')
    bDeleteInput = getOpt('bDeleteInput')
    iOutputLayer = getOpt('iOutputLayer')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def coerceBrepObject(rhObj):
        rdObj = rs.coercerhinoobject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj



    rdBs_In = [] # Flat list.
    rgBs_In = [] # Flat list. May be breps of individual faces if the latter were subobject selected.
    rgB1s_1F_PerB_In = [] # Nested lists.
    rgPlanes_PerB_In = [] # Nested lists.


    def collectTargetBrepData(rhBreps, rdBs_In, rgBs_In, rgB1s_1F_PerB_In, rgPlanes_PerB_In):
        for rhB in rhBreps:
            rgB = coerceBrep(rhB)
            if not rgB.IsValid:
                print("Invalid Brep in input, but projection to it" \
                    " will still be attempted." \
                    "  Check results!")

            rdBs_In.append(coerceBrepObject(rhB))
            rgBs_In.append(rgB)

            if not bLoose:
                rgB1s_1F_PerB_In.append([])
                rgPlanes_PerB_In.append([])

                for rgF in rgB.Faces:
                    rgB_1F = rgF.DuplicateFace(duplicateMeshes=False)

                    # Shrink face, otherwise Curve.ProjectToBrep
                    # for some toroidal RevSurfaces may return None.
                    # If this doesn't work, conversion from RevSurface
                    # to NurbsSurface may be a solution.
                    rgB_1F.Faces.ShrinkFaces()

                    rgB1s_1F_PerB_In[-1].append(rgB_1F)
                    rgSrf = rgF.UnderlyingSurface()
                    bSuccess, rgPlane = rgSrf.TryGetPlane(tolerance=2**-32)
                    if bSuccess:
                        rgPlanes_PerB_In[-1].append(rgPlane)
                    else:
                        rgPlanes_PerB_In[-1].append(None)

    collectTargetBrepData(rhBreps, rdBs_In, rgBs_In, rgB1s_1F_PerB_In, rgPlanes_PerB_In)

    rgCs_Proj_All = []

    attr_Out = rd.ObjectAttributes()

    g_Out_All = []

    if iOutputLayer == Opts.listValues['iOutputLayer'].index('Current'):
        attr_Out.LayerIndex = sc.doc.Layers.CurrentLayerIndex

    len_rhObjs0 = len(rhObjs_toProj)
    idxs_AtTenths = [int(round(0.1*i*len(rhObjs_toProj),0)) for i in range(10)]

    for iO, rhObj_toProj in enumerate(rhObjs_toProj):
        if sc.escape_test(False):
            print("Projecting interrupted by user.")
            return

        if len_rhObjs0 > 10:
            if iO in idxs_AtTenths:
                Rhino.RhinoApp.SetCommandPrompt("Processing at {:d}% of {} object ...".format(
                    int(100.0 * (iO+1) / len_rhObjs0), len_rhObjs0))
        elif len_rhObjs0 > 1:
            Rhino.RhinoApp.SetCommandPrompt(
                "Projecting {} of {} objects ({})".format(
                    iO+1, len_rhObjs0, rs.coerceguid(rhObj_toProj)))
        else:
            Rhino.RhinoApp.SetCommandPrompt("Projecting object ...")


        rdObj_In = rs.coercerhinoobject(rhObj_toProj)

        if iOutputLayer == Opts.listValues['iOutputLayer'].index('Input'):
            attr_Out.LayerIndex = rdObj_In.Attributes.LayerIndex


        # PointObjects.
        if rdObj_In.ObjectType == rd.ObjectType.Point:
            rgPt0 = rdObj_In.Geometry

            gPts_from1Pt_Out = []

            for rdB, rgB in zip(rdBs_In, rgBs_In):
                if iOutputLayer == Opts.listValues['iOutputLayer'].index('TargetObject'):
                    attr_Out.LayerIndex = rdB.Attributes.LayerIndex

                pts_Proj = rg.Intersect.Intersection.ProjectPointsToBreps(
                    breps=[rgB],
                    points=[rgPt0.Location],
                    direction=vect,
                    tolerance=fTol_Proj)
                if not pts_Proj: continue

                for pt in pts_Proj:
                    gPt_Out = sc.doc.Objects.AddPoint(pt, attributes=attr_Out)
                    if gPt_Out != Guid.Empty:
                        gPts_from1Pt_Out.append(gPt_Out)

            if gPts_from1Pt_Out:
                g_Out_All.extend(gPts_from1Pt_Out)
            if bDeleteInput:
                sc.doc.Objects.Delete(item=rdObj_In)


        if rdObj_In.ObjectType not in (rd.ObjectType.Curve, rd.ObjectType.Brep):
            raise Exception("{} is not allowed as input to project.".format(rdObj_In.GetType().Name))


        # CurveObjects and BrepEdges.
        if rdObj_In.ObjectType == rd.ObjectType.Curve:
            rgC_In = rdObj_In.Geometry
        elif rdObj_In.ObjectType == rd.ObjectType.Brep:
            # For edges.
            rgC_In = rs.coercecurve(rhObj_toProj)


        # Prepd = Preprocessed
        rgCs_Prepd_perC_In = duplicateSegments(rgC_In) if bProjectCrvSegs else [rgC_In.DuplicateCurve()]

        rgCs_Res_perC_In = []

        for rgC_Prepd_perC_In in rgCs_Prepd_perC_In:
            if bLoose:
                for rdB_In, rgB_In in zip(rdBs_In, rgBs_In):
                    if iOutputLayer == Opts.listValues['iOutputLayer'].index('TargetObject'):
                        attr_Out.LayerIndex = rdB_In.Attributes.LayerIndex

                    rgCs_Projctd_from1Prepd_1B = projectCurve_Loose(
                        rgC_Prepd_perC_In,
                        rgB_In,
                        vect,
                        fTol_Proj=fTol_Proj,
                        bDebug=bDebug)

                    if rgCs_Projctd_from1Prepd_1B:
                        rgCs_Res_perC_In.extend(rgCs_Projctd_from1Prepd_1B)
            else:
                for rdB_In, rgBs_1Face, rgPlanes in zip(rdBs_In, rgB1s_1F_PerB_In, rgPlanes_PerB_In):
                    if iOutputLayer == Opts.listValues['iOutputLayer'].index('TargetObject'):
                        attr_Out.LayerIndex = rdB_In.Attributes.LayerIndex

                    rgCs_Projctd_from1Prepd_1B = projectCurve_NotLoose(
                        rgC_Prepd_perC_In,
                        rgBs_1Face,
                        rgPlanes,
                        vect,
                        fTol_Proj=fTol_Proj,
                        bPostProcess=bPostProcess,
                        fTol_MinLength=fTol_MinLength,
                        bOnlyLinesAndCubicNurbs=bOnlyLinesAndCubicNurbs,
                        bTryGetArcs=bTryGetArcs,
                        bAcceptRational=bAcceptRational,
                        bTryRebuildOthersUniform=bTryRebuildOthersUniform,
                        bDebug=bDebug)

                    if rgCs_Projctd_from1Prepd_1B:
                        rgCs_Res_perC_In.extend(rgCs_Projctd_from1Prepd_1B)


        if not rgCs_Res_perC_In:
            continue

        if bJoinPerInputCrv:
            rgCs_Out = joinCurves(
                rgCs_Res_perC_In,
                tolerance=fTol_MinLength,
                bDebug=bDebug)
            #map(sc.doc.Objects.AddCurve, rgCs_Out); sc.doc.Views.Redraw(); 1/0
            if bPostProcess:
                # TODO: WIP
                if rgC_In.IsPeriodic:
                    print("Periodic input curve found.")
                #    if convertCrvsForPeriodicInList(rgCs_toMod):
                #        bModified = True
                #        lenList = len(rgCs_toMod)

                # Even if removeShortSegmentsInEachCrv_inList was previously run,
                # new short segments may have been produced by JoinCurves.
                # This problem also occurs using _Project.
                if removeShortSegmentsInEachCrv_inList(rgCs_Out, fTol_MinLength):
                    if bDebug:
                        print("Short segment(s) removed.")
        else:
            rgCs_Out = rgCs_Res_perC_In

        gCs_Res_fromC_In = []

        for rgC in rgCs_Out:
            gC_Out = sc.doc.Objects.AddCurve(rgC, attributes=attr_Out)
            if gC_Out != Guid.Empty:
                gCs_Res_fromC_In.append(gC_Out)

        g_Out_All.extend(gCs_Res_fromC_In)

        if bDeleteInput and (rdObj_In.ObjectType == rd.ObjectType.Curve):
            # Skips deleting BreObjects.
            sc.doc.Objects.Delete(item=rdObj_In)

        # End of loop through rhObjs_toProj.


    for rgB in rgBs_In:
        if not rgB.IsDocumentControlled:
            # Brep was created from a subobject-selected face.
            rgB.Dispose()

    if rgB1s_1F_PerB_In:
        for rgB1s_1F in rgB1s_1F_PerB_In:
            for rgB in rgB1s_1F:
                rgB.Dispose()


    return g_Out_All


def getDirectionVector(iDirection):
    if Opts.listValues['iDirection'][iDirection] == 'CPlaneX':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().XAxis
    if Opts.listValues['iDirection'][iDirection] == 'CPlaneY':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().YAxis
    if Opts.listValues['iDirection'][iDirection] == 'CPlaneZ':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    elif Opts.listValues['iDirection'][iDirection] == 'WorldX':
        return rg.Vector3d.XAxis
    elif Opts.listValues['iDirection'][iDirection] == 'WorldY':
        return rg.Vector3d.YAxis
    elif Opts.listValues['iDirection'][iDirection] == 'WorldZ':
        return rg.Vector3d.ZAxis
    elif Opts.listValues['iDirection'][iDirection] == 'View':
        return sc.doc.Views.ActiveView.ActiveViewport.GetCameraFrame()[1].ZAxis
    elif Opts.listValues['iDirection'][iDirection] == 'Custom':
        return Opts.values['vectCustom']


def main():
    
    objrefs_Crvs_Pts = getInput(
        [],
        "Select curves and points to project",
        rd.ObjectType.Curve | rd.ObjectType.Point)
    if objrefs_Crvs_Pts is None: return

    #bProjectCrvSegs = Opts.values['bProjectCrvSegs']
    #fTol_MinLength = Opts.values['fTol_MinLength']
    #bLoose = Opts.values['bLoose']
    #fTol_Proj = Opts.values['fTol_Proj']
    #bPostProcess = Opts.values['bPostProcess']
    #bOnlyLinesAndCubicNurbs = Opts.values['bOnlyLinesAndCubicNurbs']
    #bTryGetArcs = Opts.values['bTryGetArcs']
    #bAcceptRational = Opts.values['bAcceptRational']
    #bTryRebuildOthersUniform = Opts.values['bTryRebuildOthersUniform']
    #bDeleteInput = Opts.values['bDeleteInput']
    #iOutputLayer = Opts.values['iOutputLayer']
    #iDirection = Opts.values['iDirection']
    #bEcho = Opts.values['bEcho']
    #bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()

    rdCrvs_toHighlight = [o.Object() for o in objrefs_Crvs_Pts]

    print("Now, select breps and faces ...")

    objrefs_Breps = getInput(
        rdCrvs_toHighlight,
        "Select surfaces and polysurfaces to project onto",
        rd.ObjectType.Brep)
    if objrefs_Breps is None: return

    # Determine vector before main routine.
    vect = getDirectionVector(Opts.values['iDirection'])
    #if Opts.listValues['iDirection'][iDirection] == 'CPlaneZ':
    #    vect = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane().ZAxis
    #elif Opts.listValues['iDirection'][iDirection] == 'WorldZ':
    #    vect = rg.Vector3d.ZAxis
    #elif Opts.listValues['iDirection'][iDirection] == 'View':
    #    rc = sc.doc.Views.ActiveView.ActiveViewport.GetCameraFrame()
    #    if not rc[0]: return
    #    vect = rc[1].ZAxis
    #elif Opts.listValues['iDirection'][iDirection] == 'Custom':
    #    rc = rs.GetLine(
    #        mode=1, point=None,
    #        message1='Projection direction',
    #        message3='Second direction point',
    #        )
    #    if not rc: return
    #    vect = rg.Vector3d(rc[1] - rc[0])


    bProjectCrvSegs = Opts.values['bProjectCrvSegs']
    iDirection = Opts.values['iDirection']
    fTol_Proj = Opts.values['fTol_Proj']
    bLoose = Opts.values['bLoose']
    bPostProcess = Opts.values['bPostProcess']
    fTol_MinLength = Opts.values['fTol_MinLength']
    bOnlyLinesAndCubicNurbs = Opts.values['bOnlyLinesAndCubicNurbs']
    bTryGetArcs = Opts.values['bTryGetArcs']
    bAcceptRational = Opts.values['bAcceptRational']
    bTryRebuildOthersUniform = Opts.values['bTryRebuildOthersUniform']
    bJoinPerInputCrv = Opts.values['bJoinPerInputCrv']
    bDeleteInput = Opts.values['bDeleteInput']
    iOutputLayer = Opts.values['iOutputLayer']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if not Opts.values['bDebug']:
        sc.doc.Views.RedrawEnabled = False

    gs_Res = processDocObjects(
        objrefs_Crvs_Pts,
        objrefs_Breps,
        vect,
        bProjectCrvSegs=bProjectCrvSegs,
        fTol_Proj=fTol_Proj,
        bLoose=bLoose,
        bPostProcess=bPostProcess,
        fTol_MinLength=fTol_MinLength,
        bOnlyLinesAndCubicNurbs=bOnlyLinesAndCubicNurbs,
        bTryGetArcs=bTryGetArcs,
        bAcceptRational=bAcceptRational,
        bTryRebuildOthersUniform=bTryRebuildOthersUniform,
        bJoinPerInputCrv=bJoinPerInputCrv,
        bDeleteInput=bDeleteInput,
        iOutputLayer=iOutputLayer,
        iDirection=iDirection,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if Opts.values['bEcho']:
        if len(objrefs_Crvs_Pts) == len(gs_Res) == 1:
            print("Object projected to one object.".format(
                len(objrefs_Crvs_Pts)))
        elif len(objrefs_Crvs_Pts) == len(gs_Res):
            print("{} objects projected to same number of objects.".format(
                len(objrefs_Crvs_Pts)))
        elif len(gs_Res) == 0:
            print("The projection missed the selected objects.".format(
                len(objrefs_Crvs_Pts), len(gs_Res)))
        else:
            print("{} objects projected to {} objects.".format(
                len(objrefs_Crvs_Pts), len(gs_Res)))

        iCt_Crvs = 0
        iCt_Closed = 0

        for gRes in gs_Res:
            if rs.IsCurve(gRes):
                iCt_Crvs += 1
                if rs.IsCurveClosed(gRes):
                    iCt_Closed += 1

        if iCt_Crvs:
            if iCt_Closed == iCt_Crvs:
                print("All curves are closed.")
            elif iCt_Closed == 0:
                print("All curves are open.")
            else:
                print("{} curves are open. {} are closed.".format(iCt_Crvs-iCt_Closed, iCt_Closed))


    if gs_Res:
        sc.doc.Objects.UnselectAll()
        sc.doc.Objects.Select(objectIds=gs_Res)
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()