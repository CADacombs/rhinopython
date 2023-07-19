"""
This script is an alternative to ExtrudeCrvTapered.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
181120: Created.
190528: Added Opts and support for edge input.
190529: Added bSimplifySrf.
190712-14: Refactored.  Polished the UX.
190724: Bug fixes in UI.
190924-200121, 220328: Import-related updates.
220809: Added CreateFromTaperedExtrudeWithRef.  Refactored.
230425-27: Branched from spb_Brep_createFromTaperedExtrude.py
    Recreated from with different input scenario and now previews via DrawConduit instead of actual DocObjects.
230719: Minor efficency improvements.

TODO: Add my variable taper routine from spb_TaperFromPathOnVerticalFrames.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
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


    key = 'bExplodePolyCrv'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitCrvAtHighKnotMulties'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSimplifyCrv'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSimplifyCrvTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bExtrudeToPlane_Loose'; keys.append(key)
    values[key] = False
    names[key] = 'Loose'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bNumInputForDist_NotAngle'; keys.append(key)
    values[key] = True
    names[key] = 'NumberEntry'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Angle', 'Dist')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistance'; keys.append(key)
    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        values[key] = 1.0
    else:
        values[key] = 10.0 * Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fAngle_Start_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'DraftAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCPlane'; keys.append(key)
    values[key] = True
    names[key] = 'PosDirPerZAxisOf'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iCornerType'; keys.append(key)
    values[key] = 2
    listValues[key] = Enum.GetNames(rg.ExtrudeCornerType)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bJoinForCorners'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fAngleTol_Deg'; keys.append(key)
    values[key] = sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__, sc.doc.Name)

    key = 'bSimplifySrf'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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

        if key == 'iGContinuity':
            cls.values[key] = cls.riOpts[key].CurrentValue
        else:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            elif key in cls.listValues:
                cls.values[key] = idxList
            else:
                return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def _getInput_Click():
    """
    Click to toggle angle and/or direction with optional input.

    Returns
        True: To recalculate and reloop
        False: To not recalculate and break out of loop with current output.
        None: To not recalculate and return without output.
    """

    gp = ri.Custom.GetPoint()

    #TODO
    bDefaultsForDir_NotAngle = True

    gp.SetCommandPrompt("Pick to flip angle")

    gp.AcceptNumber(True, acceptZero=True)
    gp.AcceptNothing(True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(gp, key)

    #if not bDefaultsForDir_NotAngle and Opts.values['bVariableAngle']:
    #    key = 'FlipAngle'; idxs_Opts[key] = go.AddOption(key)
    #    key = 'SwapAngles'; idxs_Opts[key] = go.AddOption(key)

    addOption('bExplodePolyCrv')
    addOption('bSimplifyCrv')
    if Opts.values['bSimplifyCrv']:
        addOption('fSimplifyCrvTol')
    addOption('bSplitCrvAtHighKnotMulties')
    addOption('bExtrudeToPlane_Loose')
    addOption('bNumInputForDist_NotAngle')
    addOption('fDistance')
    addOption('fAngle_Start_Deg')
    key = 'FlipDir'; idxs_Opts[key] = gp.AddOption(key)
    idxs_Opts['FlipAngle'] = gp.AddOption('FlipAngle')
    addOption('bCPlane')
    if not Opts.values['bExtrudeToPlane_Loose']:
        addOption('iCornerType')
    addOption('fDistTol')
    addOption('fAngleTol_Deg')
    addOption('bSimplifySrf')
    addOption('bEcho')
    addOption('bDebug')


    while True:

        res = gp.Get()

        if res == ri.GetResult.Cancel:
            gp.Dispose()
            return

        if res == ri.GetResult.Nothing:
            gp.Dispose()
            return False

        if res == ri.GetResult.Point:
            Opts.riOpts['fAngle_Start_Deg'].CurrentValue = -Opts.riOpts['fAngle_Start_Deg'].CurrentValue
            Opts.setValue('fAngle_Start_Deg')
            gp.Dispose()
            return True

        if res == ri.GetResult.Number:
            if bDefaultsForDir_NotAngle:
                key = 'fDistance'
            else:
                key = 'fAngle_Start_Deg'
            key = 'fDistance' if Opts.values['bNumInputForDist_NotAngle'] else 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = gp.Number()
            Opts.setValue(key)
            gp.Dispose()
            return True

        # An option was selected.

        if gp.OptionIndex() == idxs_Opts['FlipAngle']:
            key = 'fAngle_Start_Deg'
            Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
            Opts.setValue(key)
            return True

        if gp.OptionIndex() == idxs_Opts['FlipDir']:
            key = 'fDistance'
            Opts.riOpts[key].CurrentValue = -Opts.riOpts[key].CurrentValue
            Opts.setValue(key)
            return True

        key = 'iCornerType'
        if key in idxs_Opts and gp.OptionIndex() == idxs_Opts[key]:
            Opts.setValue(key, gp.Option().CurrentListOptionIndex)
            return True

        key = 'fDistTol'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            Opts.setValue(key)
            return True

        key = 'fAngleTol_Deg'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            Opts.setValue(key)
            return True

        for key in idxs_Opts:
            if gp.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, gp.Option().CurrentListOptionIndex)
                return True

        raise Exception("L312: What happened?")

        #go.Dispose()
        #return True


def _hasHighKnotMultiplicities(nurbsObj, iDir):

    if isinstance(nurbsObj, rg.NurbsCurve):
        nc = nurbsObj
        bIsPeriodic = nc.IsPeriodic
        bIsClosed = nc.IsClosed
        iDegree = nc.Degree
        knots = nc.Knots
    elif isinstance(nurbsObj, rg.NurbsSurface):
        if iDir not in (0,1):
            raise ValueError("{} passed as surface direction.".format(iDir))
        ns = nurbsObj
        bIsPeriodic = ns.IsPeriodic(iDir)
        bIsClosed = ns.IsClosed(iDir)
        iDegree = ns.Degree(iDir)
        knots = ns.KnotsV if iDir else ns.KnotsU
    else:
        raise Exception("{} is not valid input.  Need a NurbsCurve or NurbsSurface.".format(
            nurbsObj.GetType().Name))

    if bIsPeriodic:
        iKs = range(knots.Count)
    elif bIsClosed:
        iKs = range(knots.Count - iDegree)
    else:
        iKs = range(iDegree, knots.Count - iDegree)

    for iK in iKs:
        if knots.KnotMultiplicity(iK) > iDegree-2:
            return True

    return False


def _getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None

    return rc[1]


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


def _simplifyCrv(rgCrv_In, fTol, bDebug=False):
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


    # Check if non-uniform with only simple internal knots.
    ms = knotMultiplicityList(nc_WIP.Knots)
    if not(nc_WIP.IsClosed and nc_WIP.IsPeriodic):
        ms = ms[1:-1]
    if all([m == 1 for m in ms]):
        return nc_WIP


    # "_MakeUniform"
    if nc_WIP.Degree == 3:

        if nc_WIP.IsPeriodic:
            nc_WIP.Knots.CreatePeriodicKnots(knotSpacing=1.0) # Modifies existing Knots.
        else:
            nc_WIP.Knots.CreateUniformKnots(knotSpacing=1.0) # Modifies existing Knots.

        nc_WIP.Domain = rgCrv_In.Domain

        dev = _getDistancesBetweenCurves(rgCrv_In, nc_WIP)
        if dev <= fTol:
            return nc_WIP
        if bDebug: print("'MakeUniform' routine result is not within {}.".format(fTol))

    nc_WIP.Dispose()

    return _rebuildCrv(rgCrv_In, fTol)


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
            if nc.Knots.KnotMultiplicity(iK) > nc.Degree-2:
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


def _processCrvInput(rgCrv_In, bSimplifyCrv, bExplodePolyCrv, bSplitCrvAtHighKnotMulties, bMakeDeformable, fTol, bDebug=False):
    
    # Explode PolyCrvs before simplification, etc.
    if bExplodePolyCrv and isinstance(rgCrv_In, rg.PolyCurve):
        ncs_WIP = [_.ToNurbsCurve() for _ in rgCrv_In.Explode()]
    else:
        ncs_WIP = [rgCrv_In.ToNurbsCurve()]

    # Simplify to Bezier or simpler NURBS before splitting at high knot multiplicities.
    if bSimplifyCrv:
        for i in range(len(ncs_WIP)):
            rc = _simplifyCrv(ncs_WIP[i], fTol, bDebug)
            if rc:
                ncs_WIP[i].Dispose()
                ncs_WIP[i] = rc

    if bSplitCrvAtHighKnotMulties:
        rc = _splitNurbsCrvsAtPolyKnots(ncs_WIP)
        for _ in ncs_WIP: _.Dispose()
        ncs_WIP = rc

    if bMakeDeformable:
        for i in range(len(ncs_WIP)):
            if ncs_WIP[i].Degree < 3 and ncs_WIP[i].Points.Count < 4:
                ncs_WIP[i].IncreaseDegree(3)

    return ncs_WIP


def _fitSrf(rgSrf_In, fTol_FitSrf):
    ns_FromIn = rgSrf_In.ToNurbsSurface()

    sEval = "_hasHighKnotMultiplicities(ns_FromIn, 0)"; print("{}: {}".format(sEval, eval(sEval)))
    sEval = "_hasHighKnotMultiplicities(ns_FromIn, 1)"; print("{}: {}".format(sEval, eval(sEval)))

    bU_HighMs = _hasHighKnotMultiplicities(ns_FromIn, 0)
    bV_HighMs = _hasHighKnotMultiplicities(ns_FromIn, 1)
    if not (bU_HighMs or bV_HighMs):
        ns_FromIn.Dispose()
        return

    ns_Fit = ns_FromIn.Fit(
        uDegree=ns_FromIn.Degree(0),
        vDegree=ns_FromIn.Degree(1),
        fitTolerance=fTol_FitSrf)

    if ns_Fit is None:
        print("A surface could not be Fit within {}.".format(fTol_FitSrf))
        ns_FromIn.Dispose()
        return
    
    if ns_Fit.EpsilonEquals(ns_FromIn, epsilon=1e-6):
        ns_FromIn.Dispose()
        ns_Fit.Dispose()
        print("A surface Fit didn't produce a new surface.")
        return

    ns_FromIn.Dispose()

    return ns_Fit


def createBreps(rgCrvs_ToExtrude, **kwargs):
    """
    rgCrvs_ToExtrude: Multiple curves should be passed so that output Breps can be joined.
    bExtrudeToPlane_Loose = getOpt('bExtrudeToPlane_Loose')
    fAngle_Start_Deg = getOpt('fAngle_Start_Deg')
    bCPlane = getOpt('bCPlane')
    fDistance = getOpt('fDistance')
    iCornerType = getOpt('iCornerType')
    fDistTol = getOpt('fDistTol')
    bSimplifySrf = getOpt('bSimplifySrf')
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bExtrudeToPlane_Loose = getOpt('bExtrudeToPlane_Loose')
    fAngle_Start_Deg = getOpt('fAngle_Start_Deg')
    bCPlane = getOpt('bCPlane')
    fDistance = getOpt('fDistance')
    iCornerType = getOpt('iCornerType')
    fDistTol = getOpt('fDistTol')
    bSimplifySrf = getOpt('bSimplifySrf')
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if bCPlane:
        view = sc.doc.Views.ActiveView
        plane = view.ActiveViewport.ConstructionPlane()
        direction = plane.Normal
    else:
        plane = rg.Plane.WorldXY
        direction = rg.Vector3d.ZAxis

    basePoint = rg.Point3d.Origin

    if bSimplifySrf:
        fTol_Extrude = fTol_FitSrf = fDistTol / 2.0
    else:
        fTol_Extrude = fDistTol


    rgBreps2_Extrds_ToJoin = []

    for rgCrv_ToExtrude in rgCrvs_ToExtrude:
        if bExtrudeToPlane_Loose:
            rgCrv_ToExtrude.Reverse()
            plane_Temp = rg.Plane(plane)
            plane_Temp.Translate(rgCrv_ToExtrude.PointAtStart - plane.Origin)
            rgBreps1_Extrds_Raw = rg.Brep.CreateFromTaperedExtrudeWithRef(
                curve=rgCrv_ToExtrude,
                direction=direction,
                distance=fDistance,
                draftAngle=Rhino.RhinoMath.ToRadians(fAngle_Start_Deg),
                plane=plane_Temp,
                tolerance=fTol_Extrude)
        else:
            rgBreps1_Extrds_Raw = rg.Brep.CreateFromTaperedExtrude(
                    curveToExtrude=rgCrv_ToExtrude,
                    distance=fDistance,
                    direction=direction,
                    basePoint=basePoint,
                    draftAngleRadians=Rhino.RhinoMath.ToRadians(fAngle_Start_Deg),
                    cornerType=Enum.ToObject(rg.ExtrudeCornerType, iCornerType),
                    tolerance=fTol_Extrude,
                    angleToleranceRadians=Rhino.RhinoMath.ToRadians(fAngleTol_Deg))
        if not rgBreps1_Extrds_Raw:
            print("An Extrude could not be created.")
            continue

        if not bSimplifySrf:
            rgBreps2_Extrds_ToJoin.extend(rgBreps1_Extrds_Raw)
        else:
            for rgBrep1_Extrd_Raw in rgBreps1_Extrds_Raw:
                for srf_ToFit in rgBrep1_Extrd_Raw.Surfaces:
                    rc = _fitSrf(srf_ToFit, 0.1*sc.doc.ModelAbsoluteTolerance)
                    if rc:
                        rgBreps2_Extrds_ToJoin.append(rc.ToBrep())
                        rc.Dispose()
                    else:
                        rgBreps2_Extrds_ToJoin.append(srf_ToFit.ToBrep())
                rgBrep1_Extrd_Raw.Dispose()

    if len(rgBreps2_Extrds_ToJoin) == 1:
        return rgBreps2_Extrds_ToJoin

    rgBreps_Extrds_Joined = rg.Brep.JoinBreps(
        rgBreps2_Extrds_ToJoin,
        tolerance=0.1*sc.doc.ModelAbsoluteTolerance) # A tighter tolerance can help identify problem areas.

    if rgBreps_Extrds_Joined:
        iCt_F = sum([rgB.Faces.Count for rgB in rgBreps_Extrds_Joined])

    return rgBreps_Extrds_Joined


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

    rc, objrefs = ri.RhinoGet.GetMultipleObjects(
        prompt="Select curves to Extrude",
        acceptNothing=False,
        filter=rd.ObjectType.Curve)

    if rc != Rhino.Commands.Result.Success: return

    

    sk_conduit = 'conduit({})'.format(__file__) # StickyKey
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit

    rgCs_In_NoNesting = []
    for o in objrefs:
        rgC = o.Curve()
        if isinstance(rgC, rg.BrepEdge):
            rgC = rgC.DuplicateCurve()
        if isinstance(rgC, rg.PolyCurve):
            rgC.RemoveNesting()
        rgCs_In_NoNesting.append(rgC)


    while True:
        sc.escape_test()

        bExplodePolyCrv = Opts.values['bExplodePolyCrv']
        bSplitCrvAtHighKnotMulties = Opts.values['bSplitCrvAtHighKnotMulties']
        bSimplifyCrv = Opts.values['bSimplifyCrv']
        fSimplifyCrvTol = Opts.values['fSimplifyCrvTol']
        bExtrudeToPlane_Loose = Opts.values['bExtrudeToPlane_Loose']
        fAngle_Start_Deg = Opts.values['fAngle_Start_Deg']
        bCPlane = Opts.values['bCPlane']
        fDistance = Opts.values['fDistance']
        iCornerType = Opts.values['iCornerType']
        fDistTol = Opts.values['fDistTol']
        fAngleTol_Deg = Opts.values['fAngleTol_Deg']
        bSimplifySrf = Opts.values['bSimplifySrf']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        #TODO
        #bMakeDeformable = bAlignEndDirs or (fAngle_End_Deg and fAngle_End_Deg != fAngle_Start_Deg)
        bMakeDeformable = False

        if bExplodePolyCrv or bSimplifyCrv or bSplitCrvAtHighKnotMulties:
            ncs_ToExtrude = []
            for i, rgC in enumerate(rgCs_In_NoNesting):
                rc = _processCrvInput(
                    rgC,
                    bSimplifyCrv,
                    bExplodePolyCrv,
                    bSplitCrvAtHighKnotMulties,
                    bMakeDeformable,
                    fSimplifyCrvTol,
                    bDebug)
                if rc:
                    ncs_ToExtrude.extend(rc)
                else:
                    ncs_ToExtrude.append(rgCs_In_NoNesting)
        else:
            ncs_ToExtrude = rgCs_In_NoNesting[:]


        if not bExtrudeToPlane_Loose and iCornerType:
            rc = rg.Curve.JoinCurves(ncs_ToExtrude, sc.doc.ModelAbsoluteTolerance, preserveDirection=True)
            if rc and len(rc) < len(ncs_ToExtrude):
                ncs_ToExtrude = rc

        rgBs_Joined = []


        #for nc in ncs_ToExtrude: sc.doc.Objects.AddCurve(nc)
        #sc.doc.Views.Redraw()
        #return

        rgBs_Joined = createBreps(
            rgCrvs_ToExtrude=ncs_ToExtrude,
            bExtrudeToPlane_Loose=bExtrudeToPlane_Loose,
            fDistance=fDistance,
            bCPlane=bCPlane,
            fAngle_Start_Deg=fAngle_Start_Deg,
            iCornerType=iCornerType,
            fDistTol=fDistTol,
            fAngleTol_Deg=fAngleTol_Deg,
            bEcho=bEcho,
            bDebug=bDebug)

        if rgBs_Joined:
            conduit.breps = rgBs_Joined
            conduit.Enabled = True
            sc.doc.Views.Redraw()

        rc = _getInput_Click()

        conduit.Enabled = False

        if rc is None:
            for _ in rgBs_Joined: _.Dispose()
            return

        if not rc:
            return (
                rgBs_Joined,
                bEcho)

        for _ in rgBs_Joined: _.Dispose()





def main():
    """
    Gets rg.Breps and adds to document.
    """

    rc = _createGeometryInteractively()
    if rc is None: return

    rgBs, bEcho = rc

    gBs = []
    for rgB in rgBs:
        gB = sc.doc.Objects.AddBrep(rgB)
        rgB.Dispose()
        if gB != gB.Empty:
            gBs.append(gB)

    if bEcho:
        sOut = []
        if gBs:
            iCt_Fs = sum(sc.doc.Objects.FindId(g).BrepGeometry.Faces.Count for g in gBs)
            sOut.append(("{} brep(s) with {} faces".format(
                len(gBs), iCt_Fs)))
        print("Added {}".format(", ".join(sOut)))


if __name__ == '__main__': main()