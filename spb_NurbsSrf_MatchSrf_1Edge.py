"""
Complement to _MatchSrf.
This script, when provided only inputs of natural (underlying surface) edges,
produces a surfaces that follow the input knot and point structures closely,
with continuity priority G0 then G1 then G2.


To reduce command options, a MultipleMatches option will be implemented in another script.

Options not found in _MatchSrf:
    UseUnderlyingIsoCrvs: When 'Yes', entire isocurve of the reference underlying surface is matched.
    IncreaseAsNeeded: When set to 'SpanCt', degree is maintained and knots are added as needed.
    Replace: When set to 'No', a new surface is added instead of replacing the original.

Limitations:
    Does not consider control point weights.
    Does not G2 match to non-isocurve trimmed edges except for planar surfaces.
    None of the following _MatchSrf options:
        ChainEdges, OnSurface, MatchByClosestPoints, RefineMatch,
        AverageSurfaces, Tolerances, IsocurveDirection

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

from __future__ import print_function

"""
200626: Started and 1st working version.
201207: Added G2 matching.  Various UI improvements.
201208: Now will increase degree and point counts as needed for matching.
201216: Now will add only needed control points in continuity direction when surface begins single-spanned.
210217: Fixed bug in simple knot multiplier calculation.
210402: Now, similar to _MatchSrf, will match according to pick point locations,
        not auto-align sides.  Auto-align is still available in createSurface.
        Now matches to surface of lesser degree.
        Bug fixes.
211003: Changed verbage of a command option.
211008-25: Bug fix when transferring knot vector.
        Refactored for easier use as a module by other scripts.
211102: Now, simplification is attempted on non-uniform input curves in createTanSrfFromEdge.
211106: Now, SumSurfaces are allowed as input for either surface.
        Now, full Surfaces but with split edges are allowed for the surface to modify.
211117: Minor modifications in createTanSrfFromEdge.  Added print_function from __future__.
211208: Bug fix for when PolyCurves are present during reference selection.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


W = rg.IsoStatus.West
S = rg.IsoStatus.South
E = rg.IsoStatus.East
N = rg.IsoStatus.North


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iContinuity'; keys.append(key)
    values[key] = 1
    names[key] = 'Mode'
    listValues[key] = 'G0', 'G1', 'G2'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iPreserveOtherEnd'; keys.append(key)
    values[key] = 2
    listValues[key] = 'None', 'G0', 'G1', 'G2'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseUnderlyingIsoCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMaintainDegree'; keys.append(key)
    values[key] = True
    names[key] = 'IncreaseAsNeeded'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Degree', 'SpanCt')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
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

    key = 'bAddRefs'; keys.append(key)
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList):
        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def addGeoms(geoms, bRedraw=True):
    """ For debugging. """

    if not hasattr(geoms, '__iter__'):
        geoms = [geoms]

    for geom in geoms:
        if isinstance(geom, tuple):
            gOut = sc.doc.Objects.AddSurface(geom[1])
        elif isinstance(geom, rg.Curve):
            gOut = sc.doc.Objects.AddCurve(geom)
        elif isinstance(geom, rg.Surface):
            gOut = sc.doc.Objects.AddSurface(geom)
        elif isinstance(geom, rg.Point3d):
            gOut = sc.doc.Objects.AddPoint(geom)
        elif isinstance(geom, rg.Plane):
            intrvl = rg.Interval(-sc.doc.ModelAbsoluteTolerance*1000.0, sc.doc.ModelAbsoluteTolerance*1000.0)
            psrf = rg.PlaneSurface(geom, intrvl, intrvl)
            gOut = sc.doc.Objects.AddSurface(psrf)
        else:
            raise ValueError("Method to add {} missing from addGeoms.".format(geom))
    if bRedraw: sc.doc.Views.Redraw()
    return gOut


def containsShortOrCollapsedBorders(srf, tol_short=None):
    """ Variation on checking for singular BrepTrims. """
    if tol_short is None: tol_short = sc.doc.ModelAbsoluteTolerance
    pts = []
    pts.append(srf.PointAt(srf.Domain(0).T0, srf.Domain(1).T0))
    pts.append(srf.PointAt(srf.Domain(0).T0, srf.Domain(1).T1))
    pts.append(srf.PointAt(srf.Domain(0).T1, srf.Domain(1).T1))
    pts.append(srf.PointAt(srf.Domain(0).T1, srf.Domain(1).T0))

    for i in 0,1,2:
        for j in range(i+1, 3+1):
            if pts[i].DistanceTo(pts[j]) < sc.doc.ModelAbsoluteTolerance:
                return True

    return False


def isFaceSupported(rgF, bEcho):
    rgB = rgF.Brep

    if rgB.Faces.Count > 1:
        print("Face is not from a monoface brep.")
        return False

    rgF = rgB.Faces[0]
    rgS = rgF.UnderlyingSurface()

    #if not isinstance(rgS, rg.NurbsSurface):
    #    print("Underlying surface is a {}.".format(rgS.GetType().Name))
    #    return False

    ns = rgS

    if ns.IsClosed(0) or ns.IsClosed(1):
        print("Closed surfaces are not yet supported.")
        return False

    if not rgF.IsSurface:
        print("Face is not a full surface.")
        return False

    if not rgB.Trims.Count == 4:
        print("There are {} trims on this face." \
            "  Only natural edges will remain.".format(
                rgB.Trims.Count))
    #    return False

    # Reject any NS with singular trims.
    iCt_Singulars = 0
    for iT in 0,1,2,3:
        if rgB.Trims[iT].TrimType == rg.BrepTrimType.Singular:
            iCt_Singulars += 1
    if iCt_Singulars > 0:
        print("There are {} singular trims on this surface.".format(
            iCt_Singulars))
        return False

    return True


def printObjRefObjs(objref):
    for sMethod in (
        'Brep', 'Curve', 'CurveParameter', 'Edge', 'Face', 'Geometry',
        'Object', 'Point', 'SelectionMethod', 'SelectionPoint', 'Surface',
        'SurfaceParameter', 'Trim'
    ):
        print("  .{}(): {}".format(sMethod, eval("objref."+sMethod+"()")))


def getInput_ToModify():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select untrimmed surface edge to change")

    go.GeometryFilter = rd.ObjectType.EdgeFilter


    def customGeometryFilter(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index)

        if not isinstance(geom, rg.BrepTrim):
            return False

        rgT = geom

        if rgT.IsoStatus not in (W, S, E, N):
            return False

        rgB = rgT.Brep
        if rgB.Faces.Count > 1:
            return False

        rgF = rgB.Faces[0]
        rgS = rgF.UnderlyingSurface()

        #if not isinstance(rgS, rg.NurbsSurface):
        #    return False

        ns = rgS

        if ns.IsClosed(0) or ns.IsClosed(1):
            return False

        if not rgF.IsSurface:
            return False

        #if not rgB.Trims.Count == 4:
        #    return False

        if containsShortOrCollapsedBorders(ns):
            return False

        return True

    go.SetCustomGeometryFilter(customGeometryFilter)

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iContinuity')
        addOption('iPreserveOtherEnd')
        addOption('bUseUnderlyingIsoCrvs')
        addOption('bMaintainDegree')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')
        if Opts.values['bDebug']:
            addOption('bAddRefs')

        res = go.Get()
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            #print(objref.Face())
            #print(objref.Surface())
            #print(objref.Trim())

            rgF = objref.Face()
            if isFaceSupported(rgF, bEcho=True):
                return (
                    objref,
                    Opts.values['bEcho'],
                    Opts.values['bDebug'],
                    )

            sc.doc.Objects.UnselectAll()
            sc.doc.Views.Redraw()
            break # out of this while loop.

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Ref(objref_SrfToMod):
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve or edge with which to match")

    go.GeometryFilter = rd.ObjectType.Curve


    def customGeometryFilter(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index)

        if rdObj.Id == objref_SrfToMod.ObjectId:
            return False

        if isinstance(geom, rg.BrepEdge):
            edge = geom
            iFs = edge.AdjacentFaces()
            if len(iFs) != 1:
                return False
            brep = edge.Brep
            srf = brep.Faces[iFs[0]].UnderlyingSurface()
            if isinstance(srf, rg.PlaneSurface):
                return True

            if isinstance(srf, rg.RevSurface):
                return False
            if srf.IsClosed(0) or srf.IsClosed(1):
                return False
            if isinstance(srf, rg.NurbsSurface) and srf.IsRational:
                return False
            if containsShortOrCollapsedBorders(srf):
                return False
            return True
        elif isinstance(geom, rg.Curve):
            crv = geom
            if crv.IsClosed:
                return False
            if isinstance(geom, rg.PolyCurve):
                return True
            if crv.IsRational:
                return False
            return True
        else:
            print("What happened?")
            return False

    go.SetCustomGeometryFilter(customGeometryFilter)


    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iContinuity')
        addOption('iPreserveOtherEnd')
        addOption('bUseUnderlyingIsoCrvs')
        addOption('bMaintainDegree')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')
        if Opts.values['bDebug']:
            addOption('bAddRefs')

        res = go.Get()
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return (
                objref,
                Opts.values['iContinuity'],
                Opts.values['iPreserveOtherEnd'] - 1,
                Opts.values['bUseUnderlyingIsoCrvs'],
                Opts.values['bMaintainDegree'],
                Opts.values['bReplace'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                Opts.values['bAddRefs'],
                )

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def makeNonRational(nurbsGeom, tol=1e-9):
    weights = []
    if isinstance(nurbsGeom, rg.NurbsSurface):
        ns = nurbsGeom
        for iU in range(ns.Points.CountU):
            for iV in range(ns.Points.CountV):
                weights.append(ns.Points.GetWeight(iU, iV))
        if (
            (abs(1.0 - max(weights)) <= tol) and
            (abs(1.0 - min(weights)) <= tol)
        ):
            ns.MakeNonRational()
            return True
        return False
    elif isinstance(nurbsGeom, rg.NurbsCurve):
        nc = nurbsGeom
        for i in range(nc.Points.Count):
            weights.append(nc.Points.GetWeight(i))
        if abs(1.0 - max(weights)) <= tol:
            for i in range(nc.Points.Count):
                if not ns.Points[i].SetWeight(i, 1.0):
                    return False
            return True
        return False
    else:
        return


def areParamsAlignedPerPickPts(objref_A, objref_B):
    """
    For ObjRefs with non-planar Surfaces,
        use U parameterization direction per South and North picks, V for West and East.
    Otherwise, use Curve's parameterization direction.
    """

    def getSidePoints_AnySrfType(srf, side):
        """ Returns tuple(Point3d at start of NS side, Point3d at end) """
        pt_SW = srf.PointAt(srf.Domain(0).T0, srf.Domain(1).T0)
        pt_SE = srf.PointAt(srf.Domain(0).T1, srf.Domain(1).T0)
        pt_NW = srf.PointAt(srf.Domain(0).T0, srf.Domain(1).T1)
        pt_NE = srf.PointAt(srf.Domain(0).T1, srf.Domain(1).T1)
        if side == E:
            return pt_SE, pt_NE
        elif side == W:
            return pt_SW, pt_NW
        elif side == S:
            return pt_SW, pt_SE
        elif side == N:
            return pt_NW, pt_NE
        else:
            raise ValueError("{} is not a valid IsoStatus.".format(side))


    def getSidePoints_NurbsSrf(ns, side):
        """ Returns tuple(Point3d at start of NS side, Point3d at end) """
        cps = ns.Points
        pt_SW = cps.GetControlPoint(           0,            0).Location
        pt_SE = cps.GetControlPoint(cps.CountU-1,            0).Location
        pt_NW = cps.GetControlPoint(           0, cps.CountV-1).Location
        pt_NE = cps.GetControlPoint(cps.CountU-1, cps.CountV-1).Location
        if side == E:
            return pt_SE, pt_NE
        elif side == W:
            return pt_SW, pt_NW
        elif side == S:
            return pt_SW, pt_SE
        elif side == N:
            return pt_NW, pt_NE
        else:
            raise ValueError("{} IsoStatus should have been caught earlier.".format(side))


    def isPickCloserToT0(objref, bUseFullIsoCrv=False):
        pt_Sel = objref.SelectionPoint()
        rdObj = objref.Object()

        if isinstance(rdObj, rd.CurveObject):
            crv = objref.Curve()
            return crv.PointAtStart.DistanceTo(pt_Sel) < crv.PointAtEnd.DistanceTo(pt_Sel)

        trim = objref.Trim()
        if trim is None:
            edge = objref.Edge()
            if edge is None:
                printObjRefObjs(objref)
                raise ValueError("{} not valid.".format(objref.Geometry()))
            iTs = edge.TrimIndices()
            if len(iTs) != 1:
                raise ValueError("Edge has {} trims.".format(len(iTs)))
            trim = edge.Brep.Trims[iTs[0]]

        if trim.IsoStatus == rg.IsoStatus.None:
            crv = objref.Curve()
            return crv.PointAtStart.DistanceTo(pt_Sel) < crv.PointAtEnd.DistanceTo(pt_Sel)

        # Trim is W, S, E, N, X, or Y.

        crv2D = trim.DuplicateCurve()

        if bUseFullIsoCrv:
            crv3D = getIsoCurveOfSide(trim.IsoStatus, trim.Face.UnderlyingSurface())
        else:
            crv3D = trim.Edge.DuplicateCurve()

            if trim.IsReversed():
                crv3D.Reverse()

            if trim.IsoStatus in (S, E):
                pass
            elif trim.IsoStatus in (W, N):
                crv2D.Reverse()
                crv3D.Reverse()
            elif trim.IsoStatus == rg.IsoStatus.X:
                if trim.PointAt(trim.Domain.T1).Y < trim.PointAt(trim.Domain.T0).Y:
                    crv2D.Reverse()
                    crv3D.Reverse()
            elif trim.IsoStatus == rg.IsoStatus.Y:
                if trim.PointAt(trim.Domain.T1).X < trim.PointAt(trim.Domain.T0).X:
                    crv2D.Reverse()
                    crv3D.Reverse()

        #sc.doc.Objects.AddCurve(crv3D); sc.doc.Views.Redraw()#; 1/0

        return crv3D.PointAtStart.DistanceTo(pt_Sel) < crv3D.PointAtEnd.DistanceTo(pt_Sel)

        #pts_Side = getSidePoints_AnySrfType(
        #    trim.Face.UnderlyingSurface(),
        #    trim.IsoStatus)
        #return pts_Side[0].DistanceTo(pt_Sel) < pts_Side[1].DistanceTo(pt_Sel)

    if None in (objref_A.SelectionPoint(), objref_B.SelectionPoint()): return

    bPickedAtSideStart_A = isPickCloserToT0(objref_A, bUseFullIsoCrv=True)
    if bPickedAtSideStart_A is None: return

    bPickedAtSideStart_R = isPickCloserToT0(objref_B)
    if bPickedAtSideStart_R is None: return

    return bPickedAtSideStart_A == bPickedAtSideStart_R


def createTanSrfFromEdge(rgTrim, bDebug=False):
    """
    """

    rgTrim.Brep.Faces.ShrinkFaces()
    rgEdge = rgTrim.Edge
    ns = rgTrim.Face.UnderlyingSurface()

    #ncA_Start = rgEdge.ToNurbsCurve()
    ncA_Start = ns.Pushup(rgTrim, tolerance=0.5*sc.doc.ModelAbsoluteTolerance)

    if bDebug: print(ncA_Start.Points.Count)

    def simplifyCrv(ncA_Start):
        # Try to make Bezier.
        for d in 1, 2, 3, 5:
            nc_WIP = ncA_Start.Rebuild(pointCount=d+1, degree=d, preserveTangents=True)
            rc = rg.Curve.GetDistancesBetweenCurves(
                nc_WIP, ncA_Start, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
            if not rc[0]: continue
            if rc[1] > 0.5*sc.doc.ModelAbsoluteTolerance: continue
            return nc_WIP
        if bDebug: print(ncA_Start.Knots.KnotStyle)

        if ncA_Start.Knots.KnotStyle == rg.KnotStyle.QuasiUniform:
            return

        # QuasiUniform ~ Open and uniform
        # Try to at least make uniform.
        for p in range(3, ncA_Start.Points.Count+1):
            for d in 2, 3, 5:
                if (p - d) < 1: continue
                nc_WIP = ncA_Start.Rebuild(pointCount=p, degree=d, preserveTangents=True)
                rc = rg.Curve.GetDistancesBetweenCurves(
                    nc_WIP, ncA_Start, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if not rc[0]: continue
                dev = rc[1]
                if dev > 0.5*sc.doc.ModelAbsoluteTolerance: continue
                return nc_WIP

    if ncA_Start.SpanCount > 1:
        rc = simplifyCrv(ncA_Start)
        if rc is not None:
            if bDebug:
                print("Simplified starting curve for tangent surface from D{}P{} to D{}P{}.".format(
                    ncA_Start.Degree, ncA_Start.Points.Count,
                    rc.Degree, rc.Points.Count))
            ncA_Start = rc


    ncA = ncA_Start.DuplicateCurve()

    # Do this in order to have a Greville point other than only at the ends.
    if ncA.Degree == 1:
        bDegreeWasIncrFrom1To2 = ncA.IncreaseDegree(2)
    else:
        bDegreeWasIncrFrom1To2 = False

    ts = ncA.GrevilleParameters()

    pts_on_ncA = ncA.GrevillePoints()
    pts_on_ncB = []

    angle = Rhino.RhinoMath.ToRadians(90.0)

    dist = 10.0 * sc.doc.ModelAbsoluteTolerance


    def createEndPt(i, angle):
        t = ts[i]
        #pt_Start = ncA.PointAt(t)
        pt_Start = pts_on_ncA[i]

        bSuccess, u, v = ns.ClosestPoint(pt_Start)
        if not bSuccess:
            print("ClosestPoint failed." \
                "  G1+ will not occur on an edge.")
            return

        vNormal = ns.NormalAt(u, v)
        bSuccess, frame = ncA.PerpendicularFrameAt(t)
        if not bSuccess:
            print("PerpendicularFrameAt failed." \
                "  G1+ will not occur on an edge.")
            return
        
        pt_Normal = pt_Start + vNormal * dist

        pt_End = rg.Point3d(pt_Normal)

        xform_Rotation = rg.Transform.Rotation(
                angleRadians=angle,
                rotationAxis=frame.ZAxis,
                rotationCenter=frame.Origin)

        pt_End.Transform(xform_Rotation)

        return pt_End

    # Create points that are over the Face (-angel vs. angle).
    pt_Mid_Neg = createEndPt(len(ts)//2, -angle)
    pt_Mid_Pos = createEndPt(len(ts)//2,  angle)

    rc = ns.ClosestPoint(pt_Mid_Neg)
    if not rc[0]:
        raise Exception("ClosestPoint failed.")
    pfrel_Neg = rgTrim.Face.IsPointOnFace(*rc[1:], tolerance=0.1*sc.doc.ModelAbsoluteTolerance)

    rc = ns.ClosestPoint(pt_Mid_Pos)
    if not rc[0]:
        raise Exception("ClosestPoint failed.")
    pfrel_Pos = rgTrim.Face.IsPointOnFace(*rc[1:], tolerance=0.1*sc.doc.ModelAbsoluteTolerance)

    if pfrel_Neg == rg.PointFaceRelation.Interior and pfrel_Pos != rg.PointFaceRelation.Interior:
        pts_on_ncB = [createEndPt(i, -angle) for i in range(len(ts))]
    elif pfrel_Neg != rg.PointFaceRelation.Interior and pfrel_Pos == rg.PointFaceRelation.Interior:
        pts_on_ncB = [createEndPt(i,  angle) for i in range(len(ts))]
    else:
        map(sc.doc.Objects.AddPoint, (pt_Mid_Neg, pt_Mid_Pos))
        sc.doc.Objects.AddBrep(rgTrim.Face.DuplicateFace(duplicateMeshes=False))
        raise Exception(
            "PointFaceRelations are {} and {}".format(pfrel_Neg, pfrel_Pos))

    if bDegreeWasIncrFrom1To2:
        ncA = ncA_Start
        ncB = ncA.DuplicateCurve()
        ncB.SetGrevillePoints([pts_on_ncB[0], pts_on_ncB[2]])
    else:
        ncB = ncA.DuplicateCurve()
        ncB.SetGrevillePoints(pts_on_ncB)

    rgB_Loft = rg.Brep.CreateFromLoft(
            curves=[ncA, ncB],
            start=rg.Point3d.Unset,
            end=rg.Point3d.Unset,
            loftType=rg.LoftType.Straight,
            closed=False)

    if len(rgB_Loft) != 1:
        print("{} breps resulted from one curve." \
            "  Continuity increase will not occur on an edge.")

    ns_Out = rgB_Loft[0].Surfaces[0]

    return ns_Out


def getIsoCurveOfSide(isoStatus, srf):

    if isoStatus == W:
        return srf.IsoCurve(direction=1, constantParameter=srf.Domain(0).T0)
    elif isoStatus == S:
        return srf.IsoCurve(direction=0, constantParameter=srf.Domain(1).T0)
    elif isoStatus == E:
        return srf.IsoCurve(direction=1, constantParameter=srf.Domain(0).T1)
    elif isoStatus == N:
        return srf.IsoCurve(direction=0, constantParameter=srf.Domain(1).T1)
    raise ValueError(
        "Isostatus {} not allowed.  Needs to be West, South, East, or North".format(
            isoStatus))


def findMatchingCurveByEndPoints(curvesA, curvesB, bEcho=True):
    """ Returns list(int(Indices of B per order of A)) """

    idxBs_per_As = []

    fTol = 2.0 * sc.doc.ModelAbsoluteTolerance

    for iA, cA in enumerate(curvesA):

        ptsA = [cA.PointAtStart, cA.PointAtEnd]

        for iB, cB in enumerate(curvesB):
            if iB in idxBs_per_As:
                continue

            ptsB = [cB.PointAtStart, cB.PointAtEnd]

            if (
                cA.PointAtStart.DistanceTo(cB.PointAtStart) <= fTol
                and
                cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= fTol
            ):
                idxBs_per_As.append(iB)
                break # to next curveA.

            if (
                cA.PointAtStart.DistanceTo(cB.PointAtEnd) <= fTol
                and
                cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= fTol
            ):
                idxBs_per_As.append(iB)
                break # to next curveA.
        else:
            if len(curvesA) == len(curvesB):
                if bEcho:
                    print("Matching curve not found.")
                return
            idxBs_per_As.append(None)

    return idxBs_per_As


def isFaceOnT1SideOfXYIsoCrv(rgTrim):
    """
    Return:
        True if face is on surface domain's T1 of the X or Y.
        False if face is on T0 side.
    """

    if rgTrim.IsoStatus not in (rg.IsoStatus.X, rg.IsoStatus.Y):
        return

    max_edge_tol = max(edge.Tolerance for edge in rgTrim.Brep.Edges)
    tol = max((max_edge_tol, sc.doc.ModelAbsoluteTolerance))

    b, frame = rgTrim.Edge.PerpendicularFrameAt(rgTrim.Edge.Domain.Mid)
    circle = rg.Circle(plane=frame, radius=10.0*tol)
    arcCrv = rg.ArcCurve(circle=circle)
    #sc.doc.Objects.AddCurve(arcCrv)

    b, crvs, pts = rg.Intersect.Intersection.CurveBrepFace(
        arcCrv, rgTrim.Face, tolerance=0.1*tol)
    if len(pts) == 0:
        return

    #for pt in pts:
    #    sc.doc.Objects.AddPoint(pt)

    if len(pts) > 1:
        for i in range(len(pts)-1):
            ptA = pts[i]
            for j in range(i+1, len(pts)):
                ptB = pts[j]
                if ptA.DistanceTo(ptB) > max_edge_tol:
                    print("Points non-coincident within edge tolerance found.","Check results.")

    pt = pts[0]
    pt_Mid = rgTrim.Edge.PointAt(rgTrim.Edge.Domain.Mid)
    b, uOnCrv, vOnCrv = rgTrim.Face.UnderlyingSurface().ClosestPoint(pt_Mid)

    b, uOnCircle, vOnCircle = rgTrim.Face.UnderlyingSurface().ClosestPoint(pt)

    if rgTrim.IsoStatus == rg.IsoStatus.X:
        # X is 'parallel' with West and East.
        return uOnCircle > uOnCrv
    else:
        return vOnCircle > vOnCrv


def getSrfDomainsWithinIsoCrvEdgeEnds(rgTrim):
    """
    Returns (IsoStatus, NurbsSurface) or None
    """

    edge = rgTrim.Edge
    srf = rgTrim.Face.UnderlyingSurface()

    tol0 = Rhino.RhinoMath.ZeroTolerance

    if rgTrim.IsoStatus in (W, E, rg.IsoStatus.X):
        # X is with West and East because u parameter is constant.
        b, u, v = srf.ClosestPoint(edge.PointAtStart)
        t_Start = v
        if abs(t_Start-srf.Domain(1).T0) <= tol0:
            t_Start = srf.Domain(1).T0 # For Trim if it occurs.
            b, u, v = srf.ClosestPoint(edge.PointAtEnd)
            t_End = v
            if abs(t_End-srf.Domain(1).T1) <= tol0:
                srf.Dispose()
                return
        else:
            b, u, v = srf.ClosestPoint(edge.PointAtEnd)
            t_End = v
            if abs(t_End-srf.Domain(1).T1) <= tol0:
                t_End = srf.Domain(1).T1 # For Trim.
        return srf.Domain(0), rg.Interval(t_Start, t_End)
        trimmed = srf.Trim(srf.Domain(0), rg.Interval(t_Start, t_End))

    if rgTrim.IsoStatus in (S, N, rg.IsoStatus.Y):
        # Y is with South North because v parameter is constant.
        b, u, v = srf.ClosestPoint(edge.PointAtStart)
        t_Start = u
        if abs(t_Start-srf.Domain(0).T0) <= tol0:
            t_Start = srf.Domain(0).T0 # For Trim if it occurs.
            b, u, v = srf.ClosestPoint(edge.PointAtEnd)
            t_End = u
            if abs(t_End-srf.Domain(0).T1) <= tol0:
                srf.Dispose()
                return
        else:
            b, u, v = srf.ClosestPoint(edge.PointAtEnd)
            t_End = u
            if abs(t_End-srf.Domain(0).T1) <= tol0:
                t_End = srf.Domain(0).T1 # For Trim.
        return rg.Interval(t_Start, t_End), srf.Domain(1)
        trimmed = srf.Trim(rg.Interval(t_Start, t_End), srf.Domain(1))

    return rgTrim.IsoStatus, trimmed


def getShrunkNurbsSrfFromGeom(geom_R_In, bUseUnderlyingIsoCrvs):
    """
    """
    if not isinstance(geom_R_In, rg.BrepTrim): return

    trim_In = geom_R_In

    if trim_In.IsoStatus == rg.IsoStatus.None: return


    domU, domV = None, None

    if not bUseUnderlyingIsoCrvs:
        rc = getSrfDomainsWithinIsoCrvEdgeEnds(trim_In)
        if rc is not None:
            domU, domV = rc

            #_ = trim_In.Face.UnderlyingSurface().Trim(domU, domV)
            #sc.doc.Objects.AddSurface(_); sc.doc.Views.Redraw(); 1/0

    #sc.doc.Objects.AddCurve(trim_In)

    iso_Out = trim_In.IsoStatus

    if trim_In.IsoStatus == rg.IsoStatus.X:
        rc = isFaceOnT1SideOfXYIsoCrv(trim_In)
        if rc is not None:
            if rc:
                domU = rg.Interval(
                    trim_In.PointAt(trim_In.Domain.Mid).X,
                    trim_In.Face.UnderlyingSurface().Domain(0).T1,
                    )
                iso_Out = W
            else:
                domU = rg.Interval(
                    trim_In.Face.UnderlyingSurface().Domain(0).T0,
                    trim_In.PointAt(trim_In.Domain.Mid).X,
                    )
                iso_Out = E
    if trim_In.IsoStatus == rg.IsoStatus.Y:
        rc = isFaceOnT1SideOfXYIsoCrv(trim_In)
        if rc is not None:
            if rc:
                domV = rg.Interval(
                    trim_In.PointAt(trim_In.Domain.Mid).Y,
                    trim_In.Face.UnderlyingSurface().Domain(1).T1,
                    )
                iso_Out = S
            else:
                domV = rg.Interval(
                    trim_In.Face.UnderlyingSurface().Domain(1).T0,
                    trim_In.PointAt(trim_In.Domain.Mid).Y,
                    )
                iso_Out = N

    if domU is None and domV is None:
        # No trimming.
        return

    if domU is None:
        domU = trim_In.Face.UnderlyingSurface().Domain(0)
    elif domV is None:
        domV = trim_In.Face.UnderlyingSurface().Domain(1)

    trimmed = trim_In.Face.UnderlyingSurface().Trim(domU, domV)

    return iso_Out, trimmed.ToNurbsSurface()


def getNurbsGeomFromGeom(In, iContinuity, bEcho=True, bDebug=False):
    """
    Returns 4-item list containing any combination of tuples of:
        (IsoStatus, NurbsSurface),
        (NurbsCurve, None)
    """

    if not isinstance(In, rg.GeometryBase):
        raise ValueError("{} not accepted.".format(In.GetType().Name))

    if not isinstance(In, rg.BrepTrim):
        return In.ToNurbsCurve(), None

    trim = In
    srf = trim.Face.UnderlyingSurface()

    if trim.IsoStatus == rg.IsoStatus.None:
        ns_Tan = createTanSrfFromEdge(trim, bDebug)
        if ns_Tan is None:
            return In.ToNurbsCurve(), None
        ncs_TanNS = [getIsoCurveOfSide(side, ns_Tan) for side in (W,S,E,N)]
        idxB = findMatchingCurveByEndPoints([trim.Edge], ncs_TanNS, bEcho=bEcho)[0]
        if bEcho:
            if not srf.IsPlanar(1e-9):
                s = "Non-isocurve trim of a non-planar surface picked."
                if iContinuity == 2:
                    s += "  G2 will probably not be achieved"
                if iContinuity > 0:
                    s += "  G1 may only be approximately achieved."
                    s += "  Use _EdgeContinuity to check results."
                print(s)
        return (W,S,E,N)[idxB], ns_Tan

    return trim.IsoStatus, srf.ToNurbsSurface()


def createAlignedRefGeom(geom_R_In, ns_M, side_M, bMatchWithParamsAligned):
    """
    Returns either
        NurbsSurface with parameterizations aligned with ns_M along side_M or
        NurbsCurve with parameterization aligned with ns_M
    """


    def createAlignedRefSrfs(side_M, ns_R, side_R):
        """
        Notes:
        ns_R.Reverse(int direction, bool inPlace)
        Transpose(bool inPlace)
        """

        if side_M == side_R:
            pass
        elif side_M in (S, N) and side_R in (S, N):
            ns_R.Reverse(1, True)
        elif side_M in (W, E) and side_R in (W, E):
            ns_R.Reverse(0, True)
        else:
            ns_R.Transpose(True)

            if side_M == S and side_R == E:
                ns_R.Reverse(1, True)
            elif side_M == S and side_R == W:
                pass
            elif side_M == E and side_R == N:
                pass
            elif side_M == E and side_R == S:
                ns_R.Reverse(0, True)
            elif side_M == N and side_R == W:
                ns_R.Reverse(1, True)
            elif side_M == N and side_R == E:
                pass
            elif side_M == W and side_R == S:
                pass
            elif side_M == W and side_R == N:
                ns_R.Reverse(0, True)
            else:
                raise ValueError("What happened?")

        if not bMatchWithParamsAligned:
            ns_R.Reverse(0 if side_M in (S, N) else 1, True)


    if isinstance(geom_R_In[1], rg.NurbsSurface):
        side_R, ns_R = geom_R_In
        rc = createAlignedRefSrfs(side_M, ns_R, side_R)
        # Not returning IsoStatus since R M-facing side now matches that of M.
        return ns_R

    # Reference alignment is per Edge.
    c_M = getIsoCurveOfSide(side_M, ns_M)
    c_R = geom_R_In[0]

    #map(sc.doc.Objects.AddCurve, (c_A, c_R))

    if not bMatchWithParamsAligned:
        c_R.Reverse()

     # Don't bother returning tuple since none is needed for NS.
    return c_R


def transferHigherDegree(nFrom, nTo, sideFrom, sideTo):
    if isinstance(nFrom, rg.NurbsSurface):
        degreeFrom = nFrom.Degree(0 if sideFrom in (S, N) else 1)
    else:
        degreeFrom = nFrom.Degree

    if isinstance(nTo, rg.NurbsSurface):
        if sideTo in (S,N):
            if nTo.Degree(0) >= degreeFrom: return False
            return nTo.IncreaseDegreeU(degreeFrom)
        else:
            if nTo.Degree(1) >= degreeFrom: return False
            return nTo.IncreaseDegreeV(degreeFrom)
    else:
        if nTo.Degree >= degreeFrom: return False
        return nTo.IncreaseDegree(degreeFrom)


def transferDomain(nFrom, nTo, sideFrom, sideTo):
    tol0 = Rhino.RhinoMath.ZeroTolerance
    if isinstance(nFrom, rg.NurbsSurface):
        domainFrom = nFrom.Domain(0 if sideFrom in (S, N) else 1)
    else:
        domainFrom = nFrom.Domain

    if isinstance(nTo, rg.NurbsSurface):
        if sideTo in (S,N):
            return nTo.SetDomain(0, domainFrom)
        else:
            return nTo.SetDomain(1, domainFrom)
    else:
        nTo.Domain = domainFrom
        return (nTo.Domain.EpsilonEquals(domainFrom, epsilon=tol0))


def transferUniqueKnotVector(nurbsA, nurbsB, sideA=None, sideB=None, paramTol=None):
    """
    Union the knot vectors of nurbsA and nurbsB and apply to nurbsB (modify reference).

    nurbsA and nurbsB: Any combination of NurbsSurface or NurbsCurve.
    Return True if nsB was modified.
    """

    if paramTol is None: paramTol = Rhino.RhinoMath.ZeroTolerance


    def knotCount(geom, side):
        if isinstance(geom, rg.NurbsSurface):
            if side in (S, N): return geom.KnotsU.Count
            else: return geom.KnotsV.Count
        else:
            return geom.Knots.Count

    def degree(geom, side):
        if isinstance(geom, rg.NurbsSurface):
            iDir = 0 if side in (S, N) else 1
            return geom.Degree(iDir)
        else:
            return geom.Degree

    def knotMultiplicity(geom, side, iK):
        if isinstance(geom, rg.NurbsSurface):
            if side in (S, N): return geom.KnotsU.KnotMultiplicity(iK)
            else: return geom.KnotsV.KnotMultiplicity(iK)
        else:
            return geom.Knots.KnotMultiplicity(iK)

    def getKnot(geom, side, iK):
        if isinstance(geom, rg.NurbsSurface):
            if side in (S, N): return geom.KnotsU[iK]
            else: return geom.KnotsV[iK]
        else:
            return geom.Knots[iK]

    def insertKnot(geom, side, t, m):
        if isinstance(geom, rg.NurbsSurface):
            if side in (S, N): return geom.KnotsU.InsertKnot(t, m)
            else: return geom.KnotsV.InsertKnot(t, m)
        else:
            return geom.Knots.InsertKnot(t, m)


    knotCt_In = knotCount(nurbsB, sideB)

    iK = degree(nurbsA, sideA)
    while iK < knotCount(nurbsA, sideA)-degree(nurbsA, sideA):
        sc.escape_test()
        t_A = getKnot(nurbsA, sideA, iK)
        m = knotMultiplicity(nurbsA, sideA, iK)
        if abs(t_A - getKnot(nurbsB, sideB, iK)) > paramTol:
            insertKnot(nurbsB, sideB, t_A, m)
        iK += m

    return knotCount(nurbsB, sideB) > knotCt_In


def getPtRowIndicesPerG(nurbsGeom, side, iG):
    if isinstance(nurbsGeom, rg.NurbsCurve):
        if iG == 0:
            return range(nurbsGeom.Points.Count)
        return

    # NurbsSurface.
    cps = list(nurbsGeom.Points)
    ctV = nurbsGeom.Points.CountV
    ctAll = len(cps)
    if side == W:
        return range(iG*ctV, (iG+1)*ctV)
    elif side == E:
        return range(ctAll-(iG+1)*ctV, ctAll-iG*ctV)
    elif side == S:
        return range(iG, ctAll, ctV)
    elif side == N:
        return range(ctV-(iG+1), ctAll, ctV)


def getUvIdxFromNsPoint1dIdx(ns, idxFlat):
    return idxFlat // ns.Points.CountV, idxFlat % ns.Points.CountV


def project_C_Colinear_with_AB(ptA, ptB, ptC):
    return rg.Line(ptA, ptB).ClosestPoint(ptC, limitToFiniteSegment=False)


def setContinuity_G0(**kwargs):
    """
    Parameters:
        ns_M_In: NurbsSurface
        side_M: IsoStatus
        nurbs_R: NurbsSurface or NurbsCurve
        side_R: IsoStatus
    Returns
    """

    def get_kwarg(key):
        if key not in kwargs:
            raise ("{} missing from provided parameters.".format(key))
        return kwargs[key] if key in kwargs else None

    ns_M_In = get_kwarg('ns_M_In')
    side_M = get_kwarg('side_M')
    nurbs_R = get_kwarg('nurbs_R')
    side_R = get_kwarg('side_R')


    bValid, sLog = ns_M_In.IsValidWithLog()
    if not bValid: print(sLog); return


    ns_M_Out = ns_M_In.Duplicate()
    ptsR = [cp.Location for cp in nurbs_R.Points]

    # Points are ordered (u0,v0), (u0,v1), ..., (un,vn).


    idxPts = {} # Key is tuple(str('M' or 'R'), int(G continuity))
    idxPts['M',0] = getPtRowIndicesPerG(ns_M_In, side_M, 0)
    idxPts['R',0] = getPtRowIndicesPerG(nurbs_R, side_R, 0)


    for i in range(len(idxPts['M',0])):
        iU_G0_A, iV_G0_A = getUvIdxFromNsPoint1dIdx(ns_M_Out, idxPts['M',0][i])
        iR0 = idxPts['R',0][i]
        if isinstance(nurbs_R, rg.NurbsSurface):
            iU_0_R, iV_0_R = getUvIdxFromNsPoint1dIdx(nurbs_R, iR0)
            ns_M_Out.Points.SetPoint(
                iU_G0_A, iV_G0_A,
                point=ptsR[iR0],
                weight=nurbs_R.Points.GetWeight(iU_0_R, iV_0_R))
        else: # NurbsCurve.
            ns_M_Out.Points.SetPoint(
                iU_G0_A, iV_G0_A,
                point=ptsR[iR0],
                weight=nurbs_R.Points.GetWeight(iR0))

    return ns_M_Out


def setContinuity_G1(**kwargs):
    """
    Parameters:
        ns_M_BeforeAnyMatching: NurbsSurface
        ns_M_In: NurbsSurface
        side_M: IsoStatus
        nurbs_R: NurbsSurface or NurbsCurve (for planar surface)
        side_R: IsoStatus
        bModifyRowEnd_T0: bool
        bModifyRowEnd_T1: bool
    Returns
    """

    def get_kwarg(key):
        if key not in kwargs:
            raise Exception("{} missing from provided parameters.".format(key))
        return kwargs[key] if key in kwargs else None

    ns_M_BeforeAnyMatching = get_kwarg('ns_M_BeforeAnyMatching')
    ns_M_In = get_kwarg('ns_M_In')
    side_M = get_kwarg('side_M')
    nurbs_R = get_kwarg('nurbs_R')
    side_R = get_kwarg('side_R')
    bModifyRowEnd_T0 = get_kwarg('bModifyRowEnd_T0')
    bModifyRowEnd_T1 = get_kwarg('bModifyRowEnd_T1')


    if isinstance(nurbs_R, rg.NurbsCurve): return

    ns_R = nurbs_R

    bValid, sLog = ns_M_In.IsValidWithLog()
    if not bValid: print(sLog); return


    ns_M_Out = ns_M_In.Duplicate()
    ptsM_BeforeAnyMatching = [cp.Location for cp in ns_M_BeforeAnyMatching.Points]
    ptsM_Out = [cp.Location for cp in ns_M_In.Points]
    ptsR = [cp.Location for cp in ns_R.Points]

    # Points are ordered (u0,v0), (u0,v1), ..., (un,vn).


    idxPts = {} # Key is tuple(str('M' or 'R'), int(G continuity))
    for i in 0,1:
        idxPts['M',i] = getPtRowIndicesPerG(ns_M_In, side_M, iG=i)
        idxPts['R',i] = getPtRowIndicesPerG(ns_R, side_R, iG=i)

    iRowLn = len(idxPts['M',0])


    # Set G1 row colinear with reference.

    range_start = 0 if bModifyRowEnd_T0 else 1
    range_stop = iRowLn-1+1 if bModifyRowEnd_T1 else iRowLn-2+1
    Range = range(range_start, range_stop)

    for i in Range:
        iM0 = idxPts['M',0][i]
        iM1 = idxPts['M',1][i]
        iR1 = idxPts['R',1][i]
        pt_To = project_C_Colinear_with_AB(ptsM_Out[iM0], ptsR[iR1], ptsM_Out[iM1])
        iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_M_Out, iM1)
        ns_M_Out.Points.SetPoint(
            iUT_A, iVT_A, pt_To,
            weight=ns_M_BeforeAnyMatching.Points.GetWeight(iUT_A, iVT_A))


    # Update.
    ptsM_Out = [cp.Location for cp in ns_M_Out.Points]


    # To obtain G1 continuity along edge, set the tangential row at a scale
    # of the G0-G1 point distances of the reference surface.

    # Determine the average ratio of the G0-G1 point distances of surface
    # to match to the same of reference surface.
    fDistRatios = []

    if bModifyRowEnd_T0 != bModifyRowEnd_T1:
        sample = (0,) if not bModifyRowEnd_T0 else (iRowLn-1,)
    else:
        sample = 0, iRowLn-1

    for i in sample:
        iM0 = idxPts['M',0][i]
        iM1 = idxPts['M',1][i]
        iR1 = idxPts['R',1][i]
        fDist_M = ptsM_BeforeAnyMatching[iM1].DistanceTo(ptsM_BeforeAnyMatching[iM0])
        fDist_R = ptsR[iR1].DistanceTo(ptsM_Out[iM0])
        fRatio = fDist_M / fDist_R
        fDistRatios.append(fRatio)
    fAvgDistRatio = sum(fDistRatios) / float(len(fDistRatios))

    # Set tangential row using the multiple of the average G0-G1 point distance(s) ratio
    for i in Range:
        iM0 = idxPts['M',0][i]
        iM1 = idxPts['M',1][i]
        iR1 = idxPts['R',1][i]
        vG0G1_R = ptsM_Out[iM0] - ptsR[iR1]
        vG0G1_A_toSet = fAvgDistRatio * vG0G1_R
        pt_To = ptsM_Out[iM0] + vG0G1_A_toSet
        iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_M_Out, iM1)
        ns_M_Out.Points.SetPoint(
            iUT_A, iVT_A, pt_To,
            weight=ns_M_BeforeAnyMatching.Points.GetWeight(iUT_A, iVT_A))

    return ns_M_Out


def setContinuity_G2(**kwargs):
    """
    Parameters:
        ns_M_BeforeAnyMatching: NurbsSurface
        ns_M_In: NurbsSurface
        side_M: IsoStatus
        nurbs_R: NurbsSurface or NurbsCurve
        side_R: IsoStatus
        bModifyRowEnd_T0: bool
        bModifyRowEnd_T1: bool
        bDebug: bool
        bAddRefs: bool
    Returns
    """

    def get_kwarg(key):
        if key not in kwargs:
            raise ValueError("{} missing from provided parameters.".format(key))
        return kwargs[key] if key in kwargs else None

    ns_M_BeforeAnyMatching = get_kwarg('ns_M_BeforeAnyMatching')
    ns_M_In = get_kwarg('ns_M_In')
    side_M = get_kwarg('side_M')
    nurbs_R = get_kwarg('nurbs_R')
    side_R = get_kwarg('side_R')
    bModifyRowEnd_T0 = get_kwarg('bModifyRowEnd_T0')
    bModifyRowEnd_T1 = get_kwarg('bModifyRowEnd_T1')
    bDebug = get_kwarg('bDebug')
    bAddRefs = get_kwarg('bAddRefs')

    if isinstance(nurbs_R, rg.NurbsCurve): return

    ns_R = nurbs_R

    bValid, sLog = ns_M_In.IsValidWithLog()
    if not bValid: print(sLog); return


    ns_M_Out = ns_M_In.Duplicate()
    ptsM_BeforeAnyMatching = [cp.Location for cp in ns_M_BeforeAnyMatching.Points]
    ptsA_Out = [cp.Location for cp in ns_M_In.Points]
    ptsR = [cp.Location for cp in ns_R.Points]

    # Points are ordered (u0,v0), (u0,v1), ..., (un,vn).


    idxPts = {} # Key is tuple(str('M' or 'R'), int(G continuity))
    for i in 0,1,2:
        idxPts['M',i] = getPtRowIndicesPerG(ns_M_In, side_M, iG=i)
        idxPts['R',i] = getPtRowIndicesPerG(ns_R, side_R, iG=i)

    iRowLn = len(idxPts['M',0])


    if ns_R.Degree(int(side_R in (S,N))) == 1:
        # Move G2 points of A to be colinear with G0 and G1 points.

        for i in range(iRowLn):
            iA0 = idxPts['M',0][i]
            iA1 = idxPts['M',1][i]
            iA2 = idxPts['M',2][i]
            pt_To = project_C_Colinear_with_AB(ptsA_Out[iA0], ptsA_Out[iA1], ptsA_Out[iA2])
            iUT, iVT = getUvIdxFromNsPoint1dIdx(ns_M_Out, iA2)
            ns_M_Out.Points.SetPoint(
                iUT, iVT, pt_To,
                weight=ns_M_BeforeAnyMatching.Points.GetWeight(iUT, iVT))

        return ns_M_Out


    def multiplierForDegree(ns, side_Match):
        degree = ns.Degree(int(side_Match in (S,N)))
        return float(degree - 1) / float(degree)


    mA_Deg = multiplierForDegree(ns_M_Out, side_M)

    mR_Deg = multiplierForDegree(ns_R, side_R)


    def multiplierForAdjacentSimpleKnots(ns, side_Match):

        iDir_FromMatch = int(side_Match in (S,N))

        iCt_Span = ns.SpanCount(iDir_FromMatch)
        if iCt_Span == 1:
            m_ForKnots = 1.0
        else:
            knots = ns.KnotsV if iDir_FromMatch else ns.KnotsU
            iDeg_FromMSide_A = ns.Degree(iDir_FromMatch)

            if side_Match in (E, N):
                iKnot_NextToEnd = knots.Count - iDeg_FromMSide_A - 1
            else:
                iKnot_NextToEnd = iDeg_FromMSide_A

            multy_AdjKnot = knots.KnotMultiplicity(iKnot_NextToEnd)

            if multy_AdjKnot > 1:
                m_ForKnots = 1.0
            else:
                t_Knots_Span = ns.GetSpanVector(iDir_FromMatch)

                if side_Match in (E, N):
                    fLen_spanEnd = t_Knots_Span[iCt_Span] - t_Knots_Span[iCt_Span-1]
                    fLen_spanAdj = t_Knots_Span[iCt_Span-1] - t_Knots_Span[iCt_Span-2]
                else:
                    fLen_spanEnd = t_Knots_Span[1] - t_Knots_Span[0]
                    fLen_spanAdj = t_Knots_Span[2] - t_Knots_Span[1]

                m_ForKnots = (fLen_spanEnd + fLen_spanAdj) / fLen_spanEnd
        return m_ForKnots

    mA_Knots = multiplierForAdjacentSimpleKnots(ns_M_Out, side_M)
    mR_Knots = multiplierForAdjacentSimpleKnots(ns_R, side_R)

    range_start = 0 if bModifyRowEnd_T0 else 2
    range_stop = iRowLn-1+1 if bModifyRowEnd_T1 else iRowLn-3+1
    Range = range(range_start, range_stop)

    a2s = {} # Collect initial a2 locations to use for smooth projections.

    b_trans_a2_thru_a1 = False

    for i in Range:
        a0 = r0 = ptsR[idxPts['R',0][i]]
        a1 = ptsA_Out[idxPts['M',1][i]]
        r1 = ptsR[idxPts['R',1][i]]
        r2 = ptsR[idxPts['R',2][i]]

        m = ((a1 - r0).Length/(r1 - r0).Length)**2.0

        M = m * (mA_Knots / mR_Knots) * (mR_Deg / mA_Deg)

        a2 = 2.0*a1 + -a0 + M*(-2.0*r1 + r2 + r0)

        a2s[i] = a2

        plane = rg.Plane(origin=a0, normal=a1-a0)
        #addGeoms(plane, False)
        p2 = plane.ClosestPoint(a2)

        if ((a2-p2) * (a1-a0)) < 0.0:
            # a2 is on the wrong side of a1.
            b_trans_a2_thru_a1 = True

        iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_M_Out, idxPts['M',2][i])
        ns_M_Out.Points.SetPoint(
            iUT_A, iVT_A, a2,
            weight=ns_M_BeforeAnyMatching.Points.GetWeight(iUT_A, iVT_A))

        if bDebug and bAddRefs:
            #addGeoms([r0, a1, a2, r1, r2], False)
            addGeoms(a2, False)
            #print(a2)


    if not b_trans_a2_thru_a1:
        return ns_M_Out


    if bModifyRowEnd_T0 and bModifyRowEnd_T1:
        iRefs = range(iRowLn)
    elif not (bModifyRowEnd_T0 and bModifyRowEnd_T1):
        iRefs = 0, iRowLn-1
    elif not bModifyRowEnd_T0:
        iRefs = 0,
    elif not bModifyRowEnd_T1:
        iRefs = iRowLn-1,


    m0m1_ratios_sum = 0.0

    for i in iRefs:
        a0 = ptsM_BeforeAnyMatching[idxPts['M',0][i]]
        a1 = ptsM_BeforeAnyMatching[idxPts['M',1][i]]
        a2 = ptsM_BeforeAnyMatching[idxPts['M',2][i]]
        plane = rg.Plane(origin=a0, normal=a1-a0)
        p2 = plane.ClosestPoint(a2)
        p2m2 = (a2 - p2)
        ratio = p2m2.Length / (a1 - a0).Length
        m0m1_ratios_sum += ratio

    m0m1_avg_ratio = m0m1_ratios_sum / len(iRefs)


    for i in Range:
        a0 = ptsA_Out[idxPts['M',0][i]]
        a1 = ptsA_Out[idxPts['M',1][i]]

        plane = rg.Plane(origin=a0, normal=a1-a0)
        plane.Translate(m0m1_avg_ratio * (a1 - a0))
        a2 = plane.ClosestPoint(a2s[i])


        iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_M_Out, idxPts['M',2][i])
        ns_M_Out.Points.SetPoint(
            iUT_A, iVT_A, a2,
            weight=ns_M_BeforeAnyMatching.Points.GetWeight(iUT_A, iVT_A))


    if bDebug and bAddRefs: addGeoms(ns_R, False)


    return ns_M_Out


def createSurface(ns_M_In, side_M, geom_R_In, bMatchWithParamsAligned=True, **kwargs):
    """
    Parameters:
        ns_M_In: NurbsSurface to modify (duplicate).
        side_M=side_M,
        geom_R: BrepTrim or Curve for wire=geom_R_In,
        bMatchWithParamsAligned: bool
        iContinuity=iContinuity,
        iPreserveOtherEnd=iPreserveOtherEnd,
        bUseUnderlyingIsoCrvs=bUseUnderlyingIsoCrvs,
        bMaintainDegree=bMaintainDegree,
        bDebug=bDebug,

    Returns on success: New rg.NurbsSurface, int(continuity processed, not necessarily fully achieved)
    Returns on fail: None
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iContinuity = getOpt('iContinuity')
    iPreserveOtherEnd = getOpt('iPreserveOtherEnd')
    bUseUnderlyingIsoCrvs = getOpt('bUseUnderlyingIsoCrvs')
    bMaintainDegree = getOpt('bMaintainDegree')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    bAddRefs = getOpt('bAddRefs')


    if iContinuity == 2:
        if isinstance(geom_R_In, rg.BrepTrim):
           if geom_R_In.IsoStatus == rg.IsoStatus.None:
                iContinuity = 1

    bValid, sLog = ns_M_In.IsValidWithLog()
    if not bValid: print(sLog); return


    if side_M not in (W,S,E,N):
        return


    geom_R_Nurbs = getShrunkNurbsSrfFromGeom(geom_R_In, bUseUnderlyingIsoCrvs)

    if geom_R_Nurbs is None:
        geom_R_Nurbs = getNurbsGeomFromGeom(geom_R_In, iContinuity, bEcho, bDebug)
        if not geom_R_Nurbs: return

        if isinstance(geom_R_Nurbs[1], rg.NurbsSurface):
            ns_R = geom_R_Nurbs[1]
            if ns_R.IsRational:
                makeNonRational(ns_R)
                if ns_R.IsRational:
                    print("Warning... Reference surface is rational.  Check results.")
        elif isinstance(geom_R_Nurbs[0], rg.NurbsCurve):
            nc_R = geom_R_Nurbs[0]
            if nc_R.IsRational:
                makeNonRational(nc_R)
                if nc_R.IsRational:
                    print("Warning... Reference curve is rational.  Check results.")



    #if isinstance(geom_R_Nurbs[1], rg.NurbsSurface): sc.doc.Objects.AddSurface(geom_R_Nurbs[1])
    #else: sc.doc.Objects.AddCurve(geom_R_Nurbs[0])
    #sc.doc.Views.Redraw(); 1/0


    rc = createAlignedRefGeom(
        geom_R_Nurbs, ns_M_In, side_M, bMatchWithParamsAligned)
    if rc is None: return
    nurbs_R = rc

    #if isinstance(nurbs_R, rg.NurbsSurface): sc.doc.Objects.AddSurface(nurbs_R)
    #else: sc.doc.Objects.AddCurve(nurbs_R)
    #sc.doc.Views.Redraw(); 1/0


    bValid, sLog = nurbs_R.IsValidWithLog()
    if not bValid: print(sLog); return

    side_R = side_M
    side_Both = side_M

    #sc.doc.Objects.AddSurface(nurbs_R); sc.doc.Views.Redraw(); return
    #sEval = "bMatchWithParamsAligned"; print(sEval+':',eval(sEval))
    #print(side_M, side_R)



    ns_WIP = ns_M_In.Duplicate()

    ns_R_In = nurbs_R.Duplicate()



    # Along matching side of surface, match reference's degree if larger.
    # Increase degree if necessary.

    transferHigherDegree(
            nurbs_R, ns_WIP, side_Both, side_Both)

    transferHigherDegree(
            ns_WIP, nurbs_R, side_Both, side_Both)

    transferDomain(
            ns_WIP, nurbs_R, side_Both, side_Both)





    iDir_MatchSide_R = int(side_Both in (E, W))
    #iDeg_MatchSide_R = ns_R.Degree(iDir_MatchSide_R)

    iDir_MatchSide_A = int(side_Both in (E, W))
    #iDeg_MatchSide_A = ns_WIP.Degree(iDir_MatchSide_A)



    def increaseDegree(ns, iDir, iDeg):
        if iDir:
            return ns.IncreaseDegreeV(iDeg)
        else:
            return ns.IncreaseDegreeU(iDeg)

    def getDirString(iDir):
        return 'V' if iDir else 'U'


    bValid, sLog = ns_WIP.IsValidWithLog()
    if not bValid: print(sLog); return


    # Along match side of surface, match reference's knot vector since both surface degrees are the same.

    # First, transfer unique knots to Coons so that the latter can contain all unqiue knots
    # to transfer to reference surfaces.
    transferUniqueKnotVector(nurbs_R, ns_WIP, side_Both, side_Both)

    #sc.doc.Objects.AddSurface(ns_WIP)#; sc.doc.Views.Redraw(); return


    # Back to reference surfaces.
    transferUniqueKnotVector(ns_WIP, nurbs_R, side_Both, side_Both)

    #sc.doc.Objects.AddSurface(ns_R)
    #sc.doc.Views.Redraw(); return




    def getKnots(ns ,iDir):
        return ns.KnotsV if iDir else ns.KnotsU


    def addKnotsUniformly(coll_Knots, iCt_Knot_Target):
        degree = coll_Knots.KnotMultiplicity(0)

        # Single-span.
        if coll_Knots.Count == 2 * degree:
            tsK_ToAdd = []
            iCt_toAdd = iCt_Knot_Target - coll_Knots.Count
            for iK in range(iCt_Knot_Target - coll_Knots.Count):
                n = float(iK+1) * (coll_Knots[coll_Knots.Count-1] - coll_Knots[0])
                d = float(iCt_toAdd + 1)
                tsK_ToAdd.append(n/d + coll_Knots[0])
            for tK in tsK_ToAdd:
                coll_Knots.InsertKnot(tK)
            return

        # Add knots between those in multiple-spanned vector.
        while True:
            sc.escape_test()
            rangeKnotsToTraverse = range(degree-1, coll_Knots.Count - degree)
            tsK_ToAdd = []
            for iK in rangeKnotsToTraverse:
                tsK_ToAdd.append(0.5*coll_Knots[iK] + 0.5*coll_Knots[iK+1])
            for tK in tsK_ToAdd:
                coll_Knots.InsertKnot(tK)
            if coll_Knots.Count >= iCt_Knot_Target:
                return


    bValid, sLog = ns_WIP.IsValidWithLog()
    if not bValid:
        print(sLog); return
    #sc.doc.Objects.AddSurface(ns_WIP); sc.doc.Views.Redraw(); return


    # From side of surface, increase point count as needed.
    # Points needed for continuity change: iContinuity + 1
    # Points needed to preserve other side: iPreserveOtherEnd
    iDir_FromMSide_A = int(not iDir_MatchSide_A)
    iDeg_FromMSide_A = ns_WIP.Degree(iDir_FromMSide_A)


    def getPtCt(ns, iDir):
        return ns.Points.CountV if iDir else ns.Points.CountU

    iPtCt_FromMSide_A = getPtCt(ns_WIP, iDir_FromMSide_A)


    iPtCt_FromMSide_A_Min = iContinuity + 1 + iPreserveOtherEnd + 1

    if not bMaintainDegree:# and ns_WIP.SpanCount(iDir_FromMSide_A) == 1:
        # Create single-spanned at minimal degree.
        iPtCt_FromMSide_A_ToAdd = iPtCt_FromMSide_A_Min - iPtCt_FromMSide_A
        iDeg_FromSide_A_Min = iPtCt_FromMSide_A_Min - 1
        if iDeg_FromMSide_A < iDeg_FromSide_A_Min:
            if increaseDegree(ns_WIP, iDir_FromMSide_A, iDeg_FromSide_A_Min):
                print("Increased degree of surface in {} direction from {} to {}.".format(
                    getDirString(iDir_FromMSide_A),
                    iDeg_FromMSide_A,
                    iDeg_FromSide_A_Min))

                # Update variables since degree was increased.
                iDeg_FromMSide_A = ns_WIP.Degree(iDir_FromMSide_A)
                iPtCt_FromMSide_A = getPtCt(ns_WIP, iDir_FromMSide_A)

    else:
        # Maintain degree, increasing span count.
        iPtCt_FromMSide_A_ToAdd = iPtCt_FromMSide_A_Min - iPtCt_FromMSide_A

        if iPtCt_FromMSide_A_ToAdd > 0:
            
            iKnotCt_FromSide_A_ToAdd = iPtCt_FromMSide_A_ToAdd

            coll_Knots = getKnots(ns_WIP, iDir_FromMSide_A) # Reference to collection.  Do not convert to list.

            iKnotCt_FromSide_A_Required = iPtCt_FromMSide_A_ToAdd + coll_Knots.Count

            # TODO: Use the following only to maintain knot uniformity.  Otherwise, only add knots as needed.
            tsK_New = []
            iCt_Knot_Start = coll_Knots.Count

            addKnotsUniformly(
                coll_Knots=coll_Knots,
                iCt_Knot_Target=iKnotCt_FromSide_A_Required
                )

            print("Increased knot count from {} to {}.".format(
                    iCt_Knot_Start, coll_Knots.Count))
    
            #print(coll_Knots.CreateUniformKnots(knotSpacing=1.0))
    
            ## Only increase degree to an odd value.
            #if iDegIncr % 2 == 0:
            #    iDegIncr += 1
            #if ns_WIP.IncreaseDegreeV(iDegIncr):
            #    print("Increase degree of surface in V direction from {} to {}.".format(
            #        iDeg_FromMSide_A, iDegIncr))

    ns_M_BeforeAnyMatching = ns_WIP.Duplicate()


    def areInOutEpsilonEqual(In, Out):
        epsilon = 1e-6
        if Out.EpsilonEquals(In, epsilon):
            while True:
                eps_prev = epsilon
                epsilon /= 10.0
                if not Out.EpsilonEquals(In, epsilon):
                    sc.escape_test()
                    break
            if bEcho:
                print("Input and output NURBS surfaces are EpsilonEqual within {}." \
                    "  No change.".format(eps_prev))
            Out.Dispose()
            return True
        return False


    ns_M_G0 = setContinuity_G0(
        ns_M_In=ns_M_BeforeAnyMatching,
        side_M=side_Both,
        nurbs_R=nurbs_R,
        side_R=side_Both,
        )

    if ns_M_G0 is None:
        return None, None

    if iContinuity == 0:
        if bEcho: print("Modified surface toward G0.")
        if areInOutEpsilonEqual(ns_M_BeforeAnyMatching, ns_M_G0): return
        return ns_M_G0, 0

    ns_M_G1 = setContinuity_G1(
        ns_M_BeforeAnyMatching=ns_M_BeforeAnyMatching,
        ns_M_In=ns_M_G0,
        side_M=side_Both,
        nurbs_R=nurbs_R,
        side_R=side_Both,
        bModifyRowEnd_T0=True,
        bModifyRowEnd_T1=True,
        )

    if ns_M_G1 is None:
        if areInOutEpsilonEqual(ns_M_BeforeAnyMatching, ns_M_G0): return
        return ns_M_G0, 0

    if iContinuity == 1:
        if areInOutEpsilonEqual(ns_M_BeforeAnyMatching, ns_M_G1): return
        return ns_M_G1, 1


    ns_M_G2 = setContinuity_G2(
        ns_M_BeforeAnyMatching=ns_M_BeforeAnyMatching,
        ns_M_In=ns_M_G1,
        side_M=side_Both,
        nurbs_R=nurbs_R,
        side_R=side_Both,
        bModifyRowEnd_T0=True,
        bModifyRowEnd_T1=True,
        bDebug=bDebug,
        bAddRefs=bAddRefs,
        )

    if ns_M_G2 is None:
        if areInOutEpsilonEqual(ns_M_BeforeAnyMatching, ns_M_G1): return
        return ns_M_G1, 1

    if areInOutEpsilonEqual(ns_M_BeforeAnyMatching, ns_M_G1): return
    return ns_M_G2, 2


def processObjRefs(objref_M, objref_R, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iContinuity = getOpt('iContinuity')
    iPreserveOtherEnd = getOpt('iPreserveOtherEnd')
    bUseUnderlyingIsoCrvs = getOpt('bUseUnderlyingIsoCrvs')
    bMaintainDegree = getOpt('bMaintainDegree')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    bAddRefs = getOpt('bAddRefs')


    def getGeomFromObjRef(objref):
        """ Returns BrepTrim, (wire) Curve, or None """

        if not isinstance(objref, rd.ObjRef): return

        trim = objref.Trim()

        if trim is not None:
            return trim

        edge = objref.Edge()
        if edge is not None:
            if len(edge.TrimIndices()) == 1:
                trim = edge.Brep.Trims[edge.TrimIndices()[0]]
                return trim

        return objref.Curve()


    geom_M_In = getGeomFromObjRef(objref_M)
    if geom_M_In is None:
        raise ValueError("{} is not valid for surface to modify.".format(geom_R_In.GetType().Name))
    if isinstance(geom_M_In, rg.BrepTrim):
        side_M = geom_M_In.IsoStatus
        ns_M = geom_M_In.Face.UnderlyingSurface().ToNurbsSurface()
        if ns_M.IsRational:
            makeNonRational(ns_M)
            if ns_M.IsRational:
                print("Warning... Surface to modify is rational.  Check results.")
    else:
        raise ValueError("{} not supported.".format(geom_M_In))


    geom_R_In = getGeomFromObjRef(objref_R)
    if geom_R_In is None:
        raise ValueError("{} is not valid for reference.".format(geom_R_In.GetType().Name))

    bMatchWithParamsAligned = areParamsAlignedPerPickPts(objref_M, objref_R)
    #sEval = "bMatchWithParamsAligned"; print(sEval+':',eval(sEval))


    rc = createSurface(
        ns_M_In=ns_M,
        side_M=side_M,
        geom_R_In=geom_R_In,
        bMatchWithParamsAligned=bMatchWithParamsAligned,
        iContinuity=iContinuity,
        iPreserveOtherEnd=iPreserveOtherEnd,
        bUseUnderlyingIsoCrvs=bUseUnderlyingIsoCrvs,
        bMaintainDegree=bMaintainDegree,
        bEcho=bEcho,
        bDebug=bDebug,
        bAddRefs=bAddRefs,
        )

    if rc is None:
        print("Surface could not be created.")
        return

    ns_Ret, iG = rc

    if bEcho:
        print("Continuity was modified toward G{}.".format(iG))

    if not bReplace:
        gB_Out = sc.doc.Objects.AddSurface(ns_Ret)
        if gB_Out == gB_Out.Empty:
            print("Could not add modified surface.")
        else:
            print("Surface was added.")
    else:
        rdB_In = objref_M.Object()
        rgB_In = objref_M.Brep()
        rgB_WIP = rgB_In.DuplicateBrep()
        rgB_WIP.Faces.RemoveAt(objref_M.Face().FaceIndex)
        rgBs_Out = rg.Brep.CreateBooleanUnion([rgB_WIP], sc.doc.ModelAbsoluteTolerance)
        rgB_WIP.Dispose()

        if rgB_In.Faces.Count == 1:
            rgB_Out = ns_Ret.ToBrep()
            if sc.doc.Objects.Replace(objref_M.ObjectId, rgB_Out):
                gB_Out = objref_M.ObjectId
                print("Replaced monoface brep with new surface.")
            else:
                print("Could not replace monoface brep with new surface.")
        else:
            attr = rdB_In.Attributes
            gB_Out = sc.doc.Objects.AddSurface(ns_Ret, attr)
            if gB_Out == gB_Out.Empty:
                print("Could not add modified surface.")
            else:
                if len(rgBs_Out) == 1:
                    if sc.doc.Objects.Replace(objref_M.ObjectId, rgBs_Out[0]):
                        print("Added new surface and deleted face of brep.")
                    else:
                        print("Added new surface but could not delete face of brep.")
                else:
                    gBs_Out = []
                    for rgB in rgBs_Out:
                        gBs_Out.append(sc.doc.Objects.AddBrep(rgB, attr))

                    if gBs_Out[0].Empty in gBs_Out:
                        for gB_Out in gBs_Out:
                            if gB_Out != gB_Out.Empty:
                                sc.doc.Objects.Delete(objectId=gB_Out, quiet=False)
                        print("Added new surface but could not delete face of brep.")
                    else:
                        sc.doc.Objects.Delete(rdB_In)
                        print("Added new surface and deleted face of brep." \
                            "  Remainder of brep is now {} breps.".format(len(gBs_Out)))

    return gB_Out


def main():

    rc = getInput_ToModify()
    if rc is None: return

    (
        objref_M,
        bEcho,
        bDebug,
       ) = rc


    rc = getInput_Ref(objref_M)
    if rc is None: return

    (
        objref_R,
        iContinuity,
        iPreserveOtherEnd,
        bUseUnderlyingIsoCrvs,
        bMaintainDegree,
        bReplace,
        bEcho,
        bDebug,
        bAddRefs,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBrep_Ret = processObjRefs(
        objref_M,
        objref_R,
        iContinuity=iContinuity,
        iPreserveOtherEnd=iPreserveOtherEnd,
        bUseUnderlyingIsoCrvs=bUseUnderlyingIsoCrvs,
        bMaintainDegree=bMaintainDegree,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        bAddRefs=bAddRefs,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()