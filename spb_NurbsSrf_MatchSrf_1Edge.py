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
211008-21: Bug fix when transferring knot vector.
        Refactored for easier use as module by other scripts.

TODO:
    Allow to use trimmed surface when surface side is selected.
    Check continuity before modifying surface.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum


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

    key = 'bAddPts'; keys.append(key)
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
        print "Face is not from a monoface brep."
        return False

    rgF = rgB.Faces[0]
    rgS = rgF.UnderlyingSurface()

    if not isinstance(rgS, rg.NurbsSurface):
        print "Underlying surface is a {}.".format(rgS.GetType().Name)
        return False

    ns = rgS

    if ns.IsClosed(0) or ns.IsClosed(1):
        print "Closed surfaces are not yet supported."
        return False

    if not rgF.IsSurface:
        print "Face is not a full surface."
        return False

    if not rgB.Trims.Count == 4:
        print "There are {} trims on this face.  Untrim it first.".format(
            rgB.Trims.Count)
        return False

    # Reject any NS with singular trims.
    iCt_Singulars = 0
    for iT in 0,1,2,3:
        if rgB.Trims[iT].TrimType == rg.BrepTrimType.Singular:
            iCt_Singulars += 1
    if iCt_Singulars > 0:
        print "There are {} singular trims on this surface.".format(
            iCt_Singulars)
        return False

    return True


def printObjRefObjs(objref):
    for sMethod in (
        'Brep', 'Curve', 'CurveParameter', 'Edge', 'Face', 'Geometry',
        'Object', 'Point', 'SelectionMethod', 'SelectionPoint', 'Surface',
        'SurfaceParameter', 'Trim'
    ):
        print "  .{}(): {}".format(sMethod, eval("objref."+sMethod+"()"))


def getInput_ToModify():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edge of untrimmed NURBS surface to modify")

    go.GeometryFilter = rd.ObjectType.EdgeFilter


    def geomFilter_TrimOfFull4TrimNS(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if isinstance(geom, rg.BrepTrim):
            rgT = geom
            rgB = rgT.Brep
            if rgB.Faces.Count > 1:
                return False
        else:
            return False

        rgF = rgB.Faces[0]
        rgS = rgF.UnderlyingSurface()

        if not isinstance(rgS, rg.NurbsSurface):
            return False

        ns = rgS

        if ns.IsClosed(0) or ns.IsClosed(1):
            return False

        if not rgF.IsSurface:
            return False

        if not rgB.Trims.Count == 4:
            return False

        if containsShortOrCollapsedBorders(ns):
            return False

        return True

    go.SetCustomGeometryFilter(geomFilter_TrimOfFull4TrimNS)

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
            addOption('bAddPts')

        res = go.Get()
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            #print objref.Face()
            #print objref.Surface()
            #print objref.Trim()

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


    def geomFilter_Custom(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if rdObj.Id == objref_SrfToMod.ObjectId:
            return False

        if isinstance(geom, rg.BrepEdge):
            edge = geom
            iFs = edge.AdjacentFaces()
            if len(iFs) != 1:
                return False
            brep = edge.Brep
            srf = brep.Faces[iFs[0]].UnderlyingSurface()
            if srf.IsClosed(0) or srf.IsClosed(1):
                return False
            if srf.IsRational:
                return False
            if containsShortOrCollapsedBorders(srf):
                return False
            return True
        elif isinstance(geom, rg.Curve):
            crv = geom
            if crv.IsClosed:
                return False
            if crv.IsRational:
                return False
            return True
        else:
            print "What happened?"
            return False

    go.SetCustomGeometryFilter(geomFilter_Custom)


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
            addOption('bAddPts')

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
                Opts.values['bAddPts'],
                )

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


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
            #if trim.IsoStatus in (W, S, N, E):
            #    return trim

    return objref.Curve()


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


    def isPickCloserToT0(objref):
        pt_Sel = objref.SelectionPoint()
        rdObj = objref.Object()
        if rdObj == rd.CurveObject:
            crv = objref.Curve()
            return crv.PointAtStart.DistanceTo(pt_Sel) < crv.PointAtEnd.DistanceTo(pt_Sel)
        trim = objref.Trim()
        if trim is None:
            edge = objref.Edge()
            if edge is None:
                raise ValueError("{} not valid.".format(printObjRefObjs(objref)))
            iTs = edge.TrimIndices()
            if len(iTs) != 1:
                raise ValueError("Edge has {} trims.".format(len(iTs)))
            trim = edge.Brep.Trims[iTs[0]]

        if trim.IsoStatus in (rg.IsoStatus.X, rg.IsoStatus.Y, rg.IsoStatus.None):
            crv = objref.Curve()
            return crv.PointAtStart.DistanceTo(pt_Sel) < crv.PointAtEnd.DistanceTo(pt_Sel)

        pts_Side = getSidePoints_AnySrfType(
            trim.Face.UnderlyingSurface(),
            trim.IsoStatus)
        return pts_Side[0].DistanceTo(pt_Sel) < pts_Side[1].DistanceTo(pt_Sel)

    if None in (objref_A.SelectionPoint(), objref_B.SelectionPoint()): return

    bPickedAtSideStart_A = isPickCloserToT0(objref_A)
    if bPickedAtSideStart_A is None: return

    bPickedAtSideStart_R = isPickCloserToT0(objref_B)
    if bPickedAtSideStart_R is None: return

    return bPickedAtSideStart_A == bPickedAtSideStart_R


def createTanSrfFromEdge(rgTrim):
    """
    """

    rgEdge = rgTrim.Edge
    srf = rgTrim.Face.UnderlyingSurface().ToNurbsSurface()

    ncA = rgEdge.ToNurbsCurve()

    ts = ncA.GrevilleParameters()

    pts_on_ncA = ncA.GrevillePoints()
    pts_on_ncB = []

    angle = Rhino.RhinoMath.ToRadians(90.0)

    dist = 10.0 * sc.doc.ModelAbsoluteTolerance


    def createEndPt(i, angle):
        t = ts[i]
        #pt_Start = ncA.PointAt(t)
        pt_Start = pts_on_ncA[i]

        bSuccess, u, v = srf.ClosestPoint(pt_Start)
        if not bSuccess:
            print "ClosestPoint failed." \
                "  Continuity increase will not occur on an edge."
            return

        vNormal = srf.NormalAt(u, v)
        bSuccess, frame = ncA.PerpendicularFrameAt(t)
        if not bSuccess:
            print "PerpendicularFrameAt failed." \
                "  Continuity increase will not occur on an edge."
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

    rc = srf.ClosestPoint(pt_Mid_Neg)
    if not rc[0]:
        raise Exception("ClosestPoint failed.")
    pfrel_Neg = rgTrim.Face.IsPointOnFace(*rc[1:])

    rc = srf.ClosestPoint(pt_Mid_Pos)
    if not rc[0]:
        raise Exception("ClosestPoint failed.")
    pfrel_Pos = rgTrim.Face.IsPointOnFace(*rc[1:])

    if pfrel_Neg == rg.PointFaceRelation.Interior and pfrel_Pos == rg.PointFaceRelation.Exterior:
        pts_on_ncB = [createEndPt(i, -angle) for i in range(len(ts))]
    elif pfrel_Neg == rg.PointFaceRelation.Exterior and pfrel_Pos == rg.PointFaceRelation.Interior:
        pts_on_ncB = [createEndPt(i,  angle) for i in range(len(ts))]
    else:
        raise Exception("PointFaceRelations are {} and {}".format(pfrel_Neg, pfrel_Pos))

    ncB = ncA.DuplicateCurve()
    ncB.SetGrevillePoints(pts_on_ncB)

    rgB_Loft = rg.Brep.CreateFromLoft(
            curves=[ncA, ncB],
            start=rg.Point3d.Unset,
            end=rg.Point3d.Unset,
            loftType=rg.LoftType.Straight,
            closed=False)

    if len(rgB_Loft) != 1:
        print "{} breps resulted from one curve." \
            "  Continuity increase will not occur on an edge."

    #sc.doc.Objects.AddBrep(rgB_Loft[0]); sc.doc.Views.Redraw()

    return rgB_Loft[0].Surfaces[0]


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
            if bEcho:
                print "Matching curve not found."
            return
    
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
                    print "Points non-coincident within edge tolerance found.  Check results."

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


def getNurbsGeomFromGeom(In, bEcho=True):
    """
    Returns 4-item list containing any combination of tuples of:
        (IsoStatus, NurbsSurface),
        (NurbsCurve, plane)
        (NurbsCurve, None)
    """

    if not isinstance(In, rg.GeometryBase):
        raise ValueError("{} not accepted.".format(In.GetType().Name))

    if not isinstance(In, rg.BrepTrim):
        return In.ToNurbsCurve(), None

    trim = In
    srf = trim.Face.UnderlyingSurface()

    if trim.IsoStatus not in (W,S,E,N):
        b, plane = srf.TryGetPlane(tolerance=1e-6)
        if b:
            return trim.Edge.ToNurbsCurve(), plane

    if trim.IsoStatus in (W,S,E,N):
        return trim.IsoStatus, srf.ToNurbsSurface()

    if trim.IsoStatus in (rg.IsoStatus.X, rg.IsoStatus.Y):
        print "TODO: Add routine for X or Y IsoStatus."

    # Try to create tangent surface for IsoStatus of None (and temporarily X and Y).
    ns_Tan = createTanSrfFromEdge(trim)
    if ns_Tan is None:
        return In.ToNurbsCurve(), None
    ncs_TanNS = [getIsoCurveOfSide(side, ns_Tan) for side in (W,S,E,N)]
    idxB = findMatchingCurveByEndPoints([trim.Edge], ncs_TanNS, bEcho=bEcho)[0]
    return (W,S,E,N)[idxB], ns_Tan


def createAlignedRefGeom(geom_R_In, ns_M, side_M, bMatchWithParamsAligned):
    """
    Returns either
        NurbsSurface with parameterizations aligned or
        NurbsCurve with parameterization aligned with Coons
    """


    def createAlignedRefSrfs(side_M, ns_R, side_R):
        #print side_M, side_R
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
            elif side_M == N and side_R == E:
                pass
            elif side_M == N and side_R == W:
                ns_R.Reverse(1, True)
            elif side_M == W and side_R == N:
                ns_R.Reverse(0, True)
            elif side_M == W and side_R == S:
                pass
            else:
                raise ValueError("What happened?")

        if not bMatchWithParamsAligned:
            ns_R.Reverse(0 if side_M in (S, N) else 1, True)


    geom_R = geom_R_In
    #print geom_R
    if isinstance(geom_R[1], rg.NurbsSurface):
        side_R, ns_R = geom_R
        rc = createAlignedRefSrfs(side_M, ns_R, side_R)
        return ns_R, None

    c_A = getIsoCurveOfSide(side_M, ns_M)
    c_R = geom_R[0]
    if not rg.Curve.DoDirectionsMatch(c_A, c_R):
        c_R.Reverse()

    return c_R, geom_R[1]


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


def createSurface(ns_M, side_M, geom_R, bMatchWithParamsAligned=True, **kwargs):# iContinuity=1, iPreserveOtherEnd=1, bUseUnderlyingIsoCrvs=False, bMaintainDegree=True, bDebug=False):
    """
    Parameters:
        ns_M: NurbsSurface to modify.
        side_M=side_M,
        geom_R: BrepTrim or Curve for wire=geom_R_In,
        bMatchWithParamsAligned: bool
        iContinuity=iContinuity,
        iPreserveOtherEnd=iPreserveOtherEnd,
        bUseUnderlyingIsoCrvs=bUseUnderlyingIsoCrvs,
        bMaintainDegree=bMaintainDegree,
        bDebug=bDebug,

    Returns on success: rg.NurbsSurface
    Returns on fail: None
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iContinuity = getOpt('iContinuity')
    iPreserveOtherEnd = getOpt('iPreserveOtherEnd')
    bUseUnderlyingIsoCrvs = getOpt('bUseUnderlyingIsoCrvs')
    bMaintainDegree = getOpt('bMaintainDegree')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    bAddPts = getOpt('bAddPts')


    geom_R_In = geom_R

    bValid, sLog = ns_M.IsValidWithLog()
    if not bValid: print sLog; return


    if side_M not in (W,S,E,N):
        return

    geom_R_Nurbs = None

    if isinstance(geom_R_In, rg.BrepTrim):
        trim_In = geom_R_In

        domU, domV = None, None

        if (
            not bUseUnderlyingIsoCrvs and
            trim_In.IsoStatus != rg.IsoStatus.None
        ):
            rc = getSrfDomainsWithinIsoCrvEdgeEnds(trim_In)
            if rc is not None:
                domU, domV = rc

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
                        trim_In.Face.UnderlyingSurface().Domain(0).T1,
                        )
                    iso_Out = S
                else:
                    domV = rg.Interval(
                        trim_In.Face.UnderlyingSurface().Domain(0).T0,
                        trim_In.PointAt(trim_In.Domain.Mid).Y,
                        )
                    iso_Out = N

        if domU is None and domV is None:
            # No trimming.
            geom_R_Nurbs = getNurbsGeomFromGeom(geom_R_In)
        else:
            if domU is None:
                domU = trim_In.Face.UnderlyingSurface().Domain(0)
            if domV is None:
                domV = trim_In.Face.UnderlyingSurface().Domain(1)
            trimmed = trim_In.Face.UnderlyingSurface().Trim(domU, domV)
            geom_R_Nurbs = iso_Out, trimmed.ToNurbsSurface()

    if geom_R_Nurbs is None:
        geom_R_Nurbs = getNurbsGeomFromGeom(geom_R_In)
        if not geom_R_Nurbs: return


    rc = createAlignedRefGeom(
        geom_R_Nurbs, ns_M, side_M, bMatchWithParamsAligned)
    if rc is None: return
    nurbs_R, plane = rc

    #sc.doc.Objects.AddSurface(nurbs_R) if isinstance(nurbs_R, rg.NurbsSurface) else sc.doc.Objects.AddCurve(nurbs_R); sc.doc.Views.Redraw(); return


    bValid, sLog = nurbs_R.IsValidWithLog()
    if not bValid: print sLog; return

    side_R = side_M
    side_Both = side_M

    #sc.doc.Objects.AddSurface(nurbs_R); sc.doc.Views.Redraw(); return
    #sEval = "bMatchWithParamsAligned"; print sEval+':',eval(sEval)
    #print side_M, side_R



    ns_M_In = ns_M # Need this?
    ns_M_Out = ns_M.Duplicate()

    ns_R_In = nurbs_R.Duplicate()



    # Along matching side of surface, match reference's degree if larger.
    # Increase degree if necessary.

    transferHigherDegree(
            nurbs_R, ns_M_Out, side_Both, side_Both)

    transferHigherDegree(
            ns_M_Out, nurbs_R, side_Both, side_Both)

    transferDomain(
            ns_M_Out, nurbs_R, side_Both, side_Both)





    iDir_MatchSide_R = int(side_Both in (rg.IsoStatus.East, rg.IsoStatus.West))
    #iDeg_MatchSide_R = ns_R.Degree(iDir_MatchSide_R)

    iDir_MatchSide_A = int(side_Both in (rg.IsoStatus.East, rg.IsoStatus.West))
    #iDeg_MatchSide_A = ns_M_Out.Degree(iDir_MatchSide_A)



    def increaseDegree(ns, iDir, iDeg):
        if iDir:
            return ns.IncreaseDegreeV(iDeg)
        else:
            return ns.IncreaseDegreeU(iDeg)

    def getDirString(iDir):
        return 'V' if iDir else 'U'


    #if iDeg_MatchSide_A < iDeg_MatchSide_R:
    #    if increaseDegree(ns_M_Out, iDir_MatchSide_A, iDeg_MatchSide_R):
    #        print "Increased degree of surface in {} direction from {} to {}.".format(
    #            getDirString(iDir_MatchSide_A),
    #            iDeg_MatchSide_A,
    #            iDeg_MatchSide_R)
    #    else:
    #        raise Exception("Could not raise degree of surface to modify.")

    #        # Update variables since degree was increased.
    #        iDeg_MatchSide_A = ns_M_Out.Degree(iDir_MatchSide_A)

    #elif iDeg_MatchSide_A > iDeg_MatchSide_R:
    #    ns_R = ns_R_In.Duplicate()
    #    if increaseDegree(ns_R, iDir_MatchSide_R, iDeg_MatchSide_A):
    #        if bDebug:
    #            print "Increased degree of reference surface in {} direction from {} to {}.".format(
    #                getDirString(iDir_MatchSide_R),
    #                iDeg_MatchSide_R,
    #                iDeg_MatchSide_A)
    #    else:
    #        raise Exception("Could not raise degree of reference surface.")

    #        # Update variables since degree was increased.
    #        iDeg_MatchSide_R = ns_R.Degree(iDir_MatchSide_R)

    bValid, sLog = ns_M_Out.IsValidWithLog()
    if not bValid: print sLog; return


    # Along match side of surface, match reference's knot vector since both surface degrees are the same.

    # First, transfer unique knots to Coons so that the latter can contain all unqiue knots
    # to transfer to reference surfaces.
    transferUniqueKnotVector(nurbs_R, ns_M_Out, side_Both, side_Both)

    #sc.doc.Objects.AddSurface(ns_M_Out)#; sc.doc.Views.Redraw(); return


    # Back to reference surfaces.
    transferUniqueKnotVector(ns_M_Out, nurbs_R, side_Both, side_Both)

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


    bValid, sLog = ns_M_Out.IsValidWithLog()
    if not bValid:
        print sLog; return
    #sc.doc.Objects.AddSurface(ns_M_Out); sc.doc.Views.Redraw(); return


    # From side of surface, increase point count as needed.
    # Points needed for continuity change: iContinuity + 1
    # Points needed to preserve other side: iPreserveOtherEnd
    iDir_FromMSide_R = int(not iDir_MatchSide_R)

    iDir_FromMSide_A = int(not iDir_MatchSide_A)
    iDeg_FromMSide_A = ns_M_Out.Degree(iDir_FromMSide_A)


    def getPtCt(ns, iDir):
        return ns.Points.CountV if iDir else ns.Points.CountU

    iPtCt_FromMSide_A = getPtCt(ns_M_Out, iDir_FromMSide_A)


    iPtCt_FromMSide_A_Min = iContinuity + 1 + iPreserveOtherEnd + 1

    if not bMaintainDegree:# and ns_M_Out.SpanCount(iDir_FromMSide_A) == 1:
        # Create single-spanned at minimal degree.
        iPtCt_FromMSide_A_ToAdd = iPtCt_FromMSide_A_Min - iPtCt_FromMSide_A
        iDeg_FromSide_A_Min = iPtCt_FromMSide_A_Min - 1
        if iDeg_FromMSide_A < iDeg_FromSide_A_Min:
            if increaseDegree(ns_M_Out, iDir_FromMSide_A, iDeg_FromSide_A_Min):
                print "Increased degree of surface in {} direction from {} to {}.".format(
                    getDirString(iDir_FromMSide_A),
                    iDeg_FromMSide_A,
                    iDeg_FromSide_A_Min)

                # Update variables since degree was increased.
                iDeg_FromMSide_A = ns_M_Out.Degree(iDir_FromMSide_A)
                iPtCt_FromMSide_A = getPtCt(ns_M_Out, iDir_FromMSide_A)

    else:
        # Maintain degree, increasing span count.
        iPtCt_FromMSide_A_ToAdd = iPtCt_FromMSide_A_Min - iPtCt_FromMSide_A

        if iPtCt_FromMSide_A_ToAdd > 0:
            
            iKnotCt_FromSide_A_ToAdd = iPtCt_FromMSide_A_ToAdd

            coll_Knots = getKnots(ns_M_Out, iDir_FromMSide_A) # Reference to collection.  Do not convert to list.

            iKnotCt_FromSide_A_Required = iPtCt_FromMSide_A_ToAdd + coll_Knots.Count

            # TODO: Use the following only to keep uniformity.  Otherwise, only add knots as needed.
            tsK_New = []
            iCt_Knot_Start = coll_Knots.Count

            addKnotsUniformly(
                coll_Knots=coll_Knots,
                iCt_Knot_Target=iKnotCt_FromSide_A_Required
                )

            print "Increased knot count from {} to {}.".format(
                    iCt_Knot_Start, coll_Knots.Count)
    
            #print coll_Knots.CreateUniformKnots(knotSpacing=1.0)
    
            ## Only increase degree to an odd value.
            #if iDegIncr % 2 == 0:
            #    iDegIncr += 1
            #if ns_M_Out.IncreaseDegreeV(iDegIncr):
            #    print "Increase degree of surface in V direction from {} to {}.".format(
            #        iDeg_FromMSide_A, iDegIncr)


    ptsA_Start = [cp.Location for cp in ns_M_Out.Points]
    ptsM = ptsA_Start[:]
    ptsR = [cp.Location for cp in nurbs_R.Points]

    # Points are ordered u0 v0, u0 v1, ..., un vn.

    #for pt in ptsM:
    #    sc.doc.Objects.AddPoint(pt)
    #sc.doc.Views.Redraw()
    #return



    # Check and if necessary, match point indices of reference surface.

    def getUvFrom1dList(ns, idxFlat):
        """
        Convert the 1-dimensional index of Points to the U and V indices for use .
        """
        return idxFlat // ns.Points.CountV, idxFlat % ns.Points.CountV

    bValid, sLog = ns_M_Out.IsValidWithLog()
    if not bValid: print sLog; return
    #sc.doc.Objects.AddSurface(ns_M_Out); sc.doc.Views.Redraw(); return


    # Set positional row.

    def getG0PtIndices(nurbsObj, side):
        # Position.
        if isinstance(nurbsObj, rg.NurbsSurface):
            cps = list(nurbsObj.Points)
            ctV = nurbsObj.Points.CountV
            ctAll = len(cps)
            if side == rg.IsoStatus.West:
                return range(ctV)
            elif side == rg.IsoStatus.East:
                return range(ctAll-ctV, ctAll)
            elif side == rg.IsoStatus.South:
                return range(0, ctAll, ctV)
            elif side == rg.IsoStatus.North:
                return range(ctV-1, ctAll, ctV)
        else: # NurbsCurve
            return range(nurbsObj.Points.Count)

    idxPts_G0_A = getG0PtIndices(ns_M_Out, side_Both)
    #print idxPts_G0_A
    idxPts_G0_R = getG0PtIndices(nurbs_R, side_Both)
    #print idxPts_G0_R

    for iG0_A, idxG0_R in zip(idxPts_G0_A, idxPts_G0_R):
        iU_G0_A, iV_G0_A = getUvFrom1dList(ns_M_Out, iG0_A)

        if isinstance(nurbs_R, rg.NurbsSurface):
            iU_0_R, iV_0_R = getUvFrom1dList(nurbs_R, idxG0_R)
            ns_M_Out.Points.SetPoint(
                iU_G0_A, iV_G0_A,
                point=ptsR[idxG0_R],
                weight=nurbs_R.Points.GetWeight(iU_0_R, iV_0_R))
        else:
            ns_M_Out.Points.SetPoint(
                iU_G0_A, iV_G0_A,
                point=ptsR[idxG0_R],
                weight=nurbs_R.Points.GetWeight(idxG0_R))


    #sc.doc.Objects.AddSurface(ns_M_Out); sc.doc.Views.Redraw(); return


    if iContinuity < 1:
        return ns_M_Out

    if plane is None and isinstance(nurbs_R, rg.NurbsCurve):
        return ns_M_Out


    # Update.
    ptsM = [cp.Location for cp in ns_M_Out.Points]

    #for pt in ptsM:
    #    sc.doc.Objects.AddPoint(pt)
    #sc.doc.Views.Redraw(); return


    # Set tangential row.

    ns_R = nurbs_R

    def getG1PtIndices(ns, side):
        # Tangency.
        cps = list(ns.Points)
        ctV = ns.Points.CountV
        ctAll = len(cps)
        if side == rg.IsoStatus.West:
            return range(ctV, 2*ctV)
        elif side == rg.IsoStatus.East:
            return range(ctAll-2*ctV, ctAll-ctV)
        elif side == rg.IsoStatus.South:
            return range(1, ctAll, ctV)
        elif side == rg.IsoStatus.North:
            return range(ctV-2, ctAll, ctV)


    idxPts_G1_A = getG1PtIndices(ns_M_Out, side_M)
    #print idxPts_G1_A


    if plane is not None:
        # Project tangential row to plane.
        for idx in range(len(idxPts_G1_A)):# iG1_A, iG0_A, iG1_R in zip(idxPts_G1_A, idxPts_G0_A, idxPts_G1_R):
            iG1_A = idxPts_G1_A[idx]

            pt_To = plane.ClosestPoint(ptsM[iG1_A])

            if pt_To == rg.Point3d.Unset:
                raise ValueError("Plane.ClosestPoint failed.")

            iUT_A, iVT_A = getUvFrom1dList(ns_M_Out, iG1_A)
            ns_M_Out.Points.SetPoint(
                iUT_A, iVT_A, pt_To,
                weight=ns_M_Out.Points.GetWeight(iUT_A, iVT_A))

        return ns_M_Out

    # Set tangential row colinear with reference.

    idxPts_G1_R = getG1PtIndices(ns_R, side_R) if plane is None else None
    #print idxPts_G1_R

    for iG1_A, iG0_A, iG1_R in zip(idxPts_G1_A, idxPts_G0_A, idxPts_G1_R):
        rgLine = rg.Line(ptsM[iG0_A], ptsR[iG1_R])
        pt_To = rgLine.ClosestPoint(ptsM[iG1_A], limitToFiniteSegment=False)
        iUT_A, iVT_A = getUvFrom1dList(ns_M_Out, iG1_A)
        ns_M_Out.Points.SetPoint(iUT_A, iVT_A, pt_To)


    # Update.
    ptsM = [cp.Location for cp in ns_M_Out.Points]


    # To obtain G1 continuity along edge, set the tangential row at a scale
    # of the G0-G1 point distances of the reference surface.

    # Determine the average ratio of the G0-G1 point distances of surface
    # to match to the same of reference surface.
    fDistRatios = []
    for iG1_A, iG0_A, iG1_R in zip(idxPts_G1_A, idxPts_G0_A, idxPts_G1_R):
        fDistA = ptsA_Start[iG1_A].DistanceTo(ptsA_Start[iG0_A])
        fDistR = ptsR[iG1_R].DistanceTo(ptsM[iG0_A])
        fRatio = fDistA / fDistR
        fDistRatios.append(fRatio)
    fAvgDistRatio = sum(fDistRatios) / float(len(fDistRatios))


    # Set tangential row using the multiple of the average ratio
    # to the G0-G1 point distances of reference surface.
    pts_G1_A_Out = []
    for iG1_A, iG0_A, iG1_R in zip(idxPts_G1_A, idxPts_G0_A, idxPts_G1_R):
        vG0G1_R = ptsM[iG0_A] - ptsR[iG1_R]
        vG0G1_A_toSet = fAvgDistRatio * vG0G1_R
        pt_To = ptsM[iG0_A] + vG0G1_A_toSet
        iUT_A, iVT_A = getUvFrom1dList(ns_M_Out, iG1_A)
        ns_M_Out.Points.SetPoint(iUT_A, iVT_A, pt_To)
        pts_G1_A_Out.append(pt_To)


    if iContinuity < 2:
        return ns_M_Out

    # Update.
    ptsM = [cp.Location for cp in ns_M_Out.Points]


    # G2 continuity.
    iDeg_FromMSide_R = ns_R.Degree(iDir_FromMSide_R)
    m_ForDeg_A = (iDeg_FromMSide_A - 1) / iDeg_FromMSide_A
    m_ForDeg_R = (iDeg_FromMSide_R - 1) / iDeg_FromMSide_R
    
    # Replaced by above.
    #    def multiplierForDegreeDifference(iDeg_FromMSide_A, iDeg_FromMSide_R):
    #        m_ForDegs = 1.0
    #        if iDeg_FromMSide_A < iDeg_FromMSide_R:
    #            for i in range(iDeg_FromMSide_A, iDeg_FromMSide_R):
    #                mDeg = (float(i)**2) / (float(i)**2 - 1.0)
    #                m_ForDegs *= mDeg
    #        elif iDeg_FromMSide_R < iDeg_FromMSide_A:
    #            for i in range(iDeg_FromMSide_R, iDeg_FromMSide_A):
    #                mDeg = (float(i)**2 - 1.0) / (float(i)**2)
    #                m_ForDegs *= mDeg
    #        else:
    #            # iDeg_FromMSide_A == iDeg_FromMSide_R
    #            pass
    #        return m_ForDegs
    #
    #    m_ForDegs = multiplierForDegreeDifference(iDeg_FromMSide_A, iDeg_FromMSide_R)


    def multiplierForAdjacentSimpleKnots(ns, side):
        
        iDir = bool(side in (rg.IsoStatus.South, rg.IsoStatus.North))
        
        # Calculate multiplier for knot vector.
        iCt_Span = ns.SpanCount(iDir)
        if iCt_Span == 1:
            m_ForKnots = 1.0
        else:
            knots = ns.KnotsV if iDir else ns.KnotsU
    
            if side in (rg.IsoStatus.East, rg.IsoStatus.North):
                iKnot_NextToEnd = knots.Count - iDeg_FromMSide_A - 1
            else:
                iKnot_NextToEnd = iDeg_FromMSide_A
    
            multy_AdjKnot = knots.KnotMultiplicity(iKnot_NextToEnd)
            
            if multy_AdjKnot > 1:
                m_ForKnots = 1.0
            else:
                t_Knots_Span = ns.GetSpanVector(iDir)
    
                if side in (rg.IsoStatus.East, rg.IsoStatus.North):
                    fLen_spanEnd = t_Knots_Span[iCt_Span] - t_Knots_Span[iCt_Span-1]
                    fLen_spanAdj = t_Knots_Span[iCt_Span-1] - t_Knots_Span[iCt_Span-2]
                else:
                    fLen_spanEnd = t_Knots_Span[1] - t_Knots_Span[0]
                    fLen_spanAdj = t_Knots_Span[2] - t_Knots_Span[1]
    
                m_ForKnots = (fLen_spanEnd + fLen_spanAdj) / fLen_spanEnd
        return m_ForKnots

    m_ForKnots_A = multiplierForAdjacentSimpleKnots(ns_M_Out, side_M)
    m_ForKnots_R = multiplierForAdjacentSimpleKnots(ns_R, side_R)


    def getG2PtIndices(ns, side):
        # Curvature.
        cps = list(ns.Points)
        ctV = ns.Points.CountV
        ctAll = len(cps)
        if side == rg.IsoStatus.West:
            return range(2*ctV, 3*ctV)
        elif side == rg.IsoStatus.East:
            return range(ctAll-3*ctV, ctAll-2*ctV)
        elif side == rg.IsoStatus.South:
            return range(2, ctAll, ctV)
        elif side == rg.IsoStatus.North:
            return range(ctV-3, ctAll, ctV)

    idxPts_G2_A = getG2PtIndices(ns_M_Out, side_M)
    #print idxPts_G2_R
    idxPts_G2_R = getG2PtIndices(ns_R, side_R)
    #print idxPts_G2_R

    for i in range(len(idxPts_G1_A)):
        iG1_A = idxPts_G1_A[i]
        iG0_R = idxPts_G0_R[i]
        iG1_R = idxPts_G1_R[i]
        iG2_R = idxPts_G2_R[i]
        iG2_A = idxPts_G2_A[i]

        a0 = r0 = ptsR[iG0_R]
        a1 = pts_G1_A_Out[i]
        r1 = ptsR[iG1_R]
        r2 = ptsR[iG2_R]
        m = ((a1 - r0).Length/(r1 - r0).Length)**2.0

        a2 = 2.0*a1 + -a0 + (
            m *
            (m_ForKnots_A / m_ForKnots_R) *
            (m_ForDeg_R / m_ForDeg_A) *
            (-2.0*r1 + r2 + r0)
            )

        if bDebug and bAddPts:
            sc.doc.Objects.AddPoint(r0)
            sc.doc.Objects.AddPoint(a1)
            sc.doc.Objects.AddPoint(a2)
            sc.doc.Objects.AddPoint(r1)
            sc.doc.Objects.AddPoint(r2)

        #sc.doc.Objects.AddPoint(a2)

        iUT_A, iVT_A = getUvFrom1dList(ns_M_Out, iG2_A)
        ns_M_Out.Points.SetPoint(iUT_A, iVT_A, a2)

    #sc.doc.Views.Redraw()

    return ns_M_Out


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
    bAddPts = getOpt('bAddPts')


    geom_M_In = getGeomFromObjRef(objref_M)
    if geom_M_In is None:
        print "Objref not provided for object to modify."
        return
    if isinstance(geom_M_In, rg.BrepTrim):
        side_M = geom_M_In.IsoStatus
        ns_M = geom_M_In.Face.UnderlyingSurface().ToNurbsSurface()
    else:
        print geom_M_In
        raise ValueError("{} not supported.".format(geom_M_In))


    geom_R_In = getGeomFromObjRef(objref_R)
    if geom_R_In is None:
        print "Objref not provided for object to modify."
        return
    if isinstance(geom_R_In, rg.BrepTrim):
        side_R = geom_R_In.IsoStatus
        ns_R = geom_R_In.Face.UnderlyingSurface()
    else:
        raise ValueError("{} is not valid for reference.".format(geom_R_In.GetType().Name))


    bMatchWithParamsAligned = areParamsAlignedPerPickPts(objref_M, objref_R)
    #sEval = "bMatchWithParamsAligned"; print sEval+':',eval(sEval)

    # ns_R.Reverse(int direction, bool inPlace)
    # Transpose(bool inPlace)

    ns_Ret = createSurface(
        ns_M=ns_M,
        side_M=side_M,
        geom_R=geom_R_In,
        bMatchWithParamsAligned=bMatchWithParamsAligned,
        iContinuity=iContinuity,
        iPreserveOtherEnd=iPreserveOtherEnd,
        bUseUnderlyingIsoCrvs=bUseUnderlyingIsoCrvs,
        bMaintainDegree=bMaintainDegree,
        bEcho=bEcho,
        bDebug=bDebug,
        bAddPts=bAddPts,
        )

    if ns_Ret is None:
        print "Surface could not be created."
        return

    epsilon = 1e-6
    if ns_Ret.EpsilonEquals(ns_M, epsilon):
        while True:
            eps_prev = epsilon
            epsilon /= 10.0
            if not ns_Ret.EpsilonEquals(ns_M, epsilon):
                sc.escape_test()
                break
        print "Output surface EpsilonEquals input within {}.  No change.".format(eps_prev)
        ns_Ret.Dispose()
        return

    if not bReplace:
        gB_Out = sc.doc.Objects.AddSurface(ns_Ret)
        if gB_Out == gB_Out.Empty:
            print "Could not add modified surface."
        else:
            print "Surface was added."
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
                print "Replaced monoface brep with new surface."
            else:
                print "Could not replace monoface brep with new surface."
        else:
            attr = rdB_In.Attributes
            gB_Out = sc.doc.Objects.AddSurface(ns_Ret, attr)
            if gB_Out == gB_Out.Empty:
                print "Could not add modified surface."
            else:
                if len(rgBs_Out) == 1:
                    if sc.doc.Objects.Replace(objref_M.ObjectId, rgBs_Out[0]):
                        print "Added new surface and deleted face of brep."
                    else:
                        print "Added new surface but could not delete face of brep."
                else:
                    gBs_Out = []
                    for rgB in rgBs_Out:
                        gBs_Out.append(sc.doc.Objects.AddBrep(rgB, attr))

                    if gBs_Out[0].Empty in gBs_Out:
                        for gB_Out in gBs_Out:
                            if gB_Out != gB_Out.Empty:
                                sc.doc.Objects.Delete(objectId=gB_Out, quiet=False)
                        print "Added new surface but could not delete face of brep."
                    else:
                        sc.doc.Objects.Delete(rdB_In)
                        print "Added new surface and deleted face of brep." \
                            "  Remainder of brep is now {} breps.".format(len(gBs_Out))

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
        bAddPts,
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
        bAddPts=bAddPts,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()
