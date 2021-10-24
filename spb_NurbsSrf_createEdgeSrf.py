"""
As an alternative to _EdgeSrf, this script:
1. Matches unionizes knot vectors more exactly.
2. Offers face-face continuity above G0.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""
"""
211013-23: Created.

TODO:
    Add support for 3-curve input.
    Improve continuity when input includes edges of rational surfaces.
    Allow continuity choice per edge with preview.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import spb_NurbsSrf_MatchSrf_1Edge as spb


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
    values[key] = 2
    names[key] = 'ContinuityTarget'
    listValues[key] = 'G0', 'G1', 'G2'
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

        if not idxOpt: print "Add option for {} failed.".format(key)

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        #if key == 'fPatchTol':
        #    if cls.riOpts[key].CurrentValue <= 1e-9:
        #        print "Invalid input for scale value."
        #        cls.riOpts[key].CurrentValue = cls.values[key]

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get edges with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edges or curves")

    go.GeometryFilter = rd.ObjectType.Curve
    #go.GeometryAttributeFilter = (
    #    go.GeometryAttributeFilter.SurfaceBoundaryEdge |
    #    go.GeometryAttributeFilter.WireCurve)

    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iContinuity')
        addOption('bEcho')
        addOption('bDebug')

        sc.doc.Views.Redraw()

        res = go.GetMultiple(minimumNumber=4, maximumNumber=4)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return (
                objrefs,
                Opts.values['iContinuity'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        # An option was selected.
        go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
        go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getFormattedDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


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


def getCurveOfNurbs(geom):
    if geom[1] is None:
        c = geom[0]
        if not isinstance(c, rg.NurbsCurve):
            raise ValueError("Not a NurbsCurve.")
        return c

    if isinstance(geom[1], rg.Plane):
        if not isinstance(geom[0], rg.NurbsCurve):
            raise ValueError("Not a NurbsCurve.")
        return geom[0]
    else:
        c = getIsoCurveOfSide(*geom)
        if not isinstance(c, rg.NurbsCurve):
            raise ValueError("Not a NurbsCurve.")
        return c


def isPointAtEitherCrvEnd(pt, cA):
    """
    """

    tol = 2.0 * sc.doc.ModelAbsoluteTolerance
    if cA.PointAtStart.DistanceTo(pt) <= tol:
        return True
    if cA.PointAtEnd.DistanceTo(pt) <= tol:
        return True

    return False


def doTrimsFormAClosedLoop(geoms):
    """
    Trims must
        Be of full surface edge.
        Form a closed loop.
    """

    cs = [getCurveOfNurbs(geom) for geom in geoms]
    found = [[False, False] for i in 0,1,2,3]
    tol = 2.0 * sc.doc.ModelAbsoluteTolerance

    for iA in 0,1,2:
        cA = cs[iA]
        if not found[iA][0]:
            ptA = cA.PointAtStart
            for iB in range(1+iA,4):
                cB = cs[iB]
                if ptA.DistanceTo(cB.PointAtStart) <= tol:
                    #print "Start of [{}] matches start of [{}].".format(iA, iB)
                    found[iA][0] = True
                    found[iB][0] = True
                    if cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= tol:
                        raise Exception("2 curves completely overlap.")
                    break
                if ptA.DistanceTo(cB.PointAtEnd) <= tol:
                    #print "Start of [{}] matches end of [{}].".format(iA, iB)
                    found[iA][0] = True
                    found[iB][1] = True
                    if cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= tol:
                        raise Exception("2 curves completely overlap.")
                    break
            else:
                #for c in cs: sc.doc.Objects.AddCurve(c)
                return False

        if not found[iA][1]:
            ptA = cA.PointAtEnd
            for iB in range(1+iA,4):
                cB = cs[iB]
                if ptA.DistanceTo(cB.PointAtStart) <= tol:
                    #print "End of [{}] matches start of [{}].".format(iA, iB)
                    found[iA][1] = True
                    found[iB][0] = True
                    if cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= tol:
                        raise Exception("2 curves completely overlap.")
                    break
                if ptA.DistanceTo(cB.PointAtEnd) <= tol:
                    #print "End of [{}] matches end of [{}].".format(iA, iB)
                    found[iA][1] = True
                    found[iB][1] = True
                    if cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= tol:
                        raise Exception("2 curves completely overlap.")
                    break
            else:
                #for c in cs: sc.doc.Objects.AddCurve(c)
                return False

    return True


def getTrimInterval_NurbsSrf(iso, ns, pt_Ends_ToRemain, pts_Xs):
    """
    """

    b, u, v = ns.ClosestPoint(pt_Ends_ToRemain)
    if not b:
        raise Exception("ClosestPoint failed.")

    if (abs(u - ns.Domain(0).T0) <= 1e-6) and (abs(v - ns.Domain(1).T0) <= 1e-6):
        uvs_DomainEndsToKeep = (u,v)
    elif (abs(u - ns.Domain(0).T1) <= 1e-6) and (abs(v - ns.Domain(1).T0) <= 1e-6):
        uvs_DomainEndsToKeep = (u,v)
    elif (abs(u - ns.Domain(0).T1) <= 1e-6) and (abs(v - ns.Domain(1).T1) <= 1e-6):
        uvs_DomainEndsToKeep = (u,v)
    elif (abs(u - ns.Domain(0).T0) <= 1e-6) and (abs(v - ns.Domain(1).T1) <= 1e-6):
        uvs_DomainEndsToKeep = (u,v)
    else:
        uvs_DomainEndsToKeep = None

    #sEval = "t_DomainEndToKeep"; print sEval+':',eval(sEval)

    uvs_InsideDomain = []
    for pt in pts_Xs:
        b, u, v = ns.ClosestPoint(pt)
        if not b:
            raise Exception("ClosestPoint failed.")
        uvs_InsideDomain.append((u,v))


    if iso in (N, S):
        us = [uv[0] for uv in uvs_InsideDomain]
        if uvs_DomainEndsToKeep is not None:
            us += [uvs_DomainEndsToKeep[0]]
        us.sort()
        return rg.Interval(us[0], us[1]), ns.Domain(1)
    else:
        vs = [uv[1] for uv in uvs_InsideDomain]
        if uvs_DomainEndsToKeep is not None:
            vs += [uvs_DomainEndsToKeep[1]]
        vs.sort()
        return ns.Domain(0), rg.Interval(vs[0], vs[1])


def getTrimInterval_NurbsCrv(nc, pt_Ends_ToRemain, pts_Xs):
    """
    """

    b, t = nc.ClosestPoint(pt_Ends_ToRemain)
    if not b:
        raise Exception("ClosestPoint failed.")

    if (abs(t - nc.Domain.T0) <= 1e-6):
        t_DomainEndToKeep = t
    elif (abs(t - nc.Domain.T1) <= 1e-6):
        t_DomainEndToKeep = t

    ts_InsideDomain = []
    for pt in pts_Xs:
        b, t = nc.ClosestPoint(pt)
        if not b:
            raise Exception("ClosestPoint failed.")
        ts_InsideDomain.append(t)


    ts = ts_InsideDomain
    if pt_Ends_ToRemain:
        ts += [t_DomainEndToKeep]

    ts.sort()

    return rg.Interval(ts[0], ts[1])


def trimNurbsAsNeededToCloseLoop(ngs):
    """
    Trims must
        Be of full surface edge.
        Form a closed loop.
    """

    ngs_In = ngs

    cs = [getCurveOfNurbs(geom) for geom in ngs_In]
    tol = 2.0 * sc.doc.ModelAbsoluteTolerance

    #bTrims_NeedIntersect = [False]*4

    bFounds = [[False, False] for i in 0,1,2,3]

    def findEndsWithMatches():
        for iA in 0,1,2:
            cA = cs[iA]
            if not bFounds[iA][0]:
                ptA = cA.PointAtStart
                for iB in range(1+iA,4):
                    cB = cs[iB]
                    if ptA.DistanceTo(cB.PointAtStart) <= tol:
                        #print "Start of [{}] matches start of [{}].".format(iA, iB)
                        bFounds[iA][0] = True
                        bFounds[iB][0] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= tol:
                            raise Exception("2 curves completely overlap.")
                        break
                    if ptA.DistanceTo(cB.PointAtEnd) <= tol:
                        #print "Start of [{}] matches end of [{}].".format(iA, iB)
                        bFounds[iA][0] = True
                        bFounds[iB][1] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= tol:
                            raise Exception("2 curves completely overlap.")
                        break
    
            if not bFounds[iA][1]:
                ptA = cA.PointAtEnd
                for iB in range(1+iA,4):
                    cB = cs[iB]
                    if ptA.DistanceTo(cB.PointAtStart) <= tol:
                        #print "End of [{}] matches start of [{}].".format(iA, iB)
                        bFounds[iA][1] = True
                        bFounds[iB][0] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= tol:
                            raise Exception("2 curves completely overlap.")
                        break
                    if ptA.DistanceTo(cB.PointAtEnd) <= tol:
                        #print "End of [{}] matches end of [{}].".format(iA, iB)
                        bFounds[iA][1] = True
                        bFounds[iB][1] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= tol:
                            raise Exception("2 curves completely overlap.")
                        break

    findEndsWithMatches() # Modifies bFounds.

    bTrims_NeedIntersect = [not (bFounds[i][0] and bFounds[i][1]) for i in 0,1,2,3]

    if not any(bTrims_NeedIntersect):
        raise Exception("No trims can be split to create good input for surface.")

    if sum(bTrims_NeedIntersect) == 1:
        #print bTrims_NeedIntersect
        raise Exception("Only 1 trim to split.")


    iTs_NeedSplit = []
    pts_Xs = [[],[],[],[]]

    def getIntersectsForTrimming():
        for iA in 0,1,2:
            if not bTrims_NeedIntersect[iA]:
                continue
            cA = cs[iA]
            for iB in range(1+iA,4):
                if not bTrims_NeedIntersect[iB]:
                    continue
                cB = cs[iB]
                crvinters = rg.Intersect.Intersection.CurveCurve(
                    cA, cB, tolerance=tol, overlapTolerance=0.0)
                if crvinters.Count == 0: continue # to next curve.
                for crvinter in crvinters:
                    pt = crvinter.PointA
                    if not isPointAtEitherCrvEnd(pt, cA):
                        iTs_NeedSplit.append(iA)
                        pts_Xs[iA].append(pt)
                    if not isPointAtEitherCrvEnd(pt, cB):
                        iTs_NeedSplit.append(iB)
                        pts_Xs[iB].append(pt)

                    #sc.doc.Objects.AddPoint(pt)
                    #sc.doc.Views.Redraw()
                    #addNewPointToIntersectList(pt)

    getIntersectsForTrimming() # Modifies iTs_NeedSplit and pts_Xs.

    #print found
    #print pts_Xs

    ngs_Out = []

    for i in 0,1,2,3:
        ng_In = ngs_In[i]
        if len(pts_Xs[i]) == 0:
            ngs_Out.append(ng_In)
            continue

        pt_End_ToRemain = None
        if bFounds[i][0]:
            pt_End_ToRemain = cs[i].PointAtStart
        if bFounds[i][1]:
            if pt_End_ToRemain is not None:
                sEval = "len(pts_Ends)"
                raise ValueError("{}: {}".format(sEval, eval(sEval)))

            pt_End_ToRemain = cs[i].PointAtEnd


        if isinstance(ng_In[1], rg.NurbsSurface):
            iso, ns = ng_In
            intvls_Trim = getTrimInterval_NurbsSrf(iso, ns, pt_End_ToRemain, pts_Xs[i])
            ns_Trimmed = ns.Trim(intvls_Trim[0], intvls_Trim[1])
            ngs_Out.append((iso, ns_Trimmed))
        else:
            nc = ng_In[0]
            intvl_Trim = getTrimInterval_NurbsCrv(nc, pt_End_ToRemain, pts_Xs[i])
            nc_Trimmed = nc.Trim(intvl_Trim[0], intvl_Trim[1])
            ngs_Out.append((nc_Trimmed, ng_In[1]))


    return ngs_Out


def createAlignedRefGeom_PerStart(geom_R_In, ns_M, side_M):
    """
    Variation to a function in spb_Match1 in that
    the start locations of the edges are used to determine whether
    NS or NC needs to be reversed.

    Returns either
        NurbsSurface with parameterizations aligned with ns_M along side_M or
        NurbsCurve with parameterization aligned with ns_M
    """


    def areSrfParamsAlignedAlongMatchingSides(ns_A, side_A, ns_B, side_B):
        """
        """

        def getSideEndPtsInAscendingSrfParam(ns, side):
            cps = ns.Points
            pt_SW = cps.GetControlPoint(           0,            0).Location
            pt_SE = cps.GetControlPoint(cps.CountU-1,            0).Location
            pt_NE = cps.GetControlPoint(cps.CountU-1, cps.CountV-1).Location
            pt_NW = cps.GetControlPoint(           0, cps.CountV-1).Location
            if side == E:
                return pt_SE, pt_NE
            elif side == W:
                return pt_SW, pt_NW
            elif side == S:
                return pt_SW, pt_SE
            elif side == N:
                return pt_NW, pt_NE
            else:
                raise ValueError("{} is not a valid side.".format(side))


        pts_A = getSideEndPtsInAscendingSrfParam(ns_A, side_A)
        pts_B = getSideEndPtsInAscendingSrfParam(ns_B, side_B)

        #for pt in pts_A: sc.doc.Objects.AddPoint(pt)
        #for pt in pts_B: sc.doc.Objects.AddPoint(pt)

        tol = 2.0 * sc.doc.ModelAbsoluteTolerance

        if (
            (pts_A[0].DistanceTo(pts_B[0]) < tol) and 
            (pts_A[1].DistanceTo(pts_B[1]) < tol)
        ):
            return True

        if (
            (pts_A[0].DistanceTo(pts_B[1]) < tol) and 
            (pts_A[1].DistanceTo(pts_B[0]) < tol)
        ):
            return False

        #sc.doc.Objects.AddSurface(ns_A)
        #sc.doc.Objects.AddSurface(ns_B)
        #sc.doc.Views.Redraw()


        raise Exception("The 2 sides' positions do not match.")


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
                raise Exception("What happened?")

        if (not areSrfParamsAlignedAlongMatchingSides(
            ns_M, side_M, ns_R, side_M)
        ):
            ns_R.Reverse(0 if side_M in (S, N) else 1, True)


    if isinstance(geom_R_In[1], rg.NurbsSurface):
        side_R, ns_R = geom_R_In
        rc = createAlignedRefSrfs(side_M, ns_R, side_R)
        return ns_R # Don't bother returning IsoStatus since R now matches that of M.

    # Reference alignment is per Edge.
    c_M = getIsoCurveOfSide(side_M, ns_M)
    c_R = geom_R_In[0]
    if not rg.Curve.DoDirectionsMatch(c_M, c_R):
        c_R.Reverse()

    return c_R # Don't bother returning tuple since none is needed for NS.


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


def getUvIdxFromNsPoint1dIdx(ns, idxFlat):
    return idxFlat // ns.Points.CountV, idxFlat % ns.Points.CountV


def createSurface(rhCrvs_In, **kwargs):
    """
    rhCrvs_In: Can be ObjRefs or Geometry of any Curve, including BrepEdge.
    Returns NurbsSurface on success.
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iContinuity = getOpt('iContinuity')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if len(rhCrvs_In) != 4:
        return None, "{} curves provided.  Need exactly 4.".format(len(geoms_Nurbs))


    tol0 = Rhino.RhinoMath.ZeroTolerance


    geoms_In = []
    for rhCrv_In in rhCrvs_In:
        if isinstance(rhCrv_In, rd.ObjRef):
            rc = spb.getGeomFromObjRef(rhCrv_In)
            if rc is None: return
            geoms_In.append(rc)
        else:
            geoms_In.append(rhCrv_In)


    def areAllBrepTrimsUnique(geoms):
        isos = []
        srfs = []

        for geom in geoms:
            if not isinstance(geom, rg.BrepTrim): continue

            trim = geom

            iso = trim.IsoStatus
            if iso not in (W,S,E,N): continue

            srf = trim.Face.UnderlyingSurface()

            for iso_InList, srf_InList in zip(isos, srfs):
                if rg.GeometryBase.GeometryReferenceEquals(srf_InList, srf):
                    if iso_InList == iso:
                        return False

            isos.append(iso)
            srfs.append(srf)

        return True

    if not areAllBrepTrimsUnique(geoms_In):
        return None, "2 edges of the same isostatus of the same surface were selected."


    geoms_Nurbs = [spb.getNurbsGeomFromGeom(geom) for geom in geoms_In]
    if not geoms_Nurbs: return

    #for nc in [getCurveOfNurbs(geom) for geom in geoms_As_Nurbs]:
    #    sc.doc.Objects.AddCurve(nc)
    #return


    def doAnyCurvesCompletelyOverlap(geoms):
        """
        """

        cs = [getCurveOfNurbs(geom) for geom in geoms]
        found = [[False, False] for i in 0,1,2,3]
        tol = 2.0 * sc.doc.ModelAbsoluteTolerance

        for iA in 0,1,2:
            cA = cs[iA]
            if not found[iA][0]:
                ptA = cA.PointAtStart
                for iB in range(1+iA,4):
                    cB = cs[iB]
                    if ptA.DistanceTo(cB.PointAtStart) <= tol:
                        found[iA][0] = True
                        found[iB][0] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= tol:
                            return True
                        break
                    if ptA.DistanceTo(cB.PointAtEnd) <= tol:
                        found[iA][0] = True
                        found[iB][1] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= tol:
                            return True
                        break

            if not found[iA][1]:
                ptA = cA.PointAtEnd
                for iB in range(1+iA,4):
                    cB = cs[iB]
                    if ptA.DistanceTo(cB.PointAtStart) <= tol:
                        found[iA][1] = True
                        found[iB][0] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtEnd) <= tol:
                            return True
                        break
                    if ptA.DistanceTo(cB.PointAtEnd) <= tol:
                        found[iA][1] = True
                        found[iB][1] = True
                        if cA.PointAtEnd.DistanceTo(cB.PointAtStart) <= tol:
                            return True
                        break

        return False

    if doAnyCurvesCompletelyOverlap(geoms_Nurbs):
        return None, "Some edges/trims completely overlap."

    if not doTrimsFormAClosedLoop(geoms_Nurbs):
        print "No matching endpoint of some curves found,",
        rc = trimNurbsAsNeededToCloseLoop(geoms_Nurbs)
        geoms_Nurbs_ClosedLoop = rc
        if not doTrimsFormAClosedLoop(geoms_Nurbs_ClosedLoop):
            print " and attempt ot split some trims at intersections failed."
            return
        print " but attempt ot split some trims at intersections succeeded."
    else:
        geoms_Nurbs_ClosedLoop = geoms_Nurbs



    cs = [getCurveOfNurbs(geom) for geom in geoms_Nurbs_ClosedLoop]

    #for c in cs:
    #    sc.doc.Objects.AddCurve(c)
    #sc.doc.Views.Redraw(); return


    rgB_Res = rg.Brep.CreateEdgeSurface(cs)
    if rgB_Res is None:
        raise ValueError("CreateEdgeSurface returned None.")

    ns_Coons = rgB_Res.Surfaces[0]


    if ns_Coons.IsRational:
        spb.makeNonRational(ns_Coons)


    # Brep.CreateEdgeSurface (always?) produces a surface with the
    # largest degree in chain of surfaces in each direction but
    # may have missing knots and point distributions.
    # Both it and the surfaces from which the reference trims are taken
    # need to match along the relevant direction.
    # Create 4 modified surfaces to use for remainder of script.


    # Match references to the Coons.
    cs_C = [getIsoCurveOfSide(s, ns_Coons) for s in (W,S,E,N)]
    cs_R = [getCurveOfNurbs(geom) for geom in geoms_Nurbs_ClosedLoop]

    idx_R_per_CoonsWSEN = spb.findMatchingCurveByEndPoints(
        cs_C, cs_R, bEcho)
    if idx_R_per_CoonsWSEN is None: return

    geoms_Nurbs_SideMatched = [geoms_Nurbs_ClosedLoop[i] for i in idx_R_per_CoonsWSEN]

    # Align each reference to the Coons by side and parameterization.
    # This results in opposite u x v directions between each Refernce and the Coons.
    geoms_Nurbs_AR = {}

    for i, side in enumerate((W,S,E,N)):
        rc = createAlignedRefGeom_PerStart(
            geoms_Nurbs_SideMatched[i], ns_Coons, side)
        if rc is None: return
        geoms_Nurbs_AR[side] = rc


    # Match reference curve or surface degrees to Coons patch surface in relevant direction.

    bResults = []
    for side in W,S,E,N:
        bResult = spb.transferHigherDegree(
            ns_Coons, geoms_Nurbs_AR[side], side, side)
    bResults.append(bResult)
    if bDebug: print bResults

    #for ns in geoms_ARs:
    #    sc.doc.Objects.AddSurface(ns)
    #sc.doc.Views.Redraw()


    # Match reference srf relevant domains to Coons patch srf.
    bResults = []
    for side in W,S,E,N:
        bResult = spb.transferDomain(
            ns_Coons, geoms_Nurbs_AR[side], side, side)
    bResults.append(bResult)
    if bDebug: print bResults

    #for ns in geoms_ARs:
    #    sc.doc.Objects.AddSurface(ns)
    #sc.doc.Views.Redraw()



    # Match reference srf knot vectors to Coons patch srf.

    # First, transfer unique knots to Coons so that the latter can contain all unqiue knots
    # to transfer to reference surfaces.
    [transferUniqueKnotVector(geoms_Nurbs_AR[side], ns_Coons, side, side) for side in (S,N,W,E)]

    #sc.doc.Objects.AddSurface(ns_Coons)#; sc.doc.Views.Redraw(); return


    # Back to reference surfaces.
    [transferUniqueKnotVector(ns_Coons, geoms_Nurbs_AR[side], side, side) for side in (S,N,W,E)]

    #for ns in geoms_ARs:
    #    sc.doc.Objects.AddSurface(ns)
    #sc.doc.Views.Redraw(); return


    if iContinuity == 0:
        return ns_Coons

    if ns_Coons.Points.CountU < 4 and ns_Coons.Points.CountV < 4:
        if bEcho: print "Not enough points to modify continuity."
        return ns_Coons

    if bEcho:
        iCt = sum(isinstance(geoms_Nurbs_AR[key], rg.NurbsSurface) for key in geoms_Nurbs_AR)
        print "Modifying continuity of up to {} sides.".format(iCt)



    ptsM_PreG1 = [cp.Location for cp in ns_Coons.Points]

    pts_G1_corners = {}
    pts_G2_corners = {}

    ns_WIP = ns_Coons.Duplicate()


    def areThereEnoughPtsToModify(side, iG):
        if side in (W,E):
            if ns_Coons.Points.CountU < 2*(iG+1):
                if side == W:
                    print "Not enough CPs to modify to G{} along U.".format(iG)
                return False
        else:
            if ns_Coons.Points.CountV < 2*(iG+1):
                if side == S:
                    print "Not enough CPs to modify to G{} along V.".format(iG)
                return False
        return True


    for side in W,S,E,N:

        if not areThereEnoughPtsToModify(side, 1):
            continue

        geom_AR = geoms_Nurbs_AR[side]

        if isinstance(geom_AR, rg.Curve):
            continue # Since surface is already G0.

        ns_R = geom_AR

        #sc.doc.Objects.AddSurface(ns_R)

        idxPts = {} # Key is tuple(str('M' or 'R'), int(G continuity))
        for i in range(iContinuity+1):
            idxPts['M',i] = spb.getPtRowIndicesPerG(ns_Coons, side, i)
            idxPts['R',i] = spb.getPtRowIndicesPerG(ns_R, side, i)

        iRowLn = len(idxPts['M',0])


        # G0 row is already set.

        ns_M_G1 = spb.setContinuity_G1(
            ns_M_BeforeAnyMatching=ns_Coons,
            ns_M_In=ns_WIP,
            side_M=side,
            nurbs_R=ns_R,
            side_R=side,
            bModifyRowEnds=False,
            )

        if ns_M_G1 is None:
            continue

        #if side == S:
        #    sc.doc.Objects.AddSurface(ns_M_G1); sc.doc.Views.Redraw(); return

        # Update.
        ns_WIP.Dispose()
        ns_WIP = ns_M_G1
        pts_WIP = [cp.Location for cp in ns_WIP.Points]

        # Record the G1 point from each 2nd from end for later average.
        for i in (1, iRowLn-2):
            iM1 = idxPts['M',1][i]
            if iM1 in pts_G1_corners.keys():
                pts_G1_corners[iM1].append(pts_WIP[iM1])
            else:
                pts_G1_corners[iM1] = [pts_WIP[iM1]]
        #sc.doc.Objects.AddPoint(pts_WIP[idxPts['M',1][1]])
        #sc.doc.Objects.AddPoint(pts_WIP[idxPts['M',1][len(idxPts['M',1])-2]])

        if iContinuity == 2 and areThereEnoughPtsToModify(side, 2):
            ns_M_G2 = spb.setContinuity_G2(
                ns_M_BeforeAnyMatching=ns_Coons,
                ns_M_In=ns_WIP,
                side_M=side,
                nurbs_R=ns_R,
                side_R=side,
                bModifyRowEnds=False,
                bDebug=bDebug,
                bAddPts=True,
                )

            # Update.
            ns_WIP.Dispose()
            ns_WIP = ns_M_G2
            pts_WIP = [cp.Location for cp in ns_M_G2.Points]

            # Record the G2 point from each 2nd from end for later average.
            for i in (2, iRowLn-3):
                iM2 = idxPts['M',2][i]
                if iM2 in pts_G2_corners.keys():
                    pts_G1_corners[iM2].append(pts_WIP[iM2])
                else:
                    pts_G1_corners[iM2] = [pts_WIP[iM2]]
            #sc.doc.Objects.AddPoint(pts_WIP[idxPts['M',1][1]])
            #sc.doc.Objects.AddPoint(pts_WIP[idxPts['M',1][len(idxPts['M',1])-2]])

            # Return the 3rd G2 points on each end to their original positions.
            for idx in (2, iRowLn-3):
                iM2 = idxPts['M',2][idx]
                iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_Coons, iM2)
                ns_WIP.Points.SetPoint(
                    iUT_A, iVT_A, ptsM_PreG1[iM2],
                    weight=ns_Coons.Points.GetWeight(iUT_A, iVT_A))


        # Return the 2nd G1 points on each end to their original positions.
        for idx in (1, iRowLn-2):
            iM1 = idxPts['M',1][idx]
            iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_Coons, iM1)
            ns_WIP.Points.SetPoint(
                iUT_A, iVT_A, ptsM_PreG1[iM1],
                weight=ns_Coons.Points.GetWeight(iUT_A, iVT_A))


    # Average the corner G1 locations.
    for key in pts_G1_corners:
        iM1 = key
        iUT_A, iVT_A = getUvIdxFromNsPoint1dIdx(ns_Coons, iM1)
        if len(pts_G1_corners[key]) == 1:
            ns_WIP.Points.SetPoint(
                iUT_A, iVT_A, pts_G1_corners[key][0],
                weight=ns_Coons.Points.GetWeight(iUT_A, iVT_A))
        elif len(pts_G1_corners[key]) == 2:
            pt_Avg = (pts_G1_corners[key][0] + pts_G1_corners[key][1]) / 2.0
            ns_WIP.Points.SetPoint(
                iUT_A, iVT_A, pt_Avg,
                weight=ns_Coons.Points.GetWeight(iUT_A, iVT_A))
        else:
            raise Exception(
                "{} corner points recorded at one corner.".format(
                    len(pts_G1_corners[key])))

    # TODO: Average the corner G1 locations.



    #pts_WIP = [cp.Location for cp in ns_Coons.Points]

    return ns_WIP


def main():
    
    rc = getInput()
    if rc is None: return

    (
        objrefs_Crvs,
        iContinuity,
        bEcho,
        bDebug,
        ) = rc

    #print objrefs_Edges[0].Geometry()
    #print objrefs_Edges[0].Curve()
    #print objrefs_Edges[0].Edge()
    #print objrefs_Edges[0].Trim()
    #return


    #Rhino.RhinoApp.CommandPrompt = "Working ..."

    #if not bDebug: sc.doc.Views.RedrawEnabled = False


    rc = createSurface(
        rhCrvs_In=objrefs_Crvs,
        iContinuity=iContinuity,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rc is None:
        sc.doc.Views.RedrawEnabled = True
        return
    if isinstance(rc, rg.NurbsSurface):
        ns_Res = rc
    elif isinstance(rc, tuple):
        ns_Res, sLog = rc
        print sLog
        if ns_Res is None:
            sc.doc.Views.RedrawEnabled = True
            return
    else:
        raise ValueError("Bad output: {}".format(rc))

    if ns_Res.IsRational:
        print "Resultant surface is rational.  Check results."

    gB_Res = sc.doc.Objects.AddSurface(ns_Res)
    if gB_Res == gB_Res.Empty:
        print "Brep could not be added to the document."
    #else:
    #    print gB_Res
    sc.doc.Views.Redraw()
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()
    #try: main()
    #except Exception as ex:
    #    print "{0}: {1}".format(type(ex).__name__, "\n".join(ex.args))
