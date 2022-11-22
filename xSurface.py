"""
160618: Created.
...
190617: Using NurbsSurface in place of other surfaces due to bugs in RhinoCommon stated below.
190720: A Cone-related function now check whether test point is at the Apex and returns accordingly.
190722: Split some functions.  Added a function.
190727: Corrected some output and notes.
190803: Modifications for when a value < 1e-12 is passed in as a deviation tolerance to some functions.
190810: getDescription now reports the construction curves for Rev and Sum surfaces. 
191118: Moved some functions to other modules.
191124: Import-related bug fix.
191126: Modified some printed output.
201212: Added knot information to description of NurbsSurface.
210325: getDescription now tests for primitives using more tolerances and reports the primitive's radius, etc.
210327: Bug fix.
210330: Modified some tolerances to test for 'perfect' primitive match.
210604: Bug fix in knot multiplicity routine.
220317: Added splitSurfaceIntoBrep.
220623: Bug fix.

Known Issues with RhinoCommon:
    
    RevSurface.IsPlanar and
    SumSurface.IsPlanar (5 - 6.11+?):
        Some results are erroneously False.  Use NurbsSurface instead (ToNurbsSurface).
    
    RevSurface.TryGetPlane:
        Some erroneous Plane results.  Use NurbsSurface instead (ToNurbsSurface).
    
    RevSurface.TryGetSphere (5 - 6.11+?):
        Resultant sphere IsValid == False.  Use BrepFace or NurbsSurface conversion instead.
    
    RevSurface.TryGetCylinder (5 - 6.14+?):
        Cylinder not found when it should.  Convert to NurbsSurface and analyze that instead.
    
    NurbsSurface.IsSphere(tolerance) and
    NurbsSurface.TryGetSphere(tolerance) (5.? - 6.11+?):
        Sphere not found for some sizes at tolerance < 1e-9.
    
    NurbsSurface.IsTorus(tolerance) and
    NurbsSurface.TryGetTorus(tolerance) (5.? - 6.11+?):
        Torus not found for some sizes at tolerance < 1e-9.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

import System


maxExponentFor2 = 4


def getDescription(rgSrf0):
    """
    Returns: String
    """

    typeSrf = rgSrf0.GetType()

    # Check whether face's UnderlyingSurface is already a primitive shape.
    if typeSrf == rg.PlaneSurface:
        return "PlaneSurface"

    fTols_ToTry = Rhino.RhinoMath.ZeroTolerance, 1e-9, 1e-8
    # 1e-8 can compensate for ZeroTolerance scaling from inches to millimeters
    # because 2.32830643654e-10 * 25.4 ~ 5.91389834881e-09.

    for fTol in fTols_ToTry:
        if rgSrf0.IsPlanar(fTol):
            return "Planar ({}) {}".format(fTol, typeSrf.Name)
        #    elif (
        #            (typeSrf==rg.RevSurface or typeSrf==rg.SumSurface)
        #            and rgSrf0.IsPlanar(fTol)
        #    ):
        #        return "Planar ({}) {}".format(fTol, typeSrf.Name)

    if rgSrf0.IsSolid: s = "Solid"
    elif (rgSrf0.IsClosed(0) or rgSrf0.IsClosed(1)): s = "Closed"
    else: s = "Open"

    for fTol in fTols_ToTry:
        rc = rgSrf0.TryGetCylinder(fTol)
        if rc[0]:
            cyl = rc[1]
            s += ", R{:.4f} cylindrical ({:.3e})".format(cyl.Radius, fTol)
            break
        rc = rgSrf0.TryGetCone(fTol)
        if rc[0]:
            cone = rc[1]
            s += ", {:.2f} degree conical ({:.3e})".format(cone.AngleInDegrees(), fTol)
            break
        rc = rgSrf0.TryGetSphere(fTol)
        if rc[0]:
            sphere = rc[1]
            s += ", R{:.4f} spherical ({:.3e})".format(sphere.Radius, fTol)
            break
        rc = rgSrf0.TryGetTorus(fTol)
        if rc[0]:
            torus = rc[1]
            s += ", R{:.4f},r{:.4f} toric ({:.3e})".format(torus.MajorRadius, torus.MinorRadius, fTol)
            break
    
    s += " {}".format(typeSrf.Name)
    
    if typeSrf==rg.RevSurface:
        s += " with {} revolute".format(rgSrf0.Curve.GetType().Name)
    elif typeSrf==rg.SumSurface:
        s += " composed of {} and {}".format(
                rgSrf0.IsoCurve(0, rgSrf0.Domain(1).T0).GetType().Name,
                rgSrf0.IsoCurve(1, rgSrf0.Domain(0).T0).GetType().Name,
                )
    elif typeSrf == rg.NurbsSurface:
        s += "  Rational," if rgSrf0.IsRational else "  Non-rational,"
        s += "  Degrees {} x {}".format(rgSrf0.Degree(0), rgSrf0.Degree(1))
        s += "  CpCts {} x {}".format(rgSrf0.Points.CountU, rgSrf0.Points.CountV)
        s += "  SpanCts {} x {}".format(rgSrf0.SpanCount(0), rgSrf0.SpanCount(1))
        if Rhino.RhinoApp.ExeVersion >= 7:
            s += "  KnotStyles {} x {}".format(
                rgSrf0.KnotsU.KnotStyle, rgSrf0.KnotsU.KnotStyle)


        iCt_KnotMulties_U = []
        iK = 0
        while iK < rgSrf0.KnotsU.Count:
            m = rgSrf0.KnotsU.KnotMultiplicity(iK)
            iCt_KnotMulties_U.append(m)
            iK += m

        iCt_KnotMulties_V = []
        iK = 0
        while iK < rgSrf0.KnotsV.Count:
            m = rgSrf0.KnotsV.KnotMultiplicity(iK)
            iCt_KnotMulties_V.append(m)
            iK += m
        s += "  KnotMultiplicities {} x {}".format(
            iCt_KnotMulties_U, iCt_KnotMulties_V)


    return s


def isIsoStatusAtSeam(rgSrf, isoStatus):
    if isoStatus == Rhino.Geometry.IsoStatus.South:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Mid, rgSrf.Domain(1).Min)
    elif isoStatus == Rhino.Geometry.IsoStatus.East:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Max, rgSrf.Domain(1).Mid)
    elif isoStatus == Rhino.Geometry.IsoStatus.North:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Mid, rgSrf.Domain(1).Max)
    elif isoStatus == Rhino.Geometry.IsoStatus.West:
        return rgSrf.IsAtSeam(rgSrf.Domain(0).Min, rgSrf.Domain(1).Mid)


def isSrfSideSingularPerIsoStatus(rgSrf, isoStatus):
    if isoStatus == Rhino.Geometry.IsoStatus.South: return rgSrf.IsSingular(0)
    if isoStatus == Rhino.Geometry.IsoStatus.East: return rgSrf.IsSingular(1)
    if isoStatus == Rhino.Geometry.IsoStatus.North: return rgSrf.IsSingular(2)
    if isoStatus == Rhino.Geometry.IsoStatus.West: return rgSrf.IsSingular(3)


def pointsOfSenw(rgSrf, side):
    """
    side can be an Rhino.Geometry.IsoStatus or
    integers 0 (South), 1 (East), 2 (North), & 3 (West).
    """
    if side == 0 or side == Rhino.Geometry.IsoStatus.South:
        ptList = rgSrf.Points
        numU = ptList.CountU
        pts = [ptList.GetControlPoint(u, 0).Location for u in range(numU)]
    elif side == 1 or side == Rhino.Geometry.IsoStatus.East:
        ptList = rgSrf.Points
        numU = ptList.CountU
        numV = ptList.CountV
        pts = [ptList.GetControlPoint(numU-1, v).Location for v in range(numV)]
    elif side == 2 or side == Rhino.Geometry.IsoStatus.North:
        ptList = rgSrf.Points
        numU = ptList.CountU
        numV = ptList.CountV
        pts = [ptList.GetControlPoint(u, numV-1).Location for u in range(numU)]
    elif side == 3 or side == Rhino.Geometry.IsoStatus.West:
        ptList = rgSrf.Points
        numV = ptList.CountV
        pts = [ptList.GetControlPoint(0, v).Location for v in range(numV)]
    else: return
    
    return pts


def senwCurvePerSide(rgSrf, side):
    """
    side can be an Rhino.Geometry.IsoStatus or
    integers 0 (South), 1 (East), 2 (North), & 3 (West).
    Curve direction matches that of surface, which means that
    North and West curves are in opposition to the trim loop.
    """
    if side == 0 or side == Rhino.Geometry.IsoStatus.South:
        return rgSrf.IsoCurve(0, rgSrf.Domain(1).Min)
    elif side == 1 or side == Rhino.Geometry.IsoStatus.East:
        return rgSrf.IsoCurve(1, rgSrf.Domain(0).Max)
    elif side == 2 or side == Rhino.Geometry.IsoStatus.North:
        return rgSrf.IsoCurve(0, rgSrf.Domain(1).Max)
    elif side == 3 or side == Rhino.Geometry.IsoStatus.West:
        return rgSrf.IsoCurve(1, rgSrf.Domain(0).Min)


def shortSenws(rgSrf, fSenwLen_MinAllowed, bEcho=False):
    def maximumPointSpread(pts):
        distMax = 0.
        for i in range(len(pts)):
            ptA = pts[i]
            for j in range(i+1, len(pts)):
                ptB = pts[j]
                dist = ptA.DistanceTo(ptB)
                if dist > distMax: distMax = dist
        return distMax
    
    
    iSenwsWithinPtSpread = []
    for side in (0,1,2,3):
        if not rgSrf.IsSingular(side):
            crv = senwCurvePerSide(rgSrf, side)
            length = crv.GetLength()
            if length < fSenwLen_MinAllowed:
                iSenwsWithinPtSpread.append(side)
            else:
                pts_CP = pointsOfSenw(rgSrf, side)
                if maximumPointSpread(pts_CP) < fSenwLen_MinAllowed:
                    if bEcho: print "SENW found using maximumPointSpread instead of GetLength()."
                    iSenwsWithinPtSpread.append(side)
    return iSenwsWithinPtSpread


def splitSurfaceIntoBrep(rgSrf_toSplit, rgCrvs_Splitters, **kwargs):
    """
    Parameters:
        rgSrf_toSplit: Can be rg.BrepFace or other rg.Surface.
        rgCrvs_Splitters
        fTolerance
        bTryOtherTolsOnFail
        bDebug
    Returns on success:
        rg.Brep
    Returns on fail:
        None
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    bTryOtherTolsOnFail = getOpt('bTryOtherTolsOnFail')
    bDebug = getOpt('bDebug')


    if isinstance(rgSrf_toSplit, rg.BrepFace):
        rgFace_toSplit = rgSrf_toSplit
        rgBrep_TempForUnderlyingSrf = None
    elif isinstance(rgSrf_toSplit, rg.Surface):
        rgBrep_TempForUnderlyingSrf = rgSrf_toSplit.ToBrep()
        rgFace_toSplit = rgBrep_TempForUnderlyingSrf.Faces[0]
    else:
        return


    def getFormattedDistance(fDistance):
        if fDistance is None: return "(No deviation provided)"
        if fDistance < 0.001:
            return "{:.2e}".format(fDistance)
        else:
            return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


    # Create tolerance loop.
    if bTryOtherTolsOnFail:
        fTols_toTry = []
        for tolMultiplier in 1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0, 0.0625, 16.0:
            tol = tolMultiplier * fTolerance
            fTols_toTry.append(tol)
    else:
        fTols_toTry = fTolerance,


    # Split in a tolerance loop.
    for fTol_toTry in fTols_toTry:

        #
        #
        rgB_Split = rgFace_toSplit.Split(
            curves=rgCrvs_Splitters,
            tolerance=fTol_toTry)
        #
        #


        if bDebug: sEval='rgB_Split'; print sEval+':',eval(sEval)
        
        if rgB_Split is None:
            if bDebug:
                print "  Failed at fTol_toTry=={}.".format(
                    getFormattedDistance(fTol_toTry))
                for c in crvs:
                    sc.doc.Objects.AddCurve(c)
                sc.doc.Views.Redraw()
        elif rgB_Split.Faces.Count == rgFace_toSplit.Brep.Faces.Count:
            if bDebug:
                    if rgB_Split.Faces.Count == 1:
                        print "BrepFace.Split resulted in a 1-face Brep."
                    else:
                        print "BrepFace.Split resulted in a Brep with no additional faces."

            rgB_Split = None

            if not isinstance(rgFace_toSplit.UnderlyingSurface(), rg.NurbsSurface):
                s = "{} face was not split.".format(
                    rgFace_toSplit.UnderlyingSurface().GetType().Name)
                s += "  Trying NurbsSurface equivalent ..."
                print s


                def convertToNS():
                    ns = rgFace_toSplit.UnderlyingSurface().ToNurbsSurface()
                    rgB_1F_NS = ns.ToBrep()
                    if not rgB_1F_NS.IsValid:
                        return
                    if rgB_1F_NS.Faces.Count != 1:
                        return
                    return rgB_1F_NS.Faces[0]

                rgF_NS_toSplit = convertToNS()

                if rgF_NS_toSplit is not None:
                    rgB_Split = rgF_NS_toSplit.Split(
                        curves=rgCrvs_Splitters,
                        tolerance=fTol_toTry)
                    if rgB_Split.Faces.Count == rgFace_toSplit.Brep.Faces.Count:
                        rgB_Split = None
                        print "  NurbsSurface face also didn't split."
                    else:
                        print "  NurbsSurface passed.  This means that the modified" \
                            " model has an underlying surface which was converted" \
                            " from {} to a NurbsSurface.".format(
                                rgFace_toSplit.UnderlyingSurface().GetType().Name)
        else:
            if bDebug:
                sEval='rgB_Split.IsValid'; print sEval+':',eval(sEval)
                sEval='rgB_Split.Faces.Count'; print sEval+':',eval(sEval)
                #sc.doc.Objects.AddBrep(rgB_Split); sc.doc.Views.Redraw()
            if not rgB_Split.IsValid:
                rgB_Split = None

            if bDebug or abs(fTol_toTry - fTolerance) > 1e-9:
                print "  Split successful at a tolerance of {}.".format(
                    getFormattedDistance(fTol_toTry))
            break # out of tolerance loop.


    if rgBrep_TempForUnderlyingSrf: rgBrep_TempForUnderlyingSrf.Dispose()

    return rgB_Split


