"""
V7 includes knot type reporting, so the isUniform function is more useful for V5 & V6.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
190410: Created from another script.
190510: Fixed bug in testing for invalid knot vector.
210109: Now reports knot info when this script is run on its own.
210417: Now processes SumSurface when at least one of its 2 defining curves is a NurbsCurve.
210514: Fixed bug that was reporting single-spanned degree 1 as uniform.
210714: Now supports analyzing PolyCurves with a single NurbsCurve segment.
211001: Added SpanCount to printed output and fixed a printed output bug.
211010: Full knot list is now included in report.
211026: Modified format of printed output of knot parameters.
211029: Now, the underlying curve of single-spanned PolyCurves are reported.
211102: Added NurbsCurve(/Surface)KnotList.KnotStyle.
211126-27: Bug fix in incorrect reporting of knot indices per knot parameter.
        Now processes the NurbsCurve equivalent of a multi-span PolyCurve.
211208: Minimum and maximum knot deltas are now reported.
211228: getKnotInfo now reports degree for periodic curves and surfaces.
220327: Added knot count to printed output.


TODO: Test for invalid knot vector may need more work.
"""

import Rhino
import Rhino.Geometry as rg


def isUniform(knots):
    
    knots = list(knots)
    
    len_knots = len(knots)
    
    if len_knots == 2:
        return True
    elif len_knots < 2:
        return
    
    bPeriodicOrDeg1 = not Rhino.RhinoMath.EpsilonEquals(
            knots[0], knots[1], Rhino.RhinoMath.ZeroTolerance)
    
    if not bPeriodicOrDeg1 and len_knots > 2:
        if not Rhino.RhinoMath.EpsilonEquals(
                knots[-1], knots[-2], Rhino.RhinoMath.ZeroTolerance
        ):
            print("Knot vector is invalid: {}".format(knots))
            return
    
    knotVect1stSpan = knots[1] - knots[0]
    
    if not bPeriodicOrDeg1:
        for i in xrange(1, len_knots//2):
            if not Rhino.RhinoMath.EpsilonEquals(knots[i] - knots[i+1],
                    knotVect1stSpan,
                    Rhino.RhinoMath.ZeroTolerance):
                break
        else:
            print("Degree cannot be determined.")
            return
        
        degree = i + 1
    else:
        degree = 1 # Even if periodic of degree > 1.
    
    # Set up for loop's range's start and stop based on whether knot vector is periodic.
    # iKnot_Start: Start of spans to check.
    # iKnot_End: End of spans to check.
    if bPeriodicOrDeg1:
        iKnot_Start = 0
        iKnot_End = len_knots - 1
        knotVect1stSpan_Not0 = knotVect1stSpan
    else:
        iKnot_Start = degree - 1
        iKnot_End = len_knots - degree
        if iKnot_Start == iKnot_End: return True # Is this also true when periodic?
        knotVect1stSpan_Not0 = knots[degree] - knots[degree-1]
    pass
    #  A uniform knot vector's interior parameter spans are all equal.
    for i in range(iKnot_Start, iKnot_End):
        span = knots[i+1] - knots[i]
        if not Rhino.RhinoMath.EpsilonEquals(
                span,
                knotVect1stSpan_Not0,
                Rhino.RhinoMath.ZeroTolerance):
            return False
    
    return True


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


def getKnotInfo(knots, degree=None):
    """
    Parameters:
        knots: iterable of knot parameters
        degree: int  (This is only used for periodic curves and surfaces.)
    Returns a string.
    """

    iK = 0
    ts_Unique = []
    ms = []
    iKsPerM = []
    deltas_ts = []

    while iK < knots.Count:
        k = knots[iK]
        ts_Unique.append(k)
        m = knots.KnotMultiplicity(index=iK)
        ms.append(m)
        #print("{} at {:.4f}".format(m, k),
        iKsPerM.append(range(iK, iK+m))
        if iK > 0: deltas_ts.append(ts_Unique[-1]-ts_Unique[-2])
        iK += m

    s = "Domain:[{: .6f},{: .6f}]".format(ts_Unique[0], ts_Unique[-1])
    s += "  KnotCt:{}".format(knots.Count)
    s += "  SpanCt:{}".format(len(ms)-1)
    s += "  Mults({})".format(",".join(str(i) for i in ms))

    if len(ms[1:-1]) == 0:
        # For Bezier / single span.
        pass
    elif all(i == 1 for i in ms):
        # Periodic.
        if degree:
            s += "  Deg:{}".format(degree)
        s += "  {} uniform".format("Is" if isUniform(knots) else "Not")
    elif all(i == 1 for i in ms[1:-1]):
        s += "  {} uniform".format("Is" if isUniform(knots) else "Not")
    else:
        # For other non-uniform obvious by its list of multiplicies.
        pass

    if Rhino.RhinoApp.ExeVersion >= 7:
        s += "  .KnotStyle:{}".format(knots.KnotStyle)

    # Record knot parameters to at least 3 decimal places but more if
    # required to display all unique values.
    for dec_plcs in range(15, 3-1, -1):
        LIST = [" {: .{}f}".format(t, dec_plcs) for t in ts_Unique]
        if len(set(LIST)) != len(ts_Unique):
            dec_plcs += 1
            break

    s += "\n  "

    zipped = zip(iKsPerM, ts_Unique)
    s_JoinUs = []
    for idxs, t in zipped:
        if len(idxs) == 1:
            s_JoinUs.append("{: .{}f}[{}]".format(t, dec_plcs, idxs[0]))
        else:
            s_JoinUs.append("{: .{}f}[{},{}]".format(t, dec_plcs, idxs[0], idxs[-1]))
    s += " ".join(s_JoinUs)

    if abs(min(deltas_ts) - max(deltas_ts)) <= 10.0**(-dec_plcs):
        s += "  Deltas:{0:.{1}f}".format(
            deltas_ts[0], dec_plcs)
    else:
        s += "  DeltaRange:[{0:.{2}f},{1:.{2}f}]".format(
            min(deltas_ts), max(deltas_ts), dec_plcs)

    return s


def main():
    import rhinoscriptsyntax as rs


    gObj = rs.GetObject(
        "Select curve or face",
        filter=rs.filter.curve + rs.filter.surface,
        preselect=True,
        subobjects=True)
    if gObj is None: return


    rgObj = rs.coercegeometry(gObj)
    rgCrv = rs.coercecurve(gObj)

    fTol_Shape = 1e-8

    if rgCrv is None:
        # Report surface.

        rgF = rs.coercesurface(gObj)
        rgS = rgF.UnderlyingSurface()
        if isinstance(rgS, rg.NurbsSurface):
            ns = rgS
        elif isinstance(rgS, rg.PlaneSurface):
            print("Surface is a PlaneSurface.")
            return
        elif rgS.IsPlanar(fTol_Shape):
            print("Surface is a planar surface.")
            ns = rgS.ToNurbsSurface()
        elif isinstance(rgS, rg.RevSurface):
            rgCrv = rgS.Curve
            if not isinstance(rgCrv, rg.NurbsCurve):
                print("Surface is a RevSurface with a {} revolute curve.".format(
                    rgCrv.GetType().Name))
                return
            ns = rgS.ToNurbsSurface()
        elif isinstance(rgS, rg.SumSurface):
            rgC_South = rgS.IsoCurve(0, rgS.Domain(1).Min)
            print("SumSurface")
            if isinstance(rgC_South, rg.PolyCurve):
                rgC_WIP = rgC_South.CleanUp()
                if rgC_WIP.Domain.EpsilonEquals(rgC_South.Domain, epsilon=1e-9):
                    rgC_South = rgC_WIP
            if isinstance(rgC_South, rg.NurbsCurve):
                print("Curve U  Nurbs  Knot {}".format(
                    getKnotInfo(rgC_South.Knots, rgC_South.Degree)))
            else:
                print("Curve U is a {}.".format(rgC_South.GetType().Name))

            rgC_West = rgS.IsoCurve(1, rgS.Domain(0).Min)
            if isinstance(rgC_West, rg.PolyCurve):
                rgC_WIP = rgC_West.CleanUp()
                if rgC_WIP.Domain.EpsilonEquals(rgC_West.Domain, epsilon=1e-9):
                    rgC_West = rgC_WIP
            if isinstance(rgC_West, rg.NurbsCurve):
                print("Curve V  Nurbs  Knot {}".format(
                    getKnotInfo(rgC_West.Knots, rgC_South.Degree)))
            else:
                print("Curve V is a {}.".format(rgC_West.GetType().Name))

            if rgS.IsCone(fTol_Shape):
                print("Surface is a cone-shaped SumSurface.")
                return
            elif rgS.IsCylinder(fTol_Shape):
                print("Surface is a cylindrical SumSurface.")
                return
            elif rgS.IsSphere(fTol_Shape):
                print("Surface is a spherical SumSurface.")
                return
            elif rgS.IsPlanar(fTol_Shape):
                print("Surface is a planar SumSurface.")
                return
            return
        else:
            print("{} is not yet supported.".format(rgS.GetType().Name))
            return
        print("KnotsU: {}".format(getKnotInfo(ns.KnotsU, ns.Degree(0))))
        print("KnotsV: {}".format(getKnotInfo(ns.KnotsV, ns.Degree(1))))
    else:
        # Report curve.

        if isinstance(rgCrv, rg.BrepEdge):
            rgE = rgCrv
            rgCrv = rgCrv.EdgeCurve
        else:
            rgE = None

        if isinstance(rgCrv, rg.PolyCurve):
            pc = rgCrv.DuplicateCurve()
            pc.RemoveNesting()
            if pc.SegmentCount > 1:
                print("{} is a multi-segment polycurve.".format(
                    "Edge's curve" if rgE else "Curve"),
                      "NurbsCurve equivalent is analyzed.",
                      sep="  ")
                rgCrv = rgCrv.ToNurbsCurve() # Overwrote variable.
            else:
                print("Segment from 1-segment polycurve retrieved ...")
                rgCrv = pc.SegmentCurve(0) # Overwrote variable.

        if not isinstance(rgCrv, rg.NurbsCurve):
            print("{} is a {}.".format(
                "Edge's curve" if rgE else "Curve",
                rgCrv.GetType().Name))
            return

        nc = rgCrv
        knots = nc.Knots
        print(getKnotInfo(knots, nc.Degree))


if __name__ == '__main__': main()