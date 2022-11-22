"""
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
190708: Created.
...
200304: Bug fix in printed output for G0 continuity.
200515: Modified some option defaults.
200814: Bug fix in recognizing G2 and C1 continuities.
200816,17: Bug fix in acknowledging that a colinear joint is G2.
200901: Parameterization is no longer normalized to try for C continuity matches.
200916: Bug fix.
201120,22: Now supports reporting continuities up through C11 and G3 (More if C > 3).
201204: Bug fixes.
210202: Bug fix: Now, when both G3 components are tiny, G3 is returned.
210302: Now G0 continuity check is used for C0 continuity check.
210712: Added more debug printed feedback.
220117: Bug fix when checking for G2.
220121: Bug fix.
220425: Modified printed output.
220705: Modified length tolerance used in C continuity check.

TODO: Add percent tolerance for G2 continuity.  http://catiadoc.free.fr/online/cfyugfss_C2/cfyugfssut_implicitmode_0311.htm
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fDistTol'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fAngleTol_Vector_Deg'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'VectAngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAlignCrvDirs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddDot'; keys.append(key)
    values[key] = False
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
    def setValue(cls, key, idxList=None):

        if key == 'fDistTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

        if key == 'fAngleTol_Vector_Deg':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve and parameter with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves near their ends to evaluate")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.OneByOnePostSelect = True
    go.DisablePreSelect()

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        idxs_Opt['DocTols'] = go.AddOption('DocTols')
        idxs_Opt['TightTols'] = go.AddOption('TightTols')
        addOption('fDistTol')
        addOption('fAngleTol_Vector_Deg')
        addOption('bAlignCrvDirs')
        addOption('bAddDot')
        addOption('bEcho')
        addOption('bDebug')



        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if go.Option().Index == idxs_Opt['DocTols']:
            Opts.riOpts['fDistTol'].CurrentValue =  Opts.riOpts['fDistTol'].InitialValue
            Opts.riOpts['fAngleTol_Vector_Deg'].CurrentValue =  Opts.riOpts['fAngleTol_Vector_Deg'].InitialValue
        elif go.Option().Index == idxs_Opt['TightTols']:
            Opts.riOpts['fDistTol'].CurrentValue =  1e-12
            Opts.riOpts['fAngleTol_Vector_Deg'].CurrentValue =  1e-6
        elif Opts.riOpts['fDistTol'].CurrentValue < 0.0:
            Opts.riOpts['fDistTol'].CurrentValue = Opts.riOpts['fDistTol'].InitialValue
        elif Opts.riOpts['fAngleTol_Vector_Deg'].CurrentValue < 0.0:
            Opts.riOpts['fAngleTol_Vector_Deg'].CurrentValue = Opts.riOpts['fAngleTol_Vector_Deg'].InitialValue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def testWithControlPointLocations():

    fDist_Cps_TanToPos_A = (
            rgCrv_SingleSpan_A.Points[idxCp_Tan_A].Location.DistanceTo(
                    rgCrv_SingleSpan_A.Points[idxCp_Pos_A].Location))

    fDist_Cps_TanToPos_B = (
            rgCrv_SingleSpan_B.Points[idxCp_Tan_B].Location.DistanceTo(
                    rgCrv_SingleSpan_B.Points[idxCp_Pos_B].Location))

    #    sEval='fDist_Cps_TanToPos_B'; print("{}: {}".format(sEval, eval(sEval)))
    #    sEval='fDist_Cps_TanToPos_A'; print("{}: {}".format(sEval, eval(sEval)))
    #    sEval='fDist_Cps_TanToPos_B/fDist_Cps_TanToPos_A'; print("{}: {}".format(sEval, eval(sEval)))

    pt_CrvOppTan_A = rg.Point3d(
            rgCrv_SingleSpan_A.Points[idxCp_Tan_A].Location -
            (   rgCrv_SingleSpan_A.Points[idxCp_Crv_A].Location -
                rgCrv_SingleSpan_A.Points[idxCp_Tan_A].Location)
    )
    #sc.doc.Objects.AddPoint(pt_CrvOppTan_A)
    
    pt_CrvOppTan_B = rg.Point3d(
            rgCrv_SingleSpan_B.Points[idxCp_Tan_B].Location -
            (   rgCrv_SingleSpan_B.Points[idxCp_Crv_B].Location -
                rgCrv_SingleSpan_B.Points[idxCp_Tan_B].Location)
    )
    #sc.doc.Objects.AddPoint(pt_CrvOppTan_B)

    fDist_Cps_CrvToTan_A = (
            rgCrv_SingleSpan_A.Points[idxCp_Crv_A].Location.DistanceTo(
                    rgCrv_SingleSpan_A.Points[idxCp_Tan_A].Location))

    fDist_Cps_CrvToTan_B = (
            rgCrv_SingleSpan_B.Points[idxCp_Crv_B].Location.DistanceTo(
                    rgCrv_SingleSpan_B.Points[idxCp_Tan_B].Location))

    sEval='pt_CrvOppTan_A.DistanceTo(pt_CrvOppTan_B)'; print("{}: {}".format(sEval, eval(sEval)))

    if (
            pt_CrvOppTan_A.DistanceTo(pt_CrvOppTan_B) <= fDistTol
            and
            abs(1.0 - rgCrv_SingleSpan_A.Domain.Length) <= fDistTol
            and
            abs(1.0 - rgCrv_SingleSpan_B.Domain.Length) <= fDistTol
    ):
        sContinuity_MaxFound = "C2"


def formatDistance(fDistance):
    if fDistance is None:
        return "(No deviation provided)"
    elif fDistance < 1e-6:
        return "0"
    elif fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-2)):
        return "{:.1e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def continuityVectorsAt(nc, t, side=rg.CurveEvaluationSide.Default):
    """
    Returns: Tuple of these 4 items:
        3D point as a vector,
        Unit tangent vector,
        Curvature vector,
        Vector for comparing G3 continuity, not the G3 vector itself
        
        None for any of aforementioned vectors on fail.
    """

    if not isinstance(nc, rg.NurbsCurve): return

    vs = nc.DerivativeAt(
        t,
        derivativeCount=3,
        side=side)

    # Not using rg.Curve.TangentAt since it doesn't take into account CurveEvaluationSide.
    vTangency = vs[1]/vs[1].Length
    
    cross = rg.Vector3d.CrossProduct

    # For R3
    vCurvature = (
        cross(cross(vs[1], vs[2]), vs[1])
        /
        vs[1].Length**4
        )

    #sEval='vCurvature_A_Formula'; print("{}: {}".format(sEval, eval(sEval)))
    #sc.doc.Objects.AddLine(rg.Line(start=rg.Point3d(vs[0]), span=vCurvature_A_Formula))

    # For R2
    #vCurvature_A = (
    #            cross(vs[1], vs[2]) /
    #            vs[1].Length**3)


    # From Eq. 3.5 in CAN A CUBIC SPLINE CURVE BE G3
    # in Journal of Computational Mathematics:
    vTorsion = (
        (
        -3.0*(vs[1] * vs[2])*(cross(vs[1], vs[2]))
        /
        (vs[1] * vs[1])**3.0
        )
        +
        cross(vs[1], vs[3]) / (vs[1] * vs[1])**2.0
    )


    return vs[0], vTangency, vCurvature, vTorsion


def getContinuity_C(ncA, tA, crvEvalSide_A, ncB, tB, crvEvalSide_B):
    """
    Returns:
        int(C continuity)
        or
        None
    """

    ncA_WIP = ncA.DuplicateCurve()
    ncB_WIP = ncB.DuplicateCurve()

    if ncA_WIP.Degree < ncB_WIP.Degree:
        ncA_WIP.IncreaseDegree(ncB.Degree)
    elif ncA_WIP.Degree > ncB_WIP.Degree:
        ncB_WIP.IncreaseDegree(ncA.Degree)

    degree_Max = ncA_WIP.Degree


    vsDerivs_A = ncA_WIP.DerivativeAt(
        tA,
        derivativeCount=degree_Max,
        side=crvEvalSide_A)

    if vsDerivs_A[1].IsTiny():
        ncA_WIP.Dispose()
        ncB_WIP.Dispose()
        return


    vsDerivs_B = ncB_WIP.DerivativeAt(
        tB,
        derivativeCount=degree_Max,
        side=crvEvalSide_B)

    if vsDerivs_B[1].IsTiny():
        ncA_WIP.Dispose()
        ncB_WIP.Dispose()
        return

    ncA_WIP.Dispose()
    ncB_WIP.Dispose()

    length_tol = (2**(-52))**0.5 # Rhino.RhinoMath.ZeroTolerance**0.5

    for i in range(len(vsDerivs_A)):
        vA = vsDerivs_A[i]
        vB = vsDerivs_B[i]

        if vA.IsTiny() and vB.IsTiny():
            if i == degree_Max:
                return float('inf')
            continue
        elif vA.IsTiny() or vB.IsTiny():
            return i - 1

        vDelta = vA - vB
        # Don't use Vector3d.IsTiny since its epsilon is too small for this check.

        length = vDelta.Length

        if length > length_tol:
            return i - 1

    #    if not vDelta.IsTiny():
    #        return i - 1

    return float('inf')


def processCurves(rgCrv0_A, rgCrv0_B, bEvalT1End_A, bEvalT1End_B, fDistTol, fAngleTol_Vector_Deg, bAlignCrvDirs=True, bDebug=False):
    """
    Returns tuple of string, string: (strShortContinuityDescription, strLongContinuityDescription)
    """

    nc_A = rgCrv0_A.ToNurbsCurve()
    nc_B = rgCrv0_B.ToNurbsCurve()

    degA = nc_A.Degree
    degB = nc_B.Degree

    if bAlignCrvDirs and bEvalT1End_A == bEvalT1End_B:
        if nc_B.Reverse():
            bEvalT1End_B = not bEvalT1End_B
            if bDebug:
                print("Curve B direction was reversed to match Curve A.")

    if bEvalT1End_A:
        idxSpan_A = nc_A.SpanCount - 1
        domainSpan_A = nc_A.SpanDomain(idxSpan_A)
        nc_SingleSpan_A = nc_A.Trim(domainSpan_A)
        t_A = nc_SingleSpan_A.Domain.T1
        crvEvalSide_A = rg.CurveEvaluationSide.Below
        idxCp_Pos_A = nc_SingleSpan_A.Points.Count - 1
        idxCp_Tan_A = nc_SingleSpan_A.Points.Count - 2
        idxCp_Crv_A = nc_SingleSpan_A.Points.Count - 3
    else:
        idxSpan_A = 0
        domainSpan_A = nc_A.SpanDomain(idxSpan_A)
        nc_SingleSpan_A = nc_A.Trim(domainSpan_A)
        t_A = nc_SingleSpan_A.Domain.T0
        crvEvalSide_A = rg.CurveEvaluationSide.Above
        idxCp_Pos_A = 0
        idxCp_Tan_A = 1
        idxCp_Crv_A = 2

    if bEvalT1End_B:
        # Point is closer to end of curve.
        bEvalT1End_B = True
        idxSpan_B = nc_B.SpanCount - 1
        domainSpan_B = nc_B.SpanDomain(idxSpan_B)
        nc_SingleSpan_B = nc_B.Trim(domainSpan_B)
        t_B = nc_SingleSpan_B.Domain.T1
        crvEvalSide_B = rg.CurveEvaluationSide.Below
        idxCp_Pos_B = nc_SingleSpan_B.Points.Count - 1
        idxCp_Tan_B = nc_SingleSpan_B.Points.Count - 2
        idxCp_Crv_B = nc_SingleSpan_B.Points.Count - 3
    else:
        idxSpan_B = 0
        domainSpan_B = nc_B.SpanDomain(idxSpan_B)
        nc_SingleSpan_B = nc_B.Trim(domainSpan_B)
        t_B = nc_SingleSpan_B.Domain.T0
        crvEvalSide_B = rg.CurveEvaluationSide.Above
        idxCp_Pos_B = 0
        idxCp_Tan_B = 1
        idxCp_Crv_B = 2

    # All data required for continuity evaluation have been obtained from the NurbsCurves.
    nc_A.Dispose()
    nc_B.Dispose()


    vsDerivs_A = nc_SingleSpan_A.DerivativeAt(
        t_A,
        derivativeCount=max(3, degA, degB),
        side=crvEvalSide_A)

    if vsDerivs_A[1].IsTiny():
        s  = "Cannot evaluated curve."
        s += "  Check for stacked control points at end of A."
        return None, None, s


    vsDerivs_B = nc_SingleSpan_B.DerivativeAt(
        t_B,
        derivativeCount=max(3, degA, degB),
        side=crvEvalSide_B)

    if vsDerivs_B[1].IsTiny():
        s  = "Cannot evaluated curve."
        s += "  Check for stacked control points at end of B."
        return None, None, s


    fZeroDistTol = 1e-9


    if bDebug:
        sEval='fZeroDistTol'; print("{}: {}".format(sEval, eval(sEval)))
        for i in range(max(degA, degB)+1):
            print('-'*40)
            sEval='vsDerivs_A[{}]'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='vsDerivs_B[{}]'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='vsDerivs_A[{}].Length'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='vsDerivs_B[{}].Length'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='vsDerivs_A[{0}] / vsDerivs_A[{0}].Length'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='vsDerivs_B[{0}] / vsDerivs_B[{0}].Length'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='vsDerivs_A[{0}] - vsDerivs_B[{0}]'.format(i); print("{}: {}".format(sEval, eval(sEval)))
            sEval='(vsDerivs_A[{0}] - vsDerivs_B[{0}]).Length'.format(i); print("{}: {}".format(sEval, eval(sEval)))


    # 1st derivative.

    vectTanB_toCompareWithA = vsDerivs_B[1]
    if bAlignCrvDirs and bEvalT1End_A == bEvalT1End_B:
        # Flip direction of B because curve directions do not match.
        bReversedTanForA = vectTanB_toCompareWithA.Reverse()

    fAngle_BetweenDerivs_1st_Deg = Rhino.RhinoMath.ToDegrees(
        rg.Vector3d.VectorAngle(vsDerivs_A[1], vectTanB_toCompareWithA))

    fLengthDiff_d1 = abs(vsDerivs_A[1].Length - vsDerivs_B[1].Length)


    # 2nd derivative.

    # Curve direction doesn't affect 2nd derivatives.
    pass

    fAngle_BetweenDerivs_2nd_Deg = rg.Vector3d.VectorAngle(
            vsDerivs_A[2], vsDerivs_B[2])

    if vsDerivs_A[2].IsTiny() or vsDerivs_B[2].IsTiny():
        fLengthDiff_d2 = None
        fLengthRatioDiff_d2 = None
    else:
        
        # Use the following to normalize parameters.  Any practical purpose?
        #multiWith_B_d2 = (rgCrv_SingleSpan_B.Domain.Length**2 /
        #                  rgCrv_SingleSpan_A.Domain.Length**2)
        #sEval='multiWith_B_d2'; print("{}: {}".format(sEval, eval(sEval)))
        #fLengthDiff_d2 = abs(vsDerivs_A[2].Length -
        #                     vsDerivs_B[2].Length * multiWith_B_d2)
        
        
        fLengthDiff_d2 = abs(vsDerivs_A[2].Length - vsDerivs_B[2].Length)


    # 3rd derivative.

    fAngle_BetweenDerivs_3rd_Deg = rg.Vector3d.VectorAngle(
            vsDerivs_A[3], vsDerivs_B[3])

    if vsDerivs_A[3].IsTiny() or vsDerivs_B[3].IsTiny():
        fLengthDiff_d3 = None
        fLengthRatioDiff_d3 = None
    else:
        fLengthDiff_d3 = abs(vsDerivs_A[3].Length - vsDerivs_B[3].Length)


    # Not using TangentAt and CurvatureAt since they do not take evaluation side
    # into consideration.


    (
        vPt_A,
        vUnitTan_A,
        vCurvature_A,
        vTorsion_A
        ) = continuityVectorsAt(
            nc_SingleSpan_A,
            t_A,
            crvEvalSide_A)

    #sEval='vUnitTan_A'; print("{}: {}".format(sEval, eval(sEval)))
    #sc.doc.Objects.AddLine(rg.Line(start=rg.Point3d(vsDerivs_A[0]), span=vUnitTan_A))

    #sEval='vCurvature_A'; print("{}: {}".format(sEval, eval(sEval)))
    #sc.doc.Objects.AddLine(rg.Line(start=rg.Point3d(vG0_A), span=vCurvature_A))

    (
        vPt_B,
        vUnitTan_B,
        vCurvature_B,
        vTorsion_B
        ) = continuityVectorsAt(
            nc_SingleSpan_B,
            t_B,
            crvEvalSide_B)


    if bDebug:
        print('-'*40)
        sEval='vPt_A'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vPt_B'; print("{}: {}".format(sEval, eval(sEval)))
        print('-'*40)
        sEval='vUnitTan_A'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vUnitTan_B'; print("{}: {}".format(sEval, eval(sEval)))
        print('-'*40)
        sEval='vCurvature_A'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vCurvature_B'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vCurvature_A.Length'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vCurvature_B.Length'; print("{}: {}".format(sEval, eval(sEval)))
        if vCurvature_A.Length > 0.0:
            sEval='1.0/vCurvature_A.Length'; print("{}: {}".format(sEval, eval(sEval)))
        if vCurvature_B.Length > 0.0:
            sEval='1.0/vCurvature_B.Length'; print("{}: {}".format(sEval, eval(sEval)))
        print('-'*40)
        sEval='vTorsion_A'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vTorsion_B'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vTorsion_A.Length'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='vTorsion_B.Length'; print("{}: {}".format(sEval, eval(sEval)))

    #sc.doc.Objects.AddLine(rg.Line(start=rg.Point3d(vG0_B), span=vCurvature_B))
    #sc.doc.Views.Redraw(); return

    fAngle_BetweenTans_Deg = Rhino.RhinoMath.ToDegrees(
        rg.Vector3d.VectorAngle(vUnitTan_A, vUnitTan_B))

    fCurvature_A = vCurvature_A.Length
    fCurvature_B = vCurvature_B.Length
    fAngle_BetweenCurvatures_Deg = Rhino.RhinoMath.ToDegrees(
        rg.Vector3d.VectorAngle(vCurvature_A, vCurvature_B))

    fAngle_BetweenTorsions_Deg = Rhino.RhinoMath.ToDegrees(
        rg.Vector3d.VectorAngle(vTorsion_A, vTorsion_B))

    if bDebug:
        print('-'*40)
        sEval='fAngle_BetweenTans_Deg'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='fAngle_BetweenDerivs_1st_Deg'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='fAngle_BetweenDerivs_2nd_Deg'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='fCurvature_A'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='fCurvature_B'; print("{}: {}".format(sEval, eval(sEval)))
        if fCurvature_A > 0.0 and fCurvature_B > 0.0:
            sEval='abs(fCurvature_A-fCurvature_B)/max(fCurvature_A, fCurvature_B)'; print("{}: {}".format(sEval, eval(sEval)))
        sEval='fAngle_BetweenCurvatures_Deg'; print("{}: {}".format(sEval, eval(sEval)))
        #return
    
    ptEval_A = rg.Point3d(vsDerivs_A[0])
    ptEval_B = rg.Point3d(vsDerivs_B[0])
    fDist_BetweenEvalEndPts = ptEval_A.DistanceTo(ptEval_B)


    def getContinuityG():
        """
        Checks up through G3.
        Returns: int(G continuity) or None or not G0
        """

        if fDist_BetweenEvalEndPts > fDistTol:
            # Not G0.
            return

        if vUnitTan_A.Length == 0.0 or vUnitTan_A.Length == 0.0:
            # Stacked-control points?
            return

        if fAngle_BetweenTans_Deg > fAngleTol_Vector_Deg:
            # Not G1.
            return 0

        #        iIsParallelTo_Tan = vUnitTan_A.IsParallelTo(
        #            vUnitTan_B, fAngleTol_Vector_Deg)
        #
        #        if not iIsParallelTo_Tan:
        #            return 0

        #        # When both sides are linear, continuity is G1.
        #        if vCurvature_A.IsTiny():
        #            if vCurvature_B.IsTiny():
        #                return 1


        #if vCurvature_A.Length == 0.0 or vCurvature_B.Length == 0.0:
        #    return 1
        if vCurvature_A.IsTiny() and vCurvature_B.IsTiny():
            # Colinear.
            return 3

        if vCurvature_A.IsTiny() or vCurvature_B.IsTiny():
            # Not G2.
            return 1


        fDelta_Radii = abs(1.0/vCurvature_A.Length - 1.0/vCurvature_B.Length)
        if fDelta_Radii > fDistTol:
            # Not G2.
            return 1

        if fAngle_BetweenCurvatures_Deg > fAngleTol_Vector_Deg:
            # Not G2.
            return 1


        if vTorsion_A.IsTiny() and vTorsion_B.IsTiny():
            # Co-?.
            return 3
        
        if vTorsion_A.IsTiny() or vTorsion_B.IsTiny():
            # Not G3.
            return 2

        vDelta = vTorsion_A - vTorsion_B

        if not vDelta.IsTiny():
            # Not G3.
            return 2

        # This routine does not test for continuities above G3.
        return 3


    iCont_G_MaxFound = getContinuityG()

    iCont_C_MaxFound = getContinuity_C(
        nc_SingleSpan_A,
        t_A,
        crvEvalSide_A,
        nc_SingleSpan_B,
        t_B,
        crvEvalSide_B,
        )

    nc_SingleSpan_A.Dispose()
    nc_SingleSpan_B.Dispose()

    if iCont_G_MaxFound == 3 and iCont_C_MaxFound > iCont_G_MaxFound:
        iCont_G_MaxFound = iCont_C_MaxFound

    sModelUnitSystem = sc.doc.ModelUnitSystem.ToString().lower()
    
    # This is the output format more similar to _GCon.
    #s = "Curve end difference = {:.{}f} {}".format(
    #        fDist_BetweenEvalEndPts,
    #        sc.doc.ModelDistanceDisplayPrecision,
    #        sModelUnitSystem)
    #s += "\nRadius of curvature difference = {:.{}f} {}".format(
    #        abs(1.0/fCurvature_A - 1.0/fCurvature_B),
    #        sc.doc.ModelDistanceDisplayPrecision,
    #        sModelUnitSystem)
    #s += "\nCurvature direction difference in degrees = {:.{}f}".format(
    #        fAngle_BetweenCurvatures_Deg,
    #        sc.doc.ModelDistanceDisplayPrecision)
    #s += "\nTangent difference in degrees = {:.{}f}".format(
    #        fAngle_BetweenDerivs_1st_Deg,
    #        sc.doc.ModelDistanceDisplayPrecision)
    #if iCont_G_MaxFound:
    #    s += "\nCurves are G{}.".format(iCont_G_MaxFound)
    #else:
    #    s += "\nCurve ends are out of tolerance."
    #print(s



    # My take.
    s  = "Torsion difference: "

    if (
        vTorsion_A.Length <= fZeroDistTol and
        vTorsion_B.Length <= fZeroDistTol
    ):
        s += "(Both curves have no change in curvature at evaluated ends.)"
    elif (
        vTorsion_A.Length <= fZeroDistTol or
        vTorsion_B.Length <= fZeroDistTol
    ):
        s += "(One curve has no change in curvature at evaluated end.)"
    else:
        s += "{:.{}f} degrees".format(
                fAngle_BetweenTorsions_Deg,
                sc.doc.ModelDistanceDisplayPrecision)
        s += ", {:.{}f} magnitude".format(
                abs(vTorsion_A.Length - vTorsion_B.Length),
                sc.doc.ModelDistanceDisplayPrecision)

        s += "\n  3rd derivative vector differences:"
    
        s += " {:.{}f} degrees".format(
                fAngle_BetweenDerivs_3rd_Deg,
                sc.doc.ModelDistanceDisplayPrecision)

    if fLengthDiff_d3:
        s += ", {} magnitude".format(
            formatDistance(fLengthDiff_d3))


    s += "\nCurvature difference: "

    if (
        fCurvature_A <= fZeroDistTol and
        fCurvature_B <= fZeroDistTol
    ):
        s += "(Both curves are linear at evaluated ends.)"
    elif (
        fCurvature_A <= fZeroDistTol or
        fCurvature_B <= fZeroDistTol
    ):
        s += "(One curve is linear at evaluated end.)"
    else:
        s += "{:.{}f} degrees".format(
                fAngle_BetweenCurvatures_Deg,
                sc.doc.ModelDistanceDisplayPrecision)
        s += ", R{:.{}f} {}".format(
                abs(1.0/fCurvature_A - 1.0/fCurvature_B),
                sc.doc.ModelDistanceDisplayPrecision,
                sModelUnitSystem)
        s += "\n  2nd derivative vector differences:"
    
        s += " {:.{}f} degrees".format(
                fAngle_BetweenDerivs_2nd_Deg,
                sc.doc.ModelDistanceDisplayPrecision)
    
    if fLengthDiff_d2:
        s += ", {} magnitude".format(
            formatDistance(fLengthDiff_d2))


    s += "\nTangent difference: {:.{}f} degrees".format(
            fAngle_BetweenDerivs_1st_Deg,
            sc.doc.ModelDistanceDisplayPrecision)

    s += "\n  1st derivative vector differences:"

    s += " {:.{}f} degrees".format(
            fAngle_BetweenDerivs_1st_Deg,
            sc.doc.ModelDistanceDisplayPrecision)
    s += ", {} magnitude".format(
        formatDistance(fLengthDiff_d1))

    s += "\nCurve end difference = {:.{}f} {}".format(
            fDist_BetweenEvalEndPts,
            sc.doc.ModelDistanceDisplayPrecision,
            sModelUnitSystem)

    if iCont_G_MaxFound is not None:
        s += "\nContinuities at curves' ends are G{} and C{}.".format(
            iCont_G_MaxFound, iCont_C_MaxFound)
    else:
        s += "\nNo continuity between the curves."

    return (
            iCont_G_MaxFound,
            iCont_C_MaxFound,
            s,
    )


def processCurveObjects(objrefs, fDistTol, fAngleTol_Vector_Deg, bAlignCrvDirs, bAddDot, bEcho, bDebug):
    """
    """
    rgCrv0_A = objrefs[0].Curve()

    bSuccess, tA = rgCrv0_A.ClosestPoint(objrefs[0].SelectionPoint())
    if not bSuccess:
        rgCrv0_A.Dispose()
        return
    
    rgCrv0_B = objrefs[1].Curve()
    if isinstance(rgCrv0_B, rg.NurbsCurve):
        rgNurbsCrv_B = rgCrv0_B.Duplicate()
    else:
        rgNurbsCrv_B = rgCrv0_B.ToNurbsCurve()

    bSuccess, tB = rgCrv0_B.ClosestPoint(objrefs[1].SelectionPoint())
    if not bSuccess:
        rgCrv0_B.Dispose()
        return
    
    bEvalT1End_A = tA >= rgCrv0_A.Domain.Mid
    bEvalT1End_B = tB >= rgCrv0_B.Domain.Mid

    rc = processCurves(
            rgCrv0_A,
            rgCrv0_B,
            bEvalT1End_A,
            bEvalT1End_B,
            fDistTol,
            fAngleTol_Vector_Deg,
            bAlignCrvDirs,
            bDebug,
    )
    if rc is None:
        rgCrv0_A.Dispose()
        rgCrv0_B.Dispose()
        return

    (
        iCont_G_MaxFound,
        iCont_C_MaxFound,
        sContinuityDescr,
    ) = rc

    print(sContinuityDescr)

    if bAddDot and iCont_G_MaxFound:
        rgDot = rg.TextDot(
                text="G{}/C{}".format(iCont_G_MaxFound, iCont_C_MaxFound),
                location=rgCrv0_A.PointAt(t_A))
        rgDot.FontHeight = 11
        sc.doc.Objects.AddTextDot(rgDot)
        sc.doc.Views.Redraw()

    rgCrv0_A.Dispose()
    rgCrv0_B.Dispose()

    if iCont_G_MaxFound:
        return iCont_G_MaxFound, iCont_C_MaxFound


def main():

    objrefs = getInput()
    if objrefs is None: return

    fDistTol = Opts.values['fDistTol']
    fAngleTol_Vector_Deg = Opts.values['fAngleTol_Vector_Deg']
    bAlignCrvDirs = Opts.values['bAlignCrvDirs']
    bAddDot = Opts.values['bAddDot']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    processCurveObjects(
        objrefs,
        Opts.values['fDistTol'],
        Opts.values['fAngleTol_Vector_Deg'],
        Opts.values['bAlignCrvDirs'],
        Opts.values['bAddDot'],
        Opts.values['bEcho'],
        Opts.values['bDebug'],
        )


if __name__ == '__main__': main()