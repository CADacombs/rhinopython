"""
This script does not use Curve.GetNextDiscontinuity.
This is mainly because the offered tolerance type for G2-discontinuity
is not satifactory to the script author's own design work.


At no angular difference, control point distance tolerance is about 1e-8.
At no control point distance difference, angular tolerance is about 2.0e-8 radians.
_Rotate only works at a minimum of ~1e-5 degrees.

Important if Curve.GetNextDiscontinuity is used:
Curve.GetNextDiscontinuity reports false positives for C1 and C2 probably due to
differences in parameterization intervals.  To avoid this situation,
this script converts all curves to NurbsCurves.

See page 5 in http://graphics.stanford.edu/courses/cs348a-17-winter/ReaderNotes/handout27.pdf

Per studies it was found that
Rhino's default G2 discontinuity is where at least one of the following two items is true:
    1. abs(k0-k1) / max(k0,k1) > 0.05, where k0 and k1 are magnitudes of the curvature vectors below and above the knot.
    2. The difference in curvature vectors is > 2.0 degrees.

In comparison, CATIA uses defaults of:
    G0: 0.001 mm
    G1: 0.01 degree
    G2: 0.02 for ratio of curvature |C1 -  C2| / max(C1, C2)  (Is the angle between the curvature vectors not considered?)
http://catiadoc.free.fr/online/cfyugfss_C2/cfyugfssut_implicitmode_0311.htm
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
180905-10: Created.
...
220320: Finished implementing my own discontinuity evaluation instead of Curve.GetNextDiscontinuity.
        G2 now shares the angle tolerance with G1.
        The absolute radius difference is now evaluated as the other component for G2.
220401: Now checks whether both sides are tiny, and if so, passes G2 continuity.
220717: Now, closed, but not periodic, curves are checked at start/end.
220913: Now, modified option display.
221027: Bug fix in output of edge curve split.
230407: Modified behavior when tolerance command options are modified.
230505: Numeric option input is now more robust.
230818: Now, options can be set when curves are preselected.
240106: Bug fix.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

#from System import Math
from System import Guid


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bG2_NotG1'; keys.append(key)
    values[key] = True
    names[key] = 'G'
    riOpts[key] = ri.Custom.OptionToggle(values[key], '1', '2')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fAngleTol_Deg'; keys.append(key)
    values[key] = sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fRadiusTol'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAddPts'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitCurve'; keys.append(key)
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

        if key == 'fAngleTol_Deg':
            if cls.riOpts[key].CurrentValue < 0.0:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                return

        if key == 'fRadiusTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getPreselected(objectType=rd.ObjectType.Curve):
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    if gObjs_Preselected:
        gCrvs_Preselected = []
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
            if rdRhinoObject.Id in gObjs_Preselected:
                if rdRhinoObject.ObjectType == objectType:
                    gCrvs_Preselected.append(rdRhinoObject.Id)
        return tuple(gCrvs_Preselected)


def getInput():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")
    
    go.GeometryFilter = rd.ObjectType.Curve
    
    go.AcceptNumber(True, acceptZero=True)
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    bPreselectedObjsChecked = False
    
    idxs_Opt = {}
    
    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bG2_NotG1')
        addOption('fAngleTol_Deg')
        if Opts.values['bG2_NotG1']:
            addOption('fRadiusTol')
        addOption('bAddPts')
        addOption('bSplitCurve')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            num = go.Number()
            if num == 1.0:
                key = 'bG2_NotG1'
                Opts.riOpts[key].CurrentValue = False
            elif num == 2.0:
                key = 'bG2_NotG1'
                Opts.riOpts[key].CurrentValue = True
            else:
                key = 'fAngleTol_Deg'
                Opts.riOpts[key].CurrentValue = num
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createSpanCurvesEachSideOfParameter(rgCrv0, t):

    ts_SpanDomainEnds = [rgCrv0.SpanDomain(i).T1 for i in xrange(rgCrv0.SpanCount)]

    for idxSpan in xrange(rgCrv0.SpanCount):
        if abs(t-ts_SpanDomainEnds[idxSpan]) <= Rhino.RhinoMath.ZeroTolerance:
            idxSpan_Below = idxSpan
            break
    else:
        print("{} not at span boundary?  {}".format(t, ts_SpanDomainEnds))
        return

    crv_Below = rgCrv0.Trim(rgCrv0.SpanDomain(idxSpan_Below))
    crv_Above = rgCrv0.Trim(rgCrv0.SpanDomain((idxSpan_Below + 1) % rgCrv0.SpanCount))
    
    return crv_Below, crv_Above


def getDiscontinuities_OLD(rgCrv0, sContinuity, fAngleTol_Deg, fRadiusTol=None, bDebug=False):
    """
    Returns: list(float(parameters of discontinuities))
    """

    fAngleTol_Rad = Rhino.RhinoMath.ToRadians(fAngleTol_Deg)
    if bDebug: print("fAngleTol_Rad: {}".format(fAngleTol_Rad))

    t0 = rgCrv0.Domain.Min
    t1 = rgCrv0.Domain.Max
    
    if sContinuity not in sContinuities:
        print("Continuity type, {}, not supported.".format(sContinuity))
    
    if rgCrv0.IsClosed:
        if sContinuity == 'G1' or sContinuity == 'C1':
            continuityType = rg.Continuity.C1_locus_continuous
        elif sContinuity == 'G2':
            continuityType = rg.Continuity.G2_locus_continuous
        elif sContinuity == 'C2':
            continuityType = rg.Continuity.C2_locus_continuous
    else:
        if sContinuity == 'G1' or sContinuity == 'C1':
            continuityType = rg.Continuity.C1_continuous
        elif sContinuity == 'G2':
            continuityType = rg.Continuity.G2_continuous
        elif sContinuity == 'C2':
            continuityType = rg.Continuity.C2_continuous

    ts_discontinuities = []

    get_next = True
    while get_next:
        sc.escape_test()
        get_next, t = rgCrv0.GetNextDiscontinuity(continuityType, t0, t1)
        if get_next:
            
            rc = createSpanCurvesEachSideOfParameter(rgCrv0, t)
            if not rc: return []
            
            rgCrv_SingleSpan_Below, rgCrv_SingleSpan_Above = rc
            
            rgCrv_SingleSpan_Below.Domain = rg.Interval(0.0,1.0)
            rgCrv_SingleSpan_Above.Domain = rg.Interval(0.0,1.0)

            vect1stDeriv_Below = rgCrv_SingleSpan_Below.DerivativeAt(
                    rgCrv_SingleSpan_Below.Domain.T1,
                    derivativeCount=1,
                    side=rg.CurveEvaluationSide.Below)[1]
            vect1stDeriv_Above = rgCrv_SingleSpan_Above.DerivativeAt(
                    rgCrv_SingleSpan_Above.Domain.T0,
                    derivativeCount=1,
                    side=rg.CurveEvaluationSide.Above)[1]
            if bDebug:
                angle = rg.Vector3d.VectorAngle(
                        vect1stDeriv_Below, vect1stDeriv_Above)
                print("Angle difference: {} radians [{} degrees]".format(
                        angle, Rhino.RhinoMath.ToDegrees(angle)))
            if not vect1stDeriv_Below.IsParallelTo(
                    vect1stDeriv_Above, fAngleTol_Rad
            ):
                if bDebug:
                    sEval='vect1stDeriv_Below.IsParallelTo(vect1stDeriv_Above, fAngleTol_Rad)'; print(sEval+': ',eval(sEval))
                ts_discontinuities.append(t)
            else:
                if sContinuity == 'C1' or sContinuity == 'C2':
                    fVector_magnitude_diff = abs(
                            vect1stDeriv_Below.Length -
                            vect1stDeriv_Above.Length)
                    if bDebug:
                        print("Vector magnitude difference: {}".format(
                                fVector_magnitude_diff))
                    if fVector_magnitude_diff > fRadiusTol:
                        ts_discontinuities.append(t)
                    else:
                        if sContinuity == 'C2':
                            ts_discontinuities.append(t)
                elif sContinuity == 'G2':
                        ts_discontinuities.append(t)
            t0 = t # Advance to the next parameter.
    
    return ts_discontinuities


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
        vs[1].Length**4)

    #sEval='vCurvature_A_Formula'; print(sEval+': ',eval(sEval)
    #sc.doc.Objects.AddLine(rg.Line(start=rg.Point3d(vs[0]), span=vCurvature_A_Formula))

    # For R2
    #vCurvature_A = (
    #            cross(vs[1], vs[2]) /
    #            vs[1].Length**3)


    # From Eq. 3.5 in CAN A CUBIC SPLINE CURVE BE G3
    # in Journal of Computational Mathematics:
    vG3_Condition = (
        (
            -3.0*(vs[1] * vs[2])*(cross(vs[1], vs[2]))
            /
            (vs[1] * vs[1])**3.0
        )
        +
        cross(vs[1], vs[3]) / (vs[1] * vs[1])**2.0
        )
    
    return vs[0], vTangency, vCurvature, vG3_Condition


def getDiscontinuities(rgCrv_In, bG2_NotG1, fAngleTol_Deg, fRadiusTol, bDebug=False):
    """
    """

    fAngleTol_Rad = Rhino.RhinoMath.ToRadians(fAngleTol_Deg)
    if bDebug: print("fAngleTol_Rad: {}".format(fAngleTol_Rad))

    nc = rgCrv_In.ToNurbsCurve()

    t0 = nc.Domain.T0
    t1 = nc.Domain.T1

    ts_discontinuities = []

    bG3_discontinuousFound = False

    # This will also skip simple knot overlaps of Periodic curves.
    if nc.IsClosed and not nc.IsPeriodic:
        iK = 0
    else:
        iK = nc.Degree
    iK_Stop = nc.Knots.Count - nc.Degree

    while iK < iK_Stop:
        sc.escape_test()

        m = nc.Knots.KnotMultiplicity(iK)

        if m <= nc.Degree - 3:
            # Continuity at this knot is at least G3.
            iK += m
            continue

        if bDebug:
            print('-'*20)
            sEval='nc.Knots[iK]'; print("{}: {}".format(sEval, eval(sEval)))


        if iK == 0:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                vG3_Below,
                ) = continuityVectorsAt(
                    nc,
                    nc.Knots[nc.Knots.Count - 1],
                    rg.CurveEvaluationSide.Below)
        else:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                vG3_Below,
                ) = continuityVectorsAt(
                    nc,
                    nc.Knots[iK],
                    rg.CurveEvaluationSide.Below)


        (
            vG0_Above,
            vG1_Above,
            vG2_Above,
            vG3_Above,
            ) = continuityVectorsAt(
                nc,
                nc.Knots[iK],
                rg.CurveEvaluationSide.Above)


        if bDebug:
            sEval='vG0_Below'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG0_Above'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG1_Below'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG1_Above'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG2_Below'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG2_Below.IsTiny()'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG2_Above'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG2_Above.IsTiny()'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG3_Below'; print("{}: {}".format(sEval, eval(sEval)))
            sEval='vG3_Above'; print("{}: {}".format(sEval, eval(sEval)))


        # G1.

        if m > nc.Degree - 1:

            iParallel = vG1_Below.IsParallelTo(
                other=vG1_Above, angleTolerance=fAngleTol_Rad)

            fTanDelta_Degrees = Rhino.RhinoMath.ToDegrees(
                rg.Vector3d.VectorAngle(vG1_Below, vG1_Above))

            if iParallel == 1:
                if bDebug:
                    print("Are G1 with tangent vector difference of {:.2f} degrees.".format(
                        fTanDelta_Degrees))
            else:
                print("Not G1 at {} per tangent vector difference of {:.2f} degrees.".format(
                    nc.Knots[iK],
                    fTanDelta_Degrees))
                ts_discontinuities.append(nc.Knots[iK])
                iK += m
                continue

        if not bG2_NotG1:
            # If only checking G1.
            iK += m
            continue

        #sc.doc.Objects.AddLine(rg.Line(rg.Point3d(vG0_Below), vG1_Below))
        #sc.doc.Objects.AddLine(rg.Line(rg.Point3d(vG0_Above), vG1_Above))



        # G2.

        if m > nc.Degree - 2:

            if vG2_Below.IsTiny() and vG2_Above.IsTiny():
                # Linear.
                pass
            else:

                iParallel = vG2_Below.IsParallelTo(
                    other=vG2_Above, angleTolerance=fAngleTol_Rad)

                if bDebug: sEval='iParallel'; print("{}: {}".format(sEval, eval(sEval)))

                fCrvDelta_Degrees = Rhino.RhinoMath.ToDegrees(
                    rg.Vector3d.VectorAngle(vG2_Below, vG2_Above))

                if iParallel != 1:
                    if bDebug:
                        print("Not G2 at {} per curvature vector difference of {:.2f} degrees.".format(
                            nc.Knots[iK],
                            fCrvDelta_Degrees))
                    ts_discontinuities.append(nc.Knots[iK])
                    iK += m
                    continue
                else:
                    kappa_Below = vG2_Below.Length
                    kappa_Above = vG2_Above.Length
                    ratio_of_curvature = (abs(kappa_Below-kappa_Above) /
                                        max(kappa_Below, kappa_Above))
                    if bDebug:
                        print("Ratio of curvature: {:.2f} % ".format(100.0*ratio_of_curvature))
                    delta_radius = abs(1.0/kappa_Below-1.0/kappa_Above)
                    #if ratio_of_curvature > 0.02:
                    if delta_radius > fRadiusTol:
                        #print("Not G2 at {} per relative difference of curvature vector magnitudes being {:.2f} %.  (2% is the tolerance.)".format(
                        #    nc.Knots[iK],
                        #    100.0*ratio_of_curvature))
                        if bDebug:
                            print("{:.{}f} radius difference.".format(
                                delta_radius, sc.doc.ModelDistanceDisplayPrecision))
                        ts_discontinuities.append(nc.Knots[iK])
                        iK += m
                        continue
                    else:
                        if bDebug: print("Are G2.")

                    #sc.doc.Objects.AddLine(rg.Line(rg.Point3d(vG0_Below), vG2_Below))
                    #sc.doc.Objects.AddLine(rg.Line(rg.Point3d(vG0_Above), vG2_Above))



        # G3.

        if bG3_discontinuousFound:
            iK += m
            continue

        iParallel = vG3_Below.IsParallelTo(
            other=vG3_Above, angleTolerance=Rhino.RhinoMath.ToRadians(1.0))
        fG3Delta_Degrees = Rhino.RhinoMath.ToDegrees(
            rg.Vector3d.VectorAngle(vG3_Below, vG3_Above))
        if iParallel == 1:
            if bDebug:
                print("Are G3 with G3 component vector difference of {:.2f} degrees.".format(
                    fG3Delta_Degrees))
        else:
            bG3_discontinuousFound = True
            if bDebug:
                print("Not G3 at {} per G3 component vector difference of {:.2f} degrees.".format(
                    nc.Knots[iK],
                    fG3Delta_Degrees))
    #            sc.doc.Objects.AddPoint(rg.Point3d(vG0_Below))
    #            #ts_discontinuities.append(nc.Knots[iK])


        iK += m

    nc.Dispose()

    return ts_discontinuities


    #get_next = True
    #while get_next:
    #    sc.escape_test()

    #    if Rhino.RhinoApp.ExeVersion >= 7:
    #        get_next, t = nc.GetNextDiscontinuity(
    #            continuityType=continuityType,
    #            t0=t0,
    #            t1=t1,
    #            cosAngleTolerance=Math.Cos(fAngleTol_Rad),
    #            curvatureTolerance=sc.doc.ModelAbsoluteTolerance)
    #    else:
    #        get_next, t = nc.GetNextDiscontinuity(continuityType, t0, t1)

    #    if get_next:
    #        if bDebug:
    #            print("GetNextDiscontinuity returned t of {}".format(t))
            
            
    #        rc = createSpanCurvesEachSideOfParameter(rgCrv_In, t)
    #        if not rc: return []
            
    #        rgCrv_SingleSpan_Below, rgCrv_SingleSpan_Above = rc
            
            
    #        vectDerivs_Below = nc.DerivativeAt(
    #                t,
    #                derivativeCount=2,
    #                side=rg.CurveEvaluationSide.Below)
            
    #        if t == nc.Domain.T1:
    #            vectDerivs_Above = nc.DerivativeAt(
    #                    nc.Domain.T0,
    #                    derivativeCount=2,
    #                    side=rg.CurveEvaluationSide.Above)
    #        else:
    #            vectDerivs_Above = nc.DerivativeAt(
    #                    t,
    #                    derivativeCount=2,
    #                    side=rg.CurveEvaluationSide.Above)
            
    #        if bDebug:
    #            angle = rg.Vector3d.VectorAngle(
    #                    vectDerivs_Below[1], vectDerivs_Above[1])
    #            print("Angle difference: {} radians [{} degrees]".format(
    #                    angle, Rhino.RhinoMath.ToDegrees(angle)))
    #        if not vectDerivs_Below[1].IsParallelTo(
    #                vectDerivs_Above[1],
    #                fAngleTol_Rad
    #        ):
    #            if bDebug:
    #                sEval='vectDerivs_Below[1].IsParallelTo(vectDerivs_Above[1], fAngleTol_Rad)'; print(sEval+': ',eval(sEval))
    #            ts_discontinuities.append(t)
    #        else:
    #            if sContinuity == 'C1' or sContinuity == 'C2':
    #                if vectDerivs_Below[2].IsTiny():
    #                    pass
    #                elif vectDerivs_Above[2].IsTiny():
    #                    pass
    #                else:
    #                    multiWith_Above_d1 = (rgCrv_SingleSpan_Above.Domain.Length /
    #                                          rgCrv_SingleSpan_Below.Domain.Length)
    #                    fVector_magnitude_diff = abs(
    #                            vectDerivs_Below[1].Length -
    #                            vectDerivs_Above[1].Length * multiWith_Above_d1)
    #                    if bDebug:
    #                        print("Vector magnitude difference: {}".format(
    #                                fVector_magnitude_diff))
    #                    if fVector_magnitude_diff > fRadiusTol:
    #                        ts_discontinuities.append(t)
    #                    else:
    #                        if sContinuity == 'C2':
    #                            ts_discontinuities.append(t)
    #            elif sContinuity == 'G2':
    #                    ts_discontinuities.append(t)
    #        t0 = t # Advance to the next parameter.

    #    iK += m

    #nc.Dispose()

    #return ts_discontinuities


def processCurves(curvesAndEdges0, **kwargs):
    #, sContinuity=None, fAngleTol_Deg=None, fRadiusTol=None, bG2=None, bC2=None, bAddPts=None, bSplitCurve=None, bEcho=None, bDebug=None):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bG2_NotG1 = getOpt('bG2_NotG1')
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    fRadiusTol = getOpt('fRadiusTol')
    bAddPts = getOpt('bAddPts')
    bSplitCurve = getOpt('bSplitCurve')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    gCrvs_with_discont = []
    i_Edges_with_discont_ct = 0
    ts_discontinuities_All = []
    gPts = []
    gCrvs_Split = []

    for curveOrEdge0 in curvesAndEdges0:
        type_rdCurve_or_rgEdge = curveOrEdge0.GetType()

        gCrv_In = rs.coerceguid(curveOrEdge0)
        rdObj_In = rs.coercerhinoobject(curveOrEdge0)
        rgCrv_In = rs.coercecurve(curveOrEdge0)

        if type_rdCurve_or_rgEdge == Guid:
            gCrv_In = curveOrEdge0
            rgCrv_In = rs.coercecurve(gCrv_In)
            s = str(curveOrEdge0)
        elif type_rdCurve_or_rgEdge == rg.BrepEdge:
            rgCrv_In = curveOrEdge0.DuplicateCurve()
            s = "Edge"
        if not isinstance(rgCrv_In, rg.Curve):
            if bEcho: print("Not CurveObject or BrepEdge.")
            continue
    
        ts_discontinuities = getDiscontinuities(
                rgCrv_In=rgCrv_In,
                bG2_NotG1=bG2_NotG1,
                fAngleTol_Deg=fAngleTol_Deg,
                fRadiusTol=fRadiusTol,
                bDebug=bDebug)

        if ts_discontinuities:
            ts_discontinuities_All.extend(ts_discontinuities)
            if gCrv_In:
                gCrvs_with_discont.append(gCrv_In)
            else:
                i_Edges_with_discont_ct += 1

        # Testing.
        #for t in ts_discontinuities:
        #    vect_Below = rgCrv0.DerivativeAt(t, derivativeCount=1, side=rg.CurveEvaluationSide.Below)[1]
        #    vect_Above = rgCrv0.DerivativeAt(t, derivativeCount=1, side=rg.CurveEvaluationSide.Above)[1]
        #    print("{}".format(vect_Below)
        #    print("{}".format(vect_Above)
        #    sEval='vect_Below.EpsilonEquals(vect_Above, 1e-12)'; print(sEval+': ',eval(sEval)
        #    sEval='vect_Below.IsParallelTo(vect_Above, Rhino.RhinoMath.ToRadians(1.349934664e-6))'; print(sEval+': ',eval(sEval)
        #    sEval='vect_Below.IsParallelTo(vect_Above, Rhino.RhinoMath.ToRadians(1.349934665e-6))'; print(sEval+': ',eval(sEval)
            

        if bAddPts:
            points = []
            for t in ts_discontinuities:
                pt = rgCrv_In.PointAt(t)
                if not pt: continue
                gPt = sc.doc.Objects.AddPoint(pt)
                if gPt != Guid.Empty:
                    gPts.append(gPt)
        
        if bSplitCurve:
            rgCrvs_Split = rgCrv_In.Split(ts_discontinuities)
            if rgCrvs_Split:
                attr = rdObj_In.Attributes if rdObj_In else None
                bAddFail = False
                for rgCrv_Split in rgCrvs_Split:
                    gCrv_Split = sc.doc.Objects.AddCurve(rgCrv_Split, attr)
                    if gCrv_Split != Guid.Empty:
                        gCrvs_Split.append(gCrv_Split)
                    else:
                        bAddFail = True
                if not bAddFail and isinstance(rdObj_In, rd.CurveObject):
                    sc.doc.Objects.Delete(gCrv_In, quiet=False)
#            for t in ts_discontinuities:
#                pt = rgCrv0.PointAt(t)
#                if not pt: continue
#                gPt = sc.doc.Objects.AddPoint(pt)
#                if gPt != Guid.Empty:
#                    gPts.append(gPt)

    if bEcho:
        sContinuity = 'G2' if bG2_NotG1 else 'G1'
        if ts_discontinuities_All:
            s  = "{} {} discontinuities found in".format(
                    len(ts_discontinuities_All), sContinuity)
            if gCrvs_with_discont:
                s += " {} curves".format( len(gCrvs_with_discont))
            if i_Edges_with_discont_ct:
                if gCrvs_with_discont: s += ","
                s += ", {} edges(s)".format(len(i_Edges_with_discont_ct))
            s += "."
            print(s)
        else:
            print("No {} discontinuities found.".format(sContinuity))
    
    
    if gCrvs_Split:
        return gCrvs_Split
    else:
        return gCrvs_with_discont


def main():
    
    gCrvs0_Preselected = getPreselected(objectType=rd.ObjectType.Curve)
    
    objrefs = getInput()
    if objrefs is None: return

    bG2_NotG1 = Opts.values['bG2_NotG1']
    fAngleTol_Deg = Opts.values['fAngleTol_Deg']
    fRadiusTol = Opts.values['fRadiusTol']
    bAddPts = Opts.values['bAddPts']
    bSplitCurve = Opts.values['bSplitCurve']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    gCrvs_with_discont = processCurves(
        curvesAndEdges0=objrefs,
        bG2_NotG1=bG2_NotG1,
        fAngleTol_Deg=fAngleTol_Deg,
        fRadiusTol=fRadiusTol,
        bAddPts=bAddPts,
        bSplitCurve=bSplitCurve,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    
    if gCrvs0_Preselected:
        [sc.doc.Objects.Select(objectId=_) for _ in gCrvs0_Preselected]
    else:
        sc.doc.Objects.UnselectAll()
        [sc.doc.Objects.Select(objectId=_) for _ in gCrvs_with_discont]
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()