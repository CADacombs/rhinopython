"""
This script will split edges of Breps at G1- or G2-discontinuous knots.
It can be used after a _MergeAllEdges so that G1 maximum continuities,
have vertices.

For G2 disontinuities, a choice of two routines are available.

1. Mine
   It allows input of a tolerance for absolute radius difference.
   It also uses angle tolerance for the difference in curvature vector
   directions.  ( _GCon also uses the same tolerance for both.)
2. RhinoCommmon's Curve.GetNextDiscontinuity
    Its default value to find G2 discontinuities where at least one
    of the following is true:
        1. abs(k0-k1) / max(k0,k1) > 0.05,
        where k0 and k1 are magnitudes of the curvature vectors
        before and after the knot.
        2. The difference in curvature vectors is > 2.0 degrees.
    Using a different value has a limited effect.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
211119: Created.
221027: Removed requiring trim selection during mated edge selection.
221223-24: Added options to control G2 discontinuity search.

TODO: Add absolute curvature difference tolerance at least for RC method?
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Math


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAllowBrepSelection'; keys.append(key)
    values[key] = False
    names[key] = 'SelectionMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'EdgesOnly', 'EdgesandBreps')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bIncludeMated'; keys.append(key)
    values[key] = True
    names[key] = 'IncludeMated'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fAngleTol_Deg'; keys.append(key)
    values[key] = sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bG2_NotG1'; keys.append(key)
    values[key] = True
    names[key] = 'Discont'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'G1', 'G2')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseSPB_NotRC'; keys.append(key)
    values[key] = True
    names[key] = 'G2Method'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'RC', 'SPB')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fRadiusTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fRCCrvtrTol'; keys.append(key)
    values[key] = Rhino.RhinoMath.SqrtEpsilon
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fAngleTol_Deg':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-3:
                cls.riOpts[key].CurrentValue = 1e-3
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fRadiusTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fRCCrvtrTol':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edges to split")

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        if Opts.values['bAllowBrepSelection']:
            go.SetCommandPrompt("Select breps and/or edges")
            go.SetCommandPromptDefault("All normal breps when none are selected")
            go.GeometryFilter = rd.ObjectType.Brep | rd.ObjectType.Curve
        else:
            go.SetCommandPrompt("Select edges")
            go.GeometryFilter = rd.ObjectType.Curve

        go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.EdgeCurve |
            ri.Custom.GeometryAttributeFilter.BoundaryEdge)

        if Opts.values['bIncludeMated']:
            go.GeometryAttributeFilter |= (
                ri.Custom.GeometryAttributeFilter.MatedEdge)


        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bAllowBrepSelection')
        addOption('bIncludeMated')
        addOption('fAngleTol_Deg')
        addOption('bG2_NotG1')
        if Opts.values['bG2_NotG1']:
            addOption('bUseSPB_NotRC')
            if Opts.values['bUseSPB_NotRC']:
                addOption('fRadiusTol')
            else:
                addOption('fRCCrvtrTol')
        addOption('bEcho')
        addOption('bDebug')
        if Opts.values['bDebug']:
            addOption('bAddRefs')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

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
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fAngleTol_Deg'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def groupObjrefsPerBrep_OLD(objrefs):
    """ Returns nested list per Brep. """
    gBs = []
    objrefs_PerB = []

    for objref in objrefs:
        if objref.ObjectId in gBs:
            objrefs_PerB[gBs.index(objref.ObjectId)].append(objref)
        else:
            gBs.append(objref.ObjectId)
            objrefs_PerB.append([objref])

    return objrefs_PerB


def createSortedBrepsAndEdges_FromObjRefs(objrefs, bIncludeMated):
    """
    Parameters:
        ObjRefs of Breps or Edges
    Returns on success:
        list(GUIDs of breps) and list(list(Edge indices) per brep)
    """

    rdBs = []
    gBs = []
    idxEs = []

   # Organize input into synchronized lists of brep ids and edge indices.
    for objref in objrefs:
        rdB = objref.Object()
        gB = objref.ObjectId
        rgE = objref.Edge()
        idxE = rgE.EdgeIndex if rgE else None
        rgB = rdB.BrepGeometry


        if gB not in gBs:
            # Brep not in list, so add it as well as objref's edge.
            rdBs.append(rdB)
            gBs.append(gB)
            idxEs.append([])
            iB = len(gBs) - 1
        else:
            # Brep already in list.
            iB = gBs.index(gB)


        # 4 states: idxE is None or int, bIncludeMated is True or False.

        if (idxE is None) and bIncludeMated:
            idxEs[iB] = range(rgB.Edges.Count)
        elif (idxE is None) and not bIncludeMated:
            for idxE in xrange(rgB.Edges.Count):
                rgE = rgB.Edges[idxE]
                if rgE.Valence == rg.EdgeAdjacency.Interior:
                    continue
                idxEs[iB].append(idxE)
            if len(idxEs[iB]) == 0:
                print("This brep has no valid edges.  Add code to handle this?")
        elif (idxE is not None) and bIncludeMated:
            if idxE not in idxEs[iB]:
                idxEs[iB].append(idxE)
        elif (
            (idxE is not None) and
            (rgB.Edges[idxE].Valence == rg.EdgeAdjacency.Naked)
            ):
            if idxE not in idxEs[iB]:
                idxEs[iB].append(idxE)

    return rdBs, idxEs


def _continuityVectorsAtCurveParameter(nc, t, side=rg.CurveEvaluationSide.Default):
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


def getGeometricDiscontinuities_MyTake(rgCrv_In, bG2_NotG1, fAngleTol_Deg, fRadiusTol, bDebug=False):
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
                ) = _continuityVectorsAtCurveParameter(
                    nc,
                    nc.Knots[nc.Knots.Count - 1],
                    rg.CurveEvaluationSide.Below)
        else:
            (
                vG0_Below,
                vG1_Below,
                vG2_Below,
                vG3_Below,
                ) = _continuityVectorsAtCurveParameter(
                    nc,
                    nc.Knots[iK],
                    rg.CurveEvaluationSide.Below)


        (
            vG0_Above,
            vG1_Above,
            vG2_Above,
            vG3_Above,
            ) = _continuityVectorsAtCurveParameter(
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
                            print("{:.{}f} radius difference.".format((
                                delta_radius, sc.doc.ModelDistanceDisplayPrecision)))
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


def getGeometricDiscontinuities_UsingRhinoCommonMethod(rgCrv, bG2_NotG1=True, fAngleTol_Deg=None, fRCCrvtrTol=Rhino.RhinoMath.SqrtEpsilon):
    """
    Returns: list(float(parameters of discontinuities))
    """

    if fAngleTol_Deg is None: fAngleTol_Deg = sc.doc.ModelAngleToleranceDegrees

    cosAngleTolerance = Math.Cos(Rhino.RhinoMath.ToRadians(fAngleTol_Deg))

    curvatureTolerance = fRCCrvtrTol

    t0 = rgCrv.Domain.Min
    t1 = rgCrv.Domain.Max

    if rgCrv.IsClosed:
        if bG2_NotG1:
            continuityType = rg.Continuity.G2_locus_continuous
        else:
            continuityType = rg.Continuity.G1_locus_continuous
    else:
        if bG2_NotG1:
            continuityType = rg.Continuity.G2_continuous
        else:
            continuityType = rg.Continuity.G1_continuous

    ts = []

    get_next = True

    while get_next:
        sc.escape_test()

        get_next, t = rg.Curve.GetNextDiscontinuity(
            rgCrv,
            continuityType,
            t0,
            t1,
            cosAngleTolerance,
            curvatureTolerance)

        if get_next:
            ts.append(t)
            t0 = t # Advance to the next parameter.

    return ts


def processObjRefs(objrefs, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bG2_NotG1 = getOpt('bG2_NotG1')
    bUseSPB_NotRC = getOpt('bUseSPB_NotRC')
    fRadiusTol = getOpt('fRadiusTol')
    fRCCrvtrTol = getOpt('fRCCrvtrTol')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    bAddRefs = getOpt('bAddRefs')



    objrefs_per_B = groupObjrefsPerBrep_OLD(objrefs)

    rc = createSortedBrepsAndEdges_FromObjRefs(objrefs, Opts.values['bIncludeMated'])
    if rc is None: return

    rdBs, idxEs_PerB = rc

    iCt_EdgeIncr = 0
    iCt_Breps_Mod = 0

    #for objrefs_same_B in objrefs_per_B:

    #    edges_In_SameB = [o.Edge() for o in objrefs_same_B]
    #    idxs_Es = [e.EdgeIndex for e in edges_In_SameB]

    #    # Sort in reverse EdgeIndex order.
    #    edges_In_SameB = [e for i,e in sorted(zip(idxs_Es, edges_In_SameB), reverse=True)]

    #    rgB_Start = edges_In_SameB[0].Brep

    #    rgB_WIP = rgB_Start.DuplicateBrep()

    #    edge_count_Start = rgB_WIP.Edges.Count

    #    for edge in edges_In_SameB:

    #        crv_Start = edge.DuplicateCurve()

    #        if bG2_NotG1 and bUseSPB_NotRC:
    #            ts = getGeometricDiscontinuities_MyTake(
    #                crv_Start, bG2_NotG1, fAngleTol_Deg, fRadiusTol, bDebug)
    #        else:
    #            ts = getGeometricDiscontinuities_UsingRhinoCommonMethod(
    #                crv_Start, bG2_NotG1, fAngleTol_Deg, fRCCrvtrTol)

    #        rgB_WIP.Edges.SplitEdgeAtParameters(
    #            edgeIndex=edge.EdgeIndex,
    #            edgeParameters=ts)

    #    rgB_WIP.Compact()

    #    if rgB_WIP.Edges.Count != edge_count_Start:
    #        iCt_EdgeIncr += rgB_WIP.Edges.Count - edge_count_Start

    #        bReplaced = sc.doc.Objects.Replace(
    #            objectId=objrefs_same_B[0].ObjectId,
    #            brep=rgB_WIP)

    #        iCt_Breps_Mod += 1



    for rdB, idxEs in zip(rdBs, idxEs_PerB):

        rgB_WIP = rdB.BrepGeometry.DuplicateBrep()

        edge_count_Start = rgB_WIP.Edges.Count

        for idxE in idxEs:

            edge = rgB_WIP.Edges[idxE]

            crv_Start = edge.DuplicateCurve()

            if bG2_NotG1 and bUseSPB_NotRC:
                ts = getGeometricDiscontinuities_MyTake(
                    crv_Start, bG2_NotG1, fAngleTol_Deg, fRadiusTol, bDebug)
            else:
                ts = getGeometricDiscontinuities_UsingRhinoCommonMethod(
                    crv_Start, bG2_NotG1, fAngleTol_Deg, fRCCrvtrTol)

            rgB_WIP.Edges.SplitEdgeAtParameters(
                edgeIndex=edge.EdgeIndex,
                edgeParameters=ts)

        rgB_WIP.Compact()

        if rgB_WIP.Edges.Count != edge_count_Start:
            iCt_EdgeIncr += rgB_WIP.Edges.Count - edge_count_Start

            bReplaced = sc.doc.Objects.Replace(
                objectId=rdB.Id,
                brep=rgB_WIP)

            iCt_Breps_Mod += 1

    if iCt_EdgeIncr == 0:
        print("No edges were split.")
    else:
        print("Edge count of {} breps increased by {}.".format(
            iCt_Breps_Mod, iCt_EdgeIncr))


def main():

    objrefs = getInput()
    if objrefs is None: return

    fAngleTol_Deg = Opts.values['fAngleTol_Deg']
    bG2_NotG1 = Opts.values['bG2_NotG1']
    bUseSPB_NotRC = Opts.values['bUseSPB_NotRC']
    fRadiusTol = Opts.values['fRadiusTol']
    fRCCrvtrTol = Opts.values['fRCCrvtrTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']
    bAddRefs = Opts.values['bAddRefs']


    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBrep_Ret = processObjRefs(
        objrefs,
        fAngleTol_Deg=fAngleTol_Deg,
        bG2_NotG1=bG2_NotG1,
        bUseSPB_NotRC=bUseSPB_NotRC,
        fRadiusTol=fRadiusTol,
        fRCCrvtrTol=fRCCrvtrTol,
        bEcho=bEcho,
        bDebug=bDebug,
        bAddRefs=bAddRefs,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()