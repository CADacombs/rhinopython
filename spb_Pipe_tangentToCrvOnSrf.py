"""
This script will create a pipe with its center spine an offset from a curve on a surface.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums:
https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
240105: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import spb_OffsetNormal


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fRadius'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bUseFaceOfSelNakedEdge'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRebuildInCrv'; keys.append(key)
    values[key] = True
    names[key] = 'RebuildEdgeCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Try')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol'; keys.append(key)
    values[key] = max((1.0*sc.doc.ModelAbsoluteTolerance, 1e-6))
    names[key] = 'TargetTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bOutputCenterCrv'; keys.append(key)
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

        if key == 'fRadius':
            if abs(cls.riOpts[key].CurrentValue) < 1e-3:
                print("Radius input value is too small.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fTol':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Why is key, {}, here?  Value was not set or sticky-saved.".format(key))
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def _addCommonOptions(go):
    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    addOption('fRadius')
    addOption('bRebuildInCrv')
    addOption('fTol')
    addOption('bOutputCenterCrv')
    addOption('bEcho')
    addOption('bDebug')

    return idxs_Opt


def _getInput_Curve():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve on face")

    go.GeometryFilter = rd.ObjectType.Curve

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    #bPreselectedObjsChecked = False

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()
        addOption('bUseFaceOfSelNakedEdge')
        idxs_Opt.update(_addCommonOptions(go))

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref_CrvOnFace = go.Object(0)
            go.Dispose()

            return objref_CrvOnFace

        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_Face():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select base face")

    go.GeometryFilter = rd.ObjectType.Surface

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    go.AcceptNumber(True, acceptZero=True)
    
    idxs_Opt = {}

    bPreselectedObjsChecked = False

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()
        idxs_Opt.update(_addCommonOptions(go))

        res = go.Get()

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
            objref_Face = go.Object(0)

            go.Dispose()
    
            return objref_Face


        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_Click():
    """
    Click to toggle angle and/or direction with optional input.

    Returns:
        True: To recalculate and reloop
        False: To not recalculate and break out of loop with current output.
        None: To not recalculate and return without output.
    """

    go = ri.Custom.GetPoint()

    go.SetCommandPrompt("Left click to flip direction")

    go.SetCommandPromptDefault("Accept result")

    go.AcceptNumber(True, acceptZero=True)
    go.AcceptNothing(True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    key = 'FlipDir'; idxs_Opt[key] = go.AddOption(key)

    idxs_Opt.update(_addCommonOptions(go))

    res = go.Get()

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Nothing:
        go.Dispose()
        return False

    if res == ri.GetResult.Point:
        Opts.riOpts['fRadius'].CurrentValue = -Opts.riOpts['fRadius'].CurrentValue
        Opts.setValue('fRadius')
        go.Dispose()
        return True

    if res == ri.GetResult.Number:
        key = 'fRadius'
        Opts.riOpts[key].CurrentValue = go.Number()
        Opts.setValue(key)
        go.Dispose()
        return True

    # An option was selected.

    if go.OptionIndex() == idxs_Opt['FlipDir']:
        Opts.riOpts['fRadius'].CurrentValue = -Opts.riOpts['fRadius'].CurrentValue
        Opts.setValue('fRadius')
        go.Dispose()
        return True

    for key in idxs_Opt:
        if go.Option().Index == idxs_Opt[key]:
            Opts.setValue(key, go.Option().CurrentListOptionIndex)
            break

    go.Dispose()
    return True


def getInput():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick edge near its end for start of extension")

    go.GeometryFilter = rd.ObjectType.Curve


    def geomFilter_Curve(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if isinstance(geom, rg.BrepEdge):
            # DuplicateCurve gets the edge as a curve, which may be a subset of the EdgeCurve.
            rgC = geom.DuplicateCurve()
        elif isinstance(geom, rg.Curve):
            rgC = geom
        else:
            return False

        if rgC.IsPeriodic:
            print("Periodic curves are not supported.")
            return False

        return True


    go.SetCustomGeometryFilter(geomFilter_Curve)

    go.DisablePreSelect()

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bSameDegree')
        if not Opts.values['bSameDegree']:
            addOption('iDegree')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()


        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return objref

        if res == ri.GetResult.Number:
            key = 'iDegree'
            if go.Number() == Opts.riOpts[key].CurrentValue:
                continue
            Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 1:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                go.ClearCommandOptions()
                break


def _createCenterSpineWithCorrectedEndPts(ns_Loft, nc_TanEdge, bDebug=False):
    """
    This method corrects the end points of the center spine so that the ends
    of the pipe match the end points of the tangent curve.

    Returns: NurbsCurve
    """


    ns = ns_Loft

    # When the input curve is not used verbatim, the tangent curve spine is
    # the west isocurve of the loft.
    if nc_TanEdge is None:
        nc_TanEdge = ns.IsoCurve(1, ns.Domain(0).T0)

    # Center spine is the east isocurve of the loft.
    nc_CenterSpine = ns.IsoCurve(1, ns.Domain(0).T1)

    nc_CenterSpine_WIP = nc_CenterSpine.DuplicateCurve()

    #sc.doc.Objects.AddCurve(nc_CenterSpine_WIP)
    #sc.doc.Views.Redraw()
    #return

    # Extend as necessary, the start end of the center spine to be beyond
    # the closest point projection from the start of the tangent curve.

    while True:
        sc.escape_test()
        bSuccess, t_C_Start = nc_CenterSpine_WIP.ClosestPoint(nc_TanEdge.PointAtStart)
        if not bSuccess:
            raise Exception("ClosestPoint failed.")

        pt_C_Start = nc_CenterSpine_WIP.PointAt(t_C_Start)
        dist = pt_C_Start.DistanceTo(nc_CenterSpine_WIP.PointAtStart)
        if bDebug: sEval = "dist"; print("{}: {}".format(sEval, eval(sEval)))
        if dist > 1e-6:
            break

        nc_tmp = nc_CenterSpine_WIP.Extend(-0.001, 1.0)
        nc_CenterSpine_WIP.Dispose()
        nc_CenterSpine_WIP = nc_tmp


    # Do the same at the ends (T1) of the curves.

    while True:
        sc.escape_test()
        bSuccess, t_C_End = nc_CenterSpine_WIP.ClosestPoint(nc_TanEdge.PointAtEnd)
        if not bSuccess:
            raise Exception("ClosestPoint failed.")

        pt_C_End = nc_CenterSpine_WIP.PointAt(t_C_End)
        dist = pt_C_End.DistanceTo(nc_CenterSpine_WIP.PointAtEnd)
        if bDebug: sEval = "dist"; print("{}: {}".format(sEval, eval(sEval)))
        if dist > 1e-6:
            break

        nc_tmp = nc_CenterSpine_WIP.Extend(0, 1.001)
        nc_CenterSpine_WIP.Dispose()
        nc_CenterSpine_WIP = nc_tmp


    nc_tmp = nc_CenterSpine_WIP.Trim(t_C_Start, t_C_End)
    nc_CenterSpine_WIP.Dispose()

    return nc_tmp


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.crvs = []
        self.breps = []
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        for crv in self.crvs:
            bbox = crv.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)
        for brep in self.breps:
            bbox = brep.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PreDrawObjects(self, drawEventArgs):

        color = sc.doc.Layers.CurrentLayer.Color

        for crv in self.crvs:
            drawEventArgs.Display.DrawCurve(
                curve=crv,
                color=color,
                thickness=self.crv_thk)

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


def _createGeometryInteractively():
    """
    """

    objref_CrvToOffset = _getInput_Curve()
    if objref_CrvToOffset is None: return

    bUseFaceOfSelNakedEdge = Opts.values['bUseFaceOfSelNakedEdge']
    bDebug = Opts.values['bDebug']

    rgC_In, t_atPickedPt = objref_CrvToOffset.CurveParameter()
    if bDebug: sEval = "rgC_In"; print("{}: {}".format(sEval, eval(sEval)))
    if not rgC_In: return


    if isinstance(rgC_In, rg.BrepEdge) and bUseFaceOfSelNakedEdge:
        rgE_In = rgC_In
        if rgE_In.Valence == rg.EdgeAdjacency.Naked:
            idxF = objref_CrvToOffset.Edge().AdjacentFaces()[0]
            rgF_In = rgE_In.Brep.Faces[idxF]
        else:
            rgT_In = objref_CrvToOffset.Trim()
            if bDebug: sEval = "rgT_In.TrimIndex"; print("{}: {}".format(sEval, eval(sEval)))
            rgF_In = rgT_In.Face
    else:
        sc.doc.Objects.UnselectAll()
        objref_Face = _getInput_Face()
        if objref_Face is None: return

        rgF_In = objref_Face.Face()

    gBrep = objref_CrvToOffset.ObjectId
    sc.doc.Objects.UnselectAll()


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    fRadius = Opts.values['fRadius']
    bRebuildInCrv = Opts.values['bRebuildInCrv']
    fTol = Opts.values['fTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if isinstance(rgC_In, rg.PolyCurve):
        sEval = "rgC_In.RemoveNesting()"; print("{}: {}".format(sEval, eval(sEval)))
        #rgC_In.RemoveNesting()

    rgC_In_TrimmedToFace = spb_OffsetNormal.crvWithSpansCompletelyOnFace(
        rgC_In,
        rgF_In,
        t_atPickedPt,
        fTol=1.0*fTol,
        bDebug=bDebug)
    if rgC_In_TrimmedToFace is None: return

    if (
        bDebug and
        isinstance(rgC_In, rg.NurbsCurve) and
        isinstance(rgC_In_TrimmedToFace, rg.NurbsCurve)
    ):
        sEval = "rgC_In.EpsilonEquals(rgC_In_TrimmedToFace, 1e-6)"; print("{}: {}".format(sEval, eval(sEval)))

    if rgC_In_TrimmedToFace.IsClosed and fAngle_End_Deg:
        fAngle_End_Deg = Opts.values['fAngle_End_Deg'] = sc.sticky[Opts.stickyKeys['fAngle_End_Deg']] = None
        bVariableAngle = Opts.values['bVariableAngle'] = sc.sticky[Opts.stickyKeys['bVariableAngle']] = False


    sk_conduit = 'conduit({})'.format(__file__) # StickyKey
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit

    # For debugging.  Erase me.
    conduit = None
    conduit = DrawConduit()
    sc.sticky[sk_conduit] = conduit



    fSamplingDist = 100.0*sc.doc.ModelAbsoluteTolerance
    while True:
        ncs_toOffset, fDev_fromRebuilds = spb_OffsetNormal.prepareCrvToOffset(
            rgC_In_TrimmedToFace,
            bExplodePolyCrv=True,
            bRebuild=bRebuildInCrv,
            bSplitAtNonG2Knots=True,
            bMakeDeformable=False,
            fTol=0.1*fTol,
            bDebug=bDebug)

        if bDebug:
            sEval = "fDev_fromRebuilds"; print("{}: {}".format(sEval, eval(sEval)))


        def areCurvesDifferent(cA, ncsB):
            if len(ncsB) > 1:
                True

            ncA = cA.ToNurbsCurve()
            ncB = ncsB[0]

            epsEquals = ncA.EpsilonEquals(ncB, epsilon=1e-6)

            ncA.Dispose()

            return not epsEquals


        if areCurvesDifferent(rgC_In_TrimmedToFace, ncs_toOffset):
            print("Due to preprocessing, curve to offset is different than the input curve.")
        else:
            print("Curve on surface will be used verbatim as the curve to offset to make the center spine.")



        ncs_CenterSpines = []
        rgBs_Pipes = []
        ncs_Centers = []


        bFaceNormalIsReversedToSrf = rgF_In.OrientationIsReversed

        if bDebug:
            sEval = "fRadius>0.0"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "bFaceNormalIsReversedToSrf"; print("{}: {}".format(sEval, eval(sEval)))


        for nc_toOffset in ncs_toOffset:

            nc_toOffset.Domain = rg.Interval(0.0, 1.0)

            if bDebug:
                sEval = "nc_toOffset.Degree"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_toOffset.SpanCount"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_toOffset.Domain"; print("{}: {}".format(sEval, eval(sEval)))

            rc = spb_OffsetNormal.createOffsetCurve(
                rgCrv_In=nc_toOffset,
                rgSrf=rgF_In,
                bLoose=False,
                bAlignEndDirs=False,
                fDistance=fRadius,
                fTol=fTol-fDev_fromRebuilds,
                fSamplingDist=fSamplingDist,
                bDebug=bDebug)
            if rc is None:
                if bEcho: print("No solution found.")
                return

            nc_Offset, dev_Offset = rc

            if bDebug: sEval = "fTol-fDev_fromRebuilds-dev_Offset"; print("{}: {}".format(sEval, eval(sEval)))
            if fTol-fDev_fromRebuilds-dev_Offset < 0.0:
                raise ValueError("Curve(s) deviate more than the allotted tolerance.")

            #sc.doc.Objects.AddCurve(nc_Offset); sc.doc.Views.Redraw(); return

            if bDebug:
                sEval = "nc_Offset.Degree"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_Offset.SpanCount"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_Offset.Domain"; print("{}: {}".format(sEval, eval(sEval)))

            nc_Offset.Domain = rg.Interval(0.0, 1.0)

            ncs_CenterSpines.append(nc_Offset)

            rgBs_Loft = rg.Brep.CreateFromLoft(
                [nc_toOffset, nc_Offset],
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=rg.LoftType.Straight,
                closed=nc_Offset.IsClosed)
            if not rgBs_Loft:
                continue
            if len(rgBs_Loft) > 1:
                raise ValueError("{} breps in loft.  Should only be 1.".format(len(rgB_Loft)))

            rgB_Loft = rgBs_Loft[0]

            #sc.doc.Objects.AddBrep(rgB_Loft)
            #conduit.Enabled = False
            #sc.doc.Views.Redraw()
            #return

            def refreshBrepForNormalizedDomains(brep_In):
                ns = brep_In.Surfaces[0]
                ns.SetDomain(0, rg.Interval(0.0, 1.0))
                ns.SetDomain(1, rg.Interval(0.0, 1.0))
                brep_Out = ns.ToBrep()
                brep_In.Dispose()
                return brep_Out

            rgB_Loft = refreshBrepForNormalizedDomains(rgB_Loft)
            # 2 Curves3D and their BrepEdges will have domains of [-1, 0].  Don't know why yet.

            def reportBrepDomains(brep):
                ns = brep.Surfaces[0]
                ns.SetDomain(0, rg.Interval(0.0, 1.0))
                ns.SetDomain(1, rg.Interval(0.0, 1.0))
                brep_Out = ns.ToBrep()
                for ss in "brep.Curves2D", "brep.Curves3D", "brep.Trims", "brep.Edges":
                    print("{} domains:".format(ss))
                    cs = eval(ss)
                    for i in range(eval("{}.Count".format(ss))):
                        sEval = "cs[{}].Domain".format(i); print("{}: {}".format(i, eval(sEval)))

            #reportBrepDomains(rgB_Loft)

            #sc.doc.Objects.AddCurve(nc_toOffset)
            #sc.doc.Objects.AddCurve(nc_Offset)
            #sc.doc.Objects.AddBrep(rgB_Loft)
            #conduit.Enabled = False
            #1/0

            ns_Loft = rgB_Loft.Surfaces[0]
            if bDebug:
                sEval = "ns_Loft.Domain(0)"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "ns_Loft.Domain(1)"; print("{}: {}".format(sEval, eval(sEval)))

            nc_Center = _createCenterSpineWithCorrectedEndPts(
                ns_Loft=ns_Loft,
                nc_TanEdge=None if bRebuildInCrv else nc_toOffset,
                bDebug=bDebug)
            if not nc_Center: return
            ncs_Centers.append(nc_Center)

            rgB_Loft.Dispose()

            rgBs_Pipes_Res = rg.Brep.CreatePipe(
                rail=nc_Center,
                radius=abs(fRadius),
                localBlending=False,
                cap=rg.PipeCapMode.None,
                fitRail=False,
                absoluteTolerance=fTol-fDev_fromRebuilds-dev_Offset,
                angleToleranceRadians=sc.doc.ModelAngleToleranceRadians)
            if len(rgBs_Pipes_Res) == 0:
                raise Exception("CreatePipe returned None.")
            if len(rgBs_Pipes_Res) > 1:
                raise Exception("CreatePipe returned {} breps.".format(len(rgBs_Pipes_Res)))

            rgBs_Pipes.append(rgBs_Pipes_Res[0])

        conduit.curves = ncs_Centers
        conduit.breps = rgBs_Pipes

        if bDebug:
            conduit.crvs = ncs_CenterSpines

        conduit.Enabled = True

        sc.doc.Views.Redraw()

        if bEcho:
            sOut = []
            if len(ncs_CenterSpines) > 1: sOut.append("{} curves".format(len(ncs_CenterSpines)))
            if sOut:
                print("Calculated {}.".format(", ".join(sOut)))


        rc = _getInput_Click()

        conduit.Enabled = False

        if rc is None:
            for _ in ncs_CenterSpines: _.Dispose()
            for _ in rgBs_Pipes: _.Dispose()
            return

        if not rc:

            if bDebug:
                gs = [sc.doc.Objects.AddCurve(nc) for nc in ncs_toOffset]
                [sc.doc.Objects.AddCurve(nc) for nc in ncs_CenterSpines]
                [sc.doc.Objects.AddBrep(brep) for brep in rgBs_Pipes]
                [sc.doc.Objects.AddCurve(nc) for ncs in ncs_Centers for nc in ncs]
            else:
                #for _ in ncs_toOffset: _.Dispose()
                for _ in ncs_CenterSpines: _.Dispose()

            return rgBs_Pipes, ncs_Centers

        for _ in ncs_CenterSpines: _.Dispose()


        for _ in rgBs_Pipes: _.Dispose()


        fRadius = Opts.values['fRadius']
        bRebuildInCrv = Opts.values['bRebuildInCrv']
        fTol = Opts.values['fTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


def main():

    while True:
        rc = _createGeometryInteractively()
        if rc is None: return

        (
            rgBs_Pipes,
            ncs_CenterSpines,
            ) = rc

        bOutputCenterCrv = Opts.values['bOutputCenterCrv']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        for rgB in rgBs_Pipes:
            sc.doc.Objects.AddBrep(rgB)

        if bOutputCenterCrv:
            for nc in ncs_CenterSpines:
                sc.doc.Objects.AddCurve(nc)
        else:
            for nc in ncs_CenterSpines: nc.Dispose()

        sc.doc.Objects.UnselectAll()

        sc.doc.Views.Redraw()


if __name__ == '__main__': main()