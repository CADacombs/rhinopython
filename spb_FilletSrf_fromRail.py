"""
This script (will) creates a _FilletSrf-like surface from a BrepEdge.
This type of fillet is useful as an alternate to _Fin Direction=Tangent or
_ExtendSrf when it is desired to start a fillet/blend immediatly from an edge,
not requiring the face of the edge to be trimmed.

After the fillet is created, _Silhouette can be applied to it to for creating a profile
from which to extrude in the same direction as the silhouette view.

The Using_Sweep1 method uses the method as described at
https://discourse.mcneel.com/t/wish-filletsrftorail-with-constant-radius/158321/7
except the arcs may not be perpendicular to the edge. Instead, they are aligned to
the theoretical pipe's cross-sections.

TODO
Allow any curve on face to be selected.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums:
https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
231228-29: Created.
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
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fArcAngle'; keys.append(key)
    values[key] = 90.0
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bUseFaceOfSelNakedEdge'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bApproximate_NotSweep1'; keys.append(key)
    values[key] = False
    names[key] = 'Method'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Using_Sweep1', 'Approx_RC')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol'; keys.append(key)
    values[key] = max((1.0*sc.doc.ModelAbsoluteTolerance, 1e-6))
    names[key] = 'TargetTol'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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
    addOption('fArcAngle')
    addOption('bApproximate_NotSweep1')
    addOption('fTol')
    addOption('bEcho')
    addOption('bDebug')

    return idxs_Opt


def _getInput_Curve():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edge") #("Select curve on face")

    go.GeometryFilter = rd.ObjectType.EdgeFilter # rd.ObjectType.Curve

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
        #addOption('bUseFaceOfSelNakedEdge') # TODO: Disable this when non edges are allowed.
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

def _createArc_from_loft(ns_Loft, fRadius, fArcAngle, parameter, bFaceNormalIsReversedToSrf, bEdgeIsReversedToTrim, bDebug=False):

    ns = ns_Loft

    #edge = rgBrep_Loft.Edges[1]
    #crv = edge.DuplicateCurve()
    crv = ns.IsoCurve(1, ns.Domain(1).T1)
    crv.Domain = rg.Interval(0.0, 1.0)

    t = crv.Domain.T0 if parameter is None else parameter

    pt = crv.PointAt(t)
    vEdgeTan = crv.TangentAt(t)


    vNormal = ns.NormalAt(ns.Domain(0).T1, t)

    if fRadius > 0:
        vNormal = -vNormal
    if bFaceNormalIsReversedToSrf:
        vNormal = -vNormal
    if bEdgeIsReversedToTrim:
        vNormal = -vNormal

    vCross = rg.Vector3d.CrossProduct(vNormal, vEdgeTan)

    if fRadius > 0:
        vCross = -vCross
    if bFaceNormalIsReversedToSrf:
        vCross = -vCross
    if bEdgeIsReversedToTrim:
        vCross = -vCross

    plane = rg.Plane(
        origin=pt,
        xDirection=vCross,
        yDirection=vNormal)

    arc = rg.Arc(plane,
           center=pt,
           radius=abs(fRadius),
           angleRadians=Rhino.RhinoMath.ToRadians(fArcAngle))

    if bDebug:
        sEval = "t"; print("{}: {}".format(sEval, eval(sEval)))
        #sEval = "arc.IsValid"; print("{}: {}".format(sEval, eval(sEval)))
        #sEval = "arc.AngleDegrees"; print("{}: {}".format(sEval, eval(sEval)))
        #sEval = "arc.Center"; print("{}: {}".format(sEval, eval(sEval)))
        #sEval = "arc.StartPoint"; print("{}: {}".format(sEval, eval(sEval)))
        #sEval = "arc.EndPoint"; print("{}: {}".format(sEval, eval(sEval)))
        #sEval = "arc.Radius"; print("{}: {}".format(sEval, eval(sEval)))
        #sc.doc.Objects.AddCurve(crv)
        sc.doc.Objects.AddLine(rg.Line(start=pt, span=vEdgeTan))
        sc.doc.Objects.AddLine(rg.Line(start=pt, span=vNormal))
        sc.doc.Objects.AddLine(rg.Line(start=pt, span=vCross))
        #sc.doc.Objects.AddArc(arc)
        #sc.doc.Views.Redraw()

    return arc


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.crvs = []
        self.breps = []
        self.arcs = []
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        for crv in self.crvs:
            bbox = crv.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)
        for brep in self.breps:
            bbox = brep.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)
        for arc in self.arcs:
            bbox = arc.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PreDrawObjects(self, drawEventArgs):

        color = sc.doc.Layers.CurrentLayer.Color

        for arc in self.arcs:
            drawEventArgs.Display.DrawArc(
                arc=arc,
                color=color,
                thickness=self.crv_thk)


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


    rgEdge = objref_CrvToOffset.Edge()

    if rgEdge and bUseFaceOfSelNakedEdge and rgEdge.Valence == rg.EdgeAdjacency.Naked:
        idxF = objref_CrvToOffset.Edge().AdjacentFaces()[0]
        rgF_In = rgEdge.Brep.Faces[idxF]
    else:
        sc.doc.Objects.UnselectAll()

        objref_Face = _getInput_Face()
        if objref_Face is None: return

        sc.doc.Objects.UnselectAll()


        rgF_In = objref_Face.Face()


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    fRadius = Opts.values['fRadius']
    fArcAngle = Opts.values['fArcAngle']
    bApproximate_NotSweep1 = Opts.values['bApproximate_NotSweep1']
    fTol = Opts.values['fTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rgC_In, t_Crv0_Pick = objref_CrvToOffset.CurveParameter()

    # TODO: Change this when non-BrepEdges are allowed as input.
    edge_In = rgC_In

    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()


    rgC_In_TrimmedToFace = spb_OffsetNormal._crvWithSpansCompletelyOnFace(
        rgC_In,
        rgF_In,
        t_Crv0_Pick,
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
        sc.escape_test()

        ncs_toOffset, fDev_fromRebuilds = spb_OffsetNormal._prepareCrvToOffset(
            rgC_In_TrimmedToFace,
            bExplodePolyCrv=True,
            bRebuild=True,
            bSplitAtNonG2Knots=True,
            bMakeDeformable=False,
            fTol=0.5*fTol,
            bDebug=bDebug)

        if bDebug:
            sEval = "fDev_fromRebuilds"; print("{}: {}".format(sEval, eval(sEval)))


        ncs_Offset = []
        breps_Loft = []
        arcs_perOffset = []


        bFaceNormalIsReversedToSrf = rgF_In.OrientationIsReversed

        bEdgeIsReversedToTrim = rgF_In.Brep.Trims[edge_In.TrimIndices()[0]].IsReversed()

        if bDebug:
            sEval = "bFaceNormalIsReversedToSrf"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "bEdgeIsReversedToTrim"; print("{}: {}".format(sEval, eval(sEval)))


        for nc_toOffset in ncs_toOffset:

            nc_toOffset.Domain = rg.Interval(0.0, 1.0)

            if bDebug:
                sEval = "nc_toOffset.Degree"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_toOffset.SpanCount"; print("{}: {}".format(sEval, eval(sEval)))

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

            #sc.doc.Objects.AddCurve(nc_Offset); sc.doc.Views.Redraw(); return

            if bDebug:
                sEval = "nc_Offset.Degree"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_Offset.SpanCount"; print("{}: {}".format(sEval, eval(sEval)))

            nc_Offset.Domain = rg.Interval(0.0, 1.0)

            ncs_Offset.append(nc_Offset)

            rgB_Loft = rg.Brep.CreateFromLoft(
                [nc_toOffset, nc_Offset],
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=rg.LoftType.Straight,
                closed=nc_Offset.IsClosed)
            if not rgB_Loft:
                continue
            if len(rgB_Loft) > 1:
                raise ValueError("{} breps in loft.  Should only be 1.".format(len(rgB_Loft)))

            rgB_Loft = rgB_Loft[0]
            breps_Loft.append(rgB_Loft)

            rgB_Loft.Surfaces[0].SetDomain(0, rg.Interval(0.0, 1.0))
            rgB_Loft.Surfaces[0].SetDomain(1, rg.Interval(0.0, 1.0))

            #sc.doc.Objects.AddCurve(nc_toOffset)
            #sc.doc.Objects.AddCurve(nc_Offset)
            #sc.doc.Objects.AddBrep(rgB_Loft)
            #conduit.Enabled = False
            #1/0

            ns_Loft = rgB_Loft.Surfaces[0]
            if bDebug:
                sEval = "ns_Loft.Domain(0)"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "ns_Loft.Domain(1)"; print("{}: {}".format(sEval, eval(sEval)))

            if bDebug:
                sEval = "nc_toOffset.Domain"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_Offset.Domain"; print("{}: {}".format(sEval, eval(sEval)))


            if bApproximate_NotSweep1:
                ts_Grevs = nc_Offset.GrevilleParameters()
                if bDebug:
                    sEval = "len(ts_Grevs)"; print("{}: {}".format(sEval, eval(sEval)))
                arcs_perOffset.append([])
                for t in ts_Grevs:
                    arc = _createArc_from_loft(
                        ns_Loft,
                        fRadius,
                        fArcAngle,
                        parameter=t,
                        bFaceNormalIsReversedToSrf=bFaceNormalIsReversedToSrf,
                        bEdgeIsReversedToTrim=bEdgeIsReversedToTrim,
                        bDebug=bDebug)
                    arcs_perOffset[-1].append(arc)
            else:
                arc = _createArc_from_loft(
                    ns_Loft,
                    fRadius,
                    fArcAngle,
                    parameter=None,
                    bFaceNormalIsReversedToSrf=bFaceNormalIsReversedToSrf,
                    bEdgeIsReversedToTrim=bEdgeIsReversedToTrim,
                    bDebug=bDebug)
                arcs_perOffset.append([arc])


        conduit.arcs = [arc for arcs in arcs_perOffset for arc in arcs]

        if bDebug:
            conduit.breps = breps_Loft
            conduit.crvs = ncs_Offset

        conduit.Enabled = True

        sc.doc.Views.Redraw()

        if bEcho:
            sOut = []
            if len(ncs_Offset) > 1: sOut.append("{} curves".format(len(ncs_Offset)))
            if sOut:
                print("Calculated {}.".format(", ".join(sOut)))


        rc = _getInput_Click()

        conduit.Enabled = False

        if rc is None:
            for _ in ncs_Offset: _.Dispose()
            for _ in breps_Loft: _.Dispose()
            return

        for _ in ncs_Offset: _.Dispose()


        if not rc:

            if bDebug:
                [sc.doc.Objects.AddCurve(nc) for nc in ncs_toOffset]
                [sc.doc.Objects.AddBrep(brep) for brep in breps_Loft]
                [sc.doc.Objects.AddArc(arc) for arcs in arcs_perOffset for arc in arcs]

            for _ in ncs_Offset: _.Dispose()

            for _ in ncs_toOffset: _.Dispose()

            return (
                arcs_perOffset,
                breps_Loft,
                )

        for _ in breps_Loft: _.Dispose()


        fRadius = Opts.values['fRadius']
        fArcAngle = Opts.values['fArcAngle']
        bApproximate_NotSweep1 = Opts.values['bApproximate_NotSweep1']
        fTol = Opts.values['fTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


def _addFilletApproximation(arcs, brep):
    acs = [rg.ArcCurve(arc) for arc in arcs]

    ns = brep.Surfaces[0]
    nc = ns.IsoCurve(1, ns.Domain(1).T0)

    #s1 = rg.SweepOneRail()

    #sProps = (
    #    "AngleToleranceRadians",
    #    "ClosedSweep",
    #    "GlobalShapeBlending",
    #    "IsFreeform",
    #    "IsRoadlike",
    #    "IsRoadlikeFront",
    #    "IsRoadlikeTop",
    #    "IsRoadlineRight",
    #    "MiterType",
    #    "SweepTolerance",
    #    )

    #for sProp in sProps:
    #    sEval = "s1.{}".format(sProp); print("{}: {}".format(sEval, eval(sEval)))

    #s1.ClosedSweep = False
    #s1.MiterType = 0

    #for sProp in sProps:
    #    sEval = "s1.{}".format(sProp); print("{}: {}".format(sEval, eval(sEval)))

    #sweeps = s1.PerformSweep(
    #    rail=nc,
    #    crossSections=acs)

    #if not sweeps: return
    #for sweep in sweeps:
    #    sc.doc.Objects.AddBrep(sweep)


    sweeps = rg.Brep.CreateFromSweep(
        rail=nc,
        shapes=acs,
        startPoint=rg.Point3d.Unset,
        endPoint=rg.Point3d.Unset,
        frameType=rg.SweepFrame.Freeform,
        roadlikeNormal=rg.Vector3d.Unset,
        closed=nc.IsClosed,
        blendType=rg.SweepBlend.Local,
        miterType=rg.SweepMiter.None,
        tolerance=0.5*sc.doc.ModelAbsoluteTolerance,
        rebuildType=rg.SweepRebuild.None,
        rebuildPointCount=0,
        refitTolerance=0.0)

    #sweeps = rg.Brep.CreateFromSweep(
    #    rail=nc,
    #    shapes=acs,
    #    closed=nc.IsClosed,
    #    tolerance=sc.doc.ModelAbsoluteTolerance)

    if not sweeps: return
    for sweep in sweeps:
        sc.doc.Objects.AddBrep(sweep)
    sc.doc.Views.Redraw()


def _addSweep1_with_user_interaction(arc, loft, bDebug):
    gArc = sc.doc.Objects.AddArc(arc)

    gLoft = sc.doc.Objects.AddBrep(loft)

    sc.doc.Views.Redraw()

    #rdLoft = sc.doc.Objects.FindId(gLoft)

    ## Presuming index of edge remains the same since loft is always defined the same way.
    #compIdx = rg.ComponentIndex(
    #        rg.ComponentIndexType.BrepEdge,
    #        index=1)
    ## TODO: Set syncHighlight to bDebug
    #rdLoft.SelectSubObject(
    #    componentIndex=compIdx,
    #    select=True,
    #    syncHighlight=True,
    #    persistentSelect=True)

    # TODO: Set echo to bDebug
    Rhino.RhinoApp.RunScript("_Echo", echo=True)
    script  = "_-Sweep1".format(gArc)
    script += " _Pause".format(gArc)
    script += " _SelId {} _Enter".format(gArc)
    # _ShapeBlending is ignored when there is only 1 shape.
    script += " _Style=AlignWithSurface _Simplify=None _Enter"
    Rhino.RhinoApp.RunScript(script, echo=True)

    if not bDebug:
        sc.doc.Objects.Delete(objectId=gArc, quiet=False)
        sc.doc.Objects.Delete(objectId=gLoft, quiet=False)

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def main():

    while True:
        sc.escape_test() # TODO: Erase me.

        rc = _createGeometryInteractively()
        if rc is None: return

        bApproximate_NotSweep1 = Opts.values['bApproximate_NotSweep1']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        arcs_perOffset, rgBs_Lofts = rc

        if bApproximate_NotSweep1:
            for arcs_thisOffset, loft in zip(arcs_perOffset, rgBs_Lofts):
                _addFilletApproximation(arcs_thisOffset, loft)
            sc.doc.Objects.UnselectAll()
        else:
            if bDebug:
                sc.doc.Objects.UnselectAll()
                sc.doc.Views.Redraw()
                continue

            for arcs_thisOffset, loft in zip(arcs_perOffset, rgBs_Lofts):
                _addSweep1_with_user_interaction(arcs_thisOffset[0], loft, bDebug)

        sc.doc.Views.Redraw()


if __name__ == '__main__': main()