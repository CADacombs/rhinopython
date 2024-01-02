"""
This script creates a _FilletSrf-like surface from a BrepEdge.  It is up to the user to verify
accuracy using _EdgeContinuity, etc.

The option, Method=Using_Sweep1, uses the method described at
https://discourse.mcneel.com/t/wish-filletsrftorail-with-constant-radius/158321/7
except
    1. The sweep shapes can be multiple arcs.
    2. The arcs are aligned to the theoretical pipe's cross-sections, not necessarily
perpendicular to the input edge.

Neither Brep.CreateFromSweep nor Geometry.SweepOneRail can create a simple sweep:
https://discourse.mcneel.com/t/simple-sweep/159677/2
If this changes, the process can avoid user interaction during the actual creation of the sweep.


Send any questions, comments, or script development service needs to @spb on the McNeel Forums:
https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
231228-240102: Created.

TODO
- Allow any curve on face to be selected.

Orient arc Method:
    1. Find closest point on center spine of pipe curve.
    2. Center arc on point from #1 on plane perpendicular to the same curve.
    3. Translate so the end of the arc is coincident with Greville point on (rebuilt) input edge curve.
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

    key = 'bRebuildInCrv'; keys.append(key)
    values[key] = True
    names[key] = 'RebuildEdgeCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Try')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bApproximate_NotSweep1'; keys.append(key)
    values[key] = False
    names[key] = 'Method'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Using_Sweep1', 'Approx_RC')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOnlyOneArc_NotAll'; keys.append(key)
    values[key] = True
    names[key] = 'ArcsForSweep'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'All', 'One')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseSweep1Dialog'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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
    addOption('bRebuildInCrv')
    addOption('bApproximate_NotSweep1')
    if not Opts.values['bApproximate_NotSweep1']:
        addOption('bOnlyOneArc_NotAll')
        addOption('bUseSweep1Dialog')
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
        #addOption('bUseFaceOfSelNakedEdge') # TODO: Disable this when non-edges are allowed.
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


def _createArcs(ns_Loft, nc_TanEdge, fRadius, fArcAngle, bFaceNormalIsReversedToSrf, bEdgeIsReversedToTrim, bDebug=False):
    """
    Returns:
        list(rg.Arc)

    Arcs are aligned to center curve of theoretical pipe but are set at each
    tangent edge's Greville point.
    """


    ns = ns_Loft

    if nc_TanEdge is None:
        nc_TanEdge = ns.IsoCurve(1, ns.Domain(0).T0)

    nc_E = ns.IsoCurve(1, ns.Domain(0).T1)

    nc_E_Ext = nc_E.Extend(-0.1, 1.1)

    #if bDebug:
    #    sc.doc.Objects.AddSurface(ns)
    #    sc.doc.Objects.AddCurve(nc_W)
    #    sc.doc.Objects.AddCurve(nc_E)
    #    sc.doc.Objects.AddCurve(nc_E_Ext)

    ts_Grevs_W = nc_TanEdge.GrevilleParameters()
    pts_Grevs_W = nc_TanEdge.GrevillePoints(all=False)

    ts_Closest_on_E = []
    pts_Closest_on_N = []
    arcs = []

    #if bDebug:
    #    attr = rd.ObjectAttributes()
    #    attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    #    attr.ColorSource = rd.ObjectColorSource.ColorFromObject

    for i, ptS in enumerate(pts_Grevs_W):
        bSuccess, t_E = nc_E_Ext.ClosestPoint(ptS)
        if not bSuccess:
            raise Exception("ClosestPoint failed.")

        ts_Closest_on_E.append(t_E)
        pt_E = nc_E_Ext.PointAt(t_E)
        pts_Closest_on_N.append(pt_E)

        #if bDebug and i > 0:
        #    sc.doc.Objects.AddPoint(pt_E)


        vEdgeTan_E = nc_E.TangentAt(t_E)

        iCt_Tan_Reversed = 0
        if fRadius > 0:
            vEdgeTan_E = -vEdgeTan_E
            iCt_Tan_Reversed += 1
        if bFaceNormalIsReversedToSrf:
            vEdgeTan_E = -vEdgeTan_E
            iCt_Tan_Reversed += 1
        if bEdgeIsReversedToTrim:
            vEdgeTan_E = -vEdgeTan_E
            iCt_Tan_Reversed += 1

        if bDebug and i > 0: sEval = "iCt_Tan_Reversed"; print("{}: {}".format(sEval, eval(sEval)))

        #attr.ObjectColor = attr.ObjectColor.Lime
        #sc.doc.Objects.AddLine(rg.Line(start=pt_E, span=vEdgeTan_E), attr)


        vNormal_E = ns.NormalAt(ns.Domain(0).T1, t_E)

        iCt_Normal_Reversed = 0
        if fRadius > 0:
            vNormal_E = -vNormal_E
            iCt_Normal_Reversed += 1
        if bFaceNormalIsReversedToSrf:
            vNormal_E = -vNormal_E
            iCt_Normal_Reversed += 1
        if bEdgeIsReversedToTrim:
            vNormal_E = -vNormal_E
            iCt_Normal_Reversed += 1

        if bDebug and i > 0: sEval = "iCt_Normal_Reversed"; print("{}: {}".format(sEval, eval(sEval)))

        #attr.ObjectColor = attr.ObjectColor.Red
        #sc.doc.Objects.AddLine(rg.Line(start=pt_E, span=vNormal_E), attributes=attr)

        vCross = rg.Vector3d.CrossProduct(vNormal_E, vEdgeTan_E)

        #attr.ObjectColor = attr.ObjectColor.Blue
        #sc.doc.Objects.AddLine(rg.Line(start=pt_E, span=vCross), attr)

        #return

        plane = rg.Plane(
            origin=pt_E,
            xDirection=vCross,
            yDirection=vNormal_E)

        arc = rg.Arc(plane,
               center=pt_E,
               radius=abs(fRadius),
           angleRadians=Rhino.RhinoMath.ToRadians(fArcAngle))

        # Translate arc to the west curve for any distance error between the curves.
        xform_Fix = rg.Transform.Translation(ptS - arc.StartPoint)
        arc.Transform(xform_Fix)

        arcs.append(arc)

        #if bDebug:
            #sEval = "t_E"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "arc.IsValid"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "arc.AngleDegrees"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "arc.Center"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "arc.StartPoint"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "arc.EndPoint"; print("{}: {}".format(sEval, eval(sEval)))
            #sEval = "arc.Radius"; print("{}: {}".format(sEval, eval(sEval)))
            #attr.ObjectColor = attr.ObjectColor.Red
            #sc.doc.Objects.AddLine(rg.Line(start=pt_E, span=vNormal_E), attributes=attr)
            #attr.ObjectColor = attr.ObjectColor.Lime
            #sc.doc.Objects.AddLine(rg.Line(start=pt_E, span=vEdgeTan_E), attr)
            #attr.ObjectColor = attr.ObjectColor.Blue
            #sc.doc.Objects.AddLine(rg.Line(start=pt_E, span=vCross), attr)
            #sc.doc.Objects.AddArc(arc)
            #sc.doc.Views.Redraw()

    return arcs


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

    rgE_In = objref_CrvToOffset.Edge()
    if not rgE_In: return
    sEval = "rgE_In.EdgeIndex"; print("{}: {}".format(sEval, eval(sEval)))

    rgT_In = objref_CrvToOffset.Trim()
    sEval = "rgT_In.TrimIndex"; print("{}: {}".format(sEval, eval(sEval)))
    rgF_In = rgT_In.Face

    # TODO: Enable and modify this block when non-edges are allowed as input.
    #if rgE_In and bUseFaceOfSelNakedEdge and rgE_In.Valence == rg.EdgeAdjacency.Naked:
    #    idxF = objref_CrvToOffset.Edge().AdjacentFaces()[0]
    #    rgF_In = rgE_In.Brep.Faces[idxF]
    #else:
    #    sc.doc.Objects.UnselectAll()
    #    objref_Face = _getInput_Face()
    #    if objref_Face is None: return

    #    rgF_In = objref_Face.Face()

    gBrep = objref_CrvToOffset.ObjectId
    sc.doc.Objects.UnselectAll()


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    fRadius = Opts.values['fRadius']
    fArcAngle = Opts.values['fArcAngle']
    bRebuildInCrv = Opts.values['bRebuildInCrv']
    bApproximate_NotSweep1 = Opts.values['bApproximate_NotSweep1']
    bUseSweep1Dialog = Opts.values['bUseSweep1Dialog']
    fTol = Opts.values['fTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rgC_In, t_Crv0_Pick = objref_CrvToOffset.CurveParameter()

    # TODO: Change this when non-BrepEdges are allowed as input.
    edge_In = rgC_In

    if isinstance(rgC_In, rg.PolyCurve):
        rgC_In.RemoveNesting()

    rgC_In_TrimmedToFace = rgC_In
    # TODO: Re-enable the following when non-edges are allowed as input.  This line above should then be deleted.

    #rgC_In_TrimmedToFace = spb_OffsetNormal.crvWithSpansCompletelyOnFace(
    #    rgC_In,
    #    rgF_In,
    #    t_Crv0_Pick,
    #    fTol=1.0*fTol,
    #    bDebug=bDebug)
    #if rgC_In_TrimmedToFace is None: return

    #if (
    #    bDebug and
    #    isinstance(rgC_In, rg.NurbsCurve) and
    #    isinstance(rgC_In_TrimmedToFace, rg.NurbsCurve)
    #):
    #    sEval = "rgC_In.EpsilonEquals(rgC_In_TrimmedToFace, 1e-6)"; print("{}: {}".format(sEval, eval(sEval)))

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
            fTol=0.5*fTol,
            bDebug=bDebug)

        if bDebug:
            sEval = "fDev_fromRebuilds"; print("{}: {}".format(sEval, eval(sEval)))


        ncs_Offset = []
        breps_Loft = []
        arcs_perOffset = []


        bFaceNormalIsReversedToSrf = rgF_In.OrientationIsReversed

        bEdgeIsReversedToTrim = rgT_In.IsReversed()

        if bDebug:
            sEval = "fRadius>0.0"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "bFaceNormalIsReversedToSrf"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "bEdgeIsReversedToTrim"; print("{}: {}".format(sEval, eval(sEval)))


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

            #sc.doc.Objects.AddCurve(nc_Offset); sc.doc.Views.Redraw(); return

            if bDebug:
                sEval = "nc_Offset.Degree"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_Offset.SpanCount"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "nc_Offset.Domain"; print("{}: {}".format(sEval, eval(sEval)))

            nc_Offset.Domain = rg.Interval(0.0, 1.0)

            ncs_Offset.append(nc_Offset)

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

            breps_Loft.append(rgB_Loft)

            #sc.doc.Objects.AddCurve(nc_toOffset)
            #sc.doc.Objects.AddCurve(nc_Offset)
            #sc.doc.Objects.AddBrep(rgB_Loft)
            #conduit.Enabled = False
            #1/0

            ns_Loft = rgB_Loft.Surfaces[0]
            if bDebug:
                sEval = "ns_Loft.Domain(0)"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "ns_Loft.Domain(1)"; print("{}: {}".format(sEval, eval(sEval)))

            arcs = _createArcs(
                ns_Loft=ns_Loft,
                nc_TanEdge=None if bRebuildInCrv else nc_toOffset,
                fRadius=fRadius,
                fArcAngle=fArcAngle,
                bFaceNormalIsReversedToSrf=bFaceNormalIsReversedToSrf,
                bEdgeIsReversedToTrim=bEdgeIsReversedToTrim,
                bDebug=bDebug)
            if not arcs: return
            arcs_perOffset.append(arcs)


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

        if not rc:

            if bDebug:
                gs = [sc.doc.Objects.AddCurve(nc) for nc in ncs_toOffset]
                [sc.doc.Objects.AddCurve(nc) for nc in ncs_Offset]
                [sc.doc.Objects.AddBrep(brep) for brep in breps_Loft]
                [sc.doc.Objects.AddArc(arc) for arcs in arcs_perOffset for arc in arcs]
            else:
                #for _ in ncs_toOffset: _.Dispose()
                for _ in ncs_Offset: _.Dispose()

            return (
                arcs_perOffset,
                breps_Loft,
                ncs_toOffset,
                gBrep,
                )

        for _ in ncs_Offset: _.Dispose()


        for _ in breps_Loft: _.Dispose()


        fRadius = Opts.values['fRadius']
        fArcAngle = Opts.values['fArcAngle']
        bRebuildInCrv = Opts.values['bRebuildInCrv']
        bApproximate_NotSweep1 = Opts.values['bApproximate_NotSweep1']
        bUseSweep1Dialog = Opts.values['bUseSweep1Dialog']
        fTol = Opts.values['fTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


def _addFilletApproximation_RC(arcs, nc_toOffset, brep):
    acs = [rg.ArcCurve(arc) for arc in arcs]

    ns = brep.Surfaces[0]
    nc_forRail = nc_toOffset if nc_toOffset else ns.IsoCurve(1, ns.Domain(0).T0)

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
    #    rail=nc_W,
    #    crossSections=acs)

    #if not sweeps: return
    #for sweep in sweeps:
    #    sc.doc.Objects.AddBrep(sweep)


    sweeps = rg.Brep.CreateFromSweep(
        rail=nc_forRail,
        shapes=acs,
        startPoint=rg.Point3d.Unset,
        endPoint=rg.Point3d.Unset,
        frameType=rg.SweepFrame.Freeform,
        roadlikeNormal=rg.Vector3d.Unset,
        closed=nc_forRail.IsClosed,
        blendType=rg.SweepBlend.Local,
        miterType=rg.SweepMiter.None,
        tolerance=0.5*sc.doc.ModelAbsoluteTolerance,
        rebuildType=rg.SweepRebuild.None,
        rebuildPointCount=0,
        refitTolerance=0.0)

    #sweeps = rg.Brep.CreateFromSweep(
    #    rail=nc_W,
    #    shapes=acs,
    #    closed=nc_W.IsClosed,
    #    tolerance=sc.doc.ModelAbsoluteTolerance)

    if not sweeps: return
    for sweep in sweeps:
        sc.doc.Objects.AddBrep(sweep)
    sc.doc.Views.Redraw()


def _addSweep_using_Sweep1(arcs, loft, bUseSweep1Dialog, bDebug):

    gArcs = [sc.doc.Objects.AddArc(arc) for arc in arcs]
    gLoft = sc.doc.Objects.AddBrep(loft)

    sc.doc.Views.Redraw()


    # The persistent selection may have caused the following to fail.
    ## Presuming index of edge remains the same since loft is always defined the same way.
    #gLoft = sc.doc.Objects.AddBrep(loft)
    #rdLoft = sc.doc.Objects.FindId(gLoft)
    #compIdx = rg.ComponentIndex(
    #        rg.ComponentIndexType.BrepEdge,
    #        index=1)
    ## TODO: Set syncHighlight to bDebug
    #rdLoft.SelectSubObject(
    #    componentIndex=compIdx,
    #    select=True,
    #    syncHighlight=True,
    #    persistentSelect=True)

    # So that user can see prompt to pick rail.
    #Rhino.RhinoApp.RunScript("_Echo", echo=False)

    script  = "_Sweep1 " if bUseSweep1Dialog else "_-Sweep1 "
    script += "_Pause "
    for gArc in gArcs:
        script += "_SelId {} ".format(gArc)
    script += "_Enter"
    #script += "_Enter" if bUseSweep1Dialog else "_EnterEnd"
    #script += " _Style=Freeform _ShapeBlending=Local _Simplify=None _Enter"
    script += " _Style=AlignWithSurface _Simplify=None _Enter"
    Rhino.RhinoApp.RunScript(script, echo=True)

    if not bDebug:
        [sc.doc.Objects.Delete(objectId=gArc, quiet=False) for gArc in gArcs]
        sc.doc.Objects.Delete(objectId=gLoft, quiet=False)

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()


def main():

    while True:
        rc = _createGeometryInteractively()
        if rc is None: return

        bApproximate_NotSweep1 = Opts.values['bApproximate_NotSweep1']
        bOnlyOneArc_NotAll = Opts.values['bOnlyOneArc_NotAll']
        bUseSweep1Dialog = Opts.values['bUseSweep1Dialog']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        arcs_perOffset, rgBs_Lofts, ncs_toOffset, gBrep_Picked = rc

        if bApproximate_NotSweep1:
            for arcs_thisOffset, loft, nc_toOffset in zip(arcs_perOffset, rgBs_Lofts, ncs_toOffset):
                _addFilletApproximation_RC(
                    arcs_thisOffset,
                    nc_toOffset,
                    loft)
            sc.doc.Objects.UnselectAll()
        else:
            if bDebug:
                sc.doc.Objects.UnselectAll()
                sc.doc.Views.Redraw()
                continue

            #print(sc.doc.Objects.Lock(objectId=gBrep_Picked, ignoreLayerMode=True))

            for arcs_thisOffset, loft in zip(arcs_perOffset, rgBs_Lofts):
                _addSweep_using_Sweep1(
                    [arcs_thisOffset[0]] if bOnlyOneArc_NotAll else arcs_thisOffset,
                    loft,
                    bUseSweep1Dialog,
                    bDebug)

            #print(sc.doc.Objects.Unlock(objectId=gBrep_Picked, ignoreLayerMode=True))



        sc.doc.Views.Redraw()


if __name__ == '__main__': main()