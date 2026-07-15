#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script is an alternative to _EndBulge.

TL;DR: The command should be used to really understand it. Compare with _EndBulge
and see.

Unlike _EndBulge, there is no dragging of control points in the graphics window.
Instead, number steppers in a dialog are used to specify:
    1. The tangent vector (p1 - p0) scale relative to its starting position.
    2. Where the geometry allows it, and after #1 is applied,
    the G2 (p2) tangential sliding scale from p2's starting position.
    3. Where the geometry allows it, and after #1 & #2 are applied,
    the G3 (p3) tangential sliding scale from p3's starting position.
Unlike _EndBulge, both the picked edge and the opposite edge can be modified.
Unlike _EndBulge, a curvature graph is automatically included.
Unlike _EndBulge, the continuities to maintain for the picked edge and the opposite
edge are explicitly defined and selectable by the user.

There are more options, but the command should be used to really understand it.

This script was created by Google Gemini 3.1 Pro based on the curve version.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
260712-14: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import Eto.Forms as ef

import spb_EndBulge_Kernel as ebk


def _replace_and_preserve_modes(doc, obj_id, new_geom):
    """Safely replaces document geometry without destroying active Zebra/Draft modes."""
    import System
    obj = doc.Objects.FindId(obj_id)
    if not obj: return False
    
    active_modes = []
    
    # RhinoObject lacks a "GetActiveModes" method. We must test known GUIDs explicitly.
    # We use hasattr/getattr to safely navigate known McNeel API spelling errors.
    potential_mode_names = [
        "RhinoZebraStripeAnalysisModeId",
        "RhinoEmapAnalysisModeId",
        "RhinoDraftAngleAnalysisModeId",
        "RhinoCurvatureColorAnalyisModeId",  # Native API typo
        "RhinoCurvatureColorAnalysisModeId"  # Fallback in case it is ever patched
    ]
    
    for name in potential_mode_names:
        if hasattr(Rhino.Display.VisualAnalysisMode, name):
            guid = getattr(Rhino.Display.VisualAnalysisMode, name)
            mode = Rhino.Display.VisualAnalysisMode.Find(guid)
            if mode and obj.InVisualAnalysisMode(mode):
                if mode not in active_modes:
                    active_modes.append(mode)
                
    rc = doc.Objects.Replace(obj_id, new_geom)
    
    if rc and active_modes:
        new_obj = doc.Objects.FindId(obj_id)
        if new_obj:
            for m in active_modes:
                new_obj.EnableVisualAnalysisMode(m, True)
                
    return rc


def extract_temp_curve(ns, direction, index):
    """Slices a specific U or V row of control points into a mathematically perfect 1D NurbsCurve."""
    if direction == 'U':
        nc = rg.NurbsCurve(3, ns.IsRational, ns.OrderU, ns.Points.CountU)
        for i in range(ns.KnotsU.Count): nc.Knots[i] = ns.KnotsU[i]
        for i in range(ns.Points.CountU):
            pt = ns.Points.GetControlPoint(i, index)
            nc.Points.SetPoint(i, pt.Location, pt.Weight)
        return nc
    else:
        nc = rg.NurbsCurve(3, ns.IsRational, ns.OrderV, ns.Points.CountV)
        for i in range(ns.KnotsV.Count): nc.Knots[i] = ns.KnotsV[i]
        for i in range(ns.Points.CountV):
            pt = ns.Points.GetControlPoint(index, i)
            nc.Points.SetPoint(i, pt.Location, pt.Weight)
        return nc


def get_curvature_isocurves(ns):
    """Extracts isocurves and their UV orientation for normal-aligned curvature graphs."""
    crvs_info = []
    tol = Rhino.RhinoMath.ZeroTolerance
    
    # U-direction curves (constant V parameter)
    v_dom = ns.Domain(1)
    v_params = [v_dom.Min, v_dom.Max]
    
    internal_v = [k for k in ns.KnotsV if (k > v_dom.Min + tol and k < v_dom.Max - tol)]
    internal_v = sorted(list(set(internal_v)))
    
    if not internal_v: v_params.append(v_dom.Mid)
    else: v_params.extend(internal_v)
        
    for v in v_params:
        c = ns.IsoCurve(0, v)
        if c: crvs_info.append((c, 0, v))
        
    # V-direction curves (constant U parameter)
    u_dom = ns.Domain(0)
    u_params = [u_dom.Min, u_dom.Max]
    
    internal_u = [k for k in ns.KnotsU if (k > u_dom.Min + tol and k < u_dom.Max - tol)]
    internal_u = sorted(list(set(internal_u)))
    
    if not internal_u: u_params.append(u_dom.Mid)
    else: u_params.extend(internal_u)
        
    for u in u_params:
        c = ns.IsoCurve(1, u)
        if c: crvs_info.append((c, 1, u))
        
    return crvs_info


def draw_surface_curvature_graph(display, ns, c, direction, const_param, scale, density, color):
    """
    Manually draws a curvature graph aligned strictly to the Surface Normal.
    Bypasses the curve's osculating plane to prevent 'tapered' or 'leaning' graphs.
    """
    
    # Convert 1e-6 cm into current document units for our minimum distance threshold
    unit_scale = Rhino.RhinoMath.UnitScale(Rhino.UnitSystem.Centimeters, sc.doc.ModelUnitSystem)
    min_dist = 1e-6 * unit_scale    
    
    # Skip completely collapsed curves (singular boundaries)
    if c.GetLength() < min_dist:
        return
        
    # --- RHINO's HIDDEN EXPONENTIAL SCALE FORMULA ---
    true_scale = 2.0 ** ((scale - 100.0) / 2.0)
    
    # --- DECOUPLED SPAN DENSITY LOGIC ---
    hair_steps = max(1, int(density) + 1)
    
    if int(density) == 0:
        # Special edge-case fallback: Rhino uses exactly 12 segments for density 0
        multiplier = 12
    else:
        # Rhino's target 76-segment formula for all densities >= 1
        # Integer ceiling division: (76 + hair_steps - 1) // hair_steps
        multiplier = (76 + hair_steps - 1) // hair_steps
        
    env_steps = hair_steps * multiplier
    
    hair_t_vals = []
    env_t_vals = []
    
    # Evaluate parameters strictly along the knot spans
    for i in range(c.SpanCount):
        dom = c.SpanDomain(i)
        for j in range(hair_steps):
            hair_t_vals.append(dom.Min + (dom.Length / float(hair_steps)) * j)
        for j in range(env_steps):
            env_t_vals.append(dom.Min + (dom.Length / float(env_steps)) * j)
            
    # Cap the arrays with the absolute end parameter
    hair_t_vals.append(c.Domain.Max)
    env_t_vals.append(c.Domain.Max)
    
    # 1. Evaluate and draw the high-resolution smooth envelope
    env_pts = []
    for t in env_t_vals:
        P = c.PointAt(t)
        cv = c.CurvatureAt(t)
        
        if direction == 0: norm = ns.NormalAt(t, const_param)
        else: norm = ns.NormalAt(const_param, t)
            
        if not norm.IsValid or norm.Length < min_dist or not cv.IsValid or cv.Length > 1e5:
            env_pts.append(P)
            continue
            
        kappa_n = cv * norm
        hair = norm * (kappa_n * true_scale * -1.0)
        env_pts.append(P + hair)
        
    if len(env_pts) > 1:
        display.DrawPolyline(env_pts, color, 1)

    # 2. Evaluate and draw the targeted density hairs
    for t in hair_t_vals:
        P = c.PointAt(t)
        cv = c.CurvatureAt(t)
        
        if direction == 0: norm = ns.NormalAt(t, const_param)
        else: norm = ns.NormalAt(const_param, t)
            
        if not norm.IsValid or norm.Length < min_dist or not cv.IsValid or cv.Length > 1e5:
            continue
            
        kappa_n = cv * norm
        hair = norm * (kappa_n * true_scale * -1.0)
        display.DrawLine(P, P + hair, color, 1)


def createSurface(ns_In, boundary, fScale_Picked=1.0, fSlideG2_Picked=0.0, fSlideG3_Picked=0.0, fScale_Opp=1.0, fSlideG2_Opp=0.0, fSlideG3_Opp=0.0, iG_Picked=2, iG_Opp=2, bDebug=False):
    """Loops through the orthogonal surface grid and executes the Kernel curve engine on each row."""
    ns_Out = ns_In.Duplicate()
    global_info = None

    if any(abs(_ - Rhino.RhinoMath.ZeroTolerance) != 1.0 for _ in (fScale_Picked, fScale_Opp)):
        pass
    elif any(abs(_ - Rhino.RhinoMath.ZeroTolerance) != 0.0 for _ in (fSlideG2_Picked, fSlideG3_Picked, fSlideG2_Opp, fSlideG3_Opp)):
        pass
    else:
        return None, "All scale and slide values result in no change to the geometry.", None

    if boundary in ('U0', 'U1'):
        iPickedEnd = 0 if boundary == 'U0' else 1
        for v in range(ns_In.Points.CountV):
            nc_temp = extract_temp_curve(ns_In, 'U', v)
            
            if iPickedEnd == 0:
                s0, g2_0, g3_0, i0 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked, iG_Picked
                s1, g2_1, g3_1, i1 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp, iG_Opp
            else:
                s1, g2_1, g3_1, i1 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked, iG_Picked
                s0, g2_0, g3_0, i0 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp, iG_Opp

            nc_res, sReport, info = ebk.createCurve(
                nc_In=nc_temp,
                fScale_T0=s0, fSlideG2_T0=g2_0, fSlideG3_T0=g3_0,
                fScale_T1=s1, fSlideG2_T1=g2_1, fSlideG3_T1=g3_1,
                iG_T0=i0, iG_T1=i1, iPickedEnd=iPickedEnd, bDebug=False
            )

            if nc_res is None: return None, sReport, None
            if v == 0: global_info = info 

            for u in range(ns_In.Points.CountU):
                pt = nc_res.Points[u]
                ns_Out.Points.SetControlPoint(u, v, rg.ControlPoint(pt.Location, pt.Weight))

    else:
        iPickedEnd = 0 if boundary == 'V0' else 1
        for u in range(ns_In.Points.CountU):
            nc_temp = extract_temp_curve(ns_In, 'V', u)
            
            if iPickedEnd == 0:
                s0, g2_0, g3_0, i0 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked, iG_Picked
                s1, g2_1, g3_1, i1 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp, iG_Opp
            else:
                s1, g2_1, g3_1, i1 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked, iG_Picked
                s0, g2_0, g3_0, i0 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp, iG_Opp

            nc_res, sReport, info = ebk.createCurve(
                nc_In=nc_temp,
                fScale_T0=s0, fSlideG2_T0=g2_0, fSlideG3_T0=g3_0,
                fScale_T1=s1, fSlideG2_T1=g2_1, fSlideG3_T1=g3_1,
                iG_T0=i0, iG_T1=i1, iPickedEnd=iPickedEnd, bDebug=False
            )

            if nc_res is None: return None, sReport, None
            if u == 0: global_info = info 

            for v in range(ns_In.Points.CountV):
                pt = nc_res.Points[v]
                ns_Out.Points.SetControlPoint(u, v, rg.ControlPoint(pt.Location, pt.Weight))

    return ns_Out, None, global_info


class SrfPreviewConduit(Rhino.Display.DisplayConduit):
    def __init__(self):
        super(SrfPreviewConduit, self).__init__()
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.ns = None
        self.cg_curves = []
        self.show_graph = ebk.Opts.values['bShowGraph']
        self.show_polygon = ebk.Opts.values['bShowPolygon']
        self.show_geom = ebk.Opts.values['bShowGeom']
        self.graph_scale = ebk.Opts.values['iGraphScale']
        self.graph_density = ebk.Opts.values['iGraphDensity']
        self.is_swapped = False

    def CalculateBoundingBox(self, e):
        if self.ns: 
            bbox = self.ns.GetBoundingBox(False)
            bbox.Inflate(sc.doc.ModelAbsoluteTolerance * 100)
            e.IncludeBoundingBox(bbox)

    def PostDrawObjects(self, e):
        if not self.ns: return

        if self.show_polygon:
            # Draw U-direction CV net lines
            for v in range(self.ns.Points.CountV):
                pts = [self.ns.Points.GetControlPoint(u, v).Location for u in range(self.ns.Points.CountU)]
                e.Display.DrawPolyline(pts, self.color, 1)
            
            # Draw V-direction CV net lines
            for u in range(self.ns.Points.CountU):
                pts = [self.ns.Points.GetControlPoint(u, v).Location for v in range(self.ns.Points.CountV)]
                e.Display.DrawPolyline(pts, self.color, 1)

            # Draw CP Dots
            all_pts = [self.ns.Points.GetControlPoint(u, v).Location for u in range(self.ns.Points.CountU) for v in range(self.ns.Points.CountV)]
            e.Display.DrawPoints(all_pts, Rhino.Display.PointStyle.Simple, 3, self.color)

        # Draw Base Isocurves and Custom Normal-Aligned Curvature Graphs
        for c, direction, const_param in self.cg_curves:
            if self.show_geom and not self.is_swapped:
                e.Display.DrawCurve(c, self.color, 1)
            if self.show_graph:
                draw_surface_curvature_graph(e.Display, self.ns, c, direction, const_param, self.graph_scale, self.graph_density, self.color)


class SrfEtoDialog(ebk.EtoDialog):
    """Subclasses the Kernel Dialog to intercept Surface selections and override the preview routine."""
    def __init__(self, objref_In):
        self.is_surface = True  
        self.Title = "EndBulge Surface (Side Bulge) by SPB"
        self.objref_In = objref_In
        self.dialog_ok = False

        edge = objref_In.Edge()
        face = edge.Brep.Faces[edge.AdjacentFaces()[0]]
        self.ns_In = face.ToNurbsSurface()

        t_mid = edge.Domain.Mid
        pt_mid = edge.PointAt(t_mid)
        bSuccess, u, v = face.ClosestPoint(pt_mid)
        
        domU, domV = face.Domain(0), face.Domain(1)
        dU0, dU1 = abs(u - domU.Min), abs(domU.Max - u)
        dV0, dV1 = abs(v - domV.Min), abs(domV.Max - v)
        min_d = min(dU0, dU1, dV0, dV1)

        if min_d == dU0: self.boundary = 'U0'
        elif min_d == dU1: self.boundary = 'U1'
        elif min_d == dV0: self.boundary = 'V0'
        else: self.boundary = 'V1'

        if self.boundary in ('U0', 'U1'):
            self.nc_In = extract_temp_curve(self.ns_In, 'U', 0)
        else:
            self.nc_In = extract_temp_curve(self.ns_In, 'V', 0)
            
        # GUARANTEED BACKUP: Explicitly pull the Brep geometry to prevent extracting the 1D Edge Curve
        self.original_geom = objref_In.Brep().Duplicate()

        self._exact_scale_picked = ebk.Opts.values['fScale_Picked']
        self._exact_scale_opp = ebk.Opts.values['fScale_Opp']
        self._auto_updating = False
        self._auto_updating_slider = False  
        self.active_stepper_key = None      

        self.create_controls()
        self.setup_layout()
        self.OnLinkedModeChanged(None, None)

        self.debounce_timer = ef.UITimer()
        self.debounce_timer.Interval = 0.2
        self.debounce_timer.Elapsed += self.OnDebounceTimerElapsed

        self.LoadComplete += self.OnFormLoadComplete
        self.Closed += self.OnFormClosed

    def UpdatePreview(self):
        if not hasattr(self, 'conduit') or self.conduit is None: return

        fScale_Picked = self.ParseToFloat(self.textBoxes['fScale_Picked'].Text)

        fSlideG2_Picked = self.ParseToFloat(self.textBoxes['fSlideG2_Picked'].Text)
        fSlideG2_Picked = fSlideG2_Picked if fSlideG2_Picked is not None else 0.0

        fSlideG3_Picked = self.ParseToFloat(self.textBoxes['fSlideG3_Picked'].Text)
        fSlideG3_Picked = fSlideG3_Picked if fSlideG3_Picked is not None else 0.0

        fScale_Opp = self.ParseToFloat(self.textBoxes['fScale_Opp'].Text)

        fSlideG2_Opp = self.ParseToFloat(self.textBoxes['fSlideG2_Opp'].Text)
        fSlideG2_Opp = fSlideG2_Opp if fSlideG2_Opp is not None else 0.0
        
        fSlideG3_Opp = self.ParseToFloat(self.textBoxes['fSlideG3_Opp'].Text)
        fSlideG3_Opp = fSlideG3_Opp if fSlideG3_Opp is not None else 0.0

        if (fScale_Picked is None or fScale_Picked <= Rhino.RhinoMath.ZeroTolerance or
            fScale_Opp is None or fScale_Opp <= Rhino.RhinoMath.ZeroTolerance):
            self.conduit.ns = None
            sc.doc.Views.Redraw()
            return

        idxCont_Picked = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        idxCont_Opp = self.radioButtonLists['idxCont_Opp'].SelectedIndex

        ns_Res, sReport, info = createSurface(
            self.ns_In, self.boundary,
            fScale_Picked, fSlideG2_Picked, fSlideG3_Picked,
            fScale_Opp, fSlideG2_Opp, fSlideG3_Opp,
            idxCont_Picked - 1, idxCont_Opp - 1, False
        )

        if info is not None:
            actual_T0, actual_T1, bOverlap = info
            
            if self.boundary in ('U1', 'V1'):
                actual_Picked, actual_Opp = actual_T1, actual_T0
            else:
                actual_Picked, actual_Opp = actual_T0, actual_T1
            
            if not self._auto_updating:
                self._auto_updating = True
                changed = False
                if actual_Picked != idxCont_Picked - 1:
                    self.radioButtonLists['idxCont_Picked'].SelectedIndex = actual_Picked + 1
                    changed = True
                if actual_Opp != idxCont_Opp - 1:
                    self.radioButtonLists['idxCont_Opp'].SelectedIndex = actual_Opp + 1
                    changed = True
                if (bOverlap or changed) and self.radioButtonLists['bLinkedEnds'].SelectedIndex == 1:
                    self.radioButtonLists['bLinkedEnds'].SelectedIndex = 0
                    self.OnLinkedModeChanged(None, None)
                    changed = True
                if changed: self.UpdateControlStates()
                self._auto_updating = False

        self.conduit.ns = self.ns_In.Duplicate() if ns_Res is None else ns_Res
        
        if self.conduit.ns:
            self.conduit.cg_curves = get_curvature_isocurves(self.conduit.ns)
        else:
            self.conduit.cg_curves = []
            
        # Hide the document object while actively sliding (Conduit takes over)
        sc.doc.Objects.Hide(self.objref_In.ObjectId, True)
        self.conduit.is_swapped = False
        sc.doc.Views.Redraw()

        self.debounce_timer.Stop()
        self.debounce_timer.Start()

    def OnDisplayCheckedChanged(self, sender, e):
        # 1. Run the Kernel's native logic to update the conduit flags
        super(SrfEtoDialog, self).OnDisplayCheckedChanged(sender, e)
        
        # 2. If the document object is currently swapped in, explicitly show/hide it
        if hasattr(self, 'conduit') and self.conduit and self.conduit.is_swapped:
            if self.checkBoxes['bShowGeom'].Checked:
                sc.doc.Objects.Show(self.objref_In.ObjectId, True)
            else:
                sc.doc.Objects.Hide(self.objref_In.ObjectId, True)
            
            sc.doc.Views.Redraw()

    def OnDebounceTimerElapsed(self, sender, e):
        self.debounce_timer.Stop()
        if hasattr(self.conduit, 'ns') and self.conduit.ns:
            new_brep = self.conduit.ns.ToBrep()
            if new_brep:
                # The Golden Rule: Rhino refuses to replace hidden objects.
                # We must unhide the original geometry right BEFORE we swap it.
                sc.doc.Objects.Show(self.objref_In.ObjectId, True)
                
                # Now the replace will succeed, and Zebra is preserved
                _replace_and_preserve_modes(sc.doc, self.objref_In.ObjectId, new_brep)
                
                # If the user doesn't want to see the surface, re-hide it instantly
                if not self.checkBoxes['bShowGeom'].Checked:
                    sc.doc.Objects.Hide(self.objref_In.ObjectId, True)
                
                self.conduit.is_swapped = True
                sc.doc.Views.Redraw()

    def OnFormClosed(self, sender, e):
        # We leave the document state alone here; main() will clean it up inside the UndoRecord!
        super(SrfEtoDialog, self).OnFormClosed(sender, e)


def getInput_CLI():
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Pick edge of an untrimmed surface")
    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.EdgeCurve
    go.DisablePreSelect()
    
    # Custom filter to ensure the edge belongs to an untrimmed surface
    def edge_filter(rhObj, geom, compIdx):
        if isinstance(geom, rg.BrepEdge):
            faces = geom.AdjacentFaces()
            if faces.Count > 0:
                face = geom.Brep.Faces[faces[0]]
                if face.IsSurface:
                    return True
        return False
    go.SetCustomGeometryFilter(edge_filter)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = ebk.Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()
        
        addOption('bGUI')
        if not ebk.Opts.values['bGUI']:
            addOption('bLinkedEnds')
            addOption('fScale_Picked')
            addOption('fSlideG2_Picked')
            addOption('fSlideG3_Picked')
            addOption('fScale_Opp')
            addOption('fSlideG2_Opp')
            addOption('fSlideG3_Opp')
            addOption('idxCont_Picked')
            addOption('idxCont_Opp')
            addOption('bDeleteInput')
        
        res = go.Get()
        if res == ri.GetResult.Cancel: return None
        if res == ri.GetResult.Object: return go.Object(0)

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                ebk.Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def main():
    objref_In = getInput_CLI()
    if objref_In is None: return

    if not ebk.Opts.values['bGUI']:
        print("CLI processing for surfaces is under construction. Please use the GUI mode.")
        return
    
    Rhino.RhinoApp.SetCommandPromptMessage("Continuing in dialog...")
    key = 'conduit_srf'
    stickyKey = '{}({})'.format(key, __file__)
    if stickyKey in sc.sticky:
        sc.sticky[stickyKey].Enabled = False

    parent = Rhino.UI.RhinoEtoApp.MainWindowForDocument(sc.doc)
    dialog = SrfEtoDialog(objref_In)
    dialog.conduit = SrfPreviewConduit()
    sc.sticky[stickyKey] = dialog.conduit

    dialog.UpdatePreview()
    dialog.conduit.Enabled = True
    sc.doc.Views.Redraw()

    # START UNDO RECORD: Groups all live-swapping slider drags into one clean Undo step
    undo_sn = sc.doc.BeginUndoRecord("EndBulge Srf")

    try:
        Rhino.UI.EtoExtensions.ShowSemiModal(dialog, sc.doc, parent)

        # POST-DIALOG LOGIC
        if dialog.dialog_ok and dialog.conduit.ns:
            if dialog.conduit.ns.EpsilonEquals(objref_In.Face().UnderlyingSurface(), epsilon=Rhino.RhinoMath.ZeroTolerance):
                print("Resultant surface is the same as input surface. No changes were made to the document.")
                _replace_and_preserve_modes(sc.doc, objref_In.ObjectId, dialog.original_geom)
                sc.doc.Objects.Show(objref_In.ObjectId, True)
            elif ebk.Opts.values['bDeleteInput']:
                # Replace final and keep Zebra active
                _replace_and_preserve_modes(sc.doc, objref_In.ObjectId, dialog.conduit.ns.ToBrep())
                sc.doc.Objects.Show(objref_In.ObjectId, True)
                print("Replaced surface.")
            else:
                # Restore original with Zebra, and add the new one cleanly
                _replace_and_preserve_modes(sc.doc, objref_In.ObjectId, dialog.original_geom)
                sc.doc.Objects.Show(objref_In.ObjectId, True)
                sc.doc.Objects.AddSurface(dialog.conduit.ns)
                print("Surface was added.")
        else:
            # User Cancelled: Safely restore original geometry with Zebra intact
            _replace_and_preserve_modes(sc.doc, objref_In.ObjectId, dialog.original_geom)
            sc.doc.Objects.Show(objref_In.ObjectId, True)

    except Exception as e:
        # CRASH RECOVERY: Ensure the object is recovered and unhidden if the script throws an error
        print("Script Error Encountered: {}".format(e))
        _replace_and_preserve_modes(sc.doc, objref_In.ObjectId, dialog.original_geom)
        sc.doc.Objects.Show(objref_In.ObjectId, True)
        
    finally:
        # Guarantee cleanup regardless of OK, Cancel, or Exception
        dialog.conduit.Enabled = False
        sc.doc.EndUndoRecord(undo_sn)
        sc.doc.Views.Redraw()

if __name__ == '__main__': main()
