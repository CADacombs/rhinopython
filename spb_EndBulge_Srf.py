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
260712: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc
import sys
import os

import spb_EndBulge_Kernel as ebk

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
    samples = int(50 + (density * 2)) 
    dom = c.Domain
    step = dom.Length / samples
    
    # --- RHINO's HIDDEN EXPONENTIAL SCALE FORMULA ---
    # Converts the linear integer UI slider into the massive exponential 
    # multiplier required for visual hair scaling. (e.g., 110 = 32x, 120 = 1024x)
    true_scale = 2.0 ** ((scale - 100.0) / 2.0)
    
    tips = []
    for i in range(samples + 1):
        t = dom.Min + i * step
        if i == samples: t = dom.Max
        
        P = c.PointAt(t)
        cv = c.CurvatureAt(t)
        
        if direction == 0:
            norm = ns.NormalAt(t, const_param)
        else:
            norm = ns.NormalAt(const_param, t)
            
        if not norm.IsValid or norm.Length < 1e-6 or not cv.IsValid:
            tips.append(P)
            continue
            
        # Project 3D curve curvature onto Surface Normal (Dot Product)
        kappa_n = cv * norm
        
        # Multiply by our exponential true_scale, and -1.0 to invert to Rhino's native side
        hair = norm * (kappa_n * true_scale * -1.0)
        tip = P + hair
        
        tips.append(tip)
        display.DrawLine(P, tip, color, 1)
        
    if len(tips) > 1:
        display.DrawPolyline(tips, color, 1)


def createSurface(ns_In, boundary, fScale_Picked=1.0, fFullG2_Picked=1.0, fFullG3_Picked=1.0, fScale_Opp=1.0, fFullG2_Opp=1.0, fFullG3_Opp=1.0, iG_Picked=2, iG_Opp=2, bDebug=False):
    """Loops through the orthogonal surface grid and executes the Kernel curve engine on each row."""
    ns_Out = ns_In.Duplicate()
    global_info = None

    if boundary in ('U0', 'U1'):
        iPickedEnd = 0 if boundary == 'U0' else 1
        for v in range(ns_In.Points.CountV):
            nc_temp = extract_temp_curve(ns_In, 'U', v)
            
            # Map the generic picked/opp variables to the strict T0/T1 curve ends
            if iPickedEnd == 0:
                s0, g2_0, g3_0, i0 = fScale_Picked, fFullG2_Picked, fFullG3_Picked, iG_Picked
                s1, g2_1, g3_1, i1 = fScale_Opp, fFullG2_Opp, fFullG3_Opp, iG_Opp
            else:
                s1, g2_1, g3_1, i1 = fScale_Picked, fFullG2_Picked, fFullG3_Picked, iG_Picked
                s0, g2_0, g3_0, i0 = fScale_Opp, fFullG2_Opp, fFullG3_Opp, iG_Opp

            nc_res, sReport, info = ebk.createCurve(
                nc_In=nc_temp,
                fScale_T0=s0, fFullG2_T0=g2_0, fFullG3_T0=g3_0,
                fScale_T1=s1, fFullG2_T1=g2_1, fFullG3_T1=g3_1,
                iG_T0=i0, iG_T1=i1, iPickedEnd=iPickedEnd, bDebug=False
            )

            if nc_res is None: return None, sReport, None
            if v == 0: global_info = info # Grab overlap info from the first row

            for u in range(ns_In.Points.CountU):
                pt = nc_res.Points[u]
                ns_Out.Points.SetControlPoint(u, v, rg.ControlPoint(pt.Location, pt.Weight))

    else:
        iPickedEnd = 0 if boundary == 'V0' else 1
        for u in range(ns_In.Points.CountU):
            nc_temp = extract_temp_curve(ns_In, 'V', u)
            
            if iPickedEnd == 0:
                s0, g2_0, g3_0, i0 = fScale_Picked, fFullG2_Picked, fFullG3_Picked, iG_Picked
                s1, g2_1, g3_1, i1 = fScale_Opp, fFullG2_Opp, fFullG3_Opp, iG_Opp
            else:
                s1, g2_1, g3_1, i1 = fScale_Picked, fFullG2_Picked, fFullG3_Picked, iG_Picked
                s0, g2_0, g3_0, i0 = fScale_Opp, fFullG2_Opp, fFullG3_Opp, iG_Opp

            nc_res, sReport, info = ebk.createCurve(
                nc_In=nc_temp,
                fScale_T0=s0, fFullG2_T0=g2_0, fFullG3_T0=g3_0,
                fScale_T1=s1, fFullG2_T1=g2_1, fFullG3_T1=g3_1,
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
        self.graph_scale = ebk.Opts.values['iGraphScale']
        self.graph_density = ebk.Opts.values['iGraphDensity']

    def CalculateBoundingBox(self, e):
        if self.ns: 
            bbox = self.ns.GetBoundingBox(False)
            bbox.Inflate(sc.doc.ModelAbsoluteTolerance * 100)
            e.IncludeBoundingBox(bbox)

    def PostDrawObjects(self, e):
        if not self.ns: return

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

        # Draw Custom Normal-Aligned Curvature Graphs
        for c, direction, const_param in self.cg_curves:
            draw_surface_curvature_graph(e.Display, self.ns, c, direction, const_param, self.graph_scale, self.graph_density, self.color)


class SrfEtoDialog(ebk.EtoDialog):
    """Subclasses the Kernel Dialog to intercept Surface selections and override the preview routine."""
    def __init__(self, objref_In):
        self.Title = "SideBulge Surface"
        self.objref_In = objref_In
        self.dialog_ok = False

        face = objref_In.Face()
        if face is None: return
        self.ns_In = face.ToNurbsSurface()

        # Identify which boundary the user clicked closest to
        bSuccess, u, v = face.ClosestPoint(objref_In.SelectionPoint())
        domU, domV = face.Domain(0), face.Domain(1)
        
        dU0, dU1 = abs(u - domU.Min), abs(domU.Max - u)
        dV0, dV1 = abs(v - domV.Min), abs(domV.Max - v)
        min_d = min(dU0, dU1, dV0, dV1)

        if min_d == dU0: self.boundary = 'U0'
        elif min_d == dU1: self.boundary = 'U1'
        elif min_d == dV0: self.boundary = 'V0'
        else: self.boundary = 'V1'

        # Supply a temporary 1D curve to the underlying Kernel so UI calculations (like point count) still work natively
        if self.boundary in ('U0', 'U1'):
            self.nc_In = extract_temp_curve(self.ns_In, 'U', 0)
        else:
            self.nc_In = extract_temp_curve(self.ns_In, 'V', 0)

        self._exact_scale_picked = ebk.Opts.values['fScale_Picked']
        self._exact_scale_opp = ebk.Opts.values['fScale_Opp']
        self._auto_updating = False

        self.create_controls()
        self.setup_layout()
        self.OnLinkedModeChanged(None, None)

    def UpdatePreview(self):
        if not hasattr(self, 'conduit') or self.conduit is None: return

        fScale_Picked = self.ParseToFloat(self.textBoxes['fScale_Picked'].Text)
        fFullG2_Picked = self.numericSteppers['fFullG2_Picked'].Value
        fFullG3_Picked = self.numericSteppers['fFullG3_Picked'].Value
        fScale_Opp = self.ParseToFloat(self.textBoxes['fScale_Opp'].Text)
        fFullG2_Opp = self.numericSteppers['fFullG2_Opp'].Value
        fFullG3_Opp = self.numericSteppers['fFullG3_Opp'].Value

        if (fScale_Picked is None or fScale_Picked <= Rhino.RhinoMath.ZeroTolerance or
            fScale_Opp is None or fScale_Opp <= Rhino.RhinoMath.ZeroTolerance):
            self.conduit.ns = None
            sc.doc.Views.Redraw()
            return

        idxCont_Picked = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        idxCont_Opp = self.radioButtonLists['idxCont_Opp'].SelectedIndex

        ns_Res, sReport, info = createSurface(
            self.ns_In, self.boundary,
            fScale_Picked, fFullG2_Picked, fFullG3_Picked,
            fScale_Opp, fFullG2_Opp, fFullG3_Opp,
            idxCont_Picked - 1, idxCont_Opp - 1, False
        )

        # Dynamic UI feedback for continuity downgrades
        if info is not None:
            actual_T0, actual_T1, bOverlap = info
            
            # Map T0/T1 back to Picked/Opp based on boundary
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
        
        # Extract isocurves for the curvature graph preview
        if self.conduit.ns:
            self.conduit.cg_curves = get_curvature_isocurves(self.conduit.ns)
        else:
            self.conduit.cg_curves = []
            
        sc.doc.Views.Redraw()


def getInput_CLI():
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Pick surface near an edge")
    go.GeometryFilter = rd.ObjectType.Surface | rd.ObjectType.Brep
    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = ebk.Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()
        
        addOption('bGUI')
        if not ebk.Opts.values['bGUI']:
            addOption('bLinkedEnds')
            addOption('fScale_Picked')
            addOption('fFullG2_Picked')
            addOption('fFullG3_Picked')
            addOption('fScale_Opp')
            addOption('fFullG2_Opp')
            addOption('fFullG3_Opp')
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
        # Future Expansion: Add CLI-only surface processor logic here.
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

    Rhino.UI.EtoExtensions.ShowSemiModal(dialog, sc.doc, parent)

    dialog.conduit.Enabled = False
    sc.doc.Views.Redraw()

    if dialog.dialog_ok and dialog.conduit.ns:
        if ebk.Opts.values['bDeleteInput']:
            sc.doc.Objects.Replace(objref_In.ObjectId, dialog.conduit.ns.ToBrep())
            print("Replaced surface.")
        else:
            sc.doc.Objects.AddSurface(dialog.conduit.ns)
            print("Surface was added.")
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()