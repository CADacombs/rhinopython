#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
================================================================================
CRITICAL CORE KERNEL FILE - DO NOT DELETE OR RUN DIRECTLY
================================================================================
This file houses the shared mathematical solvers, geometric point-allocation 
engines, and Eto UI frameworks for the spb_EndBulge tool suite.

DEPENDENCY NOTE: 
This script acts as a shared library module. It cannot be executed directly and 
will throw a NameError/ImportError if run on its own. It must remain in the same 
directory as its calling command-wrapper scripts.

DEVELOPER NOTES:
To utilize this kernel in a new tool, import it as a module and interface 
directly with the EtoDialog class or the createCurve geometry engine. The solvers 
are structurally decoupled from direct object selection loops to ensure multi-
geometry flexibility.
================================================================================

This script was partially developed by Google Gemini 3.1 Pro based on another script.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
260712-16: Created by extracting refactored code from another script.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import Eto.Drawing as ed
import Eto.Forms as ef
from Rhino.UI import RhinoEtoApp, EtoExtensions


class Opts:
    keys = []
    values = {}
    offValues = {}
    onValues = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'bGUI'; keys.append(key)
    values[key] = True
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'idxCont_Picked'; keys.append(key)
    listValues[key] = 'None', 'G0', 'G1', 'G2', 'G3'
    values[key] = 3
    names[key] = 'MaintainPicked'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'idxCont_Opp'; keys.append(key)
    listValues[key] = 'None', 'G0', 'G1', 'G2', 'G3'
    values[key] = 3
    names[key] = 'MaintainOpp'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLinkedEnds'; keys.append(key)
    values[key] = True
    names[key] = 'AdjustEnds'
    offValues[key] = 'Independent'
    onValues[key] = 'Linked'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fIncrement'; keys.append(key) 
    values[key] = 0.05
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iSliderSteps'; keys.append(key)
    listValues[key] = '5', '10', '20', '50', '100', '1000'
    values[key] = 1  # Defaults to index 1 ('10')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fScale_Picked'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSlideG2_Picked'; keys.append(key)
    values[key] = 0.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSlideG3_Picked'; keys.append(key)
    values[key] = 0.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fScale_Opp'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSlideG2_Opp'; keys.append(key)
    values[key] = 0.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSlideG3_Opp'; keys.append(key)
    values[key] = 0.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bShowPolygon'; keys.append(key)
    values[key] = True
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bShowGeom'; keys.append(key)
    values[key] = True
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bShowGraph'; keys.append(key)
    values[key] = True
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iGraphScale'; keys.append(key)
    values[key] = 100
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iGraphDensity'; keys.append(key)
    values[key] = 1
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteInput'; keys.append(key)
    values[key] = True
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    offValues[key] = 'No'
    onValues[key] = 'Yes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    for key in keys:
        if key not in names:
            names[key] = key[1:].replace('_', '')

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
                idxOpt = go.AddOptionToggle(cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))
        return idxOpt

    @classmethod
    def setValue(cls, key, idxList=None):
        if key == 'fScale':
            if cls.riOpts[key].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                print("Invalid input for scale value.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue == 1:
                print("A scale of 1 will not modify the curve. Scale was not modified.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def canMaintainG3(nc, bEvalT1End):
    """Evaluates if a NURBS curve can mathematically maintain G3 continuity."""
    if nc is None or not isinstance(nc, rg.NurbsCurve):
        return False
    if nc.Points.Count < 4:
        return False
    if nc.SpanCount == 1:
        return True

    knots = nc.Knots
    degree = nc.Degree
    iKnot = knots.Count - degree - 1 if bEvalT1End else degree
    return knots.KnotMultiplicity(iKnot) >= 3


def createCurve(nc_In, fScale_T0=1.0, fSlideG2_T0=0.0, fSlideG3_T0=0.0, fScale_T1=1.0, fSlideG2_T1=0.0, fSlideG3_T1=0.0, iG_T0=2, iG_T1=2, iPickedEnd=0, bDebug=False):
    if iG_T0 is None and iG_T1 is None: return None, "Both continuity inputs are None.", None
    if nc_In.IsPeriodic: return None, "Input curve is periodic.", None
    if not isinstance(nc_In, rg.NurbsCurve): return None, "Input curve is a {}".format(nc_In.GetType().Name), None
    if any(abs(_ - Rhino.RhinoMath.ZeroTolerance) != 1.0 for _ in (fScale_T0, fScale_T1)):
        pass
    elif any(abs(_ - Rhino.RhinoMath.ZeroTolerance) != 0.0 for _ in (fSlideG2_T0, fSlideG3_T0, fSlideG2_T1, fSlideG3_T1)):
        pass
    else:
        return None, "All scale and slide values result in no change to the geometry.", None

    # --- POINT ALLOCATION ENGINE ---
    N = nc_In.Points.Count
    req_0 = max(0, iG_T0 + 1)
    req_1 = max(0, iG_T1 + 1)
    
    bOverlap = False
    if req_0 + req_1 > N:
        bOverlap = True
        if req_0 > req_1:
            alloc_0 = min(req_0, N)
            alloc_1 = N - alloc_0
        elif req_1 > req_0:
            alloc_1 = min(req_1, N)
            alloc_0 = N - alloc_1
        else:
            if iPickedEnd == 0:
                alloc_0 = min(req_0, N)
                alloc_1 = N - alloc_0
            else:
                alloc_1 = min(req_1, N)
                alloc_0 = N - alloc_1
    else:
        alloc_0 = req_0
        alloc_1 = req_1

    max_mod_T0 = alloc_0 - 1
    max_mod_T1 = alloc_1 - 1
    
    free = N - alloc_0 - alloc_1
    if free > 0:
        half = free // 2
        extra = free % 2
        scale_limit_T0 = alloc_0 + half
        scale_limit_T1 = alloc_1 + half
        if extra:
            if iPickedEnd == 0: scale_limit_T0 += 1
            else: scale_limit_T1 += 1
    else:
        scale_limit_T0 = alloc_0
        scale_limit_T1 = alloc_1

    if bDebug and bOverlap:
        print("Overlap detected. T0 continuity capped at G{}, T1 capped at G{}.".format(max_mod_T0, max_mod_T1))

    # --- BASELINE CHECK ---
    def is_baseline(s, g2, g3, scale_limit):
        if scale_limit < 2: return True 
        b = abs(s - 1.0) <= Rhino.RhinoMath.ZeroTolerance
        if scale_limit > 2: b = b and (abs(g2) <= Rhino.RhinoMath.ZeroTolerance)
        if scale_limit > 3: b = b and (abs(g3) <= Rhino.RhinoMath.ZeroTolerance)
        return b

    base_T0 = is_baseline(fScale_T0, fSlideG2_T0, fSlideG3_T0, scale_limit_T0)
    base_T1 = is_baseline(fScale_T1, fSlideG2_T1, fSlideG3_T1, scale_limit_T1)

    if base_T0 and base_T1:
        return None, "Input parameters do not lead to modification of the geometry.", (max_mod_T0, max_mod_T1, bOverlap)

    pts_Prime = [pt.Location for pt in nc_In.Points]

    # Convert 1e-6 cm into current document units for our minimum distance threshold
    unit_scale = Rhino.RhinoMath.UnitScale(Rhino.UnitSystem.Centimeters, sc.doc.ModelUnitSystem)
    min_dist = 1e-6 * unit_scale

    # ----------------------------------------------------
    # SCALE T0 END
    # ----------------------------------------------------
    if not base_T0:
        p0 = nc_In.Points[0].Location
        xform_T0 = rg.Transform.Scale(p0, fScale_T0)
        
        for i in range(1, scale_limit_T0):
            pt_p = rg.Point3d(nc_In.Points[i].Location)
            pt_p.Transform(xform_T0)
            pts_Prime[i] = pt_p

        if scale_limit_T0 > 1:
            p1 = nc_In.Points[1].Location
            p1p = pts_Prime[1]
            slide_vec = p1p - p0
            orig_len_T0 = (p1 - p0).Length

            if scale_limit_T0 > 2:
                p2 = nc_In.Points[2].Location
                # Bypass tangent division if p0 and p1 are stacked (Singularity/Pole)
                if max_mod_T0 >= 2 and orig_len_T0 > min_dist:
                    m2 = ((p1p - p0).Length / orig_len_T0)**2.0
                    p2p_base = 2.0*p1p - p0 + m2*(-2.0*p1 + p2 + p0)
                else:
                    p2p_base = pts_Prime[2]
                    
                p2_slide = slide_vec * fSlideG2_T0
                pts_Prime[2] = p2p_base + p2_slide

                if scale_limit_T0 > 3:
                    p3 = nc_In.Points[3].Location
                    if max_mod_T0 >= 3 and orig_len_T0 > min_dist:
                        m3 = ((p1p - p0).Length / orig_len_T0)**3.0
                        p3p_base = 3.0*p2p_base - 3.0*p1p + p0 + m3*(p3 - 3.0*p2 + 3.0*p1 - p0)
                        p3_comp = 3.0 * p2_slide
                    else:
                        p3p_base = pts_Prime[3]
                        p3_comp = rg.Vector3d.Zero
                        
                    p3_slide = slide_vec * fSlideG3_T0
                    pts_Prime[3] = p3p_base + p3_comp + p3_slide

    # ----------------------------------------------------
    # SCALE T1 END
    # ----------------------------------------------------
    if not base_T1:
        last = N - 1
        p0 = nc_In.Points[last].Location
        xform_T1 = rg.Transform.Scale(p0, fScale_T1)
        
        for i in range(1, scale_limit_T1):
            idx = last - i
            pt_p = rg.Point3d(nc_In.Points[idx].Location)
            pt_p.Transform(xform_T1)
            pts_Prime[idx] = pt_p

        if scale_limit_T1 > 1:
            p1 = nc_In.Points[last-1].Location
            p1p = pts_Prime[last-1]
            slide_vec = p1p - p0
            orig_len_T1 = (p1 - p0).Length

            if scale_limit_T1 > 2:
                p2 = nc_In.Points[last-2].Location
                if max_mod_T1 >= 2 and orig_len_T1 > min_dist:
                    m2 = ((p1p - p0).Length / orig_len_T1)**2.0
                    p2p_base = 2.0*p1p - p0 + m2*(-2.0*p1 + p2 + p0)
                else:
                    p2p_base = pts_Prime[last-2]
                    
                p2_slide = slide_vec * fSlideG2_T1
                pts_Prime[last-2] = p2p_base + p2_slide

                if scale_limit_T1 > 3:
                    p3 = nc_In.Points[last-3].Location
                    if max_mod_T1 >= 3 and orig_len_T1 > min_dist:
                        m3 = ((p1p - p0).Length / orig_len_T1)**3.0
                        p3p_base = 3.0*p2p_base - 3.0*p1p + p0 + m3*(p3 - 3.0*p2 + 3.0*p1 - p0)
                        p3_comp = 3.0 * p2_slide
                    else:
                        p3p_base = pts_Prime[last-3]
                        p3_comp = rg.Vector3d.Zero
                        
                    p3_slide = slide_vec * fSlideG3_T1
                    pts_Prime[last-3] = p3p_base + p3_comp + p3_slide

    # Enforce minimum distance (but ignore inherently stacked singularity points)
    for i in range(len(pts_Prime) - 1):
        if pts_Prime[i].DistanceTo(pts_Prime[i+1]) < min_dist:
            orig_dist = nc_In.Points[i].Location.DistanceTo(nc_In.Points[i+1].Location)
            if orig_dist >= min_dist:
                sReport = "Minimum control point distance (1e-6 cm) violated. Is Scale too small?"
                Rhino.RhinoApp.SetCommandPromptMessage(sReport)
                if bDebug: print(sReport)
                return None, sReport, None

    nc_Out = nc_In.Duplicate()
    for i in range(nc_Out.Points.Count):
        nc_Out.Points.SetPoint(index=i, point=pts_Prime[i], weight=nc_In.Points.GetWeight(i))

    #if bDebug:
    #    for i in range(len(pts_Prime)):
    #        rgDot = rg.TextDot("{}".format(i), pts_Prime[i])
    #        rgDot.FontHeight = 11
    #        sc.doc.Objects.AddTextDot(rgDot)
    #    sc.doc.Views.Redraw()

    return nc_Out, None, (max_mod_T0, max_mod_T1, bOverlap)


class EndBulgePreviewConduit(Rhino.Display.DisplayConduit):
    def __init__(self):
        super(EndBulgePreviewConduit, self).__init__()
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.crv = None
        self.show_graph = Opts.values['bShowGraph']
        self.show_polygon = Opts.values['bShowPolygon']
        self.show_geom = Opts.values['bShowGeom']
        self.graph_scale = Opts.values['iGraphScale']
        self.graph_density = Opts.values['iGraphDensity']

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if self.crv is None: return
        bbox = self.crv.GetBoundingBox(accurate=False) 
        bbox.Inflate(sc.doc.ModelAbsoluteTolerance * 100)
        calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PostDrawObjects(self, drawEventArgs):
        if self.crv is None: return

        if self.show_geom:
            displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
            crv_thk = displayMode.DisplayAttributes.CurveThickness + 1
            drawEventArgs.Display.DrawCurve(curve=self.crv, color=self.color, thickness=crv_thk)
            
        if self.show_polygon:
            cp_locations = [pt.Location for pt in self.crv.Points]
            drawEventArgs.Display.DrawPatternedPolyline(
                points=cp_locations, color=self.color, pattern=0x00001111, thickness=1, close=False)
            drawEventArgs.Display.DrawPoints(
                points=cp_locations, style=Rhino.Display.PointStyle.Simple, radius=3, color=self.color)
                
        if self.show_graph:
            drawEventArgs.Display.DrawCurvatureGraph(
                curve=self.crv, color=self.color, hairScale=self.graph_scale, hairDensity=self.graph_density, sampleDensity=2)


class EtoDialog(ef.Dialog):
    def __init__(self, objref_In):
        self.is_surface = getattr(self, 'is_surface', False)
        self.Title = "EndBulge by SPB"
        self.objref_In = objref_In
        self.dialog_ok = False

        rgC_In = objref_In.Curve()
        self.nc_In = rgC_In.ToNurbsCurve()

        self._exact_scale_picked = Opts.values['fScale_Picked']
        self._exact_scale_opp = Opts.values['fScale_Opp']
        self._auto_updating = False
        self._auto_updating_slider = False
        self.active_stepper_key = None

        self.create_controls()
        self.setup_layout()
        self.OnLinkedModeChanged(None, None)

        # Initialize the Debounce Timer (200ms)
        self.debounce_timer = ef.UITimer()
        self.debounce_timer.Interval = 0.2
        self.debounce_timer.Elapsed += self.OnDebounceTimerElapsed

        self.LoadComplete += self.OnFormLoadComplete
        self.Closed += self.OnFormClosed

    def UpdatePreview(self):
        if not hasattr(self, 'conduit') or self.conduit is None: return

        fScale_Picked = self.ParseToFloat(self.textBoxes['fScale_Picked'].Text)
        fScale_Opp = self.ParseToFloat(self.textBoxes['fScale_Opp'].Text)
        
        if (fScale_Picked is None or fScale_Picked <= Rhino.RhinoMath.ZeroTolerance or
            fScale_Opp is None or fScale_Opp <= Rhino.RhinoMath.ZeroTolerance):
            self.conduit.crv = None
            sc.doc.Views.Redraw()
            return

        fSlideG2_Picked = self.ParseToFloat(self.textBoxes['fSlideG2_Picked'].Text)
        fSlideG2_Picked = fSlideG2_Picked if fSlideG2_Picked is not None else 0.0
        
        fSlideG3_Picked = self.ParseToFloat(self.textBoxes['fSlideG3_Picked'].Text)
        fSlideG3_Picked = fSlideG3_Picked if fSlideG3_Picked is not None else 0.0

        fSlideG2_Opp = self.ParseToFloat(self.textBoxes['fSlideG2_Opp'].Text)
        fSlideG2_Opp = fSlideG2_Opp if fSlideG2_Opp is not None else 0.0
        
        fSlideG3_Opp = self.ParseToFloat(self.textBoxes['fSlideG3_Opp'].Text)
        fSlideG3_Opp = fSlideG3_Opp if fSlideG3_Opp is not None else 0.0

        idxCont_Picked = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        idxCont_Opp = self.radioButtonLists['idxCont_Opp'].SelectedIndex
        bDebug = self.checkBoxes['bDebug'].Checked

        bSuccess, t_AtPicked = self.nc_In.ClosestPoint(self.objref_In.SelectionPoint())
        if not bSuccess: return

        if t_AtPicked > self.nc_In.Domain.Mid:
            iPickedEnd = 1
            fScale_T1, fSlideG2_T1, fSlideG3_T1 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked
            fScale_T0, fSlideG2_T0, fSlideG3_T0 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp
            iG_T1, iG_T0 = idxCont_Picked - 1, idxCont_Opp - 1
        else:
            iPickedEnd = 0
            fScale_T0, fSlideG2_T0, fSlideG3_T0 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked
            fScale_T1, fSlideG2_T1, fSlideG3_T1 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp
            iG_T0, iG_T1 = idxCont_Picked - 1, idxCont_Opp - 1

        nc_Res, sReport, info = createCurve(
            nc_In=self.nc_In,
            fScale_T0=fScale_T0, fSlideG2_T0=fSlideG2_T0, fSlideG3_T0=fSlideG3_T0,
            fScale_T1=fScale_T1, fSlideG2_T1=fSlideG2_T1, fSlideG3_T1=fSlideG3_T1,
            iG_T0=iG_T0, iG_T1=iG_T1, iPickedEnd=iPickedEnd, bDebug=bDebug
        )

        if info is not None:
            actual_T0, actual_T1, bOverlap = info
            actual_Picked = actual_T1 if (t_AtPicked > self.nc_In.Domain.Mid) else actual_T0
            actual_Opp = actual_T0 if (t_AtPicked > self.nc_In.Domain.Mid) else actual_T1
            
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
                    
                if changed:
                    self.UpdateControlStates()
                self._auto_updating = False

        if nc_Res is None:
            self.conduit.crv = self.nc_In.Duplicate()
        else:
            self.conduit.crv = nc_Res
            
        sc.doc.Views.Redraw()

        # Restart the debounce timer so it doesn't fire during active sliding
        self.debounce_timer.Stop()
        self.debounce_timer.Start()

    def OnFormLoadComplete(self, sender, e):
        loc_key = 'spb_EndBulge_WindowLoc'
        if loc_key in sc.sticky:
            saved_x, saved_y = sc.sticky[loc_key][0], sc.sticky[loc_key][1]
            if saved_x > 0 and saved_y > 0:
                self.Location = ed.Point(saved_x, saved_y)

    def create_controls(self):
        self.labels = {}
        self.checkBoxes = {}
        self.radioButtonLists = {}
        self.numericSteppers = {}
        self.textBoxes = {}
        self.dropDowns = {}
        self.sliders = {}
        self.btnUp = {}
        self.btnDown = {}

        bSuccess, t_AtPicked = self.nc_In.ClosestPoint(self.objref_In.SelectionPoint())
        bPickedIsT1 = t_AtPicked > self.nc_In.Domain.Mid
        can_G3_Picked = canMaintainG3(self.nc_In, bPickedIsT1)
        can_G3_Opp = canMaintainG3(self.nc_In, not bPickedIsT1)

        list_Picked = Opts.listValues['idxCont_Picked'] if can_G3_Picked else ('None', 'G0', 'G1', 'G2')
        list_Opp = Opts.listValues['idxCont_Opp'] if can_G3_Opp else ('None', 'G0', 'G1', 'G2')

        val_picked = int(Opts.values['idxCont_Picked'])
        if not can_G3_Picked and val_picked == 4: val_picked = 3
        val_opp = int(Opts.values['idxCont_Opp'])
        if not can_G3_Opp and val_opp == 4: val_opp = 3

        term_low = "edge" if self.is_surface else "end"
        small_font = ed.Font(ed.SystemFont.Default, 4)

        key = 'idxCont_Picked'
        self.labels[key] = ef.Label(Text = "Picked {}:".format(term_low))
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(4, 4)
        self.radioButtonLists[key].DataStore = list_Picked
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[val_picked]
        self.radioButtonLists[key].SelectedIndexChanged += self.OnContinuityChanged

        key = 'idxCont_Opp'
        self.labels[key] = ef.Label(Text = "Opp. {}:".format(term_low))
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(4, 4)
        self.radioButtonLists[key].DataStore = list_Opp
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[val_opp]
        self.radioButtonLists[key].SelectedIndexChanged += self.OnContinuityChanged

        key = 'bLinkedEnds'
        self.labels[key] = ef.Label(Text = "Adjust {}s:".format(term_low))
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Orientation = ef.Orientation.Horizontal
        self.radioButtonLists[key].Spacing = ed.Size(16, 4)
        self.radioButtonLists[key].DataStore = (Opts.offValues[key], Opts.onValues[key])
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[int(Opts.values[key])]
        self.radioButtonLists[key].SelectedIndexChanged += self.OnLinkedModeChanged

        key = 'fIncrement'
        self.labels[key] = ef.Label(Text = "Incr.:")
        self.textBoxes[key] = ef.TextBox()
        self.textBoxes[key].Text = str(Opts.values[key])
        self.textBoxes[key].TextChanged += self.OnIncrementTextChanged

        key = 'iSliderSteps'
        self.labels[key] = ef.Label(Text = "Slider steps:")
        self.dropDowns[key] = ef.DropDown()
        self.dropDowns[key].DataStore = Opts.listValues[key]
        self.dropDowns[key].SelectedIndex = Opts.values[key]
        self.dropDowns[key].SelectedIndexChanged += self.OnSliderStepsChanged

        self.hold_timer = ef.UITimer()
        self.hold_timer.Interval = 0.15 
        self.hold_timer.Elapsed += self.OnHoldTimerElapsed
        self.hold_direction = 0

        self.slider_prev_vals = {}

        # --- JOG SLIDER GENERATOR ---
        def create_jog_slider(s_key):
            slider = ef.Slider()
            slider.Width = 132  
            slider.SnapToTick = True
            slider.TickFrequency = 1
            self.slider_prev_vals[s_key] = 0
            slider.ValueChanged += lambda s, e, k=s_key: self.OnJogSliderChanged(s, k)
            slider.MouseUp += lambda s, e, k=s_key: self.ZeroSlider(s, k)
            self.sliders[s_key] = slider

        # --- HOMEMADE STEPPER GENERATOR ---
        def create_homemade_stepper(s_key, label_text, is_scale=False):
            self.labels[s_key] = ef.Label(Text = label_text)
            self.textBoxes[s_key] = ef.TextBox()
            self.textBoxes[s_key].Text = "{:.4f}".format(float(Opts.values[s_key]))
            
            if is_scale:
                self.textBoxes[s_key].TextChanged += self.OnScaleTextChanged
            else:
                self.textBoxes[s_key].TextChanged += self.OnSlideTextChanged
                
            self.textBoxes[s_key].MouseWheel += lambda s, e, k=s_key: self.OnStepperMouseWheel(s, e, k)
            self.textBoxes[s_key].KeyDown += lambda s, e, k=s_key: self.OnStepperKeyDown(s, e, k)

            create_jog_slider(s_key)

            self.btnUp[s_key] = ef.Button(Text=unichr(9650), Width=16, Height=12)
            self.btnDown[s_key] = ef.Button(Text=unichr(9660), Width=16, Height=12)
            self.btnUp[s_key].Font = small_font
            self.btnDown[s_key].Font = small_font
            self.btnUp[s_key].MinimumSize = ed.Size(16, 12)
            self.btnDown[s_key].MinimumSize = ed.Size(16, 12)
            
            self.btnUp[s_key].MouseDown += lambda s, e, k=s_key: self.StartHoldTimer(1, k)
            self.btnUp[s_key].MouseUp += self.StopHoldTimer
            self.btnUp[s_key].MouseLeave += self.StopHoldTimer
            
            self.btnDown[s_key].MouseDown += lambda s, e, k=s_key: self.StartHoldTimer(-1, k)
            self.btnDown[s_key].MouseUp += self.StopHoldTimer
            self.btnDown[s_key].MouseLeave += self.StopHoldTimer

        create_homemade_stepper('fScale_Picked', "Scale:", True)
        create_homemade_stepper('fSlideG2_Picked', "G2 slide:", False)
        create_homemade_stepper('fSlideG3_Picked', "G3 slide:", False)
        
        create_homemade_stepper('fScale_Opp', "Scale:", True)
        create_homemade_stepper('fSlideG2_Opp', "G2 slide:", False)
        create_homemade_stepper('fSlideG3_Opp', "G3 slide:", False)

        self.UpdateSliderRanges()

        key = 'bShowGeom'
        self.checkBoxes[key] = ef.CheckBox()
        self.checkBoxes[key].Text = "Surface" if self.is_surface else "Curve"
        self.checkBoxes[key].Checked = Opts.values[key]
        self.checkBoxes[key].CheckedChanged += self.OnDisplayCheckedChanged

        key = 'bShowPolygon'
        self.checkBoxes[key] = ef.CheckBox()
        self.checkBoxes[key].Text = "Control polygon"
        self.checkBoxes[key].Checked = Opts.values[key]
        self.checkBoxes[key].CheckedChanged += self.OnDisplayCheckedChanged

        key = 'bShowGraph'
        self.checkBoxes[key] = ef.CheckBox()
        self.checkBoxes[key].Text = "CGraph"
        self.checkBoxes[key].Checked = Opts.values[key]
        self.checkBoxes[key].CheckedChanged += self.OnDisplayCheckedChanged

        key = 'iGraphScale'
        self.labels[key] = ef.Label(Text = "Scale:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 0 
        self.numericSteppers[key].MinValue = 1.0
        self.numericSteppers[key].MaxValue = 10000.0
        self.numericSteppers[key].Increment = 1.0
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnGraphScaleStepperChanged

        key = 'iGraphDensity'
        self.labels[key] = ef.Label(Text = "Density:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 0
        self.numericSteppers[key].MinValue = 0.0
        self.numericSteppers[key].MaxValue = 100.0
        self.numericSteppers[key].Increment = 1.0
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnGraphDensityStepperChanged

        key = 'bDeleteInput'
        self.checkBoxes[key] = ef.CheckBox()
        self.checkBoxes[key].Checked = Opts.values[key]
        self.checkBoxes[key].Text = "Delete input"

        key = 'bEcho'
        self.checkBoxes[key] = ef.CheckBox()
        self.checkBoxes[key].Checked = Opts.values[key]
        self.checkBoxes[key].Text = Opts.names[key]

        key = 'bDebug'
        self.checkBoxes[key] = ef.CheckBox()
        self.checkBoxes[key].Checked = Opts.values[key]
        self.checkBoxes[key].Text = Opts.names[key]
        self.checkBoxes[key].CheckedChanged += lambda s, e: self.UpdatePreview()

        # --- COMPACT UI WIDTHS ---
        lbl_w = 50
        for k in ('fIncrement', 'fScale_Picked', 'fSlideG2_Picked', 'fSlideG3_Picked', 'fScale_Opp', 'fSlideG2_Opp', 'fSlideG3_Opp'):
            self.labels[k].Width = lbl_w

        self.textBoxes['fIncrement'].Width = 80 
        self.labels['iSliderSteps'].Width = 64
        self.dropDowns['iSliderSteps'].Width = 60

        # All 6 configuration textboxes share the 64px width constraint
        for k in ('fScale_Picked', 'fSlideG2_Picked', 'fSlideG3_Picked', 'fScale_Opp', 'fSlideG2_Opp', 'fSlideG3_Opp'):
            self.textBoxes[k].Width = 64

        self.numericSteppers['iGraphScale'].Width = 45
        self.numericSteppers['iGraphDensity'].Width = 45

    def setup_layout(self):
        term_cap = "Edge" if self.is_surface else "End"
        
        def wrap(ctrl):
            st = ef.StackLayout()
            st.Orientation = ef.Orientation.Horizontal
            st.Items.Add(ctrl)
            return st
            
        def gap():
            return ef.Label(Width=12)
            
        def build_combo(key):
            stepper = ef.StackLayout(Spacing=0)
            stepper.Items.Add(self.btnUp[key])
            stepper.Items.Add(self.btnDown[key])
            combo = ef.StackLayout(Orientation=ef.Orientation.Horizontal, Spacing=0)
            combo.Items.Add(self.textBoxes[key])
            combo.Items.Add(stepper)
            return combo
        
        root_layout = ef.StackLayout()
        root_layout.Padding = ed.Padding(10)
        root_layout.Spacing = 8
        root_layout.HorizontalContentAlignment = ef.HorizontalAlignment.Left

        root_layout.Items.Add(ef.Label(Text="Continuity Constraints", Font=ed.Font(ed.SystemFont.Bold, 10)))
        cont_grid = ef.DynamicLayout()
        cont_grid.Spacing = ed.Size(4, 4)
        cont_grid.AddRow(self.labels['idxCont_Picked'], self.radioButtonLists['idxCont_Picked'])
        cont_grid.AddRow(self.labels['idxCont_Opp'], self.radioButtonLists['idxCont_Opp'])
        root_layout.Items.Add(cont_grid)
        root_layout.Items.Add(ef.Label(Height=4))

        root_layout.Items.Add(ef.Label(Text="Configurations", Font=ed.Font(ed.SystemFont.Bold, 10)))
        adj_row = ef.StackLayout()
        adj_row.Orientation = ef.Orientation.Horizontal
        adj_row.Spacing = 8
        adj_row.VerticalContentAlignment = ef.VerticalAlignment.Center
        adj_row.Items.Add(self.labels['bLinkedEnds'])
        adj_row.Items.Add(self.radioButtonLists['bLinkedEnds'])
        root_layout.Items.Add(adj_row)
        root_layout.Items.Add(ef.Label(Height=2))

        incr_grid = ef.DynamicLayout()
        incr_grid.Spacing = ed.Size(4, 4)
        
        slider_steps_stack = ef.StackLayout()
        slider_steps_stack.Orientation = ef.Orientation.Horizontal
        slider_steps_stack.Spacing = 8
        slider_steps_stack.VerticalContentAlignment = ef.VerticalAlignment.Center
        slider_steps_stack.Items.Add(self.labels['iSliderSteps'])
        slider_steps_stack.Items.Add(wrap(self.dropDowns['iSliderSteps']))

        incr_grid.AddRow(self.labels['fIncrement'], wrap(self.textBoxes['fIncrement']), gap(), slider_steps_stack, None)
        root_layout.Items.Add(incr_grid)
        root_layout.Items.Add(ef.Label(Height=6))

        root_layout.Items.Add(ef.Label(Text="Picked {}".format(term_cap), Font=ed.Font(ed.SystemFont.Bold)))
        picked_grid = ef.DynamicLayout()
        picked_grid.Spacing = ed.Size(4, 4)
        picked_grid.AddRow(self.labels['fScale_Picked'], build_combo('fScale_Picked'), gap(), self.sliders['fScale_Picked'], None)
        picked_grid.AddRow(self.labels['fSlideG2_Picked'], build_combo('fSlideG2_Picked'), gap(), self.sliders['fSlideG2_Picked'], None)
        picked_grid.AddRow(self.labels['fSlideG3_Picked'], build_combo('fSlideG3_Picked'), gap(), self.sliders['fSlideG3_Picked'], None)
        
        root_layout.Items.Add(picked_grid)
        root_layout.Items.Add(ef.Label(Height=6))

        root_layout.Items.Add(ef.Label(Text="Opposite {}".format(term_cap), Font=ed.Font(ed.SystemFont.Bold)))
        opp_grid = ef.DynamicLayout()
        opp_grid.Spacing = ed.Size(4, 4)
        opp_grid.AddRow(self.labels['fScale_Opp'], build_combo('fScale_Opp'), gap(), self.sliders['fScale_Opp'], None)
        opp_grid.AddRow(self.labels['fSlideG2_Opp'], build_combo('fSlideG2_Opp'), gap(), self.sliders['fSlideG2_Opp'], None)
        opp_grid.AddRow(self.labels['fSlideG3_Opp'], build_combo('fSlideG3_Opp'), gap(), self.sliders['fSlideG3_Opp'], None)

        root_layout.Items.Add(opp_grid)
        root_layout.Items.Add(ef.Label(Height=4))

        root_layout.Items.Add(ef.Label(Text="Display", Font=ed.Font(ed.SystemFont.Bold, 10)))
        display_group = ef.StackLayout()
        display_group.Spacing = 4
        
        disp_chk_stack = ef.StackLayout()
        disp_chk_stack.Orientation = ef.Orientation.Horizontal
        disp_chk_stack.Spacing = 8
        disp_chk_stack.Items.Add(self.checkBoxes['bShowGeom'])
        disp_chk_stack.Items.Add(self.checkBoxes['bShowPolygon'])
        
        analysis_grid = ef.DynamicLayout()
        analysis_grid.Spacing = ed.Size(4, 4)
        analysis_grid.AddRow(self.checkBoxes['bShowGraph'], ef.Label(Width=4), self.labels['iGraphScale'], wrap(self.numericSteppers['iGraphScale']), ef.Label(Width=10), self.labels['iGraphDensity'], wrap(self.numericSteppers['iGraphDensity']), None)
        
        display_group.Items.Add(disp_chk_stack)
        display_group.Items.Add(analysis_grid)
        
        root_layout.Items.Add(display_group)
        root_layout.Items.Add(ef.Label(Height=4))

        chk_stack = ef.StackLayout()
        chk_stack.Orientation = ef.Orientation.Horizontal
        chk_stack.Spacing = 20
        chk_stack.Items.Add(self.checkBoxes['bDeleteInput'])
        chk_stack.Items.Add(self.checkBoxes['bEcho'])
        chk_stack.Items.Add(self.checkBoxes['bDebug'])
        
        root_layout.Items.Add(chk_stack)
        root_layout.Items.Add(ef.Label(Height=8))

        self.ok_button = ef.Button(Text = 'OK')
        self.ok_button.Click += self.OnOKButtonClick
        save_button = ef.Button(Text = 'Save Settings')
        save_button.Click += self.OnSaveSettingsButtonClick
        self.abort_button = ef.Button(Text = 'Cancel')
        self.abort_button.Click += self.OnCancelButtonClick

        self.DefaultButton = self.ok_button
        self.AbortButton = self.abort_button

        btn_stack = ef.StackLayout()
        btn_stack.Orientation = ef.Orientation.Horizontal
        btn_stack.Spacing = 8
        btn_stack.Items.Add(self.ok_button)
        btn_stack.Items.Add(save_button)
        btn_stack.Items.Add(self.abort_button)

        root_layout.Items.Add(btn_stack)

        self.Content = root_layout
        self.AutoSize = True
        self.Resizable = False

    def UpdateSliderRanges(self):
        steps = int(self.dropDowns['iSliderSteps'].SelectedValue)
        for s in self.sliders.values():
            s.MinValue = -steps
            s.MaxValue = steps

    def OnSliderStepsChanged(self, sender, e):
        self.UpdateSliderRanges()

    def OnJogSliderChanged(self, sender, target_key):
        if self._auto_updating_slider: return
        
        delta = sender.Value - self.slider_prev_vals[target_key]
        self.slider_prev_vals[target_key] = sender.Value
        if delta == 0: return

        incr_val = self.ParseToFloat(self.textBoxes['fIncrement'].Text)
        if incr_val is None: return

        change = delta * incr_val

        self._auto_updating = True
        
        if 'Scale' in target_key:
            current_val = self._exact_scale_picked if 'Picked' in target_key else self._exact_scale_opp
            new_val = max(incr_val, current_val + change)
            if 'Picked' in target_key: self._exact_scale_picked = new_val
            else: self._exact_scale_opp = new_val
        else:
            current_val = self.ParseToFloat(self.textBoxes[target_key].Text)
            if current_val is None: current_val = 0.0
            new_val = current_val + change
            
        self.textBoxes[target_key].Text = "{:.4f}".format(new_val)

        self._auto_updating = False
        self.SyncLinkedControls()
        self.UpdatePreview()

    def ZeroSlider(self, slider, target_key):
        self._auto_updating_slider = True
        slider.Value = 0
        self.slider_prev_vals[target_key] = 0
        self._auto_updating_slider = False

    def SyncLinkedControls(self, sender=None, e=None):
        if bool(self.radioButtonLists['bLinkedEnds'].SelectedIndex):
            if self.textBoxes['fScale_Opp'].Text != self.textBoxes['fScale_Picked'].Text:
                self.textBoxes['fScale_Opp'].Text = self.textBoxes['fScale_Picked'].Text
            if self.textBoxes['fSlideG2_Opp'].Text != self.textBoxes['fSlideG2_Picked'].Text:
                self.textBoxes['fSlideG2_Opp'].Text = self.textBoxes['fSlideG2_Picked'].Text
            if self.textBoxes['fSlideG3_Opp'].Text != self.textBoxes['fSlideG3_Picked'].Text:
                self.textBoxes['fSlideG3_Opp'].Text = self.textBoxes['fSlideG3_Picked'].Text

    def UpdateControlStates(self):
        is_linked = bool(self.radioButtonLists['bLinkedEnds'].SelectedIndex)
        idx_picked = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        idx_opp = self.radioButtonLists['idxCont_Opp'].SelectedIndex

        N = self.nc_In.Points.Count
        if idx_picked + idx_opp > N:
            if idx_picked > idx_opp:
                alloc_P = min(idx_picked, N)
                alloc_O = N - alloc_P
            elif idx_opp > idx_picked:
                alloc_O = min(idx_opp, N)
                alloc_P = N - alloc_O
            else: 
                alloc_P = min(idx_picked, N)
                alloc_O = N - alloc_P
        else:
            alloc_P = idx_picked
            alloc_O = idx_opp

        free = N - alloc_P - alloc_O
        if free > 0:
            half = free // 2
            extra = free % 2
            scale_limit_P = alloc_P + half + extra
            scale_limit_O = alloc_O + half
        else:
            scale_limit_P = alloc_P
            scale_limit_O = alloc_O

        self.textBoxes['fSlideG2_Picked'].Enabled = (scale_limit_P >= 3)
        self.btnUp['fSlideG2_Picked'].Enabled = (scale_limit_P >= 3)
        self.btnDown['fSlideG2_Picked'].Enabled = (scale_limit_P >= 3)
        self.sliders['fSlideG2_Picked'].Enabled = (scale_limit_P >= 3)
        
        self.textBoxes['fSlideG3_Picked'].Enabled = (scale_limit_P >= 4)
        self.btnUp['fSlideG3_Picked'].Enabled = (scale_limit_P >= 4)
        self.btnDown['fSlideG3_Picked'].Enabled = (scale_limit_P >= 4)
        self.sliders['fSlideG3_Picked'].Enabled = (scale_limit_P >= 4)

        self.textBoxes['fScale_Opp'].Enabled = not is_linked
        self.btnUp['fScale_Opp'].Enabled = not is_linked
        self.btnDown['fScale_Opp'].Enabled = not is_linked
        self.sliders['fScale_Opp'].Enabled = not is_linked
        
        if is_linked:
            self.textBoxes['fSlideG2_Opp'].Enabled = False
            self.btnUp['fSlideG2_Opp'].Enabled = False
            self.btnDown['fSlideG2_Opp'].Enabled = False
            self.sliders['fSlideG2_Opp'].Enabled = False
            
            self.textBoxes['fSlideG3_Opp'].Enabled = False
            self.btnUp['fSlideG3_Opp'].Enabled = False
            self.btnDown['fSlideG3_Opp'].Enabled = False
            self.sliders['fSlideG3_Opp'].Enabled = False
        else:
            self.textBoxes['fSlideG2_Opp'].Enabled = (scale_limit_O >= 3)
            self.btnUp['fSlideG2_Opp'].Enabled = (scale_limit_O >= 3)
            self.btnDown['fSlideG2_Opp'].Enabled = (scale_limit_O >= 3)
            self.sliders['fSlideG2_Opp'].Enabled = (scale_limit_O >= 3)
            
            self.textBoxes['fSlideG3_Opp'].Enabled = (scale_limit_O >= 4)
            self.btnUp['fSlideG3_Opp'].Enabled = (scale_limit_O >= 4)
            self.btnDown['fSlideG3_Opp'].Enabled = (scale_limit_O >= 4)
            self.sliders['fSlideG3_Opp'].Enabled = (scale_limit_O >= 4)

    def OnContinuityChanged(self, sender, e):
        if getattr(self, '_auto_updating', False): return
        self._auto_updating = True

        N = self.nc_In.Points.Count
        idx_picked = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        idx_opp = self.radioButtonLists['idxCont_Opp'].SelectedIndex
        downgraded = False

        if sender == self.radioButtonLists['idxCont_Picked']:
            if idx_picked + idx_opp > N:
                idx_opp = max(0, N - idx_picked)
                self.radioButtonLists['idxCont_Opp'].SelectedIndex = idx_opp
                downgraded = True
        elif sender == self.radioButtonLists['idxCont_Opp']:
            if idx_picked + idx_opp > N:
                idx_picked = max(0, N - idx_opp)
                self.radioButtonLists['idxCont_Picked'].SelectedIndex = idx_picked
                downgraded = True

        if downgraded and self.radioButtonLists['bLinkedEnds'].SelectedIndex == 1:
            self.radioButtonLists['bLinkedEnds'].SelectedIndex = 0
            
        self._auto_updating = False
        self.UpdateControlStates()
        self.UpdatePreview()

    def OnLinkedModeChanged(self, sender, e):
        self.UpdateControlStates()
        if bool(self.radioButtonLists['bLinkedEnds'].SelectedIndex):
            self.SyncLinkedControls()
        self.UpdatePreview()

    def OnScaleTextChanged(self, sender, e):
        val = self.ParseToFloat(sender.Text)
        
        if not self._auto_updating and val is not None:
            if sender == self.textBoxes['fScale_Picked']:
                self._exact_scale_picked = val
            elif sender == self.textBoxes['fScale_Opp']:
                self._exact_scale_opp = val
                
        if val is not None and val > Rhino.RhinoMath.ZeroTolerance:
            sender.BackgroundColor = ed.Colors.White
        else:
            sender.BackgroundColor = ed.Colors.LightPink
            
        if not self._auto_updating:
            self.SyncLinkedControls()
            self.UpdatePreview()

    def OnSlideTextChanged(self, sender, e):
        val = self.ParseToFloat(sender.Text)
        if val is not None:
            sender.BackgroundColor = ed.Colors.White
        else:
            sender.BackgroundColor = ed.Colors.LightPink
            
        if not self._auto_updating:
            self.SyncLinkedControls()
            self.UpdatePreview()

    def ParseToFloat(self, text):
        text = text.strip()
        try:
            if '/' in text:
                num, den = text.split('/')
                return float(num) / float(den)
            return float(text)
        except (ValueError, ZeroDivisionError):
            return None

    def OnStepperMouseWheel(self, sender, e, key):
        import time
        now = time.time()
        # Time Gate: Ignore event if less than 0.15 seconds have passed
        if now - getattr(self, '_last_input_time', 0) < 0.15:
            e.Handled = True
            return
        self._last_input_time = now

        direction = 1 if e.Delta.Height > 0 else -1
        self.active_stepper_key = key
        self.AdjustStepper(direction, key)
        e.Handled = True

    def OnStepperKeyDown(self, sender, e, key):
        if e.Key == ef.Keys.Up or e.Key == ef.Keys.Down:
            import time
            now = time.time()
            # Time Gate: Ignore event if less than 0.15 seconds have passed
            if now - getattr(self, '_last_input_time', 0) < 0.15:
                e.Handled = True
                return
            self._last_input_time = now

            direction = 1 if e.Key == ef.Keys.Up else -1
            self.active_stepper_key = key
            self.AdjustStepper(direction, key)
            e.Handled = True

    def OnDebounceTimerElapsed(self, sender, e):
        self.debounce_timer.Stop()

    def StartHoldTimer(self, direction, key):
        self.active_stepper_key = key
        self.hold_direction = direction
        self.AdjustStepper(direction, key)
        self.hold_timer.Start()

    def StopHoldTimer(self, sender, e):
        self.hold_timer.Stop()

    def OnHoldTimerElapsed(self, sender, e):
        self.AdjustStepper(self.hold_direction, self.active_stepper_key)

    def AdjustStepper(self, direction, key):
        incr_val = self.ParseToFloat(self.textBoxes['fIncrement'].Text)
        if incr_val is None: return

        if 'Scale' in key:
            current_val = self._exact_scale_picked if 'Picked' in key else self._exact_scale_opp
        else:
            current_val = self.ParseToFloat(self.textBoxes[key].Text)

        if current_val is None: return

        new_val = current_val + (incr_val * direction)
        
        if 'Scale' in key:
            new_val = max(incr_val, new_val)
            if 'Picked' in key:
                self._exact_scale_picked = new_val
            else:
                self._exact_scale_opp = new_val
                
        self._auto_updating = True
        self.textBoxes[key].Text = "{:.4f}".format(new_val)
        self._auto_updating = False

        self.SyncLinkedControls()
        self.UpdatePreview()

    def OnIncrementTextChanged(self, sender, e):
        text = sender.Text.strip()
        val = None
        try:
            if '/' in text:
                num, den = text.split('/')
                val = float(num) / float(den)
            else:
                val = float(text)
        except (ValueError, ZeroDivisionError): pass
        
        if val is not None and val > Rhino.RhinoMath.ZeroTolerance:
            sender.BackgroundColor = ed.Colors.White 
        else:
            sender.BackgroundColor = ed.Colors.LightPink

    def OnDisplayCheckedChanged(self, sender, e):
        if hasattr(self, 'conduit') and self.conduit is not None:
            self.conduit.show_geom = self.checkBoxes['bShowGeom'].Checked
            self.conduit.show_polygon = self.checkBoxes['bShowPolygon'].Checked
            self.conduit.show_graph = self.checkBoxes['bShowGraph'].Checked
            sc.doc.Views.Redraw()

    def OnGraphScaleStepperChanged(self, sender, e):
        if hasattr(self, 'conduit') and self.conduit is not None:
            self.conduit.graph_scale = sender.Value
            sc.doc.Views.Redraw()

    def OnGraphDensityStepperChanged(self, sender, e):
        if hasattr(self, 'conduit') and self.conduit is not None:
            self.conduit.graph_density = int(sender.Value)
            sc.doc.Views.Redraw()

    def SaveSettings(self):
        sc.sticky[Opts.stickyKeys['bLinkedEnds']] = Opts.values['bLinkedEnds'] = bool(self.radioButtonLists['bLinkedEnds'].SelectedIndex)
        sc.sticky[Opts.stickyKeys['iSliderSteps']] = Opts.values['iSliderSteps'] = self.dropDowns['iSliderSteps'].SelectedIndex
        
        parsed_incr = self.ParseToFloat(self.textBoxes['fIncrement'].Text)
        if parsed_incr is not None:
            sc.sticky[Opts.stickyKeys['fIncrement']] = Opts.values['fIncrement'] = parsed_incr

        for end in ('_Picked', '_Opp'):
            parsed_scale = self.ParseToFloat(self.textBoxes['fScale' + end].Text)
            if parsed_scale is not None: sc.sticky[Opts.stickyKeys['fScale' + end]] = Opts.values['fScale' + end] = parsed_scale
            
            p_g2 = self.ParseToFloat(self.textBoxes['fSlideG2' + end].Text)
            if p_g2 is not None: sc.sticky[Opts.stickyKeys['fSlideG2' + end]] = Opts.values['fSlideG2' + end] = p_g2
            
            p_g3 = self.ParseToFloat(self.textBoxes['fSlideG3' + end].Text)
            if p_g3 is not None: sc.sticky[Opts.stickyKeys['fSlideG3' + end]] = Opts.values['fSlideG3' + end] = p_g3

        sc.sticky[Opts.stickyKeys['idxCont_Picked']] = Opts.values['idxCont_Picked'] = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        sc.sticky[Opts.stickyKeys['idxCont_Opp']] = Opts.values['idxCont_Opp'] = self.radioButtonLists['idxCont_Opp'].SelectedIndex
        sc.sticky[Opts.stickyKeys['bShowGeom']] = Opts.values['bShowGeom'] = self.checkBoxes['bShowGeom'].Checked
        sc.sticky[Opts.stickyKeys['bShowPolygon']] = Opts.values['bShowPolygon'] = self.checkBoxes['bShowPolygon'].Checked
        sc.sticky[Opts.stickyKeys['bShowGraph']] = Opts.values['bShowGraph'] = self.checkBoxes['bShowGraph'].Checked
        sc.sticky[Opts.stickyKeys['iGraphScale']] = Opts.values['iGraphScale'] = int(self.numericSteppers['iGraphScale'].Value)
        sc.sticky[Opts.stickyKeys['iGraphDensity']] = Opts.values['iGraphDensity'] = int(self.numericSteppers['iGraphDensity'].Value)
        sc.sticky[Opts.stickyKeys['bDeleteInput']] = Opts.values['bDeleteInput'] = self.checkBoxes['bDeleteInput'].Checked
        sc.sticky[Opts.stickyKeys['bEcho']] = Opts.values['bEcho'] = self.checkBoxes['bEcho'].Checked
        sc.sticky[Opts.stickyKeys['bDebug']] = Opts.values['bDebug'] = self.checkBoxes['bDebug'].Checked

        # SYNCHRONIZE CLI OPTIONS: Update the riOpts instances so the command line matches the GUI
        for k, opt in Opts.riOpts.items():
            if k in Opts.values:
                opt.CurrentValue = Opts.values[k]

    def OnOKButtonClick(self, sender, e):
        fScale_Picked = self.ParseToFloat(self.textBoxes['fScale_Picked'].Text)
        fScale_Opp = self.ParseToFloat(self.textBoxes['fScale_Opp'].Text)
        if (fScale_Picked is None or fScale_Picked <= Rhino.RhinoMath.ZeroTolerance or 
            fScale_Opp is None or fScale_Opp <= Rhino.RhinoMath.ZeroTolerance):
            print("Invalid inputs. No changes were applied.")
            self.dialog_ok = False
            self.Close()
            return
        self.SaveSettings()
        self.dialog_ok = True
        self.Close()

    def OnSaveSettingsButtonClick(self, sender, e):
        self.SaveSettings()
        print("Settings saved as default.")

    def OnCancelButtonClick(self, sender, e):
        self.Result = ef.DialogResult.Cancel
        self.Close()

    def OnFormClosed(self, sender, e):
        # We leave the document state alone here; main() will handle Unlock and Redraw
        sc.sticky['spb_EndBulge_WindowLoc'] = (self.Location.X, self.Location.Y)