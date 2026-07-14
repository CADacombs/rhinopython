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

This script was created by Google Gemini 3.1 Pro based on another script.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
260712-14: Created by extracting refactored code from another script.
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

    key = 'iGraphScale'; keys.append(key)
    values[key] = 100
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iGraphDensity'; keys.append(key)
    values[key] = 1
    riOpts[key] = ri.Custom.OptionInteger(values[key])
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
        return None, "Nothing to modify.", (max_mod_T0, max_mod_T1, bOverlap)

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

    if bDebug:
        for i in range(len(pts_Prime)):
            rgDot = rg.TextDot("{}".format(i), pts_Prime[i])
            rgDot.FontHeight = 11
            sc.doc.Objects.AddTextDot(rgDot)
        sc.doc.Views.Redraw()

    return nc_Out, None, (max_mod_T0, max_mod_T1, bOverlap)


class EndBulgePreviewConduit(Rhino.Display.DisplayConduit):
    def __init__(self):
        super(EndBulgePreviewConduit, self).__init__()
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.crv = None
        self.graph_scale = Opts.values['iGraphScale']
        self.graph_density = Opts.values['iGraphDensity']

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if self.crv is None: return
        bbox = self.crv.GetBoundingBox(accurate=False) 
        bbox.Inflate(sc.doc.ModelAbsoluteTolerance * 100)
        calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PostDrawObjects(self, drawEventArgs):
        if self.crv is None: return

        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

        drawEventArgs.Display.DrawCurve(curve=self.crv, color=self.color, thickness=crv_thk)
        self.cp_locations = [pt.Location for pt in self.crv.Points]

        drawEventArgs.Display.DrawPatternedPolyline(
            points=self.cp_locations, color=self.color, pattern=0x00001111, thickness=1, close=False)
        drawEventArgs.Display.DrawPoints(
            points=self.cp_locations, style=Rhino.Display.PointStyle.Simple, radius=3, color=self.color)
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

        self.create_controls()
        self.setup_layout()
        self.OnLinkedModeChanged(None, None)

    def UpdatePreview(self):
        if not hasattr(self, 'conduit') or self.conduit is None: return

        fScale_Picked = self.ParseToFloat(self.textBoxes['fScale_Picked'].Text)
        fSlideG2_Picked = self.numericSteppers['fSlideG2_Picked'].Value
        fSlideG3_Picked = self.numericSteppers['fSlideG3_Picked'].Value
        
        fScale_Opp = self.ParseToFloat(self.textBoxes['fScale_Opp'].Text)
        fSlideG2_Opp = self.numericSteppers['fSlideG2_Opp'].Value
        fSlideG3_Opp = self.numericSteppers['fSlideG3_Opp'].Value

        if (fScale_Picked is None or fScale_Picked <= Rhino.RhinoMath.ZeroTolerance or
            fScale_Opp is None or fScale_Opp <= Rhino.RhinoMath.ZeroTolerance):
            self.conduit.crv = None
            sc.doc.Views.Redraw()
            return

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

    def create_controls(self):
        self.labels = {}
        self.checkBoxes = {}
        self.radioButtonLists = {}
        self.numericSteppers = {}
        self.textBoxes = {}

        term_low = "edge" if self.is_surface else "end"
        small_font = ed.Font(ed.SystemFont.Default, 4)

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

        self.hold_timer = ef.UITimer()
        self.hold_timer.Interval = 0.15 
        self.hold_timer.Elapsed += self.OnHoldTimerElapsed
        self.hold_direction = 0
        self.active_stepper_end = 'Picked'

        key = 'fScale_Picked'
        self.labels[key] = ef.Label(Text = "Scale:")
        self.textBoxes[key] = ef.TextBox()
        self.textBoxes[key].Text = str(Opts.values[key])
        self.textBoxes[key].TextChanged += self.OnScaleTextChanged

        self.btnScaleUp_Picked = ef.Button(Text=unichr(9650), Width=16, Height=12)
        self.btnScaleDown_Picked = ef.Button(Text=unichr(9660), Width=16, Height=12)
        self.btnScaleUp_Picked.Font = small_font
        self.btnScaleDown_Picked.Font = small_font
        self.btnScaleUp_Picked.MinimumSize = ed.Size(16, 12)
        self.btnScaleDown_Picked.MinimumSize = ed.Size(16, 12)
        self.btnScaleUp_Picked.MouseDown += lambda s, e: self.StartHoldTimer(1, 'Picked')
        self.btnScaleUp_Picked.MouseUp += self.StopHoldTimer
        self.btnScaleUp_Picked.MouseLeave += self.StopHoldTimer
        self.btnScaleDown_Picked.MouseDown += lambda s, e: self.StartHoldTimer(-1, 'Picked')
        self.btnScaleDown_Picked.MouseUp += self.StopHoldTimer
        self.btnScaleDown_Picked.MouseLeave += self.StopHoldTimer

        key = 'fSlideG2_Picked'
        self.labels[key] = ef.Label(Text = "G2 slide:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 2
        self.numericSteppers[key].Increment = float(Opts.values['fIncrement'])
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnSlideValueChanged

        key = 'fSlideG3_Picked'
        self.labels[key] = ef.Label(Text = "G3 slide:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 2
        self.numericSteppers[key].Increment = float(Opts.values['fIncrement'])
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnSlideValueChanged

        key = 'fScale_Opp'
        self.labels[key] = ef.Label(Text = "Scale:")
        self.textBoxes[key] = ef.TextBox()
        self.textBoxes[key].Text = str(Opts.values[key])
        self.textBoxes[key].TextChanged += self.OnScaleTextChanged

        self.btnScaleUp_Opp = ef.Button(Text=unichr(9650), Width=16, Height=12)
        self.btnScaleDown_Opp = ef.Button(Text=unichr(9660), Width=16, Height=12)
        self.btnScaleUp_Opp.Font = small_font
        self.btnScaleDown_Opp.Font = small_font
        self.btnScaleUp_Opp.MinimumSize = ed.Size(16, 12)
        self.btnScaleDown_Opp.MinimumSize = ed.Size(16, 12)
        self.btnScaleUp_Opp.MouseDown += lambda s, e: self.StartHoldTimer(1, 'Opp')
        self.btnScaleUp_Opp.MouseUp += self.StopHoldTimer
        self.btnScaleUp_Opp.MouseLeave += self.StopHoldTimer
        self.btnScaleDown_Opp.MouseDown += lambda s, e: self.StartHoldTimer(-1, 'Opp')
        self.btnScaleDown_Opp.MouseUp += self.StopHoldTimer
        self.btnScaleDown_Opp.MouseLeave += self.StopHoldTimer

        key = 'fSlideG2_Opp'
        self.labels[key] = ef.Label(Text = "G2 slide:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 2
        self.numericSteppers[key].Increment = float(Opts.values['fIncrement'])
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnSlideValueChanged

        key = 'fSlideG3_Opp'
        self.labels[key] = ef.Label(Text = "G3 slide:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 2
        self.numericSteppers[key].Increment = float(Opts.values['fIncrement'])
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnSlideValueChanged

        key = 'iGraphScale'
        self.labels[key] = ef.Label(Text = "Graph scale:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 0 
        self.numericSteppers[key].MinValue = 1.0
        self.numericSteppers[key].MaxValue = 10000.0
        self.numericSteppers[key].Increment = 1.0
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnGraphScaleStepperChanged

        key = 'iGraphDensity'
        self.labels[key] = ef.Label(Text = "Graph density:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 0
        self.numericSteppers[key].MinValue = 0.0
        self.numericSteppers[key].MaxValue = 100.0
        self.numericSteppers[key].Increment = 1.0
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnGraphDensityStepperChanged

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

        key = 'idxCont_Picked'
        self.labels[key] = ef.Label(Text = "Maint. cont. of picked {}:".format(term_low))
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(4, 4)
        self.radioButtonLists[key].DataStore = list_Picked
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[val_picked]
        self.radioButtonLists[key].SelectedIndexChanged += self.OnContinuityChanged

        key = 'idxCont_Opp'
        self.labels[key] = ef.Label(Text = "Maint. cont. of opp. {}:".format(term_low))
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(4, 4)
        self.radioButtonLists[key].DataStore = list_Opp
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[val_opp]
        self.radioButtonLists[key].SelectedIndexChanged += self.OnContinuityChanged

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

    def setup_layout(self):
        term_cap = "Edge" if self.is_surface else "End"
        
        layout = ef.DynamicLayout()
        layout.Padding = ed.Padding(10)
        layout.Spacing = ed.Size(4, 4)

        key = 'bLinkedEnds'
        layout.AddSeparateRow(None, ed.Size(4, 4), False, False, (self.labels[key], self.radioButtonLists[key]))
        layout.AddRow(None)

        layout.AddRow(ef.Label(Text="Picked {} Configuration".format(term_cap), Font=ed.Font(ed.SystemFont.Bold, 10)))
        stepper_picked = ef.DynamicLayout()
        stepper_picked.Spacing = ed.Size(0, 0)
        stepper_picked.AddRow(self.btnScaleUp_Picked)
        stepper_picked.AddRow(self.btnScaleDown_Picked)

        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['fScale_Picked'], self.textBoxes['fScale_Picked'], stepper_picked,
            ef.Label(Width=20),
            self.labels['fIncrement'], self.textBoxes['fIncrement'],
            None
        ))
        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['fSlideG2_Picked'], self.numericSteppers['fSlideG2_Picked'],
            ef.Label(Width=20),
            self.labels['fSlideG3_Picked'], self.numericSteppers['fSlideG3_Picked'],
            None
        ))
        layout.AddRow(None)

        layout.AddRow(ef.Label(Text="Opposite {} Configuration".format(term_cap), Font=ed.Font(ed.SystemFont.Bold, 10)))
        stepper_opp = ef.DynamicLayout()
        stepper_opp.Spacing = ed.Size(0, 0)
        stepper_opp.AddRow(self.btnScaleUp_Opp)
        stepper_opp.AddRow(self.btnScaleDown_Opp)

        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['fScale_Opp'], self.textBoxes['fScale_Opp'], stepper_opp,
            None
        ))
        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['fSlideG2_Opp'], self.numericSteppers['fSlideG2_Opp'],
            ef.Label(Width=20),
            self.labels['fSlideG3_Opp'], self.numericSteppers['fSlideG3_Opp'],
            None
        ))
        layout.AddRow(None)

        layout.AddRow(ef.Label(Text="Analysis Tools", Font=ed.Font(ed.SystemFont.Bold, 10)))
        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['iGraphScale'], self.numericSteppers['iGraphScale'],
            ef.Label(Width=20),
            self.labels['iGraphDensity'], self.numericSteppers['iGraphDensity'],
            None
        ))
        layout.AddRow(None)

        cont_layout = ef.DynamicLayout()
        cont_layout.Spacing = ed.Size(4, 4)
        cont_layout.AddRow(self.labels['idxCont_Picked'], self.radioButtonLists['idxCont_Picked'])
        cont_layout.AddRow(self.labels['idxCont_Opp'], self.radioButtonLists['idxCont_Opp'])
        layout.AddSeparateRow(cont_layout)
        layout.AddRow(None)

        layout.AddSeparateRow(None, ed.Size(20, 5), False, False, (
            self.checkBoxes['bDeleteInput'], self.checkBoxes['bEcho'], self.checkBoxes['bDebug']
        ))

        self.ok_button = ef.Button(Text = 'OK')
        self.ok_button.Click += self.OnOKButtonClick
        save_button = ef.Button(Text = 'Save Settings')
        save_button.Click += self.OnSaveSettingsButtonClick
        self.abort_button = ef.Button(Text = 'Cancel')
        self.abort_button.Click += self.OnCancelButtonClick

        self.DefaultButton = self.ok_button
        self.AbortButton = self.abort_button

        button_stack = ef.StackLayout()
        button_stack.Orientation = ef.Orientation.Horizontal
        button_stack.Spacing = 8
        button_stack.Items.Add(self.ok_button)
        button_stack.Items.Add(save_button)
        button_stack.Items.Add(self.abort_button)

        layout.AddRow(None) 
        layout.AddSeparateRow(None, button_stack, None)

        self.Content = layout

    def SyncLinkedControls(self, sender=None, e=None):
        if bool(self.radioButtonLists['bLinkedEnds'].SelectedIndex):
            if self.textBoxes['fScale_Opp'].Text != self.textBoxes['fScale_Picked'].Text:
                self.textBoxes['fScale_Opp'].Text = self.textBoxes['fScale_Picked'].Text
            if self.numericSteppers['fSlideG2_Opp'].Value != self.numericSteppers['fSlideG2_Picked'].Value:
                self.numericSteppers['fSlideG2_Opp'].Value = self.numericSteppers['fSlideG2_Picked'].Value
            if self.numericSteppers['fSlideG3_Opp'].Value != self.numericSteppers['fSlideG3_Picked'].Value:
                self.numericSteppers['fSlideG3_Opp'].Value = self.numericSteppers['fSlideG3_Picked'].Value

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

        self.numericSteppers['fSlideG2_Picked'].Enabled = (scale_limit_P >= 3)
        self.numericSteppers['fSlideG3_Picked'].Enabled = (scale_limit_P >= 4)

        self.textBoxes['fScale_Opp'].Enabled = not is_linked
        self.btnScaleUp_Opp.Enabled = not is_linked
        self.btnScaleDown_Opp.Enabled = not is_linked
        
        if is_linked:
            self.numericSteppers['fSlideG2_Opp'].Enabled = False
            self.numericSteppers['fSlideG3_Opp'].Enabled = False
        else:
            self.numericSteppers['fSlideG2_Opp'].Enabled = (scale_limit_O >= 3)
            self.numericSteppers['fSlideG3_Opp'].Enabled = (scale_limit_O >= 4)

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
        
        # Identify which text box triggered the change and update its exact tracker
        if not self._auto_updating and val is not None:
            if sender == self.textBoxes['fScale_Picked']:
                self._exact_scale_picked = val
            elif sender == self.textBoxes['fScale_Opp']:
                self._exact_scale_opp = val
                
        if val is not None and val > Rhino.RhinoMath.ZeroTolerance:
            sender.BackgroundColor = ed.Colors.White
        else:
            sender.BackgroundColor = ed.Colors.LightPink
            
        self.SyncLinkedControls()
        self.UpdatePreview()

    def OnSlideValueChanged(self, sender, e):
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

    def StartHoldTimer(self, direction, end_name):
        self.active_stepper_end = end_name
        self.AdjustScale(direction)
        self.hold_direction = direction
        self.hold_timer.Start()

    def StopHoldTimer(self, sender, e):
        self.hold_timer.Stop()

    def OnHoldTimerElapsed(self, sender, e):
        self.AdjustScale(self.hold_direction)

    def AdjustScale(self, direction):
        incr_val = self.ParseToFloat(self.textBoxes['fIncrement'].Text)
        if incr_val is None: return

        if self.active_stepper_end == 'Picked':
            current_val = self._exact_scale_picked
            target_key = 'fScale_Picked'
        else:
            current_val = self._exact_scale_opp
            target_key = 'fScale_Opp'

        if current_val is not None:
            new_val = current_val + (incr_val * direction)
            if new_val > Rhino.RhinoMath.ZeroTolerance:
                if self.active_stepper_end == 'Picked':
                    self._exact_scale_picked = new_val
                else:
                    self._exact_scale_opp = new_val
                self._auto_updating = True
                self.textBoxes[target_key].Text = "{:g}".format(round(new_val, 4))
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
            self.numericSteppers['fSlideG2_Picked'].Increment = val
            self.numericSteppers['fSlideG3_Picked'].Increment = val
            self.numericSteppers['fSlideG2_Opp'].Increment = val
            self.numericSteppers['fSlideG3_Opp'].Increment = val
        else:
            sender.BackgroundColor = ed.Colors.LightPink

    def OnGraphScaleStepperChanged(self, sender, e):
        if hasattr(self, 'conduit') and self.conduit is not None:
            self.conduit.graph_scale = sender.Value
            sc.doc.Views.Redraw()

    def OnGraphDensityStepperChanged(self, sender, e):
        if hasattr(self, 'conduit') and self.conduit is not None:
            self.conduit.graph_density = int(sender.Value)
            sc.doc.Views.Redraw()

    def SaveSettings(self):
        sc.sticky[Opts.stickyKeys['iGraphScale']] = Opts.values['iGraphScale'] = int(self.numericSteppers['iGraphScale'].Value)
        sc.sticky[Opts.stickyKeys['iGraphDensity']] = Opts.values['iGraphDensity'] = int(self.numericSteppers['iGraphDensity'].Value)
        sc.sticky[Opts.stickyKeys['bLinkedEnds']] = Opts.values['bLinkedEnds'] = bool(self.radioButtonLists['bLinkedEnds'].SelectedIndex)
        
        parsed_incr = self.ParseToFloat(self.textBoxes['fIncrement'].Text)
        if parsed_incr is not None:
            sc.sticky[Opts.stickyKeys['fIncrement']] = Opts.values['fIncrement'] = parsed_incr

        for end in ('_Picked', '_Opp'):
            parsed_scale = self.ParseToFloat(self.textBoxes['fScale' + end].Text)
            if parsed_scale is not None:
                sc.sticky[Opts.stickyKeys['fScale' + end]] = Opts.values['fScale' + end] = parsed_scale
            sc.sticky[Opts.stickyKeys['fSlideG2' + end]] = Opts.values['fSlideG2' + end] = self.numericSteppers['fSlideG2' + end].Value
            sc.sticky[Opts.stickyKeys['fSlideG3' + end]] = Opts.values['fSlideG3' + end] = self.numericSteppers['fSlideG3' + end].Value
        sc.sticky[Opts.stickyKeys['idxCont_Picked']] = Opts.values['idxCont_Picked'] = self.radioButtonLists['idxCont_Picked'].SelectedIndex
        sc.sticky[Opts.stickyKeys['idxCont_Opp']] = Opts.values['idxCont_Opp'] = self.radioButtonLists['idxCont_Opp'].SelectedIndex
        sc.sticky[Opts.stickyKeys['bDeleteInput']] = Opts.values['bDeleteInput'] = self.checkBoxes['bDeleteInput'].Checked
        sc.sticky[Opts.stickyKeys['bEcho']] = Opts.values['bEcho'] = self.checkBoxes['bEcho'].Checked
        sc.sticky[Opts.stickyKeys['bDebug']] = Opts.values['bDebug'] = self.checkBoxes['bDebug'].Checked

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
        self.conduit.Enabled = False
        sc.doc.Views.Redraw()