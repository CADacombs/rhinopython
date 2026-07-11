#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script will translate the G1 [1] (1st from the end) [0] from the end [0] control point
by a scale (1 is no change).
For G2 setting of MaintainPicked or MaintainOpp, the respective G2 [2] control points will also be translated.
"""

"""
210303, 0307: Created.
260420-25, 0709-10: Added an optional dialog. Added a preview for the dialog.
        Refactored.

Tangential sliding: Translating the 3rd control point, p2, parallel to the tangent vector (p1 - p0).
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

    key = 'bBothEnds'; keys.append(key)
    values[key] = True
    names[key] = 'AdjustEnds'
    offValues[key] = 'Picked'
    onValues[key] = 'Both'
    riOpts[key] = ri.Custom.OptionToggle(values[key], offValues[key], onValues[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fScale'; keys.append(key)
    values[key] = 0.5
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fScaleIncrement'; keys.append(key) # Only for dialog.
    values[key] = 0.05
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCont_Picked'; keys.append(key)
    listValues[key] = 'None', 'G0', 'G1', 'G2', 'G3' # All items must be strings.
    values[key] = 4 # 4 represents G3
    names[key] = 'MaintainPicked'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCont_Opp'; keys.append(key)
    listValues[key] = 'None', 'G0', 'G1', 'G2', 'G3' # All items must be strings.
    values[key] = 4 # 4 represents G3
    names[key] = 'MaintainOpp'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteInput'; keys.append(key)
    values[key] = False
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
            print("What happened?")
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_CLI():
    """
    Get curve with picked end and optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick curve near an end")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve


    def geomFilter_Curve(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index)

        if isinstance(geom, rg.BrepEdge):
            # DuplicateCurve gets the edge as a curve, which may be a subset of the EdgeCurve.
            rgC = geom.DuplicateCurve()
        elif isinstance(geom, rg.Curve):
            rgC = geom
        else:
            return False

        if isinstance(rgC, rg.NurbsCurve):
            nc = rgC
        elif isinstance(rgC, rg.PolyCurve):
            if rgC.SegmentCount > 1:
                print("PolyCurve with multiple segments is ignored.")
                return False
            nc = rgC.ToNurbsCurve()
        else:
            return False

        if nc.IsPeriodic:
            print("Periodic curves are not supported.")
            return False

        if nc.Degree == 1:
            print("Ignorded degree 1 NURBS curve.")
            return False

        #if nc.Degree + 1 < nc.Points.Count:
        #    print("Ignored non-Bezier curve."
        #    return


        return True


    go.SetCustomGeometryFilter(geomFilter_Curve)

    go.DisablePreSelect()

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('bGUI')
        if not Opts.values['bGUI']:
            addOption('bBothEnds')
            addOption('fScale')
            addOption('iCont_Picked')
            addOption('iCont_Opp')
            addOption('bDeleteInput')
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
            key = 'fScale'
            goNumber = go.Number()
            if goNumber == 1:
                print("A scale of 1 will not modify the curve. Scale was not modified.")
                continue
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createCurve(nc_In, fScale, iEndToScale=2, iG_T0=2, iG_T1=2, bDebug=False):
    """
    Parameters:
        nc_In: rg.NurbsCurve
        fScale: float of scale factor
        iEndToScale: int(0 for T0, 1 for T1, or 2 for both ends of nc_In)
        iG_T0: int(-1 for No continuity, 0 for G0, 1 for G1, or 2 for G2)
        iG_T1: int(-1 for No continuity, 0 for G0, 1 for G1, or 2 for G2)
        bDebug: bool
    Returns on success: rg.NurbsCurve
    Returns on fail: None
    """

    if iG_T0 is None and iG_T1 is None: return

    if nc_In.IsPeriodic: return
    if not isinstance(nc_In, rg.NurbsCurve): return
    #if nc_In.Degree + 1 < nc_In.Points.Count: return # because not Bezier.

    if nc_In.Points.Count < (iG_T0 + 1) + (iG_T1 + 1):
        print("Curve needs {} more points to maintain continuities.".format(
            (iG_T0 + 1) + (iG_T1 + 1) - nc_In.Points.Count))
        return

    pts_Prime = []


    # ----------------------------------------------------
    # SCALE T0 END
    # ----------------------------------------------------
    if iEndToScale in (0,2):
        if iG_T0 in (0,1,2,3):
            p0 = nc_In.Points[0].Location
            pts_Prime.append(p0)

            if iG_T0 in (1,2,3):
                p1 = nc_In.Points[1].Location
                xform = rg.Transform.Scale(p0, fScale)
                p1p = rg.Point3d(p1)
                p1p.Transform(xform)
                pts_Prime.append(p1p)

                if iG_T0 in (2,3):
                    p2 = nc_In.Points[2].Location
                    m2 = ((p1p - p0).Length/(p1 - p0).Length)**2.0
                    p2p = 2.0*p1p + -p0 + m2*(-2.0*p1 + p2 + p0)
                    pts_Prime.append(p2p)

                    if iG_T0 == 3:
                        p3 = nc_In.Points[3].Location
                        m3 = ((p1p - p0).Length/(p1 - p0).Length)**3.0
                        p3p = 3.0*p2p - 3.0*p1p + p0 + m3*(p3 - 3.0*p2 + 3.0*p1 - p0)
                        pts_Prime.append(p3p)


    # ----------------------------------------------------
    # DUPLICATE UNAFFECTED MIDDLE POINTS
    # ----------------------------------------------------
    iCt_PtsNotNeededByT1Scale = nc_In.Points.Count - (iG_T1 + 1)
    for i in range(len(pts_Prime), iCt_PtsNotNeededByT1Scale):
        pts_Prime.append(nc_In.Points[i].Location)


    # ----------------------------------------------------
    # SCALE T1 END
    # ----------------------------------------------------
    if iEndToScale in (1,2):
        if iG_T1 in (0,1,2,3):
            pts_New_T1 = [] # To be reversed and extended on pts_Prime.
            p0 = nc_In.Points[nc_In.Points.Count-1].Location
            pts_New_T1.append(p0)

            if iG_T1 in (1,2,3):
                p1 = nc_In.Points[nc_In.Points.Count-2].Location
                xform = rg.Transform.Scale(p0, fScale)
                p1p = rg.Point3d(p1) # G1 CP location prime.
                p1p.Transform(xform)
                pts_New_T1.append(p1p)

                if iG_T1 in (2,3):
                    p2 = nc_In.Points[nc_In.Points.Count-3].Location
                    m2 = ((p1p - p0).Length/(p1 - p0).Length)**2.0
                    p2p = 2.0*p1p + -p0 + m2*(-2.0*p1 + p2 + p0)
                    pts_New_T1.append(p2p)

                    if iG_T1 == 3:
                        p3 = nc_In.Points[nc_In.Points.Count-4].Location
                        m3 = ((p1p - p0).Length/(p1 - p0).Length)**3.0
                        p3p = 3.0*p2p - 3.0*p1p + p0 + m3*(p3 - 3.0*p2 + 3.0*p1 - p0)
                        pts_New_T1.append(p3p)

            pts_New_T1.reverse()
            pts_Prime.extend(pts_New_T1)


    # Enforce minimum distance (1e-6 cm) converted to current document units
    unit_scale = Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Centimeters,
        sc.doc.ModelUnitSystem)
    min_dist = 1e-6 * unit_scale

    for i in range(len(pts_Prime) - 1):
        if pts_Prime[i].DistanceTo(pts_Prime[i+1]) < min_dist:
            sReport = "Minimum control point distance (1e-6 cm) violated. Is Scale too small?"
            Rhino.RhinoApp.SetCommandPromptMessage(sReport)
            if bDebug: print(sReport)
            return None


    nc_Out = nc_In.Duplicate()

    for i in range(nc_Out.Points.Count):
        nc_Out.Points.SetPoint(
            index=i,
            point=pts_Prime[i],
            weight=nc_In.Points.GetWeight(i))

    if bDebug:
        for i in range(len(pts_Prime)):
            pt = pts_Prime[i]
            rgDot = rg.TextDot("{}".format(i), pt)
            rgDot.FontHeight = 11
            sc.doc.Objects.AddTextDot(rgDot)
        sc.doc.Views.Redraw()

    return nc_Out


class EndBulgePreviewConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        super(EndBulgePreviewConduit, self).__init__()
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.crv = None
        self.graph_scale = Opts.values['iGraphScale']
        self.graph_density = Opts.values['iGraphDensity']

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if self.crv is None:
            return
        bbox = self.crv.GetBoundingBox(accurate=False) # False uses control points, thus the convex hull.
        bbox.Inflate(sc.doc.ModelAbsoluteTolerance * 100)
        calculateBoundingBoxEventArgs.IncludeBoundingBox(bbox)

    def PostDrawObjects(self, drawEventArgs):
        if self.crv is None:
            return

        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

        drawEventArgs.Display.DrawCurve(
            curve=self.crv,
            color=self.color,
            thickness=crv_thk
        )

        self.cp_locations = [pt.Location for pt in self.crv.Points]

        drawEventArgs.Display.DrawPatternedPolyline(
            points=self.cp_locations,
            color=self.color,
            pattern=0x00001111,
            thickness=1,
            close=False)

        drawEventArgs.Display.DrawPoints(
            points=self.cp_locations,
            style=Rhino.Display.PointStyle.Simple,
            radius=3,
            color=self.color
            )

        drawEventArgs.Display.DrawCurvatureGraph(
            curve=self.crv,
            color=self.color,
            hairScale=self.graph_scale,
            hairDensity=self.graph_density,
            sampleDensity=2
            )


class EtoDialog(ef.Dialog):
    """
    Reference spb_Intersect_SrfSrf.py and
    https://github.com/mcneel/rhino-developer-samples/blob/3179a8386a64602ee670cc832c77c561d1b0944b/rhinopython/SampleEtoModelessForm.py
    """

    def __init__(self, objref_In):
        self.Title = "EndBulge by Scale"
        self.objref_In = objref_In
        self.dialog_ok = False

        rgC_In = objref_In.Curve()
        if isinstance(rgC_In, rg.BrepEdge):
            self.nc_In = rgC_In.ToNurbsCurve()
        else:
            self.nc_In = rgC_In.ToNurbsCurve() # TODO: Review whether non-Nurbs should be rejected. Handles Curve and PolyCurve

        self._exact_scale = Opts.values['fScale'] 
        self._auto_updating = False

        self.create_controls()
        self.setup_layout()

        # --- FIX: FORCE INITIAL BACKGROUND VALIDATION ON STARTUP ---
        # Run the validation check on startup so if fScale is 1.0, 
        # the text box immediately turns pink.
        initial_scale = self.ParseToFloat(self.textBoxes['fScale'].Text)
        if initial_scale is not None:
            if (
                initial_scale <= Rhino.RhinoMath.ZeroTolerance or 
                abs(initial_scale - 1.0) <= Rhino.RhinoMath.ZeroTolerance
            ):
                self.textBoxes['fScale'].BackgroundColor = ed.Colors.LightPink
            else:
                self.textBoxes['fScale'].BackgroundColor = ed.Colors.White


    def UpdatePreview(self):
        """Calculates the modified curve based on UI values and updates the conduit."""
        if not hasattr(self, 'conduit') or self.conduit is None:
            return

        fScale = self.ParseToFloat(self.textBoxes['fScale'].Text)

        # Handle invalid inputs (letters, empty strings, negative numbers/zero)
        if (
            fScale is None or
            fScale <= Rhino.RhinoMath.ZeroTolerance
        ):
            self.conduit.crv = None
            sc.doc.Views.Redraw()
            return

        if abs(fScale - 1.0) <= Rhino.RhinoMath.ZeroTolerance:
            self.conduit.crv = self.nc_In.Duplicate()
            sc.doc.Views.Redraw()
            return

        bBothEnds = bool(self.radioButtonLists['bBothEnds'].SelectedIndex)
        iCont_Picked = self.radioButtonLists['iCont_Picked'].SelectedIndex
        iCont_Opp = self.radioButtonLists['iCont_Opp'].SelectedIndex
        bDebug = self.checkBoxes['bDebug'].Checked


        bSuccess, t_AtPicked = self.nc_In.ClosestPoint(self.objref_In.SelectionPoint())
        if not bSuccess:
            return

        if t_AtPicked > self.nc_In.Domain.Mid:
            iEndToScale = 2 if bBothEnds else 1
            iG_T0, iG_T1 = iCont_Opp-1, iCont_Picked-1
        else:
            iEndToScale = 2 if bBothEnds else 0
            iG_T0, iG_T1 = iCont_Picked-1, iCont_Opp-1

        nc_Res = createCurve(
            nc_In=self.nc_In,
            fScale=fScale,
            iEndToScale=iEndToScale,
            iG_T0=iG_T0,
            iG_T1=iG_T1,
            bDebug=bDebug
        )

        self.conduit.crv = nc_Res
        sc.doc.Views.Redraw()


    def create_controls(self):
        self.labels = {}
        self.checkBoxes = {}
        self.radioButtonLists = {}
        self.numericSteppers = {}
        self.textBoxes = {}

        key = 'bBothEnds'
        self.labels[key] = ef.Label(Text = "Adjust end(s):")
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(16, 4)
        self.radioButtonLists[key].DataStore = (Opts.offValues[key], Opts.onValues[key])
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[int(Opts.values[key])]


        key = 'fScale'
        #        self.labels[key] = ef.Label(Text = "Scale:")
        #        self.numericSteppers[key] = ef.NumericStepper(
        #            MinValue = 0.01,
        #            DecimalPlaces = 3,
        #            MaximumDecimalPlaces = 3,
        #            Increment = 0.01,
        #            Value = Opts.values[key],
        #            )
        self.labels[key] = ef.Label(Text = "Scale:")
        self.textBoxes[key] = ef.TextBox()
        self.textBoxes[key].Text = str(Opts.values[key])
        self.textBoxes[key].TextChanged += self.OnScaleTextChanged

        self.hold_timer = ef.UITimer()
        self.hold_timer.Interval = 0.15 # Repeats every 150 milliseconds
        self.hold_timer.Elapsed += self.OnHoldTimerElapsed
        self.hold_direction = 0

        # Custom stepper buttons
        self.btnScaleUp = ef.Button(Text=unichr(9650), Width=16, Height=12)
        self.btnScaleDown = ef.Button(Text=unichr(9660), Width=16, Height=12)
        small_font = ed.Font(ed.SystemFont.Default, 4)
        self.btnScaleUp.Font = small_font
        self.btnScaleDown.Font = small_font
        self.btnScaleUp.MinimumSize = ed.Size(16, 12)
        self.btnScaleDown.MinimumSize = ed.Size(16, 12)
        
        self.btnScaleUp.MouseDown += lambda s, e: self.StartHoldTimer(1)
        self.btnScaleUp.MouseUp += self.StopHoldTimer
        self.btnScaleUp.MouseLeave += self.StopHoldTimer

        self.btnScaleDown.MouseDown += lambda s, e: self.StartHoldTimer(-1)
        self.btnScaleDown.MouseUp += self.StopHoldTimer
        self.btnScaleDown.MouseLeave += self.StopHoldTimer


        key = 'fScaleIncrement'
        self.labels[key] = ef.Label(Text = "Incr.:")
        self.textBoxes[key] = ef.TextBox()
        self.textBoxes[key].Text = str(Opts.values[key])
        self.textBoxes[key].TextChanged += self.OnIncrementTextChanged


        key = 'iGraphScale'
        self.labels[key] = ef.Label(Text = "Graph scale:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 0 # Forces integer steps to match Rhino
        self.numericSteppers[key].MinValue = 1.0
        self.numericSteppers[key].MaxValue = 10000.0
        self.numericSteppers[key].Increment = 1.0
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnGraphScaleStepperChanged

        key = 'iGraphDensity'
        self.labels[key] = ef.Label(Text = "Graph density:")
        self.numericSteppers[key] = ef.NumericStepper()
        self.numericSteppers[key].DecimalPlaces = 0
        self.numericSteppers[key].MinValue = 1.0
        self.numericSteppers[key].MaxValue = 100.0
        self.numericSteppers[key].Increment = 1.0
        self.numericSteppers[key].Value = float(Opts.values[key])
        self.numericSteppers[key].ValueChanged += self.OnGraphDensityStepperChanged


        key = 'iCont_Picked'
        self.labels[key] = ef.Label(Text = "Cont. of picked end:    ")
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(4, 4)
        self.radioButtonLists[key].DataStore = (Opts.listValues[key])
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[int(Opts.values[key])]

        key = 'iCont_Opp'
        self.labels[key] = ef.Label(Text = "Cont. of opp. end:")
        self.radioButtonLists[key] = ef.RadioButtonList()
        self.radioButtonLists[key].Spacing = ed.Size(4, 4)
        self.radioButtonLists[key].DataStore = (Opts.listValues[key])
        self.radioButtonLists[key].SelectedValue = self.radioButtonLists[key].DataStore[int(Opts.values[key])]

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

        self.radioButtonLists['bBothEnds'].SelectedIndexChanged += lambda s, e: self.UpdatePreview()
        self.radioButtonLists['iCont_Picked'].SelectedIndexChanged += lambda s, e: self.UpdatePreview()
        self.radioButtonLists['iCont_Opp'].SelectedIndexChanged += lambda s, e: self.UpdatePreview()
        self.checkBoxes['bDebug'].CheckedChanged += lambda s, e: self.UpdatePreview()


    def setup_layout(self):
        layout = ef.DynamicLayout()
        layout.Padding = ed.Padding(10)
        layout.Spacing = ed.Size(4, 4)

        key = 'bBothEnds'
        layout.AddSeparateRow(None, ed.Size(4, 4), False, False, (self.labels[key], self.radioButtonLists[key]))


        # Stack the tiny buttons vertically
        stepper_layout = ef.DynamicLayout()
        stepper_layout.Spacing = ed.Size(0, 0)
        stepper_layout.AddRow(self.btnScaleUp)
        stepper_layout.AddRow(self.btnScaleDown)

        key = 'fScale'
        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['fScale'], 
            self.textBoxes['fScale'], 
            stepper_layout,
            ef.Label(Width=20),
            self.labels['fScaleIncrement'], 
            self.textBoxes['fScaleIncrement'], 
            None
        ))

        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (
            self.labels['iGraphScale'],
            self.numericSteppers['iGraphScale'],
            ef.Label(Width=20),
            self.labels['iGraphDensity'],
            self.numericSteppers['iGraphDensity'],
            None
        ))


        cont_layout = ef.DynamicLayout()
        cont_layout.Spacing = ed.Size(4, 4)

        key = 'iCont_Picked'
        cont_layout.AddRow(self.labels[key], self.radioButtonLists[key])

        key = 'iCont_Opp'
        cont_layout.AddRow(self.labels[key], self.radioButtonLists[key])

        layout.AddSeparateRow(cont_layout)

        layout.AddRow(None)

        layout.AddSeparateRow(None, ed.Size(20, 5), False, False, (self.checkBoxes['bDeleteInput'], self.checkBoxes['bEcho'], self.checkBoxes['bDebug']))

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

        layout.AddRow(None) # Top spacing
        layout.AddSeparateRow(None, button_stack, None)

        self.Content = layout


    def ParseToFloat(self, text):
        """Helper to safely parse both fractions and decimals from a string."""
        text = text.strip()
        try:
            if '/' in text:
                num, den = text.split('/')
                return float(num) / float(den)
            return float(text)
        except (ValueError, ZeroDivisionError):
            return None


    def OnScaleTextChanged(self, sender, e):
        val = self.ParseToFloat(sender.Text)

        if not self._auto_updating and val is not None:
            self._exact_scale = val

        # Base validation: Must be a number > 0
        is_valid = val is not None and val > Rhino.RhinoMath.ZeroTolerance

        # 1.0 is an invalid scale factor
        if is_valid and abs(val - 1.0) <= Rhino.RhinoMath.ZeroTolerance:
            is_valid = False

        if is_valid:
            sender.BackgroundColor = ed.Colors.White
        else:
            sender.BackgroundColor = ed.Colors.LightPink

        self.UpdatePreview() # For conduit.


    def StartHoldTimer(self, direction):
        """Fires immediately on first click and starts the repeat timer."""
        self.AdjustScale(direction)
        self.hold_direction = direction
        self.hold_timer.Start()


    def StopHoldTimer(self, sender, e):
        """Kills the repeat timer when the mouse is released or drifts away."""
        self.hold_timer.Stop()


    def OnHoldTimerElapsed(self, sender, e):
        """Fires continuously while the timer is running."""
        self.AdjustScale(self.hold_direction)


    def AdjustScale(self, direction):
        current_val = self._exact_scale
        incr_val = self.ParseToFloat(self.textBoxes['fScaleIncrement'].Text)

        if current_val is not None and incr_val is not None:
            new_val = current_val + (incr_val * direction)
            if new_val > Rhino.RhinoMath.ZeroTolerance:

                # Save the perfect floating-point math internally
                self._exact_scale = new_val

                # Lock the text box, update the display, then unlock it
                self._auto_updating = True
                self.textBoxes['fScale'].Text = "{:g}".format(round(new_val, 4))
                self._auto_updating = False

        self.UpdatePreview() # For conduit.


    def OnIncrementTextChanged(self, sender, e):
        text = sender.Text.strip()
        val = None

        try:
            # Handle fractions (e.g., "1/8")
            if '/' in text:
                num, den = text.split('/')
                val = float(num) / float(den)
            # Handle standard decimals (e.g., "0.125")
            else:
                val = float(text)
        except (ValueError, ZeroDivisionError):
            pass # Ignore errors while the user is mid-typing

        # Validate result and apply it to the NumericStepper
        if val is not None and val > Rhino.RhinoMath.ZeroTolerance:
            #self.numericSteppers['fScale'].Increment = val
            sender.BackgroundColor = ed.Colors.White # Clear warning
        else:
            # Provide visual feedback that the current text is invalid
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

        key = 'iGraphScale'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = int(self.numericSteppers[key].Value)

        key = 'iGraphDensity'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = int(self.numericSteppers[key].Value)

        key = 'bBothEnds'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = bool(self.radioButtonLists[key].SelectedIndex)

        key = 'fScale'
        parsed_scale = self.ParseToFloat(self.textBoxes[key].Text)
        if parsed_scale is not None:
            sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = parsed_scale

        key = 'fScaleIncrement'
        parsed_incr = self.ParseToFloat(self.textBoxes[key].Text)
        if parsed_incr is not None:
            sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = parsed_incr

        key = 'iCont_Picked'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = self.radioButtonLists[key].SelectedIndex

        key = 'iCont_Opp'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = self.radioButtonLists[key].SelectedIndex

        key = 'bDeleteInput'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = self.checkBoxes[key].Checked

        key = 'bEcho'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = self.checkBoxes[key].Checked

        key = 'bDebug'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = self.checkBoxes[key].Checked

        #        for key in Opts.keys:
        #            if key == 'bGUI': continue
        #            print(sc.sticky[Opts.stickyKeys[key]])


    def OnOKButtonClick(self, sender, e):

        fScale = self.ParseToFloat(self.textBoxes['fScale'].Text)

        # If the scale is invalid or exactly 1.0, treat it as a Cancel action 
        # so we don't output or bake a redundant duplicate curve.
        if (
            fScale is None or 
            fScale <= Rhino.RhinoMath.ZeroTolerance or 
            abs(fScale - 1.0) <= Rhino.RhinoMath.ZeroTolerance
        ):
            print("Scale is 1.0 or invalid. No changes were applied.")
            self.dialog_ok = False
            self.Close()
            return

        self.SaveSettings()
        self.dialog_ok = True
        self.Close()


    def OnSaveSettingsButtonClick(self, sender, e):
            self.SaveSettings()
            # Optional: Print a status message so the user knows it worked
            print("Settings saved as default.")


    def OnCancelButtonClick(self, sender, e):
        self.Result = ef.DialogResult.Cancel
        self.Close()


    def OnFormClosed(self, sender, e):
        self.conduit.Enabled = False
        sc.doc.Views.Redraw()


def _isScaleOK(fScale):

    if fScale is None:
        return False, "Scale is not set. Script canceled."

    if fScale == 0:
        return False, "Scale set to 0, which would result in stacked control points. Script canceled."

    if fScale < 0:
        return False, "Scale set to a negative value. Script canceled.".format(fScale)

    if abs(fScale) <= Rhino.RhinoMath.ZeroTolerance:
        return False, "Scale set to {}, which would result in almost stacked control points. Script canceled.".format(fScale)

    if fScale == 1:
        return False, "Scale set to 1, so curve cannot be modified. Script canceled."

    if abs(fScale - 1.0) <= Rhino.RhinoMath.ZeroTolerance:
        return False, "Scale set within {} of 1.0. Script canceled.".format(Rhino.RhinoMath.ZeroTolerance)

    return True, None


def _createCurve_viaGUI(objref_In):

    key = 'conduit'
    stickyKey = '{}({})'.format(key, __file__)
    if stickyKey in sc.sticky:
        old_conduit = sc.sticky[stickyKey]
        old_conduit.Enabled = False
        sc.doc.Views.Redraw()

    parent = Rhino.UI.RhinoEtoApp.MainWindowForDocument(sc.doc)
    dialog = EtoDialog(objref_In)

    dialog.conduit = EndBulgePreviewConduit()
    sc.sticky[stickyKey] = dialog.conduit

    # Run an initial preview update before opening the dialog
    dialog.UpdatePreview()
    dialog.conduit.Enabled = True
    sc.doc.Views.Redraw()

    Rhino.UI.EtoExtensions.ShowSemiModal(dialog, sc.doc, parent)

    dialog.conduit.Enabled = False
    sc.doc.Views.Redraw()

    if not dialog.dialog_ok:
        return

    return dialog.conduit.crv


def processCurveObject(objref_In, nc_Precalc=None, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bBothEnds = getOpt('bBothEnds')
    fScale = getOpt('fScale')
    iCont_Picked = getOpt('iCont_Picked')
    iCont_Opp = getOpt('iCont_Opp')
    bDeleteInput = getOpt('bDeleteInput')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if iCont_Picked == iCont_Opp == 0:
        print("Both continuities are set to 0, so curve cannot be modified. Script canceled.")
        return

    bOkay, sReport = _isScaleOK(fScale)
    if not bOkay:
        print(sReport)
        return


    if nc_Precalc is not None:
        nc_Res = nc_Precalc
    else:
        rgC_In = objref_In.Curve()

        if isinstance(rgC_In, rg.BrepEdge):
            nc_In = rgC_In.ToNurbsCurve()
        elif isinstance(rgC_In, rg.NurbsCurve):
            nc_In = rgC_In
        elif isinstance(rgC_In, rg.PolyCurve):
            nc_In = rgC_In.ToNurbsCurve()
        else:
            return


        bSuccess, t_AtPicked = nc_In.ClosestPoint(objref_In.SelectionPoint())
        if not bSuccess:
            return


        if t_AtPicked > nc_In.Domain.Mid:
            iEndToScale = 2 if bBothEnds else 1
            iG_T0, iG_T1 = iCont_Opp, iCont_Picked
        else:
            iEndToScale = 2 if bBothEnds else 0
            iG_T0, iG_T1 = iCont_Picked, iCont_Opp


        nc_Res = createCurve(
            nc_In=nc_In,
            fScale=fScale,
            iEndToScale=iEndToScale,
            iG_T0=iG_T0,
            iG_T1=iG_T1,
            bDebug=bDebug,
            )

        if nc_Res is None:
            print("Curve could not be created.")
            return


    if not bDeleteInput or objref_In.Edge():
        gC_Out = sc.doc.Objects.AddCurve(nc_Res)
        if gC_Out == gC_Out.Empty:
            print("Could not add curve.")
        else:
            print("Curve was added.")
    else:
        if sc.doc.Objects.Replace(objref_In.ObjectId, nc_Res):
            gC_Out = objref_In.ObjectId
            print("Replaced curve.")
        else:
            print("Could not replace curve.")

    return gC_Out


def main():

    rv = getInput_CLI()
    if rv is None: return

    objref_In = rv

    bGUI = Opts.values['bGUI']

    if not bGUI:
        nc_Res = None
    else:
        Rhino.RhinoApp.SetCommandPromptMessage("Continuing in dialog...")
        nc_Res = _createCurve_viaGUI(objref_In)

        if nc_Res is None: return

    bBothEnds = Opts.values['bBothEnds']
    fScale = Opts.values['fScale']
    iCont_Picked = Opts.values['iCont_Picked'] - 1
    iCont_Opp = Opts.values['iCont_Opp'] - 1
    bDeleteInput = Opts.values['bDeleteInput']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if iCont_Picked < 1 and iCont_Opp < 1:
        print("Continuity of at least one end of curve must be G1 or G2. Script canceled.")
        return

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObject(
        objref_In=objref_In,
        nc_Precalc=nc_Res,
        bBothEnds=bBothEnds,
        fScale=fScale,
        iCont_Picked=iCont_Picked,
        iCont_Opp=iCont_Opp,
        bDeleteInput=bDeleteInput,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if gC_Res is None:
        return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()