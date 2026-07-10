#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script will translate the G1 [1] (1st from the end) [0] from the end [0] control point
by a scale (1 is no change).
For G2 setting of MaintainPicked or MaintainOpp, the respective G2 [2] control points will also be translated.
"""

"""
210303, 0307: Created.
260420-25, 0907: WIP Adding preview and number slider to dynamically change preview
        before accepting value.
        Refactoring.
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
    values[key] = 0.01
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCont_Picked'; keys.append(key)
    listValues[key] = 'None', 'G0', 'G1', 'G2' # All items must be strings.
    values[key] = 3
    names[key] = 'MaintainPicked'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCont_Opp'; keys.append(key)
    listValues[key] = 'None', 'G0', 'G1', 'G2' # All items must be strings.
    values[key] = 3
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


class DrawCurvesConduit_OLD(Rhino.Display.DisplayConduit):

    def __init__(self):
        #self.color = sc.doc.Layers.CurrentLayer.Color
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.crv = None
        self.pts = None

    #def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
    #    if len(self.crv) > 0:
    #        self.bbox = self.brep.GetBoundingBox(accurate=False)
    #        calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):

        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

        drawEventArgs.Display.DrawCurve(
            curve=self.crv,
            color=self.color,
            thickness=crv_thk)

        drawEventArgs.Display.DrawPoints(
            points=self.pts,
            style=Rhino.Display.PointStyle.Simple,
            radius=4,
            color=self.color)

        #for p in self.pts:
        #    drawEventArgs.Display.DrawPoint(
        #        curve=c,
        #        color=self.color,
        #        thickness=crv_thk)


class EndBulgePreviewConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        super(EndBulgePreviewConduit, self).__init__()
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.crv = None # Holds the modified preview curve

    def PreDrawObjects(self, drawEventArgs):
        if self.crv is None:
            return

        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

        drawEventArgs.Display.DrawCurve(
            curve=self.crv,
            color=self.color,
            thickness=crv_thk
        )

        cp_locations = [pt.Location for pt in self.crv.Points]
        
        drawEventArgs.Display.DrawPolyline(
            points=cp_locations,
            color=self.color,
            thickness=1,
            lineStyle=Rhino.Display.LineStyle.Dot # Dotted line looks cleaner for polygons
        )

        # For CPs.
        # drawEventArgs.Display.DrawPoints(
        #     points=cp_locations,
        #     style=Rhino.Display.PointStyle.Simple,
        #     radius=3,
        #     color=self.color
        # )


class EtoDialog(ef.Dialog):
    """
    Reference spb_Intersect_SrfSrf.py and
    https://github.com/mcneel/rhino-developer-samples/blob/3179a8386a64602ee670cc832c77c561d1b0944b/rhinopython/SampleEtoModelessForm.py
    """


    def __init__(self, objref_In):
        self.Title = "EndBulge by Scale"
        self.objref_In = objref_In
        
        rgC_In = objref_In.Curve()
        if isinstance(rgC_In, rg.BrepEdge):
            self.nc_In = rgC_In.ToNurbsCurve()
        else:
            self.nc_In = rgC_In.ToNurbsCurve() # TODO: Review whether non-Nurbs should be rejected. Handles Curve and PolyCurve

        self.create_controls()
        self.setup_layout()


    def UpdatePreview(self):
        """Calculates the modified curve based on UI values and updates the conduit."""
        if not hasattr(self, 'conduit') or self.conduit is None:
            return

        # 1. Parse current numeric inputs
        fScale = self.ParseToFloat(self.textBoxes['fScale'].Text)
        
        # If input is mid-typing or invalid (like 1.0), don't break; use an immutable fallback or skip
        if fScale is None or abs(fScale - 1.0) <= Rhino.RhinoMath.ZeroTolerance or fScale <= 0:
            self.conduit.crv = None
            sc.doc.Views.Redraw()
            return

        # 2. Extract configuration states from UI controls
        bBothEnds = bool(self.radioButtonLists['bBothEnds'].SelectedIndex)
        iCont_Picked = self.radioButtonLists['iCont_Picked'].SelectedIndex
        iCont_Opp = self.radioButtonLists['iCont_Opp'].SelectedIndex
        bDebug = self.checkBoxes['bDebug'].Checked

        # 3. Determine which end was picked (Logic adapted from your main script)
        bSuccess, t_AtPicked = self.nc_In.ClosestPoint(self.objref_In.SelectionPoint())
        if not bSuccess:
            return

        if t_AtPicked > self.nc_In.Domain.Mid:
            iEndToScale = 2 if bBothEnds else 1
            iG_T0, iG_T1 = iCont_Opp, iCont_Picked
        else:
            iEndToScale = 2 if bBothEnds else 0
            iG_T0, iG_T1 = iCont_Picked, iCont_Opp

        # 4. Generate the new preview curve using your core math function
        nc_Res = createCurve(
            nc_In=self.nc_In,
            fScale=fScale,
            iEndToScale=iEndToScale,
            iG_T0=iG_T0,
            iG_T1=iG_T1,
            bDebug=bDebug
        )

        # 5. Push to conduit and force viewport update
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

        # Custom stepper buttons
        self.btnScaleUp = ef.Button(Text=unichr(9650), Width=16, Height=12)
        self.btnScaleDown = ef.Button(Text=unichr(9660), Width=16, Height=12)
        small_font = ed.Font(ed.SystemFont.Default, 4)
        self.btnScaleUp.Font = small_font
        self.btnScaleDown.Font = small_font
        self.btnScaleUp.MinimumSize = ed.Size(16, 12)
        self.btnScaleDown.MinimumSize = ed.Size(16, 12)
        self.btnScaleUp.Click += self.OnScaleUpClick
        self.btnScaleDown.Click += self.OnScaleDownClick

        key = 'fScaleIncrement'
        self.labels[key] = ef.Label(Text = "Incr:")
        self.textBoxes[key] = ef.TextBox()
        self.textBoxes[key].Text = str(Opts.values[key])
        self.textBoxes[key].TextChanged += self.OnIncrementTextChanged

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
            stepper_layout,             # <--- Inserted fake stepper arrows
            self.labels['fScaleIncrement'], 
            self.textBoxes['fScaleIncrement'], 
            None
        ))

        #        key = 'fScale'
        #        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (self.labels['fScale'], self.numericSteppers['fScale'], self.labels['fScaleIncrement'], self.textBoxes['fScaleIncrement'], None))

#        key = 'iCont_Picked'
#        layout.AddSeparateRow(None, ed.Size(4, 4), True, False, (self.labels[key], self.radioButtonLists[key]))


        cont_layout = ef.DynamicLayout()
        cont_layout.Spacing = ed.Size(4, 4)

        key = 'iCont_Picked'
        cont_layout.AddRow(self.labels[key], self.radioButtonLists[key])

        key = 'iCont_Opp'
        cont_layout.AddRow(self.labels[key], self.radioButtonLists[key])

        layout.AddSeparateRow(cont_layout)

        layout.AddRow(None)

        layout.AddSeparateRow(None, ed.Size(20, 5), False, False, (self.checkBoxes['bDeleteInput'], self.checkBoxes['bEcho'], self.checkBoxes['bDebug']))

        ok_button = ef.Button(Text = 'OK')
        ok_button.Click += self.OnOKButtonClick

        save_button = ef.Button(Text = 'Save Settings')
        save_button.Click += self.OnSaveSettingsButtonClick

        abort_button = ef.Button(Text = 'Cancel')
        abort_button.Click += self.OnCancelButtonClick


        button_stack = ef.StackLayout()
        button_stack.Orientation = ef.Orientation.Horizontal
        button_stack.Spacing = 8

        button_stack.Items.Add(ok_button)
        button_stack.Items.Add(save_button)
        button_stack.Items.Add(abort_button)

        layout.AddRow(None) # Top spacing
        layout.AddSeparateRow(None, button_stack, None)

        self.Content = layout


    def addCurveSetListItem(self, cs, fTol):
        item = ef.ListItem()
        total_length = 0.0
        span_count = 0
        cp_count = 0
        for c in cs:
            total_length += c.GetLength()
            nc = c.ToNurbsCurve()
            span_count += nc.SpanCount
            cp_count += nc.Points.Count
            nc.Dispose()
        iPrec = sc.doc.ModelDistanceDisplayPrecision + 1
        item.Text = ""
        item.Text += "{:<10.6f}".format(fTol)
        item.Text += "{:<{:}.{:}f}".format(total_length, iPrec+8, iPrec)
        item.Text += "{:<11}".format(len(cs))
        item.Text += "{:<11}".format(span_count)
        item.Text += "{:<6}".format(cp_count)
        self.listbox.Items.Add(item)


    def fillListBox(self):
        for i in range(len(self.cs_nested)):
            self.addCurveSetListItem(self.cs_nested[i], self.fTols[i])
        self.listbox.SelectedIndex = 0

    def addDataToForm(self, d4_cs, d4_fTols):
        self.d4_cs = d4_cs
        self.d4_fTols = d4_fTols
        self.bUnder_A = self.bUnder_B = False
        self.cs_nested = self.d4_cs[(self.bUnder_A, self.bUnder_B)]
        self.fTols = self.d4_fTols[(self.bUnder_A, self.bUnder_B)]
        self.fillListBox()

    def CreateFormControls(self):
        layout = ef.DynamicLayout()
        layout.Padding = ed.Padding(10)
        layout.Spacing = ed.Size(4, 4)








        layout.Rows.Add(self.CreateCheckBoxListRow())

        layout.AddRow(None)

        iPrec = sc.doc.ModelDistanceDisplayPrecision + 1

        sLabel = "   {:<14}{:<12}{:<10}{:<8}{:<6}".format(
            'Tol', 'Length', 'Crvs', 'Spans', 'CPs')

        layout.Rows.Add(ef.Label(Text=sLabel))
        layout.Rows.Add(self.CreateListBoxRow())

        ok_button = ef.Button(Text = 'OK')
        ok_button.Click += self.OnOKButtonClick
        self.AbortButton = ef.Button(Text = 'Cancel')
        self.AbortButton.Click += self.OnCancelButtonClick

        layout.BeginVertical()
        layout.AddRow(None, ok_button, self.AbortButton, None)
        layout.EndVertical()

        self.Content = layout

    def OnUnderlyingChange(self, sender, e):
        #print('-'*10
        #for x in self.checkBoxList.SelectedValues:
        ss = [str(item) for item in self.checkBoxList.SelectedValues]
        self.bUnder_A = self.sUnder_A in ss
        self.bUnder_B = self.sUnder_B in ss
        self.conduit.Enabled = False
        sc.doc.Views.Redraw()
        self.listbox.Items.Clear()
        self.cs_nested = self.d4_cs[(self.bUnder_A, self.bUnder_B)]
        self.fTols = self.d4_fTols[(self.bUnder_A, self.bUnder_B)]
        self.fillListBox()

    def CreateCheckBoxListRow(self):
        self.checkBoxList = ef.CheckBoxList()
        self.sUnder_A = "Use face A's underlying surface"
        self.checkBoxList.Items.Add(self.sUnder_A)
        self.sUnder_B = "Use face B's underlying surface"
        self.checkBoxList.Items.Add(self.sUnder_B)
        self.checkBoxList.Orientation = ef.Orientation.Vertical
        self.checkBoxList.SelectedValuesChanged += self.OnUnderlyingChange
        return self.checkBoxList

    def OnSelectedIndexChanged(self, sender, e):
        index = self.listbox.SelectedIndex
        if index >= 0:
            self.conduit.Enabled = False
            item = self.listbox.Items[index]
            self.conduit.curves = self.cs_nested[index]
            self.conduit.Enabled = True
            sc.doc.Views.Redraw()

    def CreateListBoxRow(self):
        self.listbox = ef.ListBox()
        #self.m_listbox.Size = ed.Size(200, 100)
        self.listbox.SelectedIndexChanged += self.OnSelectedIndexChanged
        d_row = ef.DynamicRow()
        d_row.Add(self.listbox)
        return d_row


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


    def OnScaleUpClick(self, sender, e):
        self.AdjustScale(1)


    def OnScaleDownClick(self, sender, e):
        self.AdjustScale(-1)


    def AdjustScale(self, direction):
        current_val = self.ParseToFloat(self.textBoxes['fScale'].Text)
        incr_val = self.ParseToFloat(self.textBoxes['fScaleIncrement'].Text)
        
        if current_val is not None and incr_val is not None:
            new_val = current_val + (incr_val * direction)
            if new_val > Rhino.RhinoMath.ZeroTolerance:
                # Format to a maximum of 4 decimal places, stripping trailing zeros
                self.textBoxes['fScale'].Text = "{:g}".format(round(new_val, 4))

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


    def SaveSettings(self):
        key = 'bBothEnds'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = bool(self.radioButtonLists[key].SelectedIndex)

        key = 'fScale'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = float(self.textBoxes[key].Text)

        key = 'fScaleIncrement'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = float(self.textBoxes[key].Text)

        key = 'iCont_Picked'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = bool(self.radioButtonLists[key].SelectedIndex)

        key = 'iCont_Opp'
        sc.sticky[Opts.stickyKeys[key]] = Opts.values[key] = bool(self.radioButtonLists[key].SelectedIndex)

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

        self.SaveSettings()
        self.Close()

        return


        UInt32_Undo = sc.doc.BeginUndoRecord("Surface-Surface Intersect")

        for c in self.cs_nested[i]:
            sc.doc.Objects.AddCurve(c)
        sc.doc.Views.Redraw()

        if not sc.doc.EndUndoRecord(UInt32_Undo):
            print("Warning: EndUndoRecord==False")

        self.Close()


    def OnSaveSettingsButtonClick(self, sender, e):
            self.SaveSettings()
            # Optional: Print a status message so the user knows it worked
            print("Settings saved as default.")


    def OnCancelButtonClick(self, sender, e):
        self.Close()


    def OnFormClosed(self, sender, e):
        self.conduit.Enabled = False
        sc.doc.Views.Redraw()


def getInput_GUI(objref_In):

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

    dialog.DefaultButton = ef.Button()
    dialog.AbortButton = ef.Button()

    rv = Rhino.UI.EtoExtensions.ShowSemiModal(dialog, sc.doc, parent)
    print(rv)

    dialog.conduit.Enabled = False
    sc.doc.Views.Redraw()


def createCurve(nc_In, fScale, iEndToScale=2, iG_T0=2, iG_T1=2, bDebug=False):
    """
    Parameters:
        nc_In: rg.NurbsCurve
        fScale: float of scale factor
        iTEndsToScale: int(0 for T0, 1 for T1, or 3 for both ends nc_In)
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


    if iEndToScale in (0,2):
        # Scale T0 end.


        if iG_T0 in (0,1,2):
            p0 = nc_In.Points[0].Location
            pts_Prime.append(p0)

            if iG_T0 in (1,2):
                p1 = nc_In.Points[1].Location

                xform = rg.Transform.Scale(p0, fScale)

                p1p = rg.Point3d(p1)
                p1p.Transform(xform)
                pts_Prime.append(p1p)

                if iG_T0 == 2:
                    p2 = nc_In.Points[2].Location

                    m = ((p1p - p0).Length/(p1 - p0).Length)**2.0

                    p2p = 2.0*p1p + -p0 + m*(-2.0*p1 + p2 + p0)

                    pts_Prime.append(p2p)



    # Duplicate the points that are not needed by either end scale.
    iCt_PtsNotNeededByT1Scale = nc_In.Points.Count - (iG_T1 + 1)
    for i in range(len(pts_Prime), iCt_PtsNotNeededByT1Scale):
        pts_Prime.append(nc_In.Points[i].Location)


    if iEndToScale in (1,2):
        # Scale T1 end.


        if iG_T1 in (0,1,2):
            pts_New_T1 = [] # To be reversed and extended on pts_Prime.

            p0 = nc_In.Points[nc_In.Points.Count-1].Location
            pts_New_T1.append(p0)

            if iG_T1 in (1,2):
                p1 = nc_In.Points[nc_In.Points.Count-2].Location

                xform = rg.Transform.Scale(p0, fScale)

                p1p = rg.Point3d(p1) # G1 CP location prime.
                p1p.Transform(xform)
                pts_New_T1.append(p1p)

                if iG_T1 == 2:
                    p2 = nc_In.Points[nc_In.Points.Count-3].Location

                    m = ((p1p - p0).Length/(p1 - p0).Length)**2.0

                    p2p = 2.0*p1p + -p0 + m*(-2.0*p1 + p2 + p0)

                    pts_New_T1.append(p2p)

            pts_New_T1.reverse()

            pts_Prime.extend(pts_New_T1)


    #for pt in pts_Prime:
    #    sc.doc.Objects.AddPoint(pt)
    #sc.doc.Views.Redraw(); return


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


def processCurveObject(objref_In, **kwargs):
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

    if fScale == 1:
        print("Scale set to 1, so curve cannot be modified. Script canceled.")
        return

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

    if bGUI:
        Rhino.RhinoApp.SetCommandPromptMessage("Continuing in dialog...")
        getInput_GUI(objref_In)



                #Opts.values['bBothEnds'],
                #Opts.values['fScale'],
                #-1 if Opts.values['iCont_Picked'] == 3 else Opts.values['iCont_Picked'],
                #-1 if Opts.values['iCont_Opp'] == 3 else Opts.values['iCont_Opp'],
                #Opts.values['bDeleteInput'],
                #Opts.values['bEcho'],
                #Opts.values['bDebug'],
                #)

    return

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