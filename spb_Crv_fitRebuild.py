"""
An alternative to _Rebuild for curves, this script can test different degrees and control
point counts to find a rebuild result within a deviation tolerance of the input curve.

TODO:
    Limit distance between consecutive control points?
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
170624-27: Created.
...
190609-11: Now, only curve of a single, specified degree are output.  Modified printed output.
190620: Corrected count in command prompt string.
190627: Bug fix.  Again, will now will wait for more input if curves are preselected.
190629: Added an import.  Bug fixes.
190825: Removed an option.  Added bProcessOtherRat.  Modified an option default value.
        Now, correctly processes PolyCurves with 1 segment.
190831: Modified behavior of an option value entry.
190903: Added some default values to a function.
190923: Added xNurbsCurve_fitByTranslatingControlPts.
191020: Corrected a return value of a function.
191021: Import-related update.
191209: Added code for debugging.
200113: Now, 0 for fDevTol will allow any deviation.  Now, preselected curve will immediately be processed.
200115-19: Added dialog box.  Added option to process PolyCurve segments individually.
        Removed a couple of options.
200307: When degree is not specified, now iterates from degree 3 to 5 before incrementing the number of control points.
200401, 10: Improved handling of bad deviation result for closed curves.
200610: Import-related update.  Purged some of this history.
200611: Bug fix.
200622: Modified some option inputs.
210312: Modified some option default values.  Removed a print(statement used for debugging.
210412: Added filter for PolylineCurves.
220328: Added bPreserveEndG2 option.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

import Eto.Drawing as drawing
import Eto.Forms as forms

import xCurve
import xCurve_deviation
import xNurbsCurve_fitByTranslatingControlPts


class Opts:

    keys = []
    values = {}
    names = {}
    sDialogTexts = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    sDialogTexts[key] = ""
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPreserveEndG1'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    sDialogTexts[key] = "Preserve end tangent dirs."
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPreserveEndG2'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Preserve end curvatures/radii"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bFurtherTranslateCps'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Further translate CPs to minimize CP count  (May be slow.)"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitDegree'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Limit degree"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMatchDegree'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Same as input curve"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDegree'; keys.append(key)
    values[key] = 3
    names[key] = 'CrvDegree'
    sDialogTexts[key] = ""
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitMinCpCt'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Limit min. CP count"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iMinCpCt'; keys.append(key)
    values[key] = 4
    names[key] = key[1:]
    sDialogTexts[key] = "Min. CP count"
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitMaxCpCt'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    sDialogTexts[key] = "Limit max. CP count"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMatchMaxCpCt'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Same as input curve"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iMaxCpCt'; keys.append(key)
    values[key] = 32
    names[key] = key[1:]
    sDialogTexts[key] = "Max. CP count"
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessPolyCrv'; keys.append(key)
    values[key] = True
    names[key] = 'PolyCrv'
    sDialogTexts[key] = "Polycurves"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessPolyCrvSegs'; keys.append(key)
    values[key] = True
    names[key] = 'Process'
    sDialogTexts[key] = (
        None,
        "In whole  ",
        "Each segment individually")
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'WholePolyCrv', 'EachSegment')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessLinear'; keys.append(key)
    values[key] = False
    names[key] = 'Linear'
    sDialogTexts[key] = "Lines"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessPolyline'; keys.append(key)
    values[key] = False
    names[key] = 'PolylineInWhole'
    sDialogTexts[key] = "Polylines in whole"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessArc'; keys.append(key)
    values[key] = False
    names[key] = 'Arc'
    sDialogTexts[key] = "Arcs"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessEllipse'; keys.append(key)
    values[key] = True
    names[key] = 'Ellipse'
    sDialogTexts[key] = "Ellipses"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessOtherRat'; keys.append(key)
    values[key] = True
    names[key] = 'OtherRationalNurbs'
    sDialogTexts[key] = "Other rational NURBS"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessNonRatWithInternalPolyknot'; keys.append(key)
    values[key] = True
    names[key] = 'NurbsWithInternalPolyknots'
    sDialogTexts[key] = "Having internal polyknots"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessNonRatWithSomeFullPolyknot'; keys.append(key)
    values[key] = False
    names[key] = 'Multiplicity'
    sDialogTexts[key] = (
        None,
        "Any multiplicity  ",
        "Some full")
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Any', 'SomeFull')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessUniformNonRat'; keys.append(key)
    values[key] = False
    names[key] = 'Uniform'
    sDialogTexts[key] = "Uniform"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessBezierNonRat'; keys.append(key)
    values[key] = False
    names[key] = 'Bezier'
    sDialogTexts[key] = "Bezier"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bProcessOtherNonRat'; keys.append(key)
    values[key] = True
    names[key] = 'NonrationalNurbs'
    sDialogTexts[key] = "Other"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    sDialogTexts[key] = "Output:", "Add new    ", "Replace input"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    sDialogTexts[key] = "Echo"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    sDialogTexts[key] = "Debug"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)


    @classmethod
    def loadSticky(cls):
        for key in cls.stickyKeys:
            if cls.stickyKeys[key] in sc.sticky:
                if key in cls.riOpts:
                    cls.riOpts[key].CurrentValue = cls.values[key] = sc.sticky[cls.stickyKeys[key]]
                else:
                    cls.values[key] = sc.sticky[cls.stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


class FitRebuildCurveDialog(forms.Dialog[bool]):


    def __init__(self):

        Opts.loadSticky()

        # Initialize dialog box
        self.Title = 'Rebuild Curve'
        self.Padding = drawing.Padding(10)
        self.Resizable = False

        self.divider = forms.Label()
        #self.divider.Text = ''
        self.divider.Height = 2
        self.divider.BackgroundColor = drawing.Colors.Gray


        self.radioButtonList = {}
        self.checkBoxes = {}
        self.labels = {}
        self.textBoxes = {}

        for key in Opts.keys:
            if key[0] == 'b':
                if isinstance(Opts.sDialogTexts[key], str):
                    self.checkBoxes[key] = forms.CheckBox(Text=Opts.sDialogTexts[key])
                    self.checkBoxes[key].Checked = Opts.values[key]
                elif isinstance(Opts.sDialogTexts[key], tuple):
                    self.radioButtonList[key] = forms.RadioButtonList()
                    self.radioButtonList[key].DataStore = Opts.sDialogTexts[key][1:]
                    self.radioButtonList[key].SelectedIndex = int(Opts.values[key])
                    #self.radioButtonList[key].Orientation = forms.Orientation.Horizontal
            elif key[0] == 'f':
                self.labels[key] = forms.Label(Text=Opts.sDialogTexts[key])
                self.textBoxes[key] = forms.TextBox(Text=str(Opts.values[key]))
            elif key[0] == 'i':
                self.labels[key] = forms.Label(Text=Opts.sDialogTexts[key])
                self.textBoxes[key] = forms.TextBox(Text=str(int(Opts.values[key])))
            elif key[0] == 's':
                self.radioButtonList[key] = forms.RadioButtonList()
                self.radioButtonList[key].DataStore = Opts.sDialogTexts[key]
                self.radioButtonList[key].SelectedIndex = Opts.sDialogTexts[key].index(Opts.values[key])



        # Set initial state.
        self.enableTextBox_iDegree()
        self.enableTextBox_iMinCpCt()
        self.enableTextBox_iMaxCpCt()
        self.enableCheckBox_bPreserveEndG2()
        self.enableRadioButtonList_bProcessPolyCrvSegs()
        self.enableRadioButtonList_bProcessNonRatWithOnlyFullyMultiK()
        self.enableCheckBox_bProcessBezier()



        # Events.
        self.textBoxes['fDevTol'].TextChanged += self.onChange_fDevTol

        self.checkBoxes['bLimitDegree'].CheckedChanged += self.onChange_bLimitDegree
        self.checkBoxes['bMatchDegree'].CheckedChanged += self.onChange_bMatchDegree
        self.textBoxes['iDegree'].TextChanged += self.onChange_iDegree

        self.checkBoxes['bLimitMinCpCt'].CheckedChanged += self.onChange_bLimitMinCpCt
        self.textBoxes['iMinCpCt'].TextChanged += self.onChange_iMinCpCt

        self.checkBoxes['bLimitMaxCpCt'].CheckedChanged += self.onChange_bLimitMaxCpCt
        self.checkBoxes['bMatchMaxCpCt'].CheckedChanged += self.onChange_bMatchMaxCpCt
        self.textBoxes['iMaxCpCt'].TextChanged += self.onChange_iMaxCpCt

        self.checkBoxes['bPreserveEndG1'].CheckedChanged += self.onChange_bPreserveEndG1

        self.checkBoxes['bProcessPolyCrv'].CheckedChanged += self.OnChange_bProcessPolyCrv

        self.checkBoxes['bProcessUniformNonRat'].CheckedChanged += self.OnChange_bProcessUniformNonRational
        self.checkBoxes['bProcessNonRatWithInternalPolyknot'].CheckedChanged += self.OnChange_bProcessNonRatWithAnyMultiK

        self.DefaultButton = forms.Button(Text = 'OK')
        self.DefaultButton.Click += self.onOkButtonClick

        self.AbortButton = forms.Button(Text = 'Cancel')
        self.AbortButton.Click += self.onCancelButtonClick



        # Layout.
        layout = forms.DynamicLayout()
        #layout.Spacing = drawing.Size(5, 5)
        #layout.DefaultPadding = 5
        layout.DefaultSpacing = drawing.Size(5, 5)

        layout.AddSeparateRow(
            None,
            "Deviation tolerance",
            self.textBoxes['fDevTol'],
            )

        layout.BeginVertical()
        layout.AddSpace()
        layout.AddRow(
            self.checkBoxes['bLimitDegree'],
            self.checkBoxes['bMatchDegree'],
            self.textBoxes['iDegree'],
            )
        layout.AddSpace()
        layout.AddRow(
            self.checkBoxes['bLimitMinCpCt'],
            None,
            self.textBoxes['iMinCpCt'],
            )
        layout.AddSpace()
        layout.AddRow(
            self.checkBoxes['bLimitMaxCpCt'],
            self.checkBoxes['bMatchMaxCpCt'],
            self.textBoxes['iMaxCpCt'],
            )
        layout.EndVertical()

        layout.AddSpace()
        layout.AddSpace()

        layout.BeginVertical()
        layout.AddRow(
            self.checkBoxes['bPreserveEndG1'],
            self.checkBoxes['bPreserveEndG2'],
            )
        layout.EndVertical()

        layout.AddSpace()
        layout.AddSpace()

        layout.BeginVertical()
        layout.AddRow(self.checkBoxes['bFurtherTranslateCps'])
        layout.EndVertical()

        layout.AddSpace()
        layout.AddSpace()

        layout.AddRow("Process these curves:")
        layout.EndVertical()

        layout.BeginVertical()
        layout.AddSeparateRow(
            "",
            self.checkBoxes['bProcessPolyCrv'],
            self.radioButtonList['bProcessPolyCrvSegs'])
        layout.AddSeparateRow(
            "",
            self.checkBoxes['bProcessLinear'],
            self.checkBoxes['bProcessPolyline'],
            )
        layout.AddSeparateRow(
            "",
            self.checkBoxes['bProcessArc'],
            self.checkBoxes['bProcessEllipse'],
            self.checkBoxes['bProcessOtherRat'],
            )
        layout.AddSeparateRow("", "Non-rational NURBS:")
        layout.AddSeparateRow(
            "   ",
            self.checkBoxes['bProcessNonRatWithInternalPolyknot'],
            self.radioButtonList['bProcessNonRatWithSomeFullPolyknot'],
            )
        layout.AddSeparateRow(
            "   ",
            self.checkBoxes['bProcessUniformNonRat'],
            self.checkBoxes['bProcessBezierNonRat'],
            )
        layout.AddSeparateRow(
            "   ",
            self.checkBoxes['bProcessOtherNonRat'],
            None,
            )
        layout.EndVertical()

        layout.AddSeparateRow("")
        key = 'bReplace'
        layout.AddSeparateRow(Opts.sDialogTexts[key][0], self.radioButtonList[key])
        layout.BeginVertical()
        layout.AddRow("")
        layout.AddRow(
            self.checkBoxes['bEcho'],
            None,
            self.checkBoxes['bDebug'],
            )
        layout.AddRow("")
        layout.EndVertical()


        layout.AddSeparateRow(None, self.DefaultButton, None, self.AbortButton, None)
        #layout.AddRow(self.okButton, self.cancelButton)

        # Set the dialog content
        self.Content = layout



    # Functions that set properties of controls based on the values of other controls.
    # This is used in both __init__ and by individual onChange methods.
    def enableCheckBox_bMatchDegree(self):
        self.checkBoxes['bMatchDegree'].Enabled = (
            self.checkBoxes['bLimitDegree'].Checked)


    def enableTextBox_iDegree(self):
        self.textBoxes['iDegree'].Enabled = (
            self.checkBoxes['bLimitDegree'].Checked and
            not self.checkBoxes['bMatchDegree'].Checked)


    def enableTextBox_iMinCpCt(self):
        self.textBoxes['iMinCpCt'].Enabled = (
            self.checkBoxes['bLimitMinCpCt'].Checked)


    def enableCheckBox_bMatchMaxCpCt(self):
        self.checkBoxes['bMatchMaxCpCt'].Enabled = (
            self.checkBoxes['bLimitMaxCpCt'].Checked)

    def enableTextBox_iMaxCpCt(self):
        self.textBoxes['iMaxCpCt'].Enabled = (
            self.checkBoxes['bLimitMaxCpCt'].Checked and
            not self.checkBoxes['bMatchMaxCpCt'].Checked)


    def enableCheckBox_bPreserveEndG2(self):
        self.checkBoxes['bPreserveEndG2'].Enabled = (
            self.checkBoxes['bPreserveEndG1'].Checked)



    def enableRadioButtonList_bProcessPolyCrvSegs(self):
        self.radioButtonList['bProcessPolyCrvSegs'].Enabled = (
            self.checkBoxes['bProcessPolyCrv'].Checked)


    def enableRadioButtonList_bProcessNonRatWithOnlyFullyMultiK(self):
        self.radioButtonList['bProcessNonRatWithSomeFullPolyknot'].Enabled = (
            self.checkBoxes['bProcessNonRatWithInternalPolyknot'].Checked)


    def enableCheckBox_bProcessBezier(self):
        self.checkBoxes['bProcessBezierNonRat'].Enabled = (
            self.checkBoxes['bProcessUniformNonRat'].Checked)



    def onChange_bLimitCrvDev(self, sender, e):
        self.enableTextBox_fDevTol()

    def onChange_fDevTol(self, sender, e):
        try:
            float(self.textBoxes['fDevTol'].Text)
            self.textBoxes['fDevTol'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['fDevTol'].BackgroundColor = drawing.Colors.Red


    def onChange_bLimitDegree(self, sender, e):
        self.enableCheckBox_bMatchDegree()
        self.enableTextBox_iDegree()

    def onChange_bMatchDegree(self, sender, e):
        self.enableTextBox_iDegree()

    def onChange_iDegree(self, sender, e):
        try:
            iDegree = int(self.textBoxes['iDegree'].Text)
            if iDegree < 3:
                iDegree = 3
            elif iDegree > 5:
                iDegree = 5
            self.textBoxes['iDegree'].Text = str(iDegree)
            self.textBoxes['iDegree'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['iDegree'].BackgroundColor = drawing.Colors.Red


    def onChange_bLimitMinCpCt(self, sender, e):
        self.enableTextBox_iMinCpCt()

    def onChange_iMinCpCt(self, sender, e):
        try:
            iMinCpCt = int(self.textBoxes['iMinCpCt'].Text)
            self.textBoxes['iMinCpCt'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['iMinCpCt'].BackgroundColor = drawing.Colors.Red


    def onChange_bLimitMaxCpCt(self, sender, e):
        self.enableCheckBox_bMatchMaxCpCt()
        self.enableTextBox_iMaxCpCt()

    def onChange_bMatchMaxCpCt(self, sender, e):
        self.enableTextBox_iMaxCpCt()

    def onChange_iMaxCpCt(self, sender, e):
        try:
            iMaxCpCt = int(self.textBoxes['iMaxCpCt'].Text)
            self.textBoxes['iMaxCpCt'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['iMaxCpCt'].BackgroundColor = drawing.Colors.Red


    def onChange_bPreserveEndG1(self, sender, e):
        self.enableCheckBox_bPreserveEndG2()


    def OnChange_bProcessPolyCrv(self, sender, e):
        self.enableRadioButtonList_bProcessPolyCrvSegs()


    def OnChange_bProcessNonRatWithAnyMultiK(self, sender, e):
        self.enableRadioButtonList_bProcessNonRatWithOnlyFullyMultiK()


    def OnChange_bProcessUniformNonRational(self, sender, e):
        self.enableCheckBox_bProcessBezier()




    def onCancelButtonClick(self, sender, e):
        self.Close(False)


    def onOkButtonClick(self, sender, e):
        # Save values.
        for key in Opts.keys:
            if key in self.checkBoxes:
                Opts.values[key] = self.checkBoxes[key].Checked
            elif key in self.radioButtonList:
                if key[0] == 's':
                    Opts.values[key] = (
                        Opts.sDialogTexts[key][self.radioButtonList[key].SelectedIndex]
                    )
                if key[0] == 'b':
                    Opts.values[key] = bool(self.radioButtonList[key].SelectedIndex)
            elif key in self.textBoxes:
                if key[0] == 'f':
                    try:
                        Opts.values[key] = float(self.textBoxes[key].Text)
                    except:
                        s  = "Invalid input for {}.".format(key[1:])
                        s += "  {} will be used instead.".format(Opts.values[key])
                        print(s)
                if key[0] == 'i':
                    try:
                        Opts.values[key] = int(self.textBoxes[key].Text)
                    except:
                        s  = "Invalid input for {}.".format(key[1:])
                        s += "  {} will be used instead.".format(Opts.values[key])
                        print(s)


        # Save sticky.
        for key in Opts.keys:
            if key in Opts.stickyKeys:
                sc.sticky[Opts.stickyKeys[key]] = Opts.values[key]

        self.Close(True)


def getPreselectedCurves():
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    if gObjs_Preselected:
        gCrvs_Preselected = []
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
            if rdRhinoObject.Id in gObjs_Preselected:
                if rdRhinoObject.ObjectType == rd.ObjectType.Curve:
                    gCrvs_Preselected.append(rdRhinoObject.Id)
        if gCrvs_Preselected:
            if Opts.values['bEcho']:
                s  = "{} curves".format(len(gCrvs_Preselected))
                s += " were preselected and will thus be the selection set."
                print(s)
            return tuple(gCrvs_Preselected)


def getInput():
    """Get curves with optional input."""



    def setCurveFilter():
        go_CrvTypes = ri.Custom.GetOption()

        go_CrvTypes.SetCommandPrompt("Curve type(s) to process")

        idxs_Opt_CrvTypes = {}

        while True:
            key = 'bProcessPolyCrv'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessPolyCrvSegs'; idxs_Opt_CrvTypes[key] = (Opts.riAddOpts[key](go_CrvTypes)
                                           if Opts.values['bProcessPolyCrv']
                                           else None)
            key = 'bProcessLinear'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessPolyline'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessArc'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessEllipse'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessOtherRat'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessNonRatWithInternalPolyknot'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessNonRatWithSomeFullPolyknot'; idxs_Opt_CrvTypes[key] = (
                Opts.riAddOpts[key](go_CrvTypes)
                if Opts.values['bProcessNonRatWithInternalPolyknot']
                else None)
            key = 'bProcessUniformNonRat'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            key = 'bProcessBezierNonRat'; idxs_Opt_CrvTypes[key] = (Opts.riAddOpts[key](go_CrvTypes)
                                           if Opts.values['bProcessUniformNonRat']
                                           else None)
            key = 'bProcessOtherNonRat'; idxs_Opt_CrvTypes[key] = Opts.riAddOpts[key](go_CrvTypes)
            #idxs_Opt_CrvTypes['YesToAll'] = go_CrvTypes.AddOption('YesToAll')
            #idxs_Opt_CrvTypes['NoToAll'] = go_CrvTypes.AddOption('NoToAll')


            res = go_CrvTypes.Get()
            if res != ri.GetResult.Option:
                break
            #elif go_CrvTypes.OptionIndex() == idxs_Opt_CrvTypes['YesToAll']:
            #    for key in sCrvFilterOpts:
            #        Opts.riOpts[key].CurrentValue = True
            #elif go_CrvTypes.OptionIndex() == idxs_Opt_CrvTypes['NoToAll']:
            #    for key in sCrvFilterOpts:
            #        Opts.riOpts[key].CurrentValue = False

            Opts.setValues()
            Opts.saveSticky()
            go_CrvTypes.ClearCommandOptions()



    sCrvFilterOpts = [s[8:] for s in Opts.keys if s[:8] == 'bProcess']


    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")
    
    go.GeometryFilter = rd.ObjectType.Curve
    
    #    # Custom geometry filter to only select NurbsCurve wire curves.
    #    def curvesEdgesNotInBlockInstGeomFilter(rdObj, geom, compIdx):
    #        if rdObj.ObjectType == rd.ObjectType.InstanceReference: return False
    #        print(geom
    #        print(compIdx.Index
    #        return isinstance(geom, rg.Curve)
    #        return geom.GetType() == rg.Curve
    #    go.SetCustomGeometryFilter(curvesEdgesNotInBlockInstGeomFilter)    
    
    go.AcceptNumber(True, True)
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    print("\nMaxCpCt: \"0\" will allow any amount.")

    print("Curve types to process: {}".format(
        ', '.join([s for s in sCrvFilterOpts if Opts.values['bProcess'+s]])))


    #bPreselectedObjsChecked = False

    Opts.loadSticky()

    idxs_Opts = {}

    while True:
        idxs_Opts['CrvFilter'] = go.AddOption('CrvTypesToProcess')
        key = 'fDevTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bPreserveEndG1'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bPreserveEndG2'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bFurtherTranslateCps'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bLimitDegree'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bMatchDegree'; idxs_Opts[key] = (Opts.riAddOpts[key](go)[0]
                                           if Opts.values['bLimitDegree']
                                           else None)
        key = 'iDegree'; idxs_Opts[key] = (Opts.riAddOpts[key](go)[0]
                                           if (
                                               Opts.values['bLimitDegree'] and
                                               not Opts.values['bMatchDegree'])
                                           else None)
        key = 'bLimitMinCpCt'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'iMinCpCt'; idxs_Opts[key] = (Opts.riAddOpts[key](go)[0]
                                           if Opts.values['bLimitMinCpCt']
                                           else None)
        key = 'bLimitMaxCpCt'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bMatchMaxCpCt'; idxs_Opts[key] = (Opts.riAddOpts[key](go)[0]
                                           if Opts.values['bLimitMaxCpCt']
                                           else None)
        key = 'iMaxCpCt'; idxs_Opts[key] = (Opts.riAddOpts[key](go)
                                           if (
                                               Opts.values['bLimitMaxCpCt'] and
                                               not Opts.values['bMatchDegree'])
                                           else None)
        key = 'bReplace'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        ## Use bPreselectedObjsChecked so that only objects before the
        ## first call to go.GetMultiple is considered.
        #if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
        #    bPreselectedObjsChecked = True
        #    go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple(
                    [objrefs]
                    +
                    [Opts.values[key] for key in Opts.keys])
        else:
            # An option was selected or a number was entered.
            key = 'fDevTol'
            if go.Option().Index == idxs_Opts['CrvFilter']:
                setCurveFilter()
            if res == ri.GetResult.Number:
                Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 0.0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

            
            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()


def removeNesting(rgCrv0):
    """
    Parameters:
        rgCrv0: rg.PolylineCurve or rg.
    Returns:
        Success: Curves without nesting.
        Fail or N/A: None
    """


    def convertMonoSpanPolylineCrvToLineCrv(crv):
        if isinstance(crv, rg.PolylineCurve):
            if crv.PointCount == 2:
                line = rg.Line(crv.PointAtStart, crv.PointAtEnd)
                return rg.LineCurve(line)



    sType_C0 = rgCrv0.GetType().Name

    if not sType_C0 == 'PolyCurve' or sType_C0 == 'PolylineCurve':
        return None

    rgCrv_WIP = rgCrv0.DuplicateCurve()

    if sType_C0 == 'PolylineCurve':
        rc = convertMonoSpanPolylineCrvToLineCrv(rgCrv_WIP)
        if rc:
            rgCrv_WIP.Dispose()
            return rc

    elif sType_C0 == 'PolyCurve':
        bNestingRemoved = rgCrv_WIP.RemoveNesting()

        if rgCrv_WIP.SegmentCount == 1:
            exploded = rgCrv0.Explode()[0]
            rgCrv_WIP.Dispose()
            rgCrv_WIP = exploded
            rc = convertMonoSpanPolylineCrvToLineCrv(rgCrv_WIP)
            if rc:
                rgCrv_WIP.Dispose()
                return rc

        # Process each segment in polysegment PolyCurve.

        bSuccess, pline = rgCrv_WIP.TryGetPolyline()
        if bSuccess:
            rgCrv_WIP.Dispose()
            return pline.ToPolylineCurve()

        segs_forOut = []
        bNewSeg = False

        for seg in rgCrv_WIP.DuplicateSegments():
            rc = convertMonoSpanPolylineCrvToLineCrv(seg)
            if rc:
                seg.Dispose()
                segs_forOut.append(rc)
                bNewSeg = True
            else:
                segs_forOut.append(rgCrv_WIP)

        if bNewSeg:
            rgCrv_WIP.Dispose()
            joined = rg.Curve.JoinCurves(segs_forOut)
            if len(joined) != 1:
                raise (ValueError ("{} curves resulted from Curve.JoinCurves.".format(
                    len(joined))))
            return joined
        else:
            for seg in segs_forOut: seg.Dispose()
            if bNestingRemoved: return rgCrv_WIP

    rgCrv_WIP.Dispose()


def doesCurvePassTypeFilter(rgCurve0, **kwargs):
    """
    Returns: bool
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bProcessPolyCrv = getOpt('bProcessPolyCrv')
    bProcessLinear = getOpt('bProcessLinear')
    bProcessPolyline = getOpt('bProcessPolyline')
    bProcessArc = getOpt('bProcessArc')
    bProcessEllipse = getOpt('bProcessEllipse')
    bProcessOtherRat = getOpt('bProcessOtherRat')
    bProcessNonRatWithInternalPolyknot = getOpt('bProcessNonRatWithInternalPolyknot')
    bProcessNonRatWithSomeFullPolyknot = getOpt('bProcessNonRatWithSomeFullPolyknot')
    bProcessUniformNonRat = getOpt('bProcessUniformNonRat')
    bProcessBezierNonRat = getOpt('bProcessBezierNonRat')
    bProcessOtherNonRat = getOpt('bProcessOtherNonRat')


    sType_c0 = rgCurve0.GetType().Name

    if sType_c0 == 'PolyCurve':
        if bProcessPolyCrv:
            return True, None
        else:
            return False, "Skipped PolyCurve."

    if sType_c0 == 'LineCurve':
        if bProcessLinear:
            return True, None
        else:
            return False, "Skipped LineCurve."

    if sType_c0 == 'PolylineCurve':
        if bProcessPolyline:
            return True, None
        else:
            return False, "Skipped PolylineCurve."

    if sType_c0 == 'ArcCurve':
        if bProcessArc:
            return True, None
        else:
            return False, "Skipped ArcCurve."

    if rgCurve0.IsLinear(1e-6):
        if bProcessLinear:
            return True, None
        else:
            return False, "Skipped linear {}.".format(sType_c0)
    
    rc = xCurve.getArcCurve(
            rgCrv0=rgCurve0,
            bTolByRatio=False,
            fTolRatio=None,
            fDevTol=1e-9,
            fMinNewCrvLen=100.0*sc.doc.ModelAbsoluteTolerance,
            fMaxRadius=(1e6)*sc.doc.ModelAbsoluteTolerance)
    if rc is not None and rc[0] is not None:
        if bProcessArc:
            return True, None
        else:
            return False, "Skipped arc-shaped {}.".format(sType_c0)

    if sType_c0 == "NurbsCurve":
        rc = xCurve.getEllipticalNurbsCurve(
                rgCrv0=rgCurve0,
                bTolByRatio=False,
                fTolRatio=None,
                fDevTol=1e-9,
                fMinNewCrvLen=100.0*sc.doc.ModelAbsoluteTolerance,
                fMaxRadius=(1e6)*sc.doc.ModelAbsoluteTolerance)
        if rc is not None and rc[0] is not None:
            if bProcessEllipse:
                return True, None
            else:
                return False, "Skipped elliptical-shaped {}.".format(sType_c0)

        if rgCurve0.IsRational:
            if bProcessOtherRat:
                return True, None
            else:
                return False, "Skipped rational non-arc, non-elliptical NurbsCurve."


        if xCurve.Nurbs.hasInternalPolyknots(rgCurve0):
            if bProcessNonRatWithInternalPolyknot:
                if not bProcessNonRatWithSomeFullPolyknot:
                    return True, None
                else:
                    if xCurve.Nurbs.hasSomeFullyMultiplePolyknots(rgCurve0):
                        return True, None
                    else:
                        return False, "Skipped non-rational NurbsCurve containing some internal full polyknots."
            else:
                return False, "Skipped non-rational NurbsCurve containing some internal polyknots"


        if xCurve.Nurbs.isUniform(rgCurve0):
            if bProcessUniformNonRat:
                if rgCurve0.SpanCount == 1:
                    if bProcessBezierNonRat:
                        return True, None
                    else:
                        return False, "Skipped Bezier NurbsCurve."
                else:
                    return True, None
            else:
                return False, "Skipped uniform NurbsCurve."

        if bProcessOtherNonRat:
            return True, None
        else:
            return False, "Skipped other non-rational NurbsCurve."


def rebuildCurve(rgCurve0, fDevTol, iDegree=3, bPreserveEndG1=True, bPreserveEndG2=False, bFurtherTranslateCps=True, iMinCpCt=None, iMaxCpCt=20, bDebug=False):
    """
    Returns tuple of 3 values:
        On success: rg.NurbsCurve, float(Max. deviation), None
        On fail:
            None, float(Last deviation calculated (required deviation), None
            None, None, str(Feedback)
    """
    
    sCrvType = rgCurve0.GetType().Name
    
    if sCrvType == "NurbsCurve":
        nc_In = rgCurve0.Duplicate()
    else:
        nc_In = rgCurve0.ToNurbsCurve()
    
    if not nc_In:
        return None, None, "NurbsCurve could not be constructed from {}.".format(rgCurve0)
    
    if iDegree < 0:
        iDegs = nc_In.Degree,
    elif not iDegree:
        iDegs = 3, 5
    else:
        iDegs = iDegree,
    
    if bDebug: sEval = 'iDegs'; print(sEval + ':', eval(sEval))


    if 1 in iDegs:
        if nc_In.IsClosed:
            return None, None, "Skip closed curve for Degree 1 conversion."


        # Attempt to replace curve with a line.
        # Don't Rebuild a degree 1 curve with more than 2 control points
        # because smooth curves, not polylines, should be the geometry output of this function.
        nc_Out = nc_In.Rebuild(
                pointCount=2,
                degree=1,
                preserveTangents=False)
        #


        if nc_Out is None:
            return None, None, "Could not rebuild curve to a 2-point degree 1."


        rc = xCurve_deviation.isMaxClosestDistBtwn2CrvsWithinTol(
            nc_In,
            nc_Out,
            tolerance=fDevTol)
        if rc is not False:
            dev = rc
            return nc_Out, dev, None


    ct_cp = (iDegs[0] + 1 if not iMinCpCt
             else (iMinCpCt if iMinCpCt >= iDegs[0] + 1 else iDegs[0] + 1))

    if iMaxCpCt > 0:
        if iMaxCpCt < ct_cp:
            return None, None, "Minimum control point count for curve degree already exceeds maximum allowed."
    elif iMaxCpCt == -1: 
        iMaxCpCt = nc_In.Points.Count

    dev = None

    while True:
        if sc.escape_test(False):
            s  = "Script stopped at iteration {}".format(ct_cp)
            s += " with last deviation of {}.".format(dev)
            return None, dev, s
        
        for iDeg in iDegs:

            if ct_cp < iDeg + 1: continue

            if iMaxCpCt == 0 and ct_cp > nc_In.Points.Count:
                if sCrvType == "NurbsCurve" and xCurve.Nurbs.isUniform(nc_In):
                    pass
                elif sCrvType == "PolyCurve":
                    pass
                else:
                    return None, dev, None
            elif iMaxCpCt > 0 and ct_cp > iMaxCpCt:
                return None, dev, None


            if bDebug:
                sEval = 'ct_cp'; print(sEval + ':', eval(sEval),)
                sEval = 'iDeg'; print(sEval + ':', eval(sEval))

            #
            #
            nc_Out = nc_In.Rebuild(
                    pointCount=ct_cp,
                    degree=iDeg,
                    preserveTangents=bPreserveEndG1)
            #
            #


            if bPreserveEndG2:
                bSuccess = nc_Out.SetEndCondition(
                    bSetEnd=False,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                    point=nc_Out.PointAtStart,
                    tangent=nc_In.TangentAtStart,
                    curvature=nc_In.CurvatureAt(nc_In.Domain.T0))
                if not bSuccess:
                    print("SetEndCondition failed.")
                bSuccess = nc_Out.SetEndCondition(
                    bSetEnd=True,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                    point=nc_Out.PointAtEnd,
                    tangent=nc_In.TangentAtEnd,
                    curvature=nc_In.CurvatureAt(nc_In.Domain.T1))
                if not bSuccess:
                    print("SetEndCondition failed.")


            if bFurtherTranslateCps:
                rc = xNurbsCurve_fitByTranslatingControlPts.fitCurve(
                        rgNurbsCrv_toDeform=nc_Out,
                        rgNurbsCrv_forDevComp=rgCurve0,
                        bPreserveEndG1=bPreserveEndG1,
                        bDebug=bDebug)
                if rc is not None:
                    nc_Out.Dispose()
                    nc_Out = rc


            rc = xCurve_deviation.isMaxClosestDistBtwn2CrvsWithinTol(
                nc_In,
                nc_Out,
                tolerance=fDevTol)
            if rc is not False:
                dev = rc
                return nc_Out, dev, None

            # To next degree.
        
        # To next control point count increment.

        ct_cp += 1


def rebuildPolyCurveSegments(rgPolyCrv0, fDevTol, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iDegree = getOpt('iDegree')
    bPreserveEndG1 = getOpt('bPreserveEndG1')
    bPreserveEndG2 = getOpt('bPreserveEndG2')
    bFurtherTranslateCps = getOpt('bFurtherTranslateCps')
    iMinCpCt = getOpt('iMinCpCt')
    iMaxCpCt = getOpt('iMaxCpCt')

    bProcessLinear = getOpt('bProcessLinear')
    bProcessPolyline = getOpt('bProcessPolyline')
    bProcessArc = getOpt('bProcessArc')
    bProcessEllipse = getOpt('bProcessEllipse')
    bProcessOtherRat = getOpt('bProcessOtherRat')
    bProcessNonRatWithInternalPolyknot = getOpt('bProcessNonRatWithInternalPolyknot')
    bProcessNonRatWithSomeFullPolyknot = getOpt('bProcessNonRatWithSomeFullPolyknot')
    bProcessUniformNonRat = getOpt('bProcessUniformNonRat')
    bProcessBezierNonRat = getOpt('bProcessBezierNonRat')
    bProcessOtherNonRat = getOpt('bProcessOtherNonRat')

    bDebug = getOpt('bDebug')


    if not(rgPolyCrv0, rg.PolyCurve): return None, [], [], s

    rc = removeNesting(rgPolyCrv0)

    pc_WIP = rc if rc else rgPolyCrv0.DuplicateCurve()


    segs_forOut = []
    devs_Pass = []
    devs_Fail = []
    sLogs = []

    bSomeRebuilt = False # Could also be determined by comparing sFails and segment counts.

    for iSeg, seg_In in enumerate(pc_WIP.DuplicateSegments()):
        if bDebug: sEval = 'iSeg'; print(sEval+':',eval(sEval))
        if iSeg == 2:
            pass

        bPass, sLog = doesCurvePassTypeFilter(
            rgCurve0=seg_In,
            bProcessPolyCrv=True,
            bProcessLinear=bProcessLinear,
            bProcessPolyline=bProcessPolyline,
            bProcessArc=bProcessArc,
            bProcessEllipse=bProcessEllipse,
            bProcessOtherRat=bProcessOtherRat,
            bProcessNonRatWithInternalPolyknot=bProcessNonRatWithInternalPolyknot,
            bProcessNonRatWithSomeFullPolyknot=bProcessNonRatWithSomeFullPolyknot,
            bProcessUniformNonRat=bProcessUniformNonRat,
            bProcessBezierNonRat=bProcessBezierNonRat,
            bProcessOtherNonRat=bProcessOtherNonRat,
            )
        if not bPass:
            segs_forOut.append(seg_In.Duplicate())
            sLogs.append(sLog)
            continue


        ##
        seg_forOut, dev, sLog = rebuildCurve(
                rgCurve0=seg_In,
                fDevTol=fDevTol,
                iDegree=iDegree,
                bPreserveEndG1=bPreserveEndG1,
                bPreserveEndG2=bPreserveEndG2,
                bFurtherTranslateCps=bFurtherTranslateCps,
                iMinCpCt=iMinCpCt,
                iMaxCpCt=iMaxCpCt,
                bDebug=bDebug,
                )
        ##


        if seg_forOut is None:
            segs_forOut.append(seg_In.Duplicate())
            if sLog is not None:
                sLogs.append(sLog)
            elif dev is not None:
                devs_Fail.append(dev)
        else:
            # Success.
            segs_forOut.append(seg_forOut)
            devs_Pass.append(dev)
            if sLog:
                sLogs.append(sLog)
            else:
                sLogs.append("Replaced segment.")
            bSomeRebuilt = True


    if not bSomeRebuilt:
        return None, devs_Pass, devs_Fail, sLogs

    joined = rg.Curve.JoinCurves(segs_forOut)

    for c in segs_forOut: c.Dispose()

    if len(joined) != 1:
        s = "JoinCurves returned {} curves.  They are discarded".format(len(joined))
        for c in joined: c.Dispose()
        return None, [], [], [s]

    return joined[0], devs_Pass, devs_Fail, sLogs


def processCurves(curvesAndEdges0, **kwargs):
    """
    curvesAndEdges0 = (GUIDs of CurveObjects) or BrepEdges
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    iDegree = getOpt('iDegree')
    bPreserveEndG1 = getOpt('bPreserveEndG1')
    bPreserveEndG2 = getOpt('bPreserveEndG2')
    bFurtherTranslateCps = getOpt('bFurtherTranslateCps')
    iMinCpCt = getOpt('iMinCpCt')
    iMaxCpCt = getOpt('iMaxCpCt')

    bProcessPolyCrv = getOpt('bProcessPolyCrv')
    bProcessPolyCrvSegs = getOpt('bProcessPolyCrvSegs')
    bProcessLinear = getOpt('bProcessLinear')
    bProcessPolyline = getOpt('bProcessPolyline')
    bProcessArc = getOpt('bProcessArc')
    bProcessEllipse = getOpt('bProcessEllipse')
    bProcessOtherRat = getOpt('bProcessOtherRat')
    bProcessNonRatWithInternalPolyknot = getOpt('bProcessNonRatWithInternalPolyknot')
    bProcessNonRatWithSomeFullPolyknot = getOpt('bProcessNonRatWithSomeFullPolyknot')
    bProcessUniformNonRat = getOpt('bProcessUniformNonRat')
    bProcessBezierNonRat = getOpt('bProcessBezierNonRat')
    bProcessOtherNonRat = getOpt('bProcessOtherNonRat')

    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')



    def formatDistance(fDistance):
        if fDistance is None:
            return "(No deviation provided)"
        elif fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-3)):
            return "{:.2e}".format(fDistance)
        else:
            return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


    def stringOutput(rgCrv0, rgCrv_Res, crv_dev):
    
        s = "Original curve is a {}.".format(rgCrv0.GetType().Name)
    
        # Uncomment this only if needed:
        #    if not rgNurbsCrv0:
        #        s += "\nNurbsCurve could not be constructed from {}.".format(curveOrEdge0)
        #        s += "  Input vs. output will not be stated."
        #        return s
    
        if rgCrv_Res:
            sType_Crv_Res = rgCrv_Res.GetType().Name
            
            if sType_Crv_Res == 'NurbsCurve':
                nc0 = rgCrv0.ToNurbsCurve()
            
            if sType_Crv_Res == 'NurbsCurve':
                s += "  Prop:I,O"
                s += "  {}:{},{}".format("Deg", nc0.Degree, rgCrv_Res.Degree)
                s += "  {}:{},{}".format("PtCt", nc0.Points.Count, rgCrv_Res.Points.Count)
                s += "  {}:{},{}".format("IsUniform",
                        str(xCurve.Nurbs.isUniform(nc0))[0],
                        str(xCurve.Nurbs.isUniform(rgCrv_Res))[0])
                s += "  {}:{},{}".format("IsRational",
                        str(nc0.IsRational)[0],
                        str(rgCrv_Res.IsRational)[0])
                s += "  {}:{},{}".format("IsClosed",
                        str(nc0.IsClosed)[0],
                        str(rgCrv_Res.IsClosed)[0])
                if nc0.IsClosed or rgCrv_Res.IsClosed:
                    s += "  {}:{},{}".format("IsPeriodic",
                            str(nc0.IsPeriodic)[0],
                            str(rgCrv_Res.IsPeriodic)[0])
        
            nc0.Dispose()
        
            if crv_dev:
                s += "  Deviation: {}".format(formatDistance(crv_dev))
            else:
                s += "  Curve deviation cannot be calculated!"
    
        else:
            s = "Curve could not be created using entered parameters."

        return s


    rdCs_In = []
    gCrvs0 = []
    for curveOrEdge0 in curvesAndEdges0:
        rdC_In = rs.coercerhinoobject(curveOrEdge0)
        rdCs_In.append(rdC_In)
        gC_In = rdC_In.Id
        if gC_In:
            gCrvs0.append(gC_In)

    rgCrv_Res = None

    gCrvs0_Replaced = []
    gCrvs_Added = []
    devs_all = []
    
    fTols_needed = []
    sFails = []


    idxs_AtTenths = [int(round(0.1*i*len(rdCs_In),0)) for i in range(10)]
    
    for iC, rdC_In in enumerate(rdCs_In):
        if iC in idxs_AtTenths:
            Rhino.RhinoApp.SetCommandPrompt(
                "Processing curve {} ...".format(
                    "" if len(rdCs_In) == 1 else "{} of {} ".format(
                        iC+1, len(rdCs_In))))

        rgCrv0 = rdC_In.Geometry
        #rgCrv0 = rs.coercecurve(curveOrEdge0) # Returns various rg.Curves, including BrepEdge.
        if rgCrv0 is None:
            sFails.append("Geometry for {} not found!".format(rgCrv0))
            if bDebug: print(sLog)
            continue

        bCrvIsEdge = isinstance(rgCrv0, rg.BrepEdge)

        if bCrvIsEdge:
            rgCrv0 = curveOrEdge0.EdgeCurve # This is for curve type checking.

        gCrv0 = None if bCrvIsEdge else rdC_In.Id

        sType_Crv0 = rgCrv0.GetType().Name

        rc = removeNesting(rgCrv0)
        if rc:
            rgCrv0.Dispose()
            rgCrv0 = rc

        bPass, sLog = doesCurvePassTypeFilter(
            rgCurve0=rgCrv0,
            bProcessPolyCrv=bProcessPolyCrv,
            bProcessLinear=bProcessLinear,
            bProcessPolyline=bProcessPolyline,
            bProcessArc=bProcessArc,
            bProcessEllipse=bProcessEllipse,
            bProcessOtherRat=bProcessOtherRat,
            bProcessNonRatWithInternalPolyknot=bProcessNonRatWithInternalPolyknot,
            bProcessNonRatWithSomeFullPolyknot=bProcessNonRatWithSomeFullPolyknot,
            bProcessUniformNonRat=bProcessUniformNonRat,
            bProcessBezierNonRat=bProcessBezierNonRat,
            bProcessOtherNonRat=bProcessOtherNonRat,
            )
        if not bPass:
            sFails.append(sLog)
            if bDebug: print(sLog)
            continue
    
        if bDebug: print("Original curve is a {}    {}.".format(rgCrv0.GetType().Name, gCrv0))


        ##

        if sType_Crv0 == 'PolyCurve' and bProcessPolyCrvSegs:
            rc = rebuildPolyCurveSegments(
                    rgPolyCrv0=rgCrv0,
                    fDevTol=fDevTol,
                    iDegree=iDegree,
                    bPreserveEndG1=bPreserveEndG1,
                    bPreserveEndG2=bPreserveEndG2,
                    bFurtherTranslateCps=bFurtherTranslateCps,
                    iMinCpCt=iMinCpCt,
                    iMaxCpCt=iMaxCpCt,
                    bProcessLinear=bProcessLinear,
                    bProcessPolyline=bProcessPolyline,
                    bProcessArc=bProcessArc,
                    bProcessEllipse=bProcessEllipse,
                    bProcessOtherRat=bProcessOtherRat,
                    bProcessNonRatWithInternalPolyknot=bProcessNonRatWithInternalPolyknot,
                    bProcessNonRatWithSomeFullPolyknot=bProcessNonRatWithSomeFullPolyknot,
                    bProcessUniformNonRat=bProcessUniformNonRat,
                    bProcessBezierNonRat=bProcessBezierNonRat,
                    bProcessOtherNonRat=bProcessOtherNonRat,
                    bDebug=bDebug,
                    )
            rgCrv_Res, devs_Pass, devs_Fail, sLogs = rc
            sFails.extend(sLogs)
            if rgCrv_Res is None:
                if bDebug: print(sLog)
                continue
            else:
                s  = "{} segments replaced.".format(len(devs_Pass))
                s += "  {} segments not replaced.".format(len(devs_Fail))
                if devs_Fail:
                    s += "  Last deviations ranged [{}, {}].".format(
                        formatDistance(min(devs_Fail)),
                        formatDistance(max(devs_Fail)),
                        )
                devs_all.extend(devs_Pass)
                dev = max(devs_all) if devs_all else None
                fTols_needed.extend(devs_Fail)
            sSummary = None
        else:
            rgCrv_Res, dev, sLog = rebuildCurve(
                rgCurve0=rgCrv0,
                fDevTol=fDevTol,
                iDegree=iDegree,
                bPreserveEndG1=bPreserveEndG1,
                bPreserveEndG2=bPreserveEndG2,
                bFurtherTranslateCps=bFurtherTranslateCps,
                iMinCpCt=iMinCpCt,
                iMaxCpCt=iMaxCpCt,
                bDebug=bDebug,
                )

            if bDebug or len(curvesAndEdges0) == 1:
                sSummary = stringOutput(rgCrv0, rgCrv_Res, dev)

            if rgCrv_Res is None:
                if sLog is not None:
                    sFails.append(sLog)
                    if bDebug: print(sLog)
                elif dev is not None:
                    fTols_needed.append(dev)
                    if bDebug:
                        s  = "Curve could not be Rebuilt within {}.".format(formatDistance(fDevTol))
                        s += "  Last deviation was {}.".format(formatDistance(dev))
                        print(s)
                continue

            # Success.
            devs_all.append(dev)


        if bReplace and not bCrvIsEdge:
            bReplaced = sc.doc.Objects.Replace(objectId=gCrv0, curve=rgCrv_Res)
            if not bReplaced:
                s = "Curve could not be replaced."
                rgCrv_Res.Dispose()
                sFails.append(s)
                if bDebug: print(s)
            else:
                if bDebug and sSummary: print(sSummary)
                rgCrv_Res.Dispose()
                gCrvs0_Replaced.append(gCrv0)
        else:
            gCrv_Res = sc.doc.Objects.AddCurve(rgCrv_Res)
            if gCrv_Res != Guid.Empty:
                if bReplace and bCrvsAdded and not bCrvIsEdge:
                    sc.doc.Objects.Delete(curveOrEdge0, True)
                if bDebug and sSummary: print(sSummary)
                rgCrv_Res.Dispose()
                gCrvs_Added.append(gCrv_Res)
            else:
                s = "Curve could not be added."
                rgCrv_Res.Dispose()
                sFails.append(s)
                if bDebug: print(s)


    if not bEcho: return gCrvs0_Replaced, gCrvs_Added


    # Output summary.

    if len(curvesAndEdges0) == 1:
        if sType_Crv0 == 'PolyCurve' and bProcessPolyCrvSegs:
            s = "Out of {} polycurve segments:".format(rgCrv0.SegmentCount)
            for sFail in set(sFails):
                s += "\n[{}] {}".format(sFails.count(sFail), sFail)
            if fTols_needed:
                s += "\n[{}] Segments not rebuilt within tolerance.".format(len(fTols_needed))
                s += "  Tolerances needed: [{0:.{2}e} through {1:.{2}e}]".format(
                        min(fTols_needed), max(fTols_needed), 2)
            if gCrvs0_Replaced:
                s += "\nPolyCurve was replaced."
                s += "  Segment deviations: [{0:.{2}e} through {1:.{2}e}]".format(
                        min(devs_all), max(devs_all), 2)
            if gCrvs_Added:
                s += "\nPolyCurve was added."
                s += "  Segment deviations: [{0:.{2}e} through {1:.{2}e}]".format(
                        min(devs_all), max(devs_all), 2)
        else:
            if rgCrv_Res is None and sFails: print(sFails[0])
            elif bEcho and sSummary: print(sSummary)
            s = ""
            if fTols_needed:
                s += "Curve could not be rebuilt within tolerance."
                s += "  Tolerance needed: {:.{}e}".format(fTols_needed[0], 2)
            if gCrvs0_Replaced:
                s += "Curve was replaced"
                s += " at a deviation of {:.{}e}.".format(devs_all[0], 2)
            if gCrvs_Added:
                s += "Curve was added"
                s += " at a deviation of {:.{}e}.".format(devs_all[0], 2)
    else:
        s = "Out of {} total curves:".format(len(curvesAndEdges0))
        for sFail in set(sFails):
            s += "\n[{}] {}".format(sFails.count(sFail), sFail)
        if fTols_needed:
            s += "\n[{}] Curve not rebuilt within tolerance.".format(len(fTols_needed))
            s += "  Tolerances needed: [{0:.{2}e} through {1:.{2}e}]".format(
                    min(fTols_needed), max(fTols_needed), 2)
        if gCrvs0_Replaced:
            s += "\n{} curves were replaced.".format(len(gCrvs0_Replaced))
            s += "  Deviations: [{0:.{2}e} through {1:.{2}e}]".format(
                    min(devs_all), max(devs_all), 2)
        if gCrvs_Added:
            s += "\n{} curves were added.".format(len(gCrvs_Added))
            s += "  Deviations: [{0:.{2}e} through {1:.{2}e}]".format(
                    min(devs_all), max(devs_all), 2)
        if not (gCrvs0_Replaced or gCrvs_Added):
            s += "\nNo curves were added or replaced."

    print(s)


    return gCrvs0_Replaced, gCrvs_Added


def main():

    gCrvs0_Preselected = getPreselectedCurves()

    dialog = FitRebuildCurveDialog()
    if not dialog.ShowModal(Rhino.UI.RhinoEtoApp.MainWindow): return


    if gCrvs0_Preselected:
        objrefs = None
    else:
        rc = getInput()
        if rc is None: return
        objrefs = rc[0]

    
    if Opts.values['bDebug']:
        pass
    else:
        sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll()
    
    gCrvs0_Replaced, gNurbsCrvs1 = None, None

    if Opts.values['bLimitDegree']:
        if Opts.values['bMatchDegree']:
            iDegree = -1
        else:
            iDegree = Opts.values['iDegree']
    else:
        iDegree = None

    if Opts.values['bLimitMinCpCt']:
        iMinCpCt = Opts.values['iMinCpCt']
    else:
        iMinCpCt = None

    if Opts.values['bLimitMaxCpCt']:
        if Opts.values['bMatchMaxCpCt']:
            iMaxCpCt = -1
        else:
            iMaxCpCt = Opts.values['iMaxCpCt']
    else:
        iMaxCpCt = None

    rc = processCurves(
        curvesAndEdges0=gCrvs0_Preselected if gCrvs0_Preselected else objrefs,
        fDevTol=Opts.values['fDevTol'],
        iDegree=iDegree,
        iMinCpCt=iMinCpCt,
        iMaxCpCt=iMaxCpCt,
        )
    if rc is not None:
        gCrvs0_Replaced, gNurbsCrvs1 = rc
    
    if gCrvs0_Preselected:
        [sc.doc.Objects.Select(objectId=_) for _ in gCrvs0_Preselected]
    else:
        if gCrvs0_Replaced:
            [sc.doc.Objects.Select(objectId=_) for _ in gCrvs0_Replaced]
        if gNurbsCrvs1:
            [sc.doc.Objects.Select(objectId=_) for _ in gNurbsCrvs1]

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()