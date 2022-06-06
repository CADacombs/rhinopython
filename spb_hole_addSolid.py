"""
# -*- coding: utf-8 -*- # Needed to have '°' in this script.
Creates breps that can be _BooleanDifference(d) from other breps to make holes.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
150808: First useful working version developed.
150812: ListView selection now accepts keyboard input.
150816-17: Added more error checking.
150916: Nominal in Hole Sub Category is now selected by default.
181219-20: Refactored.
        On repeats, dialog box is now populated with last values.
        Changed drill point selection to radio button with 180, 118, and text box.
211102: Fixed bug that was causing a "RhinoDotNetCrash".
        Changed default drill point selection from 180 to 118.
220604: Now allows flipping of hole direction during point picking.  Refactored.
220605: Added preview.  Refactored.

TODO: 
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System.Windows import Forms

from System.Windows.Forms import (
    Form,
    Button, Panel, RadioButton, GroupBox, CheckBox,
    Label, TextBox, FormBorderStyle,
    ListView, View, ColumnHeader, ListViewItem,
    )

from System.Drawing import Size, Point, ContentAlignment, Color, SystemColors

import math


sHoleTypes_Main = "Simple", "Counterbore", "Countersink"

dialogHoleData = {}

dialogHoleData["Simple"] = (
    ("Nominal", (
        ("1/64", 0.015625),
        ("72 drill wire gauge", 0.02500),
        ("1/32", 0.03125),
        ("60 drill wire gauge)", 0.04000),
        ("3/64", 0.046875),
        ("1/16", 0.06250),
        ("5/64", 0.078125),
        ("3/32", 0.09375),
        ("7/64", 0.109375),
        ("1/8", 0.12500),
        ("3/16", 0.18750),
        ("1/4", 0.25000),
        ("5/16", 0.31250),
        ("3/8", 0.37500),
        ("13/32", 0.40625),
        ("7/16", 0.43750),
        ("1/2", 0.50000),
        ("9/16", 0.56250),
        ("5/8", 0.62500),
        ("3/4", 0.75000),
        ("7/8", 0.87500),
        ("1", 1.00000),
        )
     ),
    ("Nominal + 1/64", (
        ("#4", 0.125),
        ("#5", 0.140625),
        ("#6", 0.15625),
        ("#8", 0.171875),
        ("#10", 0.203125),
        ("1/4", 0.265625),
        ("5/16", 0.328125),
        ("3/8", 0.390625),
        ("7/16", 0.453125),
        ("1/2", 0.515625),
        ("9/16", 0.578125),
        ("5/8", 0.640625),
        ("3/4", 0.765625),
        ("7/8", 0.890625),
        ("1", 1.015625),
        )
    ),
    ("Nominal + 1/32", (
        ("#4", 0.140625),
        ("#5", 0.15625),
        ("#6", 0.171875),
        ("#8", 0.1875),
        ("#10", 0.21875),
        ("1/4", 0.28125),
        ("5/16", 0.34375),
        ("3/8", 0.40625),
        ("7/16", 0.46875),
        ("1/2", 0.53125),
        ("9/16", 0.59375),
        ("5/8", 0.65625),
        ("3/4", 0.78125),
        ("7/8", 0.90625),
        ("1", 1.03125),
        )
    ),
    ("Coarse Tap Drill", (
        ("#4-40", 0.0890),
        ("#5-40", 0.1015),
        ("#6-32", 0.10650),
        ("#8-32", 0.13600),
        ("#10-24", 0.14950),
        ("#12-24", 0.17700),
        ("1/4-20", 0.20100),
        ("5/16-18", 0.25700),
        ("3/8-16", 0.31250),
        ("7/16-14", 0.3680),
        ("1/2-13", 0.42188),
        ("5/8-11", 0.5312),
        ("3/4-10", 0.6562),
        ("1-8", 0.8750),
        )
    ),
    ("Fine Tap Drill", (
        ("#4-48", 0.0935),
        ("#5-44", 0.1040),
        ("#6-40", 0.11300),
        ("#8-36", 0.13600),
        ("#10-32", 0.15900),
        ("#12-28", 0.18200),
        ("1/4-28", 0.21300),
        ("5/16-24", 0.27200),
        ("3/8-24", 0.33200),
        ("7/16-20", 0.3906),
        ("1/2-20", 0.4531),
        ("5/8-18", 0.5781),
        ("3/4-16", 0.6875),
        ("1-12", 0.9375),
        )
    ),
    ("Press Fit Dowel", (
        ("1/32", 0.0310),
        ("1/16", 0.0620),
        ("3/32", 0.0932),
        ("1/8", 0.1245),
        ("5/32", 0.1557),
        ("3/16", 0.1870),
        ("1/4", 0.2495),
        ("5/16", 0.3120),
        ("3/8", 0.3745),
        ("7/16", 0.4370),
        ("1/2", 0.4995),
        ("5/8", 0.6245),
        ("3/4", 0.7495),
        ("7/8", 0.8745),
        ("1", 0.9995),
        )
    ),
    ("Slip Fit Dowel", (
        ("1/32", 0.0320),
        ("1/16", 0.0635),
        ("3/32", 0.0955),
        ("1/8", 0.1265),
        ("3/16", 0.1895),
        ("1/4", 0.2520),
        ("3/8", 0.3770),
        )
    ),
    ("Close Fit Screw", (
        ("#6", 0.14400),
        ("#8", 0.16950),
        ("#10", 0.19600),
        ("#12", 0.22100),
        ("1/4", 0.25700),
        ("5/16", 0.32300),
        ("3/8", 0.38600),
        ("7/16", 0.45310),
        ("1/2", 0.51563),
        ("9/16", 0.56250),
        )
    ),
    ("Free Fit Screw", (
        ("#6", 0.14950),
        ("#8", 0.17700),
        ("#10", 0.20100),
        ("#12", 0.22800),
        ("1/4", 0.26600),
        ("5/16", 0.33200),
        ("3/8", 0.39700),
        ("7/16", 0.46780),
        ("1/2", 0.53125),
        )
     )
    )

dialogHoleData["Counterbore"] = (
    ("Nominal", (
        ("#4", 0.1120, 0.1860),
        ("#5", 0.1250, 0.2050),
        ("#6", 0.1380, 0.2280),
        ("#8", 0.1640, 0.2720),
        ("#10", 0.1900, 0.3140),
        ("1/4", 0.2500, 0.3820),
        ("5/16", 0.3125, 0.4750),
        ("3/8", 0.3750, 0.5720),
        ("7/16", 0.4375, 0.6630),
        ("1/2", 0.5000, 0.7570),
        ("5/8", 0.6250, 0.9450),
        ("3/4", 0.7500, 1.1330),
        ("7/8", 0.8750, 1.3220),
        ("1", 1.0000, 1.5100),
        )
    ),
    ("Nominal + 1/64", (
        ("#4", 0.1120 + 0.0156, 0.1860 + 0.0156),
        ("#5", 0.1250 + 0.0156, 0.2050 + 0.0156),
        ("#6", 0.1380 + 0.0156, 0.2280 + 0.0156),
        ("#8", 0.1640 + 0.0156, 0.2720 + 0.0156),
        ("#10", 0.1900 + 0.0156, 0.3140 + 0.0156),
        ("1/4", 1 / 4 + 1 / 64, 0.3820 + 0.0156),
        ("5/16", 5 / 16 + 1 / 64, 0.4750 + 0.0156),
        ("3/8", 3 / 8 + 1 / 64, 0.5720 + 0.0156),
        ("7/16", 7 / 16 + 1 / 64, 0.6630 + 0.0156),
        ("1/2", 1 / 2 + 1 / 64, 0.7570 + 0.0156),
        ("5/8", 5 / 8 + 1 / 64, 0.9450 + 0.0156),
        ("3/4", 3 / 4 + 1 / 64, 1.1330 + 0.0156),
        ("7/8", 7 / 8 + 1 / 64, 1.3220 + 0.0156),
        ("1", 1.0000 + 1 / 64, 1.5100 + 0.0156),
        )
    ),
    ("Nominal + 1/32", (
    ("#4", 0.1120 + 0.0312, 0.1860 + 0.0312),
    ("#5", 0.1250 + 0.0312, 0.2050 + 0.0312),
    ("#6", 0.1380 + 0.0312, 0.2280 + 0.0312),
    ("#8", 0.1640 + 0.0312, 0.2720 + 0.0312),
    ("#10", 0.1900 + 0.0312, 0.3140 + 0.0312),
    ("1/4", 1 / 4 + 1 / 32, 0.3820 + 0.0312),
    ("5/16", 5 / 16 + 1 / 32, 0.4750 + 0.0312),
    ("3/8", 3 / 8 + 1 / 32, 0.5720 + 0.0312),
    ("7/16", 7 / 16 + 1 / 32, 0.6630 + 0.0312),
    ("1/2", 1 / 2 + 1 / 32, 0.7570 + 0.0312),
    ("5/8", 5 / 8 + 1 / 32, 0.9450 + 0.0312),
    ("3/4", 3 / 4 + 1 / 32, 1.1330 + 0.0312),
    ("7/8", 7 / 8 + 1 / 32, 1.3220 + 0.0312),
    ("1", 1.0000 + 1 / 32, 1.5100 + 0.0312),
    )
     )
    )

dialogHoleData["Countersink"] = (
    (
        "Nominal",
        (
            ("#8", 0.164, 0.35900),
            ("#10", 0.190, 0.41100),
            ("1/4", 1 / 4, 0.53125),
            ("5/16", 5 / 16, 0.65625),
            ("3/8", 3 / 8, 0.78125),
            ("7/16", 7 / 16, 0.84375),
            ("1/2", 1 / 2, 0.93750),
            ("5/8", 5 / 8, 1.18750),
            ("3/4", 3 / 4, 1.43750),
            ("7/8", 7 / 8, 1.68750),
            ("1", 1.0000, 1.9375),
        )
    ),
    (
        "Nominal + 1/64",
        (
            ("#8", 3 / 16, 0.35900),
            ("#10", 13 / 64, 0.41100),
            ("1/4", 1 / 4 + 1 / 64, 0.53125),
            ("5/16", 5 / 16 + 1 / 64, 0.65625),
            ("3/8", 3 / 8 + 1 / 64, 0.78125),
            ("7/16", 7 / 16 + 1 / 64, 0.84375),
            ("1/2", 1 / 2 + 1 / 64, 0.93750),
            ("5/8", 5 / 8 + 1 / 64, 1.18750),
            ("3/4", 3 / 4 + 1 / 64, 1.43750),
            ("7/8", 7 / 8 + 1 / 64, 1.68750),
            ("1", 1.0000 + 1 / 64, 1.9375),
        )
    ),
    ("Nominal + 1/32", (
        ("#8", 13 / 64, 0.35900),
        ("#10", 7 / 32, 0.41100),
        ("1/4", 1 / 4 + 1 / 32, 0.53125),
        ("5/16", 5 / 16 + 1 / 32, 0.65625),
        ("3/8", 3 / 8 + 1 / 32, 0.78125),
        ("7/16", 7 / 16 + 1 / 32, 0.84375),
        ("1/2", 1 / 2 + 1 / 32, 0.93750),
        ("5/8", 5 / 8 + 1 / 32, 1.18750),
        ("3/4", 3 / 4 + 1 / 32, 1.43750),
        ("7/8", 7 / 8 + 1 / 32, 1.68750),
        ("1", 1.0000 + 1 / 32, 1.9375),
        )
     )
    )


class Opts:
    
    keys = []
    values = {}
    stickyKeys = {}
    
    key = 'iHoleType'
    keys.append(key)
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iHoleSubCat'
    keys.append(key)
    values[key] = 3
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fMainDia'
    keys.append(key)
    values[key] = 0.201
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fNegZBasicLen'
    keys.append(key)
    values[key] = 1.0
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fDrillPtAngle'
    keys.append(key)
    values[key] = 118.0
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bAddPt'
    keys.append(key)
    values[key] = True
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fCBoreDia'
    keys.append(key)
    values[key] = 1.0
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fCSinkDia'
    keys.append(key)
    values[key] = 1.0
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fCSinkMaxDia'
    keys.append(key)
    values[key] = 1.0 * Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Inches,
            sc.doc.ModelUnitSystem)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fPosZLen'
    keys.append(key)
    values[key] = 1.0
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bEcho'
    keys.append(key)
    values[key] = None
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    keys.append(key)
    values[key] = None
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if sc.sticky.has_key(stickyKeys[key]):
            values[key] = sc.sticky[stickyKeys[key]]
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
        pass


class HoleForm(Form):
    
    def __init__(self):
        
        self.holeTypes_Main = sHoleTypes_Main
        self.idxHoleType_Main_Current = None
        self.sHoleTypes_Sub_Current = None

        self.sHoleTypes_Sub = (
            [dialogHoleData["Simple"][x][0] for x in range(len(dialogHoleData["Simple"]))],
            [dialogHoleData["Counterbore"][x][0] for x in range(len(dialogHoleData["Counterbore"]))],
            [dialogHoleData["Countersink"][x][0] for x in range(len(dialogHoleData["Countersink"]))],
            )

        self.Name = self.Text = 'Fastener Hole Generator'
        self.FormBorderStyle = FormBorderStyle.FixedSingle
        self.Width = 760
        self.Height = 600
        self.CenterToScreen()
        
        self.col1X = 10
        textBox_Width = 100
        label_Width = 198
        col2X = self.col1X + label_Width + 10
        col1YNow = 10
        rowSpc = 10
        
        
        # Create Main Hole Type group box and radio buttons.
        
        self._groupBox_HTMain = GroupBox()
        self._groupBox_HTMain.Text = "Hole Type"
        self._groupBox_HTMain.Size = Size(120, len(self.holeTypes_Main) * 20 + 30)
        self._groupBox_HTMain.Location = Point(self.col1X, col1YNow)
        self.Controls.Add(self._groupBox_HTMain)
        col1YNow = self._groupBox_HTMain.Bottom
        
        self._radioButton_HoleType = []
        
        for rbIndex in range(len(self.holeTypes_Main)):
            self._radioButton_HoleType.append(RadioButton())
            self._radioButton_HoleType[rbIndex].Text = self.holeTypes_Main[rbIndex]
            self._radioButton_HoleType[rbIndex].Parent = self._groupBox_HTMain
            self._radioButton_HoleType[rbIndex].Location = Point(self.col1X, 20 * (rbIndex+1))
            self._radioButton_HoleType[rbIndex].CheckedChanged += self.HoleTypeMainChanged
        
        
        # Create Hole Sub Category group box.  Radio buttons are added in HoleTypeMainChanged.
        
        self._groupBox_HoleSubCat = GroupBox()
        self._groupBox_HoleSubCat.Text = "Hole Sub Category"
        self._groupBox_HoleSubCat.Location = Point(self._groupBox_HTMain.Right + 40, self._groupBox_HTMain.Location.Y)
        self.Controls.Add(self._groupBox_HoleSubCat)
        
        
        #Create ListView
        
        self.column1 = ColumnHeader()
        self.column1.Text = 'Nominal'
        self.column1.Width = 180
        self.column2 = ColumnHeader()
        self.column2.Text = 'Main Hole Dia.'
        self.column2.Width = 100
        self.column3 = ColumnHeader()
        self.column3.Width = 100
        
        #self.SuspendLayout()
        
        self._listView = ListView()
        self._listView.Parent = self
        self._listView.Width = self.column1.Width + self.column2.Width + self.column3.Width + 4
        self._listView.Height = 488
        self._listView.Location = Point(self.Width - self._listView.Width - 15, 10)
        self._listView.FullRowSelect = True
        self._listView.MultiSelect = False
        self._listView.GridLines = True
        self._listView.AllowColumnReorder = True
        self._listView.Columns.AddRange((self.column1, self.column2, self.column3))
        
        self._listView.View = View.Details
        
        #self._listView.Click += self.HoleSelectedInLV
        self._listView.ItemSelectionChanged += self.HoleSelectedInLV
        
        #self.ResumeLayout()
        
        
        #Create Main Hole Diameter text box and label
        
        col1YNow = 220 + rowSpc
        
        self._label_fMainDia = Label()
        self._label_fMainDia.Text = "Main Hole Diameter"
        self._label_fMainDia.Location = Point(self.col1X, col1YNow)
        self._label_fMainDia.Width = label_Width
        self.Controls.Add(self._label_fMainDia)
        
        self._textBox_fMainDia = TextBox()
        self._textBox_fMainDia.Text = str(Opts.values['fMainDia']) #"0.201"
        self._textBox_fMainDia.Location = Point(col2X, col1YNow)
        self._textBox_fMainDia.Width = textBox_Width
        self.Controls.Add(self._textBox_fMainDia)
        
        col1YNow = self._label_fMainDia.Bottom
        
        
        #Create Length text box and label
        
        col1YNow += rowSpc
        
        self._label_fBasicNegZLen = Label()
        self._label_fBasicNegZLen.Text = "Simple hole depth (-Z) (w/o tip)"
        self._label_fBasicNegZLen.Location = Point(self.col1X, col1YNow)
        self._label_fBasicNegZLen.Width = label_Width
        self.Controls.Add(self._label_fBasicNegZLen)
        
        self._textBoxLength_fBasicNegZLen = TextBox()
        self._textBoxLength_fBasicNegZLen.Text = str(Opts.values['fNegZBasicLen']) #"1.0"
        self._textBoxLength_fBasicNegZLen.Location = Point(col2X, col1YNow)
        self._textBoxLength_fBasicNegZLen.Width = textBox_Width
        self.Controls.Add(self._textBoxLength_fBasicNegZLen)
        
        col1YNow = self._textBoxLength_fBasicNegZLen.Bottom
        
        
        col1YNow += rowSpc
        
        # Create Drill Point panel, label, 3 radio buttons, and text box.
        
        radioButtons_drillPointAngle = []
        
        self._panel_DrillPointAngle = Panel()
        self._label_DrillPt = Label()
        self._radioButton_118Deg = RadioButton()
        self._radioButton_180Deg = RadioButton()
        self._radioButton_CustomDeg = RadioButton()
        self._textBox_fDrillPtAngle = TextBox()
        self.Controls.Add(self._panel_DrillPointAngle)
        
        self._panel_DrillPointAngle.Controls.Add(self._label_DrillPt)
        self._panel_DrillPointAngle.Controls.Add(self._radioButton_118Deg)
        self._panel_DrillPointAngle.Controls.Add(self._radioButton_180Deg)
        self._panel_DrillPointAngle.Controls.Add(self._radioButton_CustomDeg)
        self._panel_DrillPointAngle.Controls.Add(self._textBox_fDrillPtAngle)
        self._panel_DrillPointAngle.Location = Point(self.col1X, col1YNow)
        self._panel_DrillPointAngle.Size = Size(self._textBoxLength_fBasicNegZLen.Right-self.col1X, 20)
        self._panel_DrillPointAngle.Name = "_panel_DrillPointAngle"
        
        self._label_DrillPt.Text = "Drill point angle"
        self._label_DrillPt.Location = Point(0, 2)
        self._label_DrillPt.Width = 6*int(len(self._label_DrillPt.Text))
        rightPoint = self._label_DrillPt.Right
        
        self._radioButton_118Deg.Text = "118°"
        self._radioButton_118Deg.Location = Point(rightPoint+10, 0)
        self._radioButton_118Deg.Width = 30 + 6*int(len(self._radioButton_118Deg.Text))
        rightPoint = self._radioButton_118Deg.Right
        
        self._radioButton_180Deg.Text = "180°"
        self._radioButton_180Deg.Location = Point(rightPoint+20, 0)
        self._radioButton_180Deg.Width = 30 + 6*int(len(self._radioButton_180Deg.Text))
        rightPoint = self._radioButton_180Deg.Right
        
        self._radioButton_CustomDeg.Text = ""
        self._radioButton_CustomDeg.Location = Point(rightPoint+10, 0)
        self._radioButton_CustomDeg.Width = 18 + 6*int(len(self._radioButton_CustomDeg.Text))
        self._radioButton_CustomDeg.UseVisualStyleBackColor = True
        rightPoint = self._radioButton_CustomDeg.Right
        
        self._textBox_fDrillPtAngle.Text = str(Opts.values['fDrillPtAngle'])
        self._textBox_fDrillPtAngle.Location = Point(rightPoint, 0)
        self._textBox_fDrillPtAngle.Width = self._panel_DrillPointAngle.Width - self._textBox_fDrillPtAngle.Left
        
        if Opts.values['fDrillPtAngle'] == 118.0:
            self._radioButton_118Deg.Checked = True
        elif Opts.values['fDrillPtAngle'] == 180.0:
            self._radioButton_180Deg.Checked = True
        else:
            self._radioButton_CustomDeg.Checked = True
        
        col1YNow = self._panel_DrillPointAngle.Bottom
        
        
        #Create Add Point checkbox
        
        col1YNow += rowSpc
        
        self._checkBox_bAddPt = CheckBox()
        self._checkBox_bAddPt.CheckAlign = ContentAlignment.MiddleRight
        self._checkBox_bAddPt.TextAlign = ContentAlignment.MiddleLeft
        self._checkBox_bAddPt.Parent = self
        self._checkBox_bAddPt.Text = "Add base point at each hole"
        self._checkBox_bAddPt.Location = Point(self.col1X, col1YNow)
        self._checkBox_bAddPt.Width = label_Width + 28
        self._checkBox_bAddPt.Checked = Opts.values['bAddPt'] #True
        
        col1YNow = self._checkBox_bAddPt.Bottom
        
        
        #Create Counterbore text box and label
        
        col1YNow += rowSpc
        
        self._label_CBore = Label()
        self._label_CBore.Text = "Counterbore diameter"
        self._label_CBore.Location = Point(self.col1X, col1YNow)
        self._label_CBore.Width = label_Width
        self.Controls.Add(self._label_CBore)
        
        self._textBox_fCBoreDia = TextBox()
        if Opts.values['fCBoreDia'] is not None:
            self._textBox_fCBoreDia.Text = str(Opts.values['fCBoreDia'])
        self._textBox_fCBoreDia.Location = Point(col2X, col1YNow)
        self._textBox_fCBoreDia.Width = textBox_Width
        self.Controls.Add(self._textBox_fCBoreDia)
        
        col1YNow = self._textBox_fCBoreDia.Bottom
        
        
        #Create Countersink text box and label
        
        col1YNow += rowSpc
        
        self._label_CSink = Label()
        self._label_CSink.Text = "Countersink diameter"
        self._label_CSink.Location = Point(self.col1X, col1YNow)
        self._label_CSink.Width = label_Width
        self.Controls.Add(self._label_CSink)
        
        self._textBox_fCSinkDia = TextBox()
        if Opts.values['fCSinkDia'] is not None:
            self._textBox_fCSinkDia.Text = str(Opts.values['fCSinkDia'])
        self._textBox_fCSinkDia.Location = Point(col2X, col1YNow)
        self._textBox_fCSinkDia.Width = textBox_Width
        self.Controls.Add(self._textBox_fCSinkDia)
        
        col1YNow = self._textBox_fCSinkDia.Bottom
        
        
        #Create Countersink body diameter text box and label
        
        col1YNow += rowSpc
        
        self._label_CSinkBody = Label()
        self._label_CSinkBody.Text = "Countersink maximum diameter"
        self._label_CSinkBody.Location = Point(self.col1X, col1YNow)
        self._label_CSinkBody.Width = label_Width
        self.Controls.Add(self._label_CSinkBody)
        
        self._textBox_fCSinkMaxDia = TextBox()
        if Opts.values['fCSinkMaxDia'] is not None:
            self._textBox_fCSinkMaxDia.Text = str(Opts.values['fCSinkMaxDia']) #"1.0"
        self._textBox_fCSinkMaxDia.Location = Point(col2X, col1YNow)
        self._textBox_fCSinkMaxDia.Width = textBox_Width
        self.Controls.Add(self._textBox_fCSinkMaxDia)
        
        col1YNow = self._textBox_fCSinkMaxDia.Bottom
        
        
        #Create Negative Extrusion Length label and text box
        
        col1YNow += rowSpc
        
        self._label_fPosZLen = Label()
        self._label_fPosZLen.Text = "Length above pick point (+Z)"
        self._label_fPosZLen.Location = Point(self.col1X, col1YNow)
        self._label_fPosZLen.Width = label_Width
        self.Controls.Add(self._label_fPosZLen)
        
        self._textBox_fPosZLen = TextBox()
        self._textBox_fPosZLen.Text = str(Opts.values['fPosZLen']) #"0.1"
        self._textBox_fPosZLen.Location = Point(col2X, col1YNow)
        self.Controls.Add(self._textBox_fPosZLen)
        
        
        #Create Create and Cancel buttons
        
        iButtonPadding = 10
        
        self._panel_Buttons = Panel()
        
        self._button_Create = Button()
        self._button_Create.Parent = self._panel_Buttons
        self._button_Create.Height += iButtonPadding * 2
        
        self._panel_Buttons.Height = self._button_Create.Height + 8
        self._panel_Buttons.Width = self._listView.Left - 30
        self._panel_Buttons.Location = Point(
                self.col1X,
                self.Height - self._panel_Buttons.Height - 40
        )
        
        self._button_Create.Text = "Create"
        self._button_Create.Width += iButtonPadding * 2
        self._button_Create.Location = Point(
                (self._panel_Buttons.Width - self._button_Create.Width) / 2,
                0
        )
        self._button_Create.Click += self.CreateClick
        
        self._button_Cancel = Button()
        self._button_Cancel.Parent = self._panel_Buttons
        self._button_Cancel.Text = "Cancel"
        self._button_Cancel.Width = self._button_Create.Width
        self._button_Cancel.Height = self._button_Create.Height
        self._button_Cancel.Location = Point(
                self._panel_Buttons.Width - (self._button_Cancel.Width + iButtonPadding),
                0
        )
        
        self.Controls.Add(self._panel_Buttons)
        
        self.AcceptButton = self._button_Create
        self.CancelButton = self._button_Cancel
        
        self._radioButton_HoleType[Opts.values['iHoleType']].Checked = True
        self._radioButton_HoleSubCat[Opts.values['iHoleSubCat']].Checked = True
    
    
    def HoleTypeMainChanged(self, sender, event):
        
        if sender.Checked:
            
            self._groupBox_HoleSubCat.Controls.Clear()
            self._listView.Items.Clear()
            
            self.idxHoleType_Main_Current = self.holeTypes_Main.index(sender.Text)
            
            if "Counterbore" in sender.Text:
                self.column3.Text = "C'Bore Dia."
                self._textBox_fCBoreDia.Enabled = True
                self._textBox_fCSinkDia.Enabled = False
                self._textBox_fCSinkMaxDia.Enabled = False
                self._radioButton_118Deg.Checked = True
            elif "Countersink" in sender.Text:
                self.column3.Text = "C'Sink Dia."
                self._textBox_fCSinkDia.Enabled = True
                self._textBox_fCSinkMaxDia.Enabled = True
                self._textBox_fCBoreDia.Enabled = False
                self._radioButton_118Deg.Checked = True
            else:
                self.column3.Text = ""
                self._textBox_fCBoreDia.Enabled = False
                self._textBox_fCSinkDia.Enabled = False
                self._textBox_fCSinkMaxDia.Enabled = False
            
            self.sHoleTypes_Sub_Current = self.sHoleTypes_Sub[self.idxHoleType_Main_Current]
            
            self._radioButton_HoleSubCat = []
            
            iCt_HoleSubs = len(self.sHoleTypes_Sub_Current)
            
            self._groupBox_HoleSubCat.Size = Size(140, iCt_HoleSubs * 20 + 30)
            
            for rbIndex in range(iCt_HoleSubs):
                self._radioButton_HoleSubCat.append(RadioButton())
                self._radioButton_HoleSubCat[rbIndex].Text = self.sHoleTypes_Sub_Current[rbIndex]
                self._radioButton_HoleSubCat[rbIndex].Parent = self._groupBox_HoleSubCat
                self._radioButton_HoleSubCat[rbIndex].Location = Point(self.col1X, 20*(rbIndex+1))
                self._radioButton_HoleSubCat[rbIndex].Width = 128
                self._radioButton_HoleSubCat[rbIndex].CheckedChanged += self.holeType_Sub_Changed
    
    
    def holeType_Sub_Changed(self, sender, event):
        
        if sender.Checked:
            
            self._listView.Items.Clear()
            
            tplHoleTypeData = dialogHoleData[self.holeTypes_Main[self.idxHoleType_Main_Current]][self.sHoleTypes_Sub_Current.index(sender.Text)][1]
            
            for holeDataLine in tplHoleTypeData:
                item = ListViewItem()
                item.Text = holeDataLine[0]
                if len(holeDataLine) > 1: item.SubItems.Add(str(holeDataLine[1]))
                if len(holeDataLine) == 3:
                    item.SubItems.Add(str(holeDataLine[2]))
                self._listView.Items.Add(item)
    
    
    def HoleSelectedInLV(self, sender, event):
        if sender.SelectedItems:
            self._textBox_fMainDia.Text = sender.SelectedItems[0].SubItems[1].Text
            #self._textBox_fPosZLen.Text = str(float(sender.SelectedItems[0].SubItems[1].Text) * 2)
            if self.column3.Text == "C'Bore Dia.":
                # Adjust minimum c'bore length per its diameter.
                self._textBox_fCBoreDia.Text = sender.SelectedItems[0].SubItems[2].Text
                if float(self._textBox_fPosZLen.Text) < 2.0 * float(self._textBox_fMainDia.Text):
                    self._textBox_fPosZLen.Text = str(2.0 * float(self._textBox_fMainDia.Text))
            elif self.column3.Text == "C'Sink Dia.":
                self._textBox_fCSinkDia.Text = sender.SelectedItems[0].SubItems[2].Text
    
    
    def CreateClick(self, sender, event):
        
        #Check for valid input
        bErrorsFound = False
        
        s = "Main Hole Diameter"
        try:
            self.fMainDia = float(self._textBox_fMainDia.Text)
            self._textBox_fMainDia.BackColor = SystemColors.Window
            if not (self.fMainDia > 0.0):
                print(s + " is not > 0.")
                self.fMainDia = None
                self._textBox_fMainDia.BackColor = Color.Red
                bErrorsFound = True
        except:
            print(s + " is not a float value.")
            self.fMainDia = None
            self._textBox_fMainDia.BackColor = Color.Red
            bErrorsFound = True
        
        s = "Basic Length"
        try:
            fNegZBasicLen = float(self._textBoxLength_fBasicNegZLen.Text)
            self._textBoxLength_fBasicNegZLen.BackColor = SystemColors.Window
            if not (fNegZBasicLen > 0.0):
                print(s + " is not > 0.")
                fNegZBasicLen = None
                self._textBoxLength_fBasicNegZLen.BackColor = Color.Red
                bErrorsFound = True
        except:
            print(s + " is not a float value.")
            fNegZBasicLen = None
            self._textBoxLength_fBasicNegZLen.BackColor = Color.Red
            bErrorsFound = True
        
        s = "Drill Point Angle"
        try:
            fDrillPtAngle = float(self._textBox_fDrillPtAngle.Text)
            self._textBox_fDrillPtAngle.BackColor = SystemColors.Window
            if not (fDrillPtAngle > 0.0 and fDrillPtAngle <= 180.0):
                print(s + " is not > 0° and <= 180°.")
                fDrillPtAngle = None
                self._textBox_fDrillPtAngle.BackColor = Color.Red
                bErrorsFound = True
        except:
            print(s + " is not a float value.")
            fDrillPtAngle = None
            self._textBox_fDrillPtAngle.BackColor = Color.Red
            bErrorsFound = True
        
        bAddPt = self._checkBox_bAddPt.Checked
        
        if self._textBox_fCBoreDia.Enabled:
            s = "Counterbore Diameter"
            try:
                fCBoreDia = float(self._textBox_fCBoreDia.Text)
                self._textBox_fCBoreDia.BackColor = SystemColors.Window
                if not (fCBoreDia > 0.0):
                    print(s + " is not > 0.")
                    fCBoreDia = None
                    self._textBox_fCBoreDia.BackColor = Color.Red
                    bErrorsFound = True
                if self.fMainDia:
                    if fCBoreDia <= self.fMainDia:
                        print(s + " is <= Main Hole Diameter.")
                        fCBoreDia = None
                        self._textBox_fCBoreDia.BackColor = Color.Red
                        bErrorsFound = True
            except:
                print(s + " is not a float value.")
                fCBoreDia = None
                self._textBox_fCBoreDia.BackColor = Color.Red
                bErrorsFound = True
        else: fCBoreDia = None
        
        if self._textBox_fCSinkDia.Enabled:
            s = "Countersink Diameter"
            try:
                fCSinkDia = float(self._textBox_fCSinkDia.Text)
                self._textBox_fCSinkDia.BackColor = SystemColors.Window
                if not (fCSinkDia > 0.0):
                    print(s + " is not > 0.")
                    fCSinkDia = None
                    self._textBox_fCSinkDia.BackColor = Color.Red
                    bErrorsFound = True
                if self.fMainDia:
                    if fCSinkDia <= self.fMainDia:
                        print(s + " is <= Main Hole Diameter.")
                        fCSinkDia = None
                        self._textBox_fCSinkDia.BackColor = Color.Red
                        bErrorsFound = True
            except:
                print(s + " is not a float value.")
                fCSinkDia = None
                self._textBox_fCSinkDia.BackColor = Color.Red
                bErrorsFound = True
        else: fCSinkDia = None
        
        if self._textBox_fCSinkMaxDia.Enabled:
            if self._textBox_fCSinkMaxDia.Text:
                s = "Countersink Body Diameter"
                try:
                    fCSinkMaxDia = float(self._textBox_fCSinkMaxDia.Text)
                    self._textBox_fCSinkMaxDia.BackColor = SystemColors.Window
                    if not (fCSinkMaxDia > 0.0):
                        print(s + " is not > 0.")
                        fCSinkMaxDia = None
                        self._textBox_fCSinkMaxDia.BackColor = Color.Red
                        bErrorsFound = True
                    if fCSinkDia:
                        if fCSinkMaxDia <= fCSinkDia:
                            print(s + " is <= Countersink Diameter.")
                            fCSinkMaxDia = None
                            self._textBox_fCSinkMaxDia.BackColor = Color.Red
                            bErrorsFound = True
                except:
                    print(s + " is not a float value.")
                    fCSinkMaxDia = None
                    self._textBox_fCSinkMaxDia.BackColor = Color.Red
                    bErrorsFound = True
            else: fCSinkMaxDia = None
        else: fCSinkMaxDia = None
        
        s = "Negative Extrusion Length"
        try:
            fPosZLen = float(self._textBox_fPosZLen.Text)
            self._textBox_fPosZLen.BackColor = SystemColors.Window
            if not (fPosZLen >= 0.0):
                print(s + " is not >= 0.")
                fPosZLen = None
                self._textBox_fPosZLen.BackColor = Color.Red
                bErrorsFound = True
        except:
            print(s + " is not a float value.")
            fPosZLen = None
            self._textBox_fPosZLen.BackColor = Color.Red
            bErrorsFound = True
        
        if bErrorsFound: return False
        
        #self.Hide()
        self.Close()
        
        Opts.values['iHoleType'] = [_.Checked for _ in self._radioButton_HoleType].index(True)
        Opts.values['iHoleSubCat'] = [_.Checked for _ in self._radioButton_HoleSubCat].index(True)
        Opts.values['fMainDia'] = float(self._textBox_fMainDia.Text)
        Opts.values['fNegZBasicLen'] = float(self._textBoxLength_fBasicNegZLen.Text)
        if self._radioButton_180Deg.Checked:
            Opts.values['fDrillPtAngle'] = 180.0
        elif self._radioButton_118Deg.Checked:
            Opts.values['fDrillPtAngle'] = 118.0
        else:
            Opts.values['fDrillPtAngle'] = float(self._textBox_fDrillPtAngle.Text)
        Opts.values['bAddPt'] = self._checkBox_bAddPt.Checked
        Opts.values['fCBoreDia'] = float(self._textBox_fCBoreDia.Text)
        Opts.values['fCSinkDia'] = float(self._textBox_fCSinkDia.Text)
        Opts.values['fCSinkMaxDia'] = float(self._textBox_fCSinkMaxDia.Text)
        Opts.values['fPosZLen'] = float(self._textBox_fPosZLen.Text)
        
        Opts.saveSticky()
        
        self.DialogResult = Forms.DialogResult.OK


def createPointsAlongProfile(**kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iHoleType = getOpt('iHoleType')
    fMainDia = getOpt('fMainDia')
    fNegZBasicLen = getOpt('fNegZBasicLen')
    fDrillPtAngle = getOpt('fDrillPtAngle')
    bAddPt = getOpt('bAddPt')
    fCBoreDia = getOpt('fCBoreDia')
    fCSinkDia = getOpt('fCSinkDia')
    fCSinkMaxDia = getOpt('fCSinkMaxDia')
    fPosZLen = getOpt('fPosZLen')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    pts = []
    
    fMainRad = fMainDia / 2.0

    Pt = rg.Point3d
    
    # Accumulate endpoints without drill point.
    if iHoleType == 0:
        # Simple hole.
        pts.append(Pt(fMainRad, 0.0, fPosZLen))
        pts.append(Pt(fMainRad, 0.0, -fNegZBasicLen))
    elif iHoleType == 1:
        # Counterbore with simple hole.
        fCBoreRad = fCBoreDia / 2.0
        pts.append(Pt(fCBoreRad, 0.0, fPosZLen))
        pts.append(Pt(fCBoreRad, 0.0, 0.0))
        pts.append(Pt(fMainRad, 0.0, 0.0))
        pts.append(Pt(fMainRad, 0.0, -fNegZBasicLen))
    elif iHoleType == 2:
        # Countersink with maximum body diameter with simple hole.
        fCSinkRad = fCSinkDia / 2.0
        fCSinkBodyRad = fCSinkMaxDia / 2.0
        fTheoreticalCSinkMaxRad = fCSinkRad + fPosZLen * math.tan(math.radians(41.0))
        if fTheoreticalCSinkMaxRad > fCSinkBodyRad:
            pts.append(Pt(fCSinkBodyRad, 0.0, fPosZLen))
            pts.append(Pt(
                fCSinkBodyRad,
                0.0,
                (fCSinkBodyRad - fCSinkRad) * math.tan(math.radians(49.0)) ) )
        else:
            pts.append((fTheoreticalCSinkMaxRad, 0.0, fPosZLen))
        pts.append(Pt(
            fMainRad,
            0.0,
            -(fCSinkRad - fMainRad) * math.tan(math.radians(49.0))))
        pts.append(Pt(fMainRad, 0.0, -fNegZBasicLen))
    else:
        raise ValueError("Error!  iHoleType value, {}, is not valid.".format(iHoleType))

    if fDrillPtAngle < 180.0:
        #Add endpoint of drill angle
        pts.append(Pt(
            0.0,
            0.0,
            -fNegZBasicLen - fMainRad / math.tan(math.radians(fDrillPtAngle / 2))
            ))

    return pts


def createBrepToCopy(**kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iHoleType = getOpt('iHoleType')
    fMainDia = getOpt('fMainDia')
    fNegZBasicLen = getOpt('fNegZBasicLen')
    fDrillPtAngle = getOpt('fDrillPtAngle')
    bAddPt = getOpt('bAddPt')
    fCBoreDia = getOpt('fCBoreDia')
    fCSinkDia = getOpt('fCSinkDia')
    fCSinkMaxDia = getOpt('fCSinkMaxDia')
    fPosZLen = getOpt('fPosZLen')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    pline = rg.Polyline(
        createPointsAlongProfile(
            iHoleType=iHoleType,
            fMainDia=fMainDia,
            fNegZBasicLen=fNegZBasicLen,
            fDrillPtAngle=fDrillPtAngle,
            fCBoreDia=fCBoreDia,
            fCSinkDia=fCSinkDia,
            fCSinkMaxDia=fCSinkMaxDia,
            fPosZLen=fPosZLen,
            )
        )
    if pline is None: return


    srf = rg.RevSurface.Create(
        revolutePolyline=pline,
        axisOfRevolution=rg.Line(rg.Point3d.Origin, span=rg.Vector3d.ZAxis))

    return rg.Brep.CreateFromRevSurface(srf, capStart=True, capEnd=True)


def dynamicallyAddObjects(brep_ToCopy, bAddPt):
    """
    Get Point3d with optional input.
    """

    gp = ri.Custom.GetPoint()

    #gp.SetCommandPrompt("Left click to flip direction")


    def calculateTransform():
        pt_Picked = gp.Point()

        plane = rg.Plane(sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane())
        plane.Origin = pt_Picked

        if bFlipped:
            xform = rg.Transform.Rotation(
                angleRadians=math.pi,
                rotationAxis=rg.Vector3d.YAxis,
                rotationCenter=rg.Point3d.Origin)
            xform = rg.Transform.ChangeBasis(plane, rg.Plane.WorldXY) * xform
        else:
            xform = rg.Transform.ChangeBasis(plane, rg.Plane.WorldXY)

        return xform


    def createTransformedBrep(xform):
        xform = calculateTransform()
        brep_Out = brep_ToCopy.DuplicateBrep()
        brep_Out.Transform(xform)

        return brep_Out


    def brepDraw(sender, args):
        brep_In = args.Source.Tag
        if not brep_In: return

        xform = calculateTransform()
        if not xform: return

        brep_Out = createTransformedBrep(xform)
        if not brep_Out: return

        args.Display.DrawBrepWires(brep_Out, color=Color.AliceBlue)


    gp.DynamicDraw += brepDraw
    gp.Tag = brep_ToCopy


    idxs_Opt = {}

    idxs_Opt['Flip'] = gp.AddOption('FlipDrillDir')
    bFlipped = False

    while True:
        res = gp.Get()

        if res == ri.GetResult.Cancel:
            gp.Dispose()
            return

        if res == ri.GetResult.Point:

            pt_Picked = gp.Point()

            xform = calculateTransform()
            if not xform: return

            brep_Out = createTransformedBrep(xform)
            if not brep_Out: return


            gHoleCopy = sc.doc.Objects.AddBrep(brep_Out)

            if bAddPt: sc.doc.Objects.AddPoint(pt_Picked)

        if gp.OptionIndex() == idxs_Opt['Flip']:
            print("Flipped drill direction.")
            bFlipped = not bFlipped


def main():
    
    #fMainDia = fNegZBasicLen = fDrillPtAngle = bAddPt = fCBoreDia = fCSinkDia = fCSinkMaxDia = fPosZLen = None
    
    form = HoleForm()
    
    if not form.ShowDialog() == Forms.DialogResult.OK:
        return
    
    iHoleType = Opts.values['iHoleType']
    fMainDia = Opts.values['fMainDia']
    fNegZBasicLen = Opts.values['fNegZBasicLen']
    fDrillPtAngle = Opts.values['fDrillPtAngle']
    bAddPt = Opts.values['bAddPt']
    fCBoreDia = Opts.values['fCBoreDia']
    fCSinkDia = Opts.values['fCSinkDia']
    fCSinkMaxDia = Opts.values['fCSinkMaxDia']
    fPosZLen = Opts.values['fPosZLen']

    sc.doc.Objects.UnselectAll()

    #sc.doc.Views.RedrawEnabled = False

    brep_ToCopy = createBrepToCopy(
        iHoleType=iHoleType,
        fMainDia=fMainDia,
        fNegZBasicLen=fNegZBasicLen,
        fDrillPtAngle=fDrillPtAngle,
        bAddPt=bAddPt,
        fCBoreDia=fCBoreDia,
        fCSinkDia=fCSinkDia,
        fCSinkMaxDia=fCSinkMaxDia,
        fPosZLen=fPosZLen,
        )
    if brep_ToCopy is None: return

    #bFlipped = False


    return dynamicallyAddObjects(brep_ToCopy, bAddPt)


if __name__ == '__main__': main()