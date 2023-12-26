"""
This script is an alternative to _Intersect, but for 2 BrepFaces only.

Unique features to _Intersect:
    The underlying surface can be intersected for either or both faces.
    Results for various tolerances are calculated and the user can chose one from a dialog box.

Points are not output.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
210914-21: Created.
210929: Bug fixed and streamlined getInput.
220820: Bug fix.
231224: Added options to make the curves without the dialog.  In that case, only the input tolerance is used.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import Eto.Drawing as drawing
import Eto.Forms as forms


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bUseDialog'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fMinTol'; keys.append(key)
    values[key] = 1e-5
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fMinCrvLength'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bUnderlying1st'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUnderlying2nd'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key in ('fTol', 'fMinTol'):
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Error: {} was passed to setValue.".format(key))
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Select 2 BrepFaces with options.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 brep faces")

    go.GeometryFilter = go.GeometryFilter.Surface

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bUseDialog')
        Opts.names['fTol'] = 'MaxTol' if Opts.values['bUseDialog'] else 'Tol'
        addOption('fTol')
        if Opts.values['bUseDialog']:
            addOption('fMinTol')
        addOption('fMinCrvLength')
        addOption('bUnderlying1st')
        addOption('bUnderlying2nd')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()

            return objrefs

        if res == ri.GetResult.Number:
            key = 'fTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


class DrawCurvesConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        #self.color = sc.doc.Layers.CurrentLayer.Color
        self.color = Rhino.ApplicationSettings.AppearanceSettings.FeedbackColor
        self.curves = None

    #def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
    #    if len(self.curves) > 0:
    #        self.bbox = self.brep.GetBoundingBox(accurate=False)
    #        calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):

        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

        for c in self.curves:
            drawEventArgs.Display.DrawCurve(
                curve=c,
                color=self.color,
                thickness=crv_thk)


class EtoModelessForm(forms.Form):
    """
    Started with
    https://github.com/mcneel/rhino-developer-samples/blob/3179a8386a64602ee670cc832c77c561d1b0944b/rhinopython/SampleEtoModelessForm.py
    """


    def __init__(self):
        self.m_selecting = False
        self.Initialize()
        self.CreateFormControls()
        self.conduit = DrawCurvesConduit()

    def Initialize(self):
        self.Title = 'Surface-Surface Intersection'
        self.Padding = drawing.Padding(5)
        self.Resizable = False
        self.Maximizable = False
        self.Minimizable = False
        self.ShowInTaskbar = False
        self.MinimumSize = drawing.Size(230, 150)
        self.Closed += self.OnFormClosed

    def addCurveSetListItem(self, cs, fTol):
        item = forms.ListItem()
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
        layout = forms.DynamicLayout()
        layout.Padding = drawing.Padding(10)
        layout.Spacing = drawing.Size(5, 5)

        layout.Rows.Add(self.CreateCheckBoxListRow())

        layout.AddRow(None)

        iPrec = sc.doc.ModelDistanceDisplayPrecision + 1

        sLabel = "   {:<14}{:<12}{:<10}{:<8}{:<6}".format(
            'Tol', 'Length', 'Crvs', 'Spans', 'CPs')

        layout.Rows.Add(forms.Label(Text=sLabel))
        layout.Rows.Add(self.CreateListBoxRow())

        ok_button = forms.Button(Text = 'OK')
        ok_button.Click += self.OnOKButtonClick
        self.AbortButton = forms.Button(Text = 'Cancel')
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
        self.checkBoxList = forms.CheckBoxList()
        self.sUnder_A = "Use face A's underlying surface"
        self.checkBoxList.Items.Add(self.sUnder_A)
        self.sUnder_B = "Use face B's underlying surface"
        self.checkBoxList.Items.Add(self.sUnder_B)
        self.checkBoxList.Orientation = forms.Orientation.Vertical
        self.checkBoxList.SelectedValuesChanged += self.OnUnderlyingChange
        return self.checkBoxList

    def OnSelectedIndexChanged(self, sender, e):
        index = self.listbox.SelectedIndex
        if index >= 0:
            self.conduit.Enabled = False
            self.m_selecting = True
            item = self.listbox.Items[index]
            self.conduit.curves = self.cs_nested[index]
            self.conduit.Enabled = True
            sc.doc.Views.Redraw()
            self.m_selecting = False

    def CreateListBoxRow(self):
        self.listbox = forms.ListBox()
        #self.m_listbox.Size = drawing.Size(200, 100)
        self.listbox.SelectedIndexChanged += self.OnSelectedIndexChanged
        d_row = forms.DynamicRow()
        d_row.Add(self.listbox)
        return d_row

    def OnOKButtonClick(self, sender, e):
        i = self.listbox.SelectedIndex
        if i < 0:
            self.Close()
            print("No curve set was selected.")
            return

        UInt32_Undo = sc.doc.BeginUndoRecord("Surface-Surface Intersect")

        for c in self.cs_nested[i]:
            sc.doc.Objects.AddCurve(c)
        sc.doc.Views.Redraw()

        if not sc.doc.EndUndoRecord(UInt32_Undo):
            print("Warning: EndUndoRecord==False")

        self.Close()

    def OnCancelButtonClick(self, sender, e):
        self.Close()

    def OnFormClosed(self, sender, e):
        self.conduit.Enabled = False
        sc.doc.Views.Redraw()


def _getIntersectionMethod(rgObj_A, rgObj_B):
    if isinstance(rgObj_A, rg.Brep):
        if isinstance(rgObj_B, rg.Brep):
            return rg.Intersect.Intersection.BrepBrep
        elif isinstance(rgObj_B, rg.Surface):
            return rg.Intersect.Intersection.BrepSurface
    elif isinstance(rgObj_A, rg.Surface):
        if isinstance(rgObj_B, rg.Brep):
            return lambda b, a, tol: rg.Intersect.Intersection.BrepSurface(a, b, tol)
        elif isinstance(rgObj_B, rg.Surface):
            return rg.Intersect.Intersection.SurfaceSurface
    raise ValueError("Intersection of {} and {} is not supported by this script.".format(
        rgObj_A.GetType().Name, rgObj_B.GetType().Name))


def createCurves_NoLoops(rgObj_A, rgObj_B, fTol, fMinCrvLength):
    intersectionMethod = _getIntersectionMethod(rgObj_A, rgObj_B)

    bSuccess, rgCrvs_SSX_WIP, rgPts_SSX = intersectionMethod(
        rgObj_A,
        rgObj_B,
        fTol)

    # At least in RhinoCommon of 7 SR10,
    # even when no curves or points result, the result of
    # Intersection.BrepBrep and
    # Intersection.BrepSurface
    # is / (may still be) True, but
    # Intersection.SurfaceSurface is False.
    if not bSuccess:
        return

    if len(rgPts_SSX) > 0:
        print("{} points in intersection of tolerance {}.".format(
            rgPts_SSX,
            fTol))

    if len(rgCrvs_SSX_WIP) == 0:
        return

    rgCrvs_SSX = []
    for c in rgCrvs_SSX_WIP:
        length = c.GetLength()
        if length > fMinCrvLength:
            rgCrvs_SSX.append(c)

    return rgCrvs_SSX


def _areCurvesAlreadyInList(cs, dss, fTol, bDebug):

    # epsilon = max(0.1*fTol, 1e-6)

    dss_Passing = dss

    # Compare curves counts.
    def getRefsWithSameCounts(es, fss):
        out = []
        for fs in fss:
            if len(es) == len(fs):
                out.append(fs)
        return out

    dss_Passing = getRefsWithSameCounts(cs, dss_Passing)
    if not dss_Passing:
        if bDebug: print("Curve count is new.")
        return False

    if bDebug: print("Curve counts is already in list.")

    # Compare end points.
    def getRefsWithSameEndPts(es, fss):
        out = []
        for fs in fss:
            for e, f in zip(es, fs):
                if e.PointAtStart.DistanceTo(f.PointAtStart) > fTol:
                    break # to next fs.
                if e.PointAtEnd.DistanceTo(f.PointAtEnd) > fTol:
                    break # to next fs.
            else:
                # All ends match within fTol.
                out.append(fs)
        return out

    dss_Passing = getRefsWithSameEndPts(cs, dss_Passing)
    if not dss_Passing:
        if bDebug: print("Curve end point(s) is/are new.")
        return False


    if bDebug: print("Curve(s) with the same end point match are already in list.")

    # Compare point counts.
    def getRefsWithSameCP_cts(es, fss):
        out = []
        for fs in fss:
            for e, f in zip(es, fs):
                if not isinstance(e, rg.NurbsCurve):
                    if bDebug: print("{} found.".format(e.GetType().Name))
                    e = e.ToNurbsCurve()
                if not isinstance(f, rg.NurbsCurve):
                    if bDebug: print("{} found.".format(f.GetType().Name))
                    f = f.ToNurbsCurve()
                if bDebug: print(e.Points.Count, f.Points.Count)
                if e.Points.Count != f.Points.Count:
                    break # to next fs.
            else:
                # All point counts match.
                out.append(fs)
        return out

    dss_Passing = getRefsWithSameCP_cts(cs, dss_Passing)
    if not dss_Passing:
        if bDebug: print("Control point count is new.")
        return False


    if bDebug: print("Control point count is already in list.")


    if bDebug: sEval = "len(dss_Passing)"; print("{}: {}".format(sEval, eval(sEval)))

    return True


    #Compare lengths
    #for ds in dss:
        #fLength_Cs = 0.0
        #for c in cs:
        #    f = c.GetLength()
        #    fLength_Cs += f

        #fLength_Ds = 0.0
        #for d in ds:
        #    f = d.GetLength()
        #    fLength_Ds += f

        #if abs(fLength_Cs - fLength_Ds) <= fTol:
        #    return True


        # Compare NurbsCurve forms by EpsilonEquals.
        #for c, d in zip(cs, ds):
        #    if c.GetType() != d.GetType():
        #        break

        #    nsC = c.ToNurbsCurve()
        #    nsD = d.ToNurbsCurve()

        #    epsEqual = nsC.EpsilonEquals(nsD, epsilon)

        #    if epsEqual and (nsC.SpanCount != nsD.SpanCount):
        #        raise Exception("NCs pass EpsilonEquals but have different span counts!")

        #    if epsEqual and (nsC.Points.Count != nsD.Points.Count):
        #        raise Exception("NCs pass EpsilonEquals but have different control point counts!")


        #    nsC.Dispose()
        #    nsD.Dispose()

        #    if not epsEqual:
        #        break
        #else:
        #    return True

    return False


def _createCurves_1Set_VariousTols(rgObj_A, rgObj_B, fMaxTol, fMinTol, fMinCrvLength, bDebug):

    intersectionMethod = _getIntersectionMethod(rgObj_A, rgObj_B)

    fTol_Next = fMaxTol

    rgCs_Res = []
    fTols = []

    while True:
        sc.escape_test()

        if not fTol_Next:
            return rgCs_Res, fTols

        fTol_WIP = fTol_Next
        fTol_Next = fTol_WIP / 2.0

        if fTol_WIP == fMinTol:
            fTol_Next = None
        elif fTol_Next < (fMinTol + 0.1*fMinTol):
            fTol_Next = fMinTol

        if bDebug: sEval = "fTol_WIP"; print("{}: {}".format(sEval, eval(sEval)))

        bSuccess, rgCrvs_SSX_WIP, rgPts_SSX = intersectionMethod(
            rgObj_A,
            rgObj_B,
            fTol_WIP)

        # At least in RhinoCommon of 7 SR10,
        # even when no curves or points result, the result of
        # Intersection.BrepBrep and
        # Intersection.BrepSurface
        # is / (may still be) True, but
        # Intersection.SurfaceSurface is False.
        if not bSuccess:
            continue
            #return rgCs_Res, fTols

        if len(rgPts_SSX) > 0:
            print("{} points in intersection of tolerance {}.".format(
                rgPts_SSX,
                fTol_WIP))
        if len(rgCrvs_SSX_WIP) == 0:
            continue


        rgCrvs_SSX = []
        for c in rgCrvs_SSX_WIP:
            length = c.GetLength()
            #if length < 1.0:
            #    print(length
            if length > fMinCrvLength:
                rgCrvs_SSX.append(c)

        #rgCrvs_SSX = rgCrvs_SSX_WIP

        #for c in rgCrvs_SSX_WIP:
        #    sc.doc.Objects.AddCurve(c)
        #sc.doc.Views.Redraw(); 1/0

        if not rgCs_Res or not _areCurvesAlreadyInList(rgCrvs_SSX, rgCs_Res, fMaxTol, bDebug):
            rgCs_Res.append(rgCrvs_SSX)
            fTols.append(fTol_WIP)


def _createCurves_4Sets(rgB_A, rgB_B, fMaxTol, fMinTol, fMinCrvLength, bDebug):

    d_rgCs_Res = {}
    d_fTols = {}

    for iA, iB in (0, 0), (0, 1), (1, 0), (1, 1):

        bUnder_A = bool(iA)
        bUnder_B = bool(iB)

        key = bool(iA), bool(iB)

        rgObj_A = rgB_A.Faces[0].UnderlyingSurface() if bUnder_A else rgB_A
        rgObj_B = rgB_B.Faces[0].UnderlyingSurface() if bUnder_B else rgB_B

        rc = _createCurves_1Set_VariousTols(rgObj_A, rgObj_B, fMaxTol, fMinTol, fMinCrvLength, bDebug)

        d_rgCs_Res[key] = rc[0]
        d_fTols[key] = rc[1]

    return (
        d_rgCs_Res,
        d_fTols)


def main():

    objrefs_In = getInput()
    if objrefs_In is None: return

    bUseDialog = Opts.values['bUseDialog']
    fTol = Opts.values['fTol']
    fMinTol = Opts.values['fMinTol']
    fMinCrvLength = Opts.values['fMinCrvLength']
    bUnderlying1st = Opts.values['bUnderlying1st']
    bUnderlying2nd = Opts.values['bUnderlying2nd']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if not bUseDialog:
        srfA = objrefs_In[0].Face().UnderlyingSurface() if bUnderlying1st else objrefs_In[0].Face().DuplicateFace(duplicateMeshes=False)
        srfB = objrefs_In[1].Face().UnderlyingSurface() if bUnderlying2nd else objrefs_In[1].Face().DuplicateFace(duplicateMeshes=False)
        rgCs_Res = createCurves_NoLoops(
            srfA,
            srfB,
            fTol,
            fMinCrvLength)
        if not rgCs_Res:
            print("No valid solution.")
            return
        for rgC_Res in rgCs_Res:
            sc.doc.Objects.AddCurve(rgC_Res)
        sc.doc.Views.Redraw()
        return


    rgB_ToX_A = objrefs_In[0].Face().DuplicateFace(duplicateMeshes=False)
    rgB_ToX_B = objrefs_In[1].Face().DuplicateFace(duplicateMeshes=False)



    Rhino.RhinoApp.SetCommandPrompt("Calculating curves ...")

    rc = _createCurves_4Sets(rgB_ToX_A, rgB_ToX_B, fTol, fMinTol, fMinCrvLength, bDebug)

    (
        d_rgCs_Res,
        d_fTols,
        ) = rc


    if not any(len(d_rgCs_Res[key]) for key in d_rgCs_Res):
        print("No intersections were calculated.")
        return


    form = EtoModelessForm()
    form.addDataToForm(d_rgCs_Res, d_fTols)
    form.Owner = Rhino.UI.RhinoEtoApp.MainWindow
    form.Show()


if __name__ == '__main__': main()