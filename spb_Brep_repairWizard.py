"""
160705-160728: Created.
...
200130: Now maximum edge tolerance is used for finding matching bounding boxes from split.
        Added routine to remove slits from faces.
        Now rebuilds edges of remaining brep before extracting those to be rebuilt.
200202-24: Further development.
200313-14: Printed feedback change.
200319, 0520,21, 26, 0619, 0701,29, 1211: Import-related update.
210122: Bug fixes.
210426, 220118, 220315
220317: Imported a function.  Import-related update.
220420: Repaired for Rhino 7+ due to a RhinoCommon 7.17 script-breaking change.
220914-15, 221122, 240402, 250123, 0215: Import-related updates.

TODO:
Continuing replacing xBrepFace_trimToNakedEdges with other modules.

Add variable for correctEdgeTols?
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Diagnostics import Stopwatch
from System.Drawing import Color

import Eto.Drawing as drawing
import Eto.Forms as forms

import spb_Brep_correctEdgeTolerances
import spb_Brep_edgeFromTrimDeviations
import spb_Brep_Faces_withShortBorderLength
import xBrep_findFaceByBorderBoundingBox
import xBrep_invalid
import spb_Brep_Join
import spb_Brep_nakedEdgeLoop
import xBrep_rebuildEdges
import spb_Brep_Faces_Sliver
import spb_Brep_findFaces_SmallArea
import xBrep_splitSurfaceWithCurves
import xBrepObject
import xCurve
import xSurface


stopwatch = Stopwatch() # Will just use the same instance for all tests.


class Opts:

    data = (
        ('sBrepSelPool', 'Selection',
            ('All normal srfs/polysrfs  ', 'Selection',)),
        ('bRetrim_Bad', True, 'and retrim',),
        ('bFind_SlowArea', True, 'Area compute time (seconds) > ',),
        ('fSeconds_forArea', 1.0, ''),
        ('bCorrectETols', True, 'Correct edge tol. error > (% of absolute tol.)'),
        ('fPercent_EdgeTolMismatchTol', 100.0,''),
        ('bRemoveSlits', True, 'Remove single-edged slits'),
        ('bRebuildEs_perETol', True, 'Rebuild edges whose deviation > '),
        ('fTol_EdgeForRebuild', 1.0*sc.doc.ModelAbsoluteTolerance,''),
        ('bShrinkBeforeUntrim', True, 'Shrink surfaces before trimming',),
        ('bDelSmall', True, 'Delete faces < ',),
        ('fTol_Small', 2.0*sc.doc.ModelAbsoluteTolerance, ''),
        ('bJoin', True, 'Join breps',),
        ('fTol_EdgeForJoin', 2.0*sc.doc.ModelAbsoluteTolerance, 'Tolerance  '),
        ('bRetrim_GoodMonoface', False, 'Retrim good monofaces remaining after initial join'),
        ('bMergeEdges', True, 'Merge all edges if solid result',),
        ('bDebug', False, 'Debug mode',),
    )
    # TODO: Add to list?
    # 'All in document' to 'sBrepSelPool'
    # ('bSlitEdge', True, 'Coincident non-seam edges (slits)'),

    keys = []
    values = {}
    sDialogTexts = {}
    stickyKeys = {}

    for key, value, text in data:
        if key[0] == 'b':
            values[key] = value
            sDialogTexts[key] = text
            stickyKeys[key] = '{}({})'.format(key, __file__)
        elif key[0] == 's':
            values[key] = value
            sDialogTexts[key] = text
            stickyKeys[key] = '{}({})'.format(key, __file__)
        elif key[0] == 'f':
            values[key] = value
            sDialogTexts[key] = text
            stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
        else:
            print "{} not supported.".format(key)
            continue
        keys.append(key)


class RepairBrepDialog(forms.Dialog[bool]):


    def __init__(self):

        # Load sticky.
        for key in Opts.stickyKeys:
            if Opts.stickyKeys[key] in sc.sticky:
                Opts.values[key] = sc.sticky[Opts.stickyKeys[key]]

        # Initialize dialog box
        self.Title = 'Repair Brep (Surface and/or Polysurface)'
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
                self.checkBoxes[key] = forms.CheckBox(Text=Opts.sDialogTexts[key])
                self.checkBoxes[key].Checked = Opts.values[key]
            elif key[0] == 's':
                self.radioButtonList[key] = forms.RadioButtonList()
                self.radioButtonList[key].DataStore = Opts.sDialogTexts[key]
                self.radioButtonList[key].SelectedIndex = Opts.sDialogTexts[key].index(Opts.values[key])
            if key[0] == 'f':
                self.labels[key] = forms.Label(Text=Opts.sDialogTexts[key])
                self.textBoxes[key] = forms.TextBox(Text=str(Opts.values[key]))


        self.checkBoxes['bFind_SlowArea'].CheckedChanged += self.onChange_bSlowArea
        self.textBoxes['fSeconds_forArea'].TextChanged += self.onChange_fSeconds_forArea
        self.checkBoxes['bRebuildEs_perETol'].CheckedChanged += self.onChange_bRebuildEs_perETol
        self.textBoxes['fTol_EdgeForRebuild'].TextChanged += self.onChange_fTol_EdgeForRetrim
        self.checkBoxes['bDelSmall'].CheckedChanged += self.onChange_bDelSmall
        self.textBoxes['fTol_Small'].TextChanged += self.onChange_fTol_Small
        self.checkBoxes['bJoin'].CheckedChanged += self.onChange_bJoin
        self.textBoxes['fTol_EdgeForJoin'].TextChanged += self.onChange_fTol_EdgeForJoin

        self.DefaultButton = forms.Button(Text = 'OK')
        self.DefaultButton.Click += self.onOkButtonClick

        self.AbortButton = forms.Button(Text = 'Cancel')
        self.AbortButton.Click += self.onCancelButtonClick


        layout = forms.DynamicLayout()
        layout.Spacing = drawing.Size(5, 5)

        #layout.AddSeparateRow()
        key = 'sBrepSelPool'
        self.radioButtonList[key].Orientation = forms.Orientation.Horizontal
        layout.AddSeparateRow("Process:  ", self.radioButtonList[key])
        layout.AddSeparateRow("")

        layout.AddSeparateRow("Extract    ", self.checkBoxes['bRetrim_Bad'], ":")

        dummyCB = forms.CheckBox(Text="Invalid")
        dummyCB.Checked = True
        dummyCB.Enabled = False
        layout.AddSeparateRow(dummyCB)
        dummyCB = forms.CheckBox(Text="Whose area cannot be computed")
        dummyCB.Checked = True
        dummyCB.Enabled = False
        layout.AddSeparateRow(dummyCB)
        layout.AddSeparateRow(self.checkBoxes['bFind_SlowArea'], None,
                              self.textBoxes['fSeconds_forArea'])
        layout.AddSeparateRow("")

        layout.AddSeparateRow("Additional actions:")
        layout.AddSeparateRow(self.checkBoxes['bCorrectETols'], None,
                              self.textBoxes['fPercent_EdgeTolMismatchTol'])
        layout.AddSeparateRow(self.checkBoxes['bRemoveSlits'])
        layout.AddSeparateRow(self.checkBoxes['bRebuildEs_perETol'], None,
                              self.textBoxes['fTol_EdgeForRebuild'])
        #layout.AddRow(self.checkBoxes['bShrinkBeforeUntrim'])
        layout.AddSeparateRow(self.checkBoxes['bDelSmall'], None,
                              self.textBoxes['fTol_Small'])
        layout.AddSeparateRow(self.checkBoxes['bJoin'], None,
                              "Tolerance: ",
                              self.textBoxes['fTol_EdgeForJoin'])
        layout.AddSeparateRow(None, self.checkBoxes['bRetrim_GoodMonoface'],None)
        layout.AddSeparateRow("")
        layout.AddSeparateRow(self.checkBoxes['bDebug'])

        layout.AddSeparateRow(None, self.DefaultButton, None, self.AbortButton, None)
        #layout.AddRow(self.okButton, self.cancelButton)

        # Set the dialog content
        self.Content = layout


    def onChange_bSlowArea(self, sender, e):
        self.textBoxes['fSeconds_forArea'].Enabled = self.checkBoxes['bFind_SlowArea'].Checked


    def onChange_fSeconds_forArea(self, sender, e):
        try:
            float(self.textBoxes['fSeconds_forArea'].Text)
            self.textBoxes['fSeconds_forArea'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['fSeconds_forArea'].BackgroundColor = drawing.Colors.Red



    def onChange_bRebuildEs_perETol(self, sender, e):
        self.textBoxes['fTol_EdgeForRebuild'].Enabled = self.checkBoxes['bRebuildEs_perETol'].Checked


    def onChange_fTol_EdgeForRetrim(self, sender, e):
        try:
            float(self.textBoxes['fTol_EdgeForRebuild'].Text)
            self.textBoxes['fTol_EdgeForRebuild'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['fTol_EdgeForRebuild'].BackgroundColor = drawing.Colors.Red



    def onChange_bDelSmall(self, sender, e):
        self.textBoxes['fTol_Small'].Enabled = self.checkBoxes['bDelSmall'].Checked


    def onChange_fTol_Small(self, sender, e):
        try:
            float(self.textBoxes['fTol_Small'].Text)
            self.textBoxes['fTol_Small'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['fTol_Small'].BackgroundColor = drawing.Colors.Red


    def onChange_bJoin(self, sender, e):
        self.textBoxes['fTol_EdgeForJoin'].Enabled = self.checkBoxes['bJoin'].Checked


    def onChange_fTol_EdgeForJoin(self, sender, e):
        try:
            float(self.textBoxes['fTol_EdgeForJoin'].Text)
            self.textBoxes['fTol_EdgeForJoin'].BackgroundColor = drawing.SystemColors.ControlBackground
        except:
            self.textBoxes['fTol_EdgeForJoin'].BackgroundColor = drawing.Colors.Red


    def onCancelButtonClick(self, sender, e):
        self.Close(False)


    def onOkButtonClick(self, sender, e):
        # Save values.
        for key in Opts.keys:
            if key in self.checkBoxes:
                Opts.values[key] = self.checkBoxes[key].Checked
            elif key in self.radioButtonList:
                Opts.values[key] = (
                    Opts.sDialogTexts[key][self.radioButtonList[key].SelectedIndex]
                )
            elif key in self.textBoxes:
                try:
                    Opts.values[key] = float(self.textBoxes[key].Text)
                except:
                    print "Invalid input for tolerance." \
                        "  {} will be used instead.".format(Opts.values[key])


        # Save sticky.
        for key in Opts.keys:
            if key in Opts.stickyKeys:
                sc.sticky[Opts.stickyKeys[key]] = Opts.values[key]

        self.Close(True)


def getInput():

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")

    go.GeometryFilter = rd.ObjectType.Brep
    #    optD_Overlap_MinAllowed = ri.Custom.OptionDouble(
    #            fOverlap_MinAllowed, True, Rhino.RhinoMath.ZeroTolerance)
    #    optT_AddDots = ri.Custom.OptionToggle(bAddDot, 'No', 'Yes')
    #    optI_DotSize = ri.Custom.OptionInteger(iDotHeight, True, 3)
    #    optT_Extract = ri.Custom.OptionToggle(bExtract, 'No', 'Yes')
    #    go.AddOptionDouble('MinimumEdgeOverlapAllowed', optD_Overlap_MinAllowed)
    #    go.AddOptionToggle('AddTextDots', optT_AddDots)
    #    if optT_AddDots.CurrentValue:
    #        go.AddOptionInteger('TextDotHeight', optI_DotSize)
    #    go.AddOptionToggle('ExtractFaces', optT_Extract)
    go.EnableClearObjectsOnEntry(False)
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False
    go.EnableUnselectObjectsOnExit(False)
    go.SubObjectSelect = False
    #    go.AcceptNumber(True, False)

    while True:
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs


def processBrepObjects(rdBreps, **kwargs):
    """
    """

    gBreps = {} # Placed here so it can be found more easily in debugger.


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bRetrim_Bad = getOpt('bRetrim_Bad')
    bFind_SlowArea = getOpt('bFind_SlowArea')
    fSeconds_forArea = getOpt('fSeconds_forArea')
    bCorrectETols = getOpt('bCorrectETols')
    fPercent_EdgeTolMismatchTol = getOpt('fPercent_EdgeTolMismatchTol')
    bRemoveSlits = getOpt('bRemoveSlits')
    bRebuildEs_perETol = getOpt('bRebuildEs_perETol')
    fTol_EdgeForRebuild = getOpt('fTol_EdgeForRebuild')
    fTol_EdgeForJoin = getOpt('fTol_EdgeForJoin')
    bShrinkBeforeUntrim = getOpt('bShrinkBeforeUntrim')
    bDelSmall = getOpt('bDelSmall')
    fTol_Small = getOpt('fTol_Small')
    bJoin = getOpt('bJoin')
    bRetrim_GoodMonoface = getOpt('bRetrim_GoodMonoface')
    bMergeEdges = getOpt('bMergeEdges')
    bDebug = getOpt('bDebug')



    def coerceList(obj):
        try: out = list(set(obj))
        except: out = [obj]
        return out


    def coerceSurface(rhObj):
        if isinstance(rhObj, rg.Surface):
            return rhObj
        elif isinstance(rhObj, rg.GeometryBase):
            geom = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            #print rhObj.GeometryComponentIndex.ComponentIndexType
            geom = rhObj.Geometry()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
            geom = rdObj.Geometry
        else:
            return

        srf = None
        if isinstance(geom, rg.BrepFace):
            srf = geom.UnderlyingSurface()
        elif isinstance(geom, rg.Surface):
            srf = geom
        elif isinstance(geom, rg.Brep):
            if geom.Faces.Count == 1:
                srf = geom.Faces[0].UnderlyingSurface()

        return srf


    def getBrepGeom(rhObj):
        if isinstance(rhObj, rg.GeometryBase):
            geom = rhObj
            guid = None
        elif isinstance(rhObj, rd.ObjRef):
            geom = rhObj.Geometry()
            guid = rhObj.ObjectId
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
            if rdObj is None:
                print rhObj
                pass
            geom = rdObj.Geometry
            guid = rdObj.Id
        elif isinstance(rhObj, rd.BrepObject):
            rdObj = rhObj
            geom = rdObj.Geometry
            guid = rdObj.Id
        else:
            return

        if not isinstance(geom, rg.Brep):
            print "Not a brep: {}".format(guid)
            return

        return geom


    def coerceRhinoObject(rhObj):
        """
        'Deleted objects cannot be found by id.'
        (https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_DocObjects_Tables_ObjectTable_FindId.htm)
        """
        rdObj = None
        if isinstance(rhObj, rd.RhinoObject):
            rdObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        return rdObj


    def coerceBrepObject(rhObj):
        rdObj = coerceRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def cleanBrepList(toClean, alsoRemove=None):
        """
        Recreates list of GUIDS with dict key of toClean
        without duplicates, deleted breps, breps in dict key of alsoRemove.
        """

        # Remove duplicates.
        gBreps[toClean] = list(set(gBreps[toClean]))

        # Remove GUIDs of deleted Breps.
        # None from coerceBrepObject means that the Brep was deleted.
        gBreps[toClean] = [gB for gB in gBreps[toClean] if coerceBrepObject(gB) is not None]

        # Remove any GUIDs found in gBreps[alsoRemove].
        if alsoRemove and gBreps[alsoRemove]:
            gBreps[toClean] = [gB for gB in gBreps[toClean] if gB not in gBreps[alsoRemove]]


    def getGuidsToProcess(rdBreps):
        gBs = []
        len_rhBs0 = len(rdBreps)
        idxs_AtTenths = [int(round(0.1*i*len_rhBs0,0)) for i in range(10)]

        for iB, rdObj_Brep in enumerate(rdBreps):
            if sc.escape_test(False):
                print "Searching interrupted by user."
                return

            if len_rhBs0 > 10:
               if iB in idxs_AtTenths:
                    Rhino.RhinoApp.SetCommandPrompt(
                        "Processed {:d}% of {} breps ...".format(
                            int(100.0 * (iB+1) / len_rhBs0), len_rhBs0))
            elif len_rhBs0 > 1:
                Rhino.RhinoApp.SetCommandPrompt(
                    "Processing {} of {} breps ...".format(
                        iB+1, len_rhBs0))
            else:
                Rhino.RhinoApp.SetCommandPrompt("Processing single brep ...")

            rdBrep0 = coerceBrepObject(rdObj_Brep)
            if not rdBrep0:
                print "Non-brep DocObject passed to processBrepObjects."
                continue

            gBs.append(rdBrep0.Id)

        return gBs


    def colorRed(gObjs):
        gObjs = coerceList(gObjs)
        for gObj in gObjs:
            rdObj = coerceRhinoObject(gObj)
            attr = rdObj.Attributes
            attr.ColorSource = rd.ObjectColorSource.ColorFromObject
            attr.ObjectColor = Color.Red
            rdObj.CommitChanges()


    def formatDistance(fDistance):
        if fDistance is None: return "(No deviation provided)"
        if fDistance < 0.001:
            return "{:.2e}".format(fDistance)
        else:
            return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


    def removeMicroEdges(searchUs):
        sc.escape_test()
        if bEchoCmdPrompt:
            sCmdPrompt0 = "Removing micro edges ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)
        sc.doc.Views.RedrawEnabled = False
        for gBrep in gBreps[searchUs] if isinstance(searchUs, str) else searchUs:
            sc.doc.Objects.Select(gBrep)
        Rhino.RhinoApp.RunScript("_NoEcho _RemoveAllNakedMicroEdges _Echo", echo=False)
        sc.doc.Objects.UnselectAll()
        sc.doc.Views.RedrawEnabled = True
        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def extractInvalidFaces(keySearchUs, keyExtracted):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for invalid faces ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        if keyExtracted not in gBreps: gBreps[keyExtracted] = []

        rc = xBrep_invalid.extractBadFaces(
            gBreps[keySearchUs],
            bEcho=bDebug,
            bDebug=bDebug)

        if not rc:
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")
            return

        gBs_Bad, gBs_Good = rc

        # In case a micro edge is the problem.
        for gBrep in gBs_Bad:
            removeMicroEdges(gBs_Bad)
            for i in reversed(range(len(gBs_Bad))):
                if xBrep_invalid.isValid(gBs_Bad[i]):
                    gBs_Good.append(gBs_Bad[i])
                    del gBs_Bad[i]
                    print "MicroEdge(s) found" \
                        " after extracting invalid faces." \
                        "  Moved Good brep from Bad list to Good."
                    if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if bDebug:
            print "extractInvalidFaces:" \
                  " len(extracted) len(remaining): {}, {}".format(
                      len(gBs_Bad), len(gBs_Good))
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)
        if len(gBs_Bad) or bDebug:
            print "{} invalid faces found in {} {} breps.".format(
                len(gBs_Bad), len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if gBs_Bad:
            gBreps[keyExtracted].extend(gBs_Bad)
            gBreps[keySearchUs].extend(gBs_Good)
            cleanBrepList(keyExtracted)
            cleanBrepList(keySearchUs, keyExtracted)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def correctEdgeTols(keySearchUs):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for edges" \
                " whose tolerances don't match their actual edge deviations ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        iCt_Removed = 0

        gCorrected = spb_Brep_correctEdgeTolerances.processBrepObjects(
            gBreps[keySearchUs],
            fPercentOfMAT=fPercent_EdgeTolMismatchTol,
            bAddDot=False,
            bEcho=False,
            bDebug=bDebug,
            )

        if len(gCorrected) or bDebug:
            s  = "{} out of {} {} breps".format(
                len(gCorrected), len(gBreps[keySearchUs]), keySearchUs)
            s += " had some edge tolerances that were corrected."
            print s
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def removeSlitsFromFaces(keySearchUs):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for faces containing slit trims ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        iCt_Removed = 0

        for gB in gBreps[keySearchUs]:
            rgB = getBrepGeom(gB)

            if rgB.Faces.RemoveSlits():
                xBrepObject.replaceGeometry(gB, rgB)
                iCt_Removed += 1

            rgB.Dispose()

        if iCt_Removed or bDebug:
            print "{} out of {} {} breps had slit trims removed.".format(
                iCt_Removed, len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def extractFacesWithEdgeFromTrimTolError(keySearchUs, keyExtracted):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for faces with any " \
                "actual edge deviations not" \
                " matching BrepEdge.Tolerance values ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        if keyExtracted not in gBreps: gBreps[keyExtracted] = []

        gBs_1Fs_Found = []
        gBs_Remaining_Good = []
        iCt_BrepsFound = 0

        for gB in gBreps[keySearchUs]:
            rgB = getBrepGeom(gB)
            rgB.Compact()

            devs = spb_Brep_edgeFromTrimDeviations.getDeviations(rgB)

            idx_Fs_Found = []

            for iE in range(rgB.Edges.Count):
                if devs[iE] is None:
                    # TODO: Better handling of None?
                    idx_Fs_Found.extend(rgB.Edges[iE].AdjacentFaces())
                elif (
                    devs[iE] <= (0.01*fPercent_EdgeTolMismatchTol) * sc.doc.ModelAbsoluteTolerance and
                    rgB.Edges[iE].Tolerance <= (0.01*fPercent_EdgeTolMismatchTol) * sc.doc.ModelAbsoluteTolerance
                ):
                    pass
                elif (
                    abs(devs[iE] - rgB.Edges[iE].Tolerance) >
                    (0.01*fPercent_EdgeTolMismatchTol) * sc.doc.ModelAbsoluteTolerance
                ):
                    idx_Fs_Found.extend(rgB.Edges[iE].AdjacentFaces())

            if not idx_Fs_Found:
                rgB.Dispose()
                continue

            iCt_BrepsFound += 1

            idx_Fs_Found = sorted(set(idx_Fs_Found))

            rc = xBrepObject.extractFaces(
                gB,
                idx_Fs_Found,
                bAddOnlyMonofaces=True,
                bRetainLayer=True,
                bRetainColor=True,
                bEcho=False
            )
            if not rc:
                rgB.Dispose()
                continue

            gBs_Bad_Extracted, gBs_Good_Remaining = rc

            gBs_1Fs_Found.extend(gBs_Bad_Extracted)
            gBs_Remaining_Good.extend(gBs_Good_Remaining)

            rgB.Dispose()


        if len(gBs_1Fs_Found) or bDebug:
            s  = "{} faces extracted in {} {} breps".format(
                len(gBs_1Fs_Found), iCt_BrepsFound, keySearchUs)
            s += " with at least one edge tolerance over"
            s += " {:.1f}% of tolerance per BrepEdge.".format(
                fPercent_EdgeTolMismatchTol)
            print s
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if gBs_1Fs_Found:
            gBreps[keyExtracted].extend(gBs_1Fs_Found)
            cleanBrepList(keyExtracted)
            gBreps[keySearchUs].extend(gBs_Remaining_Good)
            cleanBrepList(keySearchUs, keyExtracted)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def retrimFaceSrfsWithExistingEdges(keyRetrimUs):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Retrimming face underlying surfaces with existing edges ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()


        iCt_Replaced = 0

        gs_SingleBrepReturn = []
        gs_MultiBrepReturn = []

        for gB in gBreps[keyRetrimUs]:
            NEs = rs.coercebrep(gB).DuplicateNakedEdgeCurves(True, True) # nakedOuter, nakedInner
            rc = xBrep_splitSurfaceWithCurves.processBrepObjects(
                [gB],
                NEs,
                bSplitUnderlyingSrf=False,
                bOnlyUseCrvsOnSrf=False,
                bSplitToCrvSegs=True,
                fTolerance=sc.doc.ModelAbsoluteTolerance,
                bTryOtherTolsOnFail=True,
                bExplode=False,
                bExtract=False,
                bEcho=bDebug,
                bDebug=bDebug)
            if not rc: continue
            gB_Split = rc[0]

            rc = xBrep_findFaceByBorderBoundingBox.processBrepObject(
                gB_Split,
                NEs,
                fTolerance=sc.doc.ModelAbsoluteTolerance,
                bDelOtherFacesNotExtract=True,
                bEcho=bDebug,
                bDebug=bDebug)
            if not rc: continue

            iCt_Replaced += 1

            if isinstance(rc, Guid):
                gs_SingleBrepReturn.append(rc)
            elif isinstance(rc, list):
                if len(rc) > 1:
                    gs_MultiBrepReturn.extend(rc)


            #rc = xBrepFace_trimToNakedEdges.processObjectsWithSrfs(
            #    gB,
            #    rhBreps_UseNakedEdges=NEs,
            #    bSplitFullSurface=True,
            #    fSplitTol=sc.doc.ModelAbsoluteTolerance,
            #    bExtractModifiedFace=False,
            #    bDelTrimmingBrep=False,
            #    bEcho=bDebug,
            #    bDebug=bDebug)
            #if not rc: continue

            #iCt_Replaced += 1

            #if len(rc) > 1:
            #    gs_MultiBrepReturn.extend(rc)

        if gs_SingleBrepReturn:
            gBreps[keyRetrimUs].extend(gs_SingleBrepReturn)
        if gs_MultiBrepReturn:
            gBreps[keyRetrimUs].extend(gs_MultiBrepReturn)

        cleanBrepList(keyRetrimUs)

        if iCt_Replaced or bDebug:
            s  = "{} monofaces replaced".format(iCt_Replaced)
            s += " in {} breps".format(keyRetrimUs)
            print s
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def rebuildFacesWithLargeEdgeTol(keySearchUs):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Rebuilding faces with large edge tolerances ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()


        rc = xBrep_rebuildEdges.processBrepObjects(
            gBreps[keySearchUs],
            fRebuildNakedTol=fTol_EdgeForRebuild,
            bEcho=bDebug if len(gBreps[keySearchUs]) <= 10 else False,
            bDebug=bDebug,
            )

        if len(rc) or bDebug:
            print "{} out of {} {} breps had their edges rebuilt.".format(
                len(rc), len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def extractBadAreaFaces(keySearchUs, keyExtracted):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for bad-area faces ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()



        def findBadAreaFaces(rgBrep, fTime_MaxTol, bDebug=False):
            # Check whether brep is invalid.
            #    if not rgBrep.IsValid:
            #        print "Brep is invalid and will be skipped.  " \
            #                "_ExtractBadSrf and repair the bad faces first."
            #        if bDebug: sPrint = 'rgBrep.IsValid'; print sPrint + ':', eval(sPrint)
            #        return [], []

            idxFs_No = []
            idxFs_Slow = []

            for iF in range(rgBrep.Faces.Count):
                rgBrep_1Face = rgBrep.Faces[iF].DuplicateFace(False)

                stopwatch.Restart()

                area = rgBrep_1Face.GetArea()

                stopwatch.Stop()
                timeElapsed = stopwatch.Elapsed.TotalSeconds

                if area == 0: # Not None, for no area compute.
                    idxFs_No.append(iF)
                elif fTime_MaxTol and timeElapsed > fTime_MaxTol:
                    idxFs_Slow.append(iF)

                rgBrep_1Face.Dispose()

            return idxFs_No + idxFs_Slow



        gBs_Bad = []
        gBs_Good = []

        for gBrep in gBreps[keySearchUs]:
            rgBrep = getBrepGeom(gBrep)
            rgBrep.Compact()

            rc = findBadAreaFaces(
                    rgBrep=rgBrep,
                    fTime_MaxTol=fSeconds_forArea if bFind_SlowArea else None,
            )
            rgBrep.Dispose()
            if not rc: continue
            idx_Faces = rc

            rc = xBrepObject.extractFaces(
                gBrep,
                idx_Faces,
                bAddOnlyMonofaces=True,
                bRetainLayer=True,
                bRetainColor=True,
                bEcho=False
            )
            if not rc: continue

            gBs_Bad, gBs_Good = rc

            # In case a micro edge is the problem.
            removeMicroEdges(gBs_Bad)

            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

            for i in reversed(range(len(gBs_Bad))):
                rgBrep = getBrepGeom(gBs_Bad[i])
                findBadAreaFaces(
                    rgBrep=rgBrep,
                    fTime_MaxTol=fSeconds_forArea if bFind_SlowArea else None,
                )
                rgBrep.Dispose()
                if not rc:
                    gBs_Good.append(gBs_Bad[i])
                    del gBs_Bad[i]
                    print "Moved Good brep from Bad list to Good."
                    if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if len(gBs_Bad) or bDebug:
            print "{} bad area faces (that are otherwise valid) found in {} {} breps.".format(
                len(gBs_Bad), len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if gBs_Bad:
            gBreps[keyExtracted].extend(gBs_Bad)
            cleanBrepList(keyExtracted)
            gBreps[keySearchUs].extend(gBs_Good)
            cleanBrepList(keySearchUs, keyExtracted)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def deleteShortBorderFaces(keySearchUs):
        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for short-border faces ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        rc = spb_Brep_Faces_withShortBorderLength.processBrepObjects(
                gBreps[keySearchUs],
                fMaxLength=2.01*fTol_Small,
                bExtract=True,
                bEcho=False,
                bDebug=False)
        if not rc:
            print "spb_Brep_Faces_withShortBorderLength returned None."
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")
            return []

        gBs_Short, gBs_RemainingFromExtract = rc

        iCt_Deleted = 0

        for gB in gBs_Short:
            bSuccess = sc.doc.Objects.Delete(objectId=gB, quiet=False)
            if bSuccess: iCt_Deleted += 1

        if iCt_Deleted or bDebug:
            print "{} short border faces deleted in {} {} breps.".format(
                iCt_Deleted, len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        # Why?
        #gBs_RemainingFromExtract = [gB for gBs in gBs_RemainingFromExtract for gB in gBs]

        if iCt_Deleted:
            if gBs_RemainingFromExtract:
                gBreps[keySearchUs].extend(gBs_RemainingFromExtract)

            cleanBrepList(keySearchUs)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def deleteSliverFaces(keySearchUs):
        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for sliver faces ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        iCt_Deleted = 0
        gBs_RemainingFromExtract = []

        for gBrep in gBreps[keySearchUs]:
            rgBrep = getBrepGeom(gBrep)
            rc = spb_Brep_Faces_Sliver.getFaces(
                rgBrep,
                fMaxSliverWidth=fTol_Small,
                bSkipFacesWithShortEdges=False,
                bSkipSliverCheckOfShortEdges=True,
                fMaxShortEdgeLength=fTol_EdgeForJoin,
                bEntireFaceMustBeASliver=True
                )
            iCt_Faces = rgBrep.Faces.Count
            rgBrep.Dispose()
            idx_rgFaces_Found = rc[0]
            if not idx_rgFaces_Found: continue
            rc = xBrepObject.removeFaces(gBrep, idx_rgFaces_Found)
            if rc is None:
                print "{} sliver faces could not be deleted.".format(len(idx_rgFaces_Found))
                if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)
            else:
                iCt_Deleted += len(idx_rgFaces_Found)
                gBs_RemainingFromExtract.extend(rc)

        if iCt_Deleted or bDebug:
            print "{} sliver faces deleted in {} {} breps.".format(
                iCt_Deleted, len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if iCt_Deleted:
            if gBs_RemainingFromExtract:
                gBreps[keySearchUs].extend(gBs_RemainingFromExtract)

            cleanBrepList(keySearchUs)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def deleteSmallAreaFaces(keySearchUs):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for small area faces ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        iCt_Deleted = 0
        gBs_RemainingFromExtract = []

        for gBrep in gBreps[keySearchUs]:
            rgBrep = getBrepGeom(gBrep)
            rc = spb_Brep_findFaces_SmallArea.getFaces(
                rgBrep,
                fAreaMinTol=(1.01*fTol_Small)**2.0,
                bEcho=False,
                bDebug=False)
            iCt_Faces = rgBrep.Faces.Count
            rgBrep.Dispose()
            idx_rgFaces_Found = rc[0]
            if not idx_rgFaces_Found: continue
            rc = xBrepObject.removeFaces(gBrep, idx_rgFaces_Found)
            if rc is None:
                print "{} small area faces could not be deleted.".format(len(idx_rgFaces_Found))
                if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)
            else:
                iCt_Deleted += len(idx_rgFaces_Found)
                gBs_RemainingFromExtract.extend(rc)

        if iCt_Deleted or bDebug:
            print "{} small area faces deleted in {} {} breps.".format(
                iCt_Deleted, len(gBreps[keySearchUs]), keySearchUs)
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        if iCt_Deleted:
            if gBs_RemainingFromExtract:
                gBreps[keySearchUs].extend(gBs_RemainingFromExtract)

            cleanBrepList(keySearchUs)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def joinBrepNakedEdges(keySearchUs):
        if bEchoCmdPrompt:
            sCmdPrompt0 = "Joining naked edges of breps ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        gBs_Remaining = []

        iCt_Joined_Total = 0

        for gBrep in gBreps[keySearchUs]:
            rgBrep = getBrepGeom(gBrep)
            if rgBrep.Faces.Count == 1:
                rgBrep.Dispose()
                continue
            iCt_EdgesJoined = rgBrep.JoinNakedEdges(fTol_EdgeForJoin)

            if not rgBrep.IsValid:
                print "  Joining naked edges of brep resulted in an invalid brep."
                if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)
                rgBrep.Dispose()
                continue

            if not iCt_EdgesJoined:
                rgBrep.Dispose()
                continue

            iCt_Joined_Total += iCt_EdgesJoined

            rc = xBrepObject.replaceGeometry(gBrep, rgBrep)
            if rc is None:
                print "  Brep object's geometry could not be replaced.".format(len(idx_rgFaces_Found))
                if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

            rgBrep.Dispose()

        if iCt_Joined_Total or bDebug:
            print "{} edge pairs joined in {} {} breps.".format(
                iCt_Joined_Total, len(gBreps[keySearchUs]), keySearchUs)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def filterPerSolidState(keySearchUs, keySolids, keyNotSolids=None):
        for gB in gBreps[keySearchUs]:
            rgB = getBrepGeom(gB)
            if rgB.IsSolid:
                gBreps[keySolids].append(gB)
            else:
                if keyNotSolids is not None: gBreps[keyNotSolids].append(gB)
            rgB.Dispose()
        cleanBrepList(keySearchUs, keySolids)
        if keyNotSolids: cleanBrepList(keySearchUs, keyNotSolids)


    def joinBreps(keysJoinUs, keySolids, keyNotSolids):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Joining breps ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()


        if isinstance(keysJoinUs, str):
            gBs = gBreps[keysJoinUs]
            keysJoinUs = [keysJoinUs]
        else:
            gBs = []
            for key in keysJoinUs:
                gBs.extend(gBreps[key])


        iCt_Joined = 0
        rc = spb_Brep_Join.joinBrepObjects(
            gBreps0=gBs,
            fJoinTol=fTol_EdgeForJoin,
            bByLayer=True,
            bByColor=True,
            bUseUiJoinCmd=False,
            bEcho=False,
            bDebug=False
        )
        # TODO: Add return GUIDs for bUseUiJoinCmd==True version of spb_Brep_Join.joinBrepObjects.
        #if not rc:
        #    # In case bUseUiJoinCmd == True can better handle the join.:
        #        rc = spb_Brep_Join.joinBrepObjects(
        #            gBreps0=gBs,
        #            fJoinTol=fTol_EdgeForJoin,
        #            bByLayer=True,
        #            bByColor=True,
        #            bUseUiJoinCmd=True,
        #            bEcho=False,
        #            bDebug=False
        #        )

        print "{} breps resulted from join of {} {} breps.".format(
            len(rc), len(gBs), ' and '.join(keysJoinUs))
        if not rc:
            if bDebug: print "No breps were joined."
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")
            return

        if bDebug: print "len(Output of joinBreps): {}".format(len(rc))

        gBreps['newly-joined'] = rc

        for key in keysJoinUs:
            cleanBrepList(key, 'newly-joined')

        filterPerSolidState('newly-joined', keySolids, keyNotSolids)

        for key in keysJoinUs:
            if key == keyNotSolids: continue
            gBreps[keyNotSolids].extend(gBreps[key])
            gBreps[key] = []

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def separateMonofaces(keysSearchUs, keyExtracted):
        """
        Returns list of GUIDs of monoface breps.
        """

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Searching for monoface breps ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()


        if isinstance(keysSearchUs, str):
            gBs = gBreps[keysSearchUs]
            keysSearchUs = [keysSearchUs]
        else:
            gBs = []
            for key in keysSearchUs:
                gBs.extend(gBreps[key])


        gBs_Mono_Out = []
        iCt_Fs_All = 0

        for gBrep in gBs:
            rgBrep = getBrepGeom(gBrep)
            rgBrep.Compact()

            if rgBrep.Faces.Count == 1:
                gBs_Mono_Out.append(gBrep)
                iCt_Fs_All += 1
            else:
                iCt_Fs_All += rgBrep.Faces.Count

            rgBrep.Dispose()

        if not gBs_Mono_Out:
            if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")
            return
        else:
            s = "{} out of {} {} breps are monoface.".format(
                len(gBs_Mono_Out), len(gBs), ' and '.join(keysSearchUs))
            #            if len(gBs_Mono_Out) == 1 and len(gBs) > 1:
            #                pass
            #            elif len(gBs_Mono_Out) > 100 and len(gBs_Mono_Out) >= 0.1 * float(iCt_Fs_All):
            #                if bDebug:
            #                    print "  Since this is {}% of all faces," \
            #                        " they will not be retrimmed ONLY due to being monofaces" \
            #                        "  during this call to separateMonofaces.".format(
            #                            100.0 * float(len(gBs_Mono_Out)) / float(iCt_Fs_All))
            #                print s + "  They will not be retrimmed at this time."
            #                return
            print s

        gBreps[keyExtracted].extend(gBs_Mono_Out)

        cleanBrepList(keyExtracted)

        for key in keysSearchUs:
            cleanBrepList(key, keyExtracted)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def shrinkFaces(keyShrinkUs):
        gBs_Shrunk = []
        for gBrep in gBreps[keyShrinkUs]:
            sc.escape_test()
            rdBrep = coerceBrepObject(gBrep)
            rgBrep = rdBrep.Geometry
            rgBrep.Compact()
            if rgBrep.Faces.ShrinkFaces():
                if sc.doc.Objects.Replace(objectId=gBrep, brep=rgBrep):
                    gBs_Shrunk.append(gBrep)
            rgBrep.Dispose()
        print "  {} faces were shrunk.".format(len(gBs_Shrunk))


    def untrimFaces(keyUntrimUs, keyUntrimmed):

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Untrimming monoface breps ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()

        gBreps_Out = []

        for gB in gBreps[keyUntrimUs]:
            rdBrep = coerceBrepObject(gB)
            rgBrep = rdBrep.Geometry
            rgBrep.Compact()
            bFail = False
            for rgSrf in rgBrep.Surfaces:
                brep = rgSrf.ToBrep()
                gBrep1 = sc.doc.Objects.AddBrep(brep, attributes=rdBrep.Attributes)
                if gBrep1 != Guid.Empty:
                    gBreps_Out.append(gBrep1)
                else:
                    print "Untrimmed Brep could not be added for {}." \
                        "  Check this.".format(gB)
            if not bFail:
                sc.doc.Objects.Delete(objectId=gB, quiet=False)
            rgBrep.Dispose()
        if bDebug: print "len(Output of untrimFaces): {}".format(len(gBreps_Out))
        print "{} faces of {} {} breps were untrimmed.".format(
            len(gBreps_Out), len(gBreps[keyUntrimUs]), keyUntrimUs)
        if gBreps_Out:
            gBreps[keyUntrimmed].extend(gBreps_Out)
            cleanBrepList(keyUntrimUs, keyUntrimmed)

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def deleteBreps(keyToDelete):
        for gB in gBreps[keyToDelete]:
            sc.escape_test()
            sc.doc.Objects.Delete(objectId=gB, quiet=False)
        cleanBrepList(keyToDelete)


    def findMatchingBrep(breps, crvs, tolerance=None, bDebug=False):
        """
        This function was copied from another script.
        
        Returns:
            Index of brep in breps.
    
        This function is used by other scripts, so leave it at the root level.
        """
    
        if tolerance is None:
            tolerance = 2.0 * sc.doc.ModelAbsoluteTolerance
    
    
    
        def createBB_ofGeometry(geoms):
            try: geoms = list(geoms)
            except: geoms = [geoms]
            bbox = rg.BoundingBox.Empty
            for geom in geoms: bbox.Union(geom.GetBoundingBox(True))
            if bbox is None: raise ValueError("Bounding box not created.")
            return bbox
    
    
        def findMatchingBoundingBox(bboxes, bbox_toMatch):
            """
            bboxes: Iterable of BoundingBoxes from which to search.
            bbox_toMatch: BoundingBox to match.
            Returns: Index of single matching bounding box on success.
            """
            
            if bDebug: print 'findMatchingBoundingBox()...'
            
            #tolerance = (3.0 * (tol_Split**2.0))**0.5 # == Total distance of tol_Split for each component.
            
            idx_Found = None
    
            epsilon = tolerance
    
            for i in range(len(bboxes)):
                if (
                        bbox_toMatch.Min.EpsilonEquals(bboxes[i].Min, epsilon=epsilon) and
                        bbox_toMatch.Max.EpsilonEquals(bboxes[i].Max, epsilon=epsilon)
                ):
                    if idx_Found is not None:
                        print "More than one match found.  Check results."
                        return
                    idx_Found = i
            
            if idx_Found is None:
                if bDebug:
                    print "Matching bounding box NOT found."
                    #sc.doc.Objects.AddBrep(bbox_toMatch.ToBrep())
                    #for bb in bboxes:
                        #sc.doc.Objects.AddBrep(bb.ToBrep())
    
            return idx_Found
    
    
        crvs_Closed = []
        segs_ofOpen = []
        for crv in crvs:
            if crv.IsClosed:
                crvs_Closed.append(crv)
            else:
                segs_ofOpen.append(crv)
    
        segs_Joined = rg.Curve.JoinCurves(
                segs_ofOpen, joinTolerance=sc.doc.ModelAbsoluteTolerance)
        segs_Joined = list(segs_Joined) # Convert from Array.
    
        bbox_of_crvs = createBB_ofGeometry(crvs_Closed + segs_Joined)
    
        #sc.doc.Objects.AddBox(rg.Box(bbox_of_crvs)); #sc.doc.Views.Redraw(); 1/0
    
        bboxes_ofNEs = []
        for brep in breps:
            NEs = [edge for edge in brep.Edges
                   if edge.Valence == rg.EdgeAdjacency.Naked]
            bboxes_ofNEs.append(createBB_ofGeometry(NEs))
    
        #for bbox in bboxes_ofNEs:
            #sc.doc.Objects.AddBox(rg.Box(bbox))
        #sc.doc.Views.Redraw(); 1/0
    
        #bboxes_ofNEs = [
        #        createBB_ofGeometry(
        #            [edge
        #                for edge in brep.Edges
        #                if edge.Valence == rg.EdgeAdjacency.Naked]
        #        )
        #        for brep in breps]
    
        #map(sc.doc.Objects.AddBox, map(rg.Box, bboxes_ofNEs))
    
        return findMatchingBoundingBox(bboxes_ofNEs, bbox_of_crvs)


    def trimUnderlyingSrfsToBrepNEs(keyToTrim, keysOfTrimming):
        """
        Parameters:
            keyToTrim: Breps of any rg.Surface, not just rg.BrepFace.
            keyOfTrimming: Can be one or more GUIDs of Breps.
        Returns:
            4 lists of GUIDs of Breps on success:
                New trimmed.
                Old that were trimmed.
                New split (trim fails).
                Old that were split (trim fails).
        """

        if bEchoCmdPrompt:
            sCmdPrompt0 = "Trimming surfaces to brep naked edges ..."
            Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

        sc.escape_test()



        NBs_All = []
        for gB0 in gBreps[keysOfTrimming]:
            for border in (
                spb_Brep_nakedEdgeLoop.createClosedCrvs_2EdgesPerVertexOnly(getBrepGeom(gB0))
            ):
                NBs_All.append(border)

        #for NB in NBs_All:
        #    sc.doc.Objects.AddCurve(NB)
        #sc.doc.Views.Redraw(); 1/0

        iCt_Trimmed = 0
        rgBs_Trimmed_to1F = []
        gBs_Trimmed_to1F = []
        gBs_beforeTrim = []
        gBs_Split = []
        gBs_beforeSplit = []


        for gB0 in gBreps[keyToTrim]:
            rdBrep0 = coerceBrepObject(gB0)
            rgBrep0 = rdBrep0.Geometry
            rgBrep0.Compact()
            rgSrf01 = rgBrep0.Surfaces[0]

            rc = xCurve.filterCurvesOnSurface(
                NBs_All,
                coerceSurface(gB0),
                bDebug=bDebug)
            if not rc: continue
            (
                rgCrvs_CompletelyOnSrf,
                rgCrvs_PartiallyOnSrf,
                ) = rc


            # TODO?: Add curve trimming here.


            # Surface will be split only with segments, not joined curves.
            segs_CompletelyOnSrf_Closed = []
            rgCrvs_CompletelyOnSrf_Closed = []
            for c in rgCrvs_CompletelyOnSrf:
                if c.IsClosed:
                    rgCrvs_CompletelyOnSrf_Closed.append(c)
            segs_CompletelyOnSrf_Closed = xCurve.duplicateSegments(
                rgCrvs_CompletelyOnSrf,
                bExplodePolyCrvs=True)
            #segs_CompletelyOnSrf_Open = xCurve.duplicateSegments(
            #    rgCrvs_CompletelyOnSrf_Open,
            #    bExplodePolyCrvs=True)
            #segs_PartiallyOnSrf = xCurve.duplicateSegments(
            #    rgCrvs_PartiallyOnSrf,
            #    bExplodePolyCrvs=True)

            rgCrvs_Splitters = segs_CompletelyOnSrf_Closed


            # This uses a loop to try a few tolerances before giving up.
            # bOnlyUseCrvsOnFace=False because curves were already analyzed.
            rgB_fromSplit = xSurface.splitSurfaceIntoBrep(
                rgSrf01,
                rgCrvs_Splitters,
                bOnlyUseCrvsOnFace=False,
                fTolerance=sc.doc.ModelAbsoluteTolerance,
                bTryOtherTolsOnFail=True,
                bShrinkSplitMonofaces=False,
                bDebug=bDebug)
            if not rgB_fromSplit: continue


            rgBs_fromSplit = []
            for rgF in rgB_fromSplit.Faces:
                rgB = rgF.DuplicateFace(False)
                rgBs_fromSplit.append(rgB)
            #sc.doc.Views.Redraw(); 1/0

            # Get maximum edge tolerance.  TODO: Make this a function.
            maximaEdgeTols = []
            for brep in rgBs_fromSplit:
                maximaEdgeTols.append(max([edge.Tolerance for edge in brep.Edges]))
            maxEdgeTol = max(maximaEdgeTols)

            idxBrep_Split_Match = findMatchingBrep(
                rgBs_fromSplit,
                segs_CompletelyOnSrf_Closed,
                tolerance=max(1.1*maxEdgeTol, fTol_EdgeForJoin),
                bDebug=bDebug)

            for c in segs_CompletelyOnSrf_Closed: c.Dispose()

            if idxBrep_Split_Match is None:
                if rgB_fromSplit.Faces.Count > 1:
                    bBrep_Split = sc.doc.Objects.AddBrep(
                        rgB_fromSplit, attributes=rdBrep0.Attributes)
                    if bBrep_Split == Guid.Empty:
                        raise ValueError("Untrimmed Brep could not be added for {}." \
                            "  Check this.".format(gB0))
                    gBs_Split.append(bBrep_Split)
                    gBs_beforeSplit.append(gB0)
            else:
                rgB_fromTrim = rgBs_fromSplit[idxBrep_Split_Match]
                gBrep_Trimmed = sc.doc.Objects.AddBrep(
                    rgB_fromTrim, attributes=rdBrep0.Attributes)
                if gBrep_Trimmed != Guid.Empty:
                    gBs_Trimmed_to1F.append(gBrep_Trimmed)
                    gBs_beforeTrim.append(gB0)
                else:
                    raise ValueError("Trimmed Brep could not be added for {}." \
                        "  Check this.".format(gB0))

            rgB_fromSplit.Dispose()
            for rgB in rgBs_fromSplit: rgB.Dispose()

        if len(gBs_Trimmed_to1F) or bDebug:
            print "{} faces were trimmed to a single face from {} {} breps.".format(
                len(gBs_Trimmed_to1F), len(gBreps[keyToTrim]), keyToTrim)

        if len(gBs_Split) or bDebug:
            print "{} other faces were split to polyface breps.".format(
                len(gBs_Split))

        if gBs_Trimmed_to1F:
            gBreps['retrimmed'].extend(gBs_Trimmed_to1F)
            gBreps['before-trim'] = gBs_beforeTrim
            cleanBrepList(keyToTrim, 'before-trim')
            deleteBreps('before-trim')

        if gBs_Split:
            gBreps['split-only'].extend(gBs_Split)
            gBreps['before-split'] = gBs_beforeSplit
            # gBs_beforeSplit should == gBreps[keyToTrim].
            pass

        if bEchoCmdPrompt: Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0 + " done.")


    def resultReport():
        s  = "Final result:"
        s += " {} solid".format(len(gBreps['solid']))
        s += ","
        s += " {} good-but-open".format(len(gBreps['good-but-open']))
        s += ","
        s += " {} split-only".format(len(gBreps['split-only']))
        s += ", and"
        s += " {} untrimmed".format(len(gBreps['untrimmed']))
        s += " breps"
        return s


    def resultReturn():
        return (
            gBreps['solid'],
            gBreps['good-but-open'],
            gBreps['split-only'],
            gBreps['untrimmed'],
            )



    gBreps['initial'] = getGuidsToProcess(rdBreps)

    bEcho = True

    bEchoCmdPrompt = True

    if bEchoCmdPrompt:
        sCmdPrompt0 = "Working ..."
        Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

    #bEchoCmdPrompt = len(gBreps['initial']) > 10


    sc.doc.Objects.UnselectAll()


    # Remove micro edges before all other actions.
    removeMicroEdges('initial')

    # Monofaces are not added to 'monoface-to-retrim'
    # at this time so that other breps can be joined,
    # resulting in naked edge borders better suited for trimming other faces.

    # Bad faces must be extracted now because
    # RhinoCommon does not allow Objects.AddBrep of some invalid breps.
    extractInvalidFaces('initial', 'monoface-to-retrim')

    gBreps['remaining'] = gBreps['initial'][:]

    if gBreps['monoface-to-retrim']:
        removeMicroEdges('remaining')

    if bCorrectETols:
        correctEdgeTols('remaining')

    if bRemoveSlits:
        removeSlitsFromFaces('remaining')

    if bCorrectETols:
        extractFacesWithEdgeFromTrimTolError('remaining', 'monoface-to-retrim-to-itself')
        retrimFaceSrfsWithExistingEdges('monoface-to-retrim-to-itself')
        gBreps['remaining'].extend(gBreps['monoface-to-retrim-to-itself'])

    if bRebuildEs_perETol:
        rebuildFacesWithLargeEdgeTol('remaining')

    # This is checked last since calculating area takes longer than preceding tests.
    extractBadAreaFaces('remaining', 'monoface-to-retrim')


    if bDelSmall:
        if gBreps['remaining']:
            deleteShortBorderFaces('remaining')
            if gBreps['remaining']:
                deleteSliverFaces('remaining')
                if gBreps['remaining']:
                    deleteSmallAreaFaces('remaining')

        if gBreps['monoface-to-retrim']:
            deleteShortBorderFaces('monoface-to-retrim')
            if gBreps['monoface-to-retrim']:
                deleteSliverFaces('monoface-to-retrim')
                if gBreps['monoface-to-retrim']:
                    deleteSmallAreaFaces('monoface-to-retrim')


    # This can reduce connecting naked edge borders.
    joinBrepNakedEdges('remaining')


    gBreps['solid'] = []
    gBreps['good-but-open'] = []
    gBreps['split-only'] = []
    gBreps['untrimmed'] = []
    gBreps['retrimmed'] = []
    gBreps['split-only'] = []
    gBreps['before-split'] = []


    filterPerSolidState('remaining', 'solid', 'good-but-open')

    if not gBreps['good-but-open'] and gBreps['solid']:
        s  = "{} breps were already closed.".format(len(gBreps['solid']))
        s += "  {} 'bad' extracted faces remain.".format(len(gBreps['monoface-to-retrim']))
        if not gBreps['monoface-to-retrim']:
            s += "  Model was not modified."
        print s
        print resultReport()
        return resultReturn()


    filterPerSolidState('good-but-open', 'solid')

    if not gBreps['monoface-to-retrim']:
        if not gBreps['good-but-open'] and gBreps['solid']:
            print "Only {} closed breps remain.".format(len(gBreps['solid']))
            print resultReport()
            return resultReturn()
        elif len(gBreps['good-but-open']) == 1:
            print "Only 1 open and {} closed breps remain.".format(len(gBreps['solid']))
            print resultReport()
            return resultReturn()


    # Initial Join.
    # Join remaining good-but-open breps with each other only.
    # The edge borders after the join will be more useful
    # as trimming curves for the other faces.

    if len(gBreps['good-but-open']) == 0:
        print "No good-but-open breps to join."
    elif len(gBreps['good-but-open']) == 1:
        print "Not enough good-but-open breps to join."
    else:
        joinBreps('good-but-open', 'solid', 'good-but-open')
        joinBrepNakedEdges('good-but-open')
        filterPerSolidState('good-but-open', 'solid')

        # If no monofaces need to be retrimmed, then the current result is final.
        if not gBreps['monoface-to-retrim'] and not bRetrim_GoodMonoface:
            if gBreps['good-but-open'] or gBreps['solid']:
                print "After a single iteration of joining," \
                    " only {} open and {} closed breps remain.".format(
                        len(gBreps['good-but-open']), len(gBreps['solid']))
                print resultReport()
                return resultReturn()

    if bRetrim_GoodMonoface:
        if gBreps['monoface-to-retrim']:
            print "Since bad breps were" \
                " found," \
                " good monoface breps will" \
                " not" \
                " be collected to be retrimmed at this time."
        else:
            print "Since bad breps were" \
                " not" \
                " found," \
                " good monoface breps will" \
                " be collected to be retrimmed at this time."
            separateMonofaces('good-but-open', 'monoface-to-retrim')



    print "__Retrim and join:"

    if gBreps['monoface-to-retrim'] and gBreps['good-but-open']:
        untrimFaces('monoface-to-retrim', 'untrimmed')

        trimUnderlyingSrfsToBrepNEs('untrimmed', 'good-but-open')


        if gBreps['retrimmed']:
            deleteBreps('split-only') # Will not use splits at this time.

            if len(gBreps['good-but-open'] + gBreps['retrimmed']) >= 2:
                joinBreps(['good-but-open', 'retrimmed'], 'solid', 'good-but-open')
                joinBrepNakedEdges('good-but-open')
                filterPerSolidState('good-but-open', 'solid')

                if not gBreps['untrimmed'] and not bRetrim_GoodMonoface:
                    if gBreps['good-but-open'] or gBreps['solid']:
                        print resultReport()
                        return resultReturn()


            if bRetrim_GoodMonoface and len(gBreps['good-but-open'] + gBreps['untrimmed']) >= 2:
                print "__Find monofaces to retrim:"
                separateMonofaces(['good-but-open', 'untrimmed'], 'monoface-to-retrim')

            if gBreps['monoface-to-retrim'] and gBreps['good-but-open']:
                untrimFaces('monoface-to-retrim', 'untrimmed')

                trimUnderlyingSrfsToBrepNEs('untrimmed', 'good-but-open')

                if gBreps['split-only']:
                    # Leave split breps in model
                    # for the user to decide what to do with them,
                    # but delete untrimmed breps.
                    deleteBreps('before-split')
                    print "{} breps that were split from {} are in model.".format(
                        len(gBreps['split-only']), len(gBreps['before-split']))

                if gBreps['retrimmed']:
                    if len(gBreps['good-but-open'] + gBreps['retrimmed']) >= 2:
                        joinBreps(['good-but-open', 'retrimmed'], 'solid', 'good-but-open')
                        joinBrepNakedEdges('good-but-open')
                        filterPerSolidState('good-but-open', 'solid')

                        if not gBreps['untrimmed'] and not bRetrim_GoodMonoface:
                            if gBreps['good-but-open'] or gBreps['solid']:
                                print resultReport()
                                return resultReturn()


    if len(gBreps['good-but-open'] + gBreps['retrimmed']) >= 2:
        print "__Join newly trimmed with each other and other good, open breps:"
        joinBreps(['good-but-open', 'retrimmed'], 'solid', 'good-but-open')
        joinBrepNakedEdges('good-but-open')
        filterPerSolidState('good-but-open', 'solid')

    print resultReport()
    return resultReturn()


def main():

    dialog = RepairBrepDialog()
    if not dialog.ShowModal(Rhino.UI.RhinoEtoApp.MainWindow): return


    if Opts.values['sBrepSelPool'] == Opts.sDialogTexts['sBrepSelPool'][0]:
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False

        rdObjs = []
        for rdObj in sc.doc.Objects.GetObjectList(iter):
            if rdObj.ObjectType == rd.ObjectType.Brep:
                rdObjs.append(rdObj)
    elif Opts.values['sBrepSelPool'] == Opts.sDialogTexts['sBrepSelPool'][1]:
        rc = getInput()
        if rc is None: return
        #oRefs, fOverlap_MinAllowed, bAddDot, iDotHeight, bExtract = ret
        rdObjs = rc
    else:
        print "Not supported yet."
        return


    if not Opts.values['bDebug']: sc.doc.Views.RedrawEnabled = False

    processBrepObjects(rdObjs)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()