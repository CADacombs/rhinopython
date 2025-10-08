"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180730: Created.
...
241101: Replaced rs.IsLayer with sc.doc.Layers.FindByFullPath("Default", Rhino.RhinoMath.UnsetIntIndex) >= 0
        since the former isn't the full path.
        Removed an import.
        Expanded printed report.
251007: Bug fix in creating the Iso view.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import System

#from System.Drawing import Color


class Opts():
    
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    
    keys = []
    
    key = 'bHideGrids'
    keys.append(key)
    values[key] = True
    names[key] = 'HideGrids'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAllWireframe'
    keys.append(key)
    values[key] = True
    names[key] = 'AllVPsToWireframe'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bChangePerspectiveToIso'
    keys.append(key)
    values[key] = True
    names[key] = 'IsoFromPerspective'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSetAngleTol'
    keys.append(key)
    values[key] = True
    names[key] = 'SetAngleTol'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fAngleTol'
    keys.append(key)
    values[key] = 0.1
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSetPrec'
    keys.append(key)
    values[key] = True
    names[key] = 'SetPrecision'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iPrecision'
    keys.append(key)
    values[key] = 5 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches else 4
    names[key] = 'Precision'
    riOpts[key] = ri.Custom.OptionInteger(values[key], True, 0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bSetBlackLayersToWhite'
    keys.append(key)
    values[key] = True
    names[key] = 'BlackLayersToWhite'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bRenameLayer01ToDefault'
    keys.append(key)
    values[key] = True
    names[key] = 'RenameLayer01ToDefault'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddDefaultLayer'
    keys.append(key)
    values[key] = True
    names[key] = 'AddDefaultLayer'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'
    keys.append(key)
    values[key] = True
    names[key] = 'Echo'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    keys.append(key)
    values[key] = False
    names[key] = 'Debug'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if sc.sticky.has_key(stickyKeys[key]):
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    @classmethod
    def addOptions(cls, go, sInclude=None):
        if sInclude is None:
            sInclude = cls.keys[:]
        
        def addOptionDouble(key):
            if key in sInclude:
                go.AddOptionDouble(cls.names[key], cls.riOpts[key])
        
        def addOptionInteger(key):
            if key in sInclude:
                go.AddOptionInteger(cls.names[key], cls.riOpts[key])
        
        def addOptionToggle(key):
            if key in sInclude:
                go.AddOptionToggle(cls.names[key], cls.riOpts[key])
        
        addOptionToggle('bHideGrids')
        addOptionToggle('bAllWireframe')
        addOptionToggle('bChangePerspectiveToIso')
        addOptionToggle('bSetAngleTol')
        if cls.values['bSetAngleTol']:
            addOptionDouble('fAngleTol')
        addOptionToggle('bSetPrec')
        if cls.values['iPrecision']:
            addOptionInteger('iPrecision')
        addOptionToggle('bSetBlackLayersToWhite')
        addOptionToggle('bRenameLayer01ToDefault')
        addOptionToggle('bAddDefaultLayer')
        addOptionToggle('bEcho')
        addOptionToggle('bDebug')
    
    @classmethod
    def setValues(cls):
        for key in cls.stickyKeys:
            if cls.riOpts[key]:
                cls.values[key] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                pass
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if cls.riOpts[key]:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
    
    @classmethod
    def processInput(cls, go):
        res = go.Result()
        
        cls.setValues()
        cls.saveSticky()
        
        # Clear and add options regardless if a number was entered or options were modified in another way.
        go.ClearCommandOptions()
        Opts.addOptions(go)


def getInput():
    
    go = ri.Custom.GetOption()
    
    go.SetCommandPrompt("Options")
    
    go.AcceptNothing(True)
    
    Opts.addOptions(go)
    
    while True:
        sc.escape_test()
        
        res = go.Get()
        
        if res == ri.GetResult.Nothing:
            break
        elif res == ri.GetResult.Cancel:
            return
        
        Opts.processInput(go)
        
        # Clear and add options regardless if a number was entered or options were modified in another way.
        go.ClearCommandOptions()
        Opts.addOptions(go)
    
    return True


def setLayerColors_BlackToWhite(bDebug=False):

    sLayers = rs.LayerNames()
    if sLayers is None or len(sLayers) < 1:
        raise Exception("sLayers is None or len(sLayers) < 1")
    
    iCt_LayersToWhite_NonRef = 0
    iCt_LayersToWhite_Ref = 0
    
    sc.doc.Views.RedrawEnabled = False
    
    for sLayer in sLayers:
        if rs.IsLayerReference(sLayer) and '::' not in sLayer:
            if bDebug:
                print("Skipped main reference layer of insert.")
            continue

        color = rs.LayerColor(sLayer)
        if (
            color.R == color.Black.R and
            color.G == color.Black.G and
            color.B == color.Black.B
        ):
            if sc.doc.Views.RedrawEnabled:
                sc.doc.Views.RedrawEnabled = False
            rs.LayerColor(sLayer, color=color.White)
            if rs.IsLayerReference(sLayer):
                iCt_LayersToWhite_Ref += 1
            else:
                iCt_LayersToWhite_NonRef += 1
            if bDebug:
                print("{} changed from black to white.".format(sLayer))
    if not sc.doc.Views.RedrawEnabled: sc.doc.Views.RedrawEnabled = True
    
    return iCt_LayersToWhite_NonRef, iCt_LayersToWhite_Ref


def run(bHideGrids=None, bAllWireframe=None, bChangePerspectiveToIso=None, bSetAngleTol=None, fAngleTol=None, bSetPrec=None, iPrecision=None, bSetBlackLayersToWhite=None, bRenameLayer01ToDefault=None, bAddDefaultLayer=None, bEcho=None, bDebug=None):
    """
    """
    
    if bHideGrids is None: bHideGrids = Opts.values['bHideGrids']
    if bAllWireframe is None: bAllWireframe = Opts.values['bAllWireframe']
    if bChangePerspectiveToIso is None: bChangePerspectiveToIso = Opts.values['bChangePerspectiveToIso']
    if bSetAngleTol is None: bSetAngleTol = Opts.values['bSetAngleTol']
    if fAngleTol is None: fAngleTol = Opts.values['fAngleTol']
    if bSetPrec is None: bSetPrec = Opts.values['bSetPrec']
    if iPrecision is None: iPrecision = Opts.values['iPrecision']
    if bSetBlackLayersToWhite is None: bSetBlackLayersToWhite = Opts.values['bSetBlackLayersToWhite']
    if bRenameLayer01ToDefault is None: bRenameLayer01ToDefault = Opts.values['bRenameLayer01ToDefault']
    if bAddDefaultLayer is None: bAddDefaultLayer = Opts.values['bAddDefaultLayer']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    sLog_Change = []
    sLog_NoChange = []
    
    if bHideGrids:
        sViews = rs.ViewNames(return_names=True, view_type=0)
        
        bSomeGridsWerePreviouslyEnabled = False
        for sView in sViews:
            # 0 = standard model views
            if rs.ShowGrid(view=sView, show=None):
                bSomeGridsWerePreviouslyEnabled = True
                if not rs.ShowGrid(view=sView, show=False):
                    sLog_Change.append("Grid disabled in {}.".format(sView))
                else:
                    print("Grid could not be disabled in {}.".format(sView))
        if not bSomeGridsWerePreviouslyEnabled:
            sLog_NoChange.append("Grids in all viewports were already disabled.")
        # Replaced:
        #rs.Command("_-Grid _ApplyTo=AllViewports _ShowGrid=No _Enter", echo=bEcho)

    if bAllWireframe:
        iCt_ViewportDisplayModeChanges = 0
        #print(Rhino.Display.DisplayModeDescription.WireframeId)
        wf = None
        for view in sc.doc.Views:
            if view.ActiveViewport.DisplayMode.Id != Rhino.Display.DisplayModeDescription.WireframeId:
                if wf is None:
                    wf = wf = Rhino.Display.DisplayModeDescription.FindByName("Wireframe")
                view.ActiveViewport.DisplayMode = wf
                if view.ActiveViewport.DisplayMode.Id != Rhino.Display.DisplayModeDescription.WireframeId:
                    raise Exception("Could not make a viewport wireframe.")
                iCt_ViewportDisplayModeChanges += 1
        #        Rhino.RhinoApp.RunScript(
        #            script="_-SetDisplayMode _Viewport=All _Mode=Wireframe",
        #            echo=bDebug)
        if iCt_ViewportDisplayModeChanges:
            sLog_Change.append(
                "{} viewports' display modes were set to Wireframe.".format(
                    iCt_ViewportDisplayModeChanges))
        else:
            sLog_NoChange.append(
                "All viewports' display modes were already set to Wireframe.")

    if bChangePerspectiveToIso:
        sViews = rs.ViewNames(return_names=True, view_type=0)

        if "Iso" in sViews:
            sLog_NoChange.append("Iso already exists.")
        elif "Perspective" in sViews:
            if "Perspective" == rs.RenameView("Perspective", "Iso"):
                sLog_Change.append("Perspective viewport changed to Iso")
        else:
            new_view = sc.doc.Views.Add(
                title="Iso",
                projection=Rhino.Display.DefinedViewportProjection.Top,
                position=System.Drawing.Rectangle(0, 0, 600, 600),
                floating=False)

        if rs.ViewProjection(view="Iso", mode=None) == 1:
            # 1 = parallel
            sLog_NoChange.append("View projection of Iso is already Parallel.")
        else:
            rs.ViewProjection(view="Iso", mode=1)
            if rs.ViewProjection(view="Iso", mode=None) == 1:
                sLog_Change.append("View projection of Iso now set to Parallel")
    
    if bSetAngleTol:
        fAngleTol_Deg_Old = rs.UnitAngleTolerance(angle_tolerance_degrees=None, in_model_units=True)
        if fAngleTol_Deg_Old == 0.1:
            sLog_NoChange.append("Model's angle tolerance is already 0.1{}.".format(chr(176)))
        else:
            rs.UnitAngleTolerance(angle_tolerance_degrees=0.1, in_model_units=True)
            fAngleTol_Deg_New = rs.UnitAngleTolerance(angle_tolerance_degrees=None, in_model_units=True)
            sLog_Change.append("Model's angle tolerance changed from {} to {}".format(
                    fAngleTol_Deg_Old, fAngleTol_Deg_New))
    
    if bSetPrec:
        if sc.doc.ModelDistanceDisplayPrecision == iPrecision:
            sLog_NoChange.append("Model distance display precision remains at {}.".format(
                    sc.doc.ModelDistanceDisplayPrecision))
        else:
            iPrec_Old = sc.doc.ModelDistanceDisplayPrecision
            sc.doc.ModelDistanceDisplayPrecision = iPrecision
            sLog_Change.append("Model distance display precision changed from {} to {}".format(
                    iPrec_Old,
                    sc.doc.ModelDistanceDisplayPrecision))
    
    # Layer modifications.
    
    sLayers = rs.LayerNames()
    
    if sLayers is None or len(sLayers) == 0:
        raise Exception("Where are the layers?")
    else:
        if bSetBlackLayersToWhite or bRenameLayer01ToDefault:
            if sc.doc.Layers.FindByFullPath("Default", Rhino.RhinoMath.UnsetIntIndex) >= 0:#rs.IsLayer('Default'):
                sLog_NoChange.append("Layer 'Default' already exists.")
            else:
                if bRenameLayer01ToDefault:
                    sLayerName_ToChange = "Layer 01"
                    sLayerName_Desired = "Default"
                    if sLayerName_ToChange in sLayers:
                        nLayers = rs.LayerCount()
                        if nLayers == 1:
                            rs.RenameLayer(sLayerName_ToChange, sLayerName_Desired)
                        else:
                            sLog_NoChange.append("{} layers exist, so {} not changed to {}.".format(
                                    nLayers,
                                    sLayerName_ToChange,
                                    sLayerName_Desired))
                
                if bAddDefaultLayer and not (sc.doc.Layers.FindByFullPath("Default", Rhino.RhinoMath.UnsetIntIndex) >= 0):#rs.IsLayer('Default'):
                    rs.AddLayer('Default')
                    if sc.doc.Layers.FindByFullPath("Default", Rhino.RhinoMath.UnsetIntIndex) >= 0:#rs.IsLayer('Default'):
                        rs.CurrentLayer('Default')
                        sLog_Change.append("'Default' layer added.")
        
        if bSetBlackLayersToWhite:
            iColors_changed_NonRef, iColors_changed_Ref = setLayerColors_BlackToWhite(bDebug)
            if iColors_changed_NonRef:
                sLog_Change.append("{} non-reference layers were changed from black to white.".format(iColors_changed_NonRef))
            else:
                sLog_NoChange.append("No non-reference layers were changed from black to white.")
            if iColors_changed_Ref:
                sLog_Change.append("{} reference layers were changed from black to white.".format(iColors_changed_Ref))
            else:
                sLog_NoChange.append("No reference layers were changed from black to white.")
    
    return sLog_Change, sLog_NoChange


def main(bEcho=None, bDebug=None):
    """
_NoEcho
_-Grid
_ApplyTo=AllViewports
_ShowGrid=No
_Enter
;
_-Runscript
(
Rhino.RenameView "Perspective", "Iso"
Rhino.ViewProjection "Iso", 1
)
;
;Set angle tolerance
_-DocumentProperties
_Units
_AngleTolerance=0.1
_EnterEnd
;
_Echo
_-Runscript
(
If Rhino.UnitSystem = 8 Then Rhino.Command "xModelDistanceDisplayPrecision 4", True
)
;
xSetBlackLayersToWhite
    """
    
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    if not getInput(): return
    
    if not bDebug: rs.EnableRedraw(False)
    
    sLog_Change, sLog_NoChange = run()
    
    rs.EnableRedraw()
    
    if not bEcho:
        return
    
    if sLog_NoChange and not sLog_Change:
        print("No changes to document.")
        return
    
    print("Not changed:")
    for s in sLog_NoChange:
        print("  " + s)
    print("Changed:")
    for s in sLog_Change:
        print("  " + s)


if __name__ == '__main__': main()