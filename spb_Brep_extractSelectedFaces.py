"""
180612: Created, starting with extractFilletFacesByRadius.py.
181109: Added Opts.  Added bRetainColor option.
181127: Enabled bCopy functionality.
181202, 190325: Updated an import name.
190330: Fixed bug concerning single extracted face not being selected.
190429: Added function to handle preselected faces.
191010: Import-related update.
191022: Feedback change.
191101: Import-related update.
250604: Modified a printed feedback.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import xBrepObject


sOpts = (
        'bAddOnlyMonofaces',
        'bRetainLayer',
        'bRetainColor',
        'bCopy',
        'bEcho',
        'bDebug',
)


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    
    @classmethod
    def init(cls):
        
        for key in sOpts:
            cls.keys.append(key)
            cls.names[key] = key[1:] # Overwrite as wanted in the following.
            
        key = 'bAddOnlyMonofaces'
        cls.keys.append(key)
        cls.values[key] = False
        cls.names[key] = key[1:]
        cls.riOpts[key] = ri.Custom.OptionToggle(cls.values[key], 'No', 'Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bRetainLayer'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = key[1:]
        cls.riOpts[key] = ri.Custom.OptionToggle(cls.values[key], 'No', 'Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bRetainColor'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = key[1:]
        cls.riOpts[key] = ri.Custom.OptionToggle(cls.values[key], 'No', 'Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bCopy'
        cls.keys.append(key)
        cls.values[key] = False
        cls.names[key] = 'Copy'
        cls.riOpts[key] = ri.Custom.OptionToggle(cls.values[key], 'No', 'Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bEcho'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = 'Echo'
        cls.riOpts[key] = ri.Custom.OptionToggle(cls.values[key], 'No', 'Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bDebug'
        cls.keys.append(key)
        cls.values[key] = False
        cls.names[key] = 'Debug'
        cls.riOpts[key] = ri.Custom.OptionToggle(cls.values[key], 'No', 'Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        # Load sticky.
        for key in cls.stickyKeys:
            if sc.sticky.has_key(cls.stickyKeys[key]):
                if cls.riOpts[key]:
                    cls.values[key] = cls.riOpts[key].CurrentValue = sc.sticky[cls.stickyKeys[key]]
                else:
                    # For OptionList.
                    cls.values[key] = sc.sticky[cls.stickyKeys[key]]
    
    @classmethod
    def addOptions(cls, go, sInclude=None):
        if sInclude is None:
            sInclude = cls.keys[:]
        if not sInclude: return
        
        def addOptionDouble(key, lower=None, upper=None):
            if key in sInclude:
                go.AddOptionDouble(cls.names[key], cls.riOpts[key])
        
        def addOptionInteger(key):
            if key in sInclude:
                go.AddOptionInteger(cls.names[key], cls.riOpts[key])
        
        def addOptionToggle(key):
            if key in sInclude:
                go.AddOptionToggle(cls.names[key], cls.riOpts[key])
        
        addOptionToggle('bAddOnlyMonofaces')
        addOptionToggle('bRetainLayer')
        addOptionToggle('bRetainColor')
        addOptionToggle('bCopy')
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
Opts.init()


def getPreselectedFaceInput():
    """
    """
    
    go = ri.Custom.GetObject()
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface
    go.EnablePreSelect(enable=True, ignoreUnacceptablePreselectedObjects=False)
    go.EnablePostSelect(enable=False)
    res = go.GetMultiple(1,0)
    
    if res != ri.GetResult.Object: return
    
    gBreps = []
    idx_Faces_AllBreps = []
    
    # Prepare lists of brep GUIDs and face indices.
    for objref in go.Objects():
        gBrep = objref.ObjectId
        if gBrep in gBreps:
            idx_Faces_AllBreps[gBreps.index(gBrep)].append(objref.Face().FaceIndex)
        else:
            gBreps.append(gBrep)
            idx_Faces_AllBreps.append([objref.Face().FaceIndex])
            rdBrep = sc.doc.Objects.Find(gBrep)
    
    sc.doc.Objects.UnselectAll()
    
    gopt = ri.Custom.GetOption()
    gopt.SetCommandPrompt("Set options for preselected faces")
    
    gopt.AcceptNothing(enable=True)
    
    Opts.addOptions(gopt)
    
    while True:
        res = gopt.Get()
        
        if res == ri.GetResult.Nothing:
            break
        elif res == ri.GetResult.Cancel:
            return
        else:
            # An option was selected or a number was entered.
            Opts.processInput(gopt)
    
    return (
            gBreps,
            idx_Faces_AllBreps,
            Opts.values['bAddOnlyMonofaces'],
            Opts.values['bRetainLayer'],
            Opts.values['bRetainColor'],
            Opts.values['bCopy'],
            Opts.values['bEcho'],
            Opts.values['bDebug'])


def getInput():
    """
    """
    
    # Get breps with optional input
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select faces to extract")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    #go.EnableUnselectObjectsOnExit(False)
    
    #sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()
    
    Opts.addOptions(go)
    
    while True:
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return
        else:
            # An option was selected or a number was entered.
            Opts.processInput(go)
        
        # Rehighlight faces already selected.
        for objref in go.Objects():
            objref.Object().HighlightSubObject(
                    componentIndex=objref.GeometryComponentIndex,
                    highlight=True)
        sc.doc.Views.Redraw()
    
    gBreps = []
    idx_Faces_AllBreps = []
    
    # Prepare lists of brep GUIDs and face indices.
    for objref in go.Objects():
        gBrep = objref.ObjectId
        if gBrep in gBreps:
            idx_Faces_AllBreps[gBreps.index(gBrep)].append(objref.Face().FaceIndex)
        else:
            gBreps.append(gBrep)
            idx_Faces_AllBreps.append([objref.Face().FaceIndex])
            rdBrep = sc.doc.Objects.Find(gBrep)
            #rdBrep.Highlight(enable=False)
    
    # Unselect any non-surface preselections.
    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()
    
    return (
            gBreps,
            idx_Faces_AllBreps,
            Opts.values['bAddOnlyMonofaces'],
            Opts.values['bRetainLayer'],
            Opts.values['bRetainColor'],
            Opts.values['bCopy'],
            Opts.values['bEcho'],
            Opts.values['bDebug'])


def processBrepObjects(gBreps0, idx_Faces_AllBreps, bAddOnlyMonofaces=None, bRetainLayer=None, bRetainColor=None, bCopy=None, bEcho=None, bDebug=None):
    """
    """
    
    if bAddOnlyMonofaces is None: bAddOnlyMonofaces = Opts.values['bAddOnlyMonofaces']
    if bRetainLayer is None: bRetainLayer = Opts.values['bRetainLayer']
    if bRetainColor is None: bRetainColor = Opts.values['bRetainColor']
    if bCopy is None: bCopy = Opts.values['bCopy']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    try:
        gBreps0 = list(gBreps0)
        try:
            dummy = list(idx_Faces_AllBreps[0])
            idx_Faces_AllBreps = list(idx_Faces_AllBreps)
        except:
            idx_Faces_AllBreps = [idx_Faces_AllBreps]
    except:
        gBreps0 = [gBreps0]
        idx_Faces_AllBreps = [idx_Faces_AllBreps] # Based on an assumption.  (Not fail-safe.)
    
    gExtracted_All = [] # Accumulation of duplicated faces (breps)
    
    for gBrep0, idx_rgFaces in zip(gBreps0, idx_Faces_AllBreps):
        if bCopy:
            gExtracted = xBrepObject.addFromSubsetOfFaces(
                    gBrep0,
                    idx_rgFaces,
                    bAddOnlyMonofaces=bAddOnlyMonofaces,
                    bRetainLayer=bRetainLayer,
                    bRetainColor=bRetainColor,
            )
            if gExtracted is None:
                print "Faces could not be duplicated."
                gExtracted = []
        else:
            # Extract faces.
            rc = xBrepObject.extractFaces(
                    gBrep0,
                    idx_rgFaces,
                    bAddOnlyMonofaces=bAddOnlyMonofaces,
                    bRetainLayer=bRetainLayer,
                    bRetainColor=bRetainColor,
            )
            if rc is None:
                print "Faces could not be extracted."
                gExtracted = []
            else: gExtracted = rc[0]
                
        
        if gExtracted is not None:
            gExtracted_All.extend(gExtracted)
    
    count_Selected = sc.doc.Objects.Select(List[Guid](gExtracted_All))
    
    if not count_Selected:
        print "Face(s) could not be extracted.  Is the polyfaced brep invalid?"
    else:
        print "{} brep(s) are selected.".format(count_Selected)
    
    return gExtracted_All


def main(bEcho=True, bDebug=False):
    
    rc = getPreselectedFaceInput() if Rhino.RhinoApp.ExeVersion >= 6 else None
    if rc is not None:
        (
                gBreps0,
                idx_Faces_AllBreps,
                bAddOnlyMonofaces,
                bRetainLayer,
                bRetainColor,
                bCopy,
                bEcho,
                bDebug) = rc
    else:
        rc = getInput()
        if rc is None: return
        (
                gBreps0,
                idx_Faces_AllBreps,
                bAddOnlyMonofaces,
                bRetainLayer,
                bRetainColor,
                bCopy,
                bEcho,
                bDebug) = rc
    
    if bDebug:
        reload(xBrepObject)
    
    processBrepObjects(
            gBreps0=gBreps0,
            idx_Faces_AllBreps=idx_Faces_AllBreps)
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main(bEcho=bool(1), bDebug=bool(0))