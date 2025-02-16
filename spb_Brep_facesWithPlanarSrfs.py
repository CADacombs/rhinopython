"""
180619: Created.
180630: Fixed bug in findFaces by converting complex conditionals into more simplers ones.
...
191215: Split from another script.
...
210221: Now will select all normal breps when input is Nothing.
210316, 220629, 0810, 0929: Modified an option default value.
230701: Fixed bug that didn't allow selection of planar faces closed in one direction.
250215: Modified an option default value.

TODO: Search only in any selected faces, not the entire brep.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import xBrepFace_tryGetPrimitiveShape
if Rhino.RhinoApp.ExeVersion == 5: import xBrepObject


sOpts = (
    'fGetShapeTol',
    'bEcho',
    'bDebug',
    )


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    

    for key in sOpts:
        keys.append(key)
        names[key] = key[1:] # Overwrite as wanted in the following.
        
    key = 'fGetShapeTol'
    #values[key] = 1e-9
    value = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    values[key] = float(format(value, '.0e'))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)
    
    key = 'bEcho'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if sc.sticky.has_key(stickyKeys[key]):
            if riOpts[key]:
                riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                #values[key] = sc.sticky[stickyKeys[key]]
                pass


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if cls.riOpts[key]:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                #sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                pass


def getInput():
    """Get open breps with optional input."""

    """
    #def getPreselectedBrepsIncludingSubObjects():

    #    # Check selected objects for breps.
    #    gBreps0_Preselected = []
    #    if sc.doc.Objects.GetSelectedObjectTypes() & rd.ObjectType.Brep:
    #        for rdRhinoObject in sc.doc.Objects.GetSelectedObjects(False,False):
    #            if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
    #                gBreps0_Preselected.append(rdRhinoObject.Id)

    #    # Check normal objects for subobject selections.


    #    iter = Rhino.DocObjects.ObjectEnumeratorSettings()
    #    iter.NormalObjects = True
    #    iter.LockedObjects = False
    #    iter.IncludeLights = False
    #    iter.IncludeGrips = False
    #    for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
    #        if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
    #            compIdxs = rdRhinoObject.GetSelectedSubObjects()
    #            if compIdxs is not None:
    #                # Add this brep to gBreps0_Preselected no matter which subobjects are selected.
    #                gBreps0_Preselected.append(rdRhinoObject.Id)
    #                rdRhinoObject.UnselectAllSubObjects()
    #                rdRhinoObject.UnhighlightAllSubObjects()
    #    return gBreps0_Preselected
    
    
    #gBreps0_Preselected = getPreselectedBrepsIncludingSubObjects()
    #if gBreps0_Preselected:
    #    for g in gBreps0_Preselected:
    #        sc.doc.Objects.Select(g)
    #    sc.doc.Views.Redraw()
    #    print "{} brep(s) are already selected.".format(len(gBreps0_Preselected))
    #else:
    #    print "No breps are selected."

    #sOptions = 'UnselectAll', 'UseOnlyPreselected', 'ContinueSelecting'

    #go = ri.Custom.GetOption()
    #go.SetCommandPrompt("Options")
    #if gBreps0_Preselected:
    #    go.SetCommandPromptDefault('PostSelectMore')
    #else:
    #    go.SetCommandPromptDefault('PostSelect')
    #go.AcceptNothing(True)
    
    #idxOpt_SetGeometryOpts = go.AddOption('SetGeometryOpts')
    #idxOpt_UsePreselected = go.AddOption('UseOnlyPreselected') if gBreps0_Preselected else None
    #if gBreps0_Preselected:
    #    idxOpt_PostSelect = go.AddOption('SelectMore')
    #else:
    #    idxOpt_PostSelect = go.AddOption('Select')
    

    #res = go.Get()
        
    #if res == ri.GetResult.Nothing:
    #    if gBreps0_Preselected:
    #        sOption = 'UseOnlyPreselected'
    #    else:
    #        sOption = 'AllNormal'
    #elif res == ri.GetResult.Cancel:
    #    return
    #else:
    #    if go.OptionIndex() == idxOpt_UsePreselected: sOption = 'UseOnlyPreselected'
    #    elif go.OptionIndex() == idxOpt_UsePreselected: sOption = 'AllNormalMonofaces'
    #    else: sOption = 'AllNormal'

    #print sOption

    #if sOption == 'AllNormal':
    #    pass
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False

    optIndices = {}

    while True:
        key = 'fGetShapeTol'
        go.AddOptionDouble(Opts.names[key], Opts.riOpts[key])
        idxOpt_SelNone = go.AddOption('UnselectAll')
        key = 'bEcho'
        go.AddOptionToggle(Opts.names[key], Opts.riOpts[key])
        key = 'bDebug'
        go.AddOptionToggle(Opts.names[key], Opts.riOpts[key])
    
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
    
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(
                    enable=False, ignoreUnacceptablePreselectedObjects=True)
            continue
    
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            
            gBreps = [objref.ObjectId for objref in go.Objects()]
            go.Dispose()
            return (
                    gBreps,
                    Opts.riOpts['fGetShapeTol'].CurrentValue,
                    Opts.riOpts['bEcho'].CurrentValue,
                    Opts.riOpts['bDebug'].CurrentValue,
            )

        if res == ri.GetResult.Nothing:
            iter = rd.ObjectEnumeratorSettings()
            iter.NormalObjects = True
            iter.LockedObjects = False
            iter.IncludeLights = False
            iter.IncludeGrips = False
            rdBrepObjects = []
            gBreps = []
            for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
                if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
                    rdBrepObjects.append(rdRhinoObject)
                    gBreps.append(rdRhinoObject.Id)
            go.Dispose()
            return (
                gBreps,
                Opts.riOpts['fGetShapeTol'].CurrentValue,
                Opts.riOpts['bEcho'].CurrentValue,
                Opts.riOpts['bDebug'].CurrentValue,
            )

        if res == ri.GetResult.Number:
            Opts.riOpts['fGetShapeTol'].CurrentValue = go.Number()
        elif go.OptionIndex() == idxOpt_SelNone:
            for objref in go.Objects():
                sc.doc.Objects.Select(objref, False)
            sc.doc.Objects.UnselectAll()
            sc.doc.Views.Redraw()
            return False
        
        key = 'fGetShapeTol'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        elif Opts.riOpts[key].CurrentValue < 2**(-52):
            Opts.riOpts[key].CurrentValue = 2**(-52) # Machine epsilon.
    
        Opts.saveSticky()
        go.ClearCommandOptions()


def getFaces(rgBrep0, idxs_rgFaces_ToSearch=None, fGetShapeTol=1e-9, bDebug=False):
    """
    """
    
    if not isinstance(rgBrep0, rg.Brep):
        print "rgBrep0 is not a Brep."
        return
    if not rgBrep0.IsValid:
        s  = "Brep is invalid and will be skipped.  "
        s += "_ExtractBadSrf and repair the bad faces first."
        print s
        return
    
    if idxs_rgFaces_ToSearch is None:
        idxs_rgFaces_ToSearch = range(rgBrep0.Faces.Count)
    
    idx_rgFaces_Found = []
    sShapes = []
    fTols_Used = []
    
    for idx_rgFace in idxs_rgFaces_ToSearch:
        rgFace = rgBrep0.Faces[idx_rgFace]
        
        if rgFace.IsClosed(0) and rgFace.IsClosed(1): continue
        
        rc = xBrepFace_tryGetPrimitiveShape.tryGetPrimitiveShape(
                rgFace0=rgFace,
                bPlane=True,
                bCylinder=False,
                bCone=False,
                bSphere=False,
                bTorus=False,
                fTolerance=fGetShapeTol,
                bDebug=bDebug
        )
        if rc[0]:
            shape, fTol_Used, sShrunkOrNot = rc[0]
            sShape = shape.GetType().Name

            if bDebug:
                print "{} found at tolerance {} from {} face."\
                        .format(sShape, fTol_Used, sShrunkOrNot)
            
            idx_rgFaces_Found.append(idx_rgFace)
            sShapes.append(sShape)
            fTols_Used.append(fTol_Used)
    
    rgBrep0.Dispose()
    
    return idx_rgFaces_Found, sShapes, fTols_Used


def processBrepObject(gBrep0, fGetShapeTol=None, bEcho=None, bDebug=None):
    """
    """
    
    if fGetShapeTol is None: fGetShapeTol = Opts.riOpts['fGetShapeTol'].CurrentValue
    if bEcho is None: bEcho = Opts.riOpts['bEcho'].CurrentValue
    if bDebug is None: bDebug = Opts.riOpts['bDebug'].CurrentValue
    
    rgBrep0 = rs.coercebrep(gBrep0)
    if rgBrep0 is None:
        print "rs.coercebrep returned None."
        return
    if not rgBrep0.IsValid:
        print "Brep {} is invalid and will be skipped.  " \
                "_ExtractBadSrf and repair the bad faces first.".format(
                gBrep0)
        return
    
    rc = getFaces(
            rgBrep0=rgBrep0,
            idxs_rgFaces_ToSearch=None,
            fGetShapeTol=fGetShapeTol,
            bDebug=bDebug,
    )
    if rc is None: return
    idx_rgFaces_Found, sShapes, fTols_Used = rc

    rgBrep0.Dispose()

    if bEcho:
        print "{} face(s) found in {}.".format(len(idx_rgFaces_Found), gBrep0)
        for sShape in set(sShapes):
            print "[{}] {}".format(sShapes.count(sShape), sShape),
        print

    if not idx_rgFaces_Found: return
    
    if rgBrep0.Faces.Count == 1:
        sc.doc.Objects.Select(gBrep0)
    else:
        if Rhino.RhinoApp.ExeVersion >= 6:
            rdBrep0 = rs.coercerhinoobject(gBrep0)
            for idx_rgFace in idx_rgFaces_Found:
                compIdx = rg.ComponentIndex(
                        rg.ComponentIndexType.BrepFace,
                        index=idx_rgFace)
                rdBrep0.SelectSubObject(
                        componentIndex=compIdx,
                        select=True, syncHighlight=True, persistentSelect=True)
        else:
            # Rhino V5.
            xBrepObject.extractFaces(
                    gBrep0,
                    idx_rgFaces_Found,
                    bCurrentLayer=False,
                    bByLayerColor=False,
                    bAddOnlyMonofaces=True,
                    bEcho=bEcho,
                    bDebug=bDebug)
    
    return rc


def processBrepObjects(gBreps0, fGetShapeTol=None, bEcho=None, bDebug=None):
    """Selects (extract for V5) faces and returns tuple of
    (idx_rgFaces_Found_PerBrep, fTols_Used_PerBrep)
    """
    
    if fGetShapeTol is None: fGetShapeTol = Opts.riOpts['fGetShapeTol'].CurrentValue
    if bEcho is None: bEcho = Opts.riOpts['bEcho'].CurrentValue
    if bDebug is None: bDebug = Opts.riOpts['bDebug'].CurrentValue
    
    gBreps_WithFacesFound = []
    idx_rgFaces_Found_PerBrep = []
    sShapes_All = []
    fTols_Used_PerBrep = []
    ct_BrepsWithFacesFound = 0
    ct_FacesFound_All = 0
    fTol_Used_Min = None
    fTol_Used_Max = None
    
    for iB, gBrep0 in enumerate(gBreps0):
        
        Rhino.RhinoApp.SetCommandPrompt(
                prompt="Searching brep {} of {} ...".format(iB+1, len(gBreps0)))
        
        rc = processBrepObject(
                gBrep0=gBrep0,
                fGetShapeTol=fGetShapeTol,
                bEcho=True if len(gBreps0) == 1 else False,
                bDebug=bDebug,
        )
        
        if rc is None:
            continue

        idx_rgFaces_Found, sShapes_1Brep, fTols_Used = rc
        
        gBreps_WithFacesFound.append(gBrep0)
        idx_rgFaces_Found_PerBrep.append(idx_rgFaces_Found)
        sShapes_All.extend(sShapes_1Brep)
        fTols_Used_PerBrep.append(fTols_Used)

        if idx_rgFaces_Found:
            ct_BrepsWithFacesFound += 1
        
        ct_FacesFound_All += len(idx_rgFaces_Found)

        if fTol_Used_Min is None:
            if fTols_Used:
                fTol_Used_Min = min(fTols_Used)
                fTol_Used_Max = max(fTols_Used)
        elif fTols_Used:
            fTol_Used_Min = min(fTol_Used_Min, min(fTols_Used))
            fTol_Used_Max = max(fTol_Used_Max, max(fTols_Used))

    if bEcho:
        if len(gBreps0) > 1:
            if ct_FacesFound_All == 1:
                print "1 face found in {} out of {} brep(s).".format(
                    ct_BrepsWithFacesFound, len(gBreps0))
            else:
                print "{} faces found in {} out of {} brep(s).".format(
                        ct_FacesFound_All, ct_BrepsWithFacesFound, len(gBreps0))
        
        if len(gBreps0) > 1:
            for sShape in set(sShapes_All):
                print "[{}] {}".format(sShapes_All.count(sShape), sShape),
            print
        
        if fTol_Used_Min is not None:
            print "Range of tolerances used: [{},{}]".\
                format(fTol_Used_Min, fTol_Used_Max)
    
        ct_Selected = len(list(sc.doc.Objects.GetSelectedObjects(False,False)))
        if ct_Selected: print "{} monoface breps are selected.".format(ct_Selected)

    return gBreps_WithFacesFound, idx_rgFaces_Found_PerBrep, fTols_Used_PerBrep


def main():
    
    rc = getInput()
    if rc is None: return

    if not Opts.riOpts['bDebug'].CurrentValue: sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll()
    
    processBrepObjects(rc[0])
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()