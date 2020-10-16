"""
190424: Started this script as an extract from another.
190426: Improved object selection behavior in getInput.
190428: Minor bug fix.
190502: Improved geometry selection filter.  Renamed an import.
190520: Added bProcessAllFaces.  Refactored Opts and getInput.
190524: Import-related updates.
190619: Refactored 1 function into 2.
190719: Minor change while tracking a bug.
190731: Bug fix for when main brep is divided into multiple breps when its faces are extracted.
191031: Import-related update.  Added a function.
191119: Import-related update.
200401: Refactored processBreps.  Simplified printed output.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

import xBrep_mergeFace


sOpts = (
        'bProcessAllFaces',
        'bTryConvert',
        'fTolerance',
        'bShrinkSrfsOfPrimitives',
        'bExtract',
        'bEcho',
        'bDebug',
)


class Opts():
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    
    for key in sOpts:
        keys.append(key)
        names[key] = key[1:] # Overwrite as wanted in the following.
    
    key = 'bProcessAllFaces'
    values[key] = True
    names[key] = "Process"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'OnlyToLargestFace', 'AllFaces')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bTryConvert'
    values[key] = True
    names[key] = 'TryConvertToPrimitive'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTolerance'
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bShrinkSrfsOfPrimitives'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bExtract'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
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
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    
    @classmethod
    def setValues(cls):
        for key in sOpts:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                pass
    
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """Get breps with optional input"""
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select polyface breps to merge")
    
    go.GeometryFilter = rd.ObjectType.PolysrfFilter
    go.SubObjectSelect = False
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    go.AcceptNumber(True, True)
    
    bPreselectedObjsChecked = False
    
    while True:
        go.AddOptionToggle(Opts.names['bProcessAllFaces'], Opts.riOpts['bProcessAllFaces'])
        go.AddOptionToggle(Opts.names['bTryConvert'], Opts.riOpts['bTryConvert'])
        go.AddOptionDouble(Opts.names['fTolerance'], Opts.riOpts['fTolerance'])
        go.AddOptionToggle(Opts.names['bShrinkSrfsOfPrimitives'], Opts.riOpts['bShrinkSrfsOfPrimitives'])
        go.AddOptionToggle(Opts.names['bExtract'], Opts.riOpts['bExtract'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return
        else:
            # An option was selected or a number was entered.
            key = 'fTolerance'
            if res == ri.GetResult.Number:
                Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 0.0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
            Opts.setValues()
            Opts.saveSticky()
            #go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            go.ClearCommandOptions()

    gBreps = [objref.ObjectId for objref in go.Objects()]
    
    go.Dispose()
    
    return (
            gBreps,
            Opts.values['bProcessAllFaces'],
            Opts.values['bTryConvert'],
            Opts.values['fTolerance'],
            Opts.values['bShrinkSrfsOfPrimitives'],
            Opts.values['bExtract'],
            Opts.values['bEcho'],
            Opts.values['bDebug'],
    )


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.Brep):
        return rhObj
    elif isinstance(rhObj, rg.Surface):
        return rhObj.ToBrep()
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

    if isinstance(geom, rg.Brep):
        return geom


def processAllFacesOfBrepObject(gBrep0, bTryConvert=None, fTolerance=None, bShrinkSrfsOfPrimitives=None, bExtract=None, bEcho=None, bDebug=None):
    """
    """
    
    if bTryConvert is None: bTryConvert = bTryConvert
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bShrinkSrfsOfPrimitives is None: bShrinkSrfsOfPrimitives = Opts.values['bShrinkSrfsOfPrimitives']
    if bExtract is None: bExtract = Opts.values['bExtract']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    # Process all remaining faces as Face A.
    
    rgBrep0 = coerceBrep(gBrep0)
    ct_Face_Start = ct_Face_Current = rgBrep0.Faces.Count
    ct_Face_End = None
    
    iCountFromLargest = 0
    fLastAreaForA = None
    
    s = "Calculating face areas for {} faces ...".format(
            ct_Face_Current)
    Rhino.RhinoApp.SetCommandPrompt(s)
    
    idx_Faces = []
    fAreas = []
    for rgFace in rgBrep0.Faces:
        rgBrep_1Face = rgFace.DuplicateFace(duplicateMeshes=True)
        fAreas.append(rgBrep_1Face.GetArea())
        idx_Faces.append(rgFace.FaceIndex)
        rgBrep_1Face.Dispose()
    
    zipped = zip(fAreas, idx_Faces)
    zipped.sort(reverse=True)
    
    rgBrep0.Dispose()
    
    while True:
        sc.escape_test()
        
        # Using areas as a measure skip most faces already determine.
        # Most because faces of same areas may be repeated to ensure none are skipped.
        fAreaFaceA = zipped[iCountFromLargest][0]
        if fLastAreaForA is not None:
            while fAreaFaceA > fLastAreaForA:
                iCountFromLargest += 1
                if iCountFromLargest == len(zipped) - 1:
                    break
                fAreaFaceA = zipped[iCountFromLargest][0]
        
        if iCountFromLargest == len(zipped) - 1:
            break
        
        idx_rgFaceA = zipped[iCountFromLargest][1]
        fLastAreaForA = fAreaFaceA
        
        s = "Processing face {} of {} ...".format(
                iCountFromLargest+1, len(zipped))
        Rhino.RhinoApp.SetCommandPrompt(s)
        if bDebug: print s + "  FaceIndex: " + str(idx_rgFaceA)
        
        rc = xBrep_mergeFace.processBrepObject(
                gBrep0=gBrep0,
                idx_rgFaceA=idx_rgFaceA,
                idx_rgFaces_Filter=None,
                bTryConvert=bTryConvert,
                fTolerance=fTolerance,
                bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
                bExtract=bExtract,
                bEcho=False,
                bDebug=bDebug,
        )
        
        if rc is None:
            iCountFromLargest += 1
            if iCountFromLargest >= len(zipped) - 1:
                break
        else:
            if bExtract:
                gBrep_Merged, gBreps1_B0LessFaces = rc
                if not gBreps1_B0LessFaces:
                    ct_Face_End = 1 # Let's see if this always makes sense.
                    break
                elif len(gBreps1_B0LessFaces) > 1:
                    print "Original Brep has been divided into {} breps.".format(
                            len(gBreps1_B0LessFaces))
                    print "Rerun this script on all new polyface breps."
                    return
            else:
                pass
            
            iCountFromLargest = 0
            rgBrep0 = coerceBrep(gBrep0)
            if rgBrep0.Faces.Count == 1:
                rgBrep0.Dispose()
                break
            
            ct_Face_Current = rgBrep0.Faces.Count
            
            s = "Calculating face areas for {} remaining faces ...".format(
                    ct_Face_Current)
            Rhino.RhinoApp.SetCommandPrompt(s)
            idx_Faces = []
            fAreas = []
            for rgFace in rgBrep0.Faces:
                rgBrep_1Face = rgFace.DuplicateFace(duplicateMeshes=True)
                fAreas.append(rgBrep_1Face.GetArea())
                idx_Faces.append(rgFace.FaceIndex)
                rgBrep_1Face.Dispose()
            
            zipped = zip(fAreas, idx_Faces)
            zipped.sort(reverse=True)
            
            rgBrep0.Dispose()
    
    if ct_Face_End is None:
        rgBrep0 = coerceBrep(gBrep0)
        ct_Face_End = rgBrep0.Faces.Count
        rgBrep0.Dispose()
    
    if ct_Face_Start != ct_Face_End:
        if bDebug:
            print "Face count changed from {} to {}.".format(
                ct_Face_Start, ct_Face_End)
    else:
        if bDebug:
            print "No faces were merged.  Face count remains at {}.".format(
                ct_Face_Start)
    
    return ct_Face_Start - ct_Face_End


def processOnlyToLargestFaceOfBrepObject(gBrep0, bTryConvert=None, fTolerance=None, bShrinkSrfsOfPrimitives=None, bExtract=None, bEcho=None, bDebug=None):
    """
    """
    
    if bTryConvert is None: bTryConvert = bTryConvert
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bShrinkSrfsOfPrimitives is None: bShrinkSrfsOfPrimitives = Opts.values['bShrinkSrfsOfPrimitives']
    if bExtract is None: bExtract = Opts.values['bExtract']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    rgBrep0 = coerceBrep(gBrep0)
    
    # Find largest face in brep.
    fAreas = []
    for rgFace in rgBrep0.Faces:
        rgBrep_1Face = rgFace.DuplicateFace(duplicateMeshes=True)
        fAreas.append(rgBrep_1Face.GetArea())
        rgBrep_1Face.Dispose()
    
    idx_rgFaceA = fAreas.index(max(fAreas))

    rc = xBrep_mergeFace.processBrepObject(
            gBrep0=gBrep0,
            idx_rgFaceA=idx_rgFaceA,
            idx_rgFaces_Filter=None,
            bTryConvert=bTryConvert,
            fTolerance=fTolerance,
            bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
            bExtract=bExtract,
            bEcho=bEcho,
            bDebug=bDebug,
    )

    return rc


def processBreps(gBreps0, bProcessAllFaces=None, bTryConvert=None, fTolerance=None, bShrinkSrfsOfPrimitives=None, bExtract=None, bEcho=None, bDebug=None):
    """
    """
    
    if bProcessAllFaces is None: bProcessAllFaces = Opts.values['bProcessAllFaces']
    if bTryConvert is None: bTryConvert = Opts.values['bTryConvert']
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bShrinkSrfsOfPrimitives is None: bShrinkSrfsOfPrimitives = Opts.values['bShrinkSrfsOfPrimitives']
    if bExtract is None: bExtract = Opts.values['bExtract']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    if bEcho:
        print "Maximum simplification tolerance allowed: {}".format(fTolerance)
    
    
    if bProcessAllFaces:
        iCt_Breps_NoMerge = 0
        iCt_Breps_Merge = 0
        iCt_Faces_Total = 0
        
        for gBrep0 in gBreps0:
            iCt_Faces = processAllFacesOfBrepObject(
                    gBrep0=gBrep0,
                    bTryConvert=bTryConvert,
                    fTolerance=fTolerance,
                    bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
                    bExtract=bExtract,
                    bEcho=bEcho,
                    bDebug=bDebug,
            )
            if isinstance(iCt_Faces, int):
                if iCt_Faces:
                    iCt_Breps_Merge += 1
                    iCt_Faces_Total += iCt_Faces
                else:
                    iCt_Breps_NoMerge += 1
        if bEcho:
            s = "{} breps with merged faces.".format(iCt_Breps_Merge)
            if iCt_Faces_Total:
                s += "  Total face count was reduced by {}.".format(iCt_Faces_Total)
            if iCt_Breps_NoMerge:
                s += "  No faces were merged in {} breps.".format(iCt_Breps_NoMerge)
            print s
    else:
        gBreps1 = []
        fTols_used_for_primitive = []
        
        for gBrep0 in gBreps0:
            rc = processOnlyToLargestFaceOfBrepObject(
                    gBrep0=gBrep0,
                    bTryConvert=bTryConvert,
                    fTolerance=fTolerance,
                    bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
                    bExtract=bExtract,
                    bEcho=bEcho,
                    bDebug=bDebug,
            )
            if rc is None: continue
            gBrep0, fTol_PrimitiveUsedA = rc
            gBreps1.append(gBrep0)
            if fTol_PrimitiveUsedA is not None:
                fTols_used_for_primitive.append(rc[1])
        
        if bEcho and gBreps1:
            ct_gBreps0 = len(gBreps1)
            if ct_gBreps0 == 1:
                print "1 brep was merged."
            else:
                print "{} breps were merged.".format(ct_gBreps0)
            if fTols_used_for_primitive:
                print "Range of simplification tolerances used to find largest faces: {} to {}".format(
                        min(fTols_used_for_primitive),
                        max(fTols_used_for_primitive))


def main():
    
    rc = getInput()
    if rc is None: return
    
    (
            gBreps0,
            bProcessAllFaces,
            bTryConvert,
            fTolerance,
            bShrinkSrfsOfPrimitives,
            bExtract,
            bEcho,
            bDebug,
    ) = rc
    
    if bDebug:
        import sys
        for sModule in list(sys.modules):
            if sModule[0] == 'x':
                try:
                    reload(sys.modules[sModule])
                    print "{} reloaded.".format(sys.modules[sModule])
                except:
                    # This is for any module name changes.
                    print "{} NOT reloaded.".format(sys.modules[sModule])
    else:
        sc.doc.Views.RedrawEnabled = False
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")
    
    processBreps(
            gBreps0=gBreps0,
            bProcessAllFaces=bProcessAllFaces,
            bTryConvert=bTryConvert,
            fTolerance=fTolerance,
            bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
            bExtract=bExtract,
            bEcho=bEcho,
            bDebug=bDebug,
    )
        
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()