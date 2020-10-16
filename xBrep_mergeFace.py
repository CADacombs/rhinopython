"""
160709-24: Created.
160725: Now, faces are placed into either a pass or fail groups based on whether they can merge.
160728: Added shape simplification tolerance input option.  Updated debugging.
171211: File name changed from mergeFace2.py to mergeFace.py.
        Removed alias maker.
180131: Simplified getInput.
        Now, simplification tolerances used to find primitives are printed.
180204: Fixed bug in getInput.
180701: Fixed minor bug (incorrect variable name).
180823: Updated an import name.
181025: Added code from mergeToLargestFace.py and updated options accordingly.
181127: Fixed bug causing a crash when a number was directly entered in getInput.
190120: Some changes due to reference modifications, including adding an import reference.
190219: Added (some?) support for Toruses and NurbsSurfaces.
190221: Added (some?) support for Cones.  Updated an import name.
190226: Minor bug fix of incorrect parameter name.
190325: Updated some import names.
190330: Fixed bug that reported false positives for NurbsSurface matches.
        Fixed variable name bug.
190412-13: Added area checker.
        Now checks and tryies various seam locations for NurbsSurface until one is exterior of all the original brep faces.
        Updated some references.
190422: Now will add (not replace) breps with an area discrepancy in red.
190423: Replaced code to move seam for RevSurfaces from external to same code as for NurbsSurfaces.
190424: Moved merge-faces-to-largest into its own script.
190426: Improved seam-moving routine.
190427-28: Improved routine for spherical faces.
190502: Added support to new seam-moving routine for closed faces.
190513: Now Splits faces using joined curves since small segments can result in Split fail.
190519: Modified printed output.
190520: Fixed bug that was preventing non-NurbsSurfaces from being fully processed.
190521: Fixed bug that was sometimes selecting the incorrect split face for a sphere.
190522-24: Refactored.  Fixed bug in numerator loop.
190619: Corrected output for a particular fail.
190620: Modified some input interaction and printed output.
190720: Minor changes while tracking a bug in an import.
190722,24: Import-related changes.
190731: Modified processBrepObject's return data.
190810: Import-related changes.
191101: Import-related changes.  Bug fix.
191117: Non-manifold breps are no longer accepted.
191118-19, 200715: Import-related updates.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

import xBrep_contiguousCoshapedFaces
import xBrep_createFromShape
import xBrepObject
import xPrimitiveShape


sOpts = (
        'bMergeAll',
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
    
    key = 'bMergeAll'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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
        if sc.sticky.has_key(stickyKeys[key]):
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    
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


def getInput():
    """Get breps with optional input"""
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select starting face for merge")
    
    go.GeometryFilter = rd.ObjectType.Surface
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.SubSurface
    
    #go.EnableClearObjectsOnEntry(False)
    #go.AlreadySelectedObjectSelect = True
    #go.DeselectAllBeforePostSelect = False
    #go.EnableUnselectObjectsOnExit(False)
    #go.SubObjectSelect = False
    go.AcceptNumber(True, True)
    go.EnableHighlight(False)
    
    go.DisablePreSelect()

    print "Shape for merged brep will be obtained from first face picked."
    
    while True:
        go.AddOptionToggle(Opts.names['bMergeAll'], Opts.riOpts['bMergeAll'])
        go.AddOptionToggle(Opts.names['bTryConvert'], Opts.riOpts['bTryConvert'])
        go.AddOptionDouble(Opts.names['fTolerance'], Opts.riOpts['fTolerance'])
        go.AddOptionToggle(Opts.names['bShrinkSrfsOfPrimitives'], Opts.riOpts['bShrinkSrfsOfPrimitives'])
        go.AddOptionToggle(Opts.names['bExtract'], Opts.riOpts['bExtract'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return
        
        # An option was selected or a number was entered.
            
        key = 'fTolerance'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()
    
    if sc.doc.Objects.UnselectAll() > 0: sc.doc.Views.Redraw() # Allows for face highlighting.
    
    rdObjRef = go.Object(0)
    gBrep = rdObjRef.ObjectId
    rgBrep = rdObjRef.Brep()
    idx_rgFaceA = rdObjRef.GeometryComponentIndex.Index
    
    bMergeAll = Opts.values['bMergeAll']
    
    idxFaces_AdjA = rgBrep.Faces[idx_rgFaceA].AdjacentFaces()
        #idxFaces_AdjA = None
    
    idx_rgFaceB = None
    
    if not bMergeAll:
        if len(idxFaces_AdjA) == 1:
            # Face has only one adjacent face, so just use the face not picked.
            idx_rgFaceB = idxFaces_AdjA[0]
        elif len(idxFaces_AdjA) > 1:
            
            rdBrep = rdObjRef.Object()
            rdBrep.HighlightSubObject(rdObjRef.GeometryComponentIndex, True)
            sc.doc.Views.Redraw()
            
            #go.DisablePreSelect()
            go.SetCommandPrompt("Select a second face")
            
            def sameBrepDiffFaceGeomFilter (rhObject, geom, compIdx):
                return rhObject.Id == gBrep and compIdx.Index != idx_rgFaceA
            
            go.SetCustomGeometryFilter(sameBrepDiffFaceGeomFilter)
            
            while True:
                if sc.escape_test(False):
                    rdBrep.HighlightSubObject(rdObjRef.GeometryComponentIndex, False)
                    sc.doc.Views.Redraw()
                    return
                
                res = go.Get()
                
                if res == ri.GetResult.Object:
                    rdObjRef = go.Object(0)
                    idx_rgFaceB = rdObjRef.GeometryComponentIndex.Index
                    rdBrep.HighlightSubObject(rdObjRef.GeometryComponentIndex, False)
                    rdBrep.HighlightSubObject(rdObjRef.GeometryComponentIndex, False)
                    sc.doc.Views.Redraw()
                    break
                elif res == ri.GetResult.Cancel:
                    rdBrep.HighlightSubObject(rdObjRef.GeometryComponentIndex, False)
                    sc.doc.Views.Redraw()
                    return
                else:
                    # An option was selected or a number was entered.
                    Opts.processInput(go)
                    
                    # Check whether MergeAll was enabled during second face pick.
                    if Opts.values['bMergeAll']:
                        idx_rgFaceB = None
                        rdBrep.HighlightSubObject(rdObjRef.GeometryComponentIndex, False)
                        sc.doc.Views.Redraw()
                        break
    
    return (
            gBrep,
            idx_rgFaceA,
            idx_rgFaceB,
            Opts.values['bTryConvert'],
            Opts.values['fTolerance'],
            Opts.values['bShrinkSrfsOfPrimitives'],
            Opts.values['bExtract'],
            Opts.values['bEcho'],
            Opts.values['bDebug'],
    )


def processBrep(rgBrep0, idx_rgFaceA, idx_rgFaces_Filter=None, bTryConvert=None, fTolerance=None, bShrinkSrfsOfPrimitives=None, bEcho=None, bDebug=None):
    """
    """
    
    if bTryConvert is None: bTryConvert = bTryConvert
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bShrinkSrfsOfPrimitives is None: bShrinkSrfsOfPrimitives = Opts.values['bShrinkSrfsOfPrimitives']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    
    if not rgBrep0.IsManifold:
        if bEcho: print "Brep is non-manifold and will be skipped."
        return
    
    
    # Get shape or NurbsSurface of FaceA.
    rgFaceA = rgBrep0.Faces[idx_rgFaceA]
    rgPrimitiveShapeA = None
    rgSrfA = None
    if bTryConvert:
        rc = xPrimitiveShape.BrepFace.tryGetPrimitiveShape(
                rgFace0=rgFaceA,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc[0] is not None:
            rgPrimitiveShapeA, fTol_Used, sShrunkOrNot = rc[0]
            shapeA = rgPrimitiveShapeA
            idx_rgFs_B0_Merge = xBrep_contiguousCoshapedFaces.getPrimitiveShapedFaces(
                    rgFaceA=rgFaceA,
                    rgPrimitiveShapeA=rgPrimitiveShapeA,
                    idx_rgFaces_Filter=idx_rgFaces_Filter,
                    fTolerance=fTolerance,
                    bDebug=bDebug)

    if rgPrimitiveShapeA is None:
        rgSrfA = rgFaceA.UnderlyingSurface()
        if not isinstance(rgSrfA, rg.NurbsSurface):
            rgSrfA.Dispose()
            return
        shapeA = rgSrfA
        idx_rgFs_B0_Merge = xBrep_contiguousCoshapedFaces.getNurbsSrfFaces(
                rgFaceA=rgFaceA,
                rgNurbsSurfaceA=rgSrfA,
                idx_rgFaces_Filter=idx_rgFaces_Filter,
                fTolerance=fTolerance,
                bDebug=bDebug)
    
    if bDebug: sEval = 'idx_rgFs_B0_Merge'; print sEval + ':', eval(sEval)

    if len(idx_rgFs_B0_Merge) < 2:
        if bEcho: print "No matching contiguous faces found for starting face."
        return
    

    if bEcho:
        s  = "Merging to {} ...".format(shapeA.GetType().Name)
        Rhino.RhinoApp.SetCommandPrompt(prompt=s)
        if bDebug: print s
    
    rgBrep_SubOf0 = rgBrep0.DuplicateSubBrep(faceIndices=idx_rgFs_B0_Merge)
    
    area0 = rg.Brep.GetArea(rgBrep_SubOf0)
    if not area0:
        if bEcho:
            s  = "Area of starting brep {} cannot be calculated".format(gBrep0)
            s += " and will be skipped."
            print s
        return

    if bDebug: sEval = 'idx_rgFs_B0_Merge'; print sEval + ':', eval(sEval)


    rgBrep_New1F = xBrep_createFromShape.replaceShape(
            rgBrep0=rgBrep_SubOf0,
            shape=shapeA,
            fTolerance=None,
            bDebug=bDebug)
    if rgBrep_New1F is None:
        if bDebug: sEval = 'rgBrep_New1F'; print sEval + ':', eval(sEval)
        return


    # Check old vs. new areas.
    area1 = rg.Brep.GetArea(rgBrep_New1F)
    if not area1:
        if bEcho:
            s  = "Area of new brep to replace {}".format(gBrep0)
            s += " cannot be calculated, and this merge will be skipped."
            print s
        return
    percentAreaChange = 100.0 * abs(1.0 - area0 / area1)
    if percentAreaChange > 1.0:
    #if not Rhino.RhinoMath.EpsilonEquals(area0, area1, epsilon=1e-4):
        if bEcho:
            s  = "Change in area,"
            s += " {:.6f} of old vs. {:.6f} of new or {:.2g}%.".format(
                    area0, area1, percentAreaChange)
            print s
        if bDebug:
            attrRed = rd.ObjectAttributes()
            attrRed.ColorSource = rd.ObjectColorSource.ColorFromObject
            attrRed.ObjectColor = Color.Red
            gBrep1 = sc.doc.Objects.AddBrep(rgBrep_New1F, attrRed)
            if gBrep1 != Guid.Empty:
                print "New, red monoface brep added over old polysurface."
            else:
                print "New monoface brep could not be added for {}".format(gBrep0)
        return

    if bShrinkSrfsOfPrimitives and not isinstance(shapeA, rg.Surface):
        rgBrep_New1F.Faces.ShrinkFaces()
    
    return rgBrep_New1F, idx_rgFs_B0_Merge


def processBrepObject(gBrep0, idx_rgFaceA, idx_rgFaces_Filter=None, bTryConvert=None, fTolerance=None, bShrinkSrfsOfPrimitives=None, bExtract=None, bEcho=None, bDebug=None):
    """
    """
    
    if bTryConvert is None: bTryConvert = bTryConvert
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bShrinkSrfsOfPrimitives is None: bShrinkSrfsOfPrimitives = Opts.values['bShrinkSrfsOfPrimitives']
    if bExtract is None: bExtract = Opts.values['bExtract']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
    rgBrep0 = rdBrep0.Geometry
    
    rc = processBrep(
            rgBrep0=rgBrep0,
            idx_rgFaceA=idx_rgFaceA,
            idx_rgFaces_Filter=idx_rgFaces_Filter,
            bTryConvert=bTryConvert,
            fTolerance=fTolerance,
            bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    if rc is None:
        rgBrep0.Dispose()
        return
    
    rgBrep_New1F, idx_rgFs_B0_Merge = rc
    
    if bExtract: # Add new brep with same attributes of original and select.
        gBrep_Merged = sc.doc.Objects.AddBrep(rgBrep_New1F, rdBrep0.Attributes)
        rgBrep_New1F.Dispose()
        rgBrep0.Dispose()
        sc.doc.Objects.UnselectAll()
        bSelected = sc.doc.Objects.Select(objectId=gBrep_Merged, select=True)
        if bSelected and bEcho:
            print "{} faces were extracted, and the single, merged face is selected.".format(
                len(idx_rgFs_B0_Merge))
        
        # Delete or extract depending on whether all faces of original brep are changing.
        gBreps1_B0LessFaces = None
        if rgBrep0.Faces.Count == len(idx_rgFs_B0_Merge):
            # Delete original.
            sc.doc.Objects.Delete(rdBrep0, True)
        else:
            gBreps1_B0LessFaces = xBrepObject.removeFaces(rdBrep0, idx_rgFs_B0_Merge)
            if gBreps1_B0LessFaces is None and bEcho:
                print "Error in removing faces from brep."
                return
        rgBrep0.Dispose()
        
        return gBrep_Merged, gBreps1_B0LessFaces
    
    else: # Join / no extract
        # Check whether all faces of original are affected.
        if rgBrep0.Faces.Count == len(idx_rgFs_B0_Merge):
            sc.doc.Objects.Replace(gBrep0, rgBrep_New1F)
            rgBrep0.Dispose()
            if bDebug: print "{} faces were merged.".format(len(idx_rgFs_B0_Merge))
            return gBrep0
        else:
            if xBrepObject.replaceFaces(
                    rdBrep0,
                    idx_rgFs_B0_Merge,
                    [rgBrep_New1F],
                    bExtract=False,
                    fTolerance_Join=2.0*fTolerance,
                    bEcho=False,
                    bDebug=bDebug):
                rgBrep0.Dispose()
                if bDebug: print "{} faces were merged.".format(len(idx_rgFs_B0_Merge))
                return gBrep0
            elif bEcho:
                rgBrep_New1F.Dispose()
                print "Error in replacing faces of brep."


def main():
    
    rc = getInput()
    if rc is None: return
    
    (
            gBrep0,
            idx_rgFaceA,
            idx_rgFaceB,
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
    
    if bEcho: print "Maximum simplification tolerance allowed: {:.2e}".format(fTolerance)
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")

    idx_rgFaces_Filter = None if idx_rgFaceB is None else [idx_rgFaceB]
    
    processBrepObject(
            gBrep0=gBrep0,
            idx_rgFaceA=idx_rgFaceA,
            idx_rgFaces_Filter=idx_rgFaces_Filter,
            bTryConvert=bTryConvert,
            fTolerance=fTolerance,
            bShrinkSrfsOfPrimitives=bShrinkSrfsOfPrimitives,
            bExtract=bExtract,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()