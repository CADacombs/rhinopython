"""
190523-24: Created.
190721: Removed an unused import.
190722: Import-related change.
191004: Bug fix.
191101: Import-related update.
191117: Import-related bug fix.
191118: Import-related update.
241224: Bug fix.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import xBrepObject
import xPrimitiveShape


sOpts = (
        'bPrimitives',
        'bNurbsSrfs',
        'fTolerance',
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
    
    key = 'bPrimitives'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bNurbsSrfs'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTolerance'
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
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
    go.SetCommandPrompt("Select starting face for match")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.SubSurface
    
    go.EnableHighlight(False)
    
    go.DisablePreSelect()
    
    go.AcceptNumber(enable=True, acceptZero=True)

    print "Base shape will be obtained from first face picked."
    
    while True:
        go.AddOptionToggle(Opts.names['bPrimitives'], Opts.riOpts['bPrimitives'])
        go.AddOptionToggle(Opts.names['bNurbsSrfs'], Opts.riOpts['bNurbsSrfs'])
        go.AddOptionDouble(Opts.names['fTolerance'], Opts.riOpts['fTolerance'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return
        
        key = 'fTolerance'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = abs(go.Number())
        if Opts.riOpts[key].CurrentValue == 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()
    
    if sc.doc.Objects.UnselectAll() > 0: sc.doc.Views.Redraw() # Allows for face highlighting.
    
    rdObjRef = go.Object(0)
    gBrep = rdObjRef.ObjectId
    rgBrep = rdObjRef.Brep()
    idx_rgFaceA = rdObjRef.GeometryComponentIndex.Index
    
    idxFaces_AdjA = rgBrep.Faces[idx_rgFaceA].AdjacentFaces()
        #idxFaces_AdjA = None
    
    return (
            gBrep,
            idx_rgFaceA,
            Opts.values['bPrimitives'],
            Opts.values['bNurbsSrfs'],
            Opts.values['fTolerance'],
            Opts.values['bEcho'],
            Opts.values['bDebug'],
    )


def getNurbsSrfFaces(rgFaceA, rgNurbsSurfaceA=None, idx_rgFaces_Filter=None, fTolerance=1e-9, bDebug=False):
    
    if rgNurbsSurfaceA is None:
        rgSrfA = rgFaceA.UnderlyingSurface()
        if not isinstance(rgSrfA, rg.NurbsSurface):
            rgSrfA.Dispose()
            return
        nsA = rgSrfA
    else:
        nsA = rgNurbsSurfaceA

    idx_rgFaceA = rgFaceA.FaceIndex
    
    # Find matching shapes of adjacent faces.
    idxFaces_Pass = [idx_rgFaceA]
    idxFaces_Fail = []
    idxFaces_LastAdded = [idx_rgFaceA]
    
    rgBrep0 = rgFaceA.Brep

    for idxFace_LastAdded in idxFaces_LastAdded:
        sc.escape_test()
        
        for idx_Face_Adj in rgBrep0.Faces[idxFace_LastAdded].AdjacentFaces():
            if idx_rgFaces_Filter and idx_Face_Adj not in idx_rgFaces_Filter:
                continue
            elif idx_Face_Adj in idxFaces_Pass:
                continue
            elif idx_Face_Adj in idxFaces_Fail:
                continue
                
            rgFaceB = rgBrep0.Faces[idx_Face_Adj]
                
            srfB = rgFaceB.UnderlyingSurface()
            if not isinstance(srfB, rg.NurbsSurface): continue
                
            nsB = srfB
            if not nsA.EpsilonEquals(other=nsB, epsilon=fTolerance):
                idxFaces_Fail.append(idx_Face_Adj)
                continue
            # Main parameters of NurbsSurface are equal, so now check all ControlPoints.
            if nsA.Points.CountU != nsB.Points.CountU:
                idxFaces_Fail.append(idx_Face_Adj)
                continue
            if nsA.Points.CountV != nsB.Points.CountV:
                idxFaces_Fail.append(idx_Face_Adj)
                continue
            for ptA, ptB in zip(nsA.Points, nsB.Points):
                if not ptA.EpsilonEquals(other=ptB, epsilon=fTolerance):
                    idxFaces_Fail.append(idx_Face_Adj)
                    break
            else:
                if bDebug: print "NurbsSurfaces are EpsilonEqual."
                idxFaces_Pass.append(idx_Face_Adj)
                idxFaces_LastAdded.append(idx_Face_Adj)
    return idxFaces_Pass


def getPrimitiveShapedFaces(rgFaceA, rgPrimitiveShapeA=None, idx_rgFaces_Filter=None, fTolerance=1e-9, bDebug=False):
        
    rgBrep0 = rgFaceA.Brep
    idx_rgFaceA = rgFaceA.FaceIndex

    if rgPrimitiveShapeA is None:
        rc = xPrimitiveShape.BrepFace.tryGetPrimitiveShape(
                rgFace0=rgFaceA,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc[0] is None: return
        rgPrimitiveShapeA = rc[0][0]
    
    bPlane = bCylinder = bCone = bSphere = bTorus = False
    
    if isinstance(rgPrimitiveShapeA, rg.Plane):
        bPlane = True
    elif isinstance(rgPrimitiveShapeA, rg.Cylinder):
        bCylinder = True
    elif isinstance(rgPrimitiveShapeA, rg.Cone):
        bCone = True
    elif isinstance(rgPrimitiveShapeA, rg.Sphere):
        bSphere = True
    elif isinstance(rgPrimitiveShapeA, rg.Torus):
        bTorus = True
    
    idxFaces_Pass = [idx_rgFaceA]
    
    # Find matching shapes of adjacent faces.
    idxFaces_Pass = [idx_rgFaceA]
    idxFaces_Fail = []
    idxFaces_LastAdded = [idx_rgFaceA]
    
    for idxFace_LastAdded in idxFaces_LastAdded:
        sc.escape_test()
        
        for idx_Face_Adj in rgBrep0.Faces[idxFace_LastAdded].AdjacentFaces():
            if idx_rgFaces_Filter and idx_Face_Adj not in idx_rgFaces_Filter:
                continue
            elif idx_Face_Adj in idxFaces_Pass:
                continue
            elif idx_Face_Adj in idxFaces_Fail:
                continue
            
            rgFaceB = rgBrep0.Faces[idx_Face_Adj]
            
            rc = xPrimitiveShape.BrepFace.tryGetPrimitiveShape(
                    rgFace0=rgFaceB,
                    bPlane=bPlane,
                    bCylinder=bCylinder,
                    bCone=bCone,
                    bSphere=bSphere,
                    bTorus=bTorus,
                    fTolerance=fTolerance,
                    bDebug=bDebug)
            if rc[0] is None:
                idxFaces_Fail.append(idx_Face_Adj)
                continue
                    
            shapeB = rc[0][0]

            if xPrimitiveShape.AnyShape.areEqual(
                    [rgPrimitiveShapeA, shapeB],
                    epsilon=fTolerance,
                    bDebug=bDebug
            ):
                idxFaces_Pass.append(idx_Face_Adj)
                idxFaces_LastAdded.append(idx_Face_Adj)
            else:
                idxFaces_Fail.append(idx_Face_Adj)
    return idxFaces_Pass


def getFaces(rgBrep0, idx_rgFaceA, idx_rgFaces_Filter=None, bMergeAll=True, bPrimitives=None, bNurbsSrfs=None, fTolerance=None, bEcho=None, bDebug=None):
    """
    """
    
    if bPrimitives is None: bPrimitives = Opts.values['bPrimitives']
    if bNurbsSrfs is None: bNurbsSrfs = Opts.values['bNurbsSrfs']
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    # Get shape of FaceA.
    rgShapeA = None
    
    rgFaceA = rgBrep0.Faces[idx_rgFaceA]
    
    if bPrimitives:
        rc = xPrimitiveShape.BrepFace.tryGetPrimitiveShape(
                rgFace0=rgFaceA,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc is not None and rc[0] is not None:
            rgShapeA, fTol_Used, sShrunkOrNot = rc[0]
            idx_rgFs_B0_Pass = getPrimitiveShapedFaces(
                    rgFaceA=rgFaceA,
                    rgPrimitiveShapeA=rgShapeA,
                    idx_rgFaces_Filter=idx_rgFaces_Filter,
                    fTolerance=fTolerance,
                    bDebug=bDebug)

    if bNurbsSrfs and rgShapeA is None:
        rgShapeA = rgFaceA.UnderlyingSurface()
        if not isinstance(rgShapeA, rg.NurbsSurface):
            rgShapeA.Dispose()
            return
        idx_rgFs_B0_Pass = getNurbsSrfFaces(
                rgFaceA=rgFaceA,
                rgNurbsSurfaceA=rgShapeA,
                idx_rgFaces_Filter=idx_rgFaces_Filter,
                fTolerance=fTolerance,
                bDebug=bDebug)

    if bDebug: sEval='idx_rgFs_B0_Pass'; print sEval+':',eval(sEval)

    if len(idx_rgFs_B0_Pass) < 2:
        if bEcho: print "No matching contiguous faces found for starting face."
        return
    
    return idx_rgFs_B0_Pass, rgShapeA


def processBrepObject(gBrep0, idx_rgFaceA, idx_rgFaces_Filter=None, bPrimitives=None, bNurbsSrfs=None, fTolerance=None, bCopy=None, bEcho=None, bDebug=None):
    """
    """
    
    if bPrimitives is None: bPrimitives = Opts.values['bPrimitives']
    if bNurbsSrfs is None: bNurbsSrfs = Opts.values['bNurbsSrfs']
    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
    rgBrep0 = rdBrep0.Geometry
    
    rc = getFaces(
            rgBrep0=rgBrep0,
            bPrimitives=bPrimitives,
            bNurbsSrfs=bNurbsSrfs,
            idx_rgFaceA=idx_rgFaceA,
            idx_rgFaces_Filter=idx_rgFaces_Filter,
            bMergeAll=True,
            fTolerance=fTolerance,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    if rc is None: return
    idx_rgFs_B0_Pass, rgShapeA = rc
    
    if rgBrep0.Faces.Count == len(idx_rgFs_B0_Pass):
        print "Entire brep matches the face."
        bSelected = rs.SelectObject(gBrep0)
        rgBrep0.Dispose()
        return bSelected
    else:
        if Rhino.RhinoApp.ExeVersion >= 6:
            xBrepObject.selectFaces(rdBrep0, idx_rgFs_B0_Pass)
            rgBrep0.Dispose()
        else:
            # Rhino V5.
            rc = xBrepObject.extractFaces(
                    gBrep0,
                    idxFaces=idx_rgFs_B0_Pass,
                    bCurrentLayer=False,
                    bByLayerColor=False,
                    bAddOnlyMonofaces=True,
                    bEcho=bEcho,
                    bDebug=bDebug)
            rgBrep0.Dispose()
            return None if rc is None else rc[0]


def main():
    
    rc = getInput()
    if rc is None: return
    
    (
            gBrep0,
            idx_rgFaceA,
            bPrimitives,
            bNurbsSrfs,
            fTolerance,
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
    
    if bEcho: print "Maximum conversion tolerance allowed: {}".format(fTolerance)
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")

    processBrepObject(
            gBrep0=gBrep0,
            bPrimitives=bPrimitives,
            bNurbsSrfs=bNurbsSrfs,
            idx_rgFaceA=idx_rgFaceA,
            idx_rgFaces_Filter=None,
            fTolerance=fTolerance,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()