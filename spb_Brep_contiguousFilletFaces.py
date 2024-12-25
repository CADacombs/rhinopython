"""
160412-16: Created.
160601: Modularizations.
160611: Bug fixes.
160823: Various minor changes.
161022: Commented out alias maker.
170701: Script name changed from ActOnChainFilletFaces to extractChainFilletFaces.
        Removed code that allowed deletion and duplication of faces.
        Added command option for destination layer of extracted faces.  All attributes, not just the layer, are actually applied.
180305: Now, extracted faces are joined instead of being monoface breps.
        Name changed from extractChainFilletFaces.py to extractContiguousFilletFaces.py
181118: Added Opts and more options.
181204: Updated an import name.
190128: Refactored and added option for fillet concavity matching.
        Renamed from extractContiguousFilletFaces.py to Brep_ContiguousFilletFaces.py.
190129: Updated an import name.
190130: Added SetCommandPrompt.
190207, ..., 0605: Updated an import name.
190620: Added Opts and changed getInput accordingly.  Changed the display name of an option for clarity.
190810: Changed an option default.
191010,1101: Import-related update.
191118: Moved a function from another module to this script.
191205: Bug fix.
"""

import Rhino
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import xBrepObject
import xSurface_radius


sOpts = (
        'bSameConcavityType',
        'bSameRad',
        'fRadTol',
        'bTanContinuity',
        'fAngleTol',
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
    
    for key in sOpts:
        keys.append(key)
        names[key] = key[1:] # Overwrite as wanted in the following.
    
    key = 'bSameConcavityType'
    values[key] = True
    names[key] = 'ConcavityType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Either', 'Same')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSameRad'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fRadTol'
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bTanContinuity'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fAngleTol'
    values[key] = 50.0 * sc.doc.ModelAngleToleranceDegrees
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bRetainLayer'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bRetainColor'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCopy'
    values[key] = False
    names[key] = 'Copy'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'
    values[key] = True
    names[key] = 'Echo'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
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


def getInput(bDebug=False):
    """
    """
    
    # Get constant radius (rolling ball) fillet face of brep.
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select starting face")
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface
    go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.SubSurface) # Doesn't allow single surfaces to be selected.
    
    go.AcceptNumber(enable=True, acceptZero=True)
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == ri.GetResult.Object 
    
    while True:
        go.AddOptionToggle(Opts.names['bSameConcavityType'], Opts.riOpts['bSameConcavityType'])
        go.AddOptionToggle(Opts.names['bSameRad'], Opts.riOpts['bSameRad'])
        go.AddOptionDouble(Opts.names['fRadTol'], Opts.riOpts['fRadTol'])
        go.AddOptionToggle(Opts.names['bTanContinuity'], Opts.riOpts['bTanContinuity'])
        if Opts.values['bTanContinuity']:
            go.AddOptionDouble(Opts.names['fAngleTol'], Opts.riOpts['fAngleTol'])
        go.AddOptionToggle(Opts.names['bRetainLayer'], Opts.riOpts['bRetainLayer'])
        go.AddOptionToggle(Opts.names['bRetainColor'], Opts.riOpts['bRetainColor'])
        go.AddOptionToggle(Opts.names['bCopy'], Opts.riOpts['bCopy'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fRadTol'].CurrentValue = go.Number()
        
        if Opts.riOpts['fRadTol'].CurrentValue < 0.0:
            Opts.riOpts['fRadTol'].CurrentValue = Opts.riOpts['fRadTol'].InitialValue
        elif Opts.riOpts['fAngleTol'].CurrentValue < 0.0:
            Opts.riOpts['fAngleTol'].CurrentValue = Opts.riOpts['fAngleTol'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()
    
    bDebug = Opts.values['bDebug']
    
    objref = go.Object(0)
    
    rgFace = objref.Face()
    
    fPlanarTol = sc.doc.ModelAbsoluteTolerance
    if rgFace.IsPlanar(fPlanarTol):
        print "Planar face selected."
        return
    
    fFaceRad0 = xSurface_radius.constantRadiusOfSurface(rgFace, Opts.values['fRadTol'], bEcho=bDebug)
    if bDebug: sPrint = 'fFaceRad0'; print sPrint + ':', eval(sPrint)
    if fFaceRad0 is None:
        print "Non-fillet face selected."
        return
    
    rdBrep0 = objref.Object() # Parent
    
    idxFace_Sel = objref.GeometryComponentIndex.Index
    
    gBrep0 = objref.ObjectId
    
    return (
            gBrep0,
            idxFace_Sel,
            Opts.values['bSameConcavityType'],
            Opts.values['bSameRad'],
            Opts.values['fRadTol'],
            Opts.values['bTanContinuity'],
            Opts.values['fAngleTol'],
            Opts.values['bRetainLayer'],
            Opts.values['bRetainColor'],
            Opts.values['bCopy'],
            Opts.values['bEcho'],
            Opts.values['bDebug'],
    )


def curvatures(rgFace):
    """
    Returns: MajorCurvature, MinorCurvature as floats.
    Sign of curvature is per face, which may not be the same as the underlying surface.
    """
    areaMassProp = Rhino.Geometry.AreaMassProperties.Compute(rgFace)
    if areaMassProp is None:
        print "Face {} is skipped because" \
              " its centroid cannot be calculated.".format(
              rgFace.FaceIndex)
        return
    
    ptCentrdW = areaMassProp.Centroid
    
    getrc, u, v = rgFace.ClosestPoint(ptCentrdW)
    
    c = rgFace.CurvatureAt(u, v)
    if c is None: return
    
    if rgFace.OrientationIsReversed:
        return -c.Kappa(0), -c.Kappa(1)
    else:
        return c.Kappa(0), c.Kappa(1)


def indicesOfContiguousFilletFaces(rgBrep0, idxFace0, bSameConcavityType=True, bSameRad=False, fRadTol=1e-6, bTanContinuity=False, fAngleTol=5.0):
    
    fPlanarTol = 0.1 * sc.doc.ModelAbsoluteTolerance
    
    fAngleTol_Rad = Rhino.RhinoMath.ToRadians(fAngleTol)
    
    rgFace0 = rgBrep0.Faces[idxFace0]
    
    fFaceRad0 = xSurface_radius.constantRadiusOfSurface(rgFace0, fRadTol)
    if fFaceRad0 is None: return
    
    if bSameConcavityType:
        rc = curvatures(rgFace0)
        if not rc: return
        kappa_ToMatch = rc[0]
    
    idxFaces_Pass = [idxFace0]
    idxFaces_Fail = []
    idxFace_LastAdded = [idxFace0]
    
    for i in idxFace_LastAdded:
        sc.escape_test()
        
        for j in rgBrep0.Faces[i].AdjacentFaces():
            #print 'j:', j
            if not j in idxFaces_Pass and not j in idxFaces_Fail:
                rgFaceX = rgBrep0.Faces[j]
                
                if rgFaceX.IsPlanar(fPlanarTol):
                    idxFaces_Fail.append(j)
                    continue
                
                fRadFace = xSurface_radius.constantRadiusOfSurface(rgFaceX, fRadTol)
                if fRadFace is None:
                    idxFaces_Fail.append(j)
                    continue
                
                if bSameConcavityType:
                    rc = curvatures(rgFaceX)
                    if not rc:
                        idxFaces_Fail.append(j)
                        continue
                    kappa_Max, kappa_Min = rc
                    if kappa_Max == 0.0 and kappa_Min == 0.0:
                        idxFaces_Fail.append(j)
                        continue
                    elif kappa_Max == 0.0:
                        if (
                                (abs(1.0/kappa_Min - 1.0/kappa_ToMatch) > fRadTol)
                        ):
                            idxFaces_Fail.append(j)
                            continue
                    elif kappa_Min == 0.0:
                        if (
                                (abs(1.0/kappa_Max - 1.0/kappa_ToMatch) > fRadTol)
                        ):
                            idxFaces_Fail.append(j)
                            continue
                    else:
                        if (
                                (abs(1.0/kappa_Max - 1.0/kappa_ToMatch) > fRadTol)
                                and
                                (abs(1.0/kappa_Min - 1.0/kappa_ToMatch) > fRadTol)
                        ):
                            idxFaces_Fail.append(j)
                            continue
                
                if bSameRad and not abs(fRadFace - fFaceRad0) <= fRadTol:
                    idxFaces_Fail.append(j)
                    continue
                
                if bTanContinuity:
                    idxEdgesA = list(rgBrep0.Faces[i].AdjacentEdges())
                    idxEdgesB = list(rgFaceX.AdjacentEdges())
                    edges_Shared = list(set(idxEdgesA) & set(idxEdgesB))
                    bTanFound = False
                    for edge_Shared in edges_Shared:
                        if rgBrep0.Edges[edge_Shared].IsSmoothManifoldEdge(
                                fAngleTol_Rad
                        ):
                            bTanFound = True
                    if not bTanFound:
                        idxFaces_Fail.append(j)
                        continue
                
                idxFaces_Pass.append(j)
                idxFace_LastAdded.append(j)
    return idxFaces_Pass


def main(bEcho=True, bDebug=False):
    
    rc = getInput()
    if rc is None: return
    (
            gBrep0,
            idxFace_Sel,
            bSameConcavityType,
            bSameRad,
            fRadTol,
            bTanContinuity,
            fAngleTol,
            bRetainLayer,
            bRetainColor,
            bCopy,
            bEcho,
            bDebug,
    ) = rc
    if idxFace_Sel is None: return
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")
    
    rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
    rgBrep0 = rdBrep0.BrepGeometry
    
    idxFaces_Pass = indicesOfContiguousFilletFaces(
            rgBrep0,
            idxFace_Sel,
            Opts.values['bSameConcavityType'],
            Opts.values['bSameRad'],
            Opts.values['fRadTol'],
            Opts.values['bTanContinuity'],
            Opts.values['fAngleTol'],
    )
    
    rc = xBrepObject.extractFaces(
            gBrep0,
            idxFaces_Pass,
            bAddOnlyMonofaces=False,
            bRetainLayer=bRetainLayer,
            bRetainColor=bRetainColor,
   )
    if rc is None:
        print "Faces could not be extracted."
        return
    else:
        gExtracted = rc[0]
    
    #    # Add brep(s) of faces to be extracted.
    #    
    #    gExtracted = xBrepObject.addFromSubsetOfFaces(
    #            gBrep0, idxFaces_Pass, bRetainLayer,
    #            bEcho=True, bDebug=False,
    #            )
    #    if gExtracted is None:
    #        print "Faces could not be extracted."
    #        return
    #    
    #    # Remove faces from brep.
    #    
    #    rc = xBrepObject.removeFaces(rdBrep0, idxFaces_Pass)
    #    
    #    if rc is None:
    #        print "xBrepObject.removeFaces returned None!"
    #        return
    
    if not gExtracted:
        print "No face of target radius found."
        return
    
    ct_Face = 0
    for gBrep1 in gExtracted:
        rgBrep1 = rs.coercebrep(gBrep1)
        ct_Face += rgBrep1.Faces.Count
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    ct_Breps_Selected = rs.SelectObjects(gExtracted)
    
    print "{} brep{} selected with {} face{}.".format(
            ct_Breps_Selected,
            " is" if ct_Breps_Selected == 1 else "s are",
            ct_Face,
            "" if ct_Face == 1 else "s",
    )
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main(bEcho=bool(1), bDebug=bool(0))
