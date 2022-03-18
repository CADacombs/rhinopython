"""
170704-05: Created, starting with extractContiguousFacesOfDraftAngle.
171016: For faster extraction of faces, addBrepWithSubsetOfFaces was split into
        addBrepOfSubsetOfFaces_JoinBreps and addBrepOfSubsetOfFaces_RemoveAt.
...
180608: Now when all faces of brep are contiguous tangent, the brep is post-selected and layer and color are set to their relevant options.
180716: Now monoface breps are allowed as input.  This makes it easier when the brep's attributes are changed to current settings.
...
190808: Added bOnlyFirstAdj.
191010: Replaced an option with bRetainLayer and bRetainColor.
191101, 200701: Import-related update.
210429: Updated fAngleTol_Deg key for sc.sticky so the value will remain the same during this Rhino session.

WIP: bCopy option.

TODO: indicesOfContiguousTangentFaces: Check all edges of each face passing tangency so that no faces are marked as failing.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import xBrepObject

from System import Guid
from System.Collections.Generic import List
from System.Diagnostics import Stopwatch


sOpts = (
        'fAngleTol_Deg',
        'bOnlyFirstAdj',
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
    
    key = 'fAngleTol_Deg'
    values[key] = 50.0 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bOnlyFirstAdj'
    values[key] = False
    names[key] = 'ContiguousDepth'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Unlimited', 'One')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
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
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    
    @classmethod
    def setValues(cls):
        for key in cls.keys:
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
    """
    Get face with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select starting face")

    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(True, acceptZero=True)
    go.EnableHighlight(False)
    
    
    idxsOpts_Main = {}

    while True:
        go.AddOptionDouble(Opts.names['fAngleTol_Deg'], Opts.riOpts['fAngleTol_Deg'])
        go.AddOptionToggle(Opts.names['bOnlyFirstAdj'], Opts.riOpts['bOnlyFirstAdj'])
        go.AddOptionToggle(Opts.names['bRetainLayer'], Opts.riOpts['bRetainLayer'])
        go.AddOptionToggle(Opts.names['bRetainColor'], Opts.riOpts['bRetainColor'])
        go.AddOptionToggle(Opts.names['bCopy'], Opts.riOpts['bCopy'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])
            #rgFace = getBrepFaceFromGetObject(objref)
            #if rgFace is not None:
            #    go.Dispose()
            #    return tuple([objref, rgFace] + [Opts.values[key] for key in sOpts])
            #sc.doc.Objects.UnselectAll() # Prepare for repeat of go.Get().
            #sc.doc.Views.Redraw()
        elif res == ri.GetResult.Cancel:
            return # Esc key was pressed.
        
        # An option was selected or a number was entered.
        
        res = go.Result()
        
        key = 'fAngleTol_Deg'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def indicesOfContiguousTangentFaces(rgBrep, idxFace0, fAngleTol_Deg, bOnlyFirstAdj=False):
    
    idxFaces_Pass = [idxFace0]
    idxFace_LastAdded = [idxFace0]
    idxEdges_Pass = []
    idxEdges_Fail = []
    idxEdges_ToCheck = list(rgBrep.Faces[idxFace0].AdjacentEdges())
    
    for iE in idxEdges_ToCheck:
        if sc.escape_test(False):
            print "Script stopped in main edges to check loop."
            return
        
        if iE in idxEdges_Pass or iE in idxEdges_Fail:
            continue
        
        if rgBrep.Edges[iE].Valence != rg.EdgeAdjacency.Interior:
            idxEdges_Fail.append(iE)
            continue
        
        if rgBrep.Edges[iE].IsSmoothManifoldEdge(
                Rhino.RhinoMath.ToRadians(fAngleTol_Deg)):
            idxEdges_Pass.append(iE)
            for iF in rgBrep.Edges[iE].AdjacentFaces():
                if sc.escape_test(False):
                    print "Script stopped in adjacent face loop."
                    return
                
                if iF not in idxFaces_Pass:
                    idxFaces_Pass.append(iF)
                    if not bOnlyFirstAdj:
                        for iEPF in rgBrep.Faces[iF].AdjacentEdges():
                            if sc.escape_test(False):
                                print "Script stopped in adjacent edge loop."
                                return
                        
                            if iEPF not in (idxEdges_Pass + idxEdges_Fail):
                                idxEdges_ToCheck.append(iEPF)
        else:
            idxEdges_Fail.append(iE)
    
    return idxFaces_Pass


def main(bDebug=False):
    
    rc = getInput()
    if rc is None: return
    (
            objref,
            fAngleTol_Deg,
            bOnlyFirstAdj,
            bRetainLayer,
            bRetainColor,
            bCopy,
            bEcho,
            bDebug,
    ) = rc
    
    rdBrep0 = objref.Object()
    rgFace = objref.Face()
    
    stopwatch = Stopwatch()
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")
    
    rgBrep0 = rdBrep0.BrepGeometry
    
    idxFaces_Pass = indicesOfContiguousTangentFaces(
            rgBrep0,
            rgFace.FaceIndex,
            fAngleTol_Deg,
            bOnlyFirstAdj=bOnlyFirstAdj)
    if idxFaces_Pass is None: return
    
    if len(idxFaces_Pass) == 0:
        print "No faces found."
        return
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    if len(idxFaces_Pass) == rgBrep0.Faces.Count:
        print "All faces of brep are contiguous tangent."
    
    rgBrep0.Dispose()
    
    rc = xBrepObject.extractFaces(
            rdBrep0,
            idxFaces_Pass,
            bAddOnlyMonofaces=False,
            bRetainLayer=bRetainLayer,
            bRetainColor=bRetainColor,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    if rc is None:
        print "Faces could not be extracted."
        return

    gExtracted = rc[0]

    nSel = rs.SelectObjects(gExtracted)
    if nSel == 0:
        print "Error none of the breps are selected!"
    else:
        print "{} brep{} selected.".format(
                    nSel,
                    (' is', 's are')[bool(nSel-1)])
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main(bDebug=bool(0))