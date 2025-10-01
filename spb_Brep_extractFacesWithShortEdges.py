"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160507: Created.
...
180628-30: Seams are no longer checked.  Various code refactoring.
...
210526: Added option quick sets.
210630: Added another option quick set.
220623: Expanded printed output.
220810: Fixed bug where text dots were not added on the current layer.
250924: Modified an option default value.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List
from System.Drawing import Color

import xBrepObject


sOpts = (
    'fLength_MaxPos',
    'bAddDot',
    'iDotHeight',
    'bExtract',
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
    
    key = 'fLength_MaxPos'
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'MaxEdgeLengthToSel'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bAddDot'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iDotHeight'
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bExtract'
    values[key] = True
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


    @classmethod
    def setValues(cls):
        for key in sOpts:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
    
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue


def getInput():
    """
    Get breps with optional input.
    """
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")
    
    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False
    
    go.AcceptNothing(True)

    go.AcceptNumber(True, acceptZero=True)
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    bPreselectedObjsChecked = False
    
    idxs_Opt = {}

    while True:
        go.AddOptionDouble(Opts.names['fLength_MaxPos'], Opts.riOpts['fLength_MaxPos'])
        key = 'TenthMT'; idxs_Opt[key] = go.AddOption(key)
        key = 'HalfMT'; idxs_Opt[key] = go.AddOption(key)
        key = 'MT'; idxs_Opt[key] = go.AddOption(key)
        key = 'DoubleMT'; idxs_Opt[key] = go.AddOption(key)
        key = 'TenMT'; idxs_Opt[key] = go.AddOption(key)
        go.AddOptionToggle(Opts.names['bAddDot'], Opts.riOpts['bAddDot'])
        if Opts.values['bAddDot']:
            go.AddOptionInteger(Opts.names['iDotHeight'], Opts.riOpts['iDotHeight'])
        #go.AddOptionToggle(Opts.names['bExtract'], Opts.riOpts['bExtract'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in sOpts])

        if res == ri.GetResult.Nothing:
            iter = rd.ObjectEnumeratorSettings()
            iter.NormalObjects = True
            iter.LockedObjects = False
            iter.IncludeLights = False
            iter.IncludeGrips = False
            rdBrepObjects = []
            for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
                if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
                    rdBrepObjects.append(rdRhinoObject)
            go.Dispose()
            return tuple([rdBrepObjects] + [Opts.values[key] for key in sOpts])

        if res == ri.GetResult.Number:
            Opts.riOpts['fLength_MaxPos'].CurrentValue = go.Number()

        elif go.Option().Index == idxs_Opt['TenthMT']:
            key = 'fLength_MaxPos'
            Opts.riOpts[key].CurrentValue = 0.1 * sc.doc.ModelAbsoluteTolerance

        elif go.Option().Index == idxs_Opt['HalfMT']:
            key = 'fLength_MaxPos'
            Opts.riOpts[key].CurrentValue = 0.5 * sc.doc.ModelAbsoluteTolerance

        elif go.Option().Index == idxs_Opt['MT']:
            key = 'fLength_MaxPos'
            Opts.riOpts[key].CurrentValue = sc.doc.ModelAbsoluteTolerance

        elif go.Option().Index == idxs_Opt['DoubleMT']:
            key = 'fLength_MaxPos'
            Opts.riOpts[key].CurrentValue = 2.0 * sc.doc.ModelAbsoluteTolerance

        elif go.Option().Index == idxs_Opt['TenMT']:
            key = 'fLength_MaxPos'
            Opts.riOpts[key].CurrentValue = 10.0 * sc.doc.ModelAbsoluteTolerance


        key = 'fLength_MaxPos'
        if Opts.riOpts[key].CurrentValue <= 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue


        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getFormattedDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        mddp = sc.doc.ModelDistanceDisplayPrecision
        p = mddp if mddp > 5 else 5
        return "{:.{}f}".format(fDistance, p)


def main():
    
    rc = getInput()
    if rc is None: return
    objrefs = rc[0]
    for key, value in zip(sOpts, rc[1:]):
        exec("{} = {}".format(key, value))
    
    iDecPlaces = sc.doc.ModelDistanceDisplayPrecision
    rgDots = []
    fLength_Min_Found = None; fLength_Max_Found = None
    gBreps1_All = [] # Accumulation of duplicated faces (breps)
    numShortEdges = 0
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    sc.doc.Objects.UnselectAll()
    
    rhBreps0 = objrefs

    len_gBreps0 = len(rhBreps0)
    idxs_AtTenths = [int(round(0.1*i*len(rhBreps0),0)) for i in range(10)]
    
    for iB, rhBrep0 in enumerate(rhBreps0):
        if sc.escape_test(False):
            print("Search was interrupted by user.")
            return
        
        if iB in idxs_AtTenths:
            Rhino.RhinoApp.SetCommandPrompt("Analyzed {:d}% of {} breps ...".format(
                int(100.0 * (iB+1) / len_gBreps0), len_gBreps0))
        
        # Obtain GUID, RhinoObject, and geometry.
        if isinstance(rhBrep0, Guid):
            gBrep0 = rhBrep0
            rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
            if rdBrep0 is None: continue
        elif isinstance(rhBrep0, rd.ObjRef):
            gBrep0 = rhBrep0.ObjectId
            rdBrep0 = rhBrep0.Object()
        elif isinstance(rhBrep0, rd.RhinoObject):
            rdBrep0 = rhBrep0
            gBrep0 = rhBrep0.Id
        else:
            print("GUID could not be obtained.")
            continue
        rgBrep0 = rdBrep0.Geometry
        if rgBrep0 is None:
            print("Brep geometry for {} cannot be obtained!".format(gBrep0))
            continue
        if not rgBrep0.IsValid and bEcho:
            print("Brep {} is invalid and will be skipped.".format(gBrep0))
            continue
        
        idx_rgFaceMatches = set()
        for rgEdge in rgBrep0.Edges:
            if rgEdge.Valence == rg.EdgeAdjacency.Interior:
                rgFaces_Adjacent = rgEdge.AdjacentFaces()
                if rgFaces_Adjacent[0] == rgFaces_Adjacent[1]:
                    continue
                # Alternative.
    #                if len(rgFaces_Adjacent) != len(set(rgFaces_Adjacent)):
    #                    continue
            fEdgeLength = rgEdge.GetLength()
            if fEdgeLength <= fLength_MaxPos:
                if ((fLength_Min_Found is None) or
                        (fEdgeLength < fLength_Min_Found)):
                            fLength_Min_Found = fEdgeLength
                if ((fLength_Max_Found is None) or
                        (fEdgeLength > fLength_Max_Found)):
                            fLength_Max_Found = fEdgeLength
                idx_rgFaces = set(rgEdge.AdjacentFaces())
                idx_rgFaceMatches |= idx_rgFaces
                if bAddDot:
                    sDot = '{0:.{1}f}'.format(fEdgeLength, iDecPlaces+1)
                    pts = rgEdge.DivideByCount(2, False)
                    if pts is not None: pt = rgEdge.PointAt(pts[0])
                    else: pt = rgEdge.PointAtStart
                    rgDot = rg.TextDot(sDot, pt)
                    rgDot.FontHeight = iDotHeight
                    rgDots.append(rgDot)
                numShortEdges += 1
        
        iFaceMatchCt = len(idx_rgFaceMatches)
        if iFaceMatchCt == 0: continue # No matches in this brep
        
        idx_rgFaceMatches = list(idx_rgFaceMatches)
        
        # Extract faces from polyface brep.
        rc = xBrepObject.extractFaces(
                gBrep0,
                idx_rgFaceMatches,
                bCurrentLayer=False,
                bByLayerColor=False,
                bAddOnlyMonofaces=True,
                bEcho=False,
                bDebug=bDebug)
        if rc is None: continue

        gBreps1_All.extend(rc[0])
        
        rgBrep0.Dispose()
        
        # End of rdBreps0 loop.
    
    if bDebug: sPrint = 'gBreps1_All'; print(sPrint + ':', eval(sPrint))

    if len(rgDots) > 0:
        attr = rd.ObjectAttributes()
        attr.LayerIndex = sc.doc.ActiveDoc.Layers.CurrentLayerIndex
        attr.ColorSource = rd.ObjectColorSource.ColorFromObject
        attr.ObjectColor = Color.Red
        for rgDot in rgDots:
            gDot = sc.doc.Objects.AddTextDot(rgDot, attr)
            sc.doc.Objects.Select(gDot)


    if numShortEdges == 0:
        print("No short edges found.")
        sc.doc.Views.RedrawEnabled = True
        return

    print("{} short edge(s) found.".format(numShortEdges))
    print("{} monoface breps are selected.".format(
            sc.doc.Objects.Select(List[Guid](gBreps1_All))))

    if numShortEdges == 1:
        print("Length of single short edge: {}".format(
            getFormattedDistance(fLength_Min_Found)))
    else:
        print("Range of short edge lengths: [{}, {}]".format(
            getFormattedDistance(fLength_Min_Found),
            getFormattedDistance(fLength_Max_Found)))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()