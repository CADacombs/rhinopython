"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
171209: Created.
...
200202: Import-related update.
200507: Updated for V6+ using new BrepEdgeList.MergeAllEdges method.  Import-related update.
210221: Corrected printed feedback.
210426: Import-related update.
210630: Modified an option default value.
250222: Import-related update.
250225: Modified an option default value.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import xBrep_compact
import spb_Brep_rebuildEdges


sOpts = (
    'fAngleTol',
    'bRebuildEdges',
    'bRebuildOnlyOutOfTol',
    'fRebuildEdgesTol',
    'bRepeat',
    'bEcho',
    'bDebug',
    )


class Opts():
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}


    key = 'fAngleTol'; keys.append(key)
    values[key] = min((0.1, sc.doc.ModelAngleToleranceDegrees))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bRebuildEdges'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bRebuildOnlyOutOfTol'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fRebuildEdgesTol'; keys.append(key)
    values[key] = 0.5 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bRepeat'; keys.append(key)
    values[key] = True
    names[key] = 'RepeatMerge'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)


    for key in keys:
        if key not in names:
            names[key] = key[1:]


    # Load sticky.
    for key in stickyKeys:
        if sc.sticky.has_key(stickyKeys[key]):
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if cls.riOpts[key]:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select breps")
    
    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False
    
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    go.AcceptNumber(True, True)
    
    bPreselectedObjsChecked = False
    
    while True:
        go.AddOptionDouble(Opts.names['fAngleTol'], Opts.riOpts['fAngleTol'])
        go.AddOptionToggle(Opts.names['bRebuildEdges'], Opts.riOpts['bRebuildEdges'])
        if Opts.riOpts['bRebuildEdges'].CurrentValue:
            go.AddOptionToggle(Opts.names['bRebuildOnlyOutOfTol'], Opts.riOpts['bRebuildOnlyOutOfTol'])
            go.AddOptionDouble(Opts.names['fRebuildEdgesTol'], Opts.riOpts['fRebuildEdgesTol'])
            go.AddOptionToggle(Opts.names['bRepeat'], Opts.riOpts['bRepeat'])
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
        elif res == ri.GetResult.Number:
            Opts.riOpts['fAngleTol'].CurrentValue = abs(go.Number())
        
        for key in 'fRebuildEdgesTol', 'fAngleTol':
            if Opts.riOpts[key].CurrentValue < 0.0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        for key in Opts.stickyKeys:
            if Opts.riOpts[key]:
                Opts.values[key] = Opts.riOpts[key].CurrentValue
            else:
                # For OptionList.
                pass
        
        Opts.saveSticky()
        
        go.ClearCommandOptions()
    
    gBreps = [objref.ObjectId for objref in go.Objects()]
    go.Dispose()
    
    return (
            gBreps,
            Opts.values['fAngleTol'],
            Opts.values['bRebuildEdges'],
            Opts.values['bRebuildOnlyOutOfTol'],
            Opts.values['fRebuildEdgesTol'],
            Opts.values['bRepeat'],
            Opts.values['bEcho'],
            Opts.values['bDebug'],
    )


def edgeCount_All(gBreps):
    if isinstance(gBreps, Guid):
        gBreps = [gBreps]
    count = 0
    for gBrep in gBreps:
        rgBrep = rs.coercebrep(gBrep)
        count += rgBrep.Edges.Count
        rgBrep.Dispose()
    return count


def edgeCount_Naked(gBreps):
    if isinstance(gBreps, Guid):
        gBreps = [gBreps]
    count = 0
    for gBrep in gBreps:
        rgBrep = rs.coercebrep(gBrep)
        for rgEdge in rgBrep.Edges:
            if rgEdge.Valence == rg.EdgeAdjacency.Naked:
                count += 1
        rgBrep.Dispose()
    return count


def mergeEdges_V5(gBreps, fAngleTol=None):
    """
    Returns: Boolean based on whether any edges were merged.
    """
    if fAngleTol is None: fAngleTol = Opts.values['fAngleTol']
    
    cts_Edges0 = []
    for gBrep in gBreps:
        if sc.escape_test(throw_exception=False, reset=False):
            print "Break in edges count."
            break
        cts_Edges0.append(edgeCount_All(gBrep))
    
    sc.doc.Objects.UnselectAll()
    for gBrep in gBreps: sc.doc.Objects.Select(gBrep)
    
    Rhino.RhinoApp.SetCommandPrompt("Merging ...")
    
    fAngleTol0 = sc.doc.ModelAngleToleranceDegrees
    sc.doc.ModelAngleToleranceDegrees = fAngleTol
    
    rs.Command("_MergeAllEdges", echo=False)
    
    sc.doc.ModelAngleToleranceDegrees = fAngleTol0
    
    sc.doc.Objects.UnselectAll()
    
    cts_Edges1 = []
    for gBrep in gBreps:
        if sc.escape_test(throw_exception=False, reset=False):
            print "Break in edges count."
            break
        cts_Edges1.append(edgeCount_All(gBrep))
    
    gBreps_WithEdgesMerged = []
    for gBrep, ct_Edges0, ct_Edges1 in zip(gBreps, cts_Edges0, cts_Edges1):
        if ct_Edges1 < ct_Edges0:
            gBreps_WithEdgesMerged.append(gBrep)
    
    return gBreps_WithEdgesMerged


def mergeEdges(gBreps, fAngleTol=None):
    """
    Returns: Boolean based on whether any edges were merged.
    """
    if fAngleTol is None: fAngleTol = Opts.values['fAngleTol']
    
    fAngleTol_Rad = Rhino.RhinoMath.ToRadians(fAngleTol)
    
    Rhino.RhinoApp.SetCommandPrompt("Merging ...")
    
    gBreps_WithEdgesMerged = []
    
    for gBrep in gBreps:
        if sc.escape_test(throw_exception=False, reset=False):
            print "Break in edges count."
            break
        
        rgB_In = rs.coercebrep(gBrep)
        
        rgB_Out = rgB_In.Duplicate()
        
        iCt_Edges_beforeMergeAll = rgB_In.Edges.Count
        
        rgB_Out.Edges.MergeAllEdges(angleTolerance=fAngleTol_Rad)
        
        iCt_Edges_afterMergeAll = rgB_Out.Edges.Count
        
        if iCt_Edges_beforeMergeAll == iCt_Edges_afterMergeAll:
            rgB_Out.Dispose()
            continue
        
        sc.doc.Objects.Replace(objectId=gBrep, brep=rgB_Out)
        
        gBreps_WithEdgesMerged.append(gBrep)
    
    return gBreps_WithEdgesMerged


def createEdgeCountReport(nEdges0, nNakedEdges0, nEdges1, nNakedEdges1, nBreps, sPrefix=''):
    s  = "{}Edge count reduced by {} ({}->{})".format(
            sPrefix, nEdges0 - nEdges1, nEdges0, nEdges1)
    s += "  ({}Naked edge count reduced by {} ({}->{}))".format(
            sPrefix, nNakedEdges0 - nNakedEdges1, nNakedEdges0, nNakedEdges1)
    s += " in {} brep(s).".format(nBreps)
    return s


def processBrepObjects(gBreps, fAngleTol=None, bRebuildEdges=None, bRebuildOnlyOutOfTol=None, fRebuildEdgesTol=None, bRepeat=None, bEcho=None, bDebug=None):
    """For faster processing, all brep objects will be modified by _MergeAllEdges at once.
    They will also be sent in plural to spb_Brep_rebuildEdges.processBreps .
    """
    
    if fAngleTol is None: fAngleTol = Opts.values['fAngleTol']
    if bRebuildEdges is None: bRebuildEdges = Opts.values['bRebuildEdges']
    if bRebuildOnlyOutOfTol is None: bRebuildOnlyOutOfTol = Opts.values['bRebuildOnlyOutOfTol']
    if fRebuildEdgesTol is None: fRebuildEdgesTol = Opts.values['fRebuildEdgesTol']
    if bRepeat is None: bRepeat = Opts.values['bRepeat']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    iPrec = sc.doc.ModelDistanceDisplayPrecision
    s = "Merge angle tolerance: {0:.{1}f} degrees.".format(
            fAngleTol, (12, iPrec)[0.1**iPrec<fAngleTol])
            # Output at 12 decimal places if tolerance is < precision display decimals.
    if bRebuildEdges:
        s += "  Rebuild edges: {0:.{1}f}.".format(
                fRebuildEdgesTol, (12, iPrec)[0.1**iPrec<fRebuildEdgesTol])
                # Output at 12 decimal places if tolerance is < precision display decimals.
        s += "  RebuildOnlyOutOfTol: {}".format(bRebuildOnlyOutOfTol)
    else:
        s += "  Edges will not be rebuilt."
    print s
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Recording starting edge counts ...")
    nEdges0 = edgeCount_All(gBreps)
    nNakedEdges0 = edgeCount_Naked(gBreps)
    
    iStep = 1
    
    if bRebuildEdges and bRepeat:
        print "STEP {}: 1st merging of edges".format(iStep)
    #else:
    #    s = "Merging edges"
    
    gBreps_WithEdgesMerged = mergeEdges(gBreps, fAngleTol)

    if not gBreps_WithEdgesMerged:
        print "None of the {} edges / {} naked edges were merged.".format(
                nEdges0, nNakedEdges0)
        nEdges1 = nEdges0
        nNakedEdges1 = nNakedEdges0
    else:
        Rhino.RhinoApp.SetCommandPrompt(prompt="Recording edge counts ...")
        nEdges1 = edgeCount_All(gBreps)
        nNakedEdges1 = edgeCount_Naked(gBreps)
        
        print createEdgeCountReport(nEdges0, nNakedEdges0, nEdges1, nNakedEdges1, len(gBreps))
    
    if not bRebuildEdges:
        return
    
    iStep = 2
    
    if bRepeat:
        print "STEP {}: 1st rebuilding of edges".format(iStep)
    else:
        print "Rebuilding edges"
    Rhino.RhinoApp.SetCommandPrompt(prompt=s+' ...')
    
    sc.doc.Objects.UnselectAll()
    gBreps_EdgesRebuilt = spb_Brep_rebuildEdges.processBrepObjects(
            gBreps,
            fRebuildNakedTol=fRebuildEdgesTol,
            fSearchTol=2.0*fRebuildEdgesTol,
            bRebuildOnlyOutOfTol=bRebuildOnlyOutOfTol,
            bRebuildSharedEdges=True,
            bRebuildVertices=True,
            bReplace=True,
            bDot=False,
            iDotHeight=None,
            bEcho=True,
            bDebug=bDebug
    )

    if gBreps_EdgesRebuilt:
        print "  {} of {} breps have had their edges rebuilt.".format(len(gBreps_EdgesRebuilt), len(gBreps))
    else:
        print "  No edges were rebuilt, so breps will not be further processed."
        bRepeat = False
    
    if not bRepeat:
        return
    
    iStep = 3
    
    s = "STEP {}: 2nd merging of edges".format(iStep)
    Rhino.RhinoApp.SetCommandPrompt(prompt=s+' ...')
    
    gBreps_WithEdgesMerged = mergeEdges(gBreps, fAngleTol)
    print s
    if not gBreps_WithEdgesMerged:
        print "None of the {} edges / {} naked edges were merged.".format(
                nEdges1, nNakedEdges1)
        nEdges2 = nEdges1
        nNakedEdges2 = nNakedEdges1
        print "Rebuilding edges not repeated because no edges were merged after the previous rebuild."
    else:
        iStep = 4
        
        Rhino.RhinoApp.SetCommandPrompt(prompt="Recording ending edge counts ...")
        nEdges2 = edgeCount_All(gBreps)
        nNakedEdges2 = edgeCount_Naked(gBreps)
        
        print createEdgeCountReport(nEdges1, nNakedEdges1, nEdges2, nNakedEdges2, len(gBreps))
        
        # 2nd Rebuild edges.
        s = "STEP {}: 2nd rebuilding of edges".format(iStep)
        Rhino.RhinoApp.SetCommandPrompt(prompt=s+' ...')
        
        sc.doc.Objects.UnselectAll()
        gBreps_EdgesRebuilt = spb_Brep_rebuildEdges.processBrepObjects(
                gBreps,
                fRebuildNakedTol=fRebuildEdgesTol,
                fSearchTol=2.0*fRebuildEdgesTol,
                bRebuildOnlyOutOfTol=bRebuildOnlyOutOfTol,
                bRebuildSharedEdges=True,
                bRebuildVertices=True,
                bReplace=True,
                bDot=False,
                iDotHeight=None,
                bEcho=True,
                bDebug=bDebug
        )
        print s
        if gBreps_EdgesRebuilt:
            print "{} of {} breps have had their edges rebuilt.".format(len(gBreps_EdgesRebuilt), len(gBreps))
        else:
            print "No edges were rebuilt."
    
    xBrep_compact.processBrepObjects(gBreps, bEcho=False, bDebug=bDebug)
    
    print createEdgeCountReport(nEdges0, nNakedEdges0, nEdges2, nNakedEdges2, len(gBreps), sPrefix="Total ")


def main():
    """
    """
    
    gBreps_Preselected = rs.SelectedObjects()
    
    rc = getInput()
    if rc is None: return
    (
        gBreps0,
        fAngleTol,
        bRebuildEdges,
        bRebuildOnlyOutOfTol,
        fRebuildEdgesTol,
        bRepeat,
        bEcho,
        bDebug,
    ) = rc
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")


    if not bDebug:
        sc.doc.Views.RedrawEnabled = False


    processBrepObjects(
            gBreps0,
            fAngleTol=fAngleTol,
            bRebuildEdges=bRebuildEdges,
            bRebuildOnlyOutOfTol=bRebuildOnlyOutOfTol,
            fRebuildEdgesTol=fRebuildEdgesTol,
            bRepeat=bRepeat,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    
    if gBreps_Preselected: rs.SelectObjects(gBreps_Preselected)
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()