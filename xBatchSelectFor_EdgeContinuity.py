"""
210320-21: Created.
210326: Added support for naked edge, including auto-matching.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAllMatedEdgesInSelBrep'; keys.append(key)
    values[key] = False
    names[key] = 'SelectionMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Edges', 'Polysrfs')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMatchNakedPairs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fMatchTol'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        if not idxOpt: print "Add option for {} failed.".format(key)

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fMatchTol':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-9:
                cls.riOpts[key].CurrentValue = 1e-9

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get breps or edges with optional input.
    """

    go = ri.Custom.GetObject()

    #go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:

        if Opts.values['bAllMatedEdgesInSelBrep']:
            go.GeometryFilter = rd.ObjectType.Brep # Curve is used here for brep edges.
            go.SetCommandPrompt("Select polysurfaces")
        else:
            go.DisablePreSelect()
            go.GeometryFilter = rd.ObjectType.Curve # Curve is used here for brep edges.
            if Opts.values['bMatchNakedPairs']:
                go.GeometryAttributeFilter = (
                        ri.Custom.GeometryAttributeFilter.EdgeCurve)
                go.SetCommandPrompt("Select edges")
            else:
                go.GeometryAttributeFilter = (
                        ri.Custom.GeometryAttributeFilter.MatedEdge)
                go.SetCommandPrompt("Select interior edges")


        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bAllMatedEdgesInSelBrep')
        if not Opts.values['bAllMatedEdgesInSelBrep']:
            addOption('bMatchNakedPairs')
            if Opts.values['bMatchNakedPairs']:
                addOption('fMatchTol')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()

            return (
                objrefs,
                Opts.values['bAllMatedEdgesInSelBrep'],
                Opts.values['bMatchNakedPairs'],
                Opts.values['fMatchTol'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fMatchTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def processBreps_MatedEdges(rhBreps, idxEdges=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBreps: list(GUIDs of Breps or rd.BrepObjects)
        idxEdges:
            None to process all edges or
            list(list of int(EdgeIndices) per brep in rhBreps)
    """

    sc.doc.Objects.UnselectAll() # For first call to _EdgeContinuity.

    Rhino.RhinoApp.SetCommandPrompt("")
    
    for iB in xrange(len(rhBreps)):
        rhB = rhBreps[iB]
        
        if isinstance(rhB, Guid):
            gB = rhB
            rdB = sc.doc.Objects.FindId(gB)
        elif isinstance(rhB, rd.BrepObject):
            rdB = rhB
        else:
            print "Item in rhBreps is {}, not BrepObject nor GUID of one.".format(
                rhB.GetType().Name)
            continue

        rgB = rdB.Geometry
        
        compIdxs = [None, None]

        if idxEdges is None:
            idxEs = xrange(rgB.Edges.Count)
        else:
            idxEs = idxEdges[iB]
        
        for idxE in idxEs:
            edge = rgB.Edges[idxE]
            
            sPrompt = "Processing edge {} of {}".format(
                    idxE+1, rgB.Edges.Count)
            if len(rhBreps) > 1:
                sPrompt += " in polysurface {} of {}".format(
                    iB+1, len(rhBreps))
            
            sPrompt += " ..."
            
            Rhino.RhinoApp.CommandPrompt = sPrompt
            if bDebug: print sPrompt
            
            if sc.escape_test(throw_exception=False):
                print "Script interrupted by user."
                return
            
            if edge.TrimCount == 2:
                iTs = edge.TrimIndices()
                for i in 0,1:
                    compIdxs[i] = rg.ComponentIndex(
                        type=rg.ComponentIndexType.BrepTrim,
                        index=iTs[i])
                for compIdx in compIdxs:
                    rdB.SelectSubObject(
                        componentIndex=compIdx,
                        select=True,
                        syncHighlight=True,
                        persistentSelect=True)
                
                rc = Rhino.RhinoApp.RunScript("_EdgeContinuity _Enter", echo=bDebug)
                
                Rhino.RhinoApp.CommandPrompt = sPrompt
                
                for compIdx in compIdxs:
                    rdB.SelectSubObject(
                        componentIndex=compIdx,
                        select=False,
                        syncHighlight=True,
                        persistentSelect=True)


def processBreps_NakedEdges(rhBreps, idxEdges=None, fMatchTol=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBreps: list(GUIDs of Breps or rd.BrepObjects)
        idxEdges:
            None to process all edges or
            list(int(EdgeIndices) per rhBreps order)
    """

    if fMatchTol is None: fMatchTol = 10.0 * sc.doc.ModelAbsoluteTolerance

    sc.doc.Objects.UnselectAll() # For first call to _EdgeContinuity.

    Rhino.RhinoApp.SetCommandPrompt("")


    def fill_rdBs_rgEs():
        for i, (rhB, idxE) in enumerate(zip(rhBreps, idxEdges)):
            if isinstance(rhB, Guid):
                rdB = sc.doc.Objects.FindId(rhB)
            elif isinstance(rhB, rd.BrepObject):
                rdB = rhB
            else:
                print "Item in rhBreps is {}, not BrepObject nor GUID of one.".format(
                    rhB.GetType().Name)
                continue
            rgB = rdB.BrepGeometry
            rdBs.append(rdB)
            rgEs.append(rgB.Edges[idxE])

    rdBs = []
    rgEs = []
    fill_rdBs_rgEs()

    if not rgEs or len(rgEs) < 2: return


    idxs_Matches_Pairs = [] # Pairs of edges by rgEs list indices (not EdgeIndex).
    idxs_Matches_Flat = [] # Indices matched.

    for iA, rgE_A in enumerate(rgEs):
        if iA in idxs_Matches_Flat: continue
        for iB in range(iA+1, len(rgEs)):
            if iB in idxs_Matches_Flat: continue
            rgE_B = rgEs[iB]
            ptA_Mid = None
            if (
                (
                    rgE_A.PointAtStart.DistanceTo(rgE_B.PointAtStart) <= fMatchTol
                    and
                    rgE_A.PointAtEnd.DistanceTo(rgE_B.PointAtEnd) <= fMatchTol
                )
                or
                (
                    rgE_A.PointAtStart.DistanceTo(rgE_B.PointAtEnd) <= fMatchTol
                    and
                    rgE_A.PointAtEnd.DistanceTo(rgE_B.PointAtStart) <= fMatchTol
                )
            ):
                if ptA_Mid is None:
                    ptA_Mid = rgE_A.PointAt(rgE_A.Domain.Mid)
                rc = rgE_B.ClosestPoint(ptA_Mid)
                if rc[0]:
                    if rgE_B.PointAt(rc[1]).DistanceTo(ptA_Mid) <= fMatchTol:
                        idxs_Matches_Pairs.append((iA, iB))
                        idxs_Matches_Flat.extend((iA, iB))

    if bDebug:
        print idxs_Matches_Pairs
        print idxs_Matches_Flat


    for i, (iA, iB) in enumerate(idxs_Matches_Pairs):

        sPrompt = "Processing edge pair {} of {}".format(
                i+1, len(idxs_Matches_Pairs))
        sPrompt += " ..."

        Rhino.RhinoApp.CommandPrompt = sPrompt
        if bDebug: print sPrompt
            
        if sc.escape_test(throw_exception=False):
            print "Script interrupted by user."
            return
        
        rdB_A = rdBs[iA]
        rdB_B = rdBs[iB]
        rgE_A = rgEs[iA]
        rgE_B = rgEs[iB]

        compIdx_A = rg.ComponentIndex(
            type=rg.ComponentIndexType.BrepTrim,
            index=rgE_A.TrimIndices()[0])
        compIdx_B = rg.ComponentIndex(
            type=rg.ComponentIndexType.BrepTrim,
            index=rgE_B.TrimIndices()[0])

        rdB_A.SelectSubObject(
            componentIndex=compIdx_A,
            select=True,
            syncHighlight=True,
            persistentSelect=True)

        rdB_B.SelectSubObject(
            componentIndex=compIdx_B,
            select=True,
            syncHighlight=True,
            persistentSelect=True)

        rc = Rhino.RhinoApp.RunScript("_EdgeContinuity _Enter", echo=bDebug)

        Rhino.RhinoApp.CommandPrompt = sPrompt

        rdB_A.SelectSubObject(
            componentIndex=compIdx_A,
            select=False,
            syncHighlight=True,
            persistentSelect=True)

        rdB_B.SelectSubObject(
            componentIndex=compIdx_B,
            select=False,
            syncHighlight=True,
            persistentSelect=True)


def main():

    if Rhino.RhinoApp.ExeVersion < 7:
        print "This script calls _EdgeContinuity," \
            " a command that was introduced in Rhino 7."
        return

    rc = getInput()
    if rc is None: return
    (
        objrefs,
        bAllMatedEdgesInSelBrep,
        bMatchNakedPairs,
        fMatchTol,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False


    if bAllMatedEdgesInSelBrep:
        processBreps_MatedEdges(
            rhBreps=[o.ObjectId for o in objrefs],
            idxEdges=None,
            bEcho=bEcho,
            bDebug=bDebug)
    else:
        gBs_ForMatedEs = []
        iEs_Mated = []

        gBs_NakedEs = []
        iEs_Naked = []

        for objref in objrefs:
            rgE = objref.Edge()
            if rgE.Valence == rg.EdgeAdjacency.Interior:
                gB = objref.ObjectId
                if gB in gBs_ForMatedEs:
                    iB = gBs_ForMatedEs.index(gB)
                    iEs_Mated[iB].append(rgE.EdgeIndex)
                else:
                    gBs_ForMatedEs.append(gB)
                    iEs_Mated.append([rgE.EdgeIndex])
            elif bMatchNakedPairs and rgE.Valence == rg.EdgeAdjacency.Naked:
                gBs_NakedEs.append(objref.ObjectId)
                iEs_Naked.append(rgE.EdgeIndex)
            else:
                print "Skipping edge of {} valence.".format(rgE.Valence)

        if gBs_ForMatedEs:
            processBreps_MatedEdges(
                rhBreps=gBs_ForMatedEs,
                idxEdges=iEs_Mated,
                bEcho=bEcho,
                bDebug=bDebug)

        if gBs_NakedEs:
            processBreps_NakedEdges(
                rhBreps=gBs_NakedEs,
                idxEdges=iEs_Naked,
                fMatchTol=fMatchTol,
                bEcho=bEcho,
                bDebug=bDebug)


    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()