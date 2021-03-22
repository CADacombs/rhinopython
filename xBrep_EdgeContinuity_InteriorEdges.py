"""
210320-21: Created.
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


    key = 'bAllEdgesInSelBrep'; keys.append(key)
    values[key] = False
    names[key] = 'SelectionMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Edges', 'Polysrfs')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bMatchNakedPairs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fPairMatchingDistTol'; keys.append(key)
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

        if key == 'fPairMatchingDistTol':
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

        if Opts.values['bAllEdgesInSelBrep']:
            go.GeometryFilter = rd.ObjectType.Brep # Curve is used here for brep edges.
            go.SetCommandPrompt("Select polysurfaces")
        else:
            go.DisablePreSelect()
            go.GeometryFilter = rd.ObjectType.Curve # Curve is used here for brep edges.
            go.GeometryAttributeFilter = (
                    ri.Custom.GeometryAttributeFilter.MatedEdge)
            #go.GeometryAttributeFilter = (
            #        ri.Custom.GeometryAttributeFilter.MatedEdge |
            #        ri.Custom.GeometryAttributeFilter.EdgeCurve)
            go.SetCommandPrompt("Select interior edges")


        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bAllEdgesInSelBrep')
        #addOption('bMatchNakedPairs')
        #addOption('fPairMatchingDistTol')
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
                Opts.values['bAllEdgesInSelBrep'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        #if res == ri.GetResult.Number:
        #    key = 'fScale'
        #    Opts.riOpts[key].CurrentValue = go.Number()
        #    Opts.setValue(key)
        #    continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def processBreps(rhBreps, idxEdges=None, bEcho=True, bDebug=False):
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


def main():

    if Rhino.RhinoApp.ExeVersion < 7:
        print "This script calls _EdgeContinuity," \
            " a command that was introduced in Rhino 7."
        return

    rc = getInput()
    if rc is None: return
    (
        objrefs,
        bAllEdgesInSelBrep,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False


    if bAllEdgesInSelBrep:
        processBreps(
            rhBreps=[o.ObjectId for o in objrefs],
            idxEdges=None,
            bEcho=bEcho,
            bDebug=bDebug)
    else:
        gBs = []
        iEs = []
        for objref in objrefs:
            gB = objref.ObjectId
            if gB in gBs:
                iB = gBs.index(gB)
                iEs[iB].append(objref.Edge().EdgeIndex)
            else:
                gBs.append(gB)
                iEs.append([objref.Edge().EdgeIndex])
        processBreps(
            rhBreps=gBs,
            idxEdges=iEs,
            bEcho=bEcho,
            bDebug=bDebug)


    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()