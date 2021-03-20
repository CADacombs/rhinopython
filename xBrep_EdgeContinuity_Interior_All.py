"""
210320: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import scriptcontext as sc


def main():
    
    gBs = rs.GetObjects("Select polysurfaces", filter=16,
            preselect=True, select=True)
    if not gBs: return
    
    sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll() # For first call to _EdgeContinuity.
    
    Rhino.RhinoApp.SetCommandPrompt("")
    
    for iB in xrange(len(gBs)):
        gB = gBs[iB]
        
        rdB = rs.coercerhinoobject(gB)
        rgB = rdB.Geometry
        
        compIdxs = [None, None]
        
        for idxE in xrange(rgB.Edges.Count):
            edge = rgB.Edges[idxE]
            
            sPrompt = "Processing edge {} of {}".format(
                    idxE+1, rgB.Edges.Count)
            if len(gBs) > 1:
                sPrompt += " in polysurface {} of {}".format(
                    iB+1, len(gBs))
            
            sPrompt += " ..."
            
            Rhino.RhinoApp.CommandPrompt = sPrompt
            
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
                
                rc = Rhino.RhinoApp.RunScript("_EdgeContinuity _Enter", echo=False)
                
                Rhino.RhinoApp.CommandPrompt = sPrompt
                
                for compIdx in compIdxs:
                    rdB.SelectSubObject(
                        componentIndex=compIdx,
                        select=False,
                        syncHighlight=True,
                        persistentSelect=True)

    sc.doc.Views.RedrawEnabled = False


if __name__ == '__main__': main()