"""
190413: Created by extracting functions from xBrep.  Updated reference.
200110: Import-related update.
"""

import Rhino.Geometry as rg
import scriptcontext as sc

import xBrepTrim


def hasMultipleTrimsOnAnyNaturalEdges(rgLoop):
    W = sum(1 for rgTrim in rgLoop.Trims if rgTrim.IsoStatus == rg.IsoStatus.West)
    S = sum(1 for rgTrim in rgLoop.Trims if rgTrim.IsoStatus == rg.IsoStatus.South)
    E = sum(1 for rgTrim in rgLoop.Trims if rgTrim.IsoStatus == rg.IsoStatus.East)
    N = sum(1 for rgTrim in rgLoop.Trims if rgTrim.IsoStatus == rg.IsoStatus.North)
    return (W > 1 or S > 1 or E > 1 or N > 1)


def endVertexOfLastTrim(rgLoop):
    if rgLoop.Trims.Count == 0: return
    
    for t in range(rgLoop.Trims.Count - 1, -1, -1):
        if rgLoop.Trims[t].TrimType == rg.BrepTrimType.Singular:
            continue
        
        return xBrepTrim.endVertexOfEdgeTrim(rgLoop.Trims[t])
