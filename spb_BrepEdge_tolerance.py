"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
221104: Created.
250105: Modified number formatting of printed output.
250223: Added the edge index to the printed output.
"""

import Rhino
import Rhino.Input as ri
import scriptcontext as sc


def formatDistance(fDistance, iPrecision=None):
    if iPrecision is None: iPrecision = sc.doc.ModelDistanceDisplayPrecision
    if fDistance is None: return "(No value provided)"
    if fDistance == 0.0: return "0"
    if fDistance < 10**-iPrecision: return "{:.2e}".format(fDistance)
    return "{:.{}g}".format(fDistance, iPrecision)


def main():
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edge")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.EdgeCurve

    res = go.Get()

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    objref = go.Object(0)

    go.Dispose()

    edge = objref.Edge()

    print("Edge[{}] tolerance: {}".format(edge.EdgeIndex, formatDistance(edge.Tolerance)))


if __name__ == '__main__': main()