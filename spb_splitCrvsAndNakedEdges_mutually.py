"""
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220715-17: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc


class _Data:
    def __init__(self):
        self.rdCs = []
        self.rdBs = []
        self.rgNEs_perB = []
        self.pts_Cs = []
        self.pts_NEs_perB = []


def _getInput():
    
    res, objrefs = Rhino.Input.RhinoGet.GetMultipleObjects(
            "Select curves and breps <All normal when none are selected>",
            acceptNothing=True,
            filter=rd.ObjectType.Brep | rd.ObjectType.Curve)
    if res == Rhino.Commands.Result.Cancel:
        return
    elif objrefs:
        return [o.Object() for o in objrefs]
    else:
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.ObjectTypeFilter = rd.ObjectType.Brep | rd.ObjectType.Curve
        return [o for o in sc.doc.Objects.GetObjectList(iter)]


def _sortInput(rdObjs, d):
    bFound = False
    
    for rdObj in rdObjs:
        if rdObj.ObjectType == rd.ObjectType.Curve:
            d.rdCs.append(rdObj)
            bFound = True
        else:
            d.rdBs.append(rdObj)
            d.rgNEs_perB.append(
                [rgE for rgE in rdObj.Geometry.Edges
                if rgE.Valence == rg.EdgeAdjacency.Naked])
            bFound = bFound or d.rgNEs_perB[-1]
    
    return bFound


def _collectIntersectionData(d):
    
    for i, rdC in enumerate(d.rdCs):
        d.pts_Cs.append([])
    
    for i, rgNEs in enumerate(d.rgNEs_perB):
        d.pts_NEs_perB.append([[] for _ in xrange(len(rgNEs))])

    bFound = False

    for iC, rdC in enumerate(d.rdCs):
        rgC = rdC.Geometry
        for iNEs, rgNEs in enumerate(d.rgNEs_perB):
            for iNE, rgNE in enumerate(rgNEs):
                ie = rg.Intersect.Intersection.CurveCurve(
                    curveA=rgC,
                    curveB=rgNE,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance,
                    overlapTolerance=0.0)
                if ie.Count == 0:
                    continue
                for iie in xrange(ie.Count):
                    if ie[iie].IsOverlap:
                        print("Overlap interesection was created and will be ignored.")
                        continue
                    d.pts_Cs[iC].append(ie[iie].PointA)
                    d.pts_NEs_perB[iNEs][iNE].append(ie[iie].PointB)
                bFound = True

    return bFound


def _cleanedPointList(pts, rgC, tol=None):
    if tol is None: tol = sc.doc.ModelAbsoluteTolerance
    pts_Out = list(rg.Point3d.CullDuplicates(pts, tolerance=tol))
    
    for i, pt in reversed(list(enumerate(pts_Out))):
        if pt.DistanceTo(rgC.PointAtStart) <= tol:
            del pts_Out[i]
        elif pt.DistanceTo(rgC.PointAtEnd) <= tol:
            del pts_Out[i]
    
    return pts_Out


def _splitCurveObjects(d):
    
    iCt_WereSplit = 0
    iCt_Result = 0
    
    for iC, rdC in enumerate(d.rdCs):
        rgC = rdC.Geometry
        pts = d.pts_Cs[iC]
        if len(pts) == 0: continue
        pts = _cleanedPointList(pts, rgC, tol=sc.doc.ModelAbsoluteTolerance)
        if len(pts) == 0: continue
        ts = [rgC.ClosestPoint(pt)[1] for pt in pts]
        rgCs_Out = rgC.Split(ts)
        if len(rgCs_Out) < 2:
            print("Curve {} was not split.".format(rdC.Id))
        for rgC_Out in rgCs_Out:
            sc.doc.Objects.AddCurve(rgC_Out)
            iCt_Result += 1
        if not sc.doc.Objects.Remove(rdC):
            print("Curve {} could not be deleted.".format(rdC.Id))
        iCt_WereSplit += 1
    
    if iCt_WereSplit == 0:
        print("No curves were split.")
    else:
        print("{} curves split into {}.".format(iCt_WereSplit, iCt_Result))


def _splitNakedEdges(d):
    iCt_WereSplit = 0
    iCt_Result = 0
    
    for iB, rdB in enumerate(d.rdBs):
        rgB = rdB.Geometry
        iCts_Splits_ThisB = 0
        for iE, rgE in enumerate(d.rgNEs_perB[iB]):
            pts = d.pts_NEs_perB[iB][iE]
            if len(pts) == 0: continue
            pts = _cleanedPointList(pts, rgE, tol=sc.doc.ModelAbsoluteTolerance)
            if len(pts) == 0: continue
            ts = [rgE.ClosestPoint(pt)[1] for pt in pts]
            idxE = rgE.EdgeIndex
            iCts_Splits = rgB.Edges.SplitEdgeAtParameters(idxE, ts)
            if iCts_Splits == 0:
                raise Exception("Edge was not split.")
            iCt_Result += iCts_Splits + 1
            iCts_Splits_ThisB += iCt_Result
            iCt_WereSplit += 1
        if iCts_Splits_ThisB == 0: continue
        if not rdB.CommitChanges():
            print("Brep {} could not be modified.".format(rdB.Id))
    
    if iCt_WereSplit == 0:
        print("No edges were split.")
    else:
        print("{} edges split into {}.".format(iCt_WereSplit, iCt_Result))


def main():
    
    rdObjs = _getInput()
    if rdObjs is None: return
    
    d = _Data()
    
    if not _sortInput(rdObjs, d):
        print("No curve or edge data in input.")
        return
    
    if not _collectIntersectionData(d):
        print("No intersections found.")
        return
    
    _splitCurveObjects(d)
    
    _splitNakedEdges(d)
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main()