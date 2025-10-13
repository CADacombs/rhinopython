"""
Multiple edges may still be found but their lengths are within machine epsilon (2**-52)
of each other.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
210425: Created.
210916: Removed an unneeded print(statement.
221214: Bug fix for object filter.
251001: Now, epsilon used is machine epsilon. Multiple edges may still result.
251012: Removed an unused global variable.

TODO: Don't allow block instances to be selected.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color


def getBreps():

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
            "Select breps <All normal when none are selected>",
            acceptNothing=True,
            filter=rd.ObjectType.Brep)

    if res == Rhino.Commands.Result.Cancel:
        return
    
    #    for o in objrefs:
    #        gci = o.GeometryComponentIndex
    #        print(gci.ComponentIndexType
    #        print(o.Brep()
    #    return
    
    if objrefs is None:
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        rdBrepObjects = []
        for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
            if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
                rdBrepObjects.append(rdRhinoObject)
    else:
        rdBrepObjects = []
        for objref in objrefs:
            rdBrepObjects.append(objref.Object())
    return rdBrepObjects


def getMinimumLengthEdgeIndices(rgBrep0):
    """
    """
    
    fLengths = []
    
    for rgE in rgBrep0.Edges:
        fLengths.append(rgE.GetLength())

    fLength_Min = min(fLengths)
    idx_Mins = []
    
    for i, fLength in enumerate(fLengths):
        if fLength-fLength_Min <= sc.doc.ModelAbsoluteTolerance:
            idx_Mins.append(i)

    return fLength_Min, idx_Mins


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
        guid = None
    elif isinstance(rhObj, rd.ObjRef):
        geom = rhObj.Geometry()
        guid = rhObj.ObjectId
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
        guid = rdObj.Id
    elif isinstance(rhObj, rd.BrepObject):
        rdObj = rhObj
        geom = rdObj.Geometry
        guid = rdObj.Id
    else:
        return

    if isinstance(geom, rg.Brep):
        if not geom.IsValid:
            print("Brep{} is invalid and will be skipped.  " \
                    "_ExtractBadSrf and repair the bad faces first.".format(
                    '' if guid is None else " {}".format(guid)))
            return
        
        return geom


def coerceRhinoObject(rhObj):
    rdObj = None
    if isinstance(rhObj, rd.BrepObject):
        rdObj = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        rdObj = rhObj.Object()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
    return rdObj


def processBrepObjects(rdBrepObjects, bIsRational=True):
    """Selects edges and adds dots.
    """
    
    try:
        rdBrepObjects = list(rdBrepObjects)
    except:
        rdBrepObjects = [rdBrepObjects]
    
    fLengths_Min = []
    idx_Es_Mins_PerB = []
    
    rdBs_Checked = []
    
    for iB, rdBrep0 in enumerate(rdBrepObjects):
        
        Rhino.RhinoApp.CommandPrompt = "Searching ..."
        
        rgBrep0 = coerceBrep(rdBrep0)
        if rgBrep0 is None:
            continue
        
        rdBs_Checked.append(rdBrep0)
        
        rc = getMinimumLengthEdgeIndices(rgBrep0=rgBrep0)
        
        fLengths_Min.append(rc[0])
        idx_Es_Mins_PerB.append(rc[1])

    if len(rdBs_Checked) == 0:
        print("No valid breps.")
        return

    fLength_Min_All = min(fLengths_Min)

    rdBs_WithMins = []
    idx_Es_Mins_OfAll = []
    rgDots = []
    fEpsilon = 2**(-52) # Machine epsilon.

    for i, fLength_Min in enumerate(fLengths_Min):
        if fLength_Min-fLength_Min_All > fEpsilon:
            continue

        # Brep has shortest edge(s).

        rdBs_WithMins.append(rdBs_Checked[i])
        rdBrep0 = rdBs_Checked[i]
        rgBrep0 = rdBrep0.Geometry

        for idx_E_Min in idx_Es_Mins_PerB[i]:
            compIdx = rg.ComponentIndex(
                type=rg.ComponentIndexType.BrepEdge,
                index=idx_E_Min)
            rdBrep0.SelectSubObject(
                componentIndex=compIdx,
                select=True,
                syncHighlight=True,
                persistentSelect=True)


            sDot = '{0:.{1}f}'.format(fLength_Min, sc.doc.ModelDistanceDisplayPrecision+1)
            rgEdge = rgBrep0.Edges[idx_E_Min]
            pts = rgEdge.DivideByCount(2, False)
            pt = rgEdge.PointAtStart if pts is None else rgEdge.PointAt(pts[0])
            rgDot = rg.TextDot(sDot, pt)
            rgDot.FontHeight = 11
            rgDots.append(rgDot)


        idx_Es_Mins_OfAll.extend(idx_Es_Mins_PerB[i])


    if len(rgDots) > 0:
        attrib = rd.ObjectAttributes()
        attrib.ColorSource = rd.ObjectColorSource.ColorFromObject
        attrib.ObjectColor = Color.Red
        for rgDot in rgDots:
            sc.doc.Objects.Select(sc.doc.Objects.AddTextDot(rgDot, attrib))


    sMinPhrase = "within ~{0:.{2}e} of ~{1:.{2}g}".format(
        fEpsilon,
        fLength_Min_All,
        sc.doc.ModelDistanceDisplayPrecision+1)


    if len(rdBrepObjects) == 1:
            print("{} edge(s) found in brep {}.".format(
                len(idx_Es_Mins_OfAll),
                sMinPhrase))
    else:
        if len(idx_Es_Mins_OfAll) == 1:
            print("1 edge found in {} out of {} brep(s) {}.".format(
                len(rdBs_WithMins),
                len(rdBrepObjects),
                sMinPhrase))
        else:
            print("{} edges found in {} out of {} brep(s) {}.".format(
                len(idx_Es_Mins_OfAll),
                len(rdBs_WithMins),
                len(rdBrepObjects),
                sMinPhrase))


def main():
    
    rc = getBreps()
    if rc is None: return
    rdBrepObjs = rc
    
    sc.doc.Objects.UnselectAll()
    
    processBrepObjects(
        rdBrepObjs,
        bIsRational=True,
        )

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()