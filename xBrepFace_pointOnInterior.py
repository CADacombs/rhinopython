"""
If the BrepFace's centroid or mid of its Surface Domains are not on the Face's interior,
a random point is generated.

200108-10: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System import Random


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'fMinDistFromBorder'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(bDebug=False):
    """
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface
    go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.SubSurface) # Doesn't allow single surfaces to be selected.
    
    go.AcceptNumber(enable=True, acceptZero=True)
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == ri.GetResult.Object 
    
    while True:
        go.AddOptionDouble(Opts.names['fMinDistFromBorder'], Opts.riOpts['fMinDistFromBorder'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.Get()
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return [objref] + [Opts.values[key] for key in Opts.keys]
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fMinDistFromBorder'].CurrentValue = go.Number()
        
        if Opts.riOpts['fMinDistFromBorder'].CurrentValue < 0.0:
            Opts.riOpts['fMinDistFromBorder'].CurrentValue = Opts.riOpts['fMinDistFromBorder'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def createPoint3d(rgFace, fMinDistFromBorder=None, bDebug=False):
    """
    Parameters:
        rgFace: BrepFace
        fMinDistFromBorder: float
    Returns:
        Point3d on success, None on failure
    """
    
    if fMinDistFromBorder is None:
        fMinDistFromBorder = 10.0 * sc.doc.ModelAbsoluteTolerance
    
    rgBrep = rgFace.Brep
    
    iDomainU = rgFace.Domain(0)
    iDomainV = rgFace.Domain(1)

    areaMassProp = Rhino.Geometry.AreaMassProperties.Compute(rgFace)
    if areaMassProp is None:
        print "Face[{}]'s AreaMassProperties cannot be calculated.".format(
            rgFace.FaceIndex)
        u = iDomainU.Mid
        v = iDomainV.Mid
    else:
        ptCentrdW = areaMassProp.Centroid
        bSuccess, u, v = rgFace.ClosestPoint(ptCentrdW)
    
    rgEdges = [rgBrep.Edges[idxEdge] for idxEdge in rgFace.AdjacentEdges()]

    rand = Random()
    
    for i in xrange(10000):
        pointFaceRelation = rgFace.IsPointOnFace(u, v)

        if bDebug: sEval = "pointFaceRelation"; print sEval+':',eval(sEval)

        if pointFaceRelation == rg.PointFaceRelation.Interior:
            pt = rgFace.PointAt(u, v)

            # If point is not at least fMinDistFromBorder from border (all edges),
            # continue searching.
            for rgEdge in rgEdges:
                b, t = rgEdge.ClosestPoint(pt)
                if b:
                    dist = pt.DistanceTo(rgEdge.PointAt(t))
                    if dist < fMinDistFromBorder:
                        break # to get another u and v.
            else: # Good point
                map(lambda x: x.Dispose(), rgEdges)
                rgFace.Dispose()
                return pt

        # Get new parameters for point.
        u = rand.NextDouble() * iDomainU.Length + iDomainU.Min
        v = rand.NextDouble() * iDomainV.Length + iDomainV.Min


def main(bEcho=True, bDebug=False):
    
    rc = getInput()
    if rc is None: return
    (
            objref0,
            fMinDistFromBorder,
            bEcho,
            bDebug,
    ) = rc
    if objref0 is None: return
    
    gBrep0 = objref0.ObjectId
    rdBrep0 = objref0.Object()
    rgBrep0 = objref0.Brep()

    if rgBrep0.Faces.Count == 1:
        rgFace0 = rgBrep0.Faces[0]
    else:
        idx_CompIdx = objref0.GeometryComponentIndex.Index
        if idx_CompIdx == -1:
            rgFace0 = rgBrep0.Faces[0]
        else:
            rgFace0 = objref0.Face()

    pt3d = createPoint3d(
        rgFace=rgFace0,
        fMinDistFromBorder=fMinDistFromBorder,
        bDebug=bDebug)

    if not pt3d:
        print "A point on the face could not be found."
        return

    gPt = sc.doc.Objects.AddPoint(pt3d)
    if gPt != Guid.Empty:
        print "Point on face was added to document."
        sc.doc.Views.Redraw()


if __name__ == '__main__': main(bEcho=bool(1), bDebug=bool(0))