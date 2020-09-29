"""
200521: Created, starting with another script.
200611, 24, 25: Bug fix.
200729: Renamed a function.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid


class Opts():
    
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


    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
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


def getInput_Curves():
    """
    Get wires with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select wire curves")
    go.SetCommandPromptDefault("Enter for all normal wires")

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve
    
    go.AcceptNothing(True)
    
    go.AcceptNumber(enable=True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False)
    go.EnableUnselectObjectsOnExit(False)
    
    idxs_Opts = {}

    while True:
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            return
        elif res == ri.GetResult.Nothing:
            rdCrvs = []
            rgEdges = []
            settings = rd.ObjectEnumeratorSettings()
            settings.NormalObjects = True
            settings.LockedObjects = False
            for rdObj in sc.doc.Objects.GetObjectList(settings):
                if rdObj.ObjectType == rd.ObjectType.Curve:
                    rdCrvs.append(rdObj)
            return tuple([rdCrvs + rgEdges] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fTolerance'].CurrentValue = go.Number()
        
        if Opts.riOpts['fTolerance'].CurrentValue < 0.0:
            Opts.riOpts['fTolerance'].CurrentValue = Opts.riOpts['fTolerance'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInput_Face():
    """Get Brepface with optional input."""
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(enable=True, acceptZero=True)
    
    idxs_Opts = {}

    while True:
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            return
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fTolerance'].CurrentValue = go.Number()
        
        if Opts.riOpts['fTolerance'].CurrentValue < 0.0:
            Opts.riOpts['fTolerance'].CurrentValue = Opts.riOpts['fTolerance'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def formatDistance(fDistance):
    if fDistance is None:
        return "(No deviation was provided.)"
    if fDistance == 0.0:
        return "Exactly zero"
    if fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-3)):
        return "{:.1e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def splitCurvesWithSurfaceEdges(rgCrvs_In, rgSrf_In, fTolerance=None, bDebug=False):
    """
    Parameters:
        rgCrvs_In
        rgSrf_In: If rg.BrepFace, then adjacent edges are used, otherwise surface domain isocurves are used.
        fTolerance
        bDebug
    Returns:
        rgCrvs_Out_SrfInterior
        rgCrvs_Out_SrfBoundary
        rgCrvs_Out_SrfExterior
        
    Intersect.Intersection.CurveCurve sometimes returns intersections at the exact ends of curves.
    In this case, the entire curve will be the result of Curve.Split.
    """
    
    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance

    rgEdges = [] # Face edges or Surface domain extent curves.

    if isinstance(rgSrf_In, rg.BrepFace):
        idxs_Edges = rgSrf_In.AdjacentEdges()
        for idxE in idxs_Edges:
            rgE = rgSrf_In.Brep.Edges[idxE]
            rgEdges.append(rgE)
    elif isinstance(rgSrf_In, rg.Surface):
        # SENW == 0,1,2,3 for IsSingular.
        
        # South.
        if not rgSrf_In.IsSingular(0):
            rgC = rgSrf_In.IsoCurve(0, rgSrf_In.Domain(1).T0)
            rgEdges.append(rgC)
        
        # East.
        if not rgSrf_In.IsSingular(1):
            rgC = rgSrf_In.IsoCurve(1, rgSrf_In.Domain(0).T1)
            rgEdges.append(rgC)
        
        # North.
        if not rgSrf_In.IsSingular(2):
            rgC = rgSrf_In.IsoCurve(0, rgSrf_In.Domain(1).T1)
            rgEdges.append(rgC)
        
        # West.
        if not rgSrf_In.IsSingular(3):
            rgC = rgSrf_In.IsoCurve(1, rgSrf_In.Domain(1).T0)
            rgEdges.append(rgC)
    else:
        raise ValueError("{} passed to splitCurvesWithSurfaceEdges".format(
            rgSrf_In.GetType().Name))


    rgCrvs_Out_SrfInterior = []
    rgCrvs_Out_SrfBoundary = []
    rgCrvs_Out_SrfExterior = []

    for rgC_toSplit in rgCrvs_In:
        t_Pts = []
        doms_Overlaps = []
        
        for rgEdge in rgEdges:
            intrscts = rg.Intersect.Intersection.CurveCurve(
                curveA=rgC_toSplit,
                curveB=rgEdge,
                tolerance=fTolerance,
                overlapTolerance=0.1*fTolerance)

            if bDebug: print "Intersection count: {}".format(intrscts.Count)

            for intsct in intrscts:
                if intsct.IsPoint:
                    #print intsct.ParameterA
                    t_Pts.append(intsct.ParameterA)
                if intsct.IsOverlap:
                    doms_Overlaps.append(intsct.OverlapA)
                    #print intsct.OverlapA.T0, intsct.OverlapA.T1
                    t_Pts.append(intsct.OverlapA.T0)
                    t_Pts.append(intsct.OverlapA.T1)

        t_Pts = sorted(set(t_Pts))
        
        # TODO: Remove near-duplicates.
        
        rgCs_fromSplit = rgC_toSplit.Split(t=t_Pts)


        # For debugging.
        if len(rgCs_fromSplit) == 1:
            pass
            #print "1 curve resulted from split."
            #c = rgCs_fromSplit[0]
            #sc.doc.Objects.AddCurve(c)
            #continue


        for c in rgCs_fromSplit:

            fLength = c.GetLength()

            if fLength == 0.0:
                if bDebug: print "Zero length curve ignored."
                c.Dispose()
                continue
            elif fLength < fTolerance:
                if bDebug: print "Short ({}) curve ignored.  Domain: [{},{}]".format(
                    formatDistance(fLength),
                    c.Domain.T0,
                    c.Domain.T1)
                c.Dispose()
                continue

            for dom in doms_Overlaps:
                if c.Domain.IncludesParameter(dom.Mid):
                    rgCrvs_Out_SrfBoundary.append(c)
                    break # To next curve from split.
            else:
                # Curve is not part of Intersection overlaps but may still be along boundary.

                pt_MidCrv = c.PointAt(c.Domain.Mid)
                #sc.doc.Objects.AddPoint(pt_MidCrv)

                b, u, v = rgSrf_In.ClosestPoint(pt_MidCrv) # if rgSrf_In is BrepFace, UnderlyingSurface is still used.
                if not b:
                    raise ValueError("ClosestPoint failed.")

                ptAtUV = rgSrf_In.PointAt(u,v)
                #sc.doc.Objects.AddPoint(ptAtUV)


                # Check whether Closest point on surfaces is beyond tolerance from the test point.
                dist = ptAtUV.DistanceTo(pt_MidCrv)
                if dist > fTolerance:
                    rgCrvs_Out_SrfExterior.append(c)
                    continue


                # Point is either interior or along an edge.

                if isinstance(rgSrf_In, rg.BrepFace):

                    # TODO: Check distance from edge?

                    ptFaceRel = rgSrf_In.IsPointOnFace(u,v)
                    if ptFaceRel == rg.PointFaceRelation.Interior:
                        rgCrvs_Out_SrfInterior.append(c)
                    else:
                        rgCrvs_Out_SrfBoundary.append(c)
                else:
                    # rgSrf_In is a Surface but not a BrepFace.

                    # TODO: Check distance from surface side?

                    if (
                        rgSrf_In.Domain(0).IncludesParameter(u) and
                        rgSrf_In.Domain(1).IncludesParameter(v)
                    ):
                        rgCrvs_Out_SrfInterior.append(c)
                    else:
                        rgCrvs_Out_SrfBoundary.append(c)


    return rgCrvs_Out_SrfInterior, rgCrvs_Out_SrfBoundary, rgCrvs_Out_SrfExterior


def processCurveObjects(rhCrvObjs, rhSrf, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getCurve(rhObj):
        if isinstance(rhObj, rd.CurveObject):
            return rhObj, rhObj.CurveGeometry

        if isinstance(rhObj, rg.Curve):
            return None, rhObj

        if isinstance(rhObj, rg.GeometryBase):
            rdObj = None
            rgObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
            rgObj = rhObj.Geometry()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
            rgObj = rdObj.Geometry
        else:
            return

        if isinstance(rgObj, rg.Curve):
            return rdObj, rgObj


    def getSurface(rhSrf):
        """
        Parameters:
            rhFace
        Returns:
            rg.BrepFace or rg.Surface
        """
        
        if isinstance(rhSrf, rg.BrepFace):
            return rhSrf

        if isinstance(rhSrf, rg.Surface):
            return rhSrf

        if isinstance(rhSrf, rg.Brep):
            rgBrep = rhSrf
            if rgBrep.Faces.Count == 1:
                return rgBrep.Faces[0]
            else:
                return

        if isinstance(rhSrf, rd.ObjRef):
            objref = rhSrf
            compIdxType = objref.GeometryComponentIndex.ComponentIndexType

            if compIdxType == rg.ComponentIndexType.BrepFace:
                return objref.Face()

            if compIdxType == rg.ComponentIndexType.InvalidType:
                rgBrep = objref.Brep()
                if rgBrep.Faces.Count == 1:
                    return rgBrep.Faces[0]
                else:
                    return

        if isinstance(rhSrf, Guid):
            rdObj = sc.doc.Objects.FindId(rhSrf) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhSrf)
            if isinstance(rdObj, rd.BrepObject):
                brep = rdObj.Geometry
                if rhSrf.Faces.Count == 1:
                    return rhSrf.Faces[0]
                else:
                    return

        srf = None
        if isinstance(geom, rg.BrepFace):
            srf = geom.UnderlyingSurface()
        elif isinstance(geom, rg.Surface):
            srf = geom

        return srf


    rdCrvs_In = []
    rgCrvs_In = []
    for o in rhCrvObjs:
        rdCrv_In, rgCrv_In = getCurve(o)
        if rdCrv_In:
            rdCrvs_In.append(rdCrv_In)
            rgCrvs_In.append(rgCrv_In)

    rgSrf_In = getSurface(rhSrf)

    gCrvs_Out = []
    
    iCt_Split_AllCrv_In = 0
    iCt_CrvsSplit_NoError = 0
    
    for i, rgCrv in enumerate(rgCrvs_In):
        
        attrs = rdCrvs_In[i].Attributes

        rc = splitCurvesWithSurfaceEdges(
            [rgCrv],
            rgSrf_In=rgSrf_In,
            fTolerance=fTolerance)

        if not rc: continue
        
        (
            cs_Interior,
            cs_Boundary,
            cs_Exterior
            ) = rc

        iCt_Added_ThisCrv_In = 0
        
        crvs_toAdd = cs_Interior + cs_Boundary + cs_Exterior
        
        if not crvs_toAdd: continue
        
        for c in crvs_toAdd:
            gC = sc.doc.Objects.AddCurve(c, attributes=attrs)
            if gC != Guid.Empty:
                gCrvs_Out.append(gC)
                iCt_Added_ThisCrv_In += 1
        
        if iCt_Added_ThisCrv_In == len(crvs_toAdd):
            sc.doc.Objects.Delete(objectId=rdCrvs_In[i].Id, quiet=False)
            iCt_Split_AllCrv_In += iCt_Added_ThisCrv_In
            iCt_CrvsSplit_NoError += 1
        else:
            print "Added {} out of {} curves.  Original curve was not removed.".format(
                iCt_Added_ThisCrv_In, len(crvs_toAdd))

    if iCt_CrvsSplit_NoError:
        print "Split {} curves into {} curves.".format(
            iCt_CrvsSplit_NoError, iCt_Split_AllCrv_In)
    else:
        print "No curves were split."

    return gCrvs_Out


def main():


    rc = getInput_Curves()
    if rc is None: return
    rhObjects_Wires = rc[0]

    rc = getInput_Face()
    if rc is None: return
    objrefs_Face = rc[0]

    #sc.doc.Objects.UnselectAll()

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gCrvs_Ret = processCurveObjects(
        rhCrvObjs=rhObjects_Wires,
        rhSrf=objrefs_Face,
        )

    if gCrvs_Ret:
        sc.doc.Objects.UnselectAll()
        for gC in gCrvs_Ret:
            sc.doc.Objects.Select(objectId=gC)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()