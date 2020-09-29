"""
200520-23: Created, starting with another script.
200526: Bug fix.  Modified a function description.
200619: Import-related update.
200629: Now checks and rejects curves completely on surface/face border.
200810: Replaced some local code with a function from an import.

TODO: Add capability of tracking CurveObjects in main routine so that CurveObjects can be selecting instead of duplicated.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from clr import StrongBox
from System import Array
from System import Guid

import xBrepFace
import xCurve


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


    key = 'bProcessSegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUnderlyingSrf'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bCompletelyOnSrfAndClosed'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCompletelyOnSrfAndOpen'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPartiallyOnSrf'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSamplingResolution'; keys.append(key)
    values[key] = 100.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAdd'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

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


def getInput_Faces():
    """Get Brepface with optional input."""
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select faces")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(enable=True, acceptZero=True)
    
    idxs_Opts = {}

    while True:
        key = 'bProcessSegs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bUnderlyingSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bCompletelyOnSrfAndClosed'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bCompletelyOnSrfAndOpen'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bPartiallyOnSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fSamplingResolution'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bAdd'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
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


def getInput_Curves(objrefs_Face):
    """
    Get curves or edges with optional input.
    """


    def getBrep(rhBrep):
        if isinstance(rhBrep, rg.Brep):
            return None, rhBrep
        elif isinstance(rhBrep, rg.GeometryBase):
            rdObj = None
            rgObj = rhBrep
        elif isinstance(rhBrep, rd.ObjRef):
            rdObj = rhBrep.Object()
            rgObj = rhBrep.Geometry()
        elif isinstance(rhBrep, Guid):
            rdObj = sc.doc.Objects.FindId(rhBrep) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhBrep)
            rgObj = rdObj.Geometry
        else:
            return

        if isinstance(rgObj, (rg.Brep, rg.BrepFace)):
            return rdObj, rgObj


    idxs_EdgesOfFaceToSplit = []
    for objref_Face in objrefs_Face:
        rdBrep_withFaceToSplit, rgBrep_withFaceToSplit = getBrep(objref_Face)
        gBrep_withFaceToSplit = rdBrep_withFaceToSplit.Id
        idxs_EdgesOfFaceToSplit.extend(objref_Face.Face().AdjacentEdges())


    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves or edges")
    go.SetCommandPromptDefault("Enter for all normal wires and brep naked edges")

    go.GeometryFilter = rd.ObjectType.Curve
    
    def notEdgeOfFaceToSplit(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx
        if isinstance(rdObj, rd.BrepObject) and rdObj.Id == gBrep_withFaceToSplit:
            if geom.EdgeIndex in idxs_EdgesOfFaceToSplit:
                print "An edge of a face to split was picked and will not be used."
                return False
        return True
    go.SetCustomGeometryFilter(notEdgeOfFaceToSplit)
    
    
    go.AcceptNothing(True)
    
    go.AcceptNumber(enable=True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False)
    go.EnableUnselectObjectsOnExit(False)
    
    idxs_Opts = {}

    while True:
        key = 'bProcessSegs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bUnderlyingSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bCompletelyOnSrfAndClosed'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bCompletelyOnSrfAndOpen'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bPartiallyOnSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fSamplingResolution'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bAdd'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
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
                elif rdObj.ObjectType == rd.ObjectType.Brep:
                    rgBrep = rdObj.BrepGeometry
                    # Structure of following conditional is to optimize speed.
                    if rdObj.Id != gBrep_withFaceToSplit:
                        for edge in rgBrep.Edges:
                            if edge.Valence == rg.EdgeAdjacency.Naked:
                                rgEdges.append(edge)
                    else:
                        for edge in rgBrep.Edges:
                            if edge.EdgeIndex not in idxs_EdgesOfFaceToSplit:
                                if edge.Valence == rg.EdgeAdjacency.Naked:
                                    rgEdges.append(edge)
                            else:
                                if Opts.values['bDebug']:
                                    print "Skipped edge {}.".format(edge.EdgeIndex)
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


def duplicateCurvesOnSurface(rgCrvs_In, rgSrf, fSamplingResolution=None, fTolerance=None, bDebug=False):
    """
    Parameteters:
        curves
        rgSrf: Honors edges when this is a rg.BrepFace.
        fSamplingResolution: float of curve division length.
        fTolerance: float of distance from curve to surface.
        bDebug
    Returns:
        (
            list(Closed curves of rgCrvs_In completely on rgSrf),
            list(Open curves of rgCrvs_In completely on rgSrf),
            list(Curves of rgCrvs_In only partially on rgSrf)
        )
    """


    if fSamplingResolution is None:
        fSamplingResolution = 100.0*sc.doc.ModelAbsoluteTolerance
    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance


    def isCurveCompletelyOnFaceBorder(crv):


        strongBox_points = StrongBox[Array[rg.Point3d]]()

        rc = crv.DivideByLength(
            fSamplingResolution,
            includeEnds=True,
            points=strongBox_points)

        if rc:
            pts = list(strongBox_points.Value)
            #for pt in pts: sc.doc.Objects.AddPoint(pt)
            #sc.doc.Views.Redraw(); 1/0
            if len(pts) == 2:
                rc = crv.DivideByCount(
                    segmentCount=2,
                    includeEnds=True,
                    points=strongBox_points)
                pts = list(strongBox_points.Value)
        else:
            crv_GetLength = crv.GetLength()
            
            if crv_GetLength <= sc.doc.ModelAbsoluteTolerance:
                if bDebug:
                    print "No points for curve that is {} long.".format(
                        crv.GetLength())
                    #sc.doc.Objects.AddCurve(crv)
                return
            else:
                rc = crv.DivideByCount(
                    segmentCount=2,
                    includeEnds=True,
                    points=strongBox_points)
                pts = list(strongBox_points.Value)
                #for pt in pts: sc.doc.Objects.AddPoint(pt)
                #sc.doc.Views.Redraw(); 1/0

        for pt in pts:
            rc = xBrepFace.is3dPointOnFace(
                rgFace, pt, fTolerance)

            if rc != rg.PointFaceRelation.Boundary:
                return False

        return True


    def isCrvCompletelyOnFace(crv):

        # Quick check for False.
        rc = xBrepFace.is3dPointOnFace(
            rgFace, crv.PointAtStart, fTolerance)

        if not rc:
            # This includes None, False, and rg.PointFaceRelation.Exterior.
            return False


        # More thorough check.

        strongBox_points = StrongBox[Array[rg.Point3d]]()

        rc = crv.DivideByLength(
            fSamplingResolution,
            includeEnds=True,
            points=strongBox_points)

        if rc:
            pts = list(strongBox_points.Value)
            #for pt in pts: sc.doc.Objects.AddPoint(pt)
            #sc.doc.Views.Redraw(); 1/0
        else:
            crv_GetLength = crv.GetLength()
            
            if crv_GetLength <= sc.doc.ModelAbsoluteTolerance:
                if bDebug:
                    print "No points for curve that is {} long.".format(
                        crv.GetLength())
                    #sc.doc.Objects.AddCurve(crv)
                return
            else:
                rc = crv.DivideByCount(
                    segmentCount=2,
                    includeEnds=True,
                    points=strongBox_points)
                pts = list(strongBox_points.Value)
                #for pt in pts: sc.doc.Objects.AddPoint(pt)
                #sc.doc.Views.Redraw(); 1/0

        for pt in pts:
            rc = xBrepFace.is3dPointOnFace(
                rgFace, pt, fTolerance)

            if not rc:
                # This includes None, False, and rg.PointFaceRelation.Exterior.
                return False

        return True


    if isinstance(rgSrf, rg.BrepFace):
        rgFace = rgSrf
        rgB_Temp = None
    else:
        rgB_Temp = rgSrf.ToBrep()
        rgFace = rgB_Temp.Faces[0]


    # Full curves completely, not completely on Surface.
    fulls_CompletelyOn, fulls_PartiallyOn = [], []
    for c in rgCrvs_In:
        if isCrvCompletelyOnFace(c):
            if isCurveCompletelyOnFaceBorder(c):
                continue
            fulls_CompletelyOn.append(c.DuplicateCurve())
        else:
            fulls_PartiallyOn.append(c.DuplicateCurve())
    #map(sc.doc.Objects.AddCurve, NBs_AllOnSrf); return


    if rgB_Temp: rgB_Temp.Dispose()


    if not fulls_CompletelyOn:
        return [], [], fulls_PartiallyOn

    fulls_CompletelyOn_Closed = []
    fulls_CompletelyOn_Open = []
    for c in fulls_CompletelyOn:
        if c.IsClosed:
            fulls_CompletelyOn_Closed.append(c)
        else:
            fulls_CompletelyOn_Open.append(c)
    return fulls_CompletelyOn_Closed, fulls_CompletelyOn_Open, fulls_PartiallyOn


def processRhinoObjects(rhObjects_BrepFace, rhObjects_CurveOrEdge, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bProcessSegs = getOpt('bProcessSegs')
    bUnderlyingSrf = getOpt('bUnderlyingSrf')
    fTolerance = getOpt('fTolerance')
    bCompletelyOnSrfAndClosed = getOpt('bCompletelyOnSrfAndClosed')
    bCompletelyOnSrfAndOpen = getOpt('bCompletelyOnSrfAndOpen')
    bPartiallyOnSrf = getOpt('bPartiallyOnSrf')
    fSamplingResolution = getOpt('fSamplingResolution')
    bAdd = getOpt('bAdd')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getRhinoObject(rhObj):
        """
        'Deleted objects cannot be found by id.'
        (https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_DocObjects_Tables_ObjectTable_FindId.htm)
        """
        rdObj = None
        if isinstance(rhObj, rd.RhinoObject):
            rdObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        return rdObj


    def getBrepObject(rhObj):
        rdObj = getRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def getSortedBrepIdsAndFaces(objrefs):
        """
        Parameters:
            list(objrefs)
        Returns:
            list(Brep GUIDs)
            list(lists(integers of Face indices) per brep)
        """
        
        gBreps_In = []
        idxs_Faces_perBrep = []
    
        for o in objrefs:
            gBrep0 = o.ObjectId
            rdBrep_In = o.Object()
            rgBrep_In = o.Brep()
        
            if not rgBrep_In.IsValid:
                print "Brep {} is invalid.  Fix first.".format(gBrep0)
                rgBrep_In.Dispose()
                continue
        
            idx_CompIdx = o.GeometryComponentIndex.Index
            if idx_CompIdx == -1:
                if gBrep0 in gBreps_In:
                    idxs_Faces_perBrep[gBreps_In.index(gBrep0)] = range(rgBrep_In.Faces.Count)
                else:
                    gBreps_In.append(gBrep0)
                    idxs_Faces_perBrep.append(range(rgBrep_In.Faces.Count))
            else:
                rgFace_Brep0 = o.Face()
                if gBrep0 in gBreps_In:
                    if rgFace_Brep0 in idxs_Faces_perBrep[gBreps_In.index(gBrep0)]:
                        continue
                    else:
                        idxs_Faces_perBrep[gBreps_In.index(gBrep0)].append(rgFace_Brep0.FaceIndex)
                else:
                    gBreps_In.append(gBrep0)
                    idxs_Faces_perBrep.append([rgFace_Brep0.FaceIndex])

        return gBreps_In, idxs_Faces_perBrep


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


    gBreps_In, idxs_rgFace_perBrep = getSortedBrepIdsAndFaces(rhObjects_BrepFace)
    if not gBreps_In: return


    rdCrvs_In = []
    rgCrvs_In = []
    for o in rhObjects_CurveOrEdge:
        rdCrv_In, rgCrv_In = getCurve(o)
        if rdCrv_In:
            rdCrvs_In.append(rdCrv_In)
            rgCrvs_In.append(rgCrv_In)


    if bProcessSegs:
        rgCrvs_toSort = xCurve.duplicateSegments(
            rgCrvs_In,
            bExplodePolyCrvs=True)
    else:
        rgCrvs_toSort = rgCrvs_In



    gCrvs_Ret = []

    for iB, gBrep0 in enumerate(gBreps_In):
        rdBrep_In = getBrepObject(gBrep0)
        rgBrep_In = rdBrep_In.Geometry

        idxFaces_TrimSuccess = []
        rgBs_1F_Res = []

        for iF, idxFace in enumerate(idxs_rgFace_perBrep[iB]):
            if bUnderlyingSrf:
                rgSrf_Ref = rgBrep_In.Faces[idxFace].UnderlyingSurface()
            else:
                rgSrf_Ref = rgBrep_In.Faces[idxFace]

            rc = duplicateCurvesOnSurface(
                    rgCrvs_toSort,
                    rgSrf_Ref,
                    fSamplingResolution=fSamplingResolution,
                    fTolerance=fTolerance,
                    bDebug=bDebug)

            if not rc: continue
            
            (
                cs_AllOnSrf_Closed,
                cs_AllOnSrf_Open,
                cs_PartiallyOnSrf,
                ) = rc

            if bCompletelyOnSrfAndClosed:
                for c in cs_AllOnSrf_Closed:
                    gCrv_Ret = sc.doc.Objects.AddCurve(c)
                    if gCrv_Ret != Guid.Empty:
                        gCrvs_Ret.append(gCrv_Ret)

            if bCompletelyOnSrfAndOpen:
                for c in cs_AllOnSrf_Open:
                    gCrv_Ret = sc.doc.Objects.AddCurve(c)
                    if gCrv_Ret != Guid.Empty:
                        gCrvs_Ret.append(gCrv_Ret)

            if bPartiallyOnSrf:
                for c in cs_PartiallyOnSrf:
                    gCrv_Ret = sc.doc.Objects.AddCurve(c)
                    if gCrv_Ret != Guid.Empty:
                        gCrvs_Ret.append(gCrv_Ret)

    return gCrvs_Ret


def main():


    rc = getInput_Faces()
    if rc is None: return

    objrefs_Face = rc[0]

    rc = getInput_Curves(objrefs_Face)
    if rc is None: return

    rhObjects_CurveOrEdge = rc[0]

    #sc.doc.Objects.UnselectAll()

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gCrvs_Ret = processRhinoObjects(
        rhObjects_BrepFace=objrefs_Face,
        rhObjects_CurveOrEdge=rhObjects_CurveOrEdge,
        )

    if gCrvs_Ret:
        sc.doc.Objects.UnselectAll()
        for gC in gCrvs_Ret:
            sc.doc.Objects.Select(objectId=gC)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()