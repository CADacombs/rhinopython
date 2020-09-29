"""
190530: Created.
191112: Added a function.
191113: Made input of some function more BrepFace focused.
191115: Added more tolerance multipliers.
200130, 0202: Modified output of a function.
200209, 14: Modified debug output.
200303: Now accepts curves with length < resolution and >= ModelAbsoluteTolerance.
200422: Refactored Opts.  Added some options.  Modified an option default value.
200430: Added bOnlyUseCrvsOnFace and option to use all normal wires and brep edges.
200505: Bug fix.
200519-23, 0619,24,25, 0701,29: Refactored.  Exported a function to its own script.  Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

import xBrep_splitSurfaceWithCurves
import xBrepObject


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


    #key = 'bSelBrep'; keys.append(key)
    #values[key] = True
    #names[key] = 'BrepSelMode'
    #riOpts[key] = ri.Custom.OptionToggle(values[key], 'Edge', 'Brep')
    #riAddOpts[key] = addOptionToggle(key, names, riOpts)
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bScanForNakedEdges'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Off', 'On')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPickPtForKeepNotRemove'; keys.append(key)
    values[key] = False
    names[key] = 'AtPointPicked'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Remove', 'Keep')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOnlyUseCrvsOnFace'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitToCrvSegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bTryOtherTolsOnFail'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
    values[key] = False
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
        for key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_TrimmingObjects():
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


    #idxs_EdgesOfFaceToSplit = []
    #rdBrep_withFaceToSplit, rgBrep_withFaceToSplit = getBrep(objref_Face)
    #gBrep_withFaceToSplit = rdBrep_withFaceToSplit.Id
    #idxs_EdgesOfFaceToSplit.extend(objref_Face.Face().AdjacentEdges())

    go = ri.Custom.GetObject()

    #def notEdgeOfFaceToSplit(rdObj, geom, compIdx):
    #    print rdObj, geom, compIdx
    #    if isinstance(rdObj, rd.BrepObject) and rdObj.Id == gBrep_withFaceToSplit:
    #        if compIdx.ComponentIndexType == rg.ComponentIndexType.BrepEdge:
    #            if geom.EdgeIndex in idxs_EdgesOfFaceToSplit:
    #                print "An edge of a face to split was picked and will not be used."
    #                return False
    #    return True
    #go.SetCustomGeometryFilter(notEdgeOfFaceToSplit)
    
    
    go.AcceptNothing(True)
    
    go.AcceptNumber(enable=True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False)
    go.EnableUnselectObjectsOnExit(False)
    
    idxs_Opts = {}

    while True:
        if Opts.values['bScanForNakedEdges']:
            go.SetCommandPromptDefault("Enter for all normal wires and brep naked edges")
        else:
            go.SetCommandPromptDefault("Enter for all normal wires")

        key = 'bScanForNakedEdges'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]

        go.GeometryFilter = rd.ObjectType.Curve | rd.ObjectType.Brep
        go.SubObjectSelect = False
        go.SetCommandPrompt("Select curves or breps")
        #if Opts.values['bSelBrep']:
        #    go.GeometryFilter = rd.ObjectType.Curve | rd.ObjectType.Brep
        #    go.SubObjectSelect = False
        #    go.SetCommandPrompt("Select curves or breps")
        #else:
        #    go.GeometryFilter = rd.ObjectType.Curve
        #    go.SubObjectSelect = True
        #    go.SetCommandPrompt("Select curves or edges")

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            return
        elif res == ri.GetResult.Nothing:
            return (
                [],
                Opts.values['bScanForNakedEdges'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return (
                objrefs,
                Opts.values['bScanForNakedEdges'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fTolerance'].CurrentValue = go.Number()
        
        if Opts.riOpts['fTolerance'].CurrentValue < 0.0:
            Opts.riOpts['fTolerance'].CurrentValue = Opts.riOpts['fTolerance'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getAllNormalObjs(bScanForNakedEdges=True):
    rdObjs = []
    settings = rd.ObjectEnumeratorSettings()
    settings.NormalObjects = True
    settings.LockedObjects = False
    for rdObj in sc.doc.Objects.GetObjectList(settings):
        if rdObj.ObjectType == rd.ObjectType.Curve:
            rdObjs.append(rdObj)
        elif rdObj.ObjectType == rd.ObjectType.Brep:
            if bScanForNakedEdges:
                rdObjs.append(rdObj)
    return rdObjs


def getInput_Face():
    """Get Brepface with optional input."""
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face to trim")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(enable=True, acceptZero=True)
    
    idxs_Opts = {}

    while True:
        key = 'bPickPtForKeepNotRemove'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bOnlyUseCrvsOnFace'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bSplitToCrvSegs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bTryOtherTolsOnFail'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        #key = 'bExtract'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            sc.doc.Objects.UnselectAll()
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


def trimSurfaceIntoFace(rgSrf_In, rhObjects_CurveOrEdge, pt_Input=None, **kwargs):
    """
    WIP: This is a Rhino.Geometry function verses the already working Rhino.DocObjects function below.
    TODO: Complete this function.

    Parameters:
        rgSrf_In: Can be rg.BrepFace or other rg.Surface.
        rhObjects_CurveOrEdge
        pt_Input: None or (rg.Point3d, bPickPtForKeepNotRemove)
        bUseCrvsBBoxIfNoPt
        bOnlyUseCrvsOnFace
        bSplitToCrvSegs
        fTolerance
        bTryOtherTolsOnFail
        bDebug
    Returns on success:
        rg.Brep (Monoface)
    Returns on fail: None
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bOnlyUseCrvsOnFace = getOpt('bOnlyUseCrvsOnFace')
    bSplitToCrvSegs = getOpt('bSplitToCrvSegs')
    fTolerance = getOpt('fTolerance')
    bTryOtherTolsOnFail = getOpt('bTryOtherTolsOnFail')
    bDebug = getOpt('bDebug')


    rgB_fromSplit = xBrep_splitSurfaceWithCurves.splitSurfaceIntoBrep(
            rgSrf_toSplit=srf_toSplit,
            rhObjects_CurveOrEdge=rhObjects_CurveOrEdge,
            fTolerance=fTolerance,
            bTryOtherTolsOnFail=bTryOtherTolsOnFail,
            bDebug=bDebug)
    if not rgB_fromSplit: return


    if pt_Input and isinstance(pt_Input[0], rg.Point3d):
        if pt_Input[1] == True:
            for rgF in rgB_fromSplit.Faces:
                b, u, v = rgF.ClosestPoint(pt_Input)
                if b:
                    if rgF.IsPointOnFace(u, v) == rg.PointFaceRelation.Interior:
                        rgB_Out = rgB_fromSplit.Faces[rgF.FaceIndex].DuplicateFace(False)
                        rgB_fromSplit.Dispose()
                        return rgB_Out
            else:
                rgB_fromSplit.Dispose()
                return
        elif pt_Input[1] == False:
            if rgB_fromSplit.Faces.Count > 2:
                print "More than 2 faces for point trim!"
                rgB_fromSplit.Dispose()
                return

            for rgF in rgB_fromSplit.Faces:
                b, u, v = rgF.ClosestPoint(pt_Input)
                if b:
                    if rgF.IsPointOnFace(u, v) == rg.PointFaceRelation.Interior:
                        rgB_Out = rgB_fromSplit.Faces[(rgF.FaceIndex+1)%2].DuplicateFace(False)
                        rgB_fromSplit.Dispose()
                        return rgB_Out
            else:
                rgB_fromSplit.Dispose()
                return


    # Use bounding box.

    def createBrepFromMatchingFace(breps, crvs, tolerance=None, bDebug=False):
        """
        Returns:
            Index of brep in breps.

        This function is used by other scripts, so leave it at the root level.
        """

        if tolerance is None:
            tolerance = 2.0 * sc.doc.ModelAbsoluteTolerance



        def createBB_ofGeometry(geoms):
            try: geoms = list(geoms)
            except: geoms = [geoms]
            bbox = rg.BoundingBox.Empty
            for geom in geoms: bbox.Union(geom.GetBoundingBox(True))
            if bbox is None: raise ValueError("Bounding box not created.")
            return bbox


        def findMatchingBoundingBox(bboxes, bbox_toMatch):
            """
            bboxes: Iterable of BoundingBoxes from which to search.
            bbox_toMatch: BoundingBox to match.
            Returns: Index of single matching bounding box on success.
            """
        
            if bDebug: print 'findMatchingBoundingBox()...'
        
            #tolerance = (3.0 * (tol_Split**2.0))**0.5 # == Total distance of tol_Split for each component.
        
            idx_Found = None

            epsilon = tolerance

            for i in range(len(bboxes)):
                if (
                        bbox_toMatch.Min.EpsilonEquals(bboxes[i].Min, epsilon=epsilon) and
                        bbox_toMatch.Max.EpsilonEquals(bboxes[i].Max, epsilon=epsilon)
                ):
                    if idx_Found is not None:
                        print "More than one match found.  Check results."
                        return
                    idx_Found = i
        
            if idx_Found is None:
                if bDebug:
                    print "Matching bounding box NOT found."
                    #sc.doc.Objects.AddBrep(bbox_toMatch.ToBrep())
                    #for bb in bboxes:
                        #sc.doc.Objects.AddBrep(bb.ToBrep())

            return idx_Found


        rgBs_fromSplit = [rgF.DuplicateFace(False) for rgF in rgB_fromSplit.Faces]


        crvs_Closed = []
        segs_ofOpen = []
        for crv in crvs:
            if crv.IsClosed:
                crvs_Closed.append(crv)
            else:
                segs_ofOpen.append(crv)

        segs_Joined = rg.Curve.JoinCurves(
                segs_ofOpen, joinTolerance=sc.doc.ModelAbsoluteTolerance)
        segs_Joined = list(segs_Joined) # Convert from Array.

        bbox_of_crvs = createBB_ofGeometry(crvs_Closed + segs_Joined)

        #sc.doc.Objects.AddBox(rg.Box(bbox_of_crvs)); #sc.doc.Views.Redraw(); 1/0

        bboxes_ofNEs = []
        for brep in breps:
            NEs = [edge for edge in brep.Edges
                   if edge.Valence == rg.EdgeAdjacency.Naked]
            bboxes_ofNEs.append(createBB_ofGeometry(NEs))

        #for bbox in bboxes_ofNEs:
            #sc.doc.Objects.AddBox(rg.Box(bbox))
        #sc.doc.Views.Redraw(); 1/0

        #bboxes_ofNEs = [
        #        createBB_ofGeometry(
        #            [edge
        #                for edge in brep.Edges
        #                if edge.Valence == rg.EdgeAdjacency.Naked]
        #        )
        #        for brep in breps]

        #map(sc.doc.Objects.AddBox, map(rg.Box, bboxes_ofNEs))

        idx = findMatchingBoundingBox(bboxes_ofNEs, bbox_of_crvs)

        if idx is None:
            for rgB in rgBs_fromSplit: rgB.Dispose()
            return

        rgB_Found = rgBs_fromSplit[idx]

        for i, rgB in enumerate(rgBs_fromSplit):
            if i != idx:
                rgB.Dispose()

        return rgB_Found



    rgB_Found = createBrepFromMatchingFace(
        rgB_fromSplit,
        rhObjects_CurveOrEdge,
        tolerance=sc.doc.ModelAbsoluteTolerance,
        bDebug=bDebug)

    rgB_fromSplit.Dispose()

    return rgB_Found


def getRhinoObjsOfSplitters(rhObjs, rhObj_withBrepFace):
    rdObj_withFace = rs.coercerhinoobject(rhObj_withBrepFace)

    rdObjs = []
    for rhObj in rhObjs:
        rdObj = rs.coercerhinoobject(rhObj)
        if isinstance(rdObj, (rd.BrepObject, rd.CurveObject)):
            if rdObj.Id == rdObj_withFace.Id:
                #print "Object to trim ({}) removed from trimming objects.".format(rdObj)
                continue
            rdObjs.append(rdObj)
    return rdObjs


def processBrepObject(rhObj_BrepFace, rhObjs_Splitters, **kwargs):
    """
    """


    def setOpt(key, value=None):
        if key in kwargs:
            return kwargs[key]
        elif key in Opts.riOpts:
            return Opts.riOpts[key].InitialValue
        else:
            return value

    bPickPtForKeepNotRemove = setOpt('bPickPtForKeepNotRemove')
    bOnlyUseCrvsOnFace = setOpt('bOnlyUseCrvsOnFace')
    bSplitToCrvSegs = setOpt('bSplitToCrvSegs')
    fTolerance = setOpt('fTolerance')
    bTryOtherTolsOnFail = setOpt('bTryOtherTolsOnFail')
    bExtract = setOpt('bExtract')
    bEcho = setOpt('bEcho')
    bDebug = setOpt('bDebug')


    rdObjs_Splitters = getRhinoObjsOfSplitters(rhObjs_Splitters, rhObj_BrepFace)

    if not rdObjs_Splitters:
        print "No splitter objects."
        return


    pt_Picked = rhObj_BrepFace.SelectionPoint()


    gBs_Split = xBrep_splitSurfaceWithCurves.processBrepObjects(
        [rhObj_BrepFace],
        rdObjs_Splitters,
        bSplitUnderlyingSrf=False,
        bOnlyUseCrvsOnSrf=bOnlyUseCrvsOnFace,
        bSplitToCrvSegs=bSplitToCrvSegs,
        fTolerance=fTolerance,
        bTryOtherTolsOnFail=bTryOtherTolsOnFail,
        bExplode=False,
        bExtract=True,
        bEcho=bEcho,
        bDebug=bDebug)

    if not gBs_Split:
        print "Face was not split."
        return

    if len(gBs_Split) > 1:
        print "More than 1 brep resulted in split."
        return

    gB_Split = gBs_Split[0]

    rgB_Split = rs.coercebrep(gB_Split)


    def faceAtPoint(rgBrep, pt_onFace, bDebug=False):
        """
        Returns face index.
        """


        # Do not use the following because sometimes the face edges are
        # stretched beyond the underlying surface (object of ClosestPoint).
        #for idxF, rgF in enumerate(rgBrep.Faces):
        #    b, u, v = rgF.ClosestPoint(pt_onFace)
        #    if not b: raise ValueError("ClosestPoint returned None.")
        #    ptFaceRel = rgF.IsPointOnFace(u, v)
        #    if ptFaceRel != rg.PointFaceRelation.Exterior:
        #        return idxF


        # The following is more exact, although also probably more expensive.
        # Point may be slightly outside of face.
        # Find point closest to each face.
        fDists = []
        for idxF, rgF in enumerate(rgBrep.Faces):
            rgB = rgF.DuplicateFace(duplicateMeshes=False)
            pt_Closest = rgB.ClosestPoint(pt_onFace)
            rgB.Dispose()
            fDist = pt_Closest.DistanceTo(pt_onFace)
            fDists.append(fDist)

        fMinDist = min(fDists)

        if fDists.count(fMinDist) > 1:
            raise ValueError("Face with point not on its exterior not found.")

        idxF = fDists.index(fMinDist)

        return idxF


    idxF_atPick = faceAtPoint(rgB_Split, pt_Picked)
    if idxF_atPick is None:
        if bOutputSplitOnFail: sc.doc.Objects.AddBrep(rgB_Split)
        rgB_Split.Dispose()
        return


    def sortFacesByIslands(rgBrep, idxF_Start):
        """
        Returns:
            list(int(Index of start face and islands))
            list(int(Index of remaining faces))
        """
        rgB = rgBrep
        if 1 <= rgB.Faces.Count <= 2:
            return [idxF_Start], [0 if idxF_Start else 1]
        
        idxFs_Keep = []
        idxFs_Skip = []

        idxFs_toAddToSkip = []


        idxFs_Adj_toNewSkip = [idxF_Start]

        while True:
            sc.escape_test()

            # Add to keep.
            idxFs_toAddToKeep = []
            for idxF in idxFs_Adj_toNewSkip:
                if idxF not in idxFs_Keep and idxF not in idxFs_Skip:
                    idxFs_toAddToKeep.append(idxF)

            if not idxFs_toAddToKeep:
                # All faces have been tested.
                return idxFs_Keep, idxFs_Skip

            idxFs_Keep.extend(idxFs_toAddToKeep)

            idxFs_Adj_toNewKeep = []
            for idxF in idxFs_toAddToKeep:
                idxFs_Adj_toNewKeep.extend(rgB.Faces[idxF].AdjacentFaces())

            idxFs_Adj_toNewKeep = list(set(idxFs_Adj_toNewKeep))

            idxFs_toAddToSkip = []
            for idxF in idxFs_Adj_toNewKeep:
                if idxF not in idxFs_Skip and idxF not in idxFs_Keep:
                    idxFs_toAddToSkip.append(idxF)

            if not idxFs_toAddToSkip:
                # All faces have been tested.
                return idxFs_Keep, idxFs_Skip

            idxFs_Skip.extend(idxFs_toAddToSkip)

            idxFs_Adj_toNewSkip = []
            for idxF in idxFs_toAddToSkip:
                idxFs_Adj_toNewSkip.extend(rgB.Faces[idxF].AdjacentFaces())

            idxFs_Adj_toNewSkip = list(set(idxFs_Adj_toNewSkip))


    idxFs_withStart, idxFs_notWithStart = sortFacesByIslands(rgB_Split, idxF_atPick)

    if bPickPtForKeepNotRemove:
        gBs_fromTrim = xBrepObject.removeFaces(gB_Split, idxFs_notWithStart)
    else:
        gBs_fromTrim = xBrepObject.removeFaces(gB_Split, idxFs_withStart)

    if not gBs_fromTrim:
        if bEcho:
            print "Trim failed."

    gBs_1F = []

    for gB in gBs_fromTrim:
        gBs_Extracted, gBs_Remaining = xBrepObject.extractFaces(gB, None, bEcho=False)
        gBs_1F.extend(gBs_Extracted)
        if len(gBs_Extracted) > 1:
            if bEcho:
                print "Adjacent faces in result.  This may be due to some ambiguity in the input."

    if bEcho:
        print "Face was trimmed into {} faces.".format(len(gBs_1F))

    return gBs_fromTrim


def main():

    rc = getInput_TrimmingObjects()
    if rc is None: return
    (
        rhObjs_Splitters,
        bScanForNakedEdges,
        bEcho,
        bDebug,
        ) = rc
    #print len(rc)

    if not rhObjs_Splitters:
        # Get all Normal objects now so before any layer, etc., changes occur during next input.
        rhObjs_Splitters = getAllNormalObjs(bScanForNakedEdges)

    if not rhObjs_Splitters:
        print "No splitters."
        return

    while True:
        sc.escape_test()
        rc = getInput_Face()
        if rc is None: return

        (
            objref_Face,
            bScanForNakedEdges,
            bPickPtForKeepNotRemove,
            bOnlyUseCrvsOnFace,
            bSplitToCrvSegs,
            fTolerance,
            bTryOtherTolsOnFail,
            bExtract,
            bEcho,
            bDebug,
            ) = rc
        #print len(rc)

        sc.doc.Objects.UnselectAll()

        Rhino.RhinoApp.SetCommandPrompt("Working ...")

        if not bDebug: sc.doc.Views.RedrawEnabled = False

        gBs_fromTrim = processBrepObject(
            objref_Face,
            rhObjs_Splitters,
            bPickPtForKeepNotRemove=bPickPtForKeepNotRemove,
            bOnlyUseCrvsOnFace=bOnlyUseCrvsOnFace,
            bSplitToCrvSegs=bSplitToCrvSegs,
            fTolerance=fTolerance,
            bTryOtherTolsOnFail=bTryOtherTolsOnFail,
            bExtract=bExtract,
            bEcho=bEcho,
            )

        sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()