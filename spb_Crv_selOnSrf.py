"""
200520-23: Created, starting with another script.
200526: Bug fix.  Modified a function description.
200619: Import-related update.
200629: Now checks and rejects curves completely on surface/face border.
200810: Replaced some local code with a function from an import.
210113: Trialing using more tolerance for checking curves against full surface borders.  See TODO.
        Debugging bug fix.
220317: Moved main function to a main library module.

TODO:
    WIP: Allow edge selection.
    Add capability of tracking CurveObjects in main routine so that CurveObjects can be selected instead of duplicated.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

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
        key = 'bUnderlyingSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
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
        key = 'bUnderlyingSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
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


def processRhinoObjects(rhObjects_BrepFace, rhObjects_CurveOrEdge, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bUnderlyingSrf = getOpt('bUnderlyingSrf')
    fTolerance = getOpt('fTolerance')
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
        
        if isinstnace(rhObj, rd.BrepObject):
            pass
        
        if isinstance(rhObj, rd.CurveObject):
            return rhObj, rhObj.CurveGeometry

        if isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
            rgObj = rhObj.Geometry()
        elif isinstance(rhObj, Guid):
            try: rdObj = sc.doc.Objects.FindId(rhObj)
            except: rdObj = sc.doc.Objects.Find(rhObj)
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


    gCs_Found = []

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

            for rdC, rgC, in zip(rdCrvs_In, rgCrvs_In):

                if rdC.Id in gCs_Found:
                    continue

                rc = xCurve.filterCurvesOnSurface(
                    rgC,
                    rgSrf_Ref,
                    fSamplingResolution=fSamplingResolution,
                    fTolerance=fTolerance,
                    bDebug=bDebug)
    
                if not rc: continue
                
                (
                    cs_AllOnSrf,
                    cs_PartiallyOnSrf,
                    ) = rc
    
                if cs_AllOnSrf:
                    gCs_Found.append(rdC.Id)
                    continue
    
                if bPartiallyOnSrf and cs_PartiallyOnSrf:
                    gCs_Found.append(rdC.Id)
                    continue
    
    return gCs_Found


def main():

    rc = getInput_Faces()
    if rc is None: return

    objrefs_Face = rc[0]

    rc = getInput_Curves(objrefs_Face)
    if rc is None: return

    rhObjects_CurveOrEdge = rc[0]

    #sc.doc.Objects.UnselectAll()

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gCs_Ret = processRhinoObjects(
        rhObjects_BrepFace=objrefs_Face,
        rhObjects_CurveOrEdge=rhObjects_CurveOrEdge,
        )

    if gCs_Ret:
        sc.doc.Objects.UnselectAll()
        for gC in gCs_Ret:
            sc.doc.Objects.Select(objectId=gC)
        sc.doc.Views.Redraw()
        print("{} curves found.".format(len(gCs_Ret)))
    else:
        print("No curves found.")


if __name__ == '__main__': main()