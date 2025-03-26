"""
This script is an alternative to _SelBadObjects and _ExtractBadSrf.
In addtion to the breps in whoel, it optionally checks their components and geometry.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160819: Created.
181020-21:  Moved a function from another script.
            Added Opts class and more options.  Refactored.
            Changed name from findInvalidFaces.py to invalidBreps.py.
190505: History modification.
191022: Bug fix.  More general input is now allowed for extractBadFaces.
191024-25: Output of extractBadFaces changed from BrepObjects to GUIDS.
        Improved command line output for bEcho value in extractBadFaces.
191103: extractBadFaces can now be passed multiple breps.
250325: Added check for breps created from brep's surfaces. Refactored.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List
from System.Drawing import Color


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    keys_Type = []


    key = 'bFacesAsBreps'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    names[key] = 'FacesAsTheirOwnBreps'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSrfsAsBreps'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    names[key] = 'SrfsAsTheirOwnBreps'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCrvs2D'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCrvs3D'; keys.append(key)
    keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEdges'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bFaces'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bLoops'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSrfs'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bTrims'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bVertices'; keys.append(key)
    keys_Type.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAnalyzePerFace'; keys.append(key)
    values[key] = True
    names[key] = 'AnalysisMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Brep', 'Face')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddDot'; keys.append(key)
    values[key] = False
    names[key] = 'DotFaces'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iDotHeight'; keys.append(key)
    values[key] = 10
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]

    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        cls.stickyKeys[key]
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


    @classmethod
    def addOptions_BrepTypes(cls, getInput):
        
        def addOptionToggle(key):
            getInput.AddOptionToggle(cls.names[key], cls.riOpts[key])
        
        addOptionToggle('bFacesAsBreps')
        addOptionToggle('bCrvs2D')
        addOptionToggle('bCrvs3D')
        addOptionToggle('bEdges')
        addOptionToggle('bFaces')
        addOptionToggle('bLoops')
        addOptionToggle('bSrfs')
        addOptionToggle('bTrims')
        addOptionToggle('bVertices')
        getInput.AddOption('YesToAll')
        getInput.AddOption('NoToAll')
    
    @classmethod
    def chooseTypes(cls):
        
        getInputOpt = ri.Custom.GetOption()
        getInputOpt.SetCommandPrompt("Brep type(s) to check")
        
        cls.addOptions_BrepTypes(getInputOpt)
        
        while True:
            if sc.escape_test(False): break
            
            res = getInputOpt.Get()
            if res != ri.GetResult.Option:
                return
            
            if getInputOpt.OptionIndex() == 10:
                for key in cls.keys_Type:
                    cls.values[key] = cls.riOpts[key].CurrentValue = True
            elif getInputOpt.OptionIndex() == 11:
                for key in cls.keys_Type:
                    cls.values[key] = cls.riOpts[key].CurrentValue = False
            else:
                for key in cls.keys_Type:
                    cls.values[key] = cls.riOpts[key].CurrentValue
            
            cls.saveSticky()
            
            getInputOpt.ClearCommandOptions()
            cls.addOptions_BrepTypes(getInputOpt)
    
    @classmethod
    def processInput(cls, go):
        res = go.Result()
        
        if go.Option().Index == 1:
            cls.chooseTypes()
        
        cls.setValues()
        
        cls.saveSticky()
        
        # Clear and add options regardless if a number was entered or options were modified in another way.
        go.ClearCommandOptions()
        Opts.addOptions(go)


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")
    go.GeometryFilter = rd.ObjectType.Brep

    go.AcceptNothing(True)
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False
    go.SubObjectSelect = False

    go.EnableClearObjectsOnEntry(False)

    idxs_Opt = {}
    def addOption(ric, key): idxs_Opt[key] = Opts.addOption(ric, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        idxs_Opt['ObjectsToCheck'] = go.AddOption('ObjectsToCheck')
        addOption(go, 'bAnalyzePerFace')
        addOption(go, 'bExtract')
        addOption(go, 'bAddDot')
        if Opts.values['bAddDot']:
            addOption(go, 'iDotHeight')
        addOption(go, 'bEcho')
        addOption(go, 'bDebug')

        sTrue = []
        sFalse = []
        for key in Opts.keys_Type:
            if Opts.values[key]:
                sTrue.append(Opts.names[key])
            else:
                sFalse.append(Opts.names[key])
        print("Brep components that will be checked: {}".format(
            ", ".join(sTrue) if sTrue else 'None'))
        print("Not checked: {}".format(", ".join(sFalse) if sFalse else 'None'))

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Nothing:
            oes = rd.ObjectEnumeratorSettings()
            oes.LockedObjects = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep
            rdBs = list(sc.doc.Objects.GetObjectList(oes))
            go.Dispose()
            return rdBs

        if go.Option().Index == idxs_Opt['ObjectsToCheck']:
            idxs_Opt.clear()

            go_ObjFilter = ri.Custom.GetOption()
            go_ObjFilter.SetCommandPrompt("Subobjects to check")

            while True:
                go_ObjFilter.ClearCommandOptions()
                idxs_Opt.clear()

                for key in Opts.keys_Type:
                    addOption(go_ObjFilter, key)
                idxs_Opt['YesToAll'] = go_ObjFilter.AddOption('YesToAll')
                idxs_Opt['NoToAll'] = go_ObjFilter.AddOption('NoToAll')

                res = go_ObjFilter.Get()

                if res != ri.GetResult.Option:
                    break

                if go_ObjFilter.OptionIndex() == idxs_Opt['YesToAll']:
                    for key in Opts.keys_Type:
                        Opts.riOpts[key].CurrentValue = True
                        Opts.setValue(key)
                    continue

                if go_ObjFilter.OptionIndex() == idxs_Opt['NoToAll']:
                    for key in Opts.keys_Type:
                        Opts.riOpts[key].CurrentValue = False
                        Opts.setValue(key)
                    continue

                for key in idxs_Opt:
                    if go_ObjFilter.Option().Index == idxs_Opt[key]:
                        Opts.setValue(key, go_ObjFilter.Option().CurrentListOptionIndex)
                        break

            go_ObjFilter.Dispose()

            continue



        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break

        go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)


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

    if not isinstance(geom, rg.Brep):
        print("Not a brep: {}".format(guid))
        return
    
    return geom


def coerceRhinoObject(rhObj):
    rdObj = None
    if isinstance(rhObj, rd.RhinoObject):
        rdObj = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        rdObj = rhObj.Object()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
    return rdObj


def reportOfInvalidComponentCounts(rgBrep, bFacesAsBreps=True, bSrfsAsBreps=True, bCrvs2D=True, bCrvs3D=True, bEdges=True, bFaces=True, bLoops=True, bSrfs=True, bTrims=True, bVertices=True, bDebug=False):
    """
    """

    ss = []

    if bFacesAsBreps:
        idxs = []
        for iF in range(rgBrep.Faces.Count):
            rgF = rgBrep.Faces[iF]
            # Check duplicate brep from each face because face itself may show as being valid.
            rgB_fromF = rgF.DuplicateFace(duplicateMeshes=False)
            if not rgB_fromF.IsValid:
                idxs.append(iF)
            rgB_fromF.Dispose()
        if not idxs:
            ss.append("0 faces found that would be invalid breps.")
        elif len(idxs) > 10:
            ss.append("{} faces found that would be invalid brep: {}, etc.".format(len(idxs), idxs[:10]))
        else:
            ss.append("{} faces found that would be invalid breps: {}".format(len(idxs), idxs[:len(idxs)]))

    if bSrfsAsBreps:
        idxs = []
        for iF in range(rgBrep.Surfaces.Count):
            rgS = rgBrep.Surfaces[iF]
            # Check duplicate brep from each face because face itself may show as being valid.
            rgB_fromS = rgS.ToBrep()
            if not rgB_fromS.IsValid:
                idxs.append(iF)
        if not idxs:
            ss.append("0 surfaces found that would be invalid breps.")
        elif len(idxs) > 10:
            ss.append("{} surfaces found that would be invalid brep: {}, etc.".format(len(idxs), idxs[:10]))
        else:
            ss.append("{} surfaces found that would be invalid breps: {}".format(len(idxs), idxs[:len(idxs)]))

    return "\n".join(ss)



    if bSrfsAsBreps and not rgB_fromS.IsValid:
        idxsFs_InvalidFacesAsBreps.append(iF)
    if bCrvs2D and any(not rgX.IsValid for rgX in rgB_fromF.Curves2D):
        idxsFs_InvalidCrvs2d.append(iF)
    if bCrvs3D and any(not rgX.IsValid for rgX in rgB_fromF.Curves3D):
        idxsFs_InvalidCrvs3d.append(iF)
    if bEdges and any(not rgX.IsValid for rgX in rgB_fromF.Edges):
        idxsFs_InvalidEdges.append(iF)
    if bFaces and any(not rgX.IsValid for rgX in rgB_fromF.Faces):
        idxsFs_InvalidFaces.append(iF)
    if bLoops and any(not rgX.IsValid for rgX in rgB_fromF.Loops):
        idxsFs_InvalidLoops.append(iF)
    if bSrfs and any(not rgX.IsValid for rgX in rgB_fromF.Surfaces):
        idxsFs_InvalidSurfaces.append(iF)
    if bTrims and any(not rgX.IsValid for rgX in rgB_fromF.Trims):
        idxsFs_InvalidTrims.append(iF)
    if bVertices and any(not rgX.IsValid for rgX in rgB_fromF.Vertices):
        idxsFs_InvalidVertices.append(iF)

    rgB_fromF.Dispose()
    rgB_fromS.Dispose()


def indicesOfFacesOfInvalidBrepComponents(rgBrep, bFacesAsBreps=True, bSrfsAsBreps=True, bCrvs2D=True, bCrvs3D=True, bEdges=True, bFaces=True, bLoops=True, bSrfs=True, bTrims=True, bVertices=True, bDebug=False):
    """
    """

    idxsFs_InvalidFacesAsBreps = []
    idxsFs_InvalidSrfsAsBreps = []
    idxsFs_InvalidCrvs2d = []
    idxsFs_InvalidCrvs3d = []
    idxsFs_InvalidEdges = []
    idxsFs_InvalidFaces = []
    idxsFs_InvalidLoops = []
    idxsFs_InvalidSurfaces = []
    idxsFs_InvalidTrims = []
    idxsFs_InvalidVertices = []

    for iF in range(rgBrep.Faces.Count):
        rgF = rgBrep.Faces[iF]
        # Check duplicate brep from each face because face itself may show as being valid.
        rgB_fromF = rgF.DuplicateFace(duplicateMeshes=False)
        rgB_fromS = rgF.UnderlyingSurface().ToBrep()

        if bFacesAsBreps and not rgB_fromF.IsValid:
            idxsFs_InvalidFacesAsBreps.append(iF)
        if bSrfsAsBreps and not rgB_fromS.IsValid:
            idxsFs_InvalidFacesAsBreps.append(iF)
        if bCrvs2D and any(not rgX.IsValid for rgX in rgB_fromF.Curves2D):
            idxsFs_InvalidCrvs2d.append(iF)
        if bCrvs3D and any(not rgX.IsValid for rgX in rgB_fromF.Curves3D):
            idxsFs_InvalidCrvs3d.append(iF)
        if bEdges and any(not rgX.IsValid for rgX in rgB_fromF.Edges):
            idxsFs_InvalidEdges.append(iF)
        if bFaces and any(not rgX.IsValid for rgX in rgB_fromF.Faces):
            idxsFs_InvalidFaces.append(iF)
        if bLoops and any(not rgX.IsValid for rgX in rgB_fromF.Loops):
            idxsFs_InvalidLoops.append(iF)
        if bSrfs and any(not rgX.IsValid for rgX in rgB_fromF.Surfaces):
            idxsFs_InvalidSurfaces.append(iF)
        if bTrims and any(not rgX.IsValid for rgX in rgB_fromF.Trims):
            idxsFs_InvalidTrims.append(iF)
        if bVertices and any(not rgX.IsValid for rgX in rgB_fromF.Vertices):
            idxsFs_InvalidVertices.append(iF)

        rgB_fromF.Dispose()
        rgB_fromS.Dispose()

    return sorted(set(
        idxsFs_InvalidFacesAsBreps +
        idxsFs_InvalidSrfsAsBreps +
        idxsFs_InvalidCrvs2d +
        idxsFs_InvalidCrvs3d +
        idxsFs_InvalidEdges +
        idxsFs_InvalidFaces +
        idxsFs_InvalidLoops +
        idxsFs_InvalidSurfaces +
        idxsFs_InvalidTrims +
        idxsFs_InvalidVertices
        ))


def isValid(rhBrep):
    brep = coerceBrep(rhBrep)
    bValid = brep.IsValid
    brep.Dispose()
    return bValid


def extractBadFaces(rhBreps, bEcho=False, bDebug=False):
    """
    Returns: gBreps_NotValid, gBreps_Valid
    """

    try: rhBreps = list(rhBreps)
    except: rhBreps = [rhBreps]

    gBs_notValid = []
    gBs_Valid = []
    gBs_toExtract = []

    for rhBrep in rhBreps:
        rdBrep = coerceRhinoObject(rhBrep)
        rgBrep = rdBrep.BrepGeometry
        if rgBrep.IsValid:
            gBs_Valid.append(rdBrep.Id)
        elif rgBrep.Faces.Count == 1:
            gBs_notValid.append(rdBrep.Id)
        else:
            gBs_toExtract.append(rdBrep.Id)
        rgBrep.Dispose() # Geometry won't be used for remainder of function.

    if not gBs_toExtract:
        return gBs_notValid, gBs_Valid


    # Record first normal object (last created) that is not a brep passed to this function.
    gFirstObj_Start = None
    for rdObj in list(sc.doc.Objects.GetObjectList(rd.ObjectType.AnyObject)): # EnumeratorWrappe -> list
        if rdObj.Id not in gBs_notValid + gBs_Valid + gBs_toExtract:
            gFirstObj_Start = rdObj
            break

    if bDebug:
        sPrint = 'gFirstObj_Start'; print(sPrint + ':', eval(sPrint))
    else:
        sc.doc.Views.RedrawEnabled = False

    # Extract bad faces using UI command.
    sc.doc.Objects.UnselectAll()
    [sc.doc.Objects.Select(g) for g in gBs_toExtract]
    if bEcho:
        Rhino.RhinoApp.RunScript("_ExtractBadSrf", echo=True)
    else:
        Rhino.RhinoApp.RunScript("_NoEcho _ExtractBadSrf _Echo", echo=False)
    # Only invalid monoface breps are selected:
    rds_notValid_Extracted = list(sc.doc.Objects.GetSelectedObjects(False, False))
    
    # Create a list of remainder of breps.
    rdObjs = sc.doc.Objects.GetObjectList(rd.ObjectType.AnyObject)
    rdBs_Valid_fromExtr = []
    for rdObj in rdObjs:
        if (
                (gFirstObj_Start is not None) and
                (rdObj.Id == gFirstObj_Start)
        ):
            # Past modified brep.
            break
        if rdObj.IsSelected(False): continue # Brep is bad.  (False is for checkSubObjects.)
        rdBs_Valid_fromExtr.append(rdObj)
    
    sc.doc.Objects.UnselectAll()

    if not bDebug: sc.doc.Views.RedrawEnabled = True

    return (
        gBs_notValid + [br.Id for br in rds_notValid_Extracted],
        gBs_Valid + [br.Id for br in rdBs_Valid_fromExtr]
    )


def dotAtSurfaceCentroid(rgSrf, text='!', iDotHeight=14):
    ptCentrdW = (
            Rhino.Geometry.AreaMassProperties.Compute(rgSrf).Centroid)
    getrc, u, v = rgSrf.ClosestPoint(ptCentrdW)
    ptCentroid = rgSrf.PointAt(u, v)
    rgDot = rg.TextDot(text, ptCentroid)
    rgDot.FontHeight = iDotHeight
    return rgDot


def processBrepObjects(rhBreps, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bFacesAsBreps = getOpt('bFacesAsBreps')
    bSrfsAsBreps = getOpt('bSrfsAsBreps')
    bCrvs2D = getOpt('bCrvs2D')
    bCrvs3D = getOpt('bCrvs3D')
    bEdges = getOpt('bEdges')
    bFaces = getOpt('bFaces')
    bLoops = getOpt('bLoops')
    bSrfs = getOpt('bSrfs')
    bTrims = getOpt('bTrims')
    bVertices = getOpt('bVertices')
    bAnalyzePerFace = getOpt('bAnalyzePerFace')
    bExtract = getOpt('bExtract')
    bAddDot = getOpt('bAddDot')
    iDotHeight = getOpt('iDotHeight')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    
    rgDots = []
    gBreps0_Pos = []
    iCt_InvalidFaces_All = 0
    gBreps_NotValid_perBrep0 = []
    gBreps_Valid_perBrep0 = []
    
    Rhino.RhinoApp.CommandPrompt = "Searching ..."

    if bAnalyzePerFace:
        for rhBrep in rhBreps:
            rdBrep0 = coerceRhinoObject(rhBrep)
            rgBrep0 = coerceBrep(rdBrep0)
            if rgBrep0 is None: continue
        
            # Get indices of invalid brep faces.
            idx_rgFaces_Invalid = indicesOfFacesOfInvalidBrepComponents(rgBrep0)
            if not idx_rgFaces_Invalid: continue
        
            iCt_InvalidFaces_All = len(idx_rgFaces_Invalid)
            if iCt_InvalidFaces_All == 0:
                # No matches in this brep
                continue
        
            if bExtract:
                if rgBrep0.Faces.Count == 1: # Single-face brep
                    gBreps_NotValid_perBrep0.append([rdBrep0.Id])
                    continue
            
                # Extract from brep with multiple faces.  Uses Rhino command.
                rc = extractBadFaces(rdBrep0)
                if rc is None: continue
                gBreps_NotValid, gBreps_Valid = rc
                gBreps_NotValid_perBrep0.append(gBreps_NotValid)
                gBreps_Valid_perBrep0.append(gBreps_Valid)
            else:
                gBreps0_Pos.append(rdBrep0.Id)
        
            if bAddDot:
                for f in idx_rgFaces_Invalid:
                    rgFace = rgBrep0.Faces[f]
                    rgDots.append(dotAtSurfaceCentroid(
                        rgFace, iDotHeight=iDotHeight))
                    rgFace.Dispose()
        
            # End of idBreps0 loop.
    
        if len(rgDots) > 0:
            attr = rd.ObjectAttributes()
            attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
            attr.ColorSource = rd.ObjectColorSource.ColorFromObject
            attr.ObjectColor = Color.FromArgb(255,0,0)
            for rgDot in rgDots:
                sc.doc.Objects.AddTextDot(rgDot, attr)
                rgDot.Dispose()
    
        if not iCt_InvalidFaces_All:
            print("No faces with invalid content found.")
            return

        print("{} faces with invalid content found.".format(iCt_InvalidFaces_All))
        if bExtract:
            sc.doc.Objects.UnselectAll()
            print("{} extracted face(s) are selected.".format(
                sc.doc.Objects.Select(List[Guid](
                    [g for item in gBreps_NotValid_perBrep0 for g in item]))))
            return gBreps_NotValid_perBrep0, gBreps_Valid_perBrep0
        else:
            gBreps0_Pos = List[Guid](gBreps0_Pos)
            sc.doc.Objects.UnselectAll()
            print("{} brep(s) are selected.".format(
                sc.doc.Objects.Select(gBreps0_Pos)))
            return gBreps0_Pos
    else:
        for rhBrep in rhBreps:
            rdBrep0 = coerceRhinoObject(rhBrep)
            rgBrep0 = coerceBrep(rdBrep0)
            if rgBrep0 is None: continue
            sReport = reportOfInvalidComponentCounts(
                rgBrep0,
                bFacesAsBreps=bFacesAsBreps,
                bSrfsAsBreps=bSrfsAsBreps,
                bCrvs2D=bCrvs2D,
                bCrvs3D=bCrvs3D,
                bEdges=bEdges,
                bFaces=bFaces,
                bLoops=bLoops,
                bSrfs=bSrfs,
                bTrims=bTrims,
                bVertices=bVertices,
                bDebug=bDebug)
            print(sReport)


def main():

    rhBreps = getInput()
    if not rhBreps: return

    bFacesAsBreps = Opts.values['bFacesAsBreps']
    bSrfsAsBreps = Opts.values['bSrfsAsBreps']
    bCrvs2D = Opts.values['bCrvs2D']
    bCrvs3D = Opts.values['bCrvs3D']
    bEdges = Opts.values['bEdges']
    bFaces = Opts.values['bFaces']
    bLoops = Opts.values['bLoops']
    bSrfs = Opts.values['bSrfs']
    bTrims = Opts.values['bTrims']
    bVertices = Opts.values['bVertices']
    bAnalyzePerFace = Opts.values['bAnalyzePerFace']
    bExtract = Opts.values['bExtract']
    bAddDot = Opts.values['bAddDot']
    iDotHeight = Opts.values['iDotHeight']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    processBrepObjects(
        rhBreps,
        bFacesAsBreps=bFacesAsBreps,
        bSrfsAsBreps=bSrfsAsBreps,
        bCrvs2D=bCrvs2D,
        bCrvs3D=bCrvs3D,
        bEdges=bEdges,
        bFaces=bFaces,
        bLoops=bLoops,
        bSrfs=bSrfs,
        bTrims=bTrims,
        bVertices=bVertices,
        bAnalyzePerFace=bAnalyzePerFace,
        bExtract=bExtract,
        bAddDot=bAddDot,
        iDotHeight=iDotHeight,
        bEcho=bEcho,
        bDebug=bDebug)


    sc.doc.Views.Redraw()


if __name__ == '__main__': main()