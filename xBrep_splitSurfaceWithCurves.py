"""
190530: Created.
191112: Added a function.
191113: Made input of some function more BrepFace focused.
191115: Added more tolerance multipliers.
200130, 0202: Modified output of a function.
200209, 14: Modified debug output.
200303: Now accepts curves with length < resolution and >= ModelAbsoluteTolerance.
200422: Refactored Opts.  Added some options.  Modified an option default value.
200430: Added bOnlyUseCrvsOnSrf and option to use all normal wires and brep edges.
200505: Bug fix.
200519-23, 26: Refactored.  Exported a function to its own script.
200611: Bug fix.
200619: Import-related update.
200625: Now Normal objects are now found in main() instead of a subsequent function.
200629: Refactored.
200701: Bug fix.
200729: Import-related update.
210122: Now, processBrepObjects processes more object types for faces and splitters.
210125: Minor bug fix.
210325: Added bAddSplittingCrvs.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

import xCurve
import xCurve_sortCurvesOnSurface
import xCurve_splitWithSurfaceEdges
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


    key = 'bScanForNakedEdges'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Off', 'On')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSelBrep'; keys.append(key)
    values[key] = True
    names[key] = 'BrepSelMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Edge', 'Brep')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitUnderlyingSrf'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOnlyUseCrvsOnSrf'; keys.append(key)
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

    key = 'bExplode'; keys.append(key)
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

    key = 'bAddSplittingCrvs'; keys.append(key)
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

    go.SetCommandPrompt("Select faces to split")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(enable=True, acceptZero=True)
    
    idxs_Opts = {}

    while True:
        key = 'bSplitUnderlyingSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bOnlyUseCrvsOnSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bSplitToCrvSegs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bTryOtherTolsOnFail'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if not Opts.values['bSplitUnderlyingSrf']:
            key = 'bExtract'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if Opts.values['bSplitUnderlyingSrf'] or Opts.values['bExtract']:
            key = 'bExplode'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if Opts.values['bDebug']:
            key = 'bAddSplittingCrvs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
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


def getInput_TrimmingObjects(objrefs_Face):
    """
    Get curves or edges with optional input.
    """


    def getBrepIdsAndFaceEdgeIdxs(objrefs_Face):
        gBs = []
        idxs_Fs_perB = []

        for objref_Face in objrefs_Face:
            gBrep_withFaceToSplit = rs.coerceguid(objref_Face)
            if gBrep_withFaceToSplit in gBs:
                iB = gBs.index(gBrep_withFaceToSplit)
                idxs_Fs_perB[iB].extend(list(objref_Face.Face().AdjacentEdges()))
            else:
                gBs.append(gBrep_withFaceToSplit)
                idxs_Fs_perB.append(list(objref_Face.Face().AdjacentEdges()))

        return gBs, idxs_Fs_perB


    gBs_withFsToSplit, idxs_Fs_perB = getBrepIdsAndFaceEdgeIdxs(objrefs_Face)

    go = ri.Custom.GetObject()

    def notEdgeOfFaceToSplit(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx
        if compIdx.ComponentIndexType == rg.ComponentIndexType.BrepEdge:
            if isinstance(rdObj, rd.BrepObject):
                if rdObj.Id in gBs_withFsToSplit:
                    iB = gBs_withFsToSplit.index(rdObj.Id)
                    if geom.EdgeIndex in idxs_Fs_perB[iB]:
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
        if Opts.values['bScanForNakedEdges']:
            go.SetCommandPromptDefault("Enter for all normal wires and brep naked edges")
        else:
            go.SetCommandPromptDefault("Enter for all normal wires")

        key = 'bScanForNakedEdges'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bSelBrep'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bSplitUnderlyingSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bOnlyUseCrvsOnSrf'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bSplitToCrvSegs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bTryOtherTolsOnFail'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if not Opts.values['bSplitUnderlyingSrf']:
            key = 'bExtract'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if Opts.values['bSplitUnderlyingSrf'] or Opts.values['bExtract']:
            key = 'bExplode'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if Opts.values['bDebug']:
            key = 'bAddSplittingCrvs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]

        if Opts.values['bSelBrep']:
            go.GeometryFilter = rd.ObjectType.Curve | rd.ObjectType.Brep
            go.SubObjectSelect = False
            go.SetCommandPrompt("Select curves or breps")
        else:
            go.GeometryFilter = rd.ObjectType.Curve
            go.SubObjectSelect = True
            go.SetCommandPrompt("Select curves or edges")

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            return
        elif res == ri.GetResult.Nothing:
            go.Dispose()
            return tuple([[]] + [Opts.values[key] for key in Opts.keys])
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


def getFormattedDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def getNakedEdgesOfBrepObject(rdBrep):
    """Using BrepObject instead of rg.Brep to use GUID to avoid brep to split."""

    rgEs = []
    rgBrep = rdBrep.BrepGeometry


    for edge in rgBrep.Edges:
        if edge.Valence == rg.EdgeAdjacency.Naked:
            rgEs.append(edge)


    #
    # Keep this in case filtering of brep to face level is reimplemented.
    #
    ## Structure of following conditional is to optimize speed.
    #if rdBrep.Id not in gBs_withFsToSplit:
    #    # Get all naked edges without restriction.
    #    for edge in rgBrep.Edges:
    #        if edge.Valence == rg.EdgeAdjacency.Naked:
    #            rgEs.append(edge)
    #else:
    #    iB = gBs_withFsToSplit.index(rdBrep.Id)
    #    for edge in rgBrep.Edges:
    #        if edge.EdgeIndex not in idxs_Fs_perB[iB]:
    #            if edge.Valence == rg.EdgeAdjacency.Naked:
    #                rgEs.append(edge)


    return rgEs


def splitSurfaceIntoBrep(rgSrf_toSplit, rgCrvs_Splitters, **kwargs):
    """
    Parameters:
        rgSrf_toSplit: Can be rg.BrepFace or other rg.Surface.
        rgCrvs_Splitters
        fTolerance
        bTryOtherTolsOnFail
        bDebug
    Returns on success:
        rg.Brep
    Returns on fail:
        None
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    bTryOtherTolsOnFail = getOpt('bTryOtherTolsOnFail')
    bDebug = getOpt('bDebug')


    if isinstance(rgSrf_toSplit, rg.BrepFace):
        rgFace_toSplit = rgSrf_toSplit
        rgBrep_TempForUnderlyingSrf = None
    elif isinstance(rgSrf_toSplit, rg.Surface):
        rgBrep_TempForUnderlyingSrf = rgSrf_toSplit.ToBrep()
        rgFace_toSplit = rgBrep_TempForUnderlyingSrf.Faces[0]
    else:
        return


    # Create tolerance loop.
    if bTryOtherTolsOnFail:
        fTols_toTry = []
        for tolMultiplier in 1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0, 0.0625, 16.0:
            tol = tolMultiplier * fTolerance
            fTols_toTry.append(tol)
    else:
        fTols_toTry = fTolerance,


    # Split in a tolerance loop.
    for fTol_toTry in fTols_toTry:

        #
        #
        rgB_Split = rgFace_toSplit.Split(
            curves=rgCrvs_Splitters,
            tolerance=fTol_toTry)
        #
        #


        if bDebug: sEval='rgB_Split'; print sEval+':',eval(sEval)
        
        if rgB_Split is None:
            if bDebug:
                print "  Failed at fTol_toTry=={}.".format(
                    getFormattedDistance(fTol_toTry))
                for c in crvs:
                    sc.doc.Objects.AddCurve(c)
                sc.doc.Views.Redraw()
        elif rgB_Split.Faces.Count == 1:
            if bDebug:
                print "BrepFace.Split resulted in a 1-face Brep."
            rgB_Split = None
        elif rgB_Split.Faces.Count == rgFace_toSplit.Brep.Faces.Count:
            if bDebug:
                print "BrepFace.Split resulted in a Brep with no additional faces."
            rgB_Split = None
        else:
            if bDebug:
                sEval='rgB_Split.IsValid'; print sEval+':',eval(sEval)
                sEval='rgB_Split.Faces.Count'; print sEval+':',eval(sEval)
                #sc.doc.Objects.AddBrep(rgB_Split); sc.doc.Views.Redraw()
            if not rgB_Split.IsValid:
                rgB_Split = None

            if bDebug or abs(fTol_toTry - fTolerance) > 1e-9:
                print "  Split successful at a tolerance of {}.".format(
                    getFormattedDistance(fTol_toTry))
            break # out of tolerance loop.


    if rgBrep_TempForUnderlyingSrf: rgBrep_TempForUnderlyingSrf.Dispose()

    return rgB_Split


def splitFaceOfBrep(rgFace_toSplit, rgCrvs_Splitters, **kwargs):
    """
    Parameters:
        rgSrf_toSplit: Can be rg.BrepFace or other rg.Surface.
        rgCrvs_Splitters
        fTolerance
        bTryOtherTolsOnFail
        bDebug
    Returns on success:
        rg.Brep
    Returns on fail:
        None
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    bTryOtherTolsOnFail = getOpt('bTryOtherTolsOnFail')
    bDebug = getOpt('bDebug')


    if not isinstance(rgFace_toSplit, rg.BrepFace):
        return

    if not rgCrvs_Splitters:
        return

    #for c in rgCrvs_Splitters: sc.doc.Objects.AddCurve(c)

    # Create tolerance loop.
    if bTryOtherTolsOnFail:
        fTols_toTry = []
        for tolMultiplier in 1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0, 0.0625, 16.0:
            tol = tolMultiplier * fTolerance
            fTols_toTry.append(tol)
    else:
        fTols_toTry = fTolerance,


    # Split in a tolerance loop.
    for fTol_toTry in fTols_toTry:

        #
        #
        rgB_Split = rgFace_toSplit.Split(
            curves=rgCrvs_Splitters,
            tolerance=fTol_toTry)
        #
        #


        if bDebug: sEval='rgB_Split'; print sEval+':',eval(sEval)
        
        if rgB_Split is None:
            if bDebug:
                print "  Failed at fTol_toTry=={}.".format(
                    getFormattedDistance(fTol_toTry))
                for c in crvs:
                    sc.doc.Objects.AddCurve(c)
                sc.doc.Views.Redraw()
        elif rgB_Split.Faces.Count == rgFace_toSplit.Brep.Faces.Count:
            if bDebug:
                print "BrepFace.Split resulted in a Brep with no additional faces."
            rgB_Split = None
        else:
            if bDebug:
                sEval='rgB_Split.IsValid'; print sEval+':',eval(sEval)
                sEval='rgB_Split.Faces.Count'; print sEval+':',eval(sEval)
                #sc.doc.Objects.AddBrep(rgB_Split); sc.doc.Views.Redraw()
            if not rgB_Split.IsValid:
                rgB_Split = None

            if bDebug or abs(fTol_toTry - fTolerance) > 1e-9:
                print "  Split successful at a tolerance of {}.".format(
                    getFormattedDistance(fTol_toTry))
            break # out of tolerance loop.


    rgFace_toSplit.Brep.Dispose()

    return rgB_Split


def getCurvesToSplitSurface(rgCrvs_tryForSplitters, rgSrf, fTolerance, bDebug=False):
    """
    Duplicates input curves that are on surface and splits to edges.
    """
    rc = xCurve_sortCurvesOnSurface.duplicateCurvesOnSurface(
        rgCrvs_tryForSplitters,
        rgSrf,
        fTolerance=fTolerance,
        bDebug=bDebug)
    if not rc: return
    (
        crvs_CompletelyOnSrf_Closed,
        crvs_CompletelyOnSrf_Open,
        crvs_PartiallyOnSrf,
        ) = rc

    if not crvs_PartiallyOnSrf:
        return (crvs_CompletelyOnSrf_Closed + crvs_CompletelyOnSrf_Open)

    # Trim the curves that lie partially on the surface mostly to get rid of the overlapping (Boundary) portions.
    rc = xCurve_splitWithSurfaceEdges.splitCurvesWithSurfaceEdges(
        rgCrvs_In=crvs_PartiallyOnSrf,
        rgSrf_In=rgSrf,
        fTolerance=fTolerance)
    (
        crvs_Trimmed_Interior,
        crvs_Trimmed_Boundary,
        crvs_Trimmed_Exterior
        ) = rc

    for c in crvs_Trimmed_Boundary + crvs_Trimmed_Exterior: c.Dispose()
    for c in crvs_PartiallyOnSrf: c.Dispose()

    return (
        crvs_CompletelyOnSrf_Closed +
        crvs_CompletelyOnSrf_Open +
        crvs_Trimmed_Interior
        )


def processBrepObjects(rhObjs_Faces, rhObjs_Splitters, **kwargs):
    """
    Parameters:
        rhObjs_Faces,
        rhObjs_Splitters,
        bSplitUnderlyingSrf
        bOnlyUseCrvsOnSrf
        bSplitToCrvSegs
        fTolerance
        bTryOtherTolsOnFail
        bExplode
        bExtract
        bEcho
        bDebug
        bAddSplittingCrvs
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bSplitUnderlyingSrf = getOpt('bSplitUnderlyingSrf')
    bOnlyUseCrvsOnSrf = getOpt('bOnlyUseCrvsOnSrf')
    bSplitToCrvSegs = getOpt('bSplitToCrvSegs')
    fTolerance = getOpt('fTolerance')
    bTryOtherTolsOnFail = getOpt('bTryOtherTolsOnFail')
    bExplode = getOpt('bExplode')
    bExtract = getOpt('bExtract')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    bAddSplittingCrvs = getOpt('bAddSplittingCrvs')


    def getBrepObject(rhObj):
        rdObj = rs.coercerhinoobject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def getSortedBrepIdsAndFaces(rhObjs_Faces):
        """
        Parameters:
            list(rhObjs_Faces)
        Returns:
            list(Brep GUIDs)
            list(lists(integers of Face indices) per brep)
        """
        
        gBreps_In = []
        idxs_Faces_perBrep = []

        for rhObj in rhObjs_Faces:
            if isinstance(rhObj, Guid):
                gBrep_In = rhObj
                rgBrep_In = rs.coercebrep(gBrep_In)
                if gBrep_In in gBreps_In:
                    idxs_Faces_perBrep[gBreps_In.index(gBrep_In)] = range(rgBrep_In.Faces.Count)
                else:
                    gBreps_In.append(gBrep_In)
                    idxs_Faces_perBrep.append(range(rgBrep_In.Faces.Count))
                continue

            if not isinstance(rhObj, rd.ObjRef):
                raise ValueError(
                    "This {} is invalid input for processBrepObjects.".format(
                        rhObj.GetType().Name))

            o = rhObj # ObjRef.

            gBrep_In = o.ObjectId
            rdBrep_In = o.Object()
            rgBrep_In = o.Brep()
        
            if not rgBrep_In.IsValid:
                print "Brep {} is invalid.  Fix first.".format(gBrep_In)
                rgBrep_In.Dispose()
                continue
        
            idx_CompIdx = o.GeometryComponentIndex.Index
            if idx_CompIdx == -1:
                if gBrep_In in gBreps_In:
                    idxs_Faces_perBrep[gBreps_In.index(gBrep_In)] = range(rgBrep_In.Faces.Count)
                else:
                    gBreps_In.append(gBrep_In)
                    idxs_Faces_perBrep.append(range(rgBrep_In.Faces.Count))
            else:
                rgFace_Brep0 = o.Face()
                if gBrep_In in gBreps_In:
                    if rgFace_Brep0 in idxs_Faces_perBrep[gBreps_In.index(gBrep_In)]:
                        continue
                    else:
                        idxs_Faces_perBrep[gBreps_In.index(gBrep_In)].append(rgFace_Brep0.FaceIndex)
                else:
                    gBreps_In.append(gBrep_In)
                    idxs_Faces_perBrep.append([rgFace_Brep0.FaceIndex])

        return gBreps_In, idxs_Faces_perBrep


    def getAllNormalObjsForSplitters():
        rgCrvs_CrvObjs = []
        rgCrvs_Edges = []
        gBreps_ofEdges = []
        settings = rd.ObjectEnumeratorSettings()
        settings.NormalObjects = True
        settings.LockedObjects = False
        for rdObj in sc.doc.Objects.GetObjectList(settings):
            if rdObj.ObjectType == rd.ObjectType.Curve:
                rgCrvs_CrvObjs.append(rdObj.Geometry)
            elif rdObj.ObjectType == rd.ObjectType.Brep:
                rgNEs = getNakedEdgesOfBrepObject(rdObj)
                for rgNE in rgNEs:
                    rgCrvs_Edges.append(rgNE.DuplicateCurve())
                    gBreps_ofEdges.append(rdObj.Id)
        return rgCrvs_CrvObjs, rgCrvs_Edges, gBreps_ofEdges


    def getSelctdObjsForSplitters(rhObjs_CurveOrBrep):
        rgCrvs_CrvObjs = []
        rgCrvs_Edges = []
        gBreps_ofEdges = []
        for o in rhObjs_CurveOrBrep:
            if isinstance(o, rg.Curve):
                rgCrvs_CrvObjs.append(o)
                continue
            rdObj = rs.coercerhinoobject(o) # coercecurve returns various rg.Curves, including rg.BrepEdge.
            if rdObj.ObjectType == rd.ObjectType.Curve:
                rgCrvs_CrvObjs.append(rdObj.Geometry)
            elif rdObj.ObjectType == rd.ObjectType.Brep:
                rgNEs = getNakedEdgesOfBrepObject(rdObj)
                for rgNE in rgNEs:
                    rgCrvs_Edges.append(rgNE.DuplicateCurve())
                    gBreps_ofEdges.append(rdObj.Id)
        return rgCrvs_CrvObjs, rgCrvs_Edges, gBreps_ofEdges


    def splitIndividualFaceOrSurface():

        gBs_Added = []

        rgCrvs_Splitters_FilteredToBrep = (
            rgCrvs_fromCrvObjs +
            [c for c, g in zip(rgCrvs_fromEdges, gBreps_ofEdges) if g != gBrep_In])

        for iF, idxFace in enumerate(idxs_rgFace_perBrep[iB]):
            rgFace_In = rgBrep_In.Faces[idxFace]

            rgBrep_1Face_B0 = rgFace_In.DuplicateFace(duplicateMeshes=False)

            if bSplitUnderlyingSrf:
                srf_toSplit = rgBrep_1Face_B0.Faces[0].UnderlyingSurface()
            else:
                srf_toSplit = rgBrep_1Face_B0.Faces[0]


            if bOnlyUseCrvsOnSrf:
                rgCrvs_Splitters_thisFace = getCurvesToSplitSurface(
                    rgCrvs_Splitters_FilteredToBrep, srf_toSplit, fTolerance, bDebug)
                if bDebug and bAddSplittingCrvs:
                    for c in rgCrvs_Splitters_thisFace: sc.doc.Objects.AddCurve(c)
            else:
                rgCrvs_Splitters_thisFace = rgCrvs_Splitters_FilteredToBrep


            if bSplitToCrvSegs:
                cs_WIP = xCurve.duplicateSegments(
                    rgCrvs_Splitters_thisFace,
                    bExplodePolyCrvs=True)
                for c in rgCrvs_Splitters_thisFace:
                    c.Dispose()
                rgCrvs_Splitters_thisFace = cs_WIP
                if bDebug and bAddSplittingCrvs:
                    for c in rgCrvs_Splitters_thisFace: sc.doc.Objects.AddCurve(c)


            rgBrep_Ret_Split = splitSurfaceIntoBrep(
                rgSrf_toSplit=srf_toSplit,
                rgCrvs_Splitters=rgCrvs_Splitters_thisFace,
                fTolerance=fTolerance,
                bTryOtherTolsOnFail=bTryOtherTolsOnFail,
                bDebug=bDebug)
            rgBrep_1Face_B0.Dispose()

            if not rgBrep_Ret_Split: continue

            rgBs_1F_Ret.append(rgBrep_Ret_Split)
            idxFaces_SplitSuccess.append(idxFace)

        idxFaces_AddBrepSuccess = []

        for i, rgB in enumerate(rgBs_1F_Ret):
            bSomeFacesNotAdded = False
            if bExplode:
                for iF, rgF in enumerate(rgB.Faces):
                    rgB_1F = rgF.DuplicateFace(True)
                    gB = sc.doc.Objects.AddBrep(
                        rgB_1F,
                        attributes=rdBrep_In.Attributes)
                    rgB_1F.IsDocumentControlled
                    rgB_1F.Dispose()
                    if gB == Guid.Empty:
                        bSomeFacesNotAdded = True
                        continue
                    gBs_Added.append(gB)
                if bSomeFacesNotAdded:
                    print "Some split faces were not added.  Check results."
                else:
                    idxFaces_AddBrepSuccess.append(idxFaces_SplitSuccess[i])
            else:
                gB = sc.doc.Objects.AddBrep(
                    rgB,
                    attributes=rdBrep_In.Attributes)
                if gB == Guid.Empty:
                    continue
                gBs_Added.append(gB)
                idxFaces_AddBrepSuccess.append(idxFaces_SplitSuccess[i])

        if idxFaces_AddBrepSuccess:
            if rgBrep_In.Faces.Count == 1:
                sc.doc.Objects.Delete(objectId=rdBrep_In.Id, quiet=False)
            else:
                gBs_Remaining = xBrepObject.removeFaces(
                    rdBrep_In, idxFaces_AddBrepSuccess)
                # Remove orginal faces.
                if not gBs_Remaining:
                    print "Remove faces fail!"
                else:
                    print "{} breps remaining from extracted brep.".format(
                        len(gBs_Remaining))

        return gBs_Added


    def splitInContextOfEntireBrep():
        # rgBrep_Split_WIP is for a full brep, containing both
        # the new faces from Split
        # and faces that are not supposed to be split.
        rgBrep_Split_WIP = rgBrep_In.Duplicate()

        rgCrvs_Splitters_FilteredToBrep = (
            rgCrvs_fromCrvObjs +
            [c for c, g in zip(rgCrvs_fromEdges, gBreps_ofEdges) if g != gBrep_In])

        for iF, idxFace in enumerate(idxs_rgFace_perBrep[iB]):
            rgFace_toSplit = rgBrep_Split_WIP.Faces[idxFace]


            if bOnlyUseCrvsOnSrf:
                rgCrvs_Splitters_thisFace = getCurvesToSplitSurface(
                    rgCrvs_Splitters_FilteredToBrep, rgFace_toSplit, fTolerance, bDebug)
                if bDebug and bAddSplittingCrvs:
                    for c in rgCrvs_Splitters_thisFace: sc.doc.Objects.AddCurve(c)
            else:
                rgCrvs_Splitters_thisFace = rgCrvs_Splitters_FilteredToBrep


            if bSplitToCrvSegs:
                cs_WIP = xCurve.duplicateSegments(
                    rgCrvs_Splitters_thisFace,
                    bExplodePolyCrvs=True)
                rgCrvs_Splitters_thisFace = cs_WIP


            rgBrep_Ret_Split = splitFaceOfBrep(
                    rgFace_toSplit=rgFace_toSplit,
                    rgCrvs_Splitters=rgCrvs_Splitters_thisFace,
                    bOnlyUseCrvsOnSrf=bOnlyUseCrvsOnSrf,
                    fTolerance=fTolerance,
                    bTryOtherTolsOnFail=bTryOtherTolsOnFail,
                    bDebug=bDebug)

            if rgBrep_Ret_Split:
                rgBrep_Split_WIP.Dispose()
                rgBrep_Split_WIP = rgBrep_Ret_Split
                idxFaces_SplitSuccess.append(idxFace)

        if not rgBrep_Ret_Split:
            return

        # Success.
        
        #xBrepObject.replaceFaces(
        #    gBrep_In,
        #    idxFaces_SplitSuccess,
        #    

        bReplaced = sc.doc.Objects.Replace(
            objectId=gBrep_In, brep=rgBrep_Ret_Split)
        
        if not bReplaced:
            print "Brep update fail for {}.".format(gBrep_In)
            rgBrep_Ret_Split.Dispose()
            return
        
        if bEcho:
            s  = "Brep was updated."
            s += "  Face count increased by {} ({} -> {}).".format(
                rgBrep_Ret_Split.Faces.Count-ct_Face_B0,
                ct_Face_B0,
                rgBrep_Ret_Split.Faces.Count)
            print s
        
        rgBrep_Ret_Split.Dispose()
        
        return [gBrep_In]


    gBreps_In, idxs_rgFace_perBrep = getSortedBrepIdsAndFaces(rhObjs_Faces)
    if not gBreps_In: return

    if not rhObjs_Splitters:
        rc = getAllNormalObjsForSplitters()
    else:
        rc = getSelctdObjsForSplitters(rhObjs_Splitters)
    if not rc: return

    (
        rgCrvs_fromCrvObjs,
        rgCrvs_fromEdges,
        gBreps_ofEdges,
        ) = rc


    gBs_Out = []

    if bSplitUnderlyingSrf:
        bExtract = True


    for iB, gBrep_In in enumerate(gBreps_In):
        rdBrep_In = getBrepObject(gBrep_In)
        rgBrep_In = rdBrep_In.Geometry

        idxFaces_SplitSuccess = []
        rgBs_1F_Ret = []

        ct_Face_B0 = rgBrep_In.Faces.Count

        if bExtract and ct_Face_B0 > 1:
            rc = splitIndividualFaceOrSurface()
        else:
            rc = splitInContextOfEntireBrep()

        if rc:
            gBs_Out.extend(rc)


    return gBs_Out


def main():
    
    rc = getInput_Faces()
    if rc is None: return
    objrefs_Faces = rc[0]

    sc.doc.Objects.UnselectAll()

    rc = getInput_TrimmingObjects(objrefs_Faces)
    if rc is None: return
    rhObjs_Splitters = rc[0]

    if not rhObjs_Splitters:
        # Get all Normal objects now so before any layer, etc., changes occur during next input.
        rhObjs_Splitters = getAllNormalObjs(Opts.values['bScanForNakedEdges'])

    if not rhObjs_Splitters:
        print "No splitters."
        return

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gBreps_Result = processBrepObjects(
        rhObjs_Faces=objrefs_Faces,
        rhObjs_Splitters=rhObjs_Splitters,
        )

    if Opts.values['bEcho']:
        if not gBreps_Result:
            print "No breps were split."
        else:
            print "{} new/modified breps resulted from splits.".format(len(gBreps_Result))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()