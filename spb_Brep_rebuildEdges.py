"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160713-14: Created.
...
200204-05: Added bRebuildOnlyOutOfTol and routines to more selectively rebuild edges.
200213: Bug fix.
200507: Improved handling of breps with edges with error in value of Tolerance.
200625: Disabled use of Brep.Standardize due to problems with short edges.
201211: Added check to maximum allowed tolerance before increasing working tolerance.
210122: Bug fix: A variable was not set before use.
210426: Added fSearchTol.
210630: Modified an option default value.  Added fRebuildSharedTol.
250224: Disabled Standardize and Compact. Modified an option default value.

TODO:
    Investigate problem with brep as demonstrated by use of Brep.Standardize.
    Standardizing Brep sometimes corrects Edge.Tolerance.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


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


    key = 'fRebuildNakedTol'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bRebuildOnlyOutOfTol'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSearchTol'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bRebuildSharedEdges'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fRebuildSharedTol'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bRebuildVertices'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='Add', onValue='Replace')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHeight'; keys.append(key)
    values[key] = 11
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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


def _get_all_normal_breps():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))


def getInput():
    """
    Get breps and options values.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal")

    go.AcceptNothing(True)

    go.GeometryFilter = rd.ObjectType.Brep
    #go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.OpenSurface

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False
    go.SubObjectSelect = False

    go.EnableClearObjectsOnEntry(False)
    #go.EnableUnselectObjectsOnExit(False)

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    print("SearchTol is also the rebuild tolerance for shared edges.")

    while True:
        key = 'fRebuildNakedTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bRebuildOnlyOutOfTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        if Opts.values['bRebuildOnlyOutOfTol']:
            key = 'fSearchTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bRebuildSharedEdges'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        if Opts.values['bRebuildSharedEdges']:
            key = 'fRebuildSharedTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bRebuildVertices'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bReplace'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bAddDot'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'iDotHeight'; idxs_Opts[key] = (Opts.riAddOpts[key](go)
                                              if Opts.values['bAddDot']
                                              else None)
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return _get_all_normal_breps()

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected or a number was entered.
        key = 'fRebuildNakedTol'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        elif Opts.riOpts[key].CurrentValue == 0.0:
            Opts.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        key = 'fSearchTol'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        elif Opts.riOpts[key].CurrentValue == 0.0:
            Opts.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        key = 'fRebuildSharedTol'
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        elif Opts.riOpts[key].CurrentValue == 0.0:
            Opts.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def setTolerancesBoxesAndFlags(rgBrep):
    rgBrep.SetTolerancesBoxesAndFlags(
        bLazy=False,
        bSetVertexTolerances=True,
        bSetEdgeTolerances=True,
        bSetTrimTolerances=True,
        bSetTrimIsoFlags=True,
        bSetTrimTypeFlags=True,
        bSetLoopTypeFlags=True,
        bSetTrimBoxes=True)


def _formatDistance(fDistance):
    if fDistance is None:
        return "(None)"
    if fDistance == Rhino.RhinoMath.UnsetValue:
        return "(Infinite)"
    if fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-1)):
        return "{:.2e}".format(fDistance)
    return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def _compute_max_edge_tol_of_face(rgFace, fRebuildNakedTol=None, bRebuildVertices=None, bEcho=True):
    """
    Computes by comparing RebuildEdges with smaller tolerance to the current Edge.Tolerances.
    """

    rgB_1F_NotRebuilt = rgFace.DuplicateFace(False)
    rgB_1F_Rebuilt = rgB_1F_NotRebuilt.DuplicateBrep()

    tolerance = min((0.1*fRebuildNakedTol, 1e-5))

    rgB_1F_Rebuilt.Faces[0].RebuildEdges(
        tolerance=tolerance,
        rebuildSharedEdges=True,
        rebuildVertices=bRebuildVertices)

    deviations = []

    tolerance = min((0.01*fRebuildNakedTol, 1e-6))

    for iE in rgB_1F_NotRebuilt.Faces[0].AdjacentEdges():
        rc = rg.Curve.GetDistancesBetweenCurves(
            rgB_1F_NotRebuilt.Edges[iE],
            rgB_1F_Rebuilt.Edges[iE],
            tolerance=tolerance)
        if not rc[0]:
            s = "An EdgeCurve deviation could not be obtained."
            if bEcho: print(s)
            return
        else:
            deviations.append(rc[1])

    rgB_1F_NotRebuilt.Dispose()
    rgB_1F_Rebuilt.Dispose()

    return max(deviations)


def _facesWithEdgesOutOfTol(rgBrep, fEdgeDevs_perFace, fSearchTol=None):
    """
    """

    idx_Fs_WithOutOfTolEs = [] # Core faces to be rebuilt.
    for iF, fEdgeTol in enumerate(fEdgeDevs_perFace):
        if fEdgeTol is None or fEdgeTol > fSearchTol:
            idx_Fs_WithOutOfTolEs.append(iF)

    idx_Fs_Adj = []
    for iF in idx_Fs_WithOutOfTolEs:
        idx_Fs_Adj.extend(rgBrep.Faces[iF].AdjacentFaces())
    idx_Fs_Adj = sorted(set(idx_Fs_Adj))

    return sorted(set(idx_Fs_WithOutOfTolEs + idx_Fs_Adj))


def _faceEdgesControlPointCounts(rgFace, bDebug=None):
    """
    Returns: List of control point counts of each edge of face.
    """
    
    if bDebug is None: bDebug=Opts.values['bDebug']
    
    rgBrep = rgFace.Loops[0].Brep
    
    iCt_Cp_Face_PreRebuild = []
    for i, iE in enumerate(rgFace.AdjacentEdges()):
        rgEdge = rgBrep.Edges[iE]
        #dotCurveMidpoint(rgEdge, str(i))
        rgNurbsCrv = rgEdge.ToNurbsCurve()
        if rgNurbsCrv is None:
            if bDebug: sEval = 'rgNurbsCrv'; print(sEval+':', eval(sEval))
            continue
        numPts = rgNurbsCrv.Points.Count
        rgNurbsCrv.Dispose()
        iCt_Cp_Face_PreRebuild.append(numPts)
    return iCt_Cp_Face_PreRebuild


def _dotCurveMidpoint(rgCrv, text='', rgb=(255, 255, 255)):
    pt = rgCrv.PointAtNormalizedLength(0.5)
    if pt.X == Rhino.RhinoMath.UnsetValue:
        pt = rgCrv.PointAtStart
    rs.ObjectColor(rs.AddTextDot(text, pt), rgb)


def processBrep(rgBrep, fRebuildNakedTol=None, bRebuildOnlyOutOfTol=True, fSearchTol=None, bRebuildSharedEdges=False, fRebuildSharedTol=None, bRebuildVertices=None, bAddDot=False, iDotHeight=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep
        fRebuildNakedTol
        bRebuildOnlyOutOfTol
        bRebuildSharedEdges
        bRebuildVertices
        bAddDot
        iDotHeight
        bEcho
        bDebug

    Returns:
        rgBrep, None
        or
        None, str(Explanation of why there is no returned brep)
    """

    if not rgBrep.IsValid:
        return None, "Invalid brep skipped."

    rgB_In = rgBrep

    iCt_Cps_In = iCt_Cps_Out = 0

    iCt_Cps_In = iCt_LineCs_Out = iCt_ArcCs_In = iCt_PolyCs_In = iCt_PLCs_In = 0
    iCt_Cps_Out = iCt_LineCs_B1 = iCt_ArcCs_Out = iCt_PolyCs_Out = iCt_PLCs_Out = 0

    for e in rgB_In.Edges:
        c = e.EdgeCurve
        if isinstance(c, rg.NurbsCurve):
            iCt_Cps_In += c.Points.Count
        elif isinstance(c, rg.LineCurve):
            iCt_LineCs_Out += 1
        elif isinstance(c, rg.ArcCurve):
            iCt_ArcCs_In += 1
        elif isinstance(c, rg.PolyCurve):
            iCt_PolyCs_In += 1
        elif isinstance(c, rg.PolylineCurve):
            iCt_PLCs_In += 1
        else:
            1/0


    fEdgeTols_In = [rgE.Tolerance for rgE in rgB_In.Edges]
    if bDebug: sEval = "max(fEdgeTols_In)"; print(sEval,'=',eval(sEval))

    rgB_Out = rgB_In.DuplicateBrep()

    rgB_Out.Standardize() # See SDK for best documentation for Standardize.
    if bDebug:
        fEdgeTols_after_Standardize = [rgE.Tolerance for rgE in rgB_In.Edges]
        sEval = "max(fEdgeTols_after_Standardize)"; print(sEval,'=',eval(sEval))

    rgB_Out.Compact()
    if bDebug:
        fEdgeTols_after_Compact = [rgE.Tolerance for rgE in rgB_In.Edges]
        sEval = "max(fEdgeTols_after_Compact)"; print(sEval,'=',eval(sEval))

    setTolerancesBoxesAndFlags(rgB_Out)

    fEdgeTols_after_SetTolBoxesAndFlags = [rgE.Tolerance for rgE in rgB_Out.Edges]
    if bDebug: sEval = "max(fEdgeTols_after_SetTolBoxesAndFlags)"; print(sEval,'=',eval(sEval))

    maxDiscrepancyOfReportedTol = abs(
        max(fEdgeTols_In) -
        max(fEdgeTols_after_SetTolBoxesAndFlags))

    bReportedTolCorrected = maxDiscrepancyOfReportedTol > 0.1 * sc.doc.ModelAbsoluteTolerance

    if bEcho:
        s  = "Max. of all input brep's BrepEdge.Tolerance"
        s += " Before=>After Brep.SetTolerancesBoxesAndFlags"
        s += ": {}=>{}".format(
            _formatDistance(max(fEdgeTols_In)),
            _formatDistance(max(fEdgeTols_after_SetTolBoxesAndFlags)))
        print(s)

    #fZeroTol = 1e-6 * Rhino.RhinoMath.UnitScale(
    #    Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)

    fEdgeDevs_perFace = []
    for rgF in rgB_In.Faces:
        fDev = _compute_max_edge_tol_of_face(
            rgF,
            fRebuildNakedTol=fRebuildNakedTol,
            bRebuildVertices=bRebuildVertices)
        fEdgeDevs_perFace.append(fDev)

    if bEcho:
        print("Computed maximum edge deviation of input brep: {}.".format(
            _formatDistance(max(fEdgeDevs_perFace))))


    if fSearchTol is None:
        # Will rebuild all.
        idx_Fs_toRE = range(rgB_Out.Faces.Count)
    else:
        if all(fEdgeDevs_perFace) and max(fEdgeDevs_perFace) <= fSearchTol:
            idx_Fs_toRE = []
        else:
            idx_Fs_toRE = _facesWithEdgesOutOfTol(
                rgB_Out,
                fEdgeDevs_perFace,
                fSearchTol=fSearchTol,
                )

        if not idx_Fs_toRE:
            s = "All curve deviations are within tolerance"
            if bReportedTolCorrected:
                return rgB_Out, s + ", but Edge.Tolerance value(s) were corrected."
            else:
                return None, s + "."

        if bDebug:
            print("Up to {} faces will have their edges rebuilt.".format(
                len(idx_Fs_toRE)))


    idx_Fs_withREs = []

    for iF in idx_Fs_toRE:

        rgFace0 = rgB_In.Faces[iF]
        rgFace1 = rgB_Out.Faces[iF]
        
        if rgB_In.Faces.Count == 1:
            iCt_Cp_Face_PreRebuild = _faceEdgesControlPointCounts(rgFace0)
            if not iCt_Cp_Face_PreRebuild :
                sEval = 'iCt_Cp_Face_PreRebuild'
                return None, sEval+':', eval(sEval)
            
            if bEcho:
                print("Control Point Count of Each Edge")
                print("Original: {}".format(iCt_Cp_Face_PreRebuild))
        
        # Rebuild to input tolerance.
        
        # TODO: Determine this:
        # Should a seam be considered a naked edge?
        # Currently, it is not
        
        
        
        if (
            fSearchTol is not None and
            None not in fEdgeDevs_perFace and
            fEdgeDevs_perFace[iF] <= fSearchTol
        ):
            # Do not rebuild.
            continue # to next face.


        # Rebuild.
        bFaceHasSomeNakedEdges = any(rgB_Out.Edges[iE].Valence == rg.EdgeAdjacency.Naked for iE in rgFace1.AdjacentEdges())
        if bFaceHasSomeNakedEdges:
            bFaceHasAllNakedEdges = all(rgB_Out.Edges[iE].Valence == rg.EdgeAdjacency.Naked for iE in rgFace1.AdjacentEdges())
        else:
            bFaceHasAllNakedEdges = False

        if bFaceHasAllNakedEdges:
            bFaceHasSomeInteriorEdges = False
            bFaceHasAllInteriorEdges = False
        else:
            bFaceHasSomeInteriorEdges = any(rgB_Out.Edges[iE].Valence == rg.EdgeAdjacency.Interior for iE in rgFace1.AdjacentEdges())
            if bFaceHasSomeInteriorEdges:
                bFaceHasAllInteriorEdges = all(rgB_Out.Edges[iE].Valence == rg.EdgeAdjacency.Interior for iE in rgFace1.AdjacentEdges())
            else:
                bFaceHasAllInteriorEdges = False
        
        if bFaceHasAllNakedEdges:
            # Using fRebuildNakedTol.
            if rgFace1.RebuildEdges(
                tolerance=fRebuildNakedTol,
                rebuildSharedEdges=False,
                rebuildVertices=bRebuildVertices
            ):
                idx_Fs_withREs.append(iF)
        elif bFaceHasAllInteriorEdges:
            if bRebuildSharedEdges:
                if rgFace1.RebuildEdges(
                    tolerance=fRebuildSharedTol,
                    rebuildSharedEdges=True,
                    rebuildVertices=bRebuildVertices
                ):
                    idx_Fs_withREs.append(iF)
        else:
            # Face has a mix of naked and interior edges.
            if bRebuildSharedEdges:
                if rgFace1.RebuildEdges(
                    tolerance=fRebuildSharedTol,
                    rebuildSharedEdges=True,
                    rebuildVertices=bRebuildVertices
                ):
                    idx_Fs_withREs.append(iF)
            
            # If bRebuildSharedEdges, this will be the 2nd RebuildEdges
            # but with different arguments.
            if rgFace1.RebuildEdges(
                tolerance=fRebuildNakedTol,
                rebuildSharedEdges=False,
                rebuildVertices=bRebuildVertices
            ):
                if iF not in idx_Fs_withREs: idx_Fs_withREs.append(iF)

        #rgB_Out.Compact()

        # TODO: Standardize can modify the brep edges if any are too short.
        #  Investigate if and how Standardize can be used.  It is now disabled.
        #rgB_Out.Standardize() # See SDK for best documentation for Standardize.
        #rgB_Out.Compact()
        #sc.doc.Objects.AddBrep(rgB_Out); sc.doc.Views.Redraw(); 1/0

        if rgB_In.Faces.Count == 1:
            # Record new control point counts and dot old and new.
            iCt_Cp_Face_PostRebuild = _faceEdgesControlPointCounts(rgFace1)
            if iCt_Cp_Face_PostRebuild is None:
                sEval = 'iCt_Cp_Face_PostRebuild'
                return None, sEval+':', eval(sEval)
            
            if bAddDot:
                for i, idxEdge in enumerate(rgFace0.AdjacentEdges()):
                    if iCt_Cp_Face_PreRebuild[i] != iCt_Cp_Face_PostRebuild[i]:
                        s = "{}: {} -> {}".format(
                                idxEdge, iCt_Cp_Face_PreRebuild[i], iCt_Cp_Face_PostRebuild[i])
                        _dotCurveMidpoint(rgB_In.Edges[idxEdge], s)
            
            if bEcho:
                print("Rebuilt:  {}".format(iCt_Cp_Face_PostRebuild))

    if bDebug:
        print("{} faces had their edges rebuilt.".format(len(set(idx_Fs_withREs))))
        fEdgeTols_after_Rebuild = [rgE.Tolerance for rgE in rgB_Out.Edges]
        sEval = "max(fEdgeTols_after_Rebuild)"; print(sEval,'=',eval(sEval))

    rgB_Out.Compact()

    if bDebug:
        fEdgeTols_after_Compact = [rgE.Tolerance for rgE in rgB_Out.Edges]
        sEval = "max(fEdgeTols_after_Compact)"; print(sEval,'=',eval(sEval))

    setTolerancesBoxesAndFlags(rgB_Out)

    fEdgeTols_after_SetTolBoxesAndFlags = [rgE.Tolerance for rgE in rgB_Out.Edges]
    if bDebug: sEval = "max(fEdgeTols_after_SetTolBoxesAndFlags)"; print(sEval,'=',eval(sEval))


    for e in rgB_Out.Edges:
        c = e.EdgeCurve
        if isinstance(c, rg.NurbsCurve):
            iCt_Cps_Out += c.Points.Count
        elif isinstance(c, rg.LineCurve):
            iCt_LineCs_B1 += 1
        elif isinstance(c, rg.ArcCurve):
            iCt_ArcCs_Out += 1
        elif isinstance(c, rg.PolyCurve):
            iCt_PolyCs_Out += 1
        elif isinstance(c, rg.PolylineCurve):
            iCt_PLCs_Out += 1
        else:
            1/0

    ct_delta_Cp = iCt_Cps_Out - iCt_Cps_In

    bChangeInCrvTypeCts = (
            iCt_LineCs_B1 != iCt_LineCs_Out
            or
            iCt_ArcCs_Out != iCt_ArcCs_In
            or
            iCt_PolyCs_Out != iCt_PolyCs_In
            or 
            iCt_PLCs_Out != iCt_PLCs_In
    )

    bChangeInNCPtCts = iCt_Cps_Out != iCt_Cps_In

    if not (bChangeInCrvTypeCts or bChangeInNCPtCts):
        s  = "No change in curve type counts or NurbsCurve point counts."
        
        for iE in xrange(rgB_In.Edges.Count):
            c0 = rgB_In.Edges[iE].EdgeCurve
            c1 = rgB_Out.Edges[iE].EdgeCurve
            rc = rg.Curve.GetDistancesBetweenCurves(c0, c1, tolerance=0.1*fRebuildNakedTol)
            if not rc[0]:
                s += "  EdgeCurve deviation could not be obtained."
                s += "  Assuming edge has changed."
                if bEcho: print(s)
                break
            if rc[1] > fRebuildNakedTol:
                tA = rc[2]
                ptA_FromGDBC = c0.PointAt(tA)
                tB = rc[3]
                ptB_FromGDBC = c1.PointAt(tB)
                rc = c0.ClosestPoint(ptB_FromGDBC)
                if rc[0]:
                    ptA_ClosestToOnB = c0.PointAt(rc[1])
                rc = c1.ClosestPoint(ptA_FromGDBC)
                if rc[0]:
                    ptB_ClosestToOnA = c1.PointAt(rc[1])
                if (
                        ptA_FromGDBC.DistanceTo(ptB_ClosestToOnA) > fRebuildNakedTol
                        or
                        ptB_FromGDBC.DistanceTo(ptA_ClosestToOnB) > fRebuildNakedTol
                ):
                    if rc[1] < 0.001:
                        s += "  But EdgeCurve deviation of {:.2e} found.".format(rc[1])
                    else:
                        s += "  But EdgeCurve deviation of {:.4f} found.".format(rc[1])
                    if bEcho: print(s)
                    break
                if bDebug: print("Bad maxDistance from GetDistancesBetweenCurves.")
        else:
            if None not in fEdgeDevs_perFace:
                rgB_Out.Dispose()
                return None, s + " All curve deviations are within search tolerance."

    if bEcho:
        s  = "Total edge control point count change:"
        s +=  " {}{} ({} -> {}).".format(
            '+' if ct_delta_Cp > 0 else '',
            ct_delta_Cp,
            iCt_Cps_In,
            iCt_Cps_Out)
        print(s)

        fEdgeDevs_perFace = []
        for rgF in rgB_Out.Faces:
            fDev = _compute_max_edge_tol_of_face(
                rgF,
                fRebuildNakedTol=fRebuildNakedTol,
                bRebuildVertices=bRebuildVertices)
            fEdgeDevs_perFace.append(fDev)

        print("Computed maximum edge tolerance of resultant brep: {}.".format(
            _formatDistance(max(fEdgeDevs_perFace))))

    if bDebug:
        print(iCt_Cps_In, iCt_LineCs_Out, iCt_ArcCs_In)
        print(iCt_Cps_Out, iCt_LineCs_B1, iCt_ArcCs_Out)


    return rgB_Out, None


def _coerceBrepObject(rhObj):
    if isinstance(rhObj, rd.BrepObject):
        return rhObj
    if isinstance(rhObj, rd.ObjRef):
        #print(rhObj.GeometryComponentIndex.ComponentIndexType
        rdObj = rhObj.Object()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
    return rdObj if isinstance(rdObj, rd.BrepObject) else None


def processBrepObjects(rhBreps, **kwargs):
    """
    rhBreps: Can be objrefs, GUIDs, 
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fRebuildNakedTol = getOpt('fRebuildNakedTol')
    bRebuildOnlyOutOfTol = getOpt('bRebuildOnlyOutOfTol')
    fSearchTol = getOpt('fSearchTol') if getOpt('bRebuildOnlyOutOfTol') else None
    bRebuildSharedEdges = getOpt('bRebuildSharedEdges')
    fRebuildSharedTol = getOpt('fRebuildSharedTol')
    bRebuildVertices = getOpt('bRebuildVertices')
    bReplace = getOpt('bReplace')
    bAddDot = getOpt('bAddDot')
    iDotHeight = getOpt('iDotHeight')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rgBreps1 = []
    sLogs = []
    gBreps1 = []

    for rhBrep0 in rhBreps:
        rdBrep_In = _coerceBrepObject(rhBrep0)

        if rdBrep_In is None:
            if isinstance(rhBrep0, rg.Brep):
                rgB0 = rhBrep0
                rdBrep_In = None
            else:
                raise ValueError(
                    "{} passed to processBrepObjects." \
                        " Only ObjRef, DocObject, Geometry, or an iterable combination" \
                        " are accepted input.".format(rhBrep0.GetType()))

        rgB0 = rdBrep_In.BrepGeometry

        rc = processBrep(
            rgBrep=rgB0,
            fRebuildNakedTol=fRebuildNakedTol,
            bRebuildOnlyOutOfTol=bRebuildOnlyOutOfTol,
            fSearchTol=fSearchTol,
            bRebuildSharedEdges=bRebuildSharedEdges,
            fRebuildSharedTol=fRebuildSharedTol,
            bRebuildVertices=bRebuildVertices,
            bAddDot=bAddDot,
            iDotHeight=iDotHeight,
            bEcho=bEcho if len(rhBreps) == 1 else False,
            bDebug=bDebug
        )
        if not rc[0]:
            sLogs.append(rc[1])
        else:
            rgB_Out = rc[0]
            sLogs.append(rc[1])
            rgBreps1.append(rgB_Out)
            if bReplace:
                bSuccess = sc.doc.Objects.Replace(rdBrep_In.Id, rgB_Out)
                if bSuccess: gBreps1.append(rdBrep_In.Id)
            else:
                gBrep1 = sc.doc.Objects.AddBrep(rgB_Out)
                bSuccess = gBrep1 != Guid.Empty
                if bSuccess: gBreps1.append(gBrep1)
        rgB0.Dispose()
    
    for b in rgBreps1: b.Dispose()
    
    iBreps0_ct = len(rhBreps)

    if bEcho:
        if len(rhBreps) == 1:
            s = ""
            if sLogs and sLogs[0]:
                s += sLogs[0]
                s += "  "
            if bReplace:
                if gBreps1:
                    s += "Brep was modified.".format(
                            len(gBreps1), iBreps0_ct)
                else:
                    s +=  "Brep was not modified."
            else:
                if gBreps1:
                    s += "Brep has been added.".format(
                            len(gBreps1), iBreps0_ct)
                else:
                    s += "Brep was not added."
            print(s)
        elif len(rhBreps) > 1:
            for sFail in set(sLogs):
                print("[{}] {}".format(sLogs.count(sFail), sFail))
            if bReplace:
                if gBreps1:
                    print("{} of {} brep(s) have been modified.".format(
                            len(gBreps1), iBreps0_ct))
                else:
                    print("No breps have been modified.")
            else:
                if gBreps1:
                    print("{} of {} brep(s) have been added.".format(
                            len(gBreps1), iBreps0_ct))
                else:
                    print("No breps have been added.")
    
    return gBreps1


def main():
    
    rhBreps = getInput()
    if rhBreps is None: return

    fRebuildNakedTol = Opts.values['fRebuildNakedTol']
    bRebuildOnlyOutOfTol = Opts.values['bRebuildOnlyOutOfTol']
    fSearchTol = Opts.values['fSearchTol']
    bRebuildSharedEdges = Opts.values['bRebuildSharedEdges']
    fRebuildSharedTol = Opts.values['fRebuildSharedTol']
    bRebuildVertices = Opts.values['bRebuildVertices']
    bReplace = Opts.values['bReplace']
    bAddDot = Opts.values['bAddDot']
    iDotHeight = Opts.values['iDotHeight']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if bDebug:
        pass
    else:
        sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    rc = processBrepObjects(
        rhBreps,
        fRebuildNakedTol=fRebuildNakedTol,
        bRebuildOnlyOutOfTol=bRebuildOnlyOutOfTol,
        fSearchTol=fSearchTol,
        bRebuildSharedEdges=bRebuildSharedEdges,
        fRebuildSharedTol=fRebuildSharedTol,
        bRebuildVertices=bRebuildVertices,
        bReplace=bReplace,
        bAddDot=bAddDot,
        iDotHeight=iDotHeight,
        bEcho=bEcho,
        bDebug=bDebug
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()