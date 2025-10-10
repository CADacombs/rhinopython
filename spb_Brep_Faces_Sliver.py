"""
This script finds faces whose edge-edge curve deviation is <= an input value.

For faces with short edges, their breps often require more repair, e.g., merging edges in adjacent faces.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160704: Created starting with FindShortEdges.py.
...
210519: Now will find faces that have all short non-seam edges.
210526: Trialing another quick select value.
210630: Added progress status of faces processed.
220420: Repaired for Rhino 7+ due to a RhinoCommon 7.17 script-breaking change.
220817: Modified an option default value.
220818: Modified options and main routine.
220915: Replaced edge count option with another option.
221214: Bug fix in options.
230405: Modified an option default value.
230826, 0915: Modified an option default value.
240402: Added an option to skip short edges in deviation checks.
240526, 250324: Modified an option default value.
250325: Added the option where the MaxSliverWidth value is also used for the MaxShortEdgeLength.
250514: Disabled check of single-edge face since the GetDistanceBetweenCurves doesn't always report correctly.
        Modified available command options per settings of other options.
250916,24: Modified some option default values.
251008-09: Added an option to define the minimum length of an edge to include for faces to skip. Modified some option default values.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import xBrepObject


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fMaxSliverWidth'; keys.append(key)
    values[key] = 2.0 * sc.doc.ModelAbsoluteTolerance # 1.8 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bSkipFacesWithShortEdges'; keys.append(key)
    values[key] = False
    names[key] = 'FacesWithShortEdges'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Include', 'Skip')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSkipSliverCheckOfShortEdges'; keys.append(key)
    values[key] = True
    names[key] = 'CheckDevsOfShortEdges'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Include', 'Skip')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseSliverTolForMaxShortEdgeLengthTol'; keys.append(key)
    values[key] = False
    names[key] = 'UseSliverTolForMaxShortLength'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fMaxEdgeLengthConsideredShort'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fIgnoreEdgesBelowThisLength'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bEntireFaceMustBeASliver'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fMaxSliverWidth':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fMaxEdgeLengthConsideredShort':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fIgnoreEdgesBelowThisLength':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getAllNormalBreps():
    oes = rd.ObjectEnumeratorSettings()
    oes.NormalObjects = True
    oes.LockedObjects = False # Default is True.
    oes.IncludeLights = False
    oes.IncludeGrips = False
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep

    go.AcceptNothing(True)

    go.AcceptNumber(True, acceptZero=True)

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    bPreselectedObjsChecked = False

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('fMaxSliverWidth')
        key = 'Tenth'; idxs_Opt[key] = go.AddOption(key)
        key = 'Half'; idxs_Opt[key] = go.AddOption(key)
        key = 'MT'; idxs_Opt[key] = go.AddOption(key)
        key = 'Double'; idxs_Opt[key] = go.AddOption(key)
        key = 'TenX'; idxs_Opt[key] = go.AddOption(key)
        addOption('bSkipFacesWithShortEdges')
        if not Opts.values['bSkipFacesWithShortEdges']:
            addOption('bSkipSliverCheckOfShortEdges')
        if Opts.values['bSkipFacesWithShortEdges'] or Opts.values['bSkipSliverCheckOfShortEdges']:
            addOption('bUseSliverTolForMaxShortEdgeLengthTol')
            if not Opts.values['bUseSliverTolForMaxShortEdgeLengthTol']:
                addOption('fMaxEdgeLengthConsideredShort')
            addOption('fIgnoreEdgesBelowThisLength')
        addOption('bEntireFaceMustBeASliver')
        addOption('bExtract')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return getAllNormalBreps()

        if res == ri.GetResult.Number:
            key = 'fMaxSliverWidth'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['Tenth']:
            key = 'fMaxSliverWidth'
            Opts.riOpts[key].CurrentValue = 0.1 * Opts.riOpts['fMaxSliverWidth'].CurrentValue
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['Half']:
            key = 'fMaxSliverWidth'
            Opts.riOpts[key].CurrentValue = 0.5 * Opts.riOpts['fMaxSliverWidth'].CurrentValue
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['MT']:
            key = 'fMaxSliverWidth'
            Opts.riOpts[key].CurrentValue = sc.doc.ModelAbsoluteTolerance
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['Double']:
            key = 'fMaxSliverWidth'
            Opts.riOpts[key].CurrentValue = 2.0 * Opts.riOpts['fMaxSliverWidth'].CurrentValue
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['TenX']:
            key = 'fMaxSliverWidth'
            Opts.riOpts[key].CurrentValue = 10.0 * Opts.riOpts['fMaxSliverWidth'].CurrentValue
            Opts.setValue(key)
            continue


        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _indexPairsOfOverlappingCurves(rgCrvs, fMaxSliverWidth, bEntireFaceMustBeASliver):

    idx_rgCrvs_OverlapPairs = []
    #fOverlap_Face_Min = None
    fOverlap_Face_MaxBelowTol = None

    if not rgCrvs:
        raise ValueError("Curve set has no curves!")

    # 180702: Added.
    for rgC in rgCrvs:
        if not rgC.IsValid:
            print("Warning: Curve is invalid, so this group of curves will be skipped because rg.Curve.GetDistancesBetweenCurves may hang.")
            return

    if len(rgCrvs) == 1:
        # TODO: Determine a better routine for single edge faces. GetDistancesBetweenCurves doesn't always work well.
        return
        
        rgC_In = rgCrvs[0]
        
        tMidLength = rgC_In.DivideByCount(
            segmentCount=2,
            includeEnds=False)[0]
        
        rgC_A = rgC_In.Trim(t0=rgC_In.Domain.T0, t1=tMidLength)
        rgC_B = rgC_In.Trim(t0=tMidLength, t1=rgC_In.Domain.T1)
        # tMidLength was rgC_In.Domain.Mid
        rvs = rg.Curve.GetDistancesBetweenCurves(
            rgC_A,
            rgC_B,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
        if not rvs[0]:
            return
        
        fOverlap_Max = rvs[1]
        sc.doc.Objects.AddCurve(rgC_A)
        sc.doc.Objects.AddCurve(rgC_B)
        return
        sc.doc.Objects.AddPoint(rgC_In.PointAt(tMidLength))
        sc.doc.Objects.AddLine(rgC_A.PointAt(rvs[2]), rgC_B.PointAt(rvs[3])); sc.doc.Views.Redraw()
        fOverlap_Min = rvs[4]
        rgC_A.Dispose()
        rgC_B.Dispose()
        if fOverlap_Max > fMaxSliverWidth:
            return
        
        return (0, 0), fOverlap_Max

    iCt_Cs = len(rgCrvs)

    fOverlap_Maxs = []
    fOverlap_Maxs_Sliver = []

    for iE_A in xrange(iCt_Cs):
        rgC_A = rgCrvs[iE_A]
        for iE_B in xrange(iE_A+1, iCt_Cs):
            rgC_B = rgCrvs[iE_B]
            rc = rg.Curve.GetDistancesBetweenCurves(
                rgC_A,
                rgC_B,
                tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
            if not rc[0]:
                return

            fOverlap_Max = rc[1]
            fOverlap_Maxs.append(fOverlap_Max)
            fOverlap_Min = rc[4]
            if fOverlap_Max <= fMaxSliverWidth:
                idx_rgCrvs_OverlapPairs.append((iE_A, iE_B))
#                    if fOverlap_Face_Min is None or fOverlap_Min < fOverlap_Face_Min:
#                        fOverlap_Face_Min = fOverlap_Min
                fOverlap_Maxs_Sliver.append(fOverlap_Max)
                #if fOverlap_Face_MaxBelowTol is None or fOverlap_Max > fOverlap_Face_MaxBelowTol:
                #    fOverlap_Face_MaxBelowTol = fOverlap_Max
            elif bEntireFaceMustBeASliver:
                # Skip any face with an overlap larger than fMaxSliverWidth.
                return

    return idx_rgCrvs_OverlapPairs, max(fOverlap_Maxs_Sliver)


def getFaces(rgBrep, fMaxSliverWidth, bSkipFacesWithShortEdges, bSkipSliverCheckOfShortEdges, fMaxEdgeLengthConsideredShort, fIgnoreEdgesBelowThisLength, bEntireFaceMustBeASliver, bEcho=True, bDebug=False):
    """
    Search all faces of brep for slivers.

    Parameters:
        rgBrep
        fMaxSliverWidth
        bSkipFacesWithShortEdges
        bSkipSliverCheckOfShortEdges
        fMaxEdgeLengthConsideredShort
        fIgnoreEdgesBelowThisLength
        bEntireFaceMustBeASliver
        bEcho
        bDebug

    Returns:
        sorted(idxs_rgFaces_Pass), nOverlapPairCt, fOverlap_Brep_MaxBelowTol
        Or None on error.
    """


    rgB = rgBrep

    if fMaxSliverWidth == None:
        fMaxSliverWidth = sc.doc.ModelAbsoluteTolerance


    # Previous sliver-finding routine.
    #    def isFaceSliver(rgFace, fSliverMaximumTolerance=2.0*sc.doc.ModelAbsoluteTolerance, bEcho=True):
    #        
    #        rgBrep_1F = rgFace.DuplicateFace(False)
    #        # Use DuplicateNakedEdgeCurves instead of AdjacentEdges to skip seams.
    #        rgCrvs = rgBrep_1F.DuplicateNakedEdgeCurves(True, True)
    #        rgBrep_1F.Dispose()
    #        
    #        # Reject faces with number of edges exceeding 2.
    #        if len(rgCrvs) != 2:
    #            map(lambda x: x.Dispose(), rgCrvs)
    #            return False
    #        
    #        rc = _indexPairsOfOverlappingCurves(rgCrvs, fSliverMaximumTolerance, bEntireFaceMustBeASliver=True)
    #        map(lambda x: x.Dispose(), rgCrvs)
    #        
    #        if rc is None or len(rc[0]) == 0: return False
    #        
    #        return True


    idxs_rgFaces_Pass = []
    #fOverlap_Min_Brep = None
    fOverlaps_Brep_MaxBelowTol = []
    nOverlapPairCt = 0


    iCt_Fs = rgBrep.Faces.Count

    idxs_AtTenths = [int(round(0.1*i*iCt_Fs,0)) for i in range(10)]

    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt


    def getCrvsToCheck(rgF, bSkipFacesWithShortEdges, bSkipSliverCheckOfShortEdges, fMaxEdgeLengthConsideredShort, fIgnoreEdgesBelowThisLength):

        #if not bSkipFacesWithShortEdges and not bSkipSliverCheckOfShortEdges:
        #    return [rgF.Brep.Edges[i] for i in rgF.AdjacentEdges()]

        rgEs_Out = []

        for iE in rgF.AdjacentEdges():
            rgE = rgF.Brep.Edges[iE]

            # Skip seams.
            if rgE.Valence == rg.EdgeAdjacency.Interior:
                rgFaces_Adj = rgE.AdjacentFaces()
                if rgFaces_Adj[0] == rgFaces_Adj[1]:
                    continue

            fLength = rgE.GetLength()

            if fIgnoreEdgesBelowThisLength and (fLength < fIgnoreEdgesBelowThisLength):
                continue
            elif fMaxEdgeLengthConsideredShort and (fLength <= fMaxEdgeLengthConsideredShort):
                if bSkipFacesWithShortEdges:
                    return
                if bSkipSliverCheckOfShortEdges:
                    continue

            rgEs_Out.append(rgE)

        return rgEs_Out


    for idxF in xrange(rgB.Faces.Count):
        if iCt_Fs == 1:
            pass
        elif iCt_Fs < 100:
            Rhino.RhinoApp.SetCommandPrompt(
                sCmdPrompt0 +
                "  Analyzing {} of {} faces".format(
                    idxF+1, iCt_Fs))
        elif idxF in idxs_AtTenths:
            Rhino.RhinoApp.SetCommandPrompt(
                sCmdPrompt0 +
                "  Analysis at {:d}% of {} faces ...".format(
                    int(100.0 * (idxF+1) / iCt_Fs), iCt_Fs))


        rgF = rgB.Faces[idxF]

        # Allowing single edge faces because they will be split and checked.


        rgCrvs = getCrvsToCheck(rgF, bSkipFacesWithShortEdges, bSkipSliverCheckOfShortEdges, fMaxEdgeLengthConsideredShort, fIgnoreEdgesBelowThisLength)

        if not rgCrvs:
            continue

        if len(rgCrvs) == 0:
            # All short non-seam edges?
            idxs_rgFaces_Pass.append(rgF.FaceIndex)
            continue

        rc = _indexPairsOfOverlappingCurves(rgCrvs, fMaxSliverWidth, bEntireFaceMustBeASliver)
        if rc is None or len(rc[0]) == 0:
            continue

        idx_rgCrvs_OverlapPairs, fOverlap_Face_MaxBelowTol = rc
        
        idxs_rgFaces_Pass.append(rgF.FaceIndex)
        
        nOverlapPairCt += len(idx_rgCrvs_OverlapPairs)
        
        fOverlaps_Brep_MaxBelowTol.append(fOverlap_Face_MaxBelowTol)
        
        #        if fOverlap_Min_Brep is None or fOverlap_Face_Min < fOverlap_Min_Brep:
        #            fOverlap_Min_Brep = fOverlap_Face_Min
        #        if fOverlap_Brep_MaxBelowTol is None or fOverlap_Face_MaxBelowTol > fOverlap_Brep_MaxBelowTol:
        #            fOverlap_Brep_MaxBelowTol = fOverlap_Face_MaxBelowTol

    Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

    return sorted(idxs_rgFaces_Pass), nOverlapPairCt, fOverlaps_Brep_MaxBelowTol


def processBrepObjects(rhBreps0, **kwargs):
    """
    Parameters:
        rhBreps0
        fMaxSliverWidth,
        bSkipFacesWithShortEdges,
        bSkipSliverCheckOfShortEdges,
        fMaxEdgeLengthConsideredShort,
        fIgnoreEdgesBelowThisLength,
        bEntireFaceMustBeASliver,
        bExtract,
        bEcho,
        bDebug,
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fMaxSliverWidth = getOpt('fMaxSliverWidth')
    bSkipFacesWithShortEdges = getOpt('bSkipFacesWithShortEdges')
    bSkipSliverCheckOfShortEdges = getOpt('bSkipSliverCheckOfShortEdges')
    bUseSliverTolForMaxShortEdgeLengthTol = getOpt('bUseSliverTolForMaxShortEdgeLengthTol')
    fMaxEdgeLengthConsideredShort = getOpt('fMaxEdgeLengthConsideredShort')
    fIgnoreEdgesBelowThisLength = getOpt('fIgnoreEdgesBelowThisLength')
    bEntireFaceMustBeASliver = getOpt('bEntireFaceMustBeASliver')
    bExtract = getOpt('bExtract')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if bUseSliverTolForMaxShortEdgeLengthTol:
        fMaxEdgeLengthConsideredShort = fMaxSliverWidth

    fOverlap_Min_All = fOverlap_Max_All = None
    fOverlaps_BelowTol_AllBreps = []
    gBreps_Extracted_All = [] # Accumulation of duplicated faces (breps)
    iOverlap_ct_total = 0
    
    len_rhBs0 = len(rhBreps0)
    idxs_AtTenths = [int(round(0.1*i*len(rhBreps0),0)) for i in range(10)]


    # escape_test is sometimes called when this script is called by a (Rhino Python Editor: New -) Command.
    if sc.escape_test(throw_exception=False, reset=True):
        print("sc.escape_test was triggered.")


    iCt_FacesSel_All = 0

    for iB, rhBrep0 in enumerate(rhBreps0):
        if sc.escape_test(False):
            print("Searching interrupted by user.")
            return
        
        if len_rhBs0 == 1:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing brep ...")
        elif len_rhBs0 <= 12:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing {} of {} breps".format(
                iB+1, len_rhBs0))
        elif iB in idxs_AtTenths:
            # May be many monoface breps.
            Rhino.RhinoApp.SetCommandPrompt("Analysis at {:d}% of {} breps ...".format(
                int(100.0 * (iB+1) / len_rhBs0), len_rhBs0))

        # Obtain GUID, RhinoObject, and geometry.
        if isinstance(rhBrep0, Guid):
            gBrep0 = rhBrep0
            rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
            if rdBrep0 is None: continue
        elif isinstance(rhBrep0, rd.ObjRef):
            gBrep0 = rhBrep0.ObjectId
            rdBrep0 = rhBrep0.Object()
        elif isinstance(rhBrep0, rd.RhinoObject):
            rdBrep0 = rhBrep0
            gBrep0 = rhBrep0.Id
        else:
            print("GUID could not be obtained.")
            continue
        
        rgBrep0 = rdBrep0.Geometry
        if rgBrep0 is None or (rgBrep0.IsSolid and rgBrep0.IsSurface): continue
        if not rgBrep0.IsValid and bEcho:
            print("Brep {} is invalid and will be skipped because Curve.GetDistancesBetweenCurves may hang.".format(gBrep0))
            continue
        
        rc = getFaces(
            rgBrep=rgBrep0,
            fMaxSliverWidth=fMaxSliverWidth,
            bSkipFacesWithShortEdges=bSkipFacesWithShortEdges,
            bSkipSliverCheckOfShortEdges=bSkipSliverCheckOfShortEdges,
            fMaxEdgeLengthConsideredShort=fMaxEdgeLengthConsideredShort,
            fIgnoreEdgesBelowThisLength=fIgnoreEdgesBelowThisLength,
            bEntireFaceMustBeASliver=bEntireFaceMustBeASliver,
            bEcho=bEcho,
            bDebug=bDebug,
            )
        rgBrep0.Dispose()
        if rc is None: continue

        (
            idxs_rgFaces_Pass,
            iCt_Overlaps_1Brep,
            fOverlaps_MaxBelowTol_1Brep,
            ) = rc
    
        if not idxs_rgFaces_Pass: continue

        if bExtract:
            rc = xBrepObject.extractFaces(
                    gBrep0,
                    idxs_rgFaces_Pass,
                    bCurrentLayer=False,
                    bByLayerColor=False,
                    bAddOnlyMonofaces=True,
                    bEcho=False,
                    bDebug=bDebug)
            if rc is None: continue
            gBreps_Extracted_All.extend(rc[0])
        else:
            iCt_FacesSel = xBrepObject.selectFaces(
                    gBrep0,
                    idxs_rgFaces_Pass,
                    bEcho=bEcho)
            iCt_FacesSel_All += iCt_FacesSel
        
        iOverlap_ct_total += iCt_Overlaps_1Brep
        
        fOverlaps_BelowTol_AllBreps.extend(fOverlaps_MaxBelowTol_1Brep)
        
        #        if fOverlap_Min_All is None or fOverlap_Min_ThisBrep < fOverlap_Min_All:
        #            fOverlap_Min_All = fOverlap_Min_ThisBrep
        #        if fOverlap_Max_All is None or fOverlaps_MaxBelowTol_1Brep > fOverlap_Max_All:
        #            fOverlap_Max_All = fOverlaps_MaxBelowTol_1Brep


    if iOverlap_ct_total == 0:
        print("No completely overlapping edges found.")
        return


    print("{} completely overlapping edge(s) found.".format(iOverlap_ct_total))

    if bExtract:
        if len(gBreps_Extracted_All) > 1:
            sc.doc.Objects.Select(List[Guid](gBreps_Extracted_All))
        else:
            sc.doc.Objects.Select(gBreps_Extracted_All[0])
        rdBreps1_Selected = [rdObj for rdObj in
                sc.doc.Objects.GetSelectedObjects(False, False)]
        print("{} monoface brep(s) selected.".format(len(rdBreps1_Selected)))
    elif iCt_FacesSel_All:
        print("{} face(s) selected.".format(iCt_FacesSel_All))

    if len(fOverlaps_BelowTol_AllBreps) == 1:
        print("Width of sliver: {0:.{1}f}".format(
            min(fOverlaps_BelowTol_AllBreps),
            sc.doc.ModelDistanceDisplayPrecision+1))
    else:
        print("Range of widths of slivers:")
        print(" [{0:.{2}f}, {1:.{2}f}]".format(
            min(fOverlaps_BelowTol_AllBreps),
            max(fOverlaps_BelowTol_AllBreps),
            sc.doc.ModelDistanceDisplayPrecision+1))


def main():

    objrefs = getInput()
    if objrefs is None: return

    fMaxSliverWidth = Opts.values['fMaxSliverWidth']
    bSkipFacesWithShortEdges = Opts.values['bSkipFacesWithShortEdges']
    bSkipSliverCheckOfShortEdges = Opts.values['bSkipSliverCheckOfShortEdges']
    bUseSliverTolForMaxShortEdgeLengthTol = Opts.values['bUseSliverTolForMaxShortEdgeLengthTol']
    fMaxEdgeLengthConsideredShort = Opts.values['fMaxEdgeLengthConsideredShort'] if Opts.values['bSkipFacesWithShortEdges'] else 0.0
    fIgnoreEdgesBelowThisLength = Opts.values['fIgnoreEdgesBelowThisLength']
    bEntireFaceMustBeASliver = Opts.values['bEntireFaceMustBeASliver']
    bExtract = Opts.values['bExtract']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    processBrepObjects(
        objrefs,
        fMaxSliverWidth=fMaxSliverWidth,
        bSkipFacesWithShortEdges=bSkipFacesWithShortEdges,
        bSkipSliverCheckOfShortEdges=bSkipSliverCheckOfShortEdges,
        fMaxEdgeLengthConsideredShort=fMaxSliverWidth if bUseSliverTolForMaxShortEdgeLengthTol else fMaxEdgeLengthConsideredShort,
        fIgnoreEdgesBelowThisLength=fIgnoreEdgesBelowThisLength,
        bEntireFaceMustBeASliver=bEntireFaceMustBeASliver,
        bExtract=bExtract,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()