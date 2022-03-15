"""
160704: Created starting with FindShortEdges.py.
...
191104: Actual edge count, not only over fOverlap_MinAllowed, is now compared with iCt_MaxEdgesOfFace.
191107: Bug fix.
200619, 0701: Import-related update.
210209: Again, only edges over a fOverlap_MinAllowed are counted.  Why was this changed on 191104?
210429: Modified an option default value.
210506: Fixed bug in main routine and added more options for quick tolerance value changes.
210519: Now will find faces that have all short non-seam edges.
210526: Trialing another quick select value.
210630: Added progress status of faces processed.
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


    key = 'fOverlap_MinAllowed'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'MinTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAnyEdgeCt'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCt_MaxEdgesOfFace'; keys.append(key)
    values[key] = 2
    names[key] = 'MaxEdgeCt'
    riOpts[key] = ri.Custom.OptionInteger(values[key], True, 1)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAllOverlapsMustBeSlivers'; keys.append(key)
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

        if not idxOpt: print "Add option for {} failed.".format(key)

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fOverlap_MinAllowed':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue == 0.0:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


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
    
    s  = "Tolerance = Minimum tolerance of the maximum edge overlap"
    s += "  0 of MaxEdgeCount means no limit"
    s += "  |  MaxEdgeCount does not include seams."
    print s

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fOverlap_MinAllowed')
        key = 'TenthMT'; idxs_Opt[key] = go.AddOption(key)
        key = 'HalfMT'; idxs_Opt[key] = go.AddOption(key)
        key = 'MT'; idxs_Opt[key] = go.AddOption(key)
        key = 'TwoMT'; idxs_Opt[key] = go.AddOption(key)
        key = 'TenMT'; idxs_Opt[key] = go.AddOption(key)
        addOption('bAnyEdgeCt')
        if not Opts.values['bAnyEdgeCt']:
            addOption('iCt_MaxEdgesOfFace')
        addOption('bAllOverlapsMustBeSlivers')
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
            return (
                objrefs,
                Opts.values['fOverlap_MinAllowed'],
                None if Opts.values['bAnyEdgeCt'] else Opts.values['iCt_MaxEdgesOfFace'],
                Opts.values['bAllOverlapsMustBeSlivers'],
                Opts.values['bExtract'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Nothing:
            iter = rd.ObjectEnumeratorSettings()
            iter.NormalObjects = True
            iter.LockedObjects = False
            iter.IncludeLights = False
            iter.IncludeGrips = False
            rdBrepObjects = []
            for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
                if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
                    rdBrepObjects.append(rdRhinoObject)
            go.Dispose()
            return (
                rdBrepObjects,
                Opts.values['fOverlap_MinAllowed'],
                None if Opts.values['bAnyEdgeCt'] else Opts.values['iCt_MaxEdgesOfFace'],
                Opts.values['bAllOverlapsMustBeSlivers'],
                Opts.values['bExtract'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fOverlap_MinAllowed'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['TenthMT']:
            key = 'fOverlap_MinAllowed'
            Opts.riOpts[key].CurrentValue = 0.1 * sc.doc.ModelAbsoluteTolerance
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['HalfMT']:
            key = 'fOverlap_MinAllowed'
            Opts.riOpts[key].CurrentValue = 0.5 * sc.doc.ModelAbsoluteTolerance
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['MT']:
            key = 'fOverlap_MinAllowed'
            Opts.riOpts[key].CurrentValue = sc.doc.ModelAbsoluteTolerance
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['TwoMT']:
            key = 'fOverlap_MinAllowed'
            Opts.riOpts[key].CurrentValue = 2.0 * sc.doc.ModelAbsoluteTolerance
            Opts.setValue(key)
            continue

        if go.Option().Index == idxs_Opt['TenMT']:
            key = 'fOverlap_MinAllowed'
            Opts.riOpts[key].CurrentValue = 10.0 * sc.doc.ModelAbsoluteTolerance
            Opts.setValue(key)
            continue


        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getFaces(rgBrep, fOverlap_MinAllowed=None, iCt_MaxEdgesOfFace=None, bAllOverlapsMustBeSlivers=True, bEcho=False, bDebug=False):
    """
    Search all faces of brep for slivers.

    Returns:
        sorted(idxs_rgFaces_Pass), nOverlapPairCt, fOverlap_Brep_MaxBelowTol
        Or None on error.
    """

    rgB = rgBrep

    if fOverlap_MinAllowed == None:
        fOverlap_MinAllowed = sc.doc.ModelAbsoluteTolerance


    # Previous sliver-finding routine.
    #    def isFaceSliver(rgFace, fSliverMaximumTolerance=2.0*sc.doc.ModelAbsoluteTolerance, bEcho=True):
    #        
    #        rgBrep_1F = rgFace.DuplicateFace(False)
    #        # Use DuplicateNakedEdgeCurves instead of AdjacentEdges to skip seams.
    #        rgCrvs = rgBrep_1F.DuplicateNakedEdgeCurves(outer=True, inner=True)
    #        rgBrep_1F.Dispose()
    #        
    #        # Reject faces with number of edges exceeding 2.
    #        if len(rgCrvs) != 2:
    #            map(lambda x: x.Dispose(), rgCrvs)
    #            return False
    #        
    #        rc = indexPairsOfOverlappingCurves(rgCrvs, fSliverMaximumTolerance, bAllOverlapsMustBeSlivers=True)
    #        map(lambda x: x.Dispose(), rgCrvs)
    #        
    #        if rc is None or len(rc[0]) == 0: return False
    #        
    #        return True


    def indexPairsOfOverlappingCurves(rgCrvs, fOverlap_MinAllowed, bAllOverlapsMustBeSlivers):

        idx_rgCrvs_OverlapPairs = []
        #fOverlap_Face_Min = None
        fOverlap_Face_MaxBelowTol = None

        # 180702: Added.
        for rgC in rgCrvs:
            if not rgC.IsValid:
                print "Warning: Curve is invalid, so this group of curves will be skipped because rg.Curve.GetDistancesBetweenCurves may hang."
                return

        if len(rgCrvs) == 1:
            rgC_In = rgCrvs[0]
            rgCs_CheckLength = [
                rgC_In.Trim(t0=rgC_In.Domain.T0, t1=rgC_In.Domain.Mid),
                rgC_In.Trim(t0=rgC_In.Domain.Mid, t1=rgC_In.Domain.T1)
            ]
            bSingle_crv_was_split = True
        elif len(rgCrvs) >= 2:
            rgCs_CheckLength = rgCrvs
            bSingle_crv_was_split = False
        else:
            raise ValueError("Curve set has no curves!")


        rgCs_LongEnough = []
        for rgC in rgCrvs:
            length = rgC.GetLength()
            if length >= fOverlap_MinAllowed:
                rgCs_LongEnough.append(rgC)


        iCt_Cs = len(rgCs_LongEnough)

        for iE_A in xrange(iCt_Cs):
            rgC_A = rgCs_LongEnough[iE_A]
            for iE_B in xrange(iE_A+1, iCt_Cs):
                rgC_B = rgCs_LongEnough[iE_B]
                rc = rg.Curve.GetDistancesBetweenCurves(
                    rgC_A,
                    rgC_B,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if rc[0]:
                    fOverlap_Max = rc[1]
                    fOverlap_Min = rc[4]
                    if fOverlap_Max < fOverlap_MinAllowed:
                        if bSingle_crv_was_split:
                            return (0, 0), fOverlap_Max

                        idx_rgCrvs_OverlapPairs.append((iE_A, iE_B))
        #                    if fOverlap_Face_Min is None or fOverlap_Min < fOverlap_Face_Min:
        #                        fOverlap_Face_Min = fOverlap_Min
                        if fOverlap_Face_MaxBelowTol is None or fOverlap_Max > fOverlap_Face_MaxBelowTol:
                            fOverlap_Face_MaxBelowTol = fOverlap_Max
                    elif bAllOverlapsMustBeSlivers:
                        # Skip any face with an overlap larger than fOverlap_MinAllowed.
                        return

        return idx_rgCrvs_OverlapPairs, fOverlap_Face_MaxBelowTol


    idxs_rgFaces_Pass = []
    #fOverlap_Min_Brep = None
    fOverlaps_Brep_MaxBelowTol = []
    nOverlapPairCt = 0


    iCt_Fs = rgBrep.Faces.Count

    idxs_AtTenths = [int(round(0.1*i*iCt_Fs,0)) for i in range(10)]

    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt


    for idxF in xrange(rgB.Faces.Count):
        
        # 180702: Commented out.
    #        rgBrep_1F = rgFace.DuplicateFace(False)
    #        # Use DuplicateNakedEdgeCurves instead of AdjacentEdges to skip seams.
    #        rgCrvs = rgBrep_1F.DuplicateNakedEdgeCurves(outer=True, inner=True)
    #        rgBrep_1F.Dispose()

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

        rgCrvs = []
        edges_face = rgF.AdjacentEdges()

        for idx_rgEdge in edges_face:
            rgEdge = rgB.Edges[idx_rgEdge]

            # Skip short edges.
            if len(edges_face) > 2:
                fLengthEdge = rgEdge.GetLength()
                if fLengthEdge < fOverlap_MinAllowed:
                    continue

            # Skip seams.
            if rgEdge.Valence == rg.EdgeAdjacency.Interior:
                rgFaces_Adj = rgEdge.AdjacentFaces()
                if rgFaces_Adj[0] == rgFaces_Adj[1]:
                    continue

            rgCrvs.append(rgEdge)

        # Skip faces with number of edges exceeding iCt_MaxEdgesOfFace.
        if iCt_MaxEdgesOfFace and len(rgCrvs) > iCt_MaxEdgesOfFace:
            continue
        
        if len(rgCrvs) == 0:
            # All short non-seam edges?
            idxs_rgFaces_Pass.append(rgF.FaceIndex)
            continue

        rc = indexPairsOfOverlappingCurves(rgCrvs, fOverlap_MinAllowed, bAllOverlapsMustBeSlivers)
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
        fOverlap_MinAllowed,
        iCt_MaxEdgesOfFace,
        bAllOverlapsMustBeSlivers,
        bExtract,
        bEcho,
        bDebug,
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fOverlap_MinAllowed = getOpt('fOverlap_MinAllowed')
    iCt_MaxEdgesOfFace = getOpt('iCt_MaxEdgesOfFace')
    bAllOverlapsMustBeSlivers = getOpt('bAllOverlapsMustBeSlivers')
    bExtract = getOpt('bExtract')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    fOverlap_Min_All = fOverlap_Max_All = None
    fOverlaps_BelowTol_AllBreps = []
    gBreps_Extracted_All = [] # Accumulation of duplicated faces (breps)
    iOverlap_ct_total = 0
    
    len_rhBs0 = len(rhBreps0)
    idxs_AtTenths = [int(round(0.1*i*len(rhBreps0),0)) for i in range(10)]


    # escape_test is sometimes called when this script is called by a (Rhino Python Editor: New -) Command.
    if sc.escape_test(throw_exception=False, reset=True):
        print "sc.escape_test was triggered."


    iCt_FacesSel_All = 0

    for iB, rhBrep0 in enumerate(rhBreps0):
        if sc.escape_test(False):
            print "Searching interrupted by user."
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
            print "GUID could not be obtained."
            continue
        
        rgBrep0 = rdBrep0.Geometry
        if rgBrep0 is None or (rgBrep0.IsSolid and rgBrep0.IsSurface): continue
        if not rgBrep0.IsValid and bEcho:
            print "Brep {} is invalid and will be skipped because rg.Curve.GetDistancesBetweenCurves may hang.".format(gBrep0)
            continue
        
        rc = getFaces(
                rgBrep=rgBrep0,
                fOverlap_MinAllowed=fOverlap_MinAllowed,
                iCt_MaxEdgesOfFace=iCt_MaxEdgesOfFace,
                bAllOverlapsMustBeSlivers=bAllOverlapsMustBeSlivers,
                bEcho=bEcho,
                bDebug=bDebug)
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
    
    if iOverlap_ct_total > 0:
        print "{} completely overlapping edge(s) found.".format(iOverlap_ct_total)
        
        if bExtract:
            if len(gBreps_Extracted_All) > 1:
                sc.doc.Objects.Select(List[Guid](gBreps_Extracted_All))
            else:
                sc.doc.Objects.Select(gBreps_Extracted_All[0])
            rdBreps1_Selected = [rdObj for rdObj in
                    sc.doc.Objects.GetSelectedObjects(False, False)]
            print "{} mono-face brep(s) selected.".format(len(rdBreps1_Selected))
        elif iCt_FacesSel_All:
            print "{} face(s) selected.".format(iCt_FacesSel_All)

        print "Range of distances of completely overlapping edges:"
        print " [{0:.{2}f}, {1:.{2}f}]".format(
                min(fOverlaps_BelowTol_AllBreps),
                max(fOverlaps_BelowTol_AllBreps),
                sc.doc.ModelDistanceDisplayPrecision+1)
    else:
        print "No completely overlapping edges found."


def main():
    
    rc = getInput()
    if rc is None: return
    (
        objrefs,
        fOverlap_MinAllowed,
        iCt_MaxEdgesOfFace,
        bAllOverlapsMustBeSlivers,
        bExtract,
        bEcho,
        bDebug,
       ) = rc
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    sc.doc.Objects.UnselectAll()

    processBrepObjects(
        rhBreps0=objrefs,
        fOverlap_MinAllowed=fOverlap_MinAllowed,
        iCt_MaxEdgesOfFace=iCt_MaxEdgesOfFace,
        bAllOverlapsMustBeSlivers=bAllOverlapsMustBeSlivers,
        bExtract=bExtract,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()