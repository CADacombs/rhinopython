"""
This script uses Brep.CutUpSurface and optionally BrepFace.RebuildEdges
(before 1st CutUpSurface) to attempt to retrim all faces of a brep individually
to simplify edge curves and BrepTrim curves.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250327-28,30: Created.
250603: Bug fixes.
251001: Bug fixes. Refactored. Added a command option. Added more printed output, debug and not debug.

useEdgeCurves
From https://developer.rhino3d.com/api/rhinocommon/rhino.geometry.brep/cutupsurface :
The 2D trimming curves are made by pulling back the 3D curves using the fitting tolerance.
If useEdgeCurves is true, the input 3D curves will be used as the edge curves in the result.
Otherwise, the edges will come from pushing up the 2D pullbacks.

Restated:
useEdgeCurves is True: Use input curves as the edge curves
useEdgeCurves is False: Generate new edge curves from the newly created BrepTrim curves.

TODO:
    Try tighter tolerances on trim fails.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fTol_Edge_and_trim'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'EdgeAndTrimTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    key = 'fTol_Join'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'JoinTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    key = 'bRebuildEdges'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bStartingCrvs_Edge_not_UV'; keys.append(key)
    values[key] = True
    names[key] = 'StartingCrvs'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'UV', 'Edge')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryShrinkOnFail'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryExtendOnFail'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryOtherTolsOnFail'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    #key = 'bUntrimFailedFaces'; keys.append(key)
    #values[key] = True
    #names[key] = 'FailedTrimFaces'
    #riOpts[key] = ri.Custom.OptionToggle(values[key], 'InputFace', 'Untrim')
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'OutputBrep'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
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
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fTol_Edge_and_trim':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList
            return

        print("Invalid key?")


def getInput(bDebug=False):
    """
    Get breps.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select breps")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Brep

    #go.AcceptNumber(True, acceptZero=False)
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == ri.GetResult.Object 
    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        addOption('fTol_Edge_and_trim')
        addOption('fTol_Join')
        addOption('bRebuildEdges')
        if Opts.values['bRebuildEdges']:
            addOption('bStartingCrvs_Edge_not_UV')
        addOption('bTryShrinkOnFail')
        addOption('bTryExtendOnFail')
        addOption('bTryOtherTolsOnFail')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fTol_Edge_and_trim'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createBrepOfUntrimmedFace(rgFace_In):
    rgB_Out = rgFace_In.DuplicateSurface().ToBrep()
    rgB_Out.Faces[0].PerFaceColor = rgFace_In.PerFaceColor
    return rgB_Out


def _areaComp(rgBrep, fAreaRef, fDeltaTol, bDebug=False):
    """
    fDeltaTol: 0.01 means 1% tolerance
    """
    fArea_New = rgBrep.GetArea() # Value is 0.0, not None, for no area compute.
    if not fArea_New:
        return False, None
    relDiff = (abs(fAreaRef-fArea_New)/max((fAreaRef, fArea_New)))
    if bDebug: sEval = "relDiff"; print(sEval,'=',eval(sEval))
    return relDiff <= fDeltaTol


def cutUpSurface(surface, curves, useEdgeCurves, tolerance, fAreaRef, bDebug=False):
    """
    Parameters:
    Returns:
        (rg.Brep, None) on success
        (None, str(Fail label)) on fail
    """

    if bDebug:
        print("cutUpSurface:")

    breps_Res = rg.Brep.CutUpSurface(
        surface=surface,
        curves=curves,
        useEdgeCurves=useEdgeCurves,
        tolerance=tolerance)

    if len(breps_Res) == 0:
        #raise Exception("No result of CutUpSurface for duplicated face.")
        return None, 'CutUpSrfFail'

    if len(breps_Res) > 1:
        for _ in breps_Res: _.Dispose()
        return None, 'MultiRes'
        #[sc.doc.Objects.AddBrep(b) for b in breps_Res]
        #print("Multiple breps from CutUpSurface. They were added to model.")
        #sc.doc.Views.Redraw()
        #CutUpSrfFail.append(iF)
        #return

    rgB_1F_Out = breps_Res[0]

    if not rgB_1F_Out.IsValid:
        rgB_1F_Out.Dispose()
        return None, 'InvalidRes'

    #sc.doc.Objects.AddBrep(rgB_1F_Out); sc.doc.Views.Redraw(); 1/0

    bAreaPass = _areaComp(
        rgB_1F_Out, fAreaRef=fAreaRef, fDeltaTol=0.05, bDebug=bDebug)

    if not bAreaPass:
        rgB_1F_Out.Dispose()
        return None, 'AreaFail'

    return rgB_1F_Out, None


def _extendSrf(srf_In, extensionLength):
    if srf_In.IsClosed(0) and srf_In.IsClosed(1):
        return

    bModified = False
    srf_WIP = srf_In.Duplicate()

    if not srf_In.IsClosed(0):
        rv = srf_WIP.Extend(
            edge=rg.IsoStatus.West,
            extensionLength=extensionLength,
            smooth=True)
        if rv:
            srf_WIP.Dispose()
            srf_WIP = rv
            bModified = True
        rv = srf_WIP.Extend(
            edge=rg.IsoStatus.East,
            extensionLength=extensionLength,
            smooth=True)
        if rv:
            srf_WIP.Dispose()
            srf_WIP = rv
            bModified = True

    if not srf_In.IsClosed(1):
        rv = srf_WIP.Extend(
            edge=rg.IsoStatus.South,
            extensionLength=extensionLength,
            smooth=True)
        if rv:
            srf_WIP.Dispose()
            srf_WIP = rv
            bModified = True
        rv = srf_WIP.Extend(
            edge=rg.IsoStatus.North,
            extensionLength=extensionLength,
            smooth=True)
        if rv:
            srf_WIP.Dispose()
            srf_WIP = rv
            bModified = True
    return srf_WIP


def ncptct(crv):
    nc = crv.ToNurbsCurve()
    ct = nc.Points.Count
    nc.Dispose()
    return ct


def _createCurvesFromEdges(brep, minLength=None, bDebug=False):
    if minLength is None: minLength = 0.1 * sc.doc.ModelAbsoluteTolerance
    curves_forCUS = []
    for iE in brep.Faces[0].AdjacentEdges():
        edge = brep.Edges[iE]
        length = rg.Curve.GetLength(edge) # Zero on failure.
        if length == 0:
            continue
        if length < minLength:
            if bDebug: print("Skipped {}-long curve.".format(length))
            continue
        curve = edge.DuplicateCurve()
        curves_forCUS.append(curve)
    return curves_forCUS


def _retrimFace(rgFace_In, fTol_Edge_and_trim, bRebuildEdges, bStartingCrvs_Edge_not_UV, bTryShrinkOnFail, bTryExtendOnFail, bTryOtherTolsOnFail, bDebug=False):
    """
    Returns tuple of 2:
        rg.Brep (Just a duplicate of face when the retrim fails.)
        None on success, str on fail
    """

    rgBs_1F = [rgFace_In.DuplicateFace(duplicateMeshes=False)]
    # rgBs_1F  [0] is DuplicateFace of rgFace_In.  The rest are mods/WIPs.

    if bDebug:
        print("_retrimFace:")
        sEval = "sum([ncptct(rgE) for rgE in rgBs_1F[0].Edges])"; print(sEval,'=',eval(sEval))
        sEval = "sum([ncptct(rgT) for rgT in rgBs_1F[0].Trims])"; print(sEval,'=',eval(sEval))

    if bRebuildEdges and not bStartingCrvs_Edge_not_UV:
        rgBs_1F.append(rgBs_1F[0].DuplicateBrep())
        bRebuilt = rgBs_1F[-1].Faces[0].RebuildEdges(
            tolerance=fTol_Edge_and_trim,
            rebuildSharedEdges=True,
            rebuildVertices=True)

        if bDebug:
            sEval = "bRebuilt"; print(sEval,'=',eval(sEval))
            sEval = "sum([ncptct(rgE) for rgE in rgBs_1F[-1].Edges])"; print(sEval,'=',eval(sEval))
            sEval = "sum([ncptct(rgT) for rgT in rgBs_1F[-1].Trims])"; print(sEval,'=',eval(sEval))

    if bDebug: print("\n## Phase 1 (P1)")

    area_In = rgBs_1F[-1].GetArea()
    if not area_In:
        for brep in rgBs_1F[:-1]: brep.Dispose()
        return rgBs_1F[-1], 'AreaFail'

    curves_forCUS = _createCurvesFromEdges(rgBs_1F[-1])
    #[sc.doc.Objects.AddCurve(c) for c in curves_forCUS]

    rgB_1F_Res_P1, sFailType = cutUpSurface(
        surface=rgBs_1F[0].Faces[0].UnderlyingSurface(),
        curves=curves_forCUS,
        useEdgeCurves=not(bRebuildEdges and bStartingCrvs_Edge_not_UV),
        tolerance=fTol_Edge_and_trim,
        fAreaRef=area_In,
        bDebug=bDebug,
        )

    # useEdgeCurves conditional explained:
    #    bRebuildEdges bStartingCrvs_Edge_not_UV : result
    #    T T (Start with edges) : F (Use edges pushed up from new BrepTrims.)
    #    T F (Start with trims) : T (Use provided curves as edges. curves are resultant edges of BrepFace.RebuildEdges of input BrepFace (not CutUpSurface output).)
    #    F T (Start with edges) : T (Use provided curves as edges. curves are edges of input BrepFace.)
    #    F F (Start with trims) : T (Use provided curves as edges. curves are edges of input BrepFace.)


    if not rgB_1F_Res_P1:
        if not (bTryShrinkOnFail or bTryExtendOnFail):
            for brep in rgBs_1F[:-1]: brep.Dispose()
            return rgBs_1F[-1], sFailType


        def try_cutUpSurface_on_modified_srf():

            rgB_1F_WIP = rgBs_1F[-1].Duplicate()

            if bTryShrinkOnFail:
                if rgB_1F_WIP.Faces[0].ShrinkFace(
                    disableSide=rg.BrepFace.ShrinkDisableSide.ShrinkAllSides
                ):
                    rgB_1F_Res_P1, _ = cutUpSurface(
                        surface=rgB_1F_WIP.Faces[0].UnderlyingSurface(),
                        curves=curves_forCUS,
                        useEdgeCurves=True,
                        tolerance=fTol_Edge_and_trim,
                        fAreaRef=area_In,
                        bDebug=bDebug,
                        )
                    if rgB_1F_Res_P1:
                        rgB_1F_WIP.Dispose()
                        return rgB_1F_Res_P1

            if not bTryExtendOnFail:
                rgB_1F_WIP.Dispose()
                return

            rgS_Extended = _extendSrf(
                rgB_1F_WIP.Faces[0].UnderlyingSurface(),
                10.0*fTol_Edge_and_trim)

            rgB_1F_WIP.Dispose()

            if rgS_Extended is None:
                return 

            rgB_1F_Res_P1, _ = cutUpSurface(
                surface=rgS_Extended,
                curves=curves_forCUS,
                useEdgeCurves=True,
                tolerance=fTol_Edge_and_trim,
                fAreaRef=area_In,
                bDebug=bDebug,
                )
            if rgB_1F_Res_P1:
                return rgB_1F_Res_P1

        rgB_1F_Res_P1 = try_cutUpSurface_on_modified_srf()
        if not rgB_1F_Res_P1:
            for brep in rgBs_1F[:-1]: brep.Dispose()
            return rgBs_1F[-1], sFailType

    if bDebug:
        sEval = "sum([ncptct(rgE) for rgE in rgB_1F_Res_P1.Edges])"; print(sEval,'=',eval(sEval))
        sEval = "sum([ncptct(rgT) for rgT in rgB_1F_Res_P1.Trims])"; print(sEval,'=',eval(sEval))

    if rgBs_1F[0].Faces[0].PerFaceColor is not None:
        rgB_1F_Res_P1.Faces[0].PerFaceColor = rgBs_1F[0].Faces[0].PerFaceColor


    # Postpone Phase 2 until further investigation seems warranted.

    for brep in rgBs_1F: brep.Dispose()
    return rgB_1F_Res_P1, None



    if not bRebuilt:
        for brep in rgBs_1F: brep.Dispose()
        return rgB_1F_Res_P1, None


    if bDebug: print("\n## Phase 2: cutUpSurface the result of Phase 1 to attempt to further simplify BrepTrim curves.")

    rgF_Res_P1 = rgB_1F_Res_P1.Faces[0]
    rgS_Res_P1 = rgF_Res_P1.UnderlyingSurface()
    rgCs_Res_P1 = [rgB_1F_Res_P1.Edges[iE].DuplicateCurve() for iE in rgF_Res_P1.AdjacentEdges()]

    rgB_1F_Res_P2, sFailType = cutUpSurface(
        surface=rgS_Res_P1,
        curves=rgCs_Res_P1,
        useEdgeCurves=True,
        tolerance=fTol_Edge_and_trim,
        fAreaRef=area_In,
        bDebug=bDebug,
        )

    if rgB_1F_Res_P2 is None:
        sc.doc.Objects.AddBrep(rgB_1F_Res_P1)
        raise Exception("Phase 2 failure. Phase 1 result was added.")
        #return rgB_1Face, 'CutUpSrfFail'

    rgB_1F_Res_P1.Dispose()

    if rgBs_1F[0].Faces[0].PerFaceColor is not None:
        rgB_1F_Res_P2.Faces[0].PerFaceColor = rgBs_1F[0].Faces[0].PerFaceColor

    if bDebug:
        sEval = "sum([ncptct(rgE) for rgE in rgB_1F_Res_P2.Edges])"; print(sEval,'=',eval(sEval))
        sEval = "sum([ncptct(rgT) for rgT in rgB_1F_Res_P2.Trims])"; print(sEval,'=',eval(sEval))

    return rgB_1F_Res_P2, None


def processBrepObject(rdBrep, fTol_Edge_and_trim, fTol_Join, bRebuildEdges, bStartingCrvs_Edge_not_UV, bTryShrinkOnFail, bTryExtendOnFail, bTryOtherTolsOnFail, bReplace, bEcho=True, bDebug=False):
    """
    """

    rdB_Full_In = rdBrep
    rgB_Full_In = rdB_Full_In.BrepGeometry

    #if rgB_Full_In.Faces.Count == 1:
    #    rdBs_1F = [rdB_Full_In]
    #else:
    #    rdBs_1F = rdB_Full_In.GetSubObjects()
    #    if not rdBs_1F:
    #        print("No GetSubObjects result for {}.".format(rdB_Full_In.Id))
    #        return

    rgBs_1F_Res = []
    rgBs_toManuallyRetrim = []

    dict_iFs = {
        'Success': [],
        'AreaFail': [],
        'CutUpSrfFail': [],
        'MultiRes': [],
        'InvalidRes': [],
        }


    #sLogs = []


    nCPs_Es_Pre = nCPs_Es_Post = nCPs_Ts_Pre = nCPs_Ts_Post = 0
    nEs_Pre = nEs_Post = nTs_Pre = nTs_Post = 0

    iDivision = 20
    iFs_atDivision = [int(round((1.0/iDivision)*i*rgB_Full_In.Faces.Count,0)) for i in range(iDivision)]
    iCounter = 0

    for iF, rgF_In in enumerate(rgB_Full_In.Faces):
        sc.escape_test()
        Rhino.RhinoApp.Wait()

        if rgB_Full_In.Faces.Count == 1:
            Rhino.RhinoApp.SetCommandPromptMessage("Processing face...")
        elif rgB_Full_In.Faces.Count < iDivision:
            Rhino.RhinoApp.SetCommandPromptMessage(
                "Processing face {} of {} ({}% complete)...".format(
                    iF+1,
                    rgB_Full_In.Faces.Count,
                    int(100*(iF+1)/rgB_Full_In.Faces.Count),
                    ))
        elif iF in iFs_atDivision:
            iCounter += 1
            Rhino.RhinoApp.SetCommandPromptMessage(
                "Processing face {} of {} ({}% complete)...".format(
                    iF+1,
                    rgB_Full_In.Faces.Count,
                    int(100*(iF+1)/rgB_Full_In.Faces.Count),
                    ))
            #Rhino.RhinoApp.SetCommandPromptMessage(
            #    "Processing face {} of {}...".format(iF+1, rgB_Full_In.Faces.Count))

        rgB_1F_Res, sFailType = _retrimFace(
            rgF_In,
            fTol_Edge_and_trim=fTol_Edge_and_trim,
            bRebuildEdges=bRebuildEdges,
            bStartingCrvs_Edge_not_UV=bStartingCrvs_Edge_not_UV,
            bTryShrinkOnFail=bTryShrinkOnFail,
            bTryExtendOnFail=bTryExtendOnFail,
            bTryOtherTolsOnFail=bTryOtherTolsOnFail,
            bDebug=bDebug,
            )

        if sFailType:
            rgBs_toManuallyRetrim.append(rgB_1F_Res)
            #rgBs_toManuallyTrim.append(createBrepOfUntrimmedFace(rgF_In))
            #rdB_1F.Dispose()
            dict_iFs[sFailType].append(iF)
            continue

        rgBs_1F_Res.append(rgB_1F_Res)
        dict_iFs['Success'].append(iF)

        nCPs_Es_Pre += sum([ncptct(rgB_Full_In.Edges[iE]) for iE in rgF_In.AdjacentEdges()])
        nCPs_Es_Post += sum([ncptct(rgE) for rgE in rgB_1F_Res.Edges])
        nCPs_Ts_Pre += sum([ncptct(rgB_Full_In.Trims[rgT.TrimIndex]) for rgL in rgF_In.Loops for rgT in rgL.Trims])
        nCPs_Ts_Post += sum([ncptct(rgT) for rgT in rgB_1F_Res.Trims])

        nEs_Pre += rgF_In.AdjacentEdges().Count
        nEs_Post += rgB_1F_Res.Edges.Count
        nTs_Pre += sum([loop.Trims.Count for loop in rgF_In.Loops])
        nTs_Post += rgB_1F_Res.Trims.Count


    if len(rgBs_1F_Res) == 0:
        pass
    elif len(rgBs_1F_Res) == 1:
        if bReplace:
            sc.doc.Objects.Replace(
                objectId=rdB_Full_In.Id,
                brep=rgBs_1F_Res[0])
        else:
            sc.doc.Objects.AddBrep(
                rgBs_1F_Res[0]
                )
    else:
        breps_Joined = rg.Brep.JoinBreps(
            brepsToJoin=rgBs_1F_Res,
            tolerance=fTol_Join)
        
        if len(breps_Joined) == 0:
            raise Exception("No result from JoinBreps.")
        if len(breps_Joined) > 1:
            print("Multiple breps output of JoinBreps.")
            gB_Outs = [
                sc.doc.Objects.AddBrep(
                    brep_Joined, attributes=rdB_Full_In.Attributes)
                    for brep_Joined in breps_Joined]
            if bReplace and not any([gB_Out == Guid.Empty for gB_Out in gB_Outs]):
                sc.doc.Objects.Delete(obj=rdB_Full_In, quiet=False)
        else:
            if bReplace:
                sc.doc.Objects.Replace(
                    objectId=rdB_Full_In.Id,
                    brep=breps_Joined[0])
            else:
                sc.doc.Objects.AddBrep(
                    breps_Joined[0]
                    )

    gB_toManuallyTrim = [
        sc.doc.Objects.AddBrep(
            rgB, attributes=rdB_Full_In.Attributes) for rgB in rgBs_toManuallyRetrim]
    if gB_toManuallyTrim:
        if not all(gB==Guid.Empty for gB in gB_toManuallyTrim):
            print("Added {} untrimmed faces. Trim this manually.".format(
                len(gB_toManuallyTrim) - gB_toManuallyTrim.count(Guid.Empty)))
        if Guid.Empty in gB_toManuallyTrim:
            print("{} surfaces could not be added. Check results.".format(
                gB_toManuallyTrim.count(Guid.Empty)))

    print("Retrimmed {} out of {} faces.".format(
        len(dict_iFs['Success']), rgB_Full_In.Faces.Count))
    for key in 'AreaFail', 'CutUpSrfFail', 'MultiRes', 'InvalidRes':
        if dict_iFs[key]:
            print("{}: {}".format(key, len(dict_iFs[key])))
    #print("Trim fail count: {}".format(
    #    len(dict_iFs['CutUpSrfFail'])))
    #print("Invalid brep result (skipped) count: {}".format(
    #    len(dict_iFs['InvalidRes'])))
    #if sLogs:
    #    print("Logs of invalid breps:")
    #    print(*sLogs)

    if nEs_Post != nEs_Pre:
        print("Edge count (for individual faces (exploded polysrf)): {} -> {}, Delta {:+}.".format(
            nEs_Pre,
            nEs_Post,
            nEs_Post-nEs_Pre))
    else:
        print("Edge CP count (for individual faces (exploded polysrf)): {} -> {}, Delta {:+}.".format(
            nCPs_Es_Pre,
            nCPs_Es_Post,
            nCPs_Es_Post-nCPs_Es_Pre))

    if nTs_Post != nTs_Pre:
        print("Trim count: {} -> {}, Delta {:+}.".format(
            nTs_Pre,
            nTs_Post,
            nTs_Post-nTs_Pre))
    else:
        print("Trim CP count: {} -> {}, Delta {:+}.".format(
            nCPs_Ts_Pre,
            nCPs_Ts_Post,
            nCPs_Ts_Post-nCPs_Ts_Pre))


def main():
    
    #res, objrefs = ri.RhinoGet.GetMultipleObjects(
    #    "Select breps",
    #    acceptNothing=False,
    #    filter=rd.ObjectType.Brep)
    #if res != Rhino.Commands.Result.Success: return

    objrefs = getInput()
    if objrefs is None: return

    fTol_Edge_and_trim = Opts.values['fTol_Edge_and_trim']
    fTol_Join = Opts.values['fTol_Join']
    bRebuildEdges = Opts.values['bRebuildEdges']
    bStartingCrvs_Edge_not_UV = Opts.values['bStartingCrvs_Edge_not_UV']
    bTryShrinkOnFail = Opts.values['bTryShrinkOnFail']
    bTryExtendOnFail = Opts.values['bTryExtendOnFail']
    bTryOtherTolsOnFail = Opts.values['bTryOtherTolsOnFail']
    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']
    
    for objref in objrefs:
        rdBrep_In = objref.Object()
        processBrepObject(
            rdBrep_In,
            fTol_Edge_and_trim=fTol_Edge_and_trim,
            fTol_Join=fTol_Join,
            bRebuildEdges=bRebuildEdges,
            bStartingCrvs_Edge_not_UV=bStartingCrvs_Edge_not_UV,
            bTryShrinkOnFail=bTryShrinkOnFail,
            bTryExtendOnFail=bTryExtendOnFail,
            bTryOtherTolsOnFail=bTryOtherTolsOnFail,
            bReplace=bReplace,
            bEcho=bEcho,
            bDebug=bDebug,
            )

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()