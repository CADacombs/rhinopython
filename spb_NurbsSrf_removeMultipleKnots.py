"""
NurbsSurfaceKnotList.RemoveMultipleKnots appears to not remove multiplicities
from full (== degree) multiplicities.
_RemoveMultiKnot does the same (at least up through version 7.25).
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221224,26: WIP: Created, starting with another script.

TODO: Filter out conical cross-section directions.
      Add support for periodic surfaces.
"""


import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import xBrep_getDistancesBetween2
import xBrepFace
import xBrepObject


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAllowBrepSelection'; keys.append(key)
    values[key] = False
    names[key] = 'SelectionMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'FacesOnly', 'FacesAndBreps')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitDev'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSrfOnFail'; keys.append(key)
    values[key] = True
    names[key] = 'AddSrfOnFail'
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

        if key == 'fDevTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

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


def _sortBrepsAndFaces(objrefs):
    """
    Parameters:
        list(objrefs)
    Returns:
        list(Brep GUIDs)
        list(lists(integers of Face indices) per brep)
    """
        
    gBs = []
    rdBs = []
    idxFs_perB = []
    
    for o in objrefs:
        gB = o.ObjectId
        rdB = o.Object()
        rgB = o.Brep()

        if not rgB.IsValid:
            print("Brep {} is invalid.  Fix first.".format(gB))
            continue

        idx_CompIdx = o.GeometryComponentIndex.Index
        if idx_CompIdx == -1:
            if gB in gBs:
                idxFs_perB[gBs.index(gB)] = range(rgB.Faces.Count)
            else:
                gBs.append(gB)
                rdBs.append(rdB)
                idxFs_perB.append(range(rgB.Faces.Count))
        else:
            rgFace_Brep0 = o.Face()
            if gB in gBs:
                if rgFace_Brep0 in idxFs_perB[gBs.index(gB)]:
                    continue
                else:
                    idxFs_perB[gBs.index(gB)].append(rgFace_Brep0.FaceIndex)
            else:
                gBs.append(gB)
                rdBs.append(rdB)
                idxFs_perB.append([rgFace_Brep0.FaceIndex])

    return rdBs, idxFs_perB


def getInput():
    """
    Get BrepFaces with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select faces")

    go.SetCommandPromptDefault("All normal breps when none are selected")

    go.AcceptNothing(True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    #go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        if Opts.values['bAllowBrepSelection']:
            go.SetCommandPrompt("Select breps and/or faces")
            go.GeometryFilter = rd.ObjectType.Brep | rd.ObjectType.Curve
        else:
            go.SetCommandPrompt("Select faces")
            go.GeometryFilter = rd.ObjectType.Surface

        go.ClearCommandOptions()

        idxs_Opt.clear()

        go.AcceptNumber(Opts.values['bLimitDev'], acceptZero=Opts.values['bLimitDev'])

        addOption('bAllowBrepSelection')
        addOption('bLimitDev')
        if Opts.values['bLimitDev']:
            addOption('fDevTol')
        addOption('bReplace')
        if Opts.values['bReplace']:
            addOption('bExtract')
        addOption('bSrfOnFail')
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
            return _sortBrepsAndFaces(objrefs)

        if res == ri.GetResult.Nothing:
            oes = rd.ObjectEnumeratorSettings()
            oes.NormalObjects = True
            oes.LockedObjects = False
            oes.IncludeLights = False
            oes.IncludeGrips = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep

            rdBs = list(sc.doc.Objects.GetObjectList(oes))
            go.Dispose()
            if len(rdBs) == 0: return

            return (
                rdBs,
                [[iF for iF in range(rdBs[iB].BrepGeometry.Faces.Count)]
                for iB in range(len(rdBs))]
                )

        if res == ri.GetResult.Number:
            if Opts.values['bLimitDev']:
                key = 'fDevTol'
                Opts.riOpts[key].CurrentValue = go.Number()
            else:
                continue

            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def knotMultiplicityList(knots):
    """Returns a list."""
    i = 0
    iMulties = []
    fKnotTs_Unique = []
    while i < knots.Count:
        sc.escape_test()
        knot = knots[i]
        fKnotTs_Unique.append(knot)
        iMulti = knots.KnotMultiplicity(index=i)
        iMulties.append(iMulti)
        #print("{} at {:.4f}".format(iMulti, knot),
        i += iMulti
    return iMulties


def multiplicityRangeOfInteriorKnots(ns, iDir):
    degree = ns.Degree(iDir)
    knots = ns.KnotsV if iDir else ns.KnotsU
    ms = []
    iK = degree
    while iK < knots.Count-degree:
        sc.escape_test()
        m = knots.KnotMultiplicity(iK)
        ms.append(m)
        iK += m
    return ms


def removeMultipleKnots(ns_In, bDebug=False):
    """
    This is a replacement for NurbsSurfaceKnotList.RemoveMultipleKnots,
    also reducing full (== degree) multiplicities to 1.
    """

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script only works in Rhino V6 and above.")
        return

    if not isinstance(ns_In, rg.NurbsSurface):
        if bDebug: print("{} skipped.".format(ns_In.GetType().Name))
        return

    if ns_In.SpanCount(0) == 1 and ns_In.SpanCount(1) == 1:
        if bDebug: print("No interior knots.")
        return

    knots_In = ns_In.KnotsU, ns_In.KnotsV

    ms_U = knotMultiplicityList(ns_In.KnotsU)
    ms_V = knotMultiplicityList(ns_In.KnotsV)

    if (
        (len(ms_U[1:-1]) == 0 or max(ms_U[1:-1]) == 1) and
        (len(ms_V[1:-1]) == 0 or max(ms_V[1:-1]) == 1)
        ):
        if bDebug: print("No interior knots with multiplicity above 1.")
        return

    ns_Out = ns_In.Duplicate()

    knots_Out = ns_Out.KnotsU, ns_Out.KnotsV

    # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Collections_NurbsSurfaceKnotList_RemoveKnots.htm
    # Remove knots from the knot vector and adjusts the remaining control points to maintain surface position as closely as possible.
    # The knots from Knots[index0] through Knots[index1 - 1] will be removed.

    for iDir in (0,1):
        knots = knots_Out[iDir]
        degree = ns_Out.Degree(iDir)
        iK = knots.Count - degree - 1
        if bDebug: sEval = 'iK'; print("{}: {}".format(sEval, eval(sEval)))
        while iK > degree:
            sc.escape_test()
            m = knots.KnotMultiplicity(iK)
            if bDebug: sEval = 'm'; print("{}: {}".format(sEval, eval(sEval)))
            index0 = iK - m + 1
            index1 = iK
            if bDebug: sEval = 'index0'; print("{}: {}".format(sEval, eval(sEval)))
            if bDebug: sEval = 'index1'; print("{}: {}".format(sEval, eval(sEval)))
            knots_Out[iDir].RemoveKnots(index0, index1)

            if bDebug: sEval = 'ns_Out.KnotsU.Count'; print("{}: {}".format(sEval, eval(sEval)))
            if bDebug: sEval = 'ns_Out.KnotsV.Count'; print("{}: {}".format(sEval, eval(sEval)))
            iK -= m
            if bDebug: sEval = 'iK'; print("{}: {}".format(sEval, eval(sEval)))

    return ns_Out


def createSurface(ns_In, fDevTol, bDebug=False):
    """
    Returns on success:
        rg.NurbsSurface
        float(surface deviation)
        None
    Returns on deviation fail:
        None
        float(surface deviation needed)
        None
    Returns on other fails:
        None
        None
        str(Reason of failure)
    """

    if Rhino.RhinoApp.ExeVersion < 6:
        raise Exception("This script only works in Rhino V6 and above.")


    if not isinstance(ns_In, rg.NurbsSurface):
        return None, None, "{} skipped.".format(ns_In.GetType().Name)

    if ns_In.SpanCount(0) == 1 and ns_In.SpanCount(1) == 1:
        return None, None, "NurbsSurface has no interior knots."

    if bDebug:
        sEval = 'ns_In.KnotsU.Count'; print("{}: {}".format(sEval, eval(sEval)))
        sEval = 'ns_In.KnotsV.Count'; print("{}: {}".format(sEval, eval(sEval)))

    ms_U = knotMultiplicityList(ns_In.KnotsU)
    ms_V = knotMultiplicityList(ns_In.KnotsV)

    if (
        (len(ms_U[1:-1]) == 0 or max(ms_U[1:-1]) == 1) and
        (len(ms_V[1:-1]) == 0 or max(ms_V[1:-1]) == 1)
        ):
        return None, None, "No interior knots with multiplicity above 1."


    ns_Out = removeMultipleKnots(ns_In, bDebug)


    rgMeshParams = rg.MeshingParameters.QualityRenderMesh

    rc = xBrep_getDistancesBetween2.getDistancesBetweenBreps(
        ns_In,
        ns_Out,
        rgMeshParams,
        bCalcBrepIntersection=False)

    if rc[0]:
        rgSrf1_LastGood = ns_Out.Duplicate()
        fDev = rc[1]


    def getNurbsSurfaceChangeDescription(rgNurbsSrf1, rgNurbsSrf2):
        s  = "  Prop:I->O"
        s += "  {}:({}->{})x({}->{})".format(
                "Deg",
                rgNurbsSrf1.OrderU-1,
                rgNurbsSrf2.OrderU-1,
                rgNurbsSrf1.OrderV-1,
                rgNurbsSrf2.OrderV-1,
        )
        s += "  {}:({}->{})x({}->{})".format(
                "PtCt",
                rgNurbsSrf1.Points.CountU,
                rgNurbsSrf2.Points.CountU,
                rgNurbsSrf1.Points.CountV,
                rgNurbsSrf2.Points.CountV,
        )
        if Rhino.RhinoApp.ExeVersion >= 7:
            s += "  {}:({}->{})x({}->{})".format(
                    "IsUniform",
                    str(isKnotVectorUniform(rgNurbsSrf1.KnotsU))[0],
                    str(isKnotVectorUniform(rgNurbsSrf2.KnotsU))[0],
                    str(isKnotVectorUniform(rgNurbsSrf1.KnotsV))[0],
                    str(isKnotVectorUniform(rgNurbsSrf2.KnotsV))[0],
            )
        s += "  {}:{}->{}".format(
                "IsRational",
                str(rgNurbsSrf1.IsRational)[0],
                str(rgNurbsSrf2.IsRational)[0],
        )
        s += "  {}:({}->{})x({}->{})".format(
                "IsClosed",
                str(rgNurbsSrf1.IsClosed(0))[0],
                str(rgNurbsSrf2.IsClosed(0))[0],
                str(rgNurbsSrf1.IsClosed(1))[0],
                str(rgNurbsSrf2.IsClosed(1))[0],
        )
        if (    rgNurbsSrf1.IsClosed(0) or rgNurbsSrf1.IsClosed(1) or
                rgNurbsSrf2.IsClosed(0) or rgNurbsSrf2.IsClosed(1)
        ):
            s += "  {}:{}->{}".format(
                    "IsPeriodic",
                    str(rgNurbsSrf1.IsPeriodic)[0],
                    str(rgNurbsSrf2.IsPeriodic)[0],
            )
        return s


    if fDevTol is None or (fDev <= fDevTol):

        if bDebug:
            rgBrep_1F_Mod = rgBs_1F_Res[0]
            s  = getNurbsSurfaceChangeDescription(ns_In, ns_Out)
            s += "  Deviation: {0:.2e}".format(fDev)
            print(s)

        return ns_Out, fDev, None

    ns_Out.Dispose()

    return None, fDev, None


def addBrepOfSrf(rgSrf, bDebug=False):
    gB_Out = sc.doc.Objects.AddSurface(rgSrf)
    if bDebug:
        if gB_Out.Guid != gB_Out.Empty:
            print("Converted underlying surface was added.")
        else:
            print("Converted underlying surface could not be added.")


def create1FaceBrepWithNewSurface(rgFace_In, **kwargs):
    """
    Parameters:
        rgFace_In
        fDevTol
        bSrfOnFail
        bDebug
    Returns on success:
        rg.Brep
        float(surface deviation)
        None
    Returns on deviation fail:
        None
        float(surface deviation needed)
        None
    Returns on other fails:
        None
        None
        str(Reason of failure)
    """

    if Rhino.RhinoApp.ExeVersion < 6:
        raise Exception("This script only works in Rhino V6 and above.")


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bExtract = getOpt('bExtract')
    bSrfOnFail = getOpt('bSrfOnFail')
    bDebug = getOpt('bDebug')


    rgSrf_In = rgFace_In.UnderlyingSurface()

    if not isinstance(rgSrf_In, rg.NurbsSurface):
        return None, None, "Skipped {}.".format(rgSrf_In.GetType().Name)


    rc = createSurface(
        rgSrf_In,
        fDevTol,
        bDebug=bDebug,
        )
    if rc[0] is None: return rc

    ns_Res, fDev, sLog = rc

    # Success in creating NurbsSurface.  Now, create correctly trimmed brep.
    
    rgB_1F_In = rgFace_In.DuplicateFace(duplicateMeshes=False)

    rgB_1F_In.Faces[0].RebuildEdges(
        tolerance=max((1e-6, 0.1*sc.doc.ModelAbsoluteTolerance)),
        rebuildSharedEdges=True,
        rebuildVertices=True)
    #rgBrep0_1Face.Faces.ShrinkFaces()

    if rgB_1F_In.IsSurface:
        rgB_1F_In.Dispose()
        rgB_1F_Out = ns_Res.ToBrep()
        ns_Res.Dispose()
        if not rgB_1F_Out.IsValid:
            rgB_1F_Out.Dispose()
            if bSrfOnFail: addBrepOfSrf(ns_Res, bDebug)
            return None, None, "Invalid brep geometry after ToBrep."
        return rgB_1F_Out, fDev, None
    
    # Test areas before trimming.
    fArea_Trimmed = rgB_1F_In.GetArea()
    if fArea_Trimmed:
        rgBrep_Untrimmed = rgB_1F_In.Faces[0].UnderlyingSurface().ToBrep()
        fArea_Untrimmed = rgBrep_Untrimmed.GetArea()
        rgBrep_Untrimmed.Dispose()
        if fArea_Untrimmed:
            if abs(fArea_Trimmed - fArea_Untrimmed) <= sc.doc.ModelAbsoluteTolerance**2:
                rgB_1F_Out = ns_Res.ToBrep()
                ns_Res.Dispose()
                if not rgB_1F_Out.IsValid:
                    rgB_1F_Out.Dispose()
                    if bSrfOnFail: addBrepOfSrf(ns_Res, bDebug)
                    return None, None, "Invalid brep geometry after ToBrep."
                return rgB_1F_Out, fDev, None

    rgB_1F_Out = xBrepFace.retrimFace(
            rgB_1F_In.Faces[0],
            rgSrf_Replacement=ns_Res,
            fSplitTol=1.0*sc.doc.ModelAbsoluteTolerance if fDevTol is None else fDevTol,
            bDebug=bDebug
    )
    rgB_1F_In.Dispose()
    ns_Res.Dispose()

    if rgB_1F_Out is None:
        if bSrfOnFail: addBrepOfSrf(ns_Res, bDebug)
        return None, None, "xBrepFace.createMonofaceBrep returned None."

    if not rgB_1F_Out.IsValid:
        rgB_1F_Out.Dispose()
        if bSrfOnFail: addBrepOfSrf(ns_Res, bDebug)
        return None, None, "An invalid brep was skipped."

    return rgB_1F_Out, fDev, None


def processBrepObject(rhBrep, idxFaces=None, **kwargs):
    """
    Parameters:
        rhBrep: rd.BrepObject or GUID of brep.
        idxFaces
        fDevTol
        bReplace
        bExtract
        bSrfOnFail
        bEcho
        bDebug
    """

    if Rhino.RhinoApp.ExeVersion < 6:
        raise Exception("This script only works in Rhino V6 and above.")


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fDevTol = getOpt('fDevTol')
    bReplace = getOpt('bReplace')
    bExtract = getOpt('bExtract')
    bSrfOnFail = getOpt('bSrfOnFail')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def isKnotVectorUniform(knots):
        return (
            (knots.KnotStyle == rg.KnotStyle.Uniform) or
            (knots.KnotStyle == rg.KnotStyle.QuasiUniform) or
            (
                (knots.KnotStyle == rg.KnotStyle.PiecewiseBezier) and
                knots.Count == knots.KnotMultiplicity(0) * 2)
            )


    rgBs_1F_Res = []
    idxsFs_Success = []
    fDevs = []
    fDevs_Needed = []
    sLogs = []


    if isinstance(rhBrep, rd.BrepObject):
        rdB_In = rhBrep
    else:
        rdB_In = rs.coercerhinoobject(rhBrep)
    rgB_In = rdB_In.BrepGeometry

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt
    
    idxs_AtTenths = [int(round(0.1*i*len(idxFaces),0)) for i in range(10)]
    
    for iF, idx_rgFace in enumerate(idxFaces):
        if sc.escape_test(False):
            raise Exception("Searching interrupted by user.")

        if iF in idxs_AtTenths:
            s = sCmdPrompt_In + ", {:d}% of {} faces in current brep ...".format(
                int(100.0 * (iF+1) / len(idxFaces)), len(idxFaces))
            
            if bDebug:
                print(s)
            else:
                Rhino.RhinoApp.SetCommandPrompt(s)

        rgFace_In = rgB_In.Faces[idx_rgFace]


        rc = create1FaceBrepWithNewSurface(
            rgFace_In,
            fDevTol=fDevTol,
            bSrfOnFail=bSrfOnFail,
            bDebug=bDebug
            )

        rgB_1F_Res, fDev, sLog = rc

        if rgB_1F_Res is None:
            if fDev is not None:
                fDevs_Needed.append(fDev)
            if sLog is not None:
                sLogs.append(sLog)
            continue

        rgBs_1F_Res.append(rgB_1F_Res)
        idxsFs_Success.append(idx_rgFace)
        fDevs.append(fDev)


    if not rgBs_1F_Res:
        return [], [], fDevs_Needed, sLogs


    gBs_NewFs = []
    gBs_Res = []


    if not bReplace:
        # Adding the new face over each old one.
        gBrep1_thisBrep = []
        for rgBrep_1F_New, idxF in zip(rgBs_1F_Res, idxsFs_Success):
            gBrep1 = sc.doc.Objects.AddBrep(rgBrep_1F_New)
            if gBrep1 != gBrep1.Empty:
                gBrep1_thisBrep.append(gBrep1)
        gBs_Res.append(gBrep1_thisBrep)
        
        if bEcho: print("{} monofaces added.".format(len(gBs_Res)))
        
        if bExtract:
            rc = xBrepObject.extractFaces(rdB_In, idxsFs_Success)
            gBs_1Fs_Extracted, gBs_RemainingGeom = rc
            if bEcho:
                print("{} monofaces extracted.".format(len(gBs_1Fs_Extracted)))
                print("{} breps of unmodified faces remain.".format(len(gBs_RemainingGeom)))
    else:
        # bReplace==True.
        if bExtract:
            rc = xBrepObject.replaceFaces(
                rdB_In,
                idxsFs_Success,
                rgBs_1F_Res,
                bExtract=True,
                fTolerance_Join=max(fDevs))
            if rc:
                gBs_NewFaces, gBs_RemainingGeom = rc
                if bEcho:
                    print("{} monofaces extracted.".format(len(gBs_NewFaces)))
                    print("{} breps of unmodified faces remain.".format(len(gBs_RemainingGeom)))
                gBs_Res = gBs_NewFaces + gBs_RemainingGeom
        else:
            fTols_Edges = [edge.Tolerance for edge in rgB_In.Edges]
            fTolerance_Join = max((
                    2.0*fDevTol if fDevTol is not None else 0.0,
                    1.1*max(fTols_Edges),
                    sc.doc.ModelAbsoluteTolerance))
            rc = xBrepObject.replaceFaces(
                    rdB_In,
                    idxsFs_Success,
                    rgBs_1F_Res,
                    bExtract=False,
                    fTolerance_Join=fTolerance_Join)
            if rc:
                gBs_Res = rc
                if bEcho:
                    print("Brep was replaced with {} revised faces.".format(
                        len(idxsFs_Success)))
            else:
                if bEcho:
                    print("Brep could not be replaced with {} revised faces.".format(
                        len(idxsFs_Success)))

    for brep in rgBs_1F_Res: brep.Dispose()
    rgB_In.Dispose()

    return gBs_Res, fDevs, fDevs_Needed, sLogs


def formatDistance(fDistance):
    """Returns: str"""
    if fDistance is None:
        return "(None)"
    elif fDistance == 0.0:
        return "exactly 0".format(fDistance)
    elif fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-2)):
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def main():

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script works only in Rhino V6 and above.")
        return

    rc = getInput()
    if rc is None: return

    rdBs_In, idxFs_PerB = rc
    if not rdBs_In: return

    fDevTol = Opts.values['fDevTol'] if Opts.values['bLimitDev'] else None
    bExtract = Opts.values['bExtract']
    bReplace = Opts.values['bReplace']
    bSrfOnFail = Opts.values['bSrfOnFail']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    bDocModified = False

    gBs_Res_perBs_In = []
    fDevs_All = []
    fDevs_Needed_All = []
    sLogs_All = []


    if len(rdBs_In) == 1:
        s = "Processing brep"
    else:
        idxs_AtTenths = [int(round(0.1*i*len(rdBs_In),0)) for i in range(10)]


    for iB, (rdB_In, idxFs) in enumerate(zip(rdBs_In, idxFs_PerB)):
        if len(rdBs_In) > 1:
            if iB in idxs_AtTenths:
                s = "Processing at {:d}% of {} breps".format(
                    int(100.0 * (iB+1) / len(rdBs_In)), len(rdBs_In))
        
            if bDebug:
                print(s)
            else:
                Rhino.RhinoApp.SetCommandPrompt(s)


        rc = processBrepObject(
            rdB_In,
            idxFs,
            fDevTol=fDevTol,
            bExtract=bExtract,
            bSrfOnFail=bSrfOnFail,
            bDebug=bDebug)
        if rc[0] is None: return

        gBs_Res, fDevs, fDevs_Needed, sLogs = rc

        gBs_Res_perBs_In.extend(gBs_Res)
        fDevs_All.extend(fDevs)
        fDevs_Needed_All.extend(fDevs_Needed)
        sLogs_All.extend(sLogs)


    sc.doc.Views.RedrawEnabled = True


    if fDevs_All:
        if len(fDevs_All) == 1:
            print("Converted surfaces with deviation {}.".format(
                formatDistance(fDevs_All[0])))
        else:
            print("Converted surfaces with deviations {} through {}.".format(
                formatDistance(min(fDevs_All)),
                formatDistance(max(fDevs_All))))

    if fDevs_Needed_All:
        if len(fDevs_Needed_All) == 1:
            s += "Need tolerance of {}.".format(
                formatDistance(fDevs_Needed_All[0]))
            s += " to convert surface."
        else:
            s = "Need tolerances of {} through {}".format(
                formatDistance(min(fDevs_Needed_All)),
                formatDistance(max(fDevs_Needed_All)))
            s += " to convert remaining convertible surfaces."
        print(s)


    for sLog in set(sLogs_All):
        print("{} of {}".format(sLogs_All.count(sLog), sLog))


    ## Select new faces.
    #gBs_1F_Added_NETList = List[Guid](gBs_1F_Added)
    #nSelected = sc.doc.Objects.Select(gBs_1F_Added_NETList)
        
    #if nSelected > 0:
    #    print("{} monoface breps of simplified surfaces are"
    #            " selected.".format(nSelected))


if __name__ == '__main__': main()