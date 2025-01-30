"""
This script will select grips (control points) of a single surface who's corresponding
Greville points are closest to input reference points and curves.

A case use of this is:
1. Intersection curves of the NURBS surface with another object is created.
2. The intersection curves, possibly one at a time, can be the reference input of this script.
3. Once selected, the grips can be transformed with visual clues to eliminate each intersection.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250127-30: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from clr import StrongBox
from System import Array


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAutoPickSrf'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClosestGripCtPerAnalysisPt'; keys.append(key)
    values[key] = 1
    riOpts[key] = ri.Custom.OptionInteger(1, setLowerLimit=True, limit=1)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddInterior'; keys.append(key)
    values[key] = True
    names[key] = 'AddInteriorGripsForCrvInput'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTolForSharedGrips'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    names[key] = 'DistTolForEquallySharedGrips'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

        if key == 'iClosestGripCtPerAnalysisPt':
            if cls.riOpts[key].CurrentValue == 0:
                print("Value must be 1 or more. Now set to 1.")
                cls.riOpts[key].CurrentValue = cls.values[key] = 1
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fTolForSharedGrips':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = cls.values[key] = Rhino.RhinoMath.ZeroTolerance
            else:
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


def getInput_Ref():
    """
    Get grips with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves and/or points to find grips")

    go.GeometryFilter = rd.ObjectType.Point | rd.ObjectType.Curve

    #go.EnablePreSelect(True, ignoreUnacceptablePreselectedObjects=False)

    go.AcceptNumber(True, acceptZero=False)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        #key = 'AllSrfGripsOn'; idxs_Opt[key] = go.AddOption(key)
        #addOption('bAutoPickSrf')
        addOption('iClosestGripCtPerAnalysisPt')
        addOption('bAddInterior')
        addOption('fTolForSharedGrips')
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
            key = 'iClosestGripCtPerAnalysisPt'
            Opts.riOpts[key].CurrentValue = abs(int(go.Number()))
            Opts.setValue(key)
            continue

        if res == ri.GetResult.Number:
            key = 'fTolForSharedGrips'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _find_indices_of_minimums(lst, levels, tol_per_level):
    zipped_lists = zip(lst, range(len(lst)))
    zipped_lists.sort()
    level = 1 # Base 1.
    idxs_Out = [zipped_lists[0][1]]
    dist_Start_of_level = zipped_lists[0][0]
    for dist, idx in zipped_lists[1:]:
        if (dist - dist_Start_of_level) <= tol_per_level:
            idxs_Out.append(idx)
        else:
            level += 1
            if level > levels:
                return idxs_Out
            dist_Start_of_level = dist
            idxs_Out.append(idx)

    return idxs_Out


def _getGripsIndices_GripsClosestToPoint(rgPoint, idxs_UVs, pt3ds_Grevs, levels, tol_per_level):
    pt = rgPoint.Location

    dists = []

    for i, (iU, iV) in enumerate(idxs_UVs):
        dist = pt3ds_Grevs[i].DistanceTo(pt)
        dists.append(dist)

    return _find_indices_of_minimums(
        dists,
        levels=levels,
        tol_per_level=tol_per_level)


def _GrevillePt3ds_NestedList(ns):
    pt3ds_Grevs = []
    for iU in range(ns.Points.CountU):
        pt3ds_Grevs.append([])
        for iV in range(ns.Points.CountV):
            pt2d = ns.Points.GetGrevillePoint(iU, iV)
            pt3d = ns.PointAt(pt2d.X, pt2d.Y)
            pt3ds_Grevs[-1].append(pt3d)
    return pt3ds_Grevs


def _findMinGrevSpan(ns):
    pt3ds_Grevs = _GrevillePt3ds_NestedList(ns)
    fMinDist = float('Inf')

    for iU in range(ns.Points.CountU):
        for iV in range(ns.Points.CountV):
            if iU < (ns.Points.CountU - 1):
                dist = pt3ds_Grevs[iU][iV].DistanceTo(pt3ds_Grevs[iU+1][iV])
                if dist < fMinDist:
                    fMinDist = dist
            if iV < (ns.Points.CountV - 1):
                dist = pt3ds_Grevs[iU][iV].DistanceTo(pt3ds_Grevs[iU][iV+1])
                if dist < fMinDist:
                    fMinDist = dist

    return fMinDist


def _curve_DivideByLength(rgCurve, segmentLength, includeEnds):
    strongBox_points = StrongBox[Array[rg.Point3d]]()

    ts = rgCurve.DivideByLength(
        segmentLength=segmentLength,
        includeEnds=True,
        points=strongBox_points)

    pts = list(strongBox_points.Value)

    return pts


def _U_V_indices_from_flat_list_index(idx_In, v_count):
    return (idx_In // v_count), (idx_In % v_count)


def _nestedList_of_NurbsSrf_CPs_from_flatList(lst):
    pass


def _get_CP_indices_between_input_indices(idxs_In, countV):
    """
    Parameters:
        idxs_In: Flat list of int.

    Returns: Flat list of int.
    """

    iVs_perU = {}
    iUs_perV = {}
    for idx in idxs_In:
        iU = idx // countV
        iV = idx % countV
        if iU not in iVs_perU:
            iVs_perU[iU] = [iV]
        else:
            iVs_perU[iU].append(iV)
        if iV not in iUs_perV:
            iUs_perV[iV] = [iU]
        else:
            iUs_perV[iV].append(iU)

    if not iVs_perU and not iUs_perV:
        return []

    idxs_Added = []

    for iU in iVs_perU:
        if len(iVs_perU[iU]) == 1:
            continue
        for iV in range(min(iVs_perU[iU])+1, max(iVs_perU[iU])):
            idx_Btwn = iU*countV + iV
            if idx_Btwn not in idxs_Added:
                idxs_Added.append(idx_Btwn)

    for iV in iUs_perV:
        if len(iUs_perV[iV]) == 1:
            continue
        for iU in range(min(iUs_perV[iV])+1, max(iUs_perV[iV])):
            idx_Btwn = iU*countV + iV
            if idx_Btwn not in idxs_Added:
                idxs_Added.append(idx_Btwn)

    return idxs_Added


def _getGripsIndices_GripsClosestToCurve(rgCurve, segmentLength, pt3ds_Grevs, levels=1, tol_per_level=None):
    pts = _curve_DivideByLength(rgCurve, segmentLength=segmentLength, includeEnds=True)

    idxs_Mins_All = []

    for pt in pts:
        dists = []
        for i in range(len(pt3ds_Grevs)):
            dist = pt3ds_Grevs[i].DistanceTo(pt)
            dists.append(dist)

        idxs_Mins_This_pt = _find_indices_of_minimums(
            dists,
            levels=levels,
            tol_per_level=tol_per_level)

        idxs_Mins_All.extend(idxs_Mins_This_pt)

    return sorted(set(idxs_Mins_All))


def selectGrips(objrefs_Ref, rdBrep_1F, iClosestGripCtPerAnalysisPt=1, bAddInterior=True, fTolForSharedGrips=None, bEcho=True, bDebug=False):
    """
    """

    rgB = rdBrep_1F.BrepGeometry
    rgS = rgB.Faces[0].UnderlyingSurface()
    if not isinstance(rgS, rg.NurbsSurface):
        if bEcho:
            print("Non-NurbsSurface passed to selectGrips. Nothing will be done.")
        return

    ns = rgS

    if fTolForSharedGrips is None:
        fTolForSharedGrips = sc.doc.ModelAbsoluteTolerance

    # Creating flat lists.

    #pt2ds_Grevs = []
    pt3ds_Grevs = []
    idxs_UVs = [] # CPs & Grevilles.

    for iU in range(ns.Points.CountU):
        for iV in range(ns.Points.CountV):
            idxs_UVs.append((iU,iV))
            pt2d = ns.Points.GetGrevillePoint(iU, iV)
            #pt2ds_Grevs.append(pt2d)
            pt3d = ns.PointAt(pt2d.X, pt2d.Y)
            pt3ds_Grevs.append(pt3d)

    fMinGrevSpan = None

    idx_Grips_to_sel = []

    for objref in objrefs_Ref:
        rgO = objref.Object().Geometry
        if isinstance(rgO, rg.Point):
            idxs_Grips = _getGripsIndices_GripsClosestToPoint(
                rgO,
                idxs_UVs=idxs_UVs,
                pt3ds_Grevs=pt3ds_Grevs,
                levels=iClosestGripCtPerAnalysisPt,
                tol_per_level=fTolForSharedGrips)
            idx_Grips_to_sel.extend(idxs_Grips)
        elif isinstance(rgO, rg.Curve):
            if fMinGrevSpan is None:
                fMinGrevSpan = _findMinGrevSpan(ns)
                #sEval = "fMinGrevSpan"; print(sEval,'=',eval(sEval))
            idxs_Grips = _getGripsIndices_GripsClosestToCurve(
                rgO,
                segmentLength=0.25*fMinGrevSpan,
                pt3ds_Grevs=pt3ds_Grevs,
                levels=iClosestGripCtPerAnalysisPt,
                tol_per_level=fTolForSharedGrips)
            idx_Grips_to_sel.extend(idxs_Grips)

            if bAddInterior:
                idxs_Grips = _get_CP_indices_between_input_indices(idxs_Grips, ns.Points.CountV)
            idx_Grips_to_sel.extend(idxs_Grips)

        else:
            if bEcho:
                print("{} is no supported as reference geometry.".format(rgO))

    if not idx_Grips_to_sel:
        return

    rdGrips = None
    rdGrips_AlreadySelected = []
    bGrips_AlreadyOn = rdBrep_1F.GripsOn

    if not bGrips_AlreadyOn:
        rdBrep_1F.GripsOn = True

    rdGrips = rdBrep_1F.GetGrips()

    if bGrips_AlreadyOn:
        if rdBrep_1F.GripsSelected:
            for i, rdG in enumerate(rdGrips):
                if rdG.IsSelected(checkSubObjects=False):
                    print(i, rdG.Index)
                    rdGrips_AlreadySelected.append(rdG.Index)

    iCt_Selected = 0

    for iG in idx_Grips_to_sel:
        rv = rdGrips[iG].Select(
            on=True,
            syncHighlight=True,
            persistentSelect=True,
            ignoreGripsState=True,
            ignoreLayerLocking=False,
            ignoreLayerVisibility=False)
        if rv:
            iCt_Selected += 1

    if bEcho:
        if not iCt_Selected:
            print("No grips were selected.")
        else:
            print("{} surface grips are selected.".format(iCt_Selected))


def _findClosestSrf(objrefs):
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    rdBs_Pass = [] # Normal, 1-face, & NurbsSurface.
    for rdB in sc.doc.Objects.GetObjectList(oes):
        if rdB.BrepGeometry.Faces.Count != 1:
            continue
        rgS = rdB.BrepGeometry.Faces[0].UnderlyingSurface()
        if not isinstance(rgS, rg.NurbsSurface):
            continue
        rdBs_Pass.append(rdB)

    pts = []
    for objref in objrefs:
        rdO = objref.Object()
        rgO = rdO.Geometry
        if isinstance(rgO, rg.Point):
            pt = rgO.Location
            pts.append(pt)
        elif isinstance(rgO, rg.Curve):
            pt = rgO.PointAtMid
            pts.append(pt)

    if not pts:
        return

    dist_per_B = []

    for rdB in rdBs_Pass:
        dists_thisB = []
        for pt in pts:
            pt_OnB = rdB.BrepGeometry.ClosestPoint(pt)
            dist = pt_OnB.DistanceTo(pt)
            dists_thisB.append(dist)
        dist_per_B.append(min(dists_thisB))

    return rdBs_Pass[dist_per_B.index(min(dist_per_B))]


def main():

    objrefs_Ref = getInput_Ref()
    if objrefs_Ref is None: return

    bAutoPickSrf = Opts.values['bAutoPickSrf']
    iClosestGripCtPerAnalysisPt = Opts.values['iClosestGripCtPerAnalysisPt']
    bAddInterior = Opts.values['bAddInterior']
    fTolForSharedGrips = Opts.values['fTolForSharedGrips']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if bAutoPickSrf:
        rdB = _findClosestSrf(objrefs_Ref)
        #rdB.Select(
        #    on=True,
        #    syncHighlight=True,
        #    persistentSelect=True,
        #    ignoreGripsState=True,
        #    ignoreLayerLocking=False,
        #    ignoreLayerVisibility=False)

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    selectGrips(
        objrefs_Ref=objrefs_Ref,
        rdBrep_1F=rdB,
        iClosestGripCtPerAnalysisPt=iClosestGripCtPerAnalysisPt,
        bAddInterior=bAddInterior,
        fTolForSharedGrips=fTolForSharedGrips,
        bEcho=bEcho,
        bDebug=bDebug,
    )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()