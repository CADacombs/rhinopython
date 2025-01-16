"""
This script is an alternative to _Drape. It uses Greville point locations for fitting to
the target object and allows selection of a starting surface.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
191018-19: Created.
...
210718: Starting surface can no longer be selected from one of the breps to fit.
210730-0808: Reduced wavy output along elevation transitions of target.
210808: Now, active CPlane's Z axis is used for drape direction.
250115: Added support for meshes as objects over which to drape. Refactored.

Starting surface's control points' X and Y are maintained.
Starting surface's Greville points are used for measurement.

WIP:
    Fix bug for highest elevation in ...HighToLow function.
    Add option to 'ledge' neighbors more than 1 point beyond target.
    Add support for mesh targets.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

#import itertools
import math
import random


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'bFlipCPlane'; keys.append(key)
    values[key] = False
    names[key] = 'DrapeDir'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'NegCPlaneZAxis', 'PosCPlaneZAxis')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUserProvidesStartingSrf'; keys.append(key)
    values[key] = False
    names[key] = 'StartingSrf'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Create', 'UserProvides')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fPointSpacing'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iExtrapolationCt'; keys.append(key)
    values[key] = 1
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=0)
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

        if key == 'fPointSpacing':
            if cls.riOpts[key].CurrentValue <= 10.0*sc.doc.ModelAbsoluteTolerance:
                print("Invalid input for tolerance.")
                cls.riOpts[key].CurrentValue = cls.values[key]

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue <= 0.0:
                print("Invalid input for tolerance.")
                cls.riOpts[key].CurrentValue = cls.values[key]

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_ObjsToDrapeOver():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select target breps or meshes")
    go.GeometryFilter = rd.ObjectType.Brep | rd.ObjectType.Mesh

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bFlipCPlane')
        addOption('bUserProvidesStartingSrf')
        if not Opts.values['bUserProvidesStartingSrf']:
            addOption('fPointSpacing')
        addOption('fTolerance')
        addOption('iExtrapolationCt')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            cPlane = go.View().ActiveViewport.ConstructionPlane()
            if Opts.values['bFlipCPlane']:
                cPlane.Flip()
            go.Dispose()
            rdBreps = [rs.coercerhinoobject(o) for o in objrefs]
            return (
                rdBreps,
                cPlane,
                Opts.values['bUserProvidesStartingSrf'],
                Opts.values['fPointSpacing'],
                Opts.values['fTolerance'],
                Opts.values['iExtrapolationCt'],
                Opts.values['bEcho'],
                Opts.values['bDebug']
                )

        if not Opts.values['bUserProvidesStartingSrf']:
            if res == ri.GetResult.Number:
                key = 'fPointSpacing'
                Opts.riOpts[key].CurrentValue = go.Number()
                Opts.setValue(key)
                continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_StartingSurface(rdObjs_toFit):
    """
    Get Surface with optional input.
    """

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select starting surface")
    go.GeometryFilter = rd.ObjectType.Surface

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fTolerance')
        addOption('iExtrapolationCt')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            if objref.ObjectId in [o.Id for o in rdObjs_toFit]:
                print("Starting surface cannot be one of the objects to fit")
                sc.doc.Objects.UnselectAll()
                go = ri.Custom.GetObject()
                go.SetCommandPrompt("Select starting surface")
                go.GeometryFilter = rd.ObjectType.Surface
                continue

            # Check that surface is Degree x and contains only simple knots.
            #rgF = objref.Surface()
            #rgS = rgF.UnderlyingSurface()

            return (
                objref,
                Opts.values['fTolerance'],
                Opts.values['iExtrapolationCt'],
                Opts.values['bEcho'],
                Opts.values['bDebug']
                )

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _prompt(sPrompt, bDebug=False):
    if bDebug: print(sPrompt)
    Rhino.RhinoApp.SetCommandPrompt(sPrompt)
    Rhino.RhinoApp.Wait()


def _promptDone(bAddWorking=True, bDebug=False):
    if bDebug: print(Rhino.RhinoApp.CommandPrompt + " done.")
    if bAddWorking:
        Rhino.RhinoApp.SetCommandPrompt("Working ...")
        Rhino.RhinoApp.Wait()


def _addKnotsToSurface(ns_In):
    """
    Returns: rg.NurbsSurface
    """

    ns_Out = ns_In.Duplicate()

    iK = ns_In.KnotsU.Count-1-ns_In.Degree(0)
    k_M = ns_In.KnotsU[iK]
    k_L = ns_In.KnotsU[iK-1]
    ns_Out.KnotsU.InsertKnot(k_L/3.0 + 2.0*k_M/3.0)

    for iK in range(ns_In.KnotsU.Count-2-ns_In.Degree(0), ns_In.Degree(0), -1):
        k_M = ns_In.KnotsU[iK]
        k_L = ns_In.KnotsU[iK-1]
        k_R = ns_In.KnotsU[iK+1]
        ns_Out.KnotsU.InsertKnot(k_R/3.0 + 2.0*k_M/3.0)
        ns_Out.KnotsU.InsertKnot(k_L/3.0 + 2.0*k_M/3.0)

    iK = ns_In.Degree(0)
    k_M = ns_In.KnotsU[iK]
    k_R = ns_In.KnotsU[iK+1]
    ns_Out.KnotsU.InsertKnot(k_R/3.0 + 2.0*k_M/3.0)

    iK = ns_In.KnotsV.Count-1-ns_In.Degree(1)
    k_M = ns_In.KnotsV[iK]
    k_L = ns_In.KnotsV[iK-1]
    ns_Out.KnotsV.InsertKnot(k_L/3.0 + 2.0*k_M/3.0)

    for iK in range(ns_In.KnotsV.Count-2-ns_In.Degree(1), ns_In.Degree(1), -1):
        k_M = ns_In.KnotsV[iK]
        k_L = ns_In.KnotsV[iK-1]
        k_R = ns_In.KnotsV[iK+1]
        ns_Out.KnotsV.InsertKnot(k_R/3.0 + 2.0*k_M/3.0)
        ns_Out.KnotsV.InsertKnot(k_L/3.0 + 2.0*k_M/3.0)

    iK = ns_In.Degree(1)
    k_M = ns_In.KnotsV[iK]
    k_R = ns_In.KnotsV[iK+1]
    ns_Out.KnotsV.InsertKnot(k_R/3.0 + 2.0*k_M/3.0)

    return ns_Out


def _getGrevillePoints(ns):
    pts_out = []

    for iU in range(ns.Points.CountU):
        pts_out.append([])
        for iV in range(ns.Points.CountV):
            u, vN = ns.Points.GetGrevillePoint(iU, iV)
            pt = ns.PointAt(u, vN)
            pts_out[-1].append(pt)

    return pts_out


def _createPointsProjectedToTargets(rgObjs_Targets, pts_toProject, bDebug=False):
    """
    rhObjects_Ref can include ObjRefs, DocObjects.RhinoObjects, GUIDS, or Geometry, but must be all the same type.
    Returns:
        list(rg.Point3ds)
    """

    _prompt("Creating target points ...", bDebug=bDebug)

    rgBreps_Targets = []
    rgMeshes_Targets = []
    for rgObj in rgObjs_Targets:
        if isinstance(rgObj, rg.Brep):
            rgBreps_Targets.append(rgObj)
        elif isinstance(rgObj, rg.Mesh):
            rgMeshes_Targets.append(rgObj)
        else:
            raise Exception()

    pts_Out = []
    
    for iU in range(len(pts_toProject)):
        pts_Out.append([])
        for iV in range(len(pts_toProject[0])):
            if sc.escape_test(throw_exception=False):
                print("User break.")
                return

            pts_Projected = []

            if rgMeshes_Targets:
                rv = rg.Intersect.Intersection.ProjectPointsToMeshes(
                    meshes=rgMeshes_Targets,
                    points=[pts_toProject[iU][iV]],
                    direction=rg.Vector3d.ZAxis,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if rv:
                    pts_Projected.extend(rv)

            if rgBreps_Targets:
                rv = rg.Intersect.Intersection.ProjectPointsToBreps(
                    breps=rgBreps_Targets,
                    points=[pts_toProject[iU][iV]],
                    direction=rg.Vector3d.ZAxis,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if rv:
                    pts_Projected.extend(rv)

            if len(pts_Projected) == 0:
                pt_to = None
            elif len(pts_Projected) == 1:
                pt_to = pts_Projected[0]
            else:
                pts = pts_Projected
                zs = []
                for pt in pts:
                    zs.append(pt.Z)
                    dist = pt.DistanceTo(pts_toProject[iU][iV])
                winning_Z = max(zs)
                pt_to = pts[zs.index(winning_Z)]

            pts_Out[-1].append(pt_to)

    _promptDone(bDebug=bDebug)

    return pts_Out


def _get_orthogonal_neighbor_count_of_point(pts, iU, iV):
    """
    Parameters:
        pts: u list of v lists of pts
        iU: int Index in u (top-level list)
        iV: int Index of v (nested list)
    Returns:
        int (0, 1, 2, 3, or 4)
    """

    idx_MaxU = len(pts)-1
    idx_MaxV = len(pts[0])-1

    ct = 0

    if iU-1 >= 0 and pts[iU-1][iV] is not None:
        ct += 1
    if iU+1 <= idx_MaxU and pts[iU+1][iV] is not None:
        ct += 1
    if iV-1 >= 0 and pts[iU][iV-1] is not None:
        ct += 1
    if iV+1 <= idx_MaxV and pts[iU][iV+1] is not None:
        ct += 1

    return ct


def _getBorderPointIndices(pts_Target):
    idxs_borderPts = []
    for iU in range(len(pts_Target)):
        for iV in range(len(pts_Target[0])):
            if (
                pts_Target[iU][iV] is None and
                _get_orthogonal_neighbor_count_of_point(pts_Target, iU, iV) > 0
            ):
                idxs_borderPts.append((iU,iV))
    return idxs_borderPts


    def _createZAxisLinesAtGrevilles(pts_Greville_toProject):
        lines = []
        for iU in range(len(pts_Greville_toProject)):
            lines.append([])
            for iV in range(len(pts_Greville_toProject[0])):
                line = rg.Line(
                        start=pts_Greville_toProject[iU][iV],
                        span=rg.Vector3d.ZAxis)
                #sc.doc.Objects.AddLine(line)
                lines[-1].append(line)
        return lines


def _hasMissingPoints(pts):
    return any([pt is None for ptsV in pts for pt in ptsV])


def _flattenNestedList(pts_in):
    """Returns: Flattened list of 1-level nested lists"""
    pts_out = []
    for iU in range(len(pts_in)):
        for iV in range(len(pts_in[iU])):
            pts_out.append(pts_in[iU][iV])
    return pts_out


def _sortAndGroupTargetsByElevation(ns, pts_In, fElevTol, bDebug=False):
    """
    Returns:
        list(lists of u and vN index tuples whose z's are within fElevTol)
    """

    _prompt("Sorting target points by elevation ...", bDebug=bDebug)

    pts_Flat = _flattenNestedList(pts_In)
    zs_Flat_FlatOrder = [pt.Z for pt in pts_Flat]
    zs_Flat_HiToLo = sorted(zs_Flat_FlatOrder, reverse=True)
    iZs_Used = []
    iZs_Sorted_Grouped_Flat = []
    z_LastTolStart = None
    for i, z in enumerate(zs_Flat_HiToLo):
        if zs_Flat_FlatOrder.index(z) in iZs_Used:
            continue
        iZs_NewElevation = [j for j, x in enumerate(zs_Flat_FlatOrder) if x == z]
        if z_LastTolStart is None or abs(z_LastTolStart-z) > fElevTol:
            iZs_Sorted_Grouped_Flat.append(iZs_NewElevation)
            z_LastTolStart = z
        else:
            iZs_Sorted_Grouped_Flat[-1].extend(iZs_NewElevation)
        iZs_Used.extend(iZs_NewElevation)

    iUiVs_Sorted = []
    for z_Group in iZs_Sorted_Grouped_Flat:
        iUiVs_Sorted.append([])
        for iZ in z_Group:
            iUiVs_Sorted[-1].append(
                (iZ // ns.Points.CountV, iZ % ns.Points.CountV))
            #sc.doc.Objects.AddPoint(pts_In[iU][iV])

    _promptDone(bDebug=False)

    return iUiVs_Sorted


def _getNeighborsPerElevationGroup(ns, uvs_Target_Groups, zs_Targets, zs_Min_AdjustedPerNeighbors, bDebug=False):
    """
    """

    _prompt("Determining neighbors of targets ...", bDebug=bDebug)

    iUiVs_Sorted_Flat = _flattenNestedList(uvs_Target_Groups)
    uvs_WithTarget = []
    iDirs = -1, 0, 1
    # 8 directions from each index location.
    neighbor_dir_deltas = [[iU_N, iV_N] for iU_N in iDirs for iV_N in iDirs]
    uvs_Neighbors_PerElevGroup = []

    for iGroup, uvs_Target_Group in enumerate(uvs_Target_Groups):
        uvs_Neighbors_PerElevGroup.append([])
        for uT, vT in uvs_Target_Group:
            #print("T:", uT, vT)
            for uD, vD in neighbor_dir_deltas:
                uN = uT + uD
                vN = vT + vD
                #print("N:", uN, vN)
                if (uN, vN) in uvs_Target_Group:
                    # Neighbor cannot be in current target group.
                    continue
                if (uN, vN) in uvs_WithTarget:
                    continue
                if not (2 < uN < (ns.Points.CountU-3)):
                    continue
                if not (2 < vN < (ns.Points.CountV-3)):
                    continue
                uvs_Neighbors_PerElevGroup[-1].append((uN, vN))
                uvs_WithTarget.append((uN, vN))
                if zs_Min_AdjustedPerNeighbors[uN][vN] < zs_Targets[uT][vT]:
                    zs_Min_AdjustedPerNeighbors[uN][vN] = zs_Targets[uT][vT]

    _promptDone(bDebug=bDebug)

    return uvs_Neighbors_PerElevGroup


def _iterativelyFit_TranslatePointsIndividually_HighToLow(pts_Target, ns_In, fTolerance, bDebug=False):
    """
    """

    ns_Out = ns_In.Duplicate()


    #map(sc.doc.Objects.AddPoint, [pt for pts_V in pts_Target for pt in pts_V]); sc.doc.Views.Redraw(); return


    uvs_Targets_InElevGroups = _sortAndGroupTargetsByElevation(
        ns_In,
        pts_Target,
        10.0*fTolerance,
        bDebug=bDebug)

    if bDebug:
        print('-'*20)
        sEval = "uvs_Targets_InElevGroups[8][22:27]"; print(sEval,'=',eval(sEval))
        sEval = "uvs_Targets_InElevGroups[9][22:29]"; print(sEval,'=',eval(sEval))
        sEval = "uvs_Targets_InElevGroups[9][22:27]"; print(sEval,'=',eval(sEval))


    if bDebug: sEval = 'len(uvs_Targets_InElevGroups)'; print(sEval+':',eval(sEval))

    zs_HighestPerElevGroup = []
    for uvsGroup in uvs_Targets_InElevGroups:
        zs = []
        for u,v in uvsGroup:
            zs.append(pts_Target[u][v].Z)
        zs_HighestPerElevGroup.append(max(zs))

    if bDebug:
        print('-'*20)
        sEval = "zs_HighestPerElevGroup"; print(sEval,'=',eval(sEval))
        sEval = "len(zs_HighestPerElevGroup)"; print(sEval,'=',eval(sEval))

        print('-'*20)
        sEval = "uvs_Targets_InElevGroups[8][22:27]"; print(sEval,'=',eval(sEval))
        sEval = "uvs_Targets_InElevGroups[9][22:29]"; print(sEval,'=',eval(sEval))
        sEval = "uvs_Targets_InElevGroups[9][22:27]"; print(sEval,'=',eval(sEval))

    #for pts_V in pts_Target:
    #    for pt in pts_V:
    #        sc.doc.Objects.AddPoint(pt)
    #    sc.doc.Views.Redraw()
    #return


    #attr = rd.ObjectAttributes()
    #attr.ColorSource = rd.ObjectColorSource.ColorFromObject
    #for iUiV_Group in uvs_Targets_InElevGroups:
    #    attr.ObjectColor = Color.FromArgb(
    #            red=random.randint(0, 255),
    #            green=random.randint(0, 255),
    #            blue=random.randint(0, 255))
    #    for iU, iV in iUiV_Group:
    #        sc.doc.Objects.AddPoint(pts_Target[iU][iV], attr)
    #sc.doc.Views.Redraw(); 1/0


    zs_Targets = [[pt.Z for pt in ptsV] for ptsV in pts_Target]

    if bDebug:
        print('-'*20)
        sEval = "zs_Targets[8][22:27]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Targets[9][22:29]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Targets[9][22:27]"; print(sEval,'=',eval(sEval))


    # Set zs_Targets in elevation groups to highest elevation.
    # This allows better control point selection in resultant surface.
    for iGroup, uvsGroup in enumerate(uvs_Targets_InElevGroups):
        for u, v in uvsGroup:
            zs_Targets[u][v] = zs_HighestPerElevGroup[iGroup]

    if bDebug:
        print('-'*20)
        sEval = "zs_Targets[8][22:27]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Targets[9][22:29]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Targets[9][22:27]"; print(sEval,'=',eval(sEval))

    zs_Min_AdjustedPerNeighbors = [zsV[:] for zsV in zs_Targets]

    if bDebug:
        print('-'*20)
        sEval = "zs_Min_AdjustedPerNeighbors[8][22:27]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Min_AdjustedPerNeighbors[9][22:29]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Min_AdjustedPerNeighbors[9][22:27]"; print(sEval,'=',eval(sEval))


    uvs_Neighbors_PerElevGroup = _getNeighborsPerElevationGroup(
        ns_In,
        uvs_Targets_InElevGroups,
        zs_Targets,
        zs_Min_AdjustedPerNeighbors,
        bDebug=bDebug)

    if bDebug:
        print('-'*20)
        sEval = "zs_Min_AdjustedPerNeighbors[8][22:27]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Min_AdjustedPerNeighbors[9][22:29]"; print(sEval,'=',eval(sEval))
        sEval = "zs_Min_AdjustedPerNeighbors[9][22:27]"; print(sEval,'=',eval(sEval))


    #print(uvs_Targets_InElevGroups[:10])
    if bDebug: sEval = 'len(uvs_Neighbors_PerElevGroup)'; print(sEval+':',eval(sEval))


    # Do not allow first neighber group to be lowered.
    for (uN, vN) in uvs_Neighbors_PerElevGroup[0]:
        cp = ns_In.Points.GetControlPoint(uN, vN)
        zs_Min_AdjustedPerNeighbors[uN][vN] = cp.Location.Z


    # Start loop at level 1, not 0, because
    # top elevations and neighbors stay at highest elevation.

    uvs_Neighbors_Cum_prev_Flat = []

    for iGroup in range(1, len(uvs_Targets_InElevGroups)):
        if iGroup > 1:
            _promptDone(False, bDebug=False)

        _prompt("Fitting elevation level {} of {} ...".format(
            iGroup+1, len(uvs_Targets_InElevGroups)), bDebug=False)
            
        uvs_Target_Group = uvs_Targets_InElevGroups[iGroup]
        uvs_NeighborsOfGroup = uvs_Neighbors_PerElevGroup[iGroup]
        uvs_Neighbors_Cum_prev_Flat.extend(uvs_Neighbors_PerElevGroup[iGroup-1])

        # Set CP locations of targets only if not already set as a neighbor.
        for uT,vT in uvs_Target_Group:
            if (uT, vT) in uvs_Neighbors_Cum_prev_Flat:
                continue
            cp = ns_Out.Points.GetControlPoint(uT, vT)
            cp.Z = pts_Target[uT][vT].Z
            ns_Out.Points.SetControlPoint(uT, vT, cp)


        ## Translate neighbors as low as possible to their target Z.

        # First, translate neighbors to their ? Z.
        for uN,vN in uvs_NeighborsOfGroup:

            # Do not modify first group.
            if (uN, vN) in uvs_Neighbors_PerElevGroup[0]:
                #print("Skipped {} {}".format(uN, vN))
                continue

            cp = ns_Out.Points.GetControlPoint(uN, vN)
            #z_Low = pts_Target[uN][vN].Z
            z_Low = zs_Min_AdjustedPerNeighbors[uN][vN]
            cp.Z = z_Low
            ns_Out.Points.SetControlPoint(uN, vN, cp)

            if (uN,vN) == (6,5):
                pass

            #    for uT,vT in uvs_Target_Group:
            #        sc.doc.Objects.AddPoint(pts_Target[uT][vT])
            #    sc.doc.Views.Redraw(); return


        # Test whether Greville points of elevation group are still on or above target.

        for uT,vT in uvs_Target_Group:
            uG,vG = ns_Out.Points.GetGrevillePoint(uT, vT)
            #print(uT, vT)
            if (uT,vT) == (6,5):
                pass
            #    sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); return
            zG = ns_Out.PointAt(uG, vG).Z
            zT = pts_Target[uT][vT].Z
            if (zG + 0.001*fTolerance) >= zT:
                continue
            #sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); 1/0
            break
        else:
            #if bDebug:
            #    _prompt("All Grevilles are on or above target.", bDebug=False)
            # No change to zs_Min_AdjustedPerNeighbors.
            #sc.doc.Objects.AddPoint(pts_Target[u][vN])
            continue # to next elevation group.
        #sc.doc.Views.Redraw(); 1/0

        ##


        #if bDebug:
        #    _prompt("Binary search the correct elevation.", bDebug=False)
        fraction_L = 0.0
        fraction_H = 1.0

        while True:
            sc.escape_test()

            #if bDebug: print("L,H 'fraction': {}, {}".format(fraction_L, fraction_H))

            fraction_M = 0.5*fraction_L + 0.5*fraction_H

            for uN,vN in uvs_NeighborsOfGroup:
                # Translate point as low as possible to its target Z.
                cp = ns_Out.Points.GetControlPoint(uN,vN)
                z_Lowest = zs_Min_AdjustedPerNeighbors[uN][vN]

                # Instead of starting surface, use elevation of closest
                # neighbor in uvs_Targets_InElevGroups.
                #def getElevationOfClosestTarget():
                #    dists = []
                #    for (uT, vT) in uvs_Targets_InElevGroups[iGroup]:
                #        dist = ((float(uN - uT))**2 + (float(vN - vT))**2)**0.5
                #        dists.append(dist)
                #        dist_Min = min(dists)
                #    zs_Winners = []
                #    for i, dist in enumerate(dists):
                #        if abs(dist_Min-dist) <= 1e-9:
                #            zs_Winners.append(zs_PerElevGroup[iGroup][i])
                #    return sum(zs_Winners) / float(len(zs_Winners))

                # Instead of starting surface, use highest elevation of
                # neighbors.
                def getHighestElevationOfNeighbors():
                    iDirs = -1, 0, 1
                    # 8 directions from each index location.
                    delta_dirs = [[u, v] for u in iDirs for v in iDirs]

                    zs_Neighbors = []

                    for uD, vD in delta_dirs:
                        uNN = uN + uD
                        vNN = vN + vD
                        zs_Neighbors.append(zs_Targets[uNN][vNN])

                    return max(zs_Neighbors)

                z_Highest = getHighestElevationOfNeighbors()


                cp.Z = z_Lowest + (z_Highest-z_Lowest)*fraction_M
                ns_Out.Points.SetControlPoint(uN,vN,cp)

            for uT,vT in uvs_Target_Group:
                uG,vG = ns_Out.Points.GetGrevillePoint(uT, vT)
                zG = ns_Out.PointAt(uG, vG).Z
                zT = pts_Target[uT][vT].Z
                if zG >= zT:
                    continue
                # Greville is too low.
                fraction_L = fraction_M
                break
            else:
                # All Grevilles are on or above target.
                fraction_H = fraction_M

            if abs(fraction_H - fraction_L) <= 0.001:
                #sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); 1/0
                break # out of while / binary search.

        for uN,vN in uvs_NeighborsOfGroup:
            # Translate point as low as possible to its target Z.
            cp = ns_Out.Points.GetControlPoint(uN,vN)
            zs_Min_AdjustedPerNeighbors[uN][vN] = cp.Location.Z

        #sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); 1/0


    if bDebug:
        print("Iteratived through {} elevation groups.".format(iGroup+1))
        print("Position points not within 3 from border nor are already translated.")

    uvs_done_Flat = uvs_Targets_InElevGroups[0] + [(uN, vN) for uvs in uvs_Neighbors_PerElevGroup for (uN, vN) in uvs]

    for uN in range(3, ns_In.Points.CountU-3):
        for vN in range(3, ns_In.Points.CountV-3):
            if (uN, vN) not in uvs_done_Flat:
                cp = ns_Out.Points.GetControlPoint(uN,vN)
                cp.Z = pts_Target[uN][vN].Z
                ns_Out.Points.SetControlPoint(uN,vN,cp)


    return ns_Out


def _iterativelyFit_TranslatePointsIndividually(pts, ns_In, fTolerance):
    """
    """

    ns_Out = ns_In.Duplicate()

    # Initially, move control points whose Grevilles are not within tolerance of their targets.
    for iU in range(ns_In.Points.CountU):
        for iV in range(ns_In.Points.CountV):
            uv_Gr = ns_Out.Points.GetGrevillePoint(iU, iV)
            pt_Gr = ns_Out.PointAt(uv_Gr[0], uv_Gr[1])
            dist = pt_Gr.DistanceTo(pts[iU][iV])
            vect = pts[iU][iV] - pt_Gr
            if vect.Length <= fTolerance:
                continue

            ns_Out.Points.SetControlPoint(iU,iV,pts[iU][iV])

    for i in xrange(200):
        bTransPts = False
        for iU in range(ns_In.Points.CountU):
            for iV in range(ns_In.Points.CountV):
                uv_Gr = ns_Out.Points.GetGrevillePoint(iU, iV)
                pt_Gr = ns_Out.PointAt(uv_Gr[0], uv_Gr[1])
                dist = pt_Gr.DistanceTo(pts[iU][iV])
                vect = pts[iU][iV] - pt_Gr
                if vect.Length <= fTolerance:
                    continue
                #print(dist, vect)
                cp = ns_Out.Points.GetControlPoint(iU,iV)
                ns_Out.Points.SetControlPoint(iU,iV,cp.Location+vect)
                bTransPts = True
        if not bTransPts:
            print("{} iterations for Grevilles to lie on target(s) within {}.".format(
                i+1, fTolerance))
            return ns_Out

    print("After {} iterations, Grevilles still do not lie on target(s) within {}.".format(
        i+1, fTolerance))

    return ns_Out


def _iterativelyFit_TranslatePoints_GroupedByDegree(pts_Target, ns_In, fTolerance):
    """
    """

    ns_Out = ns_In.Duplicate()

    # Initially, move control points whose Grevilles are not within tolerance of their targets.

    bPos_NotNeg_Z_trans_from_start = None
    iCt_Neg_Z_trans = 0
    iCt_Pos_Z_trans = 0

    for iU_Group_SW in range(ns_Out.Degree(0), ns_Out.Points.CountU-ns_Out.Degree(0), ns_Out.Degree(0)):
        if ((iU_Group_SW + ns_Out.Degree(0)) >= ns_Out.Points.CountU):
            break
        for iV_Group_SW in range(ns_Out.Degree(1), ns_Out.Points.CountV-ns_Out.Degree(1), ns_Out.Degree(1)):
            if ((iV_Group_SW + ns_Out.Degree(1)) >= ns_Out.Points.CountV):
                continue
            zs_Trans = []
            dists = []
            iUs_dists = []
            iVs_dists = []
            for iU_Sub in range(ns_Out.Degree(0)):
                for iV_Sub in range(ns_Out.Degree(1)):
                    iU = iU_Group_SW + iU_Sub
                    iV = iV_Group_SW + iV_Sub
                    cp = ns_Out.Points.GetControlPoint(iU, iV)
                    pt_cp = cp.Location
                    #sc.doc.Objects.AddPoint(pt_cp)
                    #sc.doc.Objects.AddPoint(pts_Target[iU][iV])
                    vect = pts_Target[iU][iV] - pt_cp
                    zs_Trans.append(vect.Z)
                    dist = abs(vect.Z)
                    dists.append(dist)
                    iUs_dists.append(iU)
                    iVs_dists.append(iV)
                    if dist > fTolerance:
                        if bPos_NotNeg_Z_trans_from_start is None:
                            bPos_NotNeg_Z_trans_from_start = vect.Z > 0
                        elif bPos_NotNeg_Z_trans_from_start != (vect.Z > 0):
                            print("Starting surface is not completely on one side of targets.")
                            return

            goal_dist = max(dists)

            #print(goal_dist)

            if goal_dist <= fTolerance:
                continue # to next group.

            idx_goal = dists.index(goal_dist)

            z = pts_Target[iUs_dists[idx_goal]][iVs_dists[idx_goal]].Z

            for iU_Sub in range(ns_Out.Degree(0)):
                for iV_Sub in range(ns_Out.Degree(1)):
                    iU = iU_Group_SW + iU_Sub
                    iV = iV_Group_SW + iV_Sub
                    cp_New = ns_Out.Points.GetControlPoint(iU, iV)
                    cp_New.Z = z
                    ns_Out.Points.SetControlPoint(
                        iU,
                        iV,
                        cp_New)

    #attr.ObjectColor = Color.DarkCyan
    #sc.doc.Objects.AddSurface(ns_Out, attr); sc.doc.Views.Redraw()#; return


    for i in xrange(20):
        sc.escape_test()
        print(i)
        bTransPts = False

        for iU_Group_SW in range(0, ns_Out.Points.CountU, ns_Out.Degree(0)):
            if ((iU_Group_SW + ns_Out.Degree(0)) >= ns_Out.Points.CountU):
                break
            for iV_Group_SW in range(0, ns_Out.Points.CountV, ns_Out.Degree(1)):
                if ((iV_Group_SW + ns_Out.Degree(1)) >= ns_Out.Points.CountV):
                    continue
                dists_FromStart = []
                zs_Trans = []
                iUs_dists = []
                iVs_dists = []
                for iU_Sub in range(ns_Out.Degree(0)):
                    for iV_Sub in range(ns_Out.Degree(1)):
                        iU = iU_Group_SW + iU_Sub
                        iV = iV_Group_SW + iV_Sub
                        u_Gr, v_Gr = ns_Out.Points.GetGrevillePoint(iU, iV)
                        pt_Gr = ns_Out.PointAt(u_Gr, v_Gr)
                        vect = pts_Target[iU][iV] - pt_Gr
                        zs_Trans.append(vect.Z)

                z_Trans = min(zs_Trans) if bPos_NotNeg_Z_trans_from_start else max(zs_Trans)
                if abs(z_Trans) <= fTolerance:
                    continue # to next group.

                for iU_Sub in range(ns_Out.Degree(0)):
                    for iV_Sub in range(ns_Out.Degree(1)):
                        iU = iU_Group_SW + iU_Sub
                        iV = iV_Group_SW + iV_Sub
                        cp_New = ns_Out.Points.GetControlPoint(iU, iV)
                        cp_New.Z += z_Trans
                        ns_Out.Points.SetControlPoint(
                            iU,
                            iV,
                            cp_New)

                bTransPts = True

        #sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); #return

        if not bTransPts:
            print("{} iterations for Grevilles to lie on target(s) within {}.".format(
                i+1, fTolerance))
            return ns_Out

    print("After {} iterations, Grevilles still do not lie on target(s) within {}.".format(
        i+1, fTolerance))

    return ns_Out


def _addObject(rgObj, xform):
    if xform:
        if isinstance(rgObj, rg.GeometryBase):
            rgObj_Dup = rgObj.Duplicate()
            rgObj_Dup.Transform(xform)
            sc.doc.Objects.Add(rgObj)
        elif isinstance(rgObj, rg.Point3d):
            rgPt3d_Dup = rg.Point3d(rgObj)
            rgPt3d_Dup.Transform(xform)
            sc.doc.Objects.AddPoint(rgPt3d_Dup)
        else:
            raise Exception("{} not supported yet.".format(rgObj.GetType().Name))
    else:
        if isinstance(rgObj, rg.GeometryBase):
            sc.doc.Objects.Add(rgObj)
        elif isinstance(rgObj, rg.Point3d):
            sc.doc.Objects.AddPoint(rgPt3d_Dup)


def processGeometry(rgObjs_toDrapeOver, srf_Starting, cPlane=rg.Plane.WorldXY, fTolerance=None, iExtrapolationCt=1, bDebug=False):
    """
    """

    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance

    ns_Starting = srf_Starting.ToNurbsSurface()

    #sEval = "cPlane"; print(sEval,'=',eval(sEval))
    #return

    if cPlane == rg.Plane.WorldXY:
        xform = None
    else:
        xform = rg.Transform.PlaneToPlane(cPlane, rg.Plane.WorldXY)
        [rgObjs_toDrapeOver[i].Transform(xform) for i in range(len(rgObjs_toDrapeOver))]
        ns_Starting.Transform(xform)

        # Prepare xform for output.
        xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)
        #[sc.doc.Objects.AddBrep(rgBreps_ProjectTo[i]) for i in range(len(rgBreps_ProjectTo))]
        #sc.doc.Objects.AddSurface(ns_Starting); sc.doc.Views.Redraw()#; return


    uvs_All = [(u, v)
               for u in range(ns_Starting.Points.CountU)
               for v in range(ns_Starting.Points.CountV)]


    #ns_WIP = addKnotsToSurface(ns_Starting)


    pts_Greville = _getGrevillePoints(ns_Starting)

    #[sc.doc.Objects.AddPoint(pts_Greville[u][v]) for u, v in uvs_All if pts_Greville[u][v] is not None]
    #sc.doc.Views.Redraw(); return

    pts_Target = _createPointsProjectedToTargets(
        rgObjs_Targets=rgObjs_toDrapeOver,
        pts_toProject=pts_Greville,
        bDebug=bDebug,
        )
    if not pts_Target:
        print("Projected points were not obtained.")
        return


    #[_addObject(pts_Target[u][v], xform) for u, v in uvs_All if pts_Target[u][v] is not None]
    #sc.doc.Views.Redraw(); return

    #for iU in range(len(pts_Target)):
    #    for iV in range(len(pts_Target[iU])):
    #        #print(iU, iV, pts_Target[iU][iV]))
    #        if pts_Target[iU][iV] is None:
    #            continue
    #        sc.doc.Objects.AddPoint(pts_Target[iU][iV])
            #line = rg.Line(start=pts_Target[iU][iV], span=norms_Projected[iU][iV])
            #sc.doc.Objects.AddLine(line)
    #sc.doc.Views.Redraw(); return


    def closestPointsOfNeighborsOnNormalLines(pts_In, iU, iV, iMinNeighborCt=1, bDiag=True, bLineExts=True):
        """
        bLineExts:
            When True: If neighbor's point and its neighbor's point in the same
            direction are both available, create a line through those points
            and get the ClosestPoint of that line on the normal line.
        """
        pts_Out = []

        idx_MaxU = len(pts_In)-1
        idx_MaxV = len(pts_In[0])-1

        #attr.ObjectColor = Color.FromArgb(
        #        red=random.randint(0, 255),
        #        green=random.randint(0, 255),
        #        blue=random.randint(0, 255))

        # West (Previous U).
        if iU-1 >= 0 and pts_In[iU-1][iV] is not None:
            pt = None
            if bLineExts and iU-2 >= 0 and pts_In[iU-2][iV] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU-2][iV], pts_In[iU-1][iV])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU-1][iV],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # East (Next U).
        if iU+1 <= idx_MaxU and pts_In[iU+1][iV] is not None:
            pt = None
            if bLineExts and iU+2 <= idx_MaxU and pts_In[iU+2][iV] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU+2][iV], pts_In[iU+1][iV])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU+1][iV],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # South (Previous V).
        if iV-1 >= 0 and pts_In[iU][iV-1] is not None:
            pt = None
            if bLineExts and iV-2 >= 0 and pts_In[iU][iV-2] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU][iV-2], pts_In[iU][iV-1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU][iV-1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # North (Next V).
        if iV+1 <= idx_MaxV and pts_In[iU][iV+1] is not None:
            pt = None
            if bLineExts and iV+2 <= idx_MaxV and pts_In[iU][iV+2] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU][iV+2], pts_In[iU][iV+1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU][iV+1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Southwest (Previous U, Previous V).
        if (
                bDiag and
                iU-1 >= 0 and
                iV-1 >= 0 and
                pts_In[iU-1][iV-1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU-2 >= 0 and
                    iV-2 >= 0 and
                    pts_In[iU-2][iV-2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU-2][iV-2], pts_In[iU-1][iV-1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU-1][iV-1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Southeast (Next U, Previous V).
        if (
                bDiag and
                iU+1 <= idx_MaxU and
                iV-1 >= 0 and
                pts_In[iU+1][iV-1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU+2 <= idx_MaxU and
                    iV-2 >= 0 and
                    pts_In[iU+2][iV-2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU+2][iV-2], pts_In[iU+1][iV-1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU+1][iV-1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Northwest (Previous U, Next V).
        if (
                bDiag and
                iU-1 >= 0 and
                iV+1 <= idx_MaxV and
                pts_In[iU-1][iV+1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU-2 >= 0 and
                    iV+2 <= idx_MaxV and
                    pts_In[iU-2][iV+2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU-2][iV+2], pts_In[iU-1][iV+1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU-1][iV+1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Northeast (Next U, Next V).
        if (
                bDiag and
                iU+1 <= idx_MaxU and
                iV+1 <= idx_MaxV and
                pts_In[iU+1][iV+1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU+2 <= idx_MaxU and
                    iV+2 <= idx_MaxV and
                    pts_In[iU+2][iV+2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU+2][iV+2], pts_In[iU+1][iV+1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU+1][iV+1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        if len(pts_Out) < iMinNeighborCt:
            return []

        return pts_Out


    def addMissingPointsAlongBorder(pts, idxs_pt_filter=None, iMinNeighborCt=1, bDiag=False, bLineExts=True):
        """
        Modify a copy of the list so that new points do not affect
        subsequent ones in this function call.
        """

        pts0 = pts

        pts_Out = [ptsV[:] for ptsV in pts0]

        bModificationOccured = False

        for iU in range(len(pts0)):
            for iV in range(len(pts0[0])):
                if idxs_pt_filter and not (iU, iV) in idxs_pt_filter: continue

                if pts0[iU][iV] is not None: continue

                pts = closestPointsOfNeighborsOnNormalLines(
                        pts0,
                        iU,
                        iV,
                        iMinNeighborCt=iMinNeighborCt,
                        bDiag=bDiag,
                        bLineExts=bLineExts)
                if not pts: continue

                pt_Sum = None
                for pt in pts:
                    pt_Sum = pt if pt_Sum is None else pt_Sum + pt

                pt = pt_Sum / float(len(pts))
                #sc.doc.Objects.AddPoint(pt, attr)

                pts_Out[iU][iV] = pt

                bModificationOccured = True

        # Modify original list.
        for iU in range(len(pts0)):
            for iV in range(len(pts0[0])):
                pts0[iU][iV] = pts_Out[iU][iV]

        return bModificationOccured



    if not _hasMissingPoints(pts_Target):
        if len(rgObjs_toDrapeOver) == 1 and rgObjs_toDrapeOver[0].Faces.Count == 1:
            ns1 = _iterativelyFit_TranslatePointsIndividually(
                pts_Target,
                ns_Starting,
                fTolerance,
                bDebug=bDebug)
        else:
            ns1 = _iterativelyFit_TranslatePointsIndividually_HighToLow(
                pts_Target,
                ns_Starting,
                fTolerance,
                bDebug=bDebug)

        if cPlane != rg.Plane.WorldXY:
            xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)
            ns1.Transform(xform)

        return ns1


    pts_Target_HasMissing = [ptsV[:] for ptsV in pts_Target]


    if iExtrapolationCt < 1:

        # Fill any remaining missing points with the
        # starting surface control point locations.

        def addMissingPerStartingSrf(pts, ns_Starting):
            for iU in range(len(pts_Target)):
                for iV in range(len(pts_Target[0])):
                    if pts_Target[iU][iV] is None:
                        pts_Target[iU][iV] = ns_Starting.Points.GetControlPoint(iU,iV).Location

        addMissingPerStartingSrf(pts_Target, ns_Starting)

        ns1 = _iterativelyFit_TranslatePointsIndividually_HighToLow(
                    pts_Target,
                    ns_Starting,
                    fTolerance)

        if cPlane != rg.Plane.WorldXY:
            xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)
            ns1.Transform(xform)

        return ns1



    # Extrapolate missing points.

    idxs_borderPts = _getBorderPointIndices(pts_Target)
    #print(idxs_borderPts); 1/0

    lines_thruStartingSrfNormals = _createZAxisLinesAtGrevilles(pts_Greville)
    #for line in lines_thruStartingSrfNormals: sc.doc.Objects.AddLine(line)
    #sc.doc.Views.Redraw(); 1/0

    attr = rd.ObjectAttributes()
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject

    #for bDiag1, bDiag2, bLineExts1, bLineExts2 in itertools.product(
    #            (False, True), (False, True), (False, True), (False, True)):
    for bLineExts1, bLineExts2 in ((True, False),):
    #for bLineExts1, bLineExts2 in itertools.product(
    #            (False, True), (False, True)):

        pts_Target = [ptsV[:] for ptsV in pts_Target_HasMissing]


        #for iMinNeighborCt in 4,3,2,1: #(1,): #

        #    addMissingPointsAlongBorder(
        #            pts_Target,
        #            idxs_pt_filter=idxs_borderPts,
        #            iMinNeighborCt=iMinNeighborCt,
        #            bDiag=True,
        #            bLineExts=bLineExts1)

        #    if not hasMissingPoints(pts_Target):
        #        break


        i = 2 # For 2nd iteration (not index) to be relevant with iExtrapolationCt.

        while (
            _hasMissingPoints(pts_Target) and
            ((iExtrapolationCt == 0) or (i < iExtrapolationCt))
        ):
            sc.escape_test()

            addMissingPointsAlongBorder(
                    pts_Target,
                    bDiag=False,
                    bLineExts=bLineExts2)
            i += 1



        # Fill any remaining missing points with the
        # starting surface control point locations.

        #idxs_Pts_SameAsStartingSrf = []
        #for iU in range(len(pts_Target)):
        #    for iV in range(len(pts_Target[0])):
        #        if pts_Target[iU][iV] is None:
        #            pts_Target[iU][iV] = ns_Starting.Points.GetControlPoint(iU,iV).Location
        #            idxs_Pts_SameAsStartingSrf.append((iU,iV))


        # Fill any remaining missing points with the
        # starting surface border elevation.
        for iU in range(len(pts_Target)):
            for iV in range(len(pts_Target[0])):
                if pts_Target[iU][iV] is None:
                    pts_Target[iU][iV] = rg.Point3d(
                        ns_Starting.Points.GetControlPoint(iU,iV).X,
                        ns_Starting.Points.GetControlPoint(iU,iV).Y,
                        ns_Starting.Points.GetControlPoint(0,0).Z)


        ns_Res3 = _iterativelyFit_TranslatePointsIndividually_HighToLow(
            pts_Target,
            ns_Starting,
            fTolerance)
        if not ns_Res3: continue

        if cPlane != rg.Plane.WorldXY:
            xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)
            ns_Res3.Transform(xform)

        return ns_Res3


def _coerceSurface(rhObj):
    if isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        #print(rhObj.GeometryComponentIndex.ComponentIndexType)
        geom = rhObj.Geometry()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    srf = None
    if isinstance(geom, rg.BrepFace):
        srf = geom.UnderlyingSurface()
    elif isinstance(geom, rg.Surface):
        srf = geom
    elif isinstance(geom, rg.Brep):
        if geom.Faces.Count == 1:
            srf = geom.Faces[0].UnderlyingSurface()

    return srf


def _createStartingSurface(rgObjs_Ref, cPlane=rg.Plane.WorldXY, fPointSpacing=1.0, bDebug=False):
    """
    Returns:
        rg.NurbsSurface that is degree 3 and 2 knots beyond the target on each of the 4 sides.
    """

    bb = rg.BoundingBox.Unset

    if cPlane != rg.Plane.WorldXY:
        xform = rg.Transform.PlaneToPlane(cPlane, rg.Plane.WorldXY)
        for rgObj in rgObjs_Ref:
            rgB_Xd = rgObj.Duplicate()
            rgB_Xd.Transform(xform)
            bb.Union(rgB_Xd.GetBoundingBox(accurate=True))
            rgB_Xd.Dispose()
    else:
        for rgObj in rgObjs_Ref:
            bb.Union(rgObj.GetBoundingBox(accurate=True))

    #sc.doc.Objects.AddBox(rg.Box(bb))

    degree = 3

    # The 6.0 in the following creates 4 rows of planar control points
    # on each side, thus allowing _ExtendSrf to remain planar for the
    # degree-3 surface.

    starting_srf_X_dim = math.ceil(bb.Diagonal.X + 6.0*fPointSpacing)
    uInterval = rg.Interval(0.0, starting_srf_X_dim)
    uPointCount = int(starting_srf_X_dim / fPointSpacing) + degree

    starting_srf_Y_dim = math.ceil(bb.Diagonal.Y + 6.0*fPointSpacing)
    vInterval = rg.Interval(0.0, starting_srf_Y_dim)
    vPointCount = int(starting_srf_Y_dim / fPointSpacing) + degree

    origin = rg.Point3d(
        bb.Center.X-starting_srf_X_dim/2.0,
        bb.Center.Y-starting_srf_Y_dim/2.0,
        bb.Max.Z)

    plane = rg.Plane(origin=origin, normal=rg.Vector3d.ZAxis)

    ns = rg.NurbsSurface.CreateFromPlane(
        plane=plane,
        uInterval=uInterval,
        vInterval=vInterval,
        uDegree=degree,
        vDegree=degree,
        uPointCount=uPointCount,
        vPointCount=vPointCount)

    # Set perimeter 3 points on each side to "bottom" of target.
    for iU in range(ns.Points.CountU):
        for iV in range(ns.Points.CountV):
            if (
                (degree <= iU <= ns.Points.CountU - degree - 1) and
                (degree <= iV <= ns.Points.CountV - degree - 1)
            ):
                continue
            cp_New = ns.Points.GetControlPoint(iU, iV)
            cp_New.Z = bb.Min.Z
            ns.Points.SetControlPoint(
                iU,
                iV,
                cp_New)

    #sc.doc.Objects.AddSurface(ns); sc.doc.Views.Redraw(); 1/0

    #for u in range(ns.Points.CountU):
    #    for v in range(ns.Points.CountV):
    #        ptA = ns.Points.GetControlPoint(u, v).Location
    #        ptA.Z = bb.Max.Z
    #        ptB = rg.Point3d(ptA)
    #        ptB.Z = bb.Min.Z
    #        print(sc.doc.Objects.AddLine(ptA, ptB))
    #sc.doc.Views.Redraw(); 1/0

    if cPlane != rg.Plane.WorldXY:
        xform = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)
        ns.Transform(xform)

    return ns


def processDocObject(rdObjs_toDrapeOver, objref_srf_Starting=None, cPlane=rg.Plane.WorldXY, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fPointSpacing = getOpt('fPointSpacing')
    fTolerance = getOpt('fTolerance')
    iExtrapolationCt = getOpt('iExtrapolationCt')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    rgObjs_toDrapeOver = []
    for rhObj in rdObjs_toDrapeOver:
        rgObj_toDrapeOver = rs.coercegeometry(rhObj)
        if not isinstance(rgObj_toDrapeOver, (rg.Brep, rg.Mesh)):
            print("{} is not supported for an object to drape over.".format(rgObj_toDrapeOver.GetType().Name))
            continue
        rgObjs_toDrapeOver.append(rgObj_toDrapeOver)


    if objref_srf_Starting is None:
        ns_Starting = _createStartingSurface(
            rgObjs_toDrapeOver,
            cPlane=cPlane,
            fPointSpacing=fPointSpacing,
            bDebug=bDebug)
    else:
        srf_Starting = _coerceSurface(objref_srf_Starting)
        ns_Starting = srf_Starting.ToNurbsSurface()


    if bEcho or bDebug:
        point_count = ns_Starting.Points.CountU * ns_Starting.Points.CountV
        if point_count > 10000:
            print("Start surface has {} points.".format(point_count))

    ns_Res = processGeometry(
        rgObjs_toDrapeOver=rgObjs_toDrapeOver,
        srf_Starting=ns_Starting,
        cPlane=cPlane,
        fTolerance=fTolerance,
        iExtrapolationCt=iExtrapolationCt,
        bDebug=bDebug
        )


    for brep in rgObjs_toDrapeOver: brep.Dispose()
    if objref_srf_Starting is not None: srf_Starting.Dispose()


    if not ns_Res: return

    g_ns1 = sc.doc.Objects.AddSurface(ns_Res)

    ns_Res.Dispose()

    if g_ns1 == Guid.Empty: return

    sc.doc.Views.Redraw()
    return g_ns1


def main():

    rc = getInput_ObjsToDrapeOver()
    if rc is None: return
    (
        rdObjs_toDrapeOver,
        cPlane,
        bUserProvidesStartingSrf,
        fPointSpacing,
        fTolerance,
        iExtrapolationCt,
        bEcho,
        bDebug,
        ) = rc

    if not bUserProvidesStartingSrf:
        objref_srf_Starting = None
    else:
        rc = getInput_StartingSurface(rdObjs_toDrapeOver)
        if rc is None: return
        (
            objref_srf_Starting,
            fTolerance,
            iExtrapolationCt,
            bEcho,
            bDebug,
            ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.CommandPrompt = "Working ..."

    processDocObject(
        rdObjs_toDrapeOver=rdObjs_toDrapeOver,
        objref_srf_Starting=objref_srf_Starting,
        cPlane=cPlane,
        fPointSpacing=fPointSpacing,
        fTolerance=fTolerance,
        iExtrapolationCt=iExtrapolationCt,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()