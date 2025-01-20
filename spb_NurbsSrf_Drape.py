"""
This script is an alternative to _Drape. It uses Greville point locations for fitting to
the target object and allows selection of a starting surface.

Starting surface's control points' X and Y are maintained.
Starting surface's Greville points are used for fitting analysis.

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
250116: Fixed bug for highest elevation in ...HighToLow function. Refactored.
250117: WIP: Started a rewrite of the HighToLow routine.
250118: Fixed a bug in creating the starting surface. It was previously adding an extra knot outside of the bounding box.
        Now, starting surface must be an open, degree-3 NURBS with only multiplicity-of-1 interior knots.
250119: Replaced the 2 bool options for missed targets to a 3-choice list.
        Now, negative values are allowed for SpansBeyondEachSide.

TODO:
    Create new HighToLow routine with a slightly new approach.
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
import random


class Data:
    def __init__(self):
        self


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'fTolerance'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

    key = 'fSpanSpacing'; keys.append(key)
    if sc.doc.ModelUnitSystem.Inches:
        values[key] = 1.0
    else:
        values[key] = 25.0 * Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iSpansBeyondEachSide'; keys.append(key)
    values[key] = 3
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iTargetMisses'; keys.append(key)
    values[key] = 2
    listValues[key] = (
        'FixToStartingSrf',
        'UseLowestNeighborHits',
        'LinearlyExtrapolateFromHits',
        )
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeleteStartingSrf'; keys.append(key)
    values[key] = False
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

        if key == 'fSpanSpacing':
            if cls.riOpts[key].CurrentValue <= 10.0*sc.doc.ModelAbsoluteTolerance:
                print("Invalid input for tolerance.")
                cls.riOpts[key].CurrentValue = cls.values[key]

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
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
        addOption('fTolerance')
        addOption('bUserProvidesStartingSrf')
        if not Opts.values['bUserProvidesStartingSrf']:
            addOption('fSpanSpacing')
            addOption('iSpansBeyondEachSide')
        addOption('iTargetMisses')
        if Opts.values['bUserProvidesStartingSrf']:
            addOption('bDeleteStartingSrf')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            #rdBreps = [rs.coercerhinoobject(o) for o in objrefs]
            return objrefs

        if res == ri.GetResult.Number:
            if Opts.values['bUserProvidesStartingSrf']:
                print("Numeric input ignored.")
                continue
            key = 'fSpanSpacing'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _doesNurbsSrfContainInteriorKnotsWithMultiplictyGT1(ns):
    for iDir in 0,1:
        knots = ns.KnotsV if iDir else ns.KnotsU
        degree = ns.Degree(iDir)
        iK = 0 if ns.IsPeriodic(iDir) else degree
        while iK < (knots.Count if ns.IsPeriodic(iDir) else knots.Count - degree):
            sc.escape_test()
            if knots.KnotMultiplicity(iK) > 1:
                return True
            iK += 1

    return False


def _isStartingSrfSupported(ns):
    """
    Supported means surface must be an open, degree-3 NURBS with only multiplicity-of-1 interior knots.
    """

    if not isinstance(ns, rg.NurbsSurface):
        return False

    if ns.Degree(0) != 3:
        return False

    if ns.Degree(1) != 3:
        return False

    if ns.IsClosed(0):
        return False

    if ns.IsClosed(1):
        return False

    if _doesNurbsSrfContainInteriorKnotsWithMultiplictyGT1(ns):
        return False

    return True


def getInput_StartingSurface(gObjs_toDrapeOver):
    """
    Get Surface with optional input.
    """

    if sc.doc.Objects.SelectedObjectsExist(objectType=rd.ObjectType.AnyObject, checkSubObjects=True):
        sc.doc.Objects.UnselectAll()
        sc.doc.Views.Redraw()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select starting surface")

    go.GeometryFilter = rd.ObjectType.Surface

    while True:
        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            sc.doc.Objects.UnselectAll()
            sc.doc.Views.Redraw()

            if objref.ObjectId in gObjs_toDrapeOver:
                print("Starting surface cannot be one of the objects to fit")
                go = ri.Custom.GetObject()
                go.SetCommandPrompt("Select starting surface")
                go.GeometryFilter = rd.ObjectType.Surface
                continue

            # Check that surface is Degree 3 and contains only simple knots.
            if not _isStartingSrfSupported(objref.Surface().UnderlyingSurface()):
                print("Starting surface must be an open, degree-3 NURBS with only multiplicity-of-1 interior knots.")
                go = ri.Custom.GetObject()
                go.SetCommandPrompt("Select starting surface")
                go.GeometryFilter = rd.ObjectType.Surface
                continue

            return objref


def _prompt(sPrompt, bDebug=False):
    if bDebug: print(sPrompt)
    Rhino.RhinoApp.SetCommandPrompt(sPrompt)
    Rhino.RhinoApp.Wait()


def _promptDone(bAddWorking=True, bDebug=False):
    if bDebug: print(Rhino.RhinoApp.CommandPrompt + " done.")
    if bAddWorking:
        Rhino.RhinoApp.SetCommandPrompt("Working ...")
        Rhino.RhinoApp.Wait()


def _createStartingSurface(rgObjs_Ref, cPlane=rg.Plane.WorldXY, fSpanSpacing=1.0, iSpansBeyondEachSide=4, bDebug=False):
    """
    Returns
        rg.NurbsSurface that is degree-3 and has 2 interior knots beyond the target bounding box on each of the 4 sides.
    """

    bb = rg.BoundingBox.Unset

    if cPlane == rg.Plane.WorldXY:
        xform_toW = xform_fromW = None
        for rgObj_Ref in rgObjs_Ref:
            bb.Union(rgObj_Ref.GetBoundingBox(accurate=True))
    else:
        xform_toW = rg.Transform.PlaneToPlane(cPlane, rg.Plane.WorldXY)
        for rgObj_Ref in rgObjs_Ref:
            rgObj_Ref_Dup = rgObj_Ref.Duplicate()
            rgObj_Ref_Dup.Transform(xform_toW)
            bb.Union(rgObj_Ref_Dup.GetBoundingBox(accurate=True))
            rgObj_Ref_Dup.Dispose()
        xform_fromW = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)

    #sc.doc.Objects.AddBox(rg.Box(bb))

    degree = 3

    dim_bb_x = bb.Diagonal.X
    starting_srf_X_dim = round(dim_bb_x + 2.0*float(iSpansBeyondEachSide)*fSpanSpacing, 0)
    uInterval = rg.Interval(0.0, starting_srf_X_dim)
    uPointCount = int(starting_srf_X_dim / fSpanSpacing) + degree

    dim_bb_y = bb.Diagonal.Y
    starting_srf_Y_dim = round(dim_bb_y + 2.0*float(iSpansBeyondEachSide)*fSpanSpacing, 0)
    vInterval = rg.Interval(0.0, starting_srf_Y_dim)
    vPointCount = int(starting_srf_Y_dim / fSpanSpacing) + degree

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

    #sc.doc.Objects.AddSurface(ns); sc.doc.Views.Redraw(); 1/0

    if xform_fromW:
        ns.Transform(xform_fromW)

    return ns


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


def _getGrevillePoints(ns):
    pts_out = []

    for iU in range(ns.Points.CountU):
        pts_out.append([])
        for iV in range(ns.Points.CountV):
            u, v = ns.Points.GetGrevillePoint(iU, iV)
            pt = ns.PointAt(u, v)
            pts_out[-1].append(pt)

    return pts_out


def _projectPts_toObjs(pts_In, rgObjs_Targets, bDebug=False):
    """
    Parameters
        rgObjs_Targets: rd.ObjRef, rd.RhinoObject, GUID, or Rhino.Geometry, but all must be of the same type.
        pts_In: list of rg.Point3d
        bDebug: bool

    Returns
        list of mix of rg.Point3ds and None (for misses)
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
    
    for iU in range(len(pts_In)):
        pts_Out.append([])
        for iV in range(len(pts_In[0])):
            if sc.escape_test(throw_exception=False):
                print("User break.")
                return

            pts_Projected = []

            if rgMeshes_Targets:
                rv = rg.Intersect.Intersection.ProjectPointsToMeshes(
                    meshes=rgMeshes_Targets,
                    points=[pts_In[iU][iV]],
                    direction=rg.Vector3d.ZAxis,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if rv:
                    pts_Projected.extend(rv)

            if rgBreps_Targets:
                rv = rg.Intersect.Intersection.ProjectPointsToBreps(
                    breps=rgBreps_Targets,
                    points=[pts_In[iU][iV]],
                    direction=rg.Vector3d.ZAxis,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if rv:
                    pts_Projected.extend(rv)

            if len(pts_Projected) == 0:
                pt_Out = None
            elif len(pts_Projected) == 1:
                pt_Out = pts_Projected[0]
            else:
                pts = pts_Projected
                zs = []
                for pt in pts:
                    zs.append(pt.Z)
                    dist = pt.DistanceTo(pts_In[iU][iV])
                winning_Z = max(zs)
                pt_Out = pts[zs.index(winning_Z)]

            pts_Out[-1].append(pt_Out)

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
    """
    Returns
        list of lists of lines
    """
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


def _sort2ListsPer1st(list1, list2):
    sorted_zipped = sorted(zip(list1, list2), reverse=True)
    sorted1, sorted2 = zip(*sorted_zipped)
    return sorted1, sorted2


def _sortAndGroupTargetsByElevation_High_to_low_NEW(pts, fElevTol=1e-6, bDebug=False):
    """
    Parameters
        pts: list of lists of mix of rg.Point3d or None
        fElevTol: float  Should be tighter than general tolerance to minimize tolerance stack up.

    Returns
        tuple of
            list of lists of z values rounded per fElevTol
            list of lists of u and v index tuples whose z's are within fElevTol
    """

    sEval = "fElevTol"; print(sEval,'=',eval(sEval))

    _prompt("Sorting target points by elevation ...", bDebug=bDebug)

    uvs_Flat = [(u,v) for u in range(len(pts)) for v in range(len(pts[0]))]
    pts_Flat = _flattenNestedList(pts)
    zs_Flat_perPts = [pt.Z for pt in pts_Flat]

    zs_Flat_Sorted, uvs_Flat_Sorted = _sort2ListsPer1st(zs_Flat_perPts, uvs_Flat)
    sEval = "uvs_Flat_Sorted[:9]"; print(sEval,'=',eval(sEval))
    sEval = "zs_Flat_Sorted[:9]"; print(sEval,'=',eval(sEval))

    # Group per fElevTol.
    zs_1perElevGroup = [zs_Flat_Sorted[0]]
    uvs_Sorted_and_grouped_per_z = [[]]
    for iZ, z in enumerate(zs_Flat_Sorted):
        z_Lowest_for_group = zs_1perElevGroup[-1] - fElevTol
        if z >= z_Lowest_for_group:
            uvs_Sorted_and_grouped_per_z[-1].append(uvs_Flat_Sorted[iZ])
            continue
        # New group.
        zs_1perElevGroup.append(z)
        uvs_Sorted_and_grouped_per_z.append([uvs_Flat_Sorted[iZ]])

    return zs_1perElevGroup, uvs_Sorted_and_grouped_per_z


def _sortAndGroupTargetsByElevation_High_to_low(pts, fElevTol=1e-6, bDebug=False):
    """
    Parameters
        pts: list of lists of mix of rg.Point3d or None

    Returns
        list of lists of u and v index tuples whose z's are within fElevTol
    """

    _prompt("Sorting target points by elevation ...", bDebug=bDebug)

    pts_Flat = _flattenNestedList(pts)
    zs_Flat_perPts = [pt.Z for pt in pts_Flat]
    zs_Flat_HiToLo = sorted(zs_Flat_perPts, reverse=True)
    iZs_Used = []
    iZs_Sorted_Grouped_Flat = []
    z_LastTolStart = float("inf")
    for i, z in enumerate(zs_Flat_HiToLo):
        if zs_Flat_perPts.index(z) in iZs_Used:
            continue
        iZs_NewElevation = [j for j, x in enumerate(zs_Flat_perPts) if x == z]
        if abs(z_LastTolStart-z) > fElevTol:
            iZs_Sorted_Grouped_Flat.append(iZs_NewElevation)
            z_LastTolStart = z
        else:
            iZs_Sorted_Grouped_Flat[-1].extend(iZs_NewElevation)
        iZs_Used.extend(iZs_NewElevation)

    countV = len(pts[0])
    iUiVs_Sorted = []
    for z_Group in iZs_Sorted_Grouped_Flat:
        iUiVs_Sorted.append([])
        for iZ in z_Group:
            iUiVs_Sorted[-1].append(((iZ // countV), (iZ % countV)))
            #sc.doc.Objects.AddPoint(pts[iU][iV])

    _promptDone(bDebug=False)

    return iUiVs_Sorted


def _getNeighborsPerElevationGroup(ns, uvs_Target_Groups, zs_Targets, zs_Min_AdjustedPerNeighbors, bDebug=False):
    """
    Returns list of lists of tuples of index pairs
    Also modifies zs_Min_AdjustedPerNeighbors
    """

    _prompt("Determining neighbors of targets ...", bDebug=bDebug)

    iUiVs_Sorted_Flat = _flattenNestedList(uvs_Target_Groups)
    #zs_Targets_Flat = _flattenNestedList(zs_Targets)
    uvs_Neighbors_PerElevGroup_Flat = []
    iDirs = -1, 0, 1
    # 8 directions from each index location.
    neighbor_dir_deltas = [[iU_N, iV_N] for iU_N in iDirs for iV_N in iDirs if not (iU_N == iV_N == 0) ]
    uvs_Neighbors_PerElevGroup = []

    for iGroup, uvs_Target_Group in enumerate(uvs_Target_Groups):
        sc.escape_test()
        uvs_Neighbors_PerElevGroup.append([])
        for uT, vT in uvs_Target_Group:
            #print("T:", uT, vT)
            for uD, vD in neighbor_dir_deltas:
                uN = uT + uD
                vN = vT + vD
                #print("N:", uN, vN)
                if (uN, vN) in uvs_Target_Group:
                    # Neighbor cannot be in current target group.
                    #sEval = "iGroup,uD,vD"; print(sEval,'=',eval(sEval))
                    continue
                if uD == 0 and vD == 0:
                    raise Exception("When does this happen? The previous if should have continued this. This condition should be removed from neighbor_dir_deltas.")
                if (uN, vN) in uvs_Neighbors_PerElevGroup_Flat:
                    continue
                if not (2 < uN < (ns.Points.CountU-3)):
                    continue
                if not (2 < vN < (ns.Points.CountV-3)):
                    continue
                uvs_Neighbors_PerElevGroup[-1].append((uN, vN))
                uvs_Neighbors_PerElevGroup_Flat.append((uN, vN))
                z_Min_AdjustedPerNeighbors = zs_Min_AdjustedPerNeighbors[uN][vN]
                z_Target = zs_Targets[uT][vT]
                if z_Min_AdjustedPerNeighbors < z_Target:
                    zs_Min_AdjustedPerNeighbors[uN][vN] = z_Target

    _promptDone(bDebug=bDebug)

    #sEval = "len(uvs_Neighbors_PerElevGroup)"; print(sEval,'=',eval(sEval))
    #sEval = "len(uvs_Neighbors_PerElevGroup[0])"; print(sEval,'=',eval(sEval))
    #sEval = "len(uvs_Neighbors_PerElevGroup[0][0])"; print(sEval,'=',eval(sEval))
    #sEval = "uvs_Neighbors_PerElevGroup[0]"; print(sEval,'=',eval(sEval))
    #sEval = "uvs_Neighbors_PerElevGroup[0][0]"; print(sEval,'=',eval(sEval))
    #1/0

    return uvs_Neighbors_PerElevGroup


def _highestElevation(pts):
    bb = rg.BoundingBox(points=[pt for pt in _flattenNestedList(pts) if pt is not None])
    return bb.Max.Z

    # Alternative routine.
    zMax = -float('inf')
    for us in pts:
        for v in us:
            if v is None:
                continue
            if v.Z > zMax:
                zMax = v.Z
    return zMax


def _fit_NEW_High_to_low(ns_In=None, pts_Target=None, fTolerance=None, bDebug=False):
    """
    WIP to possibly replace other High_to_low
    """

    ns_Out = ns_In.Duplicate() if ns_In.IsDocumentControlled else ns_In

    zs_1perElevGroup, uvs_Sorted_and_grouped_per_z = _sortAndGroupTargetsByElevation_High_to_low_NEW(
        pts_Target,
        fElevTol=0.1*min((sc.doc.ModelAbsoluteTolerance, fTolerance)),
        bDebug=bDebug)

    sEval = "len(zs_1perElevGroup)"; print(sEval,'=',eval(sEval))
    sEval = "len(uvs_Sorted_and_grouped_per_z)"; print(sEval,'=',eval(sEval))
    sEval = "zs_1perElevGroup[:9]"; print(sEval,'=',eval(sEval))
    sEval = "uvs_Sorted_and_grouped_per_z[:9]"; print(sEval,'=',eval(sEval))



def _fit_Iteratively_translate_pts_High_to_low(pts_Target, ns_In, fTolerance, bDebug=False):
    """
    """

    #return _fit_NEW_High_to_low(
    #    ns_In=ns_In,
    #    pts_Target=pts_Target,
    #    fTolerance=fTolerance,
    #    bDebug=bDebug)

    ns_Out = ns_In.Duplicate() if ns_In.IsDocumentControlled else ns_In


    #map(sc.doc.Objects.AddPoint, [pt for pts_V in pts_Target for pt in pts_V]); sc.doc.Views.Redraw()


    uvs_Targets_InElevGroups = _sortAndGroupTargetsByElevation_High_to_low(
        pts_Target,
        fElevTol=0.1*min((sc.doc.ModelAbsoluteTolerance, fTolerance)),
        bDebug=bDebug)

    #if bDebug:
    #    print('-'*20)
    #    sEval = "uvs_Targets_InElevGroups[8][22:27]"; print(sEval,'=',eval(sEval))
    #    sEval = "uvs_Targets_InElevGroups[9][22:29]"; print(sEval,'=',eval(sEval))
    #    sEval = "uvs_Targets_InElevGroups[9][22:27]"; print(sEval,'=',eval(sEval))


    if bDebug: sEval = 'len(uvs_Targets_InElevGroups)'; print(sEval+':',eval(sEval))
    if bDebug: sEval = 'uvs_Targets_InElevGroups'; print(sEval+':',eval(sEval))
    #return

    zs_HighestPerElevGroup = []
    for uvsGroup in uvs_Targets_InElevGroups:
        zs = []
        for u,v in uvsGroup:
            zs.append(pts_Target[u][v].Z)
        zs_HighestPerElevGroup.append(max(zs))

    #if bDebug:
    #    print('-'*20)
    #    sEval = "zs_HighestPerElevGroup"; print(sEval,'=',eval(sEval))
    #    sEval = "len(zs_HighestPerElevGroup)"; print(sEval,'=',eval(sEval))

    #    print('-'*20)
    #    sEval = "uvs_Targets_InElevGroups[8][22:27]"; print(sEval,'=',eval(sEval))
    #    sEval = "uvs_Targets_InElevGroups[9][22:29]"; print(sEval,'=',eval(sEval))
    #    sEval = "uvs_Targets_InElevGroups[9][22:27]"; print(sEval,'=',eval(sEval))

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

    #if bDebug:
    #    print('-'*20)
    #    sEval = "zs_Targets[8][22:27]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Targets[9][22:29]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Targets[9][22:27]"; print(sEval,'=',eval(sEval))


    # Set zs_Targets in elevation groups to highest elevation.
    # This allows better control point selection in resultant surface.
    for iGroup, uvsGroup in enumerate(uvs_Targets_InElevGroups):
        for u, v in uvsGroup:
            zs_Targets[u][v] = zs_HighestPerElevGroup[iGroup]

    #if bDebug:
    #    print('-'*20)
    #    sEval = "zs_Targets[8][22:27]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Targets[9][22:29]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Targets[9][22:27]"; print(sEval,'=',eval(sEval))

    zs_Min_AdjustedPerNeighbors = [zsV[:] for zsV in zs_Targets]

    #if bDebug:
    #    print('-'*20)
    #    sEval = "zs_Min_AdjustedPerNeighbors[8][22:27]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Min_AdjustedPerNeighbors[9][22:29]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Min_AdjustedPerNeighbors[9][22:27]"; print(sEval,'=',eval(sEval))


    uvs_Neighbors_PerElevGroup = _getNeighborsPerElevationGroup(
        ns_In,
        uvs_Targets_InElevGroups,
        zs_Targets,
        zs_Min_AdjustedPerNeighbors,
        bDebug=bDebug)

    #if bDebug:
    #    print('-'*20)
    #    sEval = "zs_Min_AdjustedPerNeighbors[8][22:27]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Min_AdjustedPerNeighbors[9][22:29]"; print(sEval,'=',eval(sEval))
    #    sEval = "zs_Min_AdjustedPerNeighbors[9][22:27]"; print(sEval,'=',eval(sEval))


    #print(uvs_Targets_InElevGroups[:10])
    if bDebug: sEval = 'len(uvs_Neighbors_PerElevGroup)'; print(sEval+':',eval(sEval))


    # TESTING: TODO: Remove this comment if this block of code is successful.
    #zMax = _highestElevation(pts_Target)
    #sEval = 'zMax'; print(sEval+':',eval(sEval))

    #for iU in range(ns_Out.Points.CountU):
    #    for iV in range(ns_Out.Points.CountV):
    #        cp = ns_Out.Points.GetControlPoint(iU, iV)
    #        cp.Z = zMax
    #        ns_Out.Points.SetControlPoint(iU, iV, cp)

    #sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw()


    # Set highest elevation group


    zMax = _highestElevation(pts_Target)


    # Set highest elevation group to max Z.
    for (iU, iV) in uvs_Targets_InElevGroups[0]:
        #sEval = 'iU,iV'; print(sEval+':',eval(sEval))
        cp = ns_In.Points.GetControlPoint(iU, iV)
        cp.Z = zMax
        ns_Out.Points.SetControlPoint(iU, iV, cp)


    # Set highest elevation group neighbors to max Z.
    for (iU, iV) in uvs_Neighbors_PerElevGroup[0]:
        #sEval = 'iU,iV'; print(sEval+':',eval(sEval))
        cp = ns_In.Points.GetControlPoint(iU, iV)
        #zs_Min_AdjustedPerNeighbors[iU][iV] = cp.Location.Z
        cp.Z = zMax #z_Low
        ns_Out.Points.SetControlPoint(iU, iV, cp)


    if bDebug: sEval = 'len(zs_Min_AdjustedPerNeighbors)'; print(sEval+':',eval(sEval))
    if bDebug: sEval = 'zs_Min_AdjustedPerNeighbors'; print(sEval+':',eval(sEval))
    #return


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
            #if (uT,vT) == (6,5):
            #    pass
            #    sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); return
            zG = ns_Out.PointAt(uG, vG).Z
            zT = pts_Target[uT][vT].Z
            if (zG + 0.001*fTolerance) >= zT:
                continue
            #sc.doc.Objects.AddSurface(ns_Out); sc.doc.Views.Redraw(); 1/0
            raise Exception("This still occurs. When?")
            break
        else:
            # No change to zs_Min_AdjustedPerNeighbors.

            #if bDebug:
            #    _prompt("All Grevilles are on or above target.", bDebug=False)
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
                        print(uD, vD)
                        uNN = uN + uD
                        vNN = vN + vD
                        zs_Neighbors.append(zs_Targets[uNN][vNN])
                    1/0
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


def _fit_Iteratively_translate_individual_pts(pts, ns_In, fTolerance):
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
            sc.doc.Objects.AddPoint(rgObj)


def _closestPointsOfNeighborsOnNormalLines(pts_In, pts_Greville, iU, iV, iMinNeighborCt=1, bDiag=True, bLineExts=True):
    """
    bLineExts:
        When True: If neighbor's point and its neighbor's point in the same
        direction are both available, create a line through those points
        and get the ClosestPoint of that line on the normal line.
    """

    lines_thruStartingSrfGrevilles = _createZAxisLinesAtGrevilles(pts_Greville)

    #for col in lines_thruStartingSrfGrevilles:
    #    for line in col:
    #        sc.doc.Objects.AddLine(line)
    #sc.doc.Views.Redraw(); 1/0


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
                lineA=lines_thruStartingSrfGrevilles[iU][iV],
                lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
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
                    lineA=lines_thruStartingSrfGrevilles[iU][iV],
                    lineB=line_ThruNeighbors)
            if rc[0]:
                pt = lines_thruStartingSrfGrevilles[iU][iV].PointAt(rc[1])
        else:
            pt = lines_thruStartingSrfGrevilles[iU][iV].ClosestPoint(
                    pts_In[iU+1][iV+1],
                    limitToFiniteSegment=False)
        if pt: pts_Out.append(pt)

    if len(pts_Out) < iMinNeighborCt:
        return []

    return pts_Out


def _addMissingPointsAlongBorder(pts_In, pts_Greville, idxs_pt_filter=None, iMinNeighborCt=1, bDiag=False, bLineExts=True):
    """
    Parameters
        pts_In: list of lists of Point3d

    Returns
        list of lists of Point3d if modified.
    """

    # Modify a copy of the list so that new points do not affect subsequent ones in this function call.
    pts_Out = [ptsV[:] for ptsV in pts_In]

    bModificationOccured = False

    for iU in range(len(pts_In)):
        for iV in range(len(pts_In[0])):
            if idxs_pt_filter and not (iU, iV) in idxs_pt_filter: continue

            if pts_Out[iU][iV] is not None: continue

            pts = _closestPointsOfNeighborsOnNormalLines(
                pts_Out,
                pts_Greville,
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
    #for iU in range(len(pts0)):
    #    for iV in range(len(pts0[0])):
    #        pts0[iU][iV] = pts_Out[iU][iV]

    if bModificationOccured:
        return pts_Out


def _addSrfGrevillPtsForMissing(pts_In, ns_Starting):
    """
    Returns
        list of lists of Point3d
    """
    pts_Out = []
    for iU in range(len(pts_In)):
        pts_Out.append([])
        for iV in range(len(pts_In[0])):
            if pts_In[iU][iV] is None:
                #print("Missing point found at {}, {}".format(iU, iV))
                pts_Out[iU][iV] = ns_Starting.Points.GetControlPoint(iU,iV).Location
                u, v = ns.Points.GetGrevillePoint(iU, iV)
                pt = ns.PointAt(u, v)
                pts_Out[-1].append(pt)
            else:
                pts_Out[-1].append(pts_In[iU][iV])
    return pts_Out


def _extrapolateHitsForMisses(pts_Target_In, pts_Greville):

    idxs_borderPts = _getBorderPointIndices(pts_Target_In)
    #print(idxs_borderPts); 1/0

    attr = rd.ObjectAttributes()
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject

    #for bDiag1, bDiag2, bLineExts1, bLineExts2 in itertools.product(
    #            (False, True), (False, True), (False, True), (False, True)):
    for bLineExts1, bLineExts2 in ((True, False),):
    #for bLineExts1, bLineExts2 in itertools.product(
    #            (False, True), (False, True)):

        pts_Target_Out = [ptsV[:] for ptsV in pts_Target_In]


        #for iMinNeighborCt in 4,3,2,1: #(1,): #

        #    addMissingPointsAlongBorder(
        #            pts_Target,
        #            idxs_pt_filter=idxs_borderPts,
        #            iMinNeighborCt=iMinNeighborCt,
        #            bDiag=True,
        #            bLineExts=bLineExts1)

        #    if not hasMissingPoints(pts_Target):
        #        break


        while _hasMissingPoints(pts_Target_Out):
            sc.escape_test()

            pts_Target_Out = _addMissingPointsAlongBorder(
                pts_Target_Out,
                pts_Greville,
                bDiag=False,
                bLineExts=bLineExts2)

            #sEval = "bPointsWereAdded"; print(sEval,'=',eval(sEval))


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
        for iU in range(len(pts_Target_Out)):
            for iV in range(len(pts_Target_Out[0])):
                if pts_Target_Out[iU][iV] is None:
                    raise Exception("Why are points missing?")
                    pts_Target_Out[iU][iV] = rg.Point3d(
                        ns_WIP.Points.GetControlPoint(iU,iV).X,
                        ns_WIP.Points.GetControlPoint(iU,iV).Y,
                        ns_WIP.Points.GetControlPoint(0,0).Z)

        return pts_Target_Out


def processGeometry(rgObjs_toDrapeOver, srf_Starting, cPlane=rg.Plane.WorldXY, fTolerance=None, iTargetMisses=1, bDebug=False):
    """
    Parameters
        iTargetMisses: int  0 to fix on starting surface, 1 to use lowest hit neighbor, 2 to extrapolate linearly from nearest hits
    """

    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance


    data = Data()

    ns_WIP = srf_Starting.ToNurbsSurface()

    #sEval = "cPlane"; print(sEval,'=',eval(sEval))
    #return

    if cPlane == rg.Plane.WorldXY:
        xform_toW = xform_fromW = None
    else:
        xform_toW = rg.Transform.PlaneToPlane(cPlane, rg.Plane.WorldXY)
        if not all([rgObjs_toDrapeOver[i].Transform(xform_toW) for i in range(len(rgObjs_toDrapeOver))]):
            raise Exception("No all objects to drape over can be transformed.")
        if not ns_WIP.Transform(xform_toW):
            raise Exception("Starting surface cannot be transformed.")

        # Prepare xform for output.
        xform_fromW = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, cPlane)
        #[sc.doc.Objects.AddBrep(rgBreps_ProjectTo[i]) for i in range(len(rgBreps_ProjectTo))]
        #sc.doc.Objects.AddSurface(ns_WIP); sc.doc.Views.Redraw()#; return


    uvs_All = [(u, v)
               for u in range(ns_WIP.Points.CountU)
               for v in range(ns_WIP.Points.CountV)]


    #ns_WIP = addKnotsToSurface(ns_WIP)


    pts_Greville = _getGrevillePoints(ns_WIP)

    #[sc.doc.Objects.AddPoint(pts_Greville[u][v]) for u, v in uvs_All if pts_Greville[u][v] is not None]
    #sc.doc.Views.Redraw(); return

    pts_Target = _projectPts_toObjs(
        pts_In=pts_Greville,
        rgObjs_Targets=rgObjs_toDrapeOver,
        bDebug=bDebug,
        )
    if not pts_Target:
        print("Projected points were not obtained.")
        return


    #[_addObject(pts_Target[u][v], xform_fromW) for u, v in uvs_All if pts_Target[u][v] is not None]
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



    # Set all the control points of the output surface to the highest elevation.
    # It is to simulate flattening then dropping the material onto the objects.
    zMax = _highestElevation(pts_Target)


    # This checks for ns_WIP being planar, has normal parallel to World Z axis, and origin Z at max Z of targets.
    # The starting surface created by this script is this except the Z coordinate is the top of the bounding box of the objects that are draped over.
    #def is_NS_flattened_at_top_of_target(ns, z):
    #    bSuccess, plane = ns.TryGetPlane()
    #    if not bSuccess:
    #        return False
    #    if (plane.Origin.Z - z) >= 1e-6:
    #        return False
    #    iIsParallelTo = plane.Normal.IsParallelTo(
    #        other=rg.Vector3d.ZAxis,
    #        angleTolerance=Rhino.RhinoMath.ToRadians(1e-6))
    #    return iIsParallelTo == 1


    #print(is_NS_flattened_at_top_of_target(ns_WIP, zMax))
    #return

    for iU in range(ns_WIP.Points.CountU):
        for iV in range(ns_WIP.Points.CountV):
            cp = ns_WIP.Points.GetControlPoint(iU, iV)
            cp.Z = zMax
            ns_WIP.Points.SetControlPoint(iU, iV, cp)



    if _hasMissingPoints(pts_Target):

        # Making a duplicate for debugging, etc.
        pts_Target_HasMissing = [ptsV[:] for ptsV in pts_Target]


        if iTargetMisses==0:
            # Fill any remaining missing points with the
            # starting surface control point locations.
            pts_Target = _addSrfGrevillPtsForMissing(
                pts_Target_HasMissing,
                ns_WIP)
        elif iTargetMisses==1:
            raise Exception("Need to implement lowest neighbor. Change setting to FixToStartingSrf or LinearlyExtrapolateFromHits.")
        elif iTargetMisses==2:
            pts_Target = _extrapolateHitsForMisses(
                pts_Target_HasMissing,
                pts_Greville)
            #for col in pts_Target:
            #    for pt in col:
            #        sc.doc.Objects.AddPoint(pt)
            #sc.doc.Views.Redraw()
            #return
        else:
            raise Exception("iTargetMisses must be 0, 1, or 2.")


    ns_Out = _fit_Iteratively_translate_pts_High_to_low(
        pts_Target,
        ns_WIP,
        fTolerance,
        bDebug=bDebug)

    # TODO: Reevaluate using a different approach for fitting to a single surface.
    #if (
    #    len(rgObjs_toDrapeOver) == 1 and
    #    isinstance(rgObjs_toDrapeOver, rg.Brep) and
    #    rgObjs_toDrapeOver[0].Faces.Count == 1
    #):
    #    ns_Out = _fit_Iteratively_translate_individual_pts(
    #        pts_Target,
    #        ns_WIP,
    #        fTolerance,
    #        bDebug=bDebug)
    #else:
    #    ns_Out = _fit_Iteratively_translate_pts_High_to_low(
    #        pts_Target,
    #        ns_WIP,
    #        fTolerance,
    #        bDebug=bDebug)

    if xform_fromW:
        ns_Out.Transform(xform_fromW)

    return ns_Out


def processDocObject(rhObjs_toDrapeOver, objref_srf_Starting=None, cPlane=rg.Plane.WorldXY, **kwargs):
    """
    Parameters
        rhObjs_toDrapeOver: rd.ObjRef[] or rd.RhinoObject[]
        iTargetMisses: int  0 to fix on starting surface, 1 to use lowest hit neighbor, 2 to extrapolate linearly from nearest hits
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    fSpanSpacing = getOpt('fSpanSpacing')
    iSpansBeyondEachSide = getOpt('iSpansBeyondEachSide')
    iTargetMisses = getOpt('iTargetMisses')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    rgObjs_toDrapeOver = []
    for rhObj in rhObjs_toDrapeOver:
        rgObj_toDrapeOver = rs.coercegeometry(rhObj)
        if not isinstance(rgObj_toDrapeOver, (rg.Brep, rg.Mesh)):
            print("{} is not supported for an object to drape over.".format(rgObj_toDrapeOver.GetType().Name))
            continue
        rgObjs_toDrapeOver.append(rgObj_toDrapeOver)


    if objref_srf_Starting is None:
        ns_Starting = _createStartingSurface(
            rgObjs_toDrapeOver,
            cPlane=cPlane,
            fSpanSpacing=fSpanSpacing,
            iSpansBeyondEachSide=iSpansBeyondEachSide,
            bDebug=bDebug)
        #sc.doc.Objects.AddSurface(ns_Starting); sc.doc.Views.Redraw(); return
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
        iTargetMisses=iTargetMisses,
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

    objrefs_toDrapeOver = getInput_ObjsToDrapeOver()
    if objrefs_toDrapeOver is None: return

    cPlane = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()

    bFlipCPlane = Opts.values['bFlipCPlane']
    fTolerance = Opts.values['fTolerance']
    bUserProvidesStartingSrf = Opts.values['bUserProvidesStartingSrf']
    fSpanSpacing = Opts.values['fSpanSpacing']
    iSpansBeyondEachSide = Opts.values['iSpansBeyondEachSide']
    iTargetMisses = Opts.values['iTargetMisses']
    bDeleteStartingSrf = Opts.values['bDeleteStartingSrf']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if bFlipCPlane:
        cPlane.Flip()

    if not bUserProvidesStartingSrf:
        objref_srf_Starting = None
    else:
        gObjs_toDrapeOver = [_.ObjectId for _ in objrefs_toDrapeOver]

        sc.doc.Objects.UnselectAll()
        sc.doc.Views.Redraw()

        objref_srf_Starting = getInput_StartingSurface(gObjs_toDrapeOver)
        if objref_srf_Starting is None: return

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.CommandPrompt = "Working ..."

    processDocObject(
        rhObjs_toDrapeOver=objrefs_toDrapeOver,
        objref_srf_Starting=objref_srf_Starting,
        cPlane=cPlane,
        fTolerance=fTolerance,
        fSpanSpacing=fSpanSpacing,
        iSpansBeyondEachSide=iSpansBeyondEachSide,
        iTargetMisses=iTargetMisses,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()