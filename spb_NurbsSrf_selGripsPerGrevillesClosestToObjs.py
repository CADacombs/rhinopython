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
250127: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from clr import StrongBox


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

    key = 'fDistTolForSharedGrips'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
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

        if key == 'fDistTolForSharedGrips':
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

    go.SetCommandPrompt("Select curves and/or pts for references")

    go.GeometryFilter = rd.ObjectType.Point | rd.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=False)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        #key = 'AllSrfGripsOn'; idxs_Opt[key] = go.AddOption(key)
        #addOption('bAutoPickSrf')
        addOption('iClosestGripCtPerAnalysisPt')
        addOption('fDistTolForSharedGrips')
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
            key = 'fDistTolForSharedGrips'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _printGripsAddedFeedback(rdGrips_Added, rdGrips_Ref):
    if not rdGrips_Added:
        print("No grips were added to the selection of {}.".format(
            len(rdGrips_Ref)))
    else:
        print("{} grips were added to the selection of {} for a total of {}.".format(
            len(rdGrips_Added), len(rdGrips_Ref), len(rdGrips_Added) + len(rdGrips_Ref)))


def _find_indices(lst, element):
    return [i for i, _ in enumerate(lst) if _ == element]


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

    idxs_Out = []
    for j in _find_indices(dists, min(dists)):
        if j not in idxs_Out:
            idxs_Out.append(j)

    return idxs_Out


def selectGrips(objrefs_Ref, rdBrep_1F, iClosestGripCtPerAnalysisPt=1, fDistTolForSharedGrips=None, bEcho=True, bDebug=False):
    """
    """

    rgB = rdBrep_1F.BrepGeometry
    rgS = rgB.Faces[0].UnderlyingSurface()
    if not isinstance(rgS, rg.NurbsSurface):
        if bEcho:
            print("Non-NurbsSurface passed to selectGrips. Nothing will be done.")
        return

    ns = rgS

    if fDistTolForSharedGrips is None:
        fDistTolForSharedGrips = sc.doc.ModelAbsoluteTolerance

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


    idx_Grips_to_sel = []


    for objref in objrefs_Ref:
        rgO = objref.Object().Geometry
        if isinstance(rgO, rg.Point):
            idxs_Grips = _getGripsIndices_GripsClosestToPoint(
                rgO,
                idxs_UVs=idxs_UVs,
                pt3ds_Grevs=pt3ds_Grevs,
                levels=iClosestGripCtPerAnalysisPt,
                tol_per_level=fDistTolForSharedGrips)
            idx_Grips_to_sel.extend(idxs_Grips)


    if not idx_Grips_to_sel:
        return

    rdBrep_1F.GripsOn = True

    rdGrips = rdBrep_1F.GetGrips()

    for iG in idx_Grips_to_sel:
        rdGrips[iG].Select(
            on=True,
            syncHighlight=True,
            persistentSelect=True,
            ignoreGripsState=True,
            ignoreLayerLocking=False,
            ignoreLayerVisibility=False)

    return


    for rdGrip_Ref in rdGrips_Ref:
        rdGrip_Ref.Select(True, syncHighlight=True) # For persistant selection.

    for rdGrip in rdGrips_Added_All:
        rdGrip.Select(True, syncHighlight=True) # For persistant selection.

    if bEcho:
        _printGripsAddedFeedback(rdGrips_Added_All, rdGrips_Ref)


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
    fDistTolForSharedGrips = Opts.values['fDistTolForSharedGrips']
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
        fDistTolForSharedGrips=fDistTolForSharedGrips,
        bEcho=bEcho,
        bDebug=bDebug,
    )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()