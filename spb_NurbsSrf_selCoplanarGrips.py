"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
210801-02: Created.
211106: Tolerance is now creeped up from RhinoMath.ZeroTolerance to select only grips that
        are collectively within the tolerance, not independently with reference grips.
        Bug fixed during numeric input.
250125-27: Added selection of grips per one grip and parallel plane per option.
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


    key = 'bDeterminePlaneBy3PlusGripPick'; keys.append(key)
    values[key] = False
    names[key] = 'DeterminePlaneBy'
    riOpts[key] = ri.Custom.OptionToggle(values[key], '1GripWith||Plane', 'AtLeast3Grips')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iParallelWithPlane'; keys.append(key)
    values[key] = 0
    listValues[key] = (
        'CPlaneXY',
        'CPlaneYZ',
        'CPlaneZX',
        'WorldXY',
        'WorldYZ',
        'WorldZX',
        )
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fPlanarityTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bSearchOtherSrfs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEnableGrips'; keys.append(key)
    values[key] = False
    names[key] = 'OtherSrfsGrips'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'VisibleOnly', 'EnableAsNeeded')
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

        if key == 'fPlanarityTol':
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


def _getAllNormalBreps():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))
    return [rdB for rdB in sc.doc.Objects.GetObjectList(oes) if rdB.BrepGeometry.IsSurface]


def getInput():
    """
    Get grips with optional input.
    """

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go = ri.Custom.GetObject()

        if Opts.values['bDeterminePlaneBy3PlusGripPick']:
            go.SetCommandPrompt("Select at least 3 grips for coplanar reference")
        else:
            go.SetCommandPrompt("Select 1 grip for coplanar reference")

        go.GeometryFilter = rd.ObjectType.Grip

        go.AcceptNumber(True, acceptZero=False)

        bResetGetObject = False

        while True:
            go.ClearCommandOptions()
            idxs_Opt.clear()

            key = 'AllSrfGripsOn'; idxs_Opt[key] = go.AddOption(key)
            addOption('bDeterminePlaneBy3PlusGripPick')
            if not Opts.values['bDeterminePlaneBy3PlusGripPick']:
                addOption('iParallelWithPlane')
            addOption('fPlanarityTol')
            addOption('bSearchOtherSrfs')
            if Opts.values['bSearchOtherSrfs']:
                addOption('bEnableGrips')
            addOption('bEcho')
            addOption('bDebug')

            res = go.GetMultiple(
                minimumNumber=3 if Opts.values['bDeterminePlaneBy3PlusGripPick'] else 1,
                maximumNumber=0 if Opts.values['bDeterminePlaneBy3PlusGripPick'] else 1)

            if res == ri.GetResult.Cancel:
                go.Dispose()
                return

            if res == ri.GetResult.Object:
                objrefs = go.Objects()
                go.Dispose()
                return objrefs

            if res == ri.GetResult.Number:
                key = 'fPlanarityTol'
                Opts.riOpts[key].CurrentValue = go.Number()
                Opts.setValue(key)
                continue

            if go.Option().Index == idxs_Opt['AllSrfGripsOn']:
                bRedraw = False
                for rdB in _getAllNormalBreps():
                    if not rdB.GripsOn:
                        rdB.GripsOn = True # This doesn't work on polysurface breps.
                        if rdB.GripsOn:
                            bRedraw = True
                if bRedraw:
                    sc.doc.Views.Redraw()
                continue

            # An option was selected.
            for key in idxs_Opt:
                if go.Option().Index == idxs_Opt[key]:
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)
                    if key == 'bDeterminePlaneBy3PlusGripPick':
                        bResetGetObject = True
                    break

            if bResetGetObject:
                go.Dispose()
                break


def _colinearDeviation(pts):
    line = rg.Line(pts[0], pts[-1])
    #sc.doc.Objects.AddLine(line_EndEnd)
    devs = []
    for iP in xrange(1, len(pts)-1):
        dev = line.DistanceTo(
            testPoint=pts[iP],
            limitToFiniteSegment=True)
        devs.append(dev)
    return max(devs)


def _areColinear(pts, fLinearTol):
    line = rg.Line(pts[0], pts[-1])
    #sc.doc.Objects.AddLine(line_EndEnd)
    devs = []
    for iP in xrange(1, len(pts)-1):
        fDist = line.DistanceTo(
            testPoint=pts[iP],
            limitToFiniteSegment=True)
        if fDist > fLinearTol:
            return False
    return True


def _tryGetPlane(ref_pts, fPlanarityTol):
    strongbox_plane = StrongBox[rg.Plane]()
    strongbox_double = StrongBox[float]()

    plane_fit_res = rg.Plane.FitPlaneToPoints(ref_pts, strongbox_plane, strongbox_double)
    if plane_fit_res != rg.PlaneFitResult.Success:
        return (
            False,
            None,
            "Plane could not be fit through reference point.  {}".format(plane_fit_res)
            )

    planarity = strongbox_double.Value
    if planarity > fPlanarityTol:
        return (
            False,
            None,
            "Plane fit through reference points is out of tolerance ({}).".format(planarity)
            )

    plane = strongbox_plane.Value

    if not plane.IsValid:
        return (
            False,
            None,
            "Reference points may be colinear."
            )

    fColinearTol = 10.0 * sc.doc.ModelAbsoluteTolerance

    fColinearDev = _colinearDeviation(ref_pts)
    if fColinearDev <= fColinearTol:
        return (
            False,
            None,
            "Reference points are colinear within {}.".format(fColinearDev)
            )

    # Retired: Faster routine but doesn't provide deviation distance.
    #if _areColinear(ref_pts, fColinearTol):
    #    return (
    #        False,
    #        None,
    #        "Reference points are colinear within the {} tolerance.".format(fColinearTol)
    #        )

    return True, plane, None


def _findCoplanarGrips_Per_1Grip_and_parallelPlane(rhObj, rdGrip_Ref, plane_Parallel, fPlanarityTol):

    if not rhObj.GripsOn:
        rhObj.GripsOn = True
        if not rhObj.GripsOn:
            return

    rdGrips_All_of_obj = rhObj.GetGrips()
    if rdGrips_All_of_obj is None:
        return

    plane_Ref = rg.Plane(plane_Parallel)
    plane_Ref.Origin = rdGrip_Ref.CurrentLocation

    #sc.doc.Objects.AddSurface(rg.PlaneSurface(plane_Ref))

    rdGrips_Added = []

    sEval = "fPlanarityTol"; print(sEval,'=',eval(sEval))

    for rdGrip in rdGrips_All_of_obj:

        if rdGrip.Id == rdGrip_Ref.Id:
            continue

        pt = rdGrip.CurrentLocation

        dist = abs(plane_Ref.DistanceTo(pt)) # abs because this method returns a signed distance.

        if dist <= fPlanarityTol:
            #sEval = "dist"; print(sEval,'=',eval(sEval))
            rdGrips_Added.append(rdGrip)

    return rdGrips_Added


def _findCoplanarGrips_AddingWithRefGrips(rhObj, rdGrips_Ref, fPlanarityTol):
    """
    The tolerance creeps from ZeroTolerance to fPlanarityTol so that grips are added more in the
    order of how coplanar they are will grips already in list, and less by their
    order from GetGrips().
    """
    if not rhObj.GripsOn:
        rhObj.GripsOn = True
        if not rhObj.GripsOn:
            return

    # APC = ArePointsCoplanar.
    rdGrips_Passed = []
    gGrips_Ref = [_.Id for _ in rdGrips_Ref]
    pts_Passed = [_.CurrentLocation for _ in rdGrips_Ref]

    rdGrips_toCheck = [_ for _ in rhObj.GetGrips() if _.Id not in gGrips_Ref]
    if rdGrips_toCheck is None:
        return

    iExponent = 0
    fTol = Rhino.RhinoMath.ZeroTolerance
    bMaxTolReach = fPlanarityTol < fTol

    while not bMaxTolReach:
        sc.escape_test()

        fTol = Rhino.RhinoMath.ZeroTolerance*(10**(iExponent))

        if fTol >= fPlanarityTol:
            fTol = fPlanarityTol
            bMaxTolReach = True

        rdGrips_Failed = []

        for rdGrip in rdGrips_toCheck:
            pt = rdGrip.CurrentLocation

            if rg.Point3d.ArePointsCoplanar(pts_Passed+[pt], tolerance=fTol):
                rdGrips_Passed.append(rdGrip)
                pts_Passed.append(pt)
            else:
                rdGrips_Failed.append(rdGrip)

            if not rdGrips_Failed:
                return rdGrips_Passed

        rdGrips_toCheck = rdGrips_Failed

        iExponent += 1

    return rdGrips_Passed


def _printGripsAddedFeedback(rdGrips_Added, rdGrips_Ref):
    if not rdGrips_Added:
        print("No grips were added to the selection of {}.".format(
            len(rdGrips_Ref)))
    else:
        print("{} grips were added to the selection of {} for a total of {}.".format(
            len(rdGrips_Added), len(rdGrips_Ref), len(rdGrips_Added) + len(rdGrips_Ref)))


def selectGrips_planePerMultipleGrips(objref_Grips_Ref, fPlanarityTol=Rhino.RhinoMath.ZeroTolerance, bSearchOtherSrfs=False, bEnableGrips=True, bEcho=True, bDebug=False):
    """
    """

    rdGrips_Ref = []
    pts_Grips_Refs = []
    gOwners_of_RefGrips = []

    for o in objref_Grips_Ref:
        rdGrip_Ref = o.Object()
        rdGrips_Ref.append(rdGrip_Ref)
        pts_Grips_Refs.append(rdGrip_Ref.CurrentLocation)
        gOwners_of_RefGrips.append(rdGrip_Ref.OwnerId)

    #sEval = "gOwners_of_RefGrips"; print(sEval,'=',eval(sEval))

    bSuccess, plane, sLog = _tryGetPlane(pts_Grips_Refs, fPlanarityTol)
    if not bSuccess:
        if bEcho and sLog: print(sLog)
        return

    rdGrips_Added_All = []

    for gOwner in set(gOwners_of_RefGrips):
        rdOwner = sc.doc.Objects.FindId(gOwner)
        rdGrips_Added_this_object = _findCoplanarGrips_AddingWithRefGrips(
            rdOwner,
            rdGrips_Ref=rdGrips_Ref,
            fPlanarityTol=fPlanarityTol,
            )
        rdGrips_Added_All.extend(rdGrips_Added_this_object)

    if bSearchOtherSrfs:
        iter = rd.ObjectEnumeratorSettings()
        iter.LockedObjects = False

        for rdObj in sc.doc.Objects.GetObjectList(iter):
            if rdObj.ObjectType != rd.ObjectType.Brep:
                continue
            if rdObj.BrepGeometry.Faces.Count > 1:
                continue
            if rdObj.Id in gOwners_of_RefGrips:
                continue
            bGripsWereEnabled = rdObj.GripsOn
            if not bEnableGrips and not bGripsWereEnabled:
                continue

            rdGrips_Added_this_object = _findCoplanarGrips_AddingWithRefGrips(
                rdObj,
                rdGrips_Ref=rdGrips_Ref,
                fPlanarityTol=fPlanarityTol,
                )

            if rdGrips_Added_this_object:
                rdGrips_Added_All.extend(rdGrips_Added_this_object)
            else:
                if not bGripsWereEnabled:
                    rdObj.GripsOn = False

    for rdGrip_Ref in rdGrips_Ref:
        rdGrip_Ref.Select(True, syncHighlight=True) # For persistant selection.

    for rdGrip in rdGrips_Added_All:
        rdGrip.Select(True, syncHighlight=True) # For persistant selection.

    if bEcho:
        _printGripsAddedFeedback(rdGrips_Added_All, rdGrips_Ref)


def _getParallelToPlane(iParallelWithPlane):
    if Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'CPlaneXY':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()
    if Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'CPlaneYZ':
        cplane = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()
        return rg.Plane(origin=cplane.Origin, normal=cplane.XAxis)
    if Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'CPlaneZX':
        cplane = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()
        return rg.Plane(origin=cplane.Origin, normal=-cplane.YAxis)
    elif Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'WorldXY':
        return rg.Plane.WorldXY
    elif Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'WorldYZ':
        return rg.Plane.WorldYZ
    elif Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'WorldZX':
        return rg.Plane.WorldZX
    #elif Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'View':
    #    return sc.doc.Views.ActiveView.ActiveViewport.GetCameraFrame()[1].ZAxis
    #elif Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'Custom':
    #    return Opts.values['vectCustom']
    raise Exception("What happened?")


def selectGrips_Per_1Grip_and_parallelPlane(objref_Grip_Ref, iParallelWithPlane=3, fPlanarityTol=Rhino.RhinoMath.ZeroTolerance, bSearchOtherSrfs=False, bEnableGrips=True, bEcho=True, bDebug=False):
    """
    """

    rdGrip_Ref = objref_Grip_Ref.Object()
    rdOwner = sc.doc.Objects.FindId(rdGrip_Ref.OwnerId)

    plane_Parallel = _getParallelToPlane(Opts.values['iParallelWithPlane'])

    rdGrips_Added_All = _findCoplanarGrips_Per_1Grip_and_parallelPlane(
        rdOwner,
        rdGrip_Ref=rdGrip_Ref,
        plane_Parallel=plane_Parallel,
        fPlanarityTol=fPlanarityTol,
        )

    if bSearchOtherSrfs:
        iter = rd.ObjectEnumeratorSettings()
        iter.LockedObjects = False

        for rdObj in sc.doc.Objects.GetObjectList(iter):
            if rdObj.ObjectType != rd.ObjectType.Brep:
                continue
            if rdObj.BrepGeometry.Faces.Count > 1:
                continue
            if rdObj.Id == rdGrip_Ref.OwnerId:
                continue
            bGripsWereEnabled = rdObj.GripsOn
            if not bEnableGrips and not bGripsWereEnabled:
                continue

            rdGrips_Added_this_object = _findCoplanarGrips_Per_1Grip_and_parallelPlane(
                rdObj,
                rdGrip_Ref=rdGrip_Ref,
                plane_Parallel=plane_Parallel,
                fPlanarityTol=fPlanarityTol,
                )

            if rdGrips_Added_this_object:
                rdGrips_Added_All.extend(rdGrips_Added_this_object)
            else:
                if not bGripsWereEnabled:
                    rdObj.GripsOn = False

    rdGrip_Ref.Select(True, syncHighlight=True) # For persistant selection.

    for rdGrip in rdGrips_Added_All:
        rdGrip.Select(True, syncHighlight=True) # For persistant selection.

    if bEcho:
        _printGripsAddedFeedback(rdGrips_Added_All, [rdGrip_Ref])


def main():

    objrefx = getInput()
    if objrefx is None: return

    bDeterminePlaneBy3PlusGripPick = Opts.values['bDeterminePlaneBy3PlusGripPick']
    iParallelWithPlane = Opts.values['iParallelWithPlane']
    fPlanarityTol = Opts.values['fPlanarityTol']
    bSearchOtherSrfs = Opts.values['bSearchOtherSrfs']
    bEnableGrips = Opts.values['bEnableGrips']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    if bDeterminePlaneBy3PlusGripPick:

        selectGrips_planePerMultipleGrips(
            objrefx,
            fPlanarityTol=fPlanarityTol,
            bSearchOtherSrfs=bSearchOtherSrfs,
            bEnableGrips=bEnableGrips,
            bEcho=bEcho,
            bDebug=bDebug,
        )
    else:
        selectGrips_Per_1Grip_and_parallelPlane(
            objrefx[0],
            iParallelWithPlane=iParallelWithPlane,
            fPlanarityTol=fPlanarityTol,
            bSearchOtherSrfs=bSearchOtherSrfs,
            bEnableGrips=bEnableGrips,
            bEcho=bEcho,
            bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()