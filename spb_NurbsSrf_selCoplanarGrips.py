"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
210801-02: Created.
211106: Tolerance is now creeped up from RhinoMath.ZeroTolerance to select only grips that
        are collectively within the tolerance, not independently with reference grips.
        Bug fixed during numeric input.
250125-26: WIP: Added selection of grips per one grip and parallel plane per option.

Control points at both 1 u and 1 v in from each side will be made coplanar with adjacent
control points along border of surface.

WIP:
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

            key = 'EnableAllSrfGrips'; idxs_Opt[key] = go.AddOption(key)
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

            if go.Option().Index == idxs_Opt['EnableAllSrfGrips']:
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


def _findCoplanarGripsOfRhinoObject(rhObj, rdGrips_Ref, fPlanarityTol, bEnableGrips):
    rdGrips_Out = []

    if bEnableGrips and not rhObj.GripsOn:
        rhObj.GripsOn = True
        if not rhObj.GripsOn:
            #print("Could not turn on grips for brep.")
            return

    rdGrips_All_of_obj = rhObj.GetGrips()
    if rdGrips_All_of_obj is None:
        return 0

    iTol = 0
    fTol = Rhino.RhinoMath.ZeroTolerance
    bMaxTolReach = fPlanarityTol < fTol


    for grip in rdGrips_All_of_obj:

        if grip.IsSelected(checkSubObjects=False):
            sEval = "grip.Id"; print(sEval,'=',eval(sEval))
            grip.Select(True, syncHighlight=True) # For persistant selection.

    return
    while not bMaxTolReach:
        sc.escape_test()

        ref_pts_ToAdd = []

        fTol = Rhino.RhinoMath.ZeroTolerance*(10**(iTol))

        if fTol >= fPlanarityTol:
            fTol = fPlanarityTol
            bMaxTolReach = True

        for grip in rdGrips_All_of_obj:

            if grip.IsSelected(checkSubObjects=False):
                #grip.Select(True, syncHighlight=True) # For persistant selection.
                continue

            pt = grip.CurrentLocation

            if rg.Point3d.ArePointsCoplanar(ref_pts+[pt], tolerance=fTol):
                if grip.Select(True, syncHighlight=True):
                    iCt_added_grips += 1
                    ref_pts_ToAdd.append(pt)

        ref_pts.extend(ref_pts_ToAdd)
        iTol += 1

    return iCt_added_grips


def _selectCoplanarGripsOfRhinoObject(rhObj, ref_pts, fPlanarityTol, bEnableGrips):
    iCt_added_grips = 0

    if bEnableGrips and not rhObj.GripsOn:
        rhObj.GripsOn = True
        if not rhObj.GripsOn:
            #print("Could not turn on grips for brep."
            return 0

    grips = rhObj.GetGrips()
    if grips is None:
        return 0

    iTol = 0
    fTol = Rhino.RhinoMath.ZeroTolerance
    bMaxTolReach = fPlanarityTol < fTol


    for grip in grips:

        if grip.IsSelected(checkSubObjects=False):
            grip.Select(True, syncHighlight=True) # For persistant selection.


    while not bMaxTolReach:
        sc.escape_test()

        ref_pts_ToAdd = []

        fTol = Rhino.RhinoMath.ZeroTolerance*(10**(iTol))

        if fTol >= fPlanarityTol:
            fTol = fPlanarityTol
            bMaxTolReach = True

        for grip in grips:

            if grip.IsSelected(checkSubObjects=False):
                #grip.Select(True, syncHighlight=True) # For persistant selection.
                continue

            pt = grip.CurrentLocation

            if rg.Point3d.ArePointsCoplanar(ref_pts+[pt], tolerance=fTol):
                if grip.Select(True, syncHighlight=True):
                    iCt_added_grips += 1
                    ref_pts_ToAdd.append(pt)

        ref_pts.extend(ref_pts_ToAdd)
        iTol += 1

    return iCt_added_grips


def _printGripsAddedFeedback(iCt_added_grips, ref_grips):
    if iCt_added_grips == 0:
        print("No grips were added to selection of {}.".format(
            len(ref_grips)))
    else:
        print("{} grips were added to selection of {}.".format(
            iCt_added_grips, len(ref_grips)))


def selectGrips_planePerMultipleGrips(objref_Grips_Ref, fPlanarityTol=Rhino.RhinoMath.ZeroTolerance, bSearchOtherSrfs=False, bEnableGrips=True, bEcho=True, bDebug=False):
    """
    """

    rdGrips_Ref = []
    pts_Grips = []
    gObjs_Grip_owners = []

    for o in objref_Grips_Ref:
        rdGrip = o.Object()
        rdGrips_Ref.append(rdGrip)
        pts_Grips.append(rdGrip.CurrentLocation)
        gObjs_Grip_owners.append(rdGrip.OwnerId)

    #sEval = "gObjs_Grip_owners"; print(sEval,'=',eval(sEval))

    bSuccess, plane, sLog = _tryGetPlane(pts_Grips, fPlanarityTol)
    if not bSuccess:
        if bEcho and sLog: print(sLog)
        return


    iCt_added_grips = 0


    for gBrep in set(gObjs_Grip_owners):
        rhObj = sc.doc.Objects.FindId(gBrep)
        _findCoplanarGripsOfRhinoObject(
            rhObj,
            rdGrips_Ref=rdGrips_Ref,
            fPlanarityTol=fPlanarityTol,
            bEnableGrips=bEnableGrips)
        iCt_added_grips += _selectCoplanarGripsOfRhinoObject(
            rhObj,
            ref_pts=pts_Grips,
            fPlanarityTol=fPlanarityTol,
            bEnableGrips=bEnableGrips)


    if not bSearchOtherSrfs:
        if bEcho:
            _printGripsAddedFeedback(iCt_added_grips, rdGrips_Ref)
        return

    iter = rd.ObjectEnumeratorSettings()
    iter.LockedObjects = False

    for rdObj in sc.doc.Objects.GetObjectList(iter):
        if rdObj.ObjectType == rd.ObjectType.Brep:
            if rdObj.Id in gObjs_Grip_owners:
                continue
            iCt_added_grips += _selectCoplanarGripsOfRhinoObject(
                rdObj,
                ref_pts=pts_Grips,
                fPlanarityTol=fPlanarityTol,
                bEnableGrips=bEnableGrips)

    if bEcho:
        _printGripsAddedFeedback(iCt_added_grips, rdGrips_Ref)


def _getParallelToPlane(iParallelWithPlane):
    if Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'CPlaneXY':
        return sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()
    if Opts.listValues['iParallelWithPlane'][iParallelWithPlane] == 'CPlaneYZ':
        cplane = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()
        #origin = cplane.Origin
        #xaxis = cplane.YAxis
        #yaxis = cplane.ZAxis
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


def selectGrips_planePerOneGrip(objref_Grip, iParallelWithPlane=3, fPlanarityTol=Rhino.RhinoMath.ZeroTolerance, bSearchOtherSrfs=False, bEnableGrips=True, bEcho=True, bDebug=False):
    """
    """

    rdGrip = objref_Grip.Object()

    plane = _getParallelToPlane(Opts.values['iParallelWithPlane'])

    sc.doc.Objects.AddSurface(rg.PlaneSurface(plane))
    sc.doc.Views.Redraw()


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
        selectGrips_planePerOneGrip(
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