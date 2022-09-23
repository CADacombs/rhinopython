"""
Select the 2 outer control points of linear control point region.
    For surfaces, points must share a row or column.
Option:
    Ends
        Fixed: The 2 selected control points will not be moved with those between them.
        MoveToAverageLine: the 2 selected control points will be calculated for the line
        of projection and subsequently be projected onto it. 
    DistributeEvenly=Yes: All control points between will be evenly distributed between
        selected control points.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220922-23: Created, starting as a split from another script.

TODO:
    Properly process periodic and other closed curves.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bFixedEnds'; keys.append(key)
    values[key] = True
    names[key] = 'Ends'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'MoveToAverageLine', 'Fixed')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDistributeEvenly'; keys.append(key)
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

        #if key == 'fDistance':
        #    if cls.riOpts[key].CurrentValue <= 1e-9:
        #        print("Invalid input for Distance value.")
        #        cls.riOpts[key].CurrentValue = cls.values[key]
        #        return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get control point grips.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 outer control points of span to linearize")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Grip

    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bFixedEnds')
        addOption('bDistributeEvenly')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            grips = [o.Object() for o in objrefs]
            return grips

        go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
        go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _areNeighbors(gripA, gripB):
    for directionR, directionS in ((1,0), (-1,0), (0,1), (0,-1)):
        grip_Neighbor = gripA.NeighborGrip(directionR=1, directionS=0, directionT=0, wrap=True)
        if grip_Neighbor is None: continue
        if grip_Neighbor.Index == gripB.Index:
            return True
    return False


def _getIndicesBetweenGrips(gripA, gripB):
    """
    Returns: list(int(ordered indices))
    """

    for directionR, directionS in ((1,0), (-1,0), (0,1), (0,-1)):
        idxs = []
        grip_WIP = gripA
        while True:
            sc.escape_test()
            grip_WIP = grip_WIP.NeighborGrip(
                directionR=directionR,
                directionS=directionS,
                directionT=0,
                wrap=True)

            if grip_WIP is None:
                # End of open direction.
                break

            if grip_WIP.Index == gripB.Index:
                return idxs

            if grip_WIP.Index == gripA.Index:
                # This will occur when surface is closed.
                break

            idxs.append(grip_WIP.Index)


def distributeEvenly(grips_ToMove, ptA, ptB):

    bSomePtTrans = False

    fDists_From1stRef = []
    for grip_ToMove in grips_ToMove:
        fDists_From1stRef.append(
            ptA.DistanceTo(grip_ToMove.CurrentLocation))

    grips_ToMove_Sorted = []
    for dist, grip in sorted(zip(fDists_From1stRef, grips_ToMove)):
        grips_ToMove_Sorted.append(grip)

    for i in range(len(grips_ToMove_Sorted)):
        perunum = float(i+1)/float(len(grips_ToMove_Sorted)+1)
        pt_Target = ptA + (ptB - ptA) * perunum
        grip_ToMove = grips_ToMove_Sorted[i]
        if grip_ToMove.CurrentLocation.DistanceTo(pt_Target) > Rhino.RhinoMath.ZeroTolerance:
            grip_ToMove.CurrentLocation = pt_Target
            bSomePtTrans = True

    return bSomePtTrans


def linearizeGripsBetween2(gripA, gripB, **kwargs):
    """
    """


    if gripA.Index == gripB.Index:
        return False, "The 2 grips are identical."


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]


    bFixedEnds = getOpt('bFixedEnds')
    bDistributeEvenly = getOpt('bDistributeEvenly')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if gripA.OwnerId != gripB.OwnerId:
        return False, "Grips must be from the same object."


    rdObj = rs.coercerhinoobject(gripA.OwnerId)


    idxs_Between = _getIndicesBetweenGrips(gripA, gripB)
    if idxs_Between is None:
        return False, "The grips do not share a row or column."
    elif len(idxs_Between) == 0:
        return False, "The 2 grips are consecutive."

    grips_All = rdObj.GetGrips()

    bSomePtTrans = False

    if bFixedEnds:
        idxs_ToMove = idxs_Between
        line = rg.Line(gripA.CurrentLocation, gripB.CurrentLocation)
    else:
        idxs_ToMove = [gripA.Index] + idxs_Between + [gripB.Index]
        bSuccess, line = rg.Line.TryFitLineToPoints(
            [grips_All[i].CurrentLocation for i in idxs_ToMove])


    for iT in idxs_ToMove:
        grip_ToMove = grips_All[iT]
        pt_Target = line.ClosestPoint(
            grip_ToMove.CurrentLocation, limitToFiniteSegment=False)

        if grip_ToMove.CurrentLocation.DistanceTo(pt_Target) > Rhino.RhinoMath.ZeroTolerance:
            grip_ToMove.CurrentLocation = pt_Target
            bSomePtTrans = True


    if bDistributeEvenly:
        idxs_ToMove = idxs_Between
        grips_ToMove = [grips_All[i] for i in idxs_ToMove]
        rc = distributeEvenly(grips_ToMove, gripA.CurrentLocation, gripB.CurrentLocation)
        bSomePtTrans = bSomePtTrans or rc

    if not bSomePtTrans:
        return False, "Grips are already at target locations."

    rc = sc.doc.Objects.GripUpdate(rdObj, deleteOriginal=True)

    return rc, None


def main():

    grips = getInput()
    if grips is None: return

    if not Opts.values['bDebug']: sc.doc.Views.RedrawEnabled = False

    rc = linearizeGripsBetween2(grips[0], grips[1])
    if rc is None: return
    bSuccess, sLog = rc

    if not bSuccess:
        print(sLog)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()