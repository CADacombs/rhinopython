"""
1. _PointsOn or _SolidPtOn to object to modify its grips, e.g., control points.
2. Run this script and follow prompts.
    Do not include end points with first selection, grips to move.
Options:
    DistributeEvenly=Yes: All control points will be evenly distributed between
        reference points.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

from __future__ import absolute_import, division, print_function, unicode_literals

#! python 2

"""
190904: Created.
220824: Now multiple grips will be made inline to 2 reference grips.
220922: Added DistributeEvenly option.  Refactored.
220923: Now, references are picked points, not other grips.
240829-30: Added OuterGripsAtRefPts option for when DistributeEvenly and BetweenRefPts are enabled.

TODO:
    Allow points to be scaled so each of the 2 extents are on a reference point.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bDistributeEvenly'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bBetweenRefPts'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOuterGripsAtRefPts'; keys.append(key)
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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_Grips_ToMove():
    
    #
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select grips to move")
    
    go.GeometryFilter = rd.ObjectType.Grip
    
    # go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False
    go.EnableClearObjectsOnEntry(False)
    go.EnableUnselectObjectsOnExit(False)

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)
    
    idxs_Opt = {}

    bFirstRun = True

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bDistributeEvenly')
        if Opts.values['bDistributeEvenly']:
            addOption('bBetweenRefPts')
            if Opts.values['bBetweenRefPts']:
                addOption('bOuterGripsAtRefPts')
        addOption('bEcho')
        addOption('bDebug')

        sc.doc.Views.Redraw()

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            grips_ToMove = [o.Object() for o in objrefs]
            return grips_ToMove

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break

        if bFirstRun:
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True) # Default is True.
            bFirstRun = False


def getInput_Points_Ref():
    
    pts = [None, None]

    def addOption(key): idxs_Opt[key] = Opts.addOption(gp, key)

    for i, s in zip((0,1),('1st', '2nd')):
        gp = ri.Custom.GetPoint()

        gp.SetCommandPrompt("Pick {} reference point".format(s))
    

        idxs_Opt = {}

        while True:
            gp.ClearCommandOptions()

            idxs_Opt.clear()

            addOption('bDistributeEvenly')
            if Opts.values['bDistributeEvenly']:
                addOption('bBetweenRefPts')
                if Opts.values['bBetweenRefPts']:
                    addOption('bOuterGripsAtRefPts')
            addOption('bEcho')
            addOption('bDebug')

            res = gp.Get()

            if res == ri.GetResult.Cancel:
                gp.Dispose()
                return

            if res == ri.GetResult.Point:
                pts[i] = gp.Point()
                break

            for key in idxs_Opt:
                if gp.Option().Index == idxs_Opt[key]:
                    Opts.setValue(key, gp.Option().CurrentListOptionIndex)
                    break

    return pts


def processGrips(grips_ToMove, pts_Ref, bDistributeEvenly=False, bBetweenRefPts=True, bOuterGripsAtRefPts=False, bEcho=True, bDebug=False):
    """
    """

    ptA, ptB = pts_Ref

    line = rg.Line(ptA, ptB)

    rdOwners = []
    gOwners = []

    bSomeGripsWereTranslated = False

    for grip_ToMove in grips_ToMove:

        pt_Target = line.ClosestPoint(
                grip_ToMove.CurrentLocation,
                limitToFiniteSegment=False)
    
        gOwner = grip_ToMove.OwnerId
        rdOwner = sc.doc.Objects.FindId(gOwner) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gOwner)
        if rdOwner.Id not in gOwners:
            rdOwners.append(rdOwner)
            gOwners.append(rdOwner.Id)

        if grip_ToMove.CurrentLocation.DistanceTo(pt_Target) > Rhino.RhinoMath.ZeroTolerance:
            grip_ToMove.CurrentLocation = pt_Target
            bSomeGripsWereTranslated = True


    if not bDistributeEvenly:
        if not bSomeGripsWereTranslated:
            print("Grips are already at target locations. No objects were modified.")
            return

        for rdOwner in rdOwners:
            sc.doc.Objects.GripUpdate(rdOwner, deleteOriginal=True)

        return



    # Distribute grips evenly.


    def sortGripsFromPt(grips_In, pt):
        fDists_From1stRef = []
        for grip_ToMove in grips_In:
            fDists_From1stRef.append(
                pt.DistanceTo(grip_ToMove.CurrentLocation))

        grips_ToMove_Sorted = []
        for dist, grip in sorted(zip(fDists_From1stRef, grips_In)):
            grips_ToMove_Sorted.append(grip)

        return grips_ToMove_Sorted


    def evenlyRedistributeSortedGrips(grips, ptStart, ptEnd, bOuterGripsAtExtents):
        bSomeGripsWereTranslated = False
        if bOuterGripsAtExtents:
            for i in range(len(grips)):
                perunum = float(i) / float(len(grips) - 1)
                pt_Target = ptStart + (ptEnd - ptStart) * perunum
                grip_ToMove = grips_ToMove_Sorted[i]
                if grip_ToMove.CurrentLocation.DistanceTo(pt_Target) > Rhino.RhinoMath.ZeroTolerance:
                    grip_ToMove.CurrentLocation = pt_Target
                    bSomeGripsWereTranslated = True
        else:
             for i in range(len(grips)):
                perunum = float(i + 1) / float(len(grips) + 1)
                pt_Target = ptStart + (ptEnd - ptStart) * perunum
                grip_ToMove = grips[i]
                if grip_ToMove.CurrentLocation.DistanceTo(pt_Target) > Rhino.RhinoMath.ZeroTolerance:
                    grip_ToMove.CurrentLocation = pt_Target
                    bSomeGripsWereTranslated = True

        return bSomeGripsWereTranslated


    if bBetweenRefPts:
        grips_ToMove_Sorted = sortGripsFromPt(grips_ToMove, ptA)

        bSomeGripsWereTranslated |= evenlyRedistributeSortedGrips(
            grips_ToMove_Sorted,
            ptA,
            ptB,
            bOuterGripsAtExtents=bOuterGripsAtRefPts
        )
    else:
        # Do not distribute between the 2 reference points.
        # The 2 points extents (not necessarily the reference points) remain
        # and the rest are distributed between them.
        bSuccess, line = rg.Line.TryFitLineToPoints(
            [grip.CurrentLocation for grip in grips_ToMove])
        if not bSuccess:
            print("Line could not be generated from transformed grip locations.")
        else:
            grips_ToMove_Sorted = sortGripsFromPt(grips_ToMove, line.From)

            bSomeGripsWereTranslated |= evenlyRedistributeSortedGrips(
                grips_ToMove_Sorted,
                line.From,
                line.To,
                bOuterGripsAtExtents=True
            )
            for i in range(len(grips_ToMove_Sorted)):
                perunum = float(i) / float(len(grips_ToMove_Sorted) - 1)
                pt_Target = line.From + (line.To - line.From) * perunum
                grip_ToMove = grips_ToMove_Sorted[i]
                if grip_ToMove.CurrentLocation.DistanceTo(pt_Target) > Rhino.RhinoMath.ZeroTolerance:
                    grip_ToMove.CurrentLocation = pt_Target
                    bSomeGripsWereTranslated = True

            # print(sc.doc.Objects.AddLine(line))

    if not bSomeGripsWereTranslated:
        print("Grips are already at target locations. No objects were modified.")
        return

    for rdOwner in rdOwners:
        sc.doc.Objects.GripUpdate(rdOwner, deleteOriginal=True)


def main():

    grips_ToMove = getInput_Grips_ToMove()
    if grips_ToMove is None: return

    bDebug = Opts.values['bDebug']

    pts_Ref = getInput_Points_Ref()
    if pts_Ref is None: return

    bDistributeEvenly = Opts.values['bDistributeEvenly']
    bBetweenRefPts = Opts.values['bBetweenRefPts']
    bOuterGripsAtRefPts = Opts.values['bOuterGripsAtRefPts']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    processGrips(
        grips_ToMove,
        pts_Ref,
        bDistributeEvenly=bDistributeEvenly,
        bBetweenRefPts=bBetweenRefPts,
        bOuterGripsAtRefPts=bOuterGripsAtRefPts,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()