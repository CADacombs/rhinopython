"""
1. _PointsOn or _SolidPtOn to object to modify its grips, e.g., control points.
2. Run this script and follow prompts.
    Do not include end points with first selection, grips to move.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
190904: Created.
220824: Now multiple grips will be made inline to 2 reference grips.
220922: Added bDistributeEvenly.  Refactored.
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
    
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)
    
    idxs_Opt = {}
    
    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bDistributeEvenly')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            grips_ToMove = [o.Object() for o in objrefs]
            return grips_ToMove
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Grips_Ref():
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select 2 reference grips")
    
    go.GeometryFilter = rd.ObjectType.Grip
    
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)
    
    idxs_Opt = {}
    
    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bDistributeEvenly')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)
        
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            grips_Ref = [objref.Object() for objref in objrefs]
            return grips_Ref
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def processGrips(grips_ToMove, grips_Ref, bDistributeEvenly=False, bEcho=True, bDebug=False):
    """
    """

    ptA = grips_Ref[0].CurrentLocation
    ptB = grips_Ref[1].CurrentLocation

    line = rg.Line(ptA, ptB)

    rdOwners = []
    gOwners = []

    bSomePtTrans = False

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
            bSomePtTrans = True

    if bDistributeEvenly:
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


    if not bSomePtTrans:
        print("Grips are already at target locations.")
        return

    for rdOwner in rdOwners:
        sc.doc.Objects.GripUpdate(rdOwner, deleteOriginal=True)


def main():

    grips_ToMove = getInput_Grips_ToMove()
    if grips_ToMove is None: return

    bDistributeEvenly = Opts.values['bDistributeEvenly']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    grips_Ref = getInput_Grips_Ref()
    if grips_Ref is None: return

    bDistributeEvenly = Opts.values['bDistributeEvenly']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug:
        sc.doc.Views.RedrawEnabled = False

    processGrips(
        grips_ToMove,
        grips_Ref,
        bDistributeEvenly=bDistributeEvenly,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()