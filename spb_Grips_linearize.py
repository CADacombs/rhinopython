"""
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
181012-13: Created.
201002: Added support for all other point input combinations.  Previously, only 2 special cases were supported.
211112: With input of 2 CP's with a common u or v index, the CP's between will be made colinear.

TODO:
    Properly process periodic and other closed curves.
"""

import Rhino
import Rhino.DocObjects as rd
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

        if key == 'fDistance':
            if cls.riOpts[key].CurrentValue <= 1e-9:
                print("Invalid input for Distance value.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """ Get control point grips. """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select CPs")

    go.GeometryFilter = rd.ObjectType.Grip

    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bFixedEnds')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            grips = [o.Object() for o in objrefs]
            return grips

            #object_list = Rhino.Collections.TransformObjectList()
            #object_list.AddObjects(go=go, allowGrips=True)

            #gips_All = [] # List of (Tuples of Guid of grip's parent, index of grip, Point3d).  Used by rhinoscriptsyntax.
            #for grip in object_list.GripArray():
            #    gips_All.append((grip.OwnerId, grip.Index, grip.CurrentLocation))
            #return gips_All

        go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
        go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def sortNestListOfTuples(gips_to_sort, idxTuple):
    """
    Parameters:
        gips_to_sort: list(tuple(Guid of grip's parent, index of grip, Point3d)).  Used by rhinoscriptsyntax.
        idxTuple: Index of tuple by which to sort.
    """
    
    idxList = []
    nestedList = [] # 1-level nesting of (GUID, Grip index, Point3d) per GUID.

    for gip in gips_to_sort:
        if gip[0] in idxList:
            nestedList[idxList.index(gip[idxTuple])].extend([(gip[0], gip[1], gip[2])])
        else:
            idxList.append(gip[idxTuple])
            nestedList.append([(gip[0], gip[1], gip[2])])
    
    return nestedList


def _sortGripsByParentGuid(grips):
    """
    """

    nestedList = [] # 1-level nesting of grips per parent GUID.
    gParents = []

    for grip in grips:
        if grip.OwnerId in gParents:
            nestedList[gParents.index(grip.OwnerId)].append(grip)
        else:
            gParents.append(grip.OwnerId)
            nestedList.append([grip])
    
    return nestedList


def getUvIdxFromNsPoint1dIdx(ns, idxFlat):
    return idxFlat // ns.Points.CountV, idxFlat % ns.Points.CountV


def _get_iGrip_matches(grips_perOwnerGuid):
    
    pt_Matches = []
    iGrip_per_Guid_matches = []
    
    # Find matches between first two GUID's.
    for i0 in range(len(grips_perOwnerGuid[0])):
        for i1 in range(len(grips_perOwnerGuid[1])):
            dist = grips_perOwnerGuid[0][i0].CurrentLocation.DistanceTo(
                grips_perOwnerGuid[1][i1].CurrentLocation)
            if dist <= Rhino.RhinoMath.ZeroTolerance:
                pt_Matches.append(grips_perOwnerGuid[0][i0].CurrentLocation)
                iGrip_per_Guid_matches.append([i0, i1])
    if not pt_Matches:
        if bDebug: print("No matches found.")
        return
        
    for iMatch, pt_Match in enumerate(pt_Matches):
        for iG in range(2, len(grips_perOwnerGuid)):
            for i_gip, gip in enumerate(grips_perOwnerGuid[iG]):
                if gip[2].EpsilonEquals(pt_Match, Rhino.RhinoMath.ZeroTolerance):
                    iGrip_per_Guid_matches[iMatch].append(i_gip)
                    break
            else:
                if bDebug: print("No matches found.")
                break # to pt_Matches loop.
        if len(iGrip_per_Guid_matches[iMatch]) == len(grips_perOwnerGuid):
            return tuple(iGrip_per_Guid_matches[iMatch])


def process2Objects_2GripsEach(grips_perOwnerGuid):
    """

    """

    #            for iA in 0,1:
    #                for iB in 0,1:
    #                    if gips_perGuid[0][iA][2] == gips_perGuid[1][iB][2]:
    #                        return iA, iB
    i_gip_match = _get_iGrip_matches(grips_perOwnerGuid)
    if not i_gip_match:
        print("No match")
        return
        
    ptMatch = gips_perGuid[0][i_gip_match[0]][2]
        
    idxOtherA = int(not i_gip_match[0])
    ptToMoveA = gips_perGuid[0][idxOtherA][2]
    idxOtherB = int(not i_gip_match[1])
    ptOtherB = gips_perGuid[1][idxOtherB][2]
        
    line = rg.Line(ptToMoveA, ptOtherB)
    #sc.doc.Objects.AddLine(line)
        
    pt_on_line_closest_to_match = line.ClosestPoint(testPoint=ptMatch, limitToFiniteSegment=True)
    #sc.doc.Objects.AddPoint(pt_on_line_closest_to_match)
        
    vect = ptMatch - pt_on_line_closest_to_match
    if vect.IsTiny():
        if bEcho: print("Grips are already linear.")
        return
        
    xform_moveLine = rg.Transform.Translation(vect)
    ptToMoveA.Transform(xform_moveLine)
    ptOtherB.Transform(xform_moveLine)
        
    print(gips_perGuid[0][i_gip_match[0]][1])
        
    print(rs.ObjectGripLocation(gips_perGuid[0][0][0], gips_perGuid[0][idxOtherA][1], ptToMoveA))
    rs.ObjectGripLocation(gips_perGuid[1][0][0], gips_perGuid[1][idxOtherB][1], ptOtherB)


def processGrips(grips, **kwargs):
    """bDebug
    # gip: Tuple of Guid of grip's parent, index of grip, Point3d.
    """

    grips_In = grips # So 'grips' can be used in other routines.

    for grip in grips_In:
        if not isinstance(grip, rd.GripObject):
            print("{} object passed to processGrips.  Must be GripObject.")
            return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bFixedEnds = getOpt('bFixedEnds')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')




    def minimumGipDistance():
        dists = []
        for iA in range(len(gips_perGuid)):
            for gipA in gips_perGuid[iA]:
                for iB in range(iA+1, len(gips_perGuid)):
                    for gipB in gips_perGuid[iB]:
                        dist = gipA[2].DistanceTo(gipB[2])
                        if dist > sc.doc.ModelAbsoluteTolerance:
                            dists.append(dist)
        return min(dists)


    def process8Of4():
        """
        Corner of 4?
        """
        
        i_gip_match = _get_iGrip_matches(gips_perGuid)
        if not i_gip_match:
            print("No match.")
            return
        
        ptMatch = gips_perGuid[0][i_gip_match[0]][2]
        
        gip_ToMoveA = None
        for i in range(len(gips_perGuid)):
            if len(gips_perGuid[i]) == 2:
                if gip_ToMoveA is None:
                    gip_ToMoveA = gips_perGuid[int(not i_gip_match[0])]
                    ptToMoveA = gip_ToMoveA[2]
                else:
                    gip_ToMoveB = gips_perGuid[int(not i_gip_match[0])]
                    ptToMoveB = gip_ToMoveB[2]
        
        line = rg.Line(ptToMoveA, ptOtherB)
        #sc.doc.Objects.AddLine(line)
        
        pt_on_line_closest_to_match = line.ClosestPoint(
                testPoint=ptMatch, limitToFiniteSegment=True)
        #sc.doc.Objects.AddPoint(pt_on_line_closest_to_match)
        
        vect = ptMatch - pt_on_line_closest_to_match
        if vect.IsTiny():
            if bEcho: print("Grips are already linear.")
            return
        
        xform_moveLine = rg.Transform.Translation(vect)
        ptToMoveA.Transform(xform_moveLine)
        ptOtherB.Transform(xform_moveLine)
        
        print(gips_perGuid[0][i_gip_match[0]][1])
        
        print(rs.ObjectGripLocation(gips_perGuid[0][0][0], gips_perGuid[0][idxOtherA][1], ptToMoveA))
        rs.ObjectGripLocation(gips_perGuid[1][0][0], gips_perGuid[1][idxOtherB][1], ptOtherB)


    def processOthers():
        
        #pts = [gips_All[i][2] for i in range(len(gips_All))]

        if bFixedEnds:
            pass


            ## Look for co-U or c-V pattern.
            #for gObj, ips in zip(gObjs, ips_per_Guid):
            #    nc = rs.coercecurve(gObj)
            #    face = rs.coercesurface(gObj)
            #    if nc:
            #        pass
            #    elif face:
            #        ns = face.UnderlyingSurface()
            #        iUs = []
            #        iVs = []
            #        for idxPt, pt in ips:
            #            iU,iV = getUvIdxFromNsPoint1dIdx(ns, idxPt)
            #            iUs.append(iU)
            #            iVs.append(iV)
            #        print(iUs, iVs
            #        iUs_toLinearize = []
            #        iVs_toLinearize = []
            #        for iU in set(iUs):
            #            if iUs.count(iU) > 2:
            #                iUs_toLinearize.append(iU)
            #        for iV in set(iVs):
            #            if iVs.count(iV) > 2:
            #                iVs_toLinearize.append(iV)

            #        print(iUs_toLinearize, iVs_toLinearize
        
        bSuccess, line = rg.Line.TryFitLineToPoints(
            [gips_All[i][2] for i in range(len(gips_All))])
        
        #sc.doc.Objects.AddLine(line); #sc.doc.Views.Redraw()
        
        for i in range(len(gips_All)):
            pt_projected_onto_line = line.ClosestPoint(
                    testPoint=gips_All[i][2],
                    limitToFiniteSegment=False)
            rs.ObjectGripLocation(
                gips_All[i][0],
                gips_All[i][1],
                pt_projected_onto_line)
        
        
    
    # Split gips_All into a list of (GUID, Grip index, Point3d) based on GUID order
    # and a list of the sorted GUIDs.
    #gips_perGuid = sortNestListOfTuples(gips_All, idxTuple=0)
    #if gips_perGuid is None: return

    grips_perParentGuid = _sortGripsByParentGuid(grips_In)
    if grips_perParentGuid is None: return

    if (len(grips_In) == 2) and (grips_In[0].OwnerId == grips_In[1].OwnerId):
        import spb_Grips_linearizeBetween2
        spb_Grips_linearizeBetween2.linearizeGripsBetween2(
            grips_In[0],
            grips_In[1],
            bFixedEnds=bFixedEnds,
            bDistributeEvenly=False,
            bEcho=bEcho,
            bDebug=bDebug,
            )
        return

    # Reject if not at least 3 grips.
    if len(grips_In) < 3:
        if bEcho: print("Only {} grips selected.  Nothing will be done.".format(len(grips_In)))
        return

    # If there are 2 GUIDs with 2 grips each, they may span over edge.
    if (
        len(grips_perParentGuid) == 2 and
        len(grips_perParentGuid[0]) == 2 and
        len(grips_perParentGuid[1]) == 2
    ):
        process2Objects_2GripsEach(grips_perParentGuid)
    elif (
            len(gips_perGuid) == 4 and
            all(len(gips_of_Guid) == 2 for gips_of_Guid in gips_perGuid)
    ):
       process8Of4()
    else:
        processOthers()


def main():

    grips = getInput()
    if grips is None: return

    if not Opts.values['bDebug']: sc.doc.Views.RedrawEnabled = False

    processGrips(grips)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()