"""
This script is an alternative to _FilletEdge, using RhinoCommon's
Brep.CreateFilletEdges with flexibility in tolerance.

Advantage:
    Not only can a start tolerance be specified, but when the main method,
    Brep.CreateFilletEdges, fails, it is recalled with smaller, then optionally,
    larger tolerances.

Limitations:
    All fillets are constant.
    Fillets are not editable, i.e, _FilletEdge _Edit .

To minimize variation from the starting tolerance, it is recommended to fillet
each chain of fillets separately.

Send any questions, comments, or script development service needs to @spb on
the McNeel Forums ( https://discourse.mcneel.com/ )
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190903: Created.
191209: Modified an option default value.
221105: Bug fix.  Refactored.
        Disabled bConstantRadius option until UI is updated to better accomodate it.
221108: Bug fix.
221111: Now filters interior edges in main function since GetObject.EnablePreSelect doesn't do this for BrepEdges(/subobjects?).
221112: Added option to increase tolerance when others fillet result fails with other tolerances.
        Now allows input of different size fillets.

TODO:
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum
from System import Guid
from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fRadius'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    #key = 'bConstantRadius'; keys.append(key)
    #values[key] = True
    #riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    #key = 'fRadius_End'; keys.append(key)
    #values[key] = 1.0
    #names[key] = 'EndRadius'
    #riOpts[key] = ri.Custom.OptionDouble(values[key])
    #stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iRailType'; keys.append(key)
    values[key] = 1
    listValues[key] = Enum.GetNames(rg.RailType)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_Start'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    names[key] = 'StartingTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol_Min'; keys.append(key)
    values[key] = 0.0001
    names[key] = 'MinTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bIncrTolOnFail'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_Max'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'MaxTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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

        if key == 'fRadius':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key]
                print("Invalid radius input ignored.")
                return

            value = 4.0*sc.doc.ModelAbsoluteTolerance
            if cls.riOpts[key].CurrentValue < value:
                cls.riOpts[key].CurrentValue = cls.values[key] = value
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                print("Radius input was too small and therefore adjusted.")
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key == 'fRadius_End':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key]
                print("Invalid radius input ignored.")
                return

            value = 4.0*sc.doc.ModelAbsoluteTolerance
            if cls.riOpts[key].CurrentValue < value:
                cls.riOpts[key].CurrentValue = cls.values[key] = value
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                print("Radius input was too small and therefore adjusted.")
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key == 'fTol_Start':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6
                print("Tolerance input was too small and therefore adjusted.")

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key == 'fTol_Min':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6
                print("MinTol input was too small and therefore adjusted.")

            if cls.riOpts[key].CurrentValue > cls.riOpts['fTol_Start'].CurrentValue:
                cls.riOpts[key].CurrentValue = cls.riOpts['fTol_Start'].CurrentValue
                print("MinTol input was too small and therefore adjusted.")

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key == 'fTol_Max':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6
                print("MaxTol input was too small and therefore adjusted.")

            if cls.riOpts[key].CurrentValue < cls.riOpts['fTol_Start'].CurrentValue:
                cls.riOpts[key].CurrentValue = cls.riOpts['fTol_Start'].CurrentValue
                print("MaxTol input was too small and therefore adjusted.")

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def getInput(bPrevBrepsArePresent):
    """
    Get BrepEdges with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edges to fillet")
    go.SetCommandPromptDefault("Confirm edges & radius")

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.MatedEdge |
            ri.Custom.GeometryAttributeFilter.EdgeCurve)

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go_Main on repeats of While loop.
    go.EnablePreSelect(True, ignoreUnacceptablePreselectedObjects=True) # This doesn't work for edges.  Maybe it doesn't work for subobjects.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(True, acceptZero=True)

    go.AcceptNothing(True)

    idxs_Opts = {}

    sc.doc.Objects.UnselectAll()

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        #Opts.names['fRadius'] = 'Radius' if Opts.values['bConstantRadius'] else 'StartRadius'
        addOption('fRadius')
        #addOption('bConstantRadius')
        #if not Opts.values['bConstantRadius']:
        #    addOption('fRadius_End')
        addOption('iRailType')
        addOption('fTol_Start')
        addOption('fTol_Min')
        addOption('bIncrTolOnFail')
        if Opts.values['bIncrTolOnFail']:
            addOption('fTol_Max')
        if bPrevBrepsArePresent:
            idxs_Opts['UsePrevInput'] = go.AddOption('UsePrevInput')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return []

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = abs(go.Number())
            Opts.setValue(key)
            continue

        # An option was selected.
        if go.Option().Index == idxs_Opts['UsePrevInput']:
            go.Dispose()
            return 'UsePrevInput'

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


class DrawRadiusDotsConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.prs = None
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode

    def PostDrawObjects(self, drawEventArgs):
        if not self.prs: return

        for pt, rad in self.prs:

            if rad == 0.0:
                continue

            rc = drawEventArgs.Display.DrawDot(
                worldPosition=pt,
                text="{:.{}f}".format(rad, sc.doc.ModelDistanceDisplayPrecision),
                dotColor=sc.doc.Layers.CurrentLayer.Color,
                textColor=Color.Black if sc.doc.Layers.CurrentLayer.Color != Color.Black else Color.White)


def formatDistance(fDistance, iPrecision=None):
    if iPrecision is None:
        iPrecision = sc.doc.ModelDistanceDisplayPrecision
    if fDistance is None:
        return "(No value provided)"
    if fDistance == 0.0:
        return "0"
    if fDistance < 0.01:
        return "{:.2e}".format(fDistance)
    return "{:.{}g}".format(fDistance, iPrecision)


def processBrep(rgBrep_In, idxs_rgEdges_In, fRadii_In, bConstantRadius=True, **kwargs):
    """
    Parameters:
        rgBrep_In
        idxs_rgEdges_In
        fRadii
        bConstantRadius
        fRadius_End
        iRailType
        fTol_Start
        fTol_Min
        fTol_Max
        bEcho
        bDebug


    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    #fRadius_End = getOpt('fRadius_End')
    iRailType = getOpt('iRailType')
    fTol_Start = getOpt('fTol_Start')
    fTol_Min = getOpt('fTol_Min')
    fTol_Max = getOpt('fTol_Max')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if fTol_Min > fTol_Start: fTol_Min = fTol_Start
    if fTol_Max < fTol_Start: fTol_Max = fTol_Start


    def filterInteriorEdges(idxs_rgEdges_In, fRadii_In, bEcho=True):
        idxs_rgEs_Out = []
        idxs_Ignored = []
        fRadii_Out = []
        for idxE, fRadius in zip(idxs_rgEdges_In, fRadii_In):
            if rgBrep_In.Edges[idxE].Valence == rg.EdgeAdjacency.Interior:
                idxs_rgEs_Out.append(idxE)
                fRadii_Out.append(fRadius)
            else:
                idxs_Ignored.append(idxE)
        if bEcho and idxs_Ignored:
            print("{} non-interior edges are ignored.".format(len(idxs_Ignored)))
        return idxs_rgEs_Out, fRadii_Out


    idxs_rgEdges_WIP, fRadii_WIP = filterInteriorEdges(idxs_rgEdges_In, fRadii_In, bEcho)


    def nakedBorderCount(rgBrep):
        necs = rgBrep.DuplicateNakedEdgeCurves(nakedOuter=True, nakedInner=False)
        if necs.Count == 0: return 0
        jcs = rg.Curve.JoinCurves(necs, joinTolerance=1e-7)
        for nec in necs: nec.Dispose()
        ct_Open = [jc.IsClosed for jc in jcs].count(False)
        if ct_Open:
            raise Exception("{} border curves are open after JoinCurves.".format(len(ct_Open)))
        len_jcs = len(jcs)
        for jc in jcs: jc.Dispose()
        return len_jcs

    ct_Naked_In = nakedBorderCount(rgBrep_In)


    def maxEdgeTol(rgBrep):
        return max([edge.Tolerance for edge in rgBrep.Edges])

    maxEdgeTol_In = maxEdgeTol(rgBrep_In)


    fTol_WIP = fTol_Start
    bDecreasingTol = True


    def getNextTolerance(fTol_In, fTol_Start, fTol_Min, fTol_Max, bEcho):

        if fTol_In == fTol_Min:
            if fTol_Max is None or fTol_Max == fTol_Start:
                if bEcho: print("Fillet could not be created.")
                return
            fTol_Out = fTol_Start * 2.0

            if fTol_Out > fTol_Max:
                return fTol_Max

            return fTol_Out

        if fTol_In <= fTol_Start:
            fTol_Out = fTol_In / 2.0

            if fTol_Out < fTol_Min:
                return fTol_Min

            return fTol_Out

        if fTol_In == fTol_Max:
            if bEcho: print("Fillet could not be created.")
            return

        fTol_Out = fTol_In * 2.0

        if fTol_Out > fTol_Max:
            return fTol_Max

        return fTol_Out


    while True:
        sc.escape_test()


        rgBs_Out = rg.Brep.CreateFilletEdges(
            brep=rgBrep_In,
            edgeIndices=idxs_rgEdges_WIP,
            startRadii=fRadii_WIP,
            endRadii=fRadii_WIP, # if bConstantRadius else TBD,
            blendType=rg.BlendType.Fillet,
            railType=Enum.ToObject(rg.RailType, iRailType),
            tolerance=fTol_WIP)

        if rgBs_Out.Count == 0:
            # rgBs_Out does not return None.

            print("At {} tolerance, CreateFilletEdges returned an empty array.".format(
                fTol_WIP))

            fTol_WIP = getNextTolerance(fTol_WIP, fTol_Start, fTol_Min, fTol_Max, bEcho)
            if fTol_WIP is None: return

            continue


        if rgBs_Out.Count != 1:
            raise Exception("{} breps resulted from CreateFilletEdges.".format(rgBs_Out.Count))

        rgB_Out = rgBs_Out[0]

        ct_Naked_Out = nakedBorderCount(rgB_Out)

        if ct_Naked_Out > ct_Naked_In:
            print("At {} tolerance, naked border count change: {} -> {}.".format(
                fTol_WIP,
                ct_Naked_In,
                ct_Naked_Out))

            fTol_WIP = getNextTolerance(fTol_WIP, fTol_Start, fTol_Min, fTol_Max, bEcho)
            if fTol_WIP is None: return

            rgB_Out.Dispose()

            continue

        # No increase in naked border count.

        print("Fillet created at {} tolerance.".format(fTol_WIP))

        maxEdgeTol_Out = maxEdgeTol(rgB_Out)
        if maxEdgeTol_Out == maxEdgeTol_In:
            print("No change in maximum edge tolerance, {}.".format(
                formatDistance(maxEdgeTol_Out)))
        else:
            print("Maximum edge tolerance change: {} -> {}.".format(
                formatDistance(maxEdgeTol_In),
                formatDistance(maxEdgeTol_Out)))
        if bEcho and fTol_WIP != fTol_Start:
            print("To obtain filleted brep, tolerance was {} to {}.".format(
                'increased' if fTol_WIP > fTol_Start else 'decreased',
                fTol_WIP))
        return rgB_Out


def processBrepObject(rhBrep_In, idxs_rgEdges, fRadii, bConstantRadius=True, **kwargs):
    """
    Parameters:
        rhBrep_In
        fRadius
        bConstantRadius
        fRadius_End
        iRailType
        fTol_Start
        fTol_Min
        fTol_Max
        bReplace
        bEcho
        bDebug
    Return on success: GUID of resultant brep.
    Return on fail: None
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    #fRadius_End = getOpt('fRadius_End')
    iRailType = getOpt('iRailType')
    fTol_Start = getOpt('fTol_Start')
    fTol_Min = getOpt('fTol_Min')
    fTol_Max = getOpt('fTol_Max')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rdB_In = None
    if isinstance(rhBrep_In, Guid):
        gB_In = rhBrep_In
        rdB_In = sc.doc.Objects.FindId(gB_In)
    elif isinstance(rhBrep_In, rd.BrepObject):
        rdB_In = rhBrep_In
        gB_In = rdB_In.Id
    if not rdB_In: return


    rgB_In = rdB_In.BrepGeometry

    if not rgB_In.IsValid:
        print("Brep {} is invalid.  Fix and rerun this script.".format(gB_In))
        return


    rgB_Res = processBrep(
        rgB_In,
        idxs_rgEdges,
        fRadii,
        bConstantRadius=bConstantRadius,
        #fRadius_End=fRadius_End,
        iRailType=iRailType,
        fTol_Start=fTol_Start,
        fTol_Min=fTol_Min,
        fTol_Max=fTol_Max,
        bEcho=bEcho,
        bDebug=bDebug,
        )
        
    if rgB_Res is None: return

    if bReplace:
        if sc.doc.Objects.Replace(objectId=rdB_In.Id, brep=rgB_Res):
            print("Replaced brep.")
            rgB_Res.Dispose()
            return gB_In
        else:
            if bEcho: print("Failed replacing brep.")
            rgB_Res.Dispose()
            return
    else:
        gB_Out = sc.doc.Objects.AddBrep(rgB_Res, rdB_In.Attributes)
        if gB_Out != gBrep1.Empty:
            if bEcho: print("Added brep.")
            rgB_Res.Dispose()
            return gB_Out
        else:
            if bEcho: print("Failed adding brep.  Is it invalid?")
            rgB_Res.Dispose()
            return


def main():

    gBs_In = []
    idx_rgEdges_PerBrep = []
    rads_for_idxEs_PerBrep = []


    prevSels_Present = []
    skey_prevSels = 'prevSels({})'.format(__file__)
    if skey_prevSels in sc.sticky:
        prevSels_Saved = sc.sticky[skey_prevSels]
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        iter.ObjectTypeFilter = rd.ObjectType.Brep
        gBs_Saved = [d.Id for d in sc.doc.Objects.GetObjectList(iter)]
        for gB, pts, rads in prevSels_Saved:
            if gB in  gBs_Saved:
                prevSels_Present.append((gB, pts, rads))


    skey_conduit = 'conduit({})'.format(__file__)
    if (skey_conduit in sc.sticky) and sc.sticky[skey_conduit]:
        conduit = sc.sticky[skey_conduit]
    else:
        conduit = DrawRadiusDotsConduit()
        sc.sticky[skey_conduit] = conduit

    conduit.Enabled = False # Turns off the conduit if left on from a previous execution of this script.
    sc.doc.Views.Redraw()



    def get_BE_FromInput(rhObj):
        if not hasattr(rhObj, '__iter__'):
            if not isinstance(rhObj, rd.ObjRef): return

            objref = rhObj

            gB = objref.ObjectId
            idxE = objref.Edge().EdgeIndex

            return gB, idxE

        if len(rhObj) != 2: return
        rhB, idxE = rhObj
        if not isinstance(idxE, int): return

        if isinstance(rhB, rd.BrepObject):
            rdB = rhB
            gB = rdB.Id
        elif isinstance(rhB, Guid):
            gB = rhB

        return gB, idxE


    def sortInputPerBrep(rhObjs_In):

        for rhObj in rhObjs_In:

            rc = get_BE_FromInput(rhObj)
            if not rc: return

            gB, idxE = rc

            if gB in gBs_In:
                if idxE in idx_rgEdges_PerBrep[gBs_In.index(gB)]:
                    # Modify radius value.
                    idxEperB = idx_rgEdges_PerBrep[gBs_In.index(gB)].index(idxE)
                    rads_for_idxEs_PerBrep[gBs_In.index(gB)][idxEperB] = fRadius
                    continue

                idx_rgEdges_PerBrep[gBs_In.index(gB)].append(idxE)
                rads_for_idxEs_PerBrep[gBs_In.index(gB)].append(fRadius)
                continue

            gBs_In.append(gB)
            idx_rgEdges_PerBrep.append([idxE])
            rads_for_idxEs_PerBrep.append([fRadius])


    def getEdgeIndex_MatchingMidPoint(rgBrep, midpt):
        for edge in rgBrep.Edges:
            ts = list(edge.DivideByCount(2, includeEnds=False))
            if not ts:
                print("Midpoint of edge[{}] not found.".format(edge.EdgeIndex))
                return
            if edge.PointAt(ts[0]).DistanceTo(midpt) < 0.1 * sc.doc.ModelAbsoluteTolerance:
                return edge.EdgeIndex


    def prepareDataForPreviewAndSticky():
        zipped = zip(gBs_In, idx_rgEdges_PerBrep, rads_for_idxEs_PerBrep)


        rgBs = []
        pts_Mids_All = [] # Won't include matching radius of 0.0.
        rads_All = [] # Won't include radius of 0.0.
        
        bprs = [] # list of tuples of (GUID, pts, radii) for saving input for future use of script.

        for gB_In, idxs_Es, rads_for_idxEs in zipped:
            rdB = sc.doc.Objects.FindId(gB_In)
            rgB = rdB.Geometry
            rgBs.append(rgB)
            pts_Mids_ThisB = [] # Won't include matching radius of 0.0.
            rads_ThisB = [] # Won't include radius of 0.0.

            for i, rad in enumerate(rads_for_idxEs):
                if rad == 0.0:
                    continue
                iE = idxs_Es[i]
                ts = list(rgB.Edges[iE].DivideByCount(2, includeEnds=False))
                if not ts:
                    print("Midpoint of edge[{}] not found.".format(iE))
                    return
                pt = rgB.Edges[iE].PointAt(ts[0])
                pts_Mids_ThisB.append(pt)
                rads_ThisB.append(rad)

            if not pts_Mids_ThisB: continue
            bprs.append((rdB.Id, pts_Mids_ThisB, rads_ThisB))
            pts_Mids_All.extend(pts_Mids_ThisB)
            rads_All.extend(rads_ThisB)

        sc.sticky[skey_prevSels] = bprs

        if not rgBs: return

        return pts_Mids_All, rads_All



    while True:
        rc = getInput(bool(prevSels_Present))

        if rc is None:
            conduit.Enabled = False
            del conduit
            del sc.sticky[skey_conduit]
            sc.sticky[skey_conduit] = None
            return

        if rc == 'UsePrevInput':
            sc.doc.Objects.UnselectAll()
            gBs_In = []
            idx_rgEdges_PerBrep = []
            rads_for_idxEs_PerBrep = []
            for gB, pts, rads in prevSels_Present:
                rgB = sc.doc.Objects.FindId(gB).BrepGeometry
                idx_rgEdges_ThisB = []
                rads_for_idxEs_ThisB = []
                for i, pt in enumerate(pts):
                    idxE = getEdgeIndex_MatchingMidPoint(rgB, pt)
                    if idxE is None: continue
                    idx_rgEdges_ThisB.append(idxE)
                    rads_for_idxEs_ThisB.append(rads[i])
                if idx_rgEdges_ThisB:
                    gBs_In.append(gB)
                    idx_rgEdges_PerBrep.append(idx_rgEdges_ThisB)
                    rads_for_idxEs_PerBrep.append(rads_for_idxEs_ThisB)

            rc = prepareDataForPreviewAndSticky()
            if not rc: continue
            pts_Mids_All, rads_All = rc

            conduit.prs = zip(pts_Mids_All, rads_All)

            conduit.Enabled = True
            sc.doc.Views.Redraw()

            continue


        fRadius = Opts.values['fRadius']
        #bConstantRadius = Opts.values['bConstantRadius']
        #fRadius_End = Opts.values['fRadius_End']
        iRailType = Opts.values['iRailType']
        fTol_Start = Opts.values['fTol_Start']
        fTol_Min = Opts.values['fTol_Min']
        fTol_Max = Opts.values['fTol_Max'] if Opts.values['bIncrTolOnFail'] else None
        bReplace = Opts.values['bReplace']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        sc.doc.Objects.UnselectAll()
        sc.doc.Views.Redraw()

        if not rc:
            conduit.Enabled = False
            del conduit
            break

        objrefs = rc

        sortInputPerBrep(objrefs)

        rc = prepareDataForPreviewAndSticky()
        if not rc: continue
        pts_Mids_All, rads_All = rc

        conduit.prs = zip(pts_Mids_All, rads_All)

        conduit.Enabled = True
        sc.doc.Views.Redraw()


    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBreps_Res_All = []

    zipped = zip(gBs_In, idx_rgEdges_PerBrep, rads_for_idxEs_PerBrep)

    for iB, (gB_In, idxs_Es, rads_for_idxEs) in enumerate(zipped):

        sCmdPrompt = "Processing brep {} of {}".format(iB+1, len(gBs_In))
        Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt)

        gB_Ret = processBrepObject(
            rhBrep_In=gB_In,
            idxs_rgEdges=idxs_Es,
            fRadii=rads_for_idxEs,
            bConstantRadius=True,
            #fRadius_End=fRadius_End,
            iRailType=iRailType,
            fTol_Start=fTol_Start,
            fTol_Min=fTol_Min,
            fTol_Max=fTol_Max,
            bReplace=bReplace,
            bEcho=bEcho,
            bDebug=bDebug,
            )

        if gB_Ret:
            gBreps_Res_All.append(gB_Ret)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()