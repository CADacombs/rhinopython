"""
"""

from __future__ import print_function

"""
211112: Created.
"""

import Rhino.DocObjects as rd
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fTolerance'; keys.append(key)
    values[key] = 2.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAllowWires'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAllowMatedEdges'; keys.append(key)
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

    key = 'bAddRefs'; keys.append(key)
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList):

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue <= 1e-9:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_EdgesToSplit():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edges to split")

    go.GeometryFilter = rd.ObjectType.EdgeFilter

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fTolerance')
        addOption('bEcho')
        addOption('bDebug')
        if Opts.values['bDebug']:
            addOption('bAddRefs')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()

            return (
                objrefs,
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fTolerance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Refs(objref_SrfToMod):
    """
    """

    sc.doc.Objects.UnselectAll()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select reference edges")

    go.GeometryFilter = rd.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.EdgeCurve |
            ri.Custom.GeometryAttributeFilter.BoundaryEdge)

        if Opts.values['bAllowWires']:
            go.GeometryAttributeFilter |= (
                ri.Custom.GeometryAttributeFilter.WireCurve)

        if Opts.values['bAllowMatedEdges']:
            go.GeometryAttributeFilter |= (
                ri.Custom.GeometryAttributeFilter.MatedEdge)


        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fTolerance')
        addOption('bAllowWires')
        addOption('bAllowMatedEdges')
        addOption('bEcho')
        addOption('bDebug')
        if Opts.values['bDebug']:
            addOption('bAddRefs')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return (
                objrefs,
                Opts.values['fTolerance'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                Opts.values['bAddRefs'],
                )

        if res == ri.GetResult.Number:
            key = 'fTolerance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def processObjRefs(objrefs_M, objrefs_R, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    bAddRefs = getOpt('bAddRefs')


    def groupObjrefsPerBrep(objrefs):
        gBs = []
        objrefs_PerB = []

        for objref in objrefs:
            if objref.ObjectId in gBs:
                objrefs_PerB[gBs.index(objref.ObjectId)].append(objref)
            else:
                gBs.append(objref.ObjectId)
                objrefs_PerB.append([objref])

        return objrefs_PerB


    objrefs_per_B = groupObjrefsPerBrep(objrefs_M)

    rgCs_R = [o.Curve() for o in objrefs_R]

    for objrefs_same_B in objrefs_per_B:

        edges_M_In = [o.Edge() for o in objrefs_same_B]
        idxs_Es = [e.EdgeIndex for e in edges_M_In]

        edges_M_In = [e for i,e in sorted(zip(idxs_Es, edges_M_In), reverse=True)]

        rgB_WIP = edges_M_In[0].Brep.DuplicateBrep()

        edge_count_Start = rgB_WIP.Edges.Count

        for edge in edges_M_In:

            ts = []

            for c_R in rgCs_R:
                for bSuccess, t in (edge.ClosestPoint(c_R.PointAtStart), edge.ClosestPoint(c_R.PointAtEnd)):
                    if not bSuccess: continue
                    pt_Closest = edge.PointAt(t)
                    if (
                        (pt_Closest.DistanceTo(edge.PointAtStart) > fTolerance) and
                        (pt_Closest.DistanceTo(edge.PointAtEnd) > fTolerance)
                    ):
                        ts.append(t)
                        #sc.doc.Objects.AddPoint(pt_Closest)


            rgB_WIP.Edges.SplitEdgeAtParameters(
                edgeIndex=edge.EdgeIndex,
                edgeParameters=ts)

        rgB_WIP.Compact()

        edge_count_End = rgB_WIP.Edges.Count

        if edge_count_End == edge_count_Start:
            print("No edges were split.")
            continue

        bReplaced = sc.doc.Objects.Replace(
            objectId=objrefs_same_B[0].ObjectId,
            brep=rgB_WIP)

        if bReplaced:
            print("Edge count increase: {}".format(edge_count_End-edge_count_Start))


def main():

    rc = getInput_EdgesToSplit()
    if rc is None: return

    (
        objrefs_M,
        bEcho,
        bDebug,
       ) = rc


    rc = getInput_Refs(objrefs_M)
    if rc is None: return

    (
        objrefs_R,
        fTolerance,
        bEcho,
        bDebug,
        bAddRefs,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBrep_Ret = processObjRefs(
        objrefs_M,
        objrefs_R,
        fTolerance=fTolerance,
        bEcho=bEcho,
        bDebug=bDebug,
        bAddRefs=bAddRefs,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()