"""
Knot removal methods will not deform curve when there is adequate parametric continuity
through the knots.
When removing knots one by one and the knot remains, Location of it before and after knot removal
can be combined with curve deviation to test whether knot should be removed.
Unfortunately, location is of surface parameter space,
and thus determining a tolerance is not straightforward.
Possibly, tolerance should be RhinoMath.ZeroTolerance or a fraction of a domain length of the surface.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221024-26: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAllowBrepSelection'; keys.append(key)
    values[key] = False
    names[key] = 'SelectionMode'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'EdgesOnly', 'EdgesandBreps')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDecimalTol'; keys.append(key)
    values[key] = False
    names[key] = 'ToleranceType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'AbsValueOfSrfParamSpace', 'DecimalOfTrimCrvDomain')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_Absolute'; keys.append(key)
    values[key] = 1e-9
    names[key] = 'Tol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol_Decimal'; keys.append(key)
    values[key] = 1e-9
    names[key] = 'Decimal'
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

        if key == 'fTol_Absolute':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'fTol_Decimal':
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


def getInput():
    """
    Get breps and/or edges with optional input.
    """

    go = ri.Custom.GetObject()

    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.EdgeCurve

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        if Opts.values['bAllowBrepSelection']:
            go.SetCommandPrompt("Select breps and/or edges")
            go.SetCommandPromptDefault("All normal breps when none are selected")
            go.GeometryFilter = rd.ObjectType.Brep | rd.ObjectType.Curve
        else:
            go.SetCommandPrompt("Select edges")
            go.GeometryFilter = rd.ObjectType.Curve


        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('bAllowBrepSelection')
        addOption('bDecimalTol')
        if Opts.values['bDecimalTol']:
            addOption('fTol_Decimal')
        else:
            addOption('fTol_Absolute')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Nothing:
            #print(len(go.Objects()))
            oes = rd.ObjectEnumeratorSettings()
            oes.NormalObjects = True
            oes.LockedObjects = False
            oes.IncludeLights = False
            oes.IncludeGrips = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep

            rdBrepObjs = list(sc.doc.Objects.GetObjectList(oes))
            go.Dispose()
            if len(rdBrepObjs) == 0: return

            return rdBrepObjs

        if res == ri.GetResult.Number:
            if Opts.values['bDecimalTol']:
                key = 'fTol_Decimal'
            else:
                key = 'fTol_Absolute'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def knotMultiplicities(nc):
    mp = []
    iK = 0
    while iK < nc.Knots.Count:
        m = nc.Knots.KnotMultiplicity(iK)
        mp.append(m)
        iK += m
    return mp


def getMaxDev(rgCrvA, rgCrvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            rgCrvA,
            rgCrvB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
    if rc[0]:
        return rc[1]


def removeKnots(nc_In, bDecimalTol, fTol, bEcho=True, bDebug=False):
    """
    """


    if not isinstance(nc_In, rg.NurbsCurve):
        return

    if nc_In.SpanCount == 1:
        return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if bDecimalTol is None: bDecimalTol = Opts.values['bDecimalTol']
    if fTol is None:
        if bDecimalTol:
            fTol = getOpt('fTol_Decimal')
        else:
            fTol = getOpt('fTol_Absolute')


    fTol_Actual = max((1e-9, fTol/nc_In.Domain.Length)) if bDecimalTol else fTol
    #print(fTol_Actual)

    nc_Res = nc_In.DuplicateCurve()


    if bDebug:
        sEval = "knotMultiplicities(nc_Res)"; print("{}: {}".format(sEval, eval(sEval)))


    iK = nc_Res.Knots.Count - nc_Res.Degree - 1

    nc_WIP = nc_Res.Duplicate()
    nc_Success = None


    while iK >= nc_Res.Degree:
        sc.escape_test()

        if nc_Success is None:
            nc_WIP = nc_Res.Duplicate()
        else:
            nc_WIP = nc_Success.Duplicate()

        m = nc_WIP.Knots.KnotMultiplicity(iK)
        if bDebug: sEval = "m"; print("{}: {}".format(sEval, eval(sEval)))

        for i in range(m):
            iK_Start = iK - m + i + 1
            if bDebug: sEval = "iK_Start"; print("{}: {}".format(sEval, eval(sEval)))
            iK_End = iK + 1
            if bDebug: sEval = "iK_End"; print("{}: {}".format(sEval, eval(sEval)))

            bKnotsRemoved = nc_WIP.Knots.RemoveKnots(iK_Start, iK_End)
            if bDebug: sEval = "bKnotsRemoved"; print("{}: {}".format(sEval, eval(sEval)))

            if not bKnotsRemoved:
                continue

            if getMaxDev(nc_WIP, nc_In) > fTol_Actual:
                nc_WIP.Dispose()
                if nc_Success is None:
                    nc_WIP = nc_Res.Duplicate()
                else:
                    nc_WIP = nc_Success.Duplicate()
                continue

            # Success.
            if nc_Success is not None:
                nc_Success.Dispose()
            nc_Success = nc_WIP
            nc_WIP = None
            m = nc_Success.Knots.KnotMultiplicity(iK)
            if bDebug: sEval = "m"; print("{}: {}".format(sEval, eval(sEval)))
            break

        if nc_Success is None:
            pass

        iK -= m

    return nc_Success


def processBrep(rgBrep, idx_Edges=None, bDecimalTol=None, fTol=None, bEcho=True, bDebug=False):
    """
    rgBrep: Input and modified.
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if bDecimalTol is None: bDecimalTol = Opts.values['bDecimalTol']
    if fTol is None:
        if bDecimalTol:
            fTol = getOpt('fTol_Decimal')
        else:
            fTol = getOpt('fTol_Absolute')


    iTrims_Modified = []
    iFaces_Modified = []

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt


    if idx_Edges is None:
        idxTrims_ToMod = range(rgBrep.Trims.Count)
    else:
        idxTrims_ToMod = []
        for iE in idx_Edges:
            for iT in rgBrep.Edges[iE].TrimIndices():
                idxTrims_ToMod.append(iT)


    for iT in idxTrims_ToMod:
        Rhino.RhinoApp.CommandPrompt = sCmdPrompt_In + "  Processing at {} of {} BrepTrims in this Brep.".format(
            iT+1,
            len(idxTrims_ToMod))

        rgT = rgBrep.Trims[iT]
        rgC_In = rgT.TrimCurve

        #print(rgC_In.Domain.Length)

        nc_Res = removeKnots(
            rgC_In,
            bDecimalTol=bDecimalTol,
            fTol=fTol,
            bEcho=bEcho,
            bDebug=bDebug,
            )
        if nc_Res is None: continue

        i_New = rgBrep.Curves2D.Add(nc_Res)
        bSetTrimCurve = rgT.SetTrimCurve(i_New)

        if bSetTrimCurve:
            iTrims_Modified.append(iT)
            iFaces_Modified.append(rgT.Face.FaceIndex)

    iFaces_Modified = sorted(set(iFaces_Modified))


    #sEval = "rgBrep_In.Curves2D.Count"; print("{}: {}".format(sEval, eval(sEval)))
    if bDebug:
        sEval = "len(iTrims_Modified)"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "len(iFaces_Modified)"; print("{}: {}".format(sEval, eval(sEval)))

    if len(iTrims_Modified) == 0: return

    for iF in iFaces_Modified:
        face = rgBrep.Faces[iF]
        face.RebuildEdges(
            tolerance=sc.doc.ModelAbsoluteTolerance,
            rebuildSharedEdges=True,
            rebuildVertices=True)

    rgBrep.Compact()
    bIsValid, sLog = rgBrep.IsValidWithLog()
    if not bIsValid:
        rgBrep.Repair(tolerance=sc.doc.ModelAbsoluteTolerance)
        bIsValid = rgBrep.IsValid
        if bIsValid and bEcho:
            print("Error(s) fixed by Brep.Repair: {}".format(sLog))
        else:
            print("Error(s) that was not fixed by Brep.Repair: {}".format(sLog))

    #sEval = "rgBrep_In.Curves2D.Count"; print("{}: {}".format(sEval, eval(sEval)))

    return True


def createSortedBrepsAndEdges_FromObjRefs(objrefs):
    """
    Parameters:
        objrefs of Edges
    Returns on success:
        list(GUIDs of breps) and list(list(Edge indices) per brep)
    """

    rdBreps = []
    gBreps = []
    idxs_Edges = []

   # Organize input into synchronized lists of brep ids and edge indices.
    for objref in objrefs:
        rdBrep = objref.Object()
        gBrep = objref.ObjectId
        idx_Edge = objref.Edge().EdgeIndex
        
        if gBrep not in gBreps:
            # Brep not in list, so add it as well as objref's edge.
            rdBreps.append(rdBrep)
            gBreps.append(gBrep)
            idxs_Edges.append([idx_Edge])
        else:
            # Brep already in list, so add the objref's edge to the relative index.
            iB = gBreps.index(gBrep)
            if idx_Edge not in idxs_Edges[iB]:
                idxs_Edges[iB].append(idx_Edge)

    return rdBreps, idxs_Edges


def processBrepObjects(rhObjs_ToModify, bDecimalTol=None, fTol=None, bEcho=True, bDebug=False):
    """
    Parameters:
        rhObjs_In: (rd.ObjRef of BrepEdge) or rd.BrepObject
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if bDecimalTol is None: bDecimalTol = Opts.values['bDecimalTol']
    if fTol is None:
        if bDecimalTol:
            fTol = getOpt('fTol_Decimal')
        else:
            fTol = getOpt('fTol_Absolute')


    if isinstance(rhObjs_ToModify[0], rd.BrepObject):
        rdBs = rhObjs_ToModify
        idx_Edges = None
    elif isinstance(rhObjs_ToModify[0], rd.ObjRef):
        rdBs = [o.Object() for o in rhObjs_ToModify]
        gBs = [rdB.Id for rdB in rdBs]
        edge = rhObjs_ToModify[0].Edge()
        if edge is None:
            idx_Edges = None
        else:
            rdBs, idx_Edges = createSortedBrepsAndEdges_FromObjRefs(rhObjs_ToModify)

    #print(rdBs_In)
    #print(idx_Edges)

    gB_Modified = []


    for iB, rdB in enumerate(rdBs):

        Rhino.RhinoApp.SetCommandPrompt(
                prompt="Searching brep {} of {} ...".format(iB+1, len(rdBs)))

        rgB_In = rdB.Geometry


        bModified = processBrep(
            rgBrep=rdB.Geometry,
            idx_Edges=idx_Edges[iB] if idx_Edges else None,
            bDecimalTol=bDecimalTol,
            fTol = fTol,
            bEcho=True if len(rdBs) == 1 else False,
            bDebug=bDebug,
            )

        if not bModified:
            continue

        if rdB.CommitChanges():
            gB_Modified.append(rdB.Id)
            if bEcho and len(rdBs) > 1:
                print("Brep {} was modified.".format(rdB.Id))

    if gB_Modified:
        print("{} breps were modified.".format(len(gB_Modified)))
    else:
        print("No breps were modified.")

    return gB_Modified


def main():

    rhBrepObjs_In = getInput()
    if rhBrepObjs_In is None: return

    bDecimalTol = Opts.values['bDecimalTol']
    if bDecimalTol:
        fTol = Opts.values['fTol_Decimal']
    else:
        fTol = Opts.values['fTol_Absolute']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    processBrepObjects(
        rhBrepObjs_In,
        bDecimalTol=bDecimalTol,
        fTol=fTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()