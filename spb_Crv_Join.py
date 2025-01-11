"""
This script is an alternative to _Join for curves.

Unlike _Join,
1. It can optionally join curves that share a layer and/or color.
2. It optionally does not simplify/convert curves.
3. It has a tolerance setting.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180207: Created, starting with joinBreps.
...
180316: Added support for brep edge curves.
...
190513: Fixed bug where breps of selected edges were being joined.
250109-10: Reenabled the RhinoCommon method.
        Disabled support for BrepEdges until script is working properly for just wires.
        Refactored.

TODO:
    Reenable support for BrepEdges.
    Check whether there are benefits for using RC's Curve.JoinCurves.

_Join (at least in V8.14) appears to use 1.8 * ModelAbsoluteTolerance as the
endpoint distance (gap or overlap) threshold to perform the join.

simpleJoin in Curve.JoinCurves will preserve curve types. That also means contiguous
LineCurves will connect as LineCurve segments in PolyCurves instead of converting to PolylineCurve segment.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fJoinTol'; keys.append(key)
    values[key] = max((0.1 * sc.doc.ModelAbsoluteTolerance, 1e-6)) # Default in _Join is 1.8 * sc.doc.ModelAbsoluteTolerance.
    names[key] = 'Tolerance'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bByLayer'; keys.append(key)
    values[key] = True
    names[key] = 'JoinByLayer'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bByColor'; keys.append(key)
    values[key] = True
    names[key] = 'JoinByColor'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseUiJoinCmd'; keys.append(key)
    values[key] = False
    names[key] = 'Use'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'RhinoCommon', 'JoinCommand')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSimpleJoin'; keys.append(key)
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

        if key == 'fJoinTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.values[key] = cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return
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


def _isCurveJoinable(rdCrv):
    """Joinable is valid and open."""
    rgC = rdCrv.CurveGeometry
    if not rgC.IsValid: return False
    return not rgC.IsClosed


def _getAllNormalJoinableWires():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Curve
    return [rdC for rdC in sc.doc.Objects.GetObjectList(oes) if _isCurveJoinable(rdC)]


def getInput():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = (
        ri.Custom.GeometryAttributeFilter.OpenCurve
        |
        ri.Custom.GeometryAttributeFilter.WireCurve
        )

    go.SubObjectSelect = False

    #go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.

    #go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    #bSequentialSel = False
    #print("Seq1Sel = Sequential single selection")

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()


        addOption('fJoinTol')
        idxs_Opt['Toggle'] = go.AddOption("ToggleFilters")
        addOption('bByLayer')
        addOption('bByColor')
        #addOption('bUseUiJoinCmd')
        if not Opts.values['bUseUiJoinCmd']:
            addOption('bSimpleJoin')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=2, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked:
            if (
                    go.ObjectsWerePreselected and
                    go.ObjectCount == 1
            ):
                if go.ObjectCount == 1:
                    if bDebug: print("Only 1 object was preselected.")
                    bPreselectedObjsChecked = True
                    go.EnablePreSelect(enable=False, ignoreUnacceptablePreselectedObjects=True)
                    go.AlreadySelectedObjectSelect = True
                    continue
                elif bDebug:
                    print("{} objects were preselected.".format(go.ObjectCount))
            else:
                # This is another good place to print join option settings.
                pass

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return _getAllNormalJoinableWires()

        if res == ri.GetResult.Object:
            if go.ObjectCount < 2:
                print("At least 2 breps must be selected.")
                go.Dispose()
                return

            objrefs = go.Objects()
            go.Dispose()

            return objrefs

        #if res == ri.GetResult.Option and go.OptionIndex() == 1:
        #    bSequentialSel = True
        #    break
        
        if res == ri.GetResult.Number:
            key = 'fJoinTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        if 'Toggle' in idxs_Opt and go.Option().Index == idxs_Opt['Toggle']:
            for key in 'bByLayer', 'bByColor':
                Opts.riOpts[key].CurrentValue = not Opts.riOpts[key].CurrentValue
                Opts.setValue(key)
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break

    1/0
    # Unnecessary due to go.GetMultiple' minimumNumber=2 ?
    #    # At least 2 curves must be selected.
    #    if go.ObjectCount < 2: return
    
    gCrvs_FromWire = []; gCrvs_FromEdge = []
    
    for objref in go.Objects():
        if objref.GeometryComponentIndex.Index == -1:
            gCrvs_FromWire.append(objref.ObjectId)
        else:
            rdObj = objref.Object()
            if rdObj.ObjectType == rd.ObjectType.InstanceReference:
                print("Objects in block instances are not supported.")
                continue
            
            # Add curves of edges.
            gCrv = sc.doc.Objects.AddCurve(objref.Geometry())
            if gCrv != Guid.Empty:
                gCrvs_FromEdge.append(gCrv)
    
    go.Dispose()


def join_Curve_RhinoCommon(rgCrvs_In, fJoinTol=None, bSimpleJoin=False, bDebug=False):
    """
    Returns:
        list(rg.Curves from Join)
        list(int of indices of rgCrvs_In)
        str: Error
    """

    if len(rgCrvs_In) == 0:
        return None, None, "No curves in set to join."

    if len(rgCrvs_In) < 2:
        return None, None, "Only 1 curves in set to join."

    if fJoinTol is None:
        fJoinTol = 1.8 * sc.doc.ModelAbsoluteTolerance # The same as _Join for curves (per Rhino 8.14).

    rgCs_Joined, idxs_Out_per_In = rg.Curve.JoinCurves(
        rgCrvs_In,
        joinTolerance=fJoinTol,
        preserveDirection=False,
        simpleJoin=bSimpleJoin)

    #sEval = "len(rgCs_Joined)"; print(sEval, '=', eval(sEval))
    #sEval = "max(idxs_Out_per_In) + 1"; print(sEval, '=', eval(sEval))
    #sEval = "idxs_Out_per_In"; print(sEval, '=', eval(sEval))
    #print(idxs_Out_per_In)

    #for rgC_Joined in rgCs_Joined:
    #    sc.doc.Objects.AddCurve(rgC_Joined)
    #sc.doc.Views.Redraw()
    #return

    idxs_Out_per_In = list(idxs_Out_per_In)

    rgCs_Out = [] # Actual joins, not all output from JoinCurves.
    idxs_In_perOut = []
    for idx_Joined in range(len(rgCs_Joined)):
        if idxs_Out_per_In.count(idx_Joined) > 1:
            rgCs_Out.append(rgCs_Joined[idx_Joined])
            idxs_In_perOut.append([])

    #sEval = "len(rgCs_Out)"; print(sEval, '=', eval(sEval))


    for idx_In, idx_Out_per_In in enumerate(idxs_Out_per_In):
        iCt = idxs_Out_per_In.count(idx_Out_per_In)
        if iCt > 1:
            idxs_In_perOut[idx_Out_per_In].append(idx_In)

    return rgCs_Out, idxs_In_perOut, None

    """
    For simpleJoin=True,
    Set True to use the simple joining method. In general, set this parameter to false.
    (https://developer.rhino3d.com/api/rhinocommon/rhino.geometry.curve/joincurves)
    
    RhinoMergeCurves()
    Description: Join a bunch of ON_Curves into one or more ON_Curves Parameters:
    input_curves [in] Array of pointers to ON_Curves to be joined output [out] Array
    of pointers to joined results join_tol [in] max distance between endpoints to be joined.
    If join_tol < ON_EPSILON, use CRhinoDoc::AbsoluteTolerance() bPreserveDir [in] if TRUE,
    don't reverse input curves to get them to join key [out] if non-null, curves[i] is part of output[key[i]] WARNING - key[i] may be -1 for some i, in particular if curves[i] is extremely short, key[i] will be -1 and curves[i] will not contribute to the joined results. Returns: @untitled table TRUE Success False Failure Remarks: Join as many of the input curves as have matching endpoints. If the input curve is a NURBS curve or a line, the endpoints within the specified tolerance are trued up to meet exactly. All of the input curves are copied and the caller must free the results. When curves are joined they are made into polycurves. Memory for the curves is allocated and becomes the responsibility of the caller.
    (https://developer.rhino3d.com/api/cpp/group___rhino.html#ga422fa1a9386624b9cbdb77faaeb1bc5b)
    """


def join_Cmd(gCrvs0, fJoinTol=None, bEcho=False, bDebug=False):
    if rs.SelectObjects(gCrvs0) == 0:
        print(s + "Error!  Crvs should be selected.  Command stopped.")
        return
    
    fModelTol_Start = sc.doc.ModelAbsoluteTolerance
    fModelTol_ForJoin = sc.doc.ModelAbsoluteTolerance = fJoinTol / 1.8
    print("ModelAbsoluteTolerance changed from {:.15f} to {:.15f}" \
          " to use {:.15f} ({:.15f} / 1.8) join tolerance.".format(
            fModelTol_Start,
            sc.doc.ModelAbsoluteTolerance,
            fJoinTol,
            fJoinTol))
    
    rs.Command("_Join", bEcho)
    
    sc.doc.ModelAbsoluteTolerance = fModelTol_Start
    print("ModelAbsoluteTolerance changed from {:.15f} to {:.15f}.".format(
            fModelTol_ForJoin,
            sc.doc.ModelAbsoluteTolerance))
    
    sc.doc.Objects.UnselectAll()


def _getCurveObjects(rhCrvs):
    if all((isinstance(rhC, rd.CurveObject) for rhC in rhCrvs)): return rhCrvs

    rdCs_In = []
    for rhC in rhCrvs:
        if isinstance(rhC, rd.CurveObject):
            rdCs_In.append(rhC)
        else:
            rdC = rs.coercerhinoobject(rhC)
            if rdC.ObjectType != rd.ObjectType.Curve:
                raise Exception("{} passed to separateInputIntoJoiningSets.".format(
                    rdC.GetType().Name))
            rdCs_In.append(rdC)
    return rdCs_In


def _separateInputIntoJoiningSets(rhCrvs, bByLayer=True, bByColor=True, bDebug=False):

    rdCs_In = _getCurveObjects(rhCrvs)

    if not (bByLayer or bByColor):
        return [rdCs_In] # Nested 1-level for actual separation (see below).


    if bByLayer and bByColor:
        rdCs_Per_to_join = [] # Nested 1-level.
        sLayers = [] # Full layer paths (with '::').
        colors = [] # May be System.Drawing.Color or int for rs.ObjectColorSource.
        for rdC in rdCs_In:
            sLayer = rs.ObjectLayer(rdC) # Full layer path (with '::').
            ocs = rs.ObjectColorSource(rdC)
            color = rs.ObjectColor(rdC) if ocs == 1 else rs.ObjectColorSource(rdC)

            if (sLayer not in sLayers) or (color not in colors):
                rdCs_Per_to_join.append([rdC])
                sLayers.append(sLayer)
                colors.append(color)
                continue

            idxs_Layer_match = [idx for idx, _ in enumerate(sLayers) if _ == sLayer]
            idxs_color_match = [idx for idx, _ in enumerate(colors) if _ == color]
            matches = list(set(idxs_Layer_match) & set(idxs_color_match))
            if len(matches) == 0:
                rdCs_Per_to_join.append([rdC])
                sLayers.append(sLayer)
                colors.append(color)
            elif len(matches) == 1:
                idx_Match = matches[0]
                rdCs_Per_to_join[idx_Match].append(rdC)
            else:
                raise Exception("Multiple layer and color matches.")

        if bDebug:
            sEval = "sLayers"; print(sEval, '=', eval(sEval))
            sEval = "colors"; print(sEval, '=', eval(sEval))

        return rdCs_Per_to_join


    if bByLayer:
        rdCs_Per_to_join = [] # Nested 1-level.
        sLayers = [] # Full layer paths (with '::').
        for rdC in rdCs_In:
            sLayer = rs.ObjectLayer(rdC) # Full layer path (with '::').
            if sLayer in sLayers:
                idx = sLayers.index(sLayer)
                rdCs_Per_to_join[idx].append(rdC)
            else:
                rdCs_Per_to_join.append([rdC])
                sLayers.append(sLayer)

        if bDebug: sEval = "sLayers"; print(sEval, '=', eval(sEval))

        return rdCs_Per_to_join


    if bByColor:
        rdCs_Per_to_join = [] # Nested 1-level.
        colors = [] # May be System.Drawing.Color or int for rs.ObjectColorSource.
        for rdC in rdCs_In:
            ocs = rs.ObjectColorSource(rdC)
            color = rs.ObjectColor(rdC) if ocs == 1 else rs.ObjectColorSource(rdC)

            if color in colors:
                idx = colors.index(color)
                rdCs_Per_to_join[idx].append(rdC)
            else:
                rdCs_Per_to_join.append([rdC])
                colors.append(color)

        if bDebug: sEval = "colors"; print(sEval, '=', eval(sEval))

        return rdCs_Per_to_join

    raise Exception("What happened?")


def join_CurveObjects(rhCrvs_In, bByLayer=True, bByColor=True, bUseUiJoinCmd=False, bSimpleJoin=False, fJoinTol=None, bEcho=False, bDebug=False):

    rdCrvs_Separated = _separateInputIntoJoiningSets(rhCrvs_In, bByLayer, bByColor)

    gCs_Out = []
    iCt_gCs_Closed = 0
    gCs_DeleteFails = []
    idxs_rgCrvs_perSet_Joined = []

    for rdCrvs_PerSet in rdCrvs_Separated:
        rgCrvs_perSet = [rdC.CurveGeometry for rdC in rdCrvs_PerSet]
        if bUseUiJoinCmd:
            return # TODO: Reenable.
            # Join using _Join command not considering object layer or color.
            return join_Cmd(rgCrvs_perSet, fJoinTol, bEcho, bDebug)
        else:
            # Join using RC's JoinCurves not considering object layer or color.
            rvs = join_Curve_RhinoCommon(
                rgCrvs_perSet,
                fJoinTol=fJoinTol,
                bSimpleJoin=bSimpleJoin,
                bDebug=bDebug)
            if rvs is None:
                continue

            rgCs_Joined, idxs_rgCrvs_perSet, sLog = rvs

            if not rgCs_Joined:
                if bDebug and sLog:
                    print(sLog)
                continue


            for idx_Joined, rgC_Joined in enumerate(rgCs_Joined):

                # Use the first curve that was joined regardless of layer and color filter settings.
                attr = rdCrvs_PerSet[idxs_rgCrvs_perSet[idx_Joined][0]].Attributes

                gC_Out = sc.doc.Objects.AddCurve(rgC_Joined, attr)

                if gC_Out != Guid.Empty:
                    gCs_Out.append(gC_Out)
                    for idx_rdC_In in idxs_rgCrvs_perSet[idx_Joined]:
                        bDeleted = sc.doc.Objects.Delete(rdCrvs_PerSet[idx_rdC_In], quiet=False)
                        if bDeleted:
                            idxs_rgCrvs_perSet_Joined.append(rdCrvs_PerSet[idx_rdC_In])
                        else:
                            gCs_DeleteFails.append(rdCrvs_PerSet[idx_rdC_In].Id)

                    if rgC_Joined.IsClosed:
                        iCt_gCs_Closed += 1


    if gCs_Out:
        if not gCs_DeleteFails:
            if iCt_gCs_Closed == len(gCs_Out):
                print("Replaced {} curves with {} joined, closed curves.".format(
                    len(idxs_rgCrvs_perSet_Joined),
                    len(gCs_Out),
                    )
                      )
            else:
                print("Replaced {} curves with {} joined curves. {} of which are closed.".format(
                    len(idxs_rgCrvs_perSet_Joined),
                    len(gCs_Out),
                    iCt_gCs_Closed),
                      )
        else:
            print("Added {} joined curves.".format(len(gCs_Out)))
            print("Failed to delete {} curves.".format(len(gCs_DeleteFails)))
    else:
        print("None of the {} curves were joined.".format(len(rhCrvs_In)))

    return gCs_Out


def main():

    #gCrvs_Preselected = []
    
    #for gObj in rs.SelectedObjects():
    #    if rs.ObjectType(gObj) == rs.filter.curve:
    #        gCrvs_Preselected.append(gObj)
    
    #ct_gCrvs_Preselected = len(gCrvs_Preselected)
    
    #if ct_gCrvs_Preselected == 0:
    #    rc = getInput()
    #    if rc is None: return

    rhCrvs = getInput() # ObjRefs, rd.CurveObjects, or None.
    if rhCrvs is None: return

    fJoinTol = Opts.values['fJoinTol']
    bByLayer = Opts.values['bByLayer']
    bByColor = Opts.values['bByColor']
    bUseUiJoinCmd = Opts.values['bUseUiJoinCmd']
    bSimpleJoin = Opts.values['bSimpleJoin']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    #if gCrvs0_FromEdge:
    #    # If any edge curves were selected, curves will be joined ignoring layer and color.
    #    if bByLayer:
    #        print("Edge curve in selection, so all curves' layers will be ignored.")
    #    bByLayer = False
        
    #    if bByColor:
    #        print("Edge curve in selection, so all curves' layers will be ignored.")
    #    bByColor = False

    if not bDebug: sc.doc.Views.RedrawEnabled = True

    sc.doc.Objects.UnselectAll()
    
    gCs_Res = join_CurveObjects(
        rhCrvs,
        bByLayer=bByLayer,
        bByColor=bByColor,
        bUseUiJoinCmd=bUseUiJoinCmd,
        bSimpleJoin=bSimpleJoin,
        fJoinTol=fJoinTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    # Note that one case when rgCrvs1 is None is when _Join is used.
    if gCs_Res:
        sc.doc.Views.Redraw()

    return



    if len(gCrvs0_FromWire + gCrvs0_FromEdge) == len(rgCrvs1):
        print("Unable to join any curves.")
        if len(gCrvs0_FromEdge): sc.doc.Views.Redraw()
        return
    
    isClosedCt = 0
    
    gCrvs1 = []
    
    for rgCrv1 in rgCrvs1:
        

        
        # Don't modify monoface curves.
        if rgCrv1.Faces.Count == 1:
            continue
        
        if not rgCrv1.IsValid:
            print("The geometry of at least one of the new curves is invalid."
                  "Crv will not be modified.")
            continue
        
        if gCrv0:
            gCrv0 = gCrvs0_FromMarker[0] # Get 1st to take its attribute instead of the last's.
            
            rdCrv0 = rs.coercerhinoobject(gCrv0) # Using last gCrv0 in gCrvs0_FromMarker.
            if rdCrv0 is None and Rhino.RhinoApp.ExeVersion >= 6:
                # If using Rhino V6 or newer, try using sc.doc.Objects.Find.
                print("sc.doc.Objects.FindId failed! Using sc.doc.Objects.Find instead ...")
                rdCrv0 = sc.doc.Objects.Find(gCrv0)
            if rdCrv0 is None:
                print("DocObject for {} not found!".format(gCrv0))
                continue
            
            attr = rdCrv0.Attributes.Duplicate()
        else:
            attr = None
        
        gCrv1 = sc.doc.Objects.AddCrv(rgCrv1, attr)
        if gCrv1 == Guid.Empty:
            print("GUID is empty!  Check results.")
            continue
        
        if rgCrv1.IsClosed:
            isClosedCt += 1
        
        # Remove user string marker.
        removeMarkers(gCrv1)
        
        # Delete old curves.
        for gCrv0 in gCrvs0_FromMarker:
            sc.doc.Objects.Delete(gCrv0, quiet=False)
        
        gCrvs1.append(gCrv1)
    
    #rs.MatchObjectAttributes(gCrvs1, gCrvs0[0])
    
    s = "{} wire curves and {} curves from edges joined into".format(
            len(gCrvs0_FromWire), len(gCrvs0_FromEdge))
    sPolyInfos = []
    if isClosedCt > 0:
        sPolyInfos.append(" {} closed polycurves(s)".format(isClosedCt))
    if isClosedCt < len(gCrvs1):
        sPolyInfos.append(" {} open polycurves(s).".format(
                len(gCrvs1) - isClosedCt))
    print(s + ",".join(sPolyInfos))
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()