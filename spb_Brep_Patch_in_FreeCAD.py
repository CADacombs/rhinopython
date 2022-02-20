"""
Complement to _Patch by using Surface::Filling in FreeCAD.

Although Surface::Filling is supposed to support G2 matching, tested results
show that setting G2 produces G0 results.  Also, G1 knots exist within the new
surface, making the G2 goal moot in some ways.


This script uses headless FreeCAD ( https://wiki.freecadweb.org/Headless_FreeCAD )

This script been tested on
Rhino 7, 7.13.21348.13001 & TBD on Windows
FreeCAD 0.019.3 (0.019.24267) on Windows


Requirements to run:
1. Rhino 7.  Other versions have not been tested.  If you try doing so,
    please inform me (see below) the results.

2. FreeCAD, preferrably the same version as stated above,
    https://www.freecadweb.org/downloads.php
    installed where Rhino running the script can access it.

3. Modification of the code in the "Variables to modify per OS"... section
    below as needed.  Notice that as provided, .stp and .py communication files
    are created on the current user's Desktop.


TODO:
Do not export redundant breps, as when multiple edges are selected on a face.
For previous input option, should the input be merged with any current input?


Send any questions, comments, or script development service needs to @spb on
the McNeel Forums ( https://discourse.mcneel.com/ )
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220112-17: Created, starting with "FilletEdge"-in-FreeCAD script.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

import os
from subprocess import Popen, PIPE
from threading import Timer



# Variables to modify per OS, installation paths, and desired placement of
# .stp and .py files for communication between Rhino and FreeCAD.

sFCcmd_Path = '"{}"'.format(r"C:\Program Files\FreeCAD 0.19\bin\FreeCADcmd.exe")

desktop = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')

sPath_Script_for_FC = desktop + r"\to_FC.py"

sPath_STEP_to_FC = desktop + r"\to_FC.stp"

sPath_STEP_from_FC = desktop + r"\from_FC.stp"

#



s_FC = []

def fc(s=''): s_FC.append(s)


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAutoCompleteLoop'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iFilter'; keys.append(key)
    values[key] = 0
    names[key] = 'Filter'
    listValues[key] = 'EdgesAndWires', 'EdgesOnly', 'WiresOnly'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bG1_NotG0'; keys.append(key)
    values[key] = True
    names[key] = 'Continuity'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'G0', 'G1')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fFreeCADTimeoutSecs'; keys.append(key)
    values[key] = 10.0
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True, limit=0.0)
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

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


def ptsAreClose(ptA, ptB):
    return ptA.DistanceTo(ptB) <= sc.doc.ModelAbsoluteTolerance


def sort_crvs_to_closed_loop(crvs_In):
    """
    First curve in crvs_In will be the starting curve.
    """

    def matches(pt, idxC_Skip):
        idxs_Matches = []
        bT1s = [] # False for PointAtStart, True for PointAtEnd
        for i, crv in enumerate(crvs_In):

            if i in idxC_Skip: continue

            if ptsAreClose(pt, crv.PointAtStart):
                idxs_Matches.append(i)
                bT1s.append(False)
                continue

            if ptsAreClose(pt, crv.PointAtEnd):
                idxs_Matches.append(i)
                bT1s.append(True)

        return idxs_Matches, bT1s


    idxCs_Ordered = [0] # Indices of crvs_In.
    pt = crvs_In[0].PointAtEnd

    while True:
        sc.escape_test()

        for idxC, crv in enumerate(crvs_In):
            if (idxC != 0) and (idxC in idxCs_Ordered): continue
        
            idx_Matches, bT1s = matches(pt, [idxC]+idxCs_Ordered)
        
            if len(idx_Matches) == 0:
                if not ptsAreClose(pt, crvs_In[0].PointAtStart):
                    print(
                        "Loop doesn't close with start curve.")
                    return

                # Success.
                return idxCs_Ordered


            if len(idx_Matches) == 1:
                idx_Match = idx_Matches[0]
                idxCs_Ordered.extend(idx_Matches)
                if bT1s[0]:
                    pt = crvs_In[idx_Match].PointAtStart
                else:
                    pt = crvs_In[idx_Match].PointAtEnd
                #sc.doc.Objects.AddCurve(crvs_In[idx_Match]); sc.doc.Views.Redraw()
                break # to while loop.
        
            if len(idx_Matches) > 1:
                print(
                    "More than one branch detected.",
                    "Giving up finding single closed loop.")
                return


def midPt(crv):
    t = list(crv.DivideByCount(2, includeEnds=False))[0]
    return crv.PointAt(t)


def autoSearchForOnlyClosedBoundaryLoop(objref):
    """
    """

    crv_Start = objref.Edge()
    if crv_Start is None:
        crv_Start = objref.Curve()
        bWire = True
    else:
        bWire = False


    def closestPt(crv, pt):
        return crv.PointAt(crv.ClosestPoint(pt)[1])




    oes = rd.ObjectEnumeratorSettings()
    oes.NormalObjects = True
    oes.LockedObjects = False
    oes.IncludeLights = False
    oes.IncludeGrips = False

    if Opts.values['iFilter'] == 0:
        oes.ObjectTypeFilter = rd.ObjectType.Brep | rd.ObjectType.Curve
    elif Opts.values['iFilter'] == 1:
        oes.ObjectTypeFilter = rd.ObjectType.Brep
    elif Opts.values['iFilter'] == 2:
        oes.ObjectTypeFilter = rd.ObjectType.Curve

    crvs_Filtered = []
    gCrvs_Filtered = []

    for rdO in sc.doc.Objects.GetObjectList(oes):
        if rdO.ObjectType == rd.ObjectType.Curve:
            if bWire and rdO.Id == objref.ObjectId: continue
            wire = rdO.Geometry
            if wire.IsClosed: continue
            crvs_Filtered.append(wire)
            gCrvs_Filtered.append(rdO.Id)
            continue

        brep = rdO.Geometry

        for iE, edge in enumerate(brep.Edges):
            if edge.Valence != rg.EdgeAdjacency.Naked: continue
            if not bWire and (rdO.Id == objref.ObjectId) and (iE == crv_Start.EdgeIndex):
                continue
            crvs_Filtered.append(edge)
            gCrvs_Filtered.append(rdO.Id)


    #print(len(crvs_Filtered))


    def crvsOverlap(crvA, crvB):
        if ptsAreClose(crvA.PointAtStart, crvB.PointAtStart):
            if ptsAreClose(crvA.PointAtEnd, crvB.PointAtEnd):
                if ptsAreClose(midPt(crvA), closestPt(crvB, midPt(crvA))):
                    return True
        elif ptsAreClose(crvA.PointAtStart, crvB.PointAtEnd):
            if ptsAreClose(crvA.PointAtEnd, crvB.PointAtStart):
                if ptsAreClose(midPt(crvA), closestPt(crvB, midPt(crvA))):
                    return True


    def overlaps(crvs):
        idxOverlaps = []
        for i in range(len(crvs) - 1):
            if i in idxOverlaps: continue
            for j in range(i+1, len(crvs)):
                if j in idxOverlaps: continue
                if crvsOverlap(crvs[i], crvs[j]):
                    idxOverlaps.extend([i, j])
        return idxOverlaps

    idxs_Overlaps = overlaps(crvs_Filtered)
    #print(idxs_Overlaps)

    for i in sorted(set(idxs_Overlaps), reverse=True):
        del crvs_Filtered[i]
        del gCrvs_Filtered[i]
    #print(len(crvs_Filtered))


    #print(*crvs_Filtered, sep='\n')

    crvs_to_Sort = [crv_Start] + crvs_Filtered
    gCrvs_to_Sort = [objref.ObjectId] + gCrvs_Filtered


    idxCs_Ordered = sort_crvs_to_closed_loop(crvs_to_Sort)
    if idxCs_Ordered is None: return

    #print(idxCs_Ordered)

    return (
        [crvs_Filtered[i] for i in idxCs_Ordered],
        [gCrvs_Filtered[i] for i in idxCs_Ordered])


def getInput():
    """
    Get Curves, including BrepEdges, with optional input.
    """

    sk_PrevSel = 'UsePrevSelection({})({})'.format(__file__, sc.doc.Name)


    def getClosedLoopFromObjRefs(objrefs):

        rgCs = []
        gCs = []

        for objref in objrefs:
            crv = objref.Edge()
            if crv is None:
                crv = objref.Curve()
            rgCs.append(crv)
            gCs.append(objref.ObjectId)
        idxs_Sorted = sort_crvs_to_closed_loop(rgCs)
        if idxs_Sorted is None:
            print("AInput curves do not form a single closed loop.")
            return

        if len(idxs_Sorted) < objrefs.Count:
            print("BInput curves do not form a single closed loop.")
            return

        return (
            [gCs[idx] for idx in idxs_Sorted],
            [rgCs[idx] for idx in idxs_Sorted])


    go = ri.Custom.GetObject()

    go.GeometryFilter = rd.ObjectType.Curve # Curve is also used for brep edges.
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go_Main on repeats of While loop.
    go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opt = {}

    sc.doc.Objects.UnselectAll()

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iFilter')
        addOption('bG1_NotG0')
        addOption('bAutoCompleteLoop')
        if not Opts.values['bAutoCompleteLoop']:
            if (sk_PrevSel in sc.sticky) and sc.sticky[sk_PrevSel]:
                key = 'UsePrevSelection'; idxs_Opt[key] = go.AddOption(key)
        addOption('fFreeCADTimeoutSecs')
        addOption('bEcho')
        addOption('bDebug')

        if Opts.values['bAutoCompleteLoop']:
            go.SetCommandPrompt("Select curve/edge of boundary")
            res = go.Get()
        else:
            go.SetCommandPrompt("Select curves/edges for boundary")
            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)


        if Opts.values['iFilter'] == 0:
            go.GeometryAttributeFilter = (
                    ri.Custom.GeometryAttributeFilter.BoundaryEdge |
                    ri.Custom.GeometryAttributeFilter.EdgeCurve |
                    ri.Custom.GeometryAttributeFilter.OpenCurve |
                    ri.Custom.GeometryAttributeFilter.WireCurve)
        if Opts.values['iFilter'] == 1:
            go.GeometryAttributeFilter = (
                    ri.Custom.GeometryAttributeFilter.BoundaryEdge |
                    ri.Custom.GeometryAttributeFilter.EdgeCurve |
                    ri.Custom.GeometryAttributeFilter.OpenCurve)

        if Opts.values['iFilter'] == 2:
            go.GeometryAttributeFilter = (
                    ri.Custom.GeometryAttributeFilter.OpenCurve |
                    ri.Custom.GeometryAttributeFilter.WireCurve)



        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()

            if Opts.values['bAutoCompleteLoop']:
                rc = autoSearchForOnlyClosedBoundaryLoop(objrefs[0])

                if rc is None:
                    sc.doc.Objects.UnselectAll()
                    print("AutoCompleteLoop failed.  Option disabled.")
                    key = 'bAutoCompleteLoop'
                    Opts.riOpts[key].CurrentValue = False
                    Opts.setValue(key)
                    continue

                go.Dispose()

                sc.sticky[sk_PrevSel] = None

                return rc

            rc = getClosedLoopFromObjRefs(objrefs)
            if not rc: continue

            go.Dispose()

            sc.sticky[sk_PrevSel] = objrefs

            return rc


        if res == ri.GetResult.Number:
            num_In = int(go.Number())
            if num_In not in (0,1):
                print("Numeric input is invalid.")
                continue
            key = 'bG1_NotG0'
            Opts.riOpts[key].CurrentValue = bool(num_In)
            Opts.setValue(key)
            continue


        # An option was selected.
        if (
            'UsePrevSelection' in idxs_Opt and
            go.Option().Index == idxs_Opt['UsePrevSelection']
        ):
            objrefs = sc.sticky[sk_PrevSel]
            #print(len(objrefs))
            #for objref in objrefs:
            #    go.AppendToPickList(objref)
            #    sc.doc.Views.Redraw()
            rc = getClosedLoopFromObjRefs(objrefs)
            if not rc: continue

            go.Dispose()

            sc.sticky[sk_PrevSel] = objrefs

            return rc


        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


class DrawRadiusDotsConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        print("in __init__")
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


def create_rgs_for_FC(rgCrvs):
    """
    """

    rgObjs_Out = []

    for rgC in rgCrvs:
        if not isinstance(rgC, rg.BrepEdge):
            rgObjs_Out.append(rgC)
            continue

        # BrepEdge
        edge = rgC

        if False:
            rg.BrepEdge.Brep
            rg.BrepEdge.AdjacentFaces()

        for iF in edge.AdjacentFaces():
            face = edge.Brep.Faces[iF]
            rgB_Out = face.DuplicateFace(False)
            if not rgB_Out.IsValid:
                print("New monoface brep is not valid.  Script canceled.")
                return
            rgObjs_Out.append(rgB_Out)

    return rgObjs_Out


def create_FC_code_Point3Ds_to_Vectors(pts_Mids):
    fc('V=FreeCAD.Vector')
    fc('midVs = []')
    for pt in pts_Mids:
        fc('midVs.append(V({}, {}, {}))'.format(pt.X, pt.Y, pt.Z))
    fc('')
    return True


def export_to_STEP(rgObjs, bDebug=False):
    """
    rgObjs: list(Rhino.Geometry.Brep and/or Rhino.Geometry.Curve)
    """

    rdObj_MostRecent = sc.doc.Objects.MostRecentObject()
    uInt32_MostRecent = rdObj_MostRecent.RuntimeSerialNumber

    gObjs_for_STEP = []

    for rgObj in rgObjs:
        if isinstance(rgObj, rg.Brep):
            gObj_for_STEP = sc.doc.Objects.AddBrep(rgObj)
        else:
            gObj_for_STEP = sc.doc.Objects.AddCurve(rgObj)

        if gObj_for_STEP == gObj_for_STEP.Empty:
            print("Could not add object for a {} to document.".format(
                rgObj.GetType().Name))
            sc.doc.Objects.Delete(objectIds=gObjs_for_STEP, quiet=False)
            return

    rdObjs_toFC = sc.doc.Objects.AllObjectsSince(uInt32_MostRecent)

    sc.doc.Objects.UnselectAll()
    rc = [sc.doc.Objects.Select(objectId=rdObj.Id) for rdObj in rdObjs_toFC]

    scriptToExportStp = " ".join([
        "_-Export",
        '"{}"'.format(sPath_STEP_to_FC),
        "_Schema=AP214AutomotiveDesign",
        "_ExportParameterSpaceCurves=No",
        "_SplitClosedSurfaces=No",
        "_Enter"])

    if not Rhino.RhinoApp.RunScript(scriptToExportStp, echo=bDebug): return
    rc = [sc.doc.Objects.Delete(objectId=rdObj.Id, quiet=False) for rdObj in rdObjs_toFC]
    if None in rc: return

    return True


def import_STEP_into_Rhino(bDebug=False):
    """
    """
    rdObj_MostRecent = sc.doc.Objects.MostRecentObject()
    uInt32_MostRecent = rdObj_MostRecent.RuntimeSerialNumber

    scriptToImportStp = " ".join([
        "_-Import",
        '"{}"'.format(sPath_STEP_from_FC),
        "_JoinSurfaces=Yes",
        "_LimitFaces=No",
        "_SkipInvisibles=Yes",
        "_ShowBadObjectWarning=Yes",
        "_ShowNestedBlockWarning=No",
        "_EnterEnd"])

    Rhino.RhinoApp.SetCommandPrompt("Importing STEP ...")

    Rhino.RhinoApp.RunScript(scriptToImportStp, echo=bDebug)

    return sc.doc.Objects.AllObjectsSince(uInt32_MostRecent)


def create_FC_code_Determine_curves_for_border():


    #fc('for i, face in enumerate(shape.Faces):')


    fc('def index_input_mid_match(edge):')
    fc('    try:')
    fc('        midE=edge.valueAt(edge.Curve.parameterAtDistance(edge.Length/2.0, edge.FirstParameter))')
    fc('    except:')
    fc('        print(edge.Curve)')
    fc('        Part.show(edge)')
    fc('        Part.Vertex(edge.Curve.parameterAtDistance(0.0))')
    fc('        1/0')
    fc('    for i, v in enumerate(midVs):')
    fc('        if midE.distanceToPoint(v) < 0.01:')
    fc('            return i')
    fc()
    fc()
    fc('edges_Border = []')
    fc('idx_edges_to_fillet = []')
    fc('conts_per_edge = []')
    fc()
    fc('for iE, edge in enumerate(shape.Edges):')
    #fc('    print(edge.Length)')
    fc('    if edge.Degenerated:')
    fc('        print("Degenerated Edge skipped.".format(edge.Length))')
    fc('        continue')
    fc('    if edge.Length < 1e-6:')
    fc('        print("Edge with length {} skipped.".format(edge.Length))')
    fc('        continue')
    fc('    if edge.Length < 0.01:')
    fc('        print("Edge with length {} skipped.".format(edge.Length))')
    fc('        continue')
    fc('    idxMid = index_input_mid_match(edge)')
    fc('    if idxMid is None:')
    fc('        continue')
    fc('    edges_Border.append(edge)')
    fc('    idx_edges_to_fillet.append(iE+1)') # Part Edge indices are base 1 in FreeCAD.
    fc('    conts_per_edge.append(iConts_at_pts[idxMid])')
    fc('    if len(edges_Border) == len(midVs):')
    fc('        break')
    fc('else:')
    fc('    raise Exception("Not all the input edges were matched to the edges in FreeCAD.")')
    fc()
    fc()


def subprocess_FC(sArgs, fFreeCADTimeoutSecs=5.0, bEcho=True):
    """
    Parameters:
        sArgs: str(Although Popen also accepts a list, trying this when calling FreeCADcmd.exe has so far been unnsuccesful.)
    Returns on success: True
    Returns on fail: False

    References:
    https://web.archive.org/web/20160709205942/http://www.ironpython.info/index.php?title=Launching_Sub-Processes
    https://www.blog.pythonlibrary.org/2016/05/17/python-101-how-to-timeout-a-subprocess/
    """

    kill = lambda process: process.kill()

    p = Popen(sArgs, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True) # shell=True hides console window.

    my_timer = Timer(fFreeCADTimeoutSecs, kill, [p])


    bSuccess = False

    try:
        my_timer.start()
        stdout, stderr = p.communicate()
    except:
        pass
    else:
        if stderr or p.returncode:
            if bEcho:
                print("Failed to get FreeCAD to complete the filleting.")
                print("stderr from FreeCAD: {}".format(stderr if stderr else "None"))
                if p.returncode:
                    print("returncode: {}".format(p.returncode))
        else:
            bSuccess = True
    finally:
        my_timer.cancel()

    return bSuccess


def createBrep(rgCrvs_In, iContinuities, fFreeCADTimeoutSecs=10.0, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        idxs_Edges: list(int(EdgeIndex))
        fFreeCADTimeoutSecs: float
        bEcho: bool
        bDebug: bool
    """

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt

    #if not rgBrep_In.IsValid:
    #    if bEcho:
    #        print("{}-face brep is invalid.  Fix first.  Brep is skipped.".format(
    #            rgBrep_In.Faces.Count))
    #    return False, None


    rgObjs_ToPatch = create_rgs_for_FC(rgCrvs_In)

    if not rgObjs_ToPatch:
        return


    Rhino.RhinoApp.SetCommandPrompt("Exporting STEP ...")
    if not export_to_STEP(rgObjs_ToPatch, bDebug): return

    Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt_In)


    pts_Mids = [midPt(rgC_In) for rgC_In in rgCrvs_In]


    # Internal distance unit of FreeCAD is the millimeter.
    if sc.doc.ModelUnitSystem != Rhino.UnitSystem.Millimeters:
        fScaleUnit = Rhino.RhinoMath.UnitScale(
            sc.doc.ModelUnitSystem,
            Rhino.UnitSystem.Millimeters) # from, to

        # Breps sent to FC via STEP do not need to be scaled.

        for iP, pt in enumerate(pts_Mids):
            pts_Mids[iP].X = pt.X * fScaleUnit
            pts_Mids[iP].Y = pt.Y * fScaleUnit
            pts_Mids[iP].Z = pt.Z * fScaleUnit


    fc('import FreeCAD as App')
    fc('import Part')
    fc()
    fc()


    if not create_FC_code_Point3Ds_to_Vectors(pts_Mids):
        raise Exception("Could not create_FC_code_Point3Ds_to_Vectors.") 

    fc('iConts_at_pts = []')
    for iCont in iContinuities:
        fc('iConts_at_pts.append({})'.format(iCont))
    fc()


    fc('shape = Part.Shape()')
    fc('shape.read(r"{}")'.format(sPath_STEP_to_FC))

    fc('Part.show(shape)')

    fc('doc = FreeCAD.ActiveDocument')
    fc('doc.recompute()')



    # For debugging.
    #fc()
    #script_FC = "\n".join(s_FC)
    #with open(sPath_Script_for_FC, "w") as f: 
    #    f.write(script_FC) 
    #return True, None
    #


    fc()
    fc()


    #create_FC_code_Determine_curves_for_border()


    #fc('shapes_Es = [e.Curve.toShape() for e in edges_Border]')

    # Add Part::PartFeatures.
    fc('face_objs = []')
    fc('edge_objs = []')
    fc('fevs = [] # PartFeature Face, Edge.Name, Vector (midpoint)')
    fc('for i, face in enumerate(shape.Faces):')
    fc('    face_obj = doc.addObject("Part::Feature","Face{}".format(i))')
    fc('    face_obj.Shape = face')
    fc('    face_objs.append(face_obj)')
    fc('    for j, edge in enumerate(face_obj.Shape.Edges):')
    fc('        midE=edge.valueAt(edge.Curve.parameterAtDistance(edge.Length/2.0, edge.FirstParameter))')
    fc('        fevs.append((face_obj, j+1, midE))')
    #fc('        edge_obj = doc.addObject("Part::Feature","Edge{} {}".format(i, j))')
    #fc('        edge_obj.Shape = edge')
    #fc('        edge_objs.append(edge_obj)')
    fc()
    fc('doc.recompute()')
    fc()

    fc('fevs_ordered = []')
    fc('for i, v in enumerate(midVs):')
    fc('    for j, fev in enumerate(fevs):')
    fc('        if fev[2].distanceToPoint(v) < 0.01:')
    fc('            fevs_ordered.append(fev)')
    fc()
    fc()


    fc('# Patch')
    #fc('edge_objs = []')
    #fc('for i, edge in enumerate(edges_Border):')
    #fc('    edge_obj = doc.addObject("Part::Feature","Edge{}".format(i))')
    #fc('    edge_obj.Shape = edge')
    #fc('    edge_objs.append(edge_obj)')
    #fc()
    #fc('doc.recompute()')
    fc()
    #fc('data_for_featFilling = []')
    #fc('for i, eo in enumerate(edges_Border):')
    #fc('    data_for_featFilling.append((eo, ("Edge1".format(i))))')
    #fc()
    fc('patch = doc.addObject("Surface::Filling","Patch")')
    fc()
    #fc('doc.Patch.Base = doc.Shape')
    #fc('patch.BoundaryEdges = edge_objs')
    #fc('patch.BoundaryEdges = [(face_objs[i], "Edge1".format(i)) for i in range(len(edge_objs))]')
    fc('patch.BoundaryEdges = [(fevs_ordered[i][0], "Edge{}".format(fevs_ordered[i][1])) for i in range(len(fevs_ordered))]')
    #fc('    (edge_objs[0], "Edge1"),')
    #fc('    (edge_objs[1], "Edge1"),')
    #fc('    (edge_objs[2], "Edge1"),')
    #fc('    (edge_objs[3], "Edge1"),')
    #fc('    (edge_objs[4], "Edge1"),')
    #fc('    (edge_objs[5], "Edge1"),')
    #fc('    ]')
    fc()
    #fc('patch.BoundaryEdges = [(edges_Border[0], "Edge1"), (edges_Border[1], "Edge1")]')
    #fc('patch.BoundaryFaces = [face_objs[i].Label for i in range(len(face_objs))]')
    #fc('patch.BoundaryFaces = [')
    #fc('    face_objs[0].Label,')
    #fc('    face_objs[1].Label,')
    #fc('    face_objs[2].Label,')
    #fc('    face_objs[3].Label,')
    #fc('    face_objs[4].Label,')
    #fc('    face_objs[5].Label,')
    #fc('    ]')
    fc('patch.BoundaryFaces = ["Face1" for i in range(len(face_objs))]')
    fc()
    fc('patch.BoundaryOrder = iConts_at_pts')
    fc()
    fc('doc.recompute()')
    fc('doc.Shape.Visibility = False')
    fc()
    fc()


    # Create STEP file for Rhino.
    fc('import Import')
    fc('Import.export([doc.Patch], r"{}")'.format(sPath_STEP_from_FC))
    fc()
    fc()


    # Create script file.
    fc()
    script_FC = "\n".join(s_FC)

    with open(sPath_Script_for_FC, "w") as f: 
        f.write(script_FC) 


    # Send script to headless FreeCAD (FreeCADcmd.exe) via subprocess.
    sForConsole = sFCcmd_Path + " RunMode=Script " + '"{}"'.format(sPath_Script_for_FC)


    Rhino.RhinoApp.SetCommandPrompt("Running FreeCAD ...")

    rc = subprocess_FC(sForConsole, fFreeCADTimeoutSecs=fFreeCADTimeoutSecs, bEcho=bEcho)

    Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt_In)

    return bool(rc)


def createBrepObject(rhCrvs_In, iContinuities, fFreeCADTimeoutSecs=10.0, bEcho=True, bDebug=False):
    """
    Parameters:
        rhCrvs_In: list(ObjRef or rg.Curves)
        iContinuities: list(int(Geometric continuity order))
        fFreeCADTimeoutSecs: float
        bEcho: bool
        bDebug: bool
    """

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt


    rgObjs_For_FC = []

    #for rgC in rhCrvs_In:
    #    if isinstance(rgC, rg.BrepEdge):
            
    #    if not rgC.IsValid:

    #if not rdB_In: return


    #if not rgB_In.IsValid:
    #    print("Brep {} is invalid.  Fix first.  Brep is skipped.".format(gB_In))
    #    return

    if not createBrep(
        rhCrvs_In,
        iContinuities=iContinuities,
        fFreeCADTimeoutSecs=fFreeCADTimeoutSecs,
        bEcho=bEcho,
        bDebug=bDebug
    ):
        return

    rdB_from_FC = import_STEP_into_Rhino(bDebug)
    if not rdB_from_FC:
        if bEcho:
            print("Patch was not imported.")
        return

    return rdB_from_FC


def main():


    gBs_In = [] # This is used as the index for the brep in the following routine.
    idx_rgEdges_PerBrep = []
    rads_for_idxEs_PerBrep = []


    prevSels_Present = []

    #skey_conduit = 'conduit({})'.format(__file__)
    #if (skey_conduit in sc.sticky) and sc.sticky[skey_conduit]:
    #    conduit = sc.sticky[skey_conduit]
    #else:
    #    conduit = DrawRadiusDotsConduit()
    #    sc.sticky[skey_conduit] = conduit

    #conduit.Enabled = False # Turns off the conduit if left on from a previous execution of this script.
    #sc.doc.Views.Redraw()


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


    rc = getInput()
    if rc is None: return
    print(*rc, sep='\n')
    gBs, rgCs = rc


    #while True:

    #    if rc is None:
    #        conduit.Enabled = False
    #        del conduit
    #        del sc.sticky[skey_conduit]
    #        sc.sticky[skey_conduit] = None
    #        return

    #    fRadius = Opts.values['fRadius']
    #    fFreeCADTimeoutSecs = Opts.values['fFreeCADTimeoutSecs']
    #    bEcho = Opts.values['bEcho']
    #    bDebug = Opts.values['bDebug']

    #    sc.doc.Objects.UnselectAll()
    #    sc.doc.Views.Redraw()

    #    if not rc:
    #        print("Proceeding to create fillets.")
    #        conduit.Enabled = False
    #        del conduit
    #        break

    #    objrefs = rc

    #    rc = prepareDataForPreviewAndSticky()
    #    if not rc: continue
    #    pts_Mids_All, rads_All = rc

    #    conduit.prs = zip(pts_Mids_All, rads_All)

    #    conduit.Enabled = True
    #    sc.doc.Views.Redraw()


    #if not bDebug: sc.doc.Views.RedrawEnabled = False


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    rc = createBrepObject(
        rgCs,
        iContinuities=[int(Opts.values['bG1_NotG0']) for i in range(len(rgCs))],
        fFreeCADTimeoutSecs=Opts.values['fFreeCADTimeoutSecs'],
        bEcho=Opts.values['bEcho'],
        bDebug=Opts.values['bDebug'])


    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()