"""
Complement to _FilletEdge.  

This script uses headless FreeCAD ( https://wiki.freecadweb.org/Headless_FreeCAD )

The most significant differences from Rhino's _FilletEdge are at the
intersections of non-tangent edges.

Warning: Check the resultant B-Reps for missing faces.

This script been tested on
Rhino 7, 7.13.21348.13001 on Windows
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
For previous input option, should the input be merged with any current input?
Maybe?: Allow variable radius fillet input.
    FreeCAD only creates variable fillets that have 2 radii and are G1 blended
    between the edges so that the target fillet radii can G1 match to
    an adjacent fillet at a tangent edge.  They are not linearly blended.


Send any questions, comments, or script development service needs to @spb on
the McNeel Forums ( https://discourse.mcneel.com/ )
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
211222-220103: Created.
220106: subprocess.Popen now runs in a threading.Timer.
220110: Removed code for PartDesign Fillet since it has less options than Part Fillet.
        Added multiple fillet radius input capability similar to _FilletEdge.
220111: Added preview of text dots on edges containing their assigned radius.
        Added option to use previous input.
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


    key = 'fRadius'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

        if key == 'fRadius':
            if cls.riOpts[key].CurrentValue < 2.0*sc.doc.ModelAbsoluteTolerance:
                cls.riOpts[key].CurrentValue = 0.0

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(bPrevBrepsArePresent):
    """
    Get BrepEdges with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edges to fillet")
    go.SetCommandPromptDefault("Confirm edges & radius")

    go.GeometryFilter = rd.ObjectType.Curve # Curve is also used for brep edges.
    go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.MatedEdge |
            ri.Custom.GeometryAttributeFilter.EdgeCurve)

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go_Main on repeats of While loop.
    go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(enable=True, acceptZero=True)

    go.AcceptNothing(True)

    idxs_Opt = {}

    sc.doc.Objects.UnselectAll()

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fRadius')
        if bPrevBrepsArePresent:
            idxs_Opt['UsePrevInput'] = go.AddOption('UsePrevInput')
        addOption('fFreeCADTimeoutSecs')
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
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        if go.Option().Index == idxs_Opt['UsePrevInput']:
            go.Dispose()
            return 'UsePrevInput'

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


def create_subBs_per_Es_to_fillet(rgB, idxs_Es):
    """
    """


    # Don't bother creating subsets when the face count is low.
    if rgB.Faces.Count < 13:
        return [rgB.DuplicateBrep()], range(rgB.Faces.Count)


    iFs = []

    for iE in idxs_Es:
        for iF in rgB.Edges[iE].AdjacentFaces():
            if iF in iFs: continue
            iFs.append(iF)


    # Add faces adjacent to all faces in iFs so that the faces at the ends
    # of the fillets are processed correctly.
    iEs_NoFillet = []
    for iF in iFs:
        for iE in rgB.Faces[iF].AdjacentEdges():
            if iE in idxs_Es: continue
            if iE in iEs_NoFillet: continue
            iEs_NoFillet.append(iE)

    for iE in iEs_NoFillet:
        for iF_A in rgB.Edges[iE].AdjacentFaces():
            if iF_A not in iFs:
                iFs.append(iF_A)


    # In case there are any faces missing from the sub-Brep.
    iEs_SharedVs = []
    for iE_Fillet in idxs_Es:
        v = rgB.Edges[iE_Fillet].StartVertex
        for iE in v.EdgeIndices():
            if iE in idxs_Es: continue
            if iE in iEs_NoFillet: continue
            iEs_SharedVs.append(iE)
        v = rgB.Edges[iE_Fillet].EndVertex
        for iE in v.EdgeIndices():
            if iE in idxs_Es: continue
            if iE in iEs_NoFillet: continue
            iEs_SharedVs.append(iE)

    for iE in iEs_SharedVs:
        for iF_V in rgB.Edges[iE].AdjacentFaces():
            if iF_V not in iFs:
                print("Face[{}] added via Vertex routine".format(iF_V))
                iFs.append(iF_V)


    print("Face count of Brep: {}".format(rgB.Faces.Count))
    print("Face count of sub-Brep: {}".format(len(iFs)))


    # Don't bother creating subsets in some cases.
    if len(iFs) == rgB.Faces.Count:
        return [rgB.DuplicateBrep()], range(rgB.Faces.Count)

    if rgB.Faces.Count - len(iFs) < 13:
        return [rgB.DuplicateBrep()], range(rgB.Faces.Count)


    print("Creating sub-Brep for export.")

    rgBs_ToJoin = [rgB.Faces[iF].DuplicateFace(False) for iF in iFs]

    join_tol = max(
        sc.doc.ModelAbsoluteTolerance,
        max([rgB.Edges[iE].Tolerance for iE in idxs_Es]))

    return rg.Brep.JoinBreps(rgBs_ToJoin, tolerance=join_tol), iFs


def create_FC_code_Point3Ds_to_Vectors(pts_Mids):
    fc('from FreeCAD import Base')
    fc('V=Base.Vector')
    fc('midVs = []')
    for pt in pts_Mids:
        fc('midVs.append(V({}, {}, {}))'.format(pt.X, pt.Y, pt.Z))
    fc('')
    return True


def export_to_STEP(rgBrep, bDebug=False):
    rdObj_MostRecent = sc.doc.Objects.MostRecentObject()
    uInt32_MostRecent = rdObj_MostRecent.RuntimeSerialNumber

    gBrep_for_STEP = sc.doc.Objects.AddBrep(rgBrep)
    if gBrep_for_STEP == gBrep_for_STEP.Empty: return

    sc.doc.Objects.UnselectAll()

    rdB_to_FC = sc.doc.Objects.AllObjectsSince(uInt32_MostRecent)[0]

    sc.doc.Objects.Select(objectId=rdB_to_FC.Id)

    scriptToExportStp = " ".join([
        "_-Export",
        '"{}"'.format(sPath_STEP_to_FC),
        "_Schema=AP214AutomotiveDesign",
        "_ExportParameterSpaceCurves=No",
        "_SplitClosedSurfaces=No",
        "_Enter"])

    if not Rhino.RhinoApp.RunScript(scriptToExportStp, echo=bDebug): return

    if not sc.doc.Objects.Delete(rdB_to_FC): return

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


def create_FC_code_Determine_edges_to_fillet():

    fc('def idx_input_mid_match(edge):')
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
    fc('edges_to_fillet = []')
    fc('idx_edges_to_fillet = []')
    fc('rads_per_idxE = []')
    fc()
    fc('for iE, edge in enumerate(doc.Fillet.Base.Shape.Edges):')
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
    fc('    idxMid = idx_input_mid_match(edge)')
    fc('    if idxMid is None:')
    fc('        continue')
    fc('    edges_to_fillet.append(edge)')
    fc('    idx_edges_to_fillet.append(iE+1)') # Part Edge indices are base 1 in FreeCAD.
    fc('    rads_per_idxE.append(rads_per_pts_Mids[idxMid])')
    fc('    if len(edges_to_fillet) == len(midVs):')
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


def processBrep(rgBrep_In, idxs_Edges, fRadii, fFreeCADTimeoutSecs=10.0, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        idxs_Edges: list(int(EdgeIndex))
        fRadii: list(floats in idxs_Edges order)
        fFreeCADTimeoutSecs: float
        bEcho: bool
        bDebug: bool
    """

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt

    if not rgBrep_In.IsValid:
        if bEcho:
            print("{}-face brep is invalid.  Fix first.  Brep is skipped.".format(
                rgBrep_In.Faces.Count))
        return False, None


    rgBs_ToFillet, idx_Fs_to_FC = create_subBs_per_Es_to_fillet(rgBrep_In, idxs_Edges)

    if not rgBs_ToFillet:
        if bEcho:
            print("No shells to fillet in a brep.")
        return False, None


    pts_Mids = []
    for iE in idxs_Edges:
        ts = list(rgBrep_In.Edges[iE].DivideByCount(2, includeEnds=False))
        if not ts:
            print("Midpoint of edge[{}] not found.".format(iE))
            return False, None
        pt = rgBrep_In.Edges[iE].PointAt(ts[0])
        pts_Mids.append(pt)


    # Internal distance unit of FreeCAD is millimeter.
    if sc.doc.ModelUnitSystem != Rhino.UnitSystem.Millimeters:
        fScaleUnit = Rhino.RhinoMath.UnitScale(
            sc.doc.ModelUnitSystem,
            Rhino.UnitSystem.Millimeters) # from, to

        # Breps sent to FC via STEP do not need to be scaled.

        for iP, pt in enumerate(pts_Mids):
            pts_Mids[iP].X = pt.X * fScaleUnit
            pts_Mids[iP].Y = pt.Y * fScaleUnit
            pts_Mids[iP].Z = pt.Z * fScaleUnit


        #fRadius *= fScaleUnit
        fRadii = [r*fScaleUnit for r in fRadii]


    fc('import FreeCAD as App')
    fc('import Part')
    fc()
    fc()


    if not create_FC_code_Point3Ds_to_Vectors(pts_Mids):
        raise Exception("Could not create_FC_code_Point3Ds_to_Vectors.") 

    fc('rads_per_pts_Mids = []')
    for r in fRadii:
        fc('rads_per_pts_Mids.append({})'.format(r))


    for rgB_ToFillet in rgBs_ToFillet:

        Rhino.RhinoApp.SetCommandPrompt("Exporting STEP ...")

        if not export_to_STEP(rgB_ToFillet, bDebug): return

        Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt_In)

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


        fc('doc.addObject("Part::Fillet","Fillet")')
        fc('doc.Fillet.Base = doc.Shape')
        fc('doc.recompute()')


        create_FC_code_Determine_edges_to_fillet()


        # Fillet.
        fc('data_for_featFillet = []')
        fc('for i, iE in enumerate(idx_edges_to_fillet):')
        fc('    data_for_featFillet.append((iE, rads_per_idxE[i], rads_per_idxE[i]))')
        fc()
        fc('doc.Fillet.Edges = data_for_featFillet')
        fc('doc.recompute()')
        fc('doc.Shape.Visibility = False')
        fc()
        fc()


    # Create STEP file for Rhino.
    fc('import Import')
    fc('Import.export([FreeCAD.ActiveDocument.Fillet], r"{}")'.format(sPath_STEP_from_FC))


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

    if not rc: return False, None


    def create_subBs_per_Fs_to_remove(rgB, iFs_to_Remove):
        """
        """

        rgB_Out = rgB.DuplicateBrep()

        for iF in sorted(iFs_to_Remove, reverse=True):
            rgB_Out.Faces.RemoveAt(iF)

        return rgB_Out


    if len(idx_Fs_to_FC) == rgBrep_In.Faces.Count:
        rgB_Fs_not_sent_to_FC = None
    else:
        rgB_Fs_not_sent_to_FC = create_subBs_per_Fs_to_remove(rgBrep_In, idx_Fs_to_FC)


    return True, rgB_Fs_not_sent_to_FC


def processBrepObject(rhBrep_In, idxs_Es, fRadii, fFreeCADTimeoutSecs=10.0, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        idxs_Es: list(int(EdgeIndex))
        fRadii: list(floats in idxs_Es order)
        fFreeCADTimeoutSecs: float
        bEcho: bool
        bDebug: bool
    """

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt

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
        print("Brep {} is invalid.  Fix first.  Brep is skipped.".format(gB_In))
        return

    bSuccess, rgB_Fs_not_sent_to_FC = processBrep(
        rgB_In,
        idxs_Es,
        fRadii=fRadii,
        fFreeCADTimeoutSecs=fFreeCADTimeoutSecs,
        bEcho=bEcho,
        bDebug=bDebug)
    if not bSuccess: return

    rdBs_from_FC = import_STEP_into_Rhino(bDebug)
    if not rdBs_from_FC:
        if bEcho:
            print("Filleted breps were not imported.")
        return


    Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt_In)


    if rgB_Fs_not_sent_to_FC is None:
        rgBs_forOut = [rdB.BrepGeometry for rdB in rdBs_from_FC]
    else:
        join_tol = max(
            sc.doc.ModelAbsoluteTolerance,
            max([rgB_In.Edges[iE].Tolerance for iE in idxs_Es]))

        rgBs_forOut = rg.Brep.JoinBreps(
            [rgB_Fs_not_sent_to_FC] + [rdB.BrepGeometry for rdB in rdBs_from_FC],
            tolerance=join_tol)

        if len(rgBs_forOut) > 1:
            if bEcho:
                print("More than 1 brep resulted from join.")

    if len(rgBs_forOut) == 1:
        sc.doc.Objects.Replace(gB_In, rgBs_forOut[0])

    if len(rgBs_forOut) > 1:
        for rgB_Out in rgBs_forOut:
            sc.doc.Objects.AddBrep(rgB_Out, rdB_In.Attributes)
        sc.doc.Objects.Delete(rdB_In)

    for rdB in rdBs_from_FC:
        sc.doc.Objects.Delete(rdB)


def main():


    gBs_In = [] # This is used as the index for the brep in the following routine.
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
        fFreeCADTimeoutSecs = Opts.values['fFreeCADTimeoutSecs']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        sc.doc.Objects.UnselectAll()
        sc.doc.Views.Redraw()

        if not rc:
            print("Proceeding to create fillets.")
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

        s_FC[:] = [] # So that a new FreeCAD script is made for each brep.

        rc = processBrepObject(
            gB_In,
            idxs_Es,
            fRadii=rads_for_idxEs,
            fFreeCADTimeoutSecs=fFreeCADTimeoutSecs,
            bEcho=bEcho,
            bDebug=bDebug)

        if rc:
            gBreps_Res_All.extend(rc)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()