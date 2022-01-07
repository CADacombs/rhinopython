"""
Complement to _FilletEdge.  

This script uses headless FreeCAD ( https://wiki.freecadweb.org/Headless_FreeCAD )

The most significant differences from Rhino's _FilletEdge are at the
intersections of non-tangent edges.

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
Allow multiple fillet radius input similar to _FilletEdge.
Add preview of edges selected.  Just display text dots with assigned radius?
Allow edges to be unselected.
Allow variable radius fillet input.
    FreeCAD only provides variable fillets whose surfaces are G1 blended
    between the edges so that the target fillet radii can G1 match to
    an adjacent fillet at a tangent edge.  They are not linearly blended.


Send any questions, comments, or script development service needs to @spb on
the McNeel Forums ( https://discourse.mcneel.com/ )
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
211222-220103: Created.
220106: subprocess.Popen now runs in a threading.Timer.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

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

    key = 'bDirectlyFilletShell'; keys.append(key)
    values[key] = True
    names[key] = 'ObjToFillet'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Feature', 'Shell')
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

        if key == 'fRadius':
            if cls.riOpts[key].CurrentValue < 2.0*sc.doc.ModelAbsoluteTolerance:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get BrepEdges with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edges to fillet")

    go.GeometryFilter = rd.ObjectType.Curve # Curve is also used for brep edges.
    go.GeometryAttributeFilter = (
            ri.Custom.GeometryAttributeFilter.MatedEdge |
            ri.Custom.GeometryAttributeFilter.EdgeCurve)

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go_Main on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fRadius')
        addOption('bDirectlyFilletShell')
        addOption('fFreeCADTimeoutSecs')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()

            return (
                objrefs,
                Opts.values['fRadius'],
                Opts.values['bDirectlyFilletShell'],
                Opts.values['fFreeCADTimeoutSecs'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def sortInputPerBrep(rhObjs_In):

    gBs = [] # This is used as the index for the brep in the following routine.
    idx_Es_PerBrep = []


    def getDataFromInput(rhObj):
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


    for rhObj in rhObjs_In:

        rc = getDataFromInput(rhObj)
        if not rc: return

        gB, idxE = rc

        if gB in gBs:
            if idxE in idx_Es_PerBrep[gBs.index(gB)]: continue
            idx_Es_PerBrep[gBs.index(gB)].append(idxE)
            continue

        gBs.append(gB)
        idx_Es_PerBrep.append([idxE])

    return gBs, idx_Es_PerBrep


def fc(s=''): s_FC.append(s)


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
        for iF_A in rgB.Edges[iE].AdjacentFaces():
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


def create_FC_code_Determine_edges_to_fillet(bDirectlyFilletShell):
    fc('def isEdgeMid_in_input_Mids(edge, midVs):')
    fc('    try:')
    fc('        midE=edge.valueAt(edge.Curve.parameterAtDistance(edge.Length/2.0, edge.FirstParameter))')
    fc('    except:')
    fc('        print(edge.Curve)')
    fc('        Part.show(edge)')
    fc('        Part.Vertex(edge.Curve.parameterAtDistance(0.0))')
    fc('        1/0')
    fc('    for v in midVs:')
    fc('        if midE.distanceToPoint(v) < 0.01:')
    fc('            return True')
    fc('    return False')
    fc()
    fc('edges_to_fillet = []')
    fc('idx_edges_to_fillet = []')


    # Edge indices are base 1 in FreeCAD.

    fc('for iE, edge in enumerate({}, start=1):'.format(
        'shell.Edges' if bDirectlyFilletShell else 'doc.Fillet.Base.Shape.Edges'))
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
    fc('    if isEdgeMid_in_input_Mids(edge, midVs):')
    fc('        edges_to_fillet.append(edge)')
    fc('        idx_edges_to_fillet.append(iE)')
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


def processBrepShells_FcShell(rgBs_ToFillet, fRadius, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        fRadius: float
        bEcho: bool
        bDebug: bool
    """

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt


    fc('shapes = []')


    for rgB_ToFillet in rgBs_ToFillet:

        Rhino.RhinoApp.SetCommandPrompt("Exporting STEP ...")

        if not export_to_STEP(rgB_ToFillet, bDebug): return

        Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt_In)

        fc('shell = Part.Shell()')

        fc('shell.read(r"{}")'.format(sPath_STEP_to_FC))

        fc('shell.isValid()')
        fc('if not shell.isValid():')
        fc('    shell.fix(1e-6, 1e-6, 1e-6)')
        fc()
        fc('if not shell.isValid():')
        fc('    raise Exception("Shell was not fixed.")')
        fc()

        #fc('Part.show(shell)')
        #fc('Part.show(shell)')


        # For debugging.
        #fc()
        #script_FC = "\n".join(s_FC)
        #with open(sPath_Script_for_FC, "w") as f: 
        #    f.write(script_FC) 
        #return True, None
        #


        create_FC_code_Determine_edges_to_fillet(bDirectlyFilletShell=True)



        # Fillet.
        fc('shape = shell.makeFillet({}, edges_to_fillet)'.format(fRadius))
        fc()
        fc('try:')
        fc('    shape = shell.makeFillet({}, edges_to_fillet)'.format(fRadius))
        fc('except Exception as e:')
        fc('    raise Exception(traceback.format_exc())')
        fc()
        fc('shapes.append(shape)')
        fc()
        fc()


    # Create STEP file for Rhino.
    fc('compound = Part.Compound(shapes)')
    fc('compound.exportStep(r"{}")'.format(sPath_STEP_from_FC))


def processBrepShells_FcFeature(rgBs_ToFillet, fRadius, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        fRadius: float
        bEcho: bool
        bDebug: bool
    """

    sCmdPrompt_In = Rhino.RhinoApp.CommandPrompt


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


        fc('doc.addObject("Part::Fillet","spbPartFeature_Fillet")')
        fc('doc.Fillet.Base = doc.Shape')
        fc('doc.recompute()')


        create_FC_code_Determine_edges_to_fillet(bDirectlyFilletShell=False)


        # Fillet.
        fc('data_for_featFillet = []')
        fc('for iE in idx_edges_to_fillet:')
        fc('    data_for_featFillet.append((iE,{},{}))'.format(fRadius, fRadius))
        fc()
        fc('doc.Fillet.Edges = data_for_featFillet')
        fc('doc.recompute()')
        fc('doc.Shape.Visibility = False')
        fc()
        fc()


    # Create STEP file for Rhino.
    fc('import Import')
    fc('Import.export([FreeCAD.ActiveDocument.Fillet], r"{}")'.format(sPath_STEP_from_FC))


def processBrep(rgBrep_In, idxs_Edges, fRadius, bDirectlyFilletShell=True, fFreeCADTimeoutSecs=10.0, bEcho=True, bDebug=False):
    """
    Parameters:
        rgBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        fRadius: float
        bDirectlyFilletShell: bool  False means to fillet feature.
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


        fRadius *= fScaleUnit


    fc('import FreeCAD as App')
    fc('import Part')
    fc()
    fc()


    sVar_Vs = create_FC_code_Point3Ds_to_Vectors(pts_Mids)


    if bDirectlyFilletShell:
        rc = processBrepShells_FcShell(
            rgBs_ToFillet,
            fRadius=fRadius,
            bEcho=bEcho,
            bDebug=bDebug)
    else:
        rc = processBrepShells_FcFeature(
            rgBs_ToFillet,
            fRadius=fRadius,
            bEcho=bEcho,
            bDebug=bDebug)


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


def processBrepObject(rhBrep_In, idxs_Es, fRadius, bDirectlyFilletShell=True, fFreeCADTimeoutSecs=10.0, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBrep_In: ObjRef or tuple(rd.BrepObject or GUID, float(edge index))
        idxs_Es: list(int(EdgeIndex))
        fRadius: float
            True to subprocess
            False to send code to Clipboard for FreeCAD GUI.
                Code to send result to Rhino is not prepared.
        bDirectlyFilletShell: bool
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
        fRadius=fRadius,
        bDirectlyFilletShell=bDirectlyFilletShell,
        fFreeCADTimeoutSecs=fFreeCADTimeoutSecs,
        bEcho=bEcho,
        bDebug=bDebug)

    if not bSuccess: return

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

    rdBs_from_FC = sc.doc.Objects.AllObjectsSince(uInt32_MostRecent)

    Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt_In)


    if not rdBs_from_FC:
        if bEcho:
            print("Filleted breps were not imported.")
        return

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

    rc = getInput()
    if rc is None: return

    (
        objrefs,
        fRadius,
        bDirectlyFilletShell,
        fFreeCADTimeoutSecs,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gBreps_Res_All = []

    rc = sortInputPerBrep(objrefs)
    if not rc: return

    gBs_In, idx_rgEdges_PerBrep = rc

    zipped = zip(gBs_In, idx_rgEdges_PerBrep)

    for iB, (gB_In, idxs_Es) in enumerate(zipped):

        sCmdPrompt = "Processing brep {} of {}".format(iB+1, len(gBs_In))

        Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt)

        rc = processBrepObject(
            gB_In,
            idxs_Es,
            fRadius=fRadius,
            bDirectlyFilletShell=bDirectlyFilletShell,
            fFreeCADTimeoutSecs=fFreeCADTimeoutSecs,
            bEcho=bEcho,
            bDebug=bDebug)

        if rc:
            gBreps_Res_All.extend(rc)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()