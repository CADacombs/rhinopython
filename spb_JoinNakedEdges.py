"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160303: Created.  'JoinNakedEdges' apparently is an undocumented command.
160304: Changed minimum tolerance allowed to be entered to 0.  Changed formatting of printed decimals.
160728: Modularized and added more printed output.
190625: Added Opts.  Refactored getInput.
190822: Modified an option default value.
200630: Import-related update and modified an optino default value.
230928: Bug fixes.  Added a command option.  Refactored.
231001: Refactored.  Added text (indices of bad brep components) to dots.
250215: '+' now is included in reporting of changes in edge counts.
250916: Modified an option default value.
251009-10: Added an option to join edges iteratively through tolerance values.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

#import re


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'fMaxTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bIterateToMaxTol'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fStartTol'; keys.append(key)
    values[key] = 1e-5 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMarkBadGeom'; keys.append(key)
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
        if sc.sticky.has_key(stickyKeys[key]):
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
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

        if key == 'fMaxTol':
            if cls.riOpts[key].CurrentValue < 0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            if cls.values['fStartTol'] > cls.values['fMaxTol']:
                sc.sticky[cls.stickyKeys['fStartTol']] = cls.values['fStartTol'] = cls.riOpts['fStartTol'].CurrentValue = cls.riOpts['fMaxTol'].CurrentValue
            return

        if key == 'fStartTol':
            if cls.riOpts[key].CurrentValue < 1e-6 * Rhino.RhinoMath.UnitScale(
                Rhino.UnitSystem.Millimeters, sc.doc.ModelUnitSystem):
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            if cls.values['fMaxTol'] < cls.values['fStartTol']:
                sc.sticky[cls.stickyKeys['fMaxTol']] = cls.values['fMaxTol'] = cls.riOpts['fMaxTol'].CurrentValue = cls.riOpts['fStartTol'].CurrentValue
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        if key in cls.stickyKeys:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def get_all_normal_polysrf_breps():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    rdBs_Out = []
    for rdB in sc.doc.Objects.GetObjectList(oes):
        #if not rdB.BrepGeometry.IsValid:
        #    continue
        if rdB.BrepGeometry.IsSolid:
            continue
        if rdB.BrepGeometry.Faces.Count == 1:
            continue
        rdBs_Out.append(rdB)

    return rdBs_Out


def getInput():
    """
    Get brep face.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select polysrf breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.PolysrfFilter
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.OpenPolysrf
    go.SubObjectSelect = False

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fMaxTol')
        addOption('bIterateToMaxTol')
        if Opts.values['bIterateToMaxTol']:
            addOption('fStartTol')
        addOption('bMarkBadGeom')
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
            go.Dispose()
            rdBs_Out = get_all_normal_polysrf_breps()
            if not rdBs_Out:
                print("No open polysrf breps are selectable.")
                return
            return rdBs_Out


        if res == ri.GetResult.Number:
            key = 'fMaxTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _countsOfEdgeValences(rgBrep0):
    """
    From https://developer.rhino3d.com/api/rhinocommon/rhino.geometry.edgeadjacency
    EdgeAdjacency enum
    None = 0            Edge is not used by any faces and is therefore superfluous.
    Naked = 1           Edge is used by a single face.
    Interior = 2        Edge is used by two adjacent faces.
    NonManifold = 3     Edge is used by three or more adjacent faces.
    """
    numEdgeValence = [0]*4
    for e in rgBrep0.Edges:
        numEdgeValence[e.Valence.value__] += 1
        e.Dispose()
    return numEdgeValence 


def getMaxDeviationToJoin(rgBrep0, fMaxDevToFlag):
    rgCrvs = rgBrep0.DuplicateNakedEdgeCurves(True, True)
    rgBrep0.Dispose()
    
    rs.Prompt(message="Checking deviations of naked edges ...")
    
    # Check all combinations of curves for maximum deviation.
    fMaxDev = 0.
    for i in range(len(rgCrvs)):
        for j in range(i + 1, len(rgCrvs)):
            rc = rg.Curve.GetDistancesBetweenCurves(
                    rgCrvs[i], rgCrvs[j], sc.doc.ModelAbsoluteTolerance)
            if rc[0]:
                fDev = rc[1]
                if fDev <= fMaxDevToFlag:
                    if fDev > fMaxDev: fMaxDev = fDev
        rgCrvs[i].Dispose()
    
    if fMaxDev:
        print("Largest deviation <= {:f}: {:f}".format(fMaxDevToFlag, fMaxDev))
        fMaxTol = fMaxDev
    else:
        print("No deviations <= {:f} exist.".format(fMaxDevToFlag))
        fMaxTol = 2. * fMaxDevToFlag
    
    return rs.GetReal("Enter maximum deviation to join", fMaxTol, 0.)


def _dotEdge(rgBrep, idxEdge):
    rgEdge = rgBrep.Edges[idxEdge]
    pt = rgEdge.PointAtNormalizedLength(0.5)
    if pt.X == Rhino.RhinoMath.UnsetValue:
        pt = rgEdge.PointAtStart

    rgDot = rg.TextDot(text="E{}".format(idxEdge), location=pt)
    rgDot.FontHeight = 11 if Rhino.RhinoApp.ExeVersion >= 6 else 14

    attr = rd.ObjectAttributes()
    attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex # If not done, layer of index 0 will be used.
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr.ObjectColor = Color.Red

    gDot = sc.doc.Objects.AddTextDot(rgDot, attr)
    
    if gDot != gDot.Empty: return rgDot

    print("AddTextDot failed for E[{}}.".format(idxEdge))


def _dotFace(rgBrep, idxFace):

    rgFace = rgBrep.Faces[idxFace]

    amp = rg.AreaMassProperties.Compute(rgFace)
    if amp is None:
        print("Surface not dotted because its AreaMassProperties,"
            "and thus its centroid, cannot be calculated.")
        return

    ptCentroid = amp.Centroid
    getrc, u, v = rgFace.ClosestPoint(ptCentroid)
    if rg.BrepFace.IsPointOnFace(rgFace, u, v):
        location = rgFace.PointAt(u, v)
    else:
        pts = []
        for idxE in rgFace.AdjacentEdges():
            rgE = rgBrep.Edges[idxE]
            bSuccess, t = rgE.ClosestPoint(ptCentroid)
            if bSuccess: pts.append(rgE.PointAt(t))
        pt3dlist = Rhino.Collections.Point3dList(pts)
        location = pt3dlist.ClosestPointInList(ptCentroid)

    rgDot = rg.TextDot(text="F{}".format(idxFace), location=location)
    rgDot.FontHeight = 11 if Rhino.RhinoApp.ExeVersion >= 6 else 14

    attr = rd.ObjectAttributes()
    attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex # If not done, layer of index 0 will be used.
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr.ObjectColor = Color.Red

    gDot = sc.doc.Objects.AddTextDot(rgDot, attr)

    if gDot != gDot.Empty: return rgDot

    print("AddTextDot failed for F[{}].".format(idxFace))


def _report_invalid_Brep_after_JoinNakedEdges(brep, bEcho=True, bMarkBadGeom=True):
    if bEcho:
        if bMarkBadGeom:
            print("Brep became invalid at marked area(s).")
        else:
            print("Brep became invalid at latest JoinNakedEdges.")
        sLog = brep.IsValidWithLog()[1]
        print(sLog)
    if bMarkBadGeom:
        ## Alternative
        #m = re.search(r"\[([A-Za-z0-9_]+)\]", sLog)
        #print(m.group(1))
        if 'ON_Brep.m_F[' in sLog:
            idxF = int(
                    sLog.split('ON_Brep.m_F[', 1)[1].\
                    split('] is invalid.')[0])
            _dotFace(brep, idxF)
        if 'ON_Brep.m_E[' in sLog:
            idxE = int(
                    sLog.split('ON_Brep.m_E[', 1)[1].\
                    split('] is invalid.')[0])
            _dotEdge(brep, idxE)


def _formatDistance(fDistance):
    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance == 0.0:
        return "0"

    if fDistance < 10**-(sc.doc.ModelDistanceDisplayPrecision-1):
        return "{:.4e}".format(fDistance)

    return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def _formatDelta(integer):
    sign_char = '+' if integer > 0 else '' # Empty string for zero and negative
    return "{0:{sign}d}".format(integer, sign=sign_char)


def processBrep(rgBrep_In, fMaxTol, fStartTol, bMarkBadGeom=True, bEcho=True, bDebug=False):
    """
    """

    if rgBrep_In.IsSolid:
        if bEcho:
            print("Brep is already closed.")
        return

    if not rgBrep_In.IsValid:
        if bEcho:
            print("Starting Brep is invalid and will not be processed.")
        return

    rgBrep_PrevWIP = None
    rgBrep_CurrentWIP = rgBrep_In.DuplicateBrep()
    fTol_Current = fStartTol
    ct_JoinedEdges_Total = 0

    Rhino.RhinoApp.SetCommandPrompt(prompt="Joining naked edges ...")

    while True:
        sc.escape_test()

        iJoinCt = rgBrep_CurrentWIP.JoinNakedEdges(fTol_Current)

        if bEcho and rgBrep_CurrentWIP.IsValid:
            print("{} edges were joined JoinNakedEdges(tolerance={}).".format(
                iJoinCt,
                _formatDistance(fTol_Current)))

        if not rgBrep_CurrentWIP.IsValid and (bEcho or bMarkBadGeom):
            _report_invalid_Brep_after_JoinNakedEdges(rgBrep_CurrentWIP, bEcho, bMarkBadGeom)
            rgBrep_CurrentWIP.Dispose()
            print("{} edges were joined JoinNakedEdges(tolerance={:f}) to an invalid brep.".format(
                iJoinCt,
                fTol_Current))
            if ct_JoinedEdges_Total:
                return rgBrep_PrevWIP # Will be None if no JoinNakedEdges were successful and didn't create an invalid brep.
            else:
                rgBrep_PrevWIP.Dispose()
                return

        ct_JoinedEdges_Total += iJoinCt

        if rgBrep_PrevWIP:
            rgBrep_PrevWIP.Dispose()

        if fTol_Current == fMaxTol:
            break

        rgBrep_PrevWIP = rgBrep_CurrentWIP
        rgBrep_CurrentWIP = rgBrep_CurrentWIP.DuplicateBrep()

        fTol_Current *= 10.0
        #sEval = "fMaxTol, fTol_Current, fMaxTol - fTol_Current"; print(sEval,'=',eval(sEval))
        #sEval = "abs(fMaxTol - fTol_Current)/max((fMaxTol, fTol_Current))"; print(sEval,'=',eval(sEval))
        if fTol_Current > fMaxTol:
            fTol_Current = fMaxTol
        elif ((fMaxTol - fTol_Current) / fMaxTol) <= 0.1:
            fTol_Current = fMaxTol

    if ct_JoinedEdges_Total == 0:
        rgBrep_CurrentWIP.Dispose()
        return

    return rgBrep_CurrentWIP


def coerceBrepObject(rhBrep_In):
    if isinstance(rhBrep_In, rd.BrepObject):
        return rhBrep_In
    elif isinstance(rhBrep_In, Guid):
        return sc.doc.Objects.FindId(rhBrep_In)
    else:
        raise Exception("{} passed to function. Should be BrepObject or GUID of one).".format(
            rdBrep_In.GetType().Name))


def processBrepObject(rhBrep_In, fMaxTol=None, fStartTol=None, bMarkBadGeom=True, bEcho=None, bDebug=None):
    """
    rhBrep_In: BrepObject or GUID of one
    """

    if fMaxTol is None: fMaxTol = Opts.values['fMaxTol']
    if fStartTol is None: fStartTol = Opts.values['fStartTol']
    if bMarkBadGeom is None: bMarkBadGeom = Opts.values['bMarkBadGeom']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']


    rdBrep_In = coerceBrepObject(rhBrep_In)
    rgBrep_In = rdBrep_In.BrepGeometry

    rgBrep_Res = processBrep(
        rgBrep_In,
        fMaxTol=fMaxTol,
        fStartTol=fStartTol,
        bMarkBadGeom=bMarkBadGeom,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if rgBrep_Res is None:
        return

    # Report face count and edge valence counts.
    ct_Faces_In = rgBrep_In.Faces.Count
    cts_of_EdgeValences_In = _countsOfEdgeValences(rgBrep_In)

    print("Counts:",
          "{} faces,".format(ct_Faces_In),
          "{} none edges,".format(cts_of_EdgeValences_In[0]),
          "{} naked edges,".format(cts_of_EdgeValences_In[1]),
          "{} interior edges,".format(cts_of_EdgeValences_In[2]),
          "{} non-manifold edges".format(cts_of_EdgeValences_In[3]))

    if not sc.doc.Objects.Replace(rdBrep_In.Id, rgBrep_Res):
        print("Brep was not replaced.")
        return False

    sc.doc.Views.Redraw()
    ct_Faces_Res = rgBrep_Res.Faces.Count
    cts_of_EdgeValences_Res = _countsOfEdgeValences(rgBrep_Res)

    print("Counts: Before -> After (Delta)")

    if ct_Faces_In != ct_Faces_Res:
        print("Faces: {} - {} ({})".format(
            ct_Faces_In,
            ct_Faces_Res,
            _formatDelta(ct_Faces_Res-ct_Faces_In),
            ))

    if cts_of_EdgeValences_In[0] or cts_of_EdgeValences_Res[0]:
        print("None edges: {} -> {} ({})".format(
            cts_of_EdgeValences_In[0],
            cts_of_EdgeValences_Res[0],
            _formatDelta(cts_of_EdgeValences_Res[0]-cts_of_EdgeValences_In[0]),
            ))

    print("Naked edges: {} -> {} ({})".format(
        cts_of_EdgeValences_In[1],
        cts_of_EdgeValences_Res[1],
        _formatDelta(cts_of_EdgeValences_Res[1]-cts_of_EdgeValences_In[1]),
        ))

    print("Interior edges: {} -> {} ({})".format(
        cts_of_EdgeValences_In[2],
        cts_of_EdgeValences_Res[2],
        _formatDelta(cts_of_EdgeValences_Res[2]-cts_of_EdgeValences_In[2]),
        ))

    if cts_of_EdgeValences_In[3] or cts_of_EdgeValences_Res[3]:
        print("Non-manifold edges: {} -> {} ({})".format(
            cts_of_EdgeValences_In[3],
            cts_of_EdgeValences_Res[3],
            _formatDelta(cts_of_EdgeValences_Res[3]-cts_of_EdgeValences_In[3]),
            ))

    return rdBrep_In


    #print("No edges were joined within a deviation of {:f}.".format(fMaxTol))
    return iJoinCt_Total


def main():
    rhObjs_In = getInput()
    if rhObjs_In is None: return
    
    fMaxTol = Opts.values['fMaxTol']
    bIterateToMaxTol = Opts.values['bIterateToMaxTol']
    fStartTol = Opts.values['fStartTol']
    bMarkBadGeom = Opts.values['bMarkBadGeom']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")

    gBs_Res = []

    for rhObj in rhObjs_In:
        gBrep = rs.coerceguid(rhObj)

        rc = processBrepObject(
            gBrep,
            fMaxTol=fMaxTol,
            fStartTol=fStartTol if bIterateToMaxTol else fMaxTol,
            bMarkBadGeom=bMarkBadGeom,
            bEcho=bEcho,
            bDebug=bDebug)

        if rc is None:
            continue
        elif not rc:
            print("No error in brep geometry, but no edges were joined.")
        else:
            print("Edges were joined.")
            gBs_Res.append(rc)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()