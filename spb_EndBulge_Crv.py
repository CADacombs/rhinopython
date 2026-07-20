#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script is an alternative to _EndBulge.

TL;DR: The command should be used to really understand it. Compare with _EndBulge
and see.

Unlike _EndBulge, there is no dragging of control points in the graphics window.
Instead, number steppers in a dialog are used to specify:
    1. The tangent vector (p1 - p0) scale relative to its starting position.
    2. Where the geometry allows it, and after #1 is applied,
    the G2 (p2) tangential sliding scale from p2's starting position.
    3. Where the geometry allows it, and after #1 & #2 are applied,
    the G3 (p3) tangential sliding scale from p3's starting position.
Unlike _EndBulge, both the picked end and the opposite end can be modified.
Unlike _EndBulge, a curvature graph is automatically included.
Unlike _EndBulge, the continuities to maintain for the picked end and the opposite
end are explicitly defined and selectable by the user.

There are more options, but the command should be used to really understand it.

This script was partially developed by Google Gemini 3.1 Pro.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
210303, 0307: Created.
260420-25, 0709-18: Added an optional dialog. Added a preview for the dialog.
        Refactored.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import spb_EndBulge_Kernel as ebk


def getInput_CLI():
    """
    Get curve with picked end and optional input.
    """
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Pick curve near an end")
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    def geomFilter_Curve(rdObj, geom, compIdx):
        if isinstance(geom, rg.BrepEdge):
            rgC = geom.DuplicateCurve()
        elif isinstance(geom, rg.Curve):
            rgC = geom
        else:
            return False

        if isinstance(rgC, rg.NurbsCurve):
            nc = rgC
        elif isinstance(rgC, rg.PolyCurve):
            if rgC.SegmentCount > 1:
                print("PolyCurve with multiple segments is ignored.")
                return False
            nc = rgC.ToNurbsCurve()
        else:
            return False

        if nc.IsPeriodic:
            print("Periodic curves are not supported.")
            return False
        if nc.Degree == 1:
            print("Ignored degree 1 NURBS curve.")
            return False
        return True

    go.SetCustomGeometryFilter(geomFilter_Curve)
    go.DisablePreSelect()
    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = ebk.Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        addOption('bDialog')
        if not ebk.Opts.values['bDialog']:
            addOption('idxCont_Picked')
            addOption('idxCont_Opp')
            addOption('bLinkedEnds')
            if ebk.Opts.values['bLinkedEnds']:
                ebk.Opts.names['fScale_Picked'] = 'Scale'
                ebk.Opts.names['fSlideG2_Picked'] = 'SlideG2'
                ebk.Opts.names['fSlideG3_Picked'] = 'SlideG3'
            else:
                ebk.Opts.names['fScale_Picked'] = 'Scale_Picked'
                ebk.Opts.names['fSlideG2_Picked'] = 'SlideG2Picked'
                ebk.Opts.names['fSlideG3_Picked'] = 'SlideG3Picked'
                ebk.Opts.names['fScale_Opp'] = 'fScaleOpp'
                ebk.Opts.names['fSlideG2_Opp'] = 'fSlideG2Opp'
                ebk.Opts.names['fSlideG3_Opp'] = 'fSlideG3Opp'
            addOption('fScale_Picked')
            addOption('fSlideG2_Picked')
            addOption('fSlideG3_Picked')
            if not ebk.Opts.values['bLinkedEnds']:
                addOption('fScale_Opp')
                addOption('fSlideG2_Opp')
                addOption('fSlideG3_Opp')
            addOption('bDeleteInput')
            addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return objref

        if res == ri.GetResult.Number:
            key = 'fScale_Picked'
            ebk.Opts.riOpts[key].CurrentValue = go.Number()
            ebk.Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                ebk.Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _createCurve_viaGUI(objref_In):
    key = 'conduit'
    stickyKey = '{}({})'.format(key, __file__)
    if stickyKey in sc.sticky:
        old_conduit = sc.sticky[stickyKey]
        old_conduit.Enabled = False
        sc.doc.Views.Redraw()

    parent = Rhino.UI.RhinoEtoApp.MainWindowForDocument(sc.doc)
    dialog = ebk.EtoDialog(objref_In)

    dialog.conduit = ebk.EndBulgePreviewConduit()
    sc.sticky[stickyKey] = dialog.conduit

    # Run an initial preview update before opening the dialog
    dialog.UpdatePreview()
    dialog.conduit.Enabled = True
    sc.doc.Views.Redraw()

    Rhino.UI.EtoExtensions.ShowSemiModal(dialog, sc.doc, parent)

    dialog.conduit.Enabled = False
    sc.doc.Views.Redraw()

    if not dialog.dialog_ok:
        return None

    return dialog.conduit.crv


def processCurveObject(objref_In, nc_Precalc=None, **kwargs):
    def getOpt(key): return kwargs[key] if key in kwargs else ebk.Opts.values[key]

    bLinkedEnds = getOpt('bLinkedEnds')
    fScale_Picked = getOpt('fScale_Picked')
    fSlideG2_Picked = getOpt('fSlideG2_Picked') 
    fSlideG3_Picked = getOpt('fSlideG3_Picked') 
    fScale_Opp = getOpt('fScale_Opp')
    fSlideG2_Opp = getOpt('fSlideG2_Opp')
    fSlideG3_Opp = getOpt('fSlideG3_Opp')
    idxCont_Picked = getOpt('idxCont_Picked')
    idxCont_Opp = getOpt('idxCont_Opp')
    bDeleteInput = getOpt('bDeleteInput')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if nc_Precalc is not None:
        nc_Res = nc_Precalc
    else:
        # --- CLI EXECUTION PATH ---
        rgC_In = objref_In.Curve()

        if isinstance(rgC_In, rg.BrepEdge):
            nc_In = rgC_In.ToNurbsCurve()
        elif isinstance(rgC_In, rg.NurbsCurve):
            nc_In = rgC_In
        elif isinstance(rgC_In, rg.PolyCurve):
            nc_In = rgC_In.ToNurbsCurve()
        else:
            return None

        bSuccess, t_AtPicked = nc_In.ClosestPoint(objref_In.SelectionPoint())
        if not bSuccess: return None

        # Mirror linked ends strictly for CLI
        if bLinkedEnds:
            fScale_Opp = fScale_Picked
            fSlideG2_Opp = fSlideG2_Picked
            fSlideG3_Opp = fSlideG3_Picked

        if t_AtPicked > nc_In.Domain.Mid:
            iPickedEnd = 1
            fScale_T1, fSlideG2_T1, fSlideG3_T1 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked
            fScale_T0, fSlideG2_T0, fSlideG3_T0 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp
            iG_T1, iG_T0 = idxCont_Picked - 1, idxCont_Opp - 1
            name_T1, name_T0 = "Picked end", "Opposite end"
        else:
            iPickedEnd = 0
            fScale_T0, fSlideG2_T0, fSlideG3_T0 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked
            fScale_T1, fSlideG2_T1, fSlideG3_T1 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp
            iG_T0, iG_T1 = idxCont_Picked - 1, idxCont_Opp - 1
            name_T0, name_T1 = "Picked end", "Opposite end"

        # --- STRICT CLI VALIDATION ---
        errors = []
        N = nc_In.Points.Count
        req_0 = max(0, iG_T0 + 1)
        req_1 = max(0, iG_T1 + 1)

        if req_0 + req_1 > N:
            errors.append("Continuity constraints overlap (not enough control points).")

        if iG_T0 == 3 and not ebk.canMaintainG3(nc_In, False):
            errors.append("{} continuity cannot be G3 (requires internal knot multiplicity >= 3).".format(name_T0))
        if iG_T1 == 3 and not ebk.canMaintainG3(nc_In, True):
            errors.append("{} continuity cannot be G3 (requires internal knot multiplicity >= 3).".format(name_T1))

        if not errors:
            alloc_0 = req_0
            alloc_1 = req_1
            free = N - alloc_0 - alloc_1
            if free > 0:
                half = free // 2
                extra = free % 2
                scale_limit_T0 = alloc_0 + half
                scale_limit_T1 = alloc_1 + half
                if extra:
                    if iPickedEnd == 0: scale_limit_T0 += 1
                    else: scale_limit_T1 += 1
            else:
                scale_limit_T0 = alloc_0
                scale_limit_T1 = alloc_1

            tol = Rhino.RhinoMath.ZeroTolerance
            if abs(fSlideG2_T0) > tol and scale_limit_T0 < 3:
                errors.append("{} G2 slide is not allowed due to point count constraints.".format(name_T0))
            if abs(fSlideG3_T0) > tol and scale_limit_T0 < 4:
                errors.append("{} G3 slide is not allowed due to point count constraints.".format(name_T0))
            if abs(fSlideG2_T1) > tol and scale_limit_T1 < 3:
                errors.append("{} G2 slide is not allowed due to point count constraints.".format(name_T1))
            if abs(fSlideG3_T1) > tol and scale_limit_T1 < 4:
                errors.append("{} G3 slide is not allowed due to point count constraints.".format(name_T1))

        # Abort if any violations occurred
        if errors:
            if bEcho:
                print("CLI Error: Invalid settings applied. No modifications made.")
                for err in errors:
                    print(" - " + err)
            return None

        # Execute creation if validation passes
        nc_Res, sReport, info = ebk.createCurve(
            nc_In=nc_In,
            fScale_T0=fScale_T0, fSlideG2_T0=fSlideG2_T0, fSlideG3_T0=fSlideG3_T0,
            fScale_T1=fScale_T1, fSlideG2_T1=fSlideG2_T1, fSlideG3_T1=fSlideG3_T1,
            iG_T0=iG_T0,
            iG_T1=iG_T1,
            iPickedEnd=iPickedEnd,
            bDebug=bDebug
        )

        if nc_Res is None:
            if bEcho: print("Curve was not generated. {}".format(sReport))
            return None

    if nc_Res.EpsilonEquals(objref_In.Curve().ToNurbsCurve(), epsilon=Rhino.RhinoMath.ZeroTolerance):
        print("Resultant curve is the same as input curve. No changes were made to the document.")
        return

    # Failsafe variable assignment
    gC_Out = None
    
    if not bDeleteInput or objref_In.Edge():
        gC_Out = sc.doc.Objects.AddCurve(nc_Res)
        if gC_Out == gC_Out.Empty:
            if bEcho: print("Could not add curve.")
            gC_Out = None
        else:
            if bEcho: print("Curve was added.")
    else:
        if sc.doc.Objects.Replace(objref_In.ObjectId, nc_Res):
            gC_Out = objref_In.ObjectId
            if bEcho: print("Replaced curve.")
        else:
            if bEcho: print("Could not replace curve.")

    return gC_Out


def main():
    rv = getInput_CLI()
    if rv is None: return
    objref_In = rv

    bDialog = ebk.Opts.values['bDialog']

    if not bDialog:
        # --- CLI EXECUTION PATH ---
        sc.doc.Objects.UnselectAll()
        gC_Res = processCurveObject(
            objref_In=objref_In,
            nc_Precalc=None,
            bLinkedEnds=ebk.Opts.values['bLinkedEnds'],
            fScale_Picked=ebk.Opts.values['fScale_Picked'], fSlideG2_Picked=ebk.Opts.values['fSlideG2_Picked'], fSlideG3_Picked=ebk.Opts.values['fSlideG3_Picked'],
            fScale_Opp=ebk.Opts.values['fScale_Opp'], fSlideG2_Opp=ebk.Opts.values['fSlideG2_Opp'], fSlideG3_Opp=ebk.Opts.values['fSlideG3_Opp'],
            idxCont_Picked=ebk.Opts.values['idxCont_Picked'],
            idxCont_Opp=ebk.Opts.values['idxCont_Opp'],
            bDeleteInput=ebk.Opts.values['bDeleteInput'], bEcho=ebk.Opts.values['bEcho'], bDebug=ebk.Opts.values['bDebug']
        )
        if gC_Res is not None: sc.doc.Views.Redraw()
        return

    # --- GUI EXECUTION PATH ---
    Rhino.RhinoApp.SetCommandPromptMessage("Continuing in dialog...")
    key = 'conduit_crv'
    stickyKey = '{}({})'.format(key, __file__)
    if stickyKey in sc.sticky:
        sc.sticky[stickyKey].Enabled = False

    parent = Rhino.UI.RhinoEtoApp.MainWindowForDocument(sc.doc)
    dialog = ebk.EtoDialog(objref_In)
    dialog.conduit = ebk.EndBulgePreviewConduit()
    sc.sticky[stickyKey] = dialog.conduit

    # Lock the original object so it acts as a visual reference during the command
    sc.doc.Objects.Lock(objref_In.ObjectId, True)

    dialog.UpdatePreview()
    dialog.conduit.Enabled = True
    sc.doc.Views.Redraw()

    # START UNDO RECORD
    undo_sn = sc.doc.BeginUndoRecord("EndBulge Crv")

    try:
        Rhino.UI.EtoExtensions.ShowSemiModal(dialog, sc.doc, parent)

        # Suspend redraws immediately after dialog closes to prevent unlocking flicker
        if not ebk.Opts.values['bDebug']: sc.doc.Views.RedrawEnabled = False
        sc.doc.Objects.UnselectAll()

        # UNLOCK HERE: Failsafe so doc.Objects.Replace() succeeds without error
        sc.doc.Objects.Unlock(objref_In.ObjectId, True)

        if dialog.dialog_ok and dialog.conduit.crv:
            # Execute the final curve replacement with the freshest variables from the dialog
            gC_Res = processCurveObject(
                objref_In=objref_In,
                nc_Precalc=dialog.conduit.crv,
                bLinkedEnds=ebk.Opts.values['bLinkedEnds'],
                fScale_Picked=ebk.Opts.values['fScale_Picked'], fSlideG2_Picked=ebk.Opts.values['fSlideG2_Picked'], fSlideG3_Picked=ebk.Opts.values['fSlideG3_Picked'],
                fScale_Opp=ebk.Opts.values['fScale_Opp'], fSlideG2_Opp=ebk.Opts.values['fSlideG2_Opp'], fSlideG3_Opp=ebk.Opts.values['fSlideG3_Opp'],
                idxCont_Picked=ebk.Opts.values['idxCont_Picked'],
                idxCont_Opp=ebk.Opts.values['idxCont_Opp'],
                bDeleteInput=ebk.Opts.values['bDeleteInput'], bEcho=ebk.Opts.values['bEcho'], bDebug=ebk.Opts.values['bDebug']
            )

    except Exception as e:
        import traceback
        print("Script Error Encountered: {}".format(e))
        print("Standard Traceback:")
        print(traceback.format_exc())
        
    finally:
        # Cleanly disable the conduit and re-enable redrawing
        dialog.conduit.Enabled = False
        sc.doc.Objects.Unlock(objref_In.ObjectId, True) # Secondary failsafe just in case
        sc.doc.EndUndoRecord(undo_sn)
        sc.doc.Views.RedrawEnabled = True
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()