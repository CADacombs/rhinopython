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

This script was modified by Google Gemini 3.1 Pro.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
210303, 0307: Created.
260420-25, 0709-14: Added an optional dialog. Added a preview for the dialog.
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

        addOption('bGUI')
        if not ebk.Opts.values['bGUI']:
            addOption('bLinkedEnds')
            addOption('fScale_Picked')
            addOption('fSlideG2_Picked')
            addOption('fSlideG3_Picked')
            addOption('fScale_Opp')
            addOption('fSlideG2_Opp')
            addOption('fSlideG3_Opp')
            addOption('idxCont_Picked')
            addOption('idxCont_Opp')
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

        if t_AtPicked > nc_In.Domain.Mid:
            iPickedEnd = 1
            fScale_T1, fSlideG2_T1, fSlideG3_T1 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked
            fScale_T0, fSlideG2_T0, fSlideG3_T0 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp
            iG_T1, iG_T0 = idxCont_Picked - 1, idxCont_Opp - 1
        else:
            iPickedEnd = 0
            fScale_T0, fSlideG2_T0, fSlideG3_T0 = fScale_Picked, fSlideG2_Picked, fSlideG3_Picked
            fScale_T1, fSlideG2_T1, fSlideG3_T1 = fScale_Opp, fSlideG2_Opp, fSlideG3_Opp
            iG_T0, iG_T1 = idxCont_Picked - 1, idxCont_Opp - 1

        can_G3_T0 = ebk.canMaintainG3(nc_In, False)
        can_G3_T1 = ebk.canMaintainG3(nc_In, True)

        if not can_G3_T0 and iG_T0 == 3:
            iG_T0 = 2
            if bEcho: print("T0 continuity downgraded to G2. G3 requires internal knot multiplicity >= 3.")
        if not can_G3_T1 and iG_T1 == 3:
            iG_T1 = 2
            if bEcho: print("T1 continuity downgraded to G2. G3 requires internal knot multiplicity >= 3.")

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
            if bEcho: print("Curve could not be created. {}".format(sReport))
            return None

    if nc_Res.EpsilonEquals(objref_In.Curve().ToNurbsCurve(), epsilon=Rhino.RhinoMath.ZeroTolerance):
        print("Resultant curve is the same as input curve. No changes were made to the document.")
        return

    if not bDeleteInput or objref_In.Edge():
        gC_Out = sc.doc.Objects.AddCurve(nc_Res)
        if gC_Out == gC_Out.Empty:
            if bEcho: print("Could not add curve.")
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

    bGUI = ebk.Opts.values['bGUI']
    if not bGUI:
        nc_Res = None
    else:
        Rhino.RhinoApp.SetCommandPromptMessage("Continuing in dialog...")
        nc_Res = _createCurve_viaGUI(objref_In)
        if nc_Res is None: return

    # Extract parsed variables from the Kernel's Opts dictionary
    bLinkedEnds = ebk.Opts.values['bLinkedEnds']
    fScale_Picked = ebk.Opts.values['fScale_Picked']
    fSlideG2_Picked = ebk.Opts.values['fSlideG2_Picked'] 
    fSlideG3_Picked = ebk.Opts.values['fSlideG3_Picked'] 
    fScale_Opp = ebk.Opts.values['fScale_Opp']
    fSlideG2_Opp = ebk.Opts.values['fSlideG2_Opp']
    fSlideG3_Opp = ebk.Opts.values['fSlideG3_Opp']
    
    idxCont_Picked = ebk.Opts.values['idxCont_Picked']
    idxCont_Opp = ebk.Opts.values['idxCont_Opp']
    bDeleteInput = ebk.Opts.values['bDeleteInput']
    bEcho = ebk.Opts.values['bEcho']
    bDebug = ebk.Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False
    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObject(
        objref_In=objref_In,
        nc_Precalc=nc_Res,
        bLinkedEnds=bLinkedEnds,
        fScale_Picked=fScale_Picked, fSlideG2_Picked=fSlideG2_Picked, fSlideG3_Picked=fSlideG3_Picked,
        fScale_Opp=fScale_Opp, fSlideG2_Opp=fSlideG2_Opp, fSlideG3_Opp=fSlideG3_Opp,
        idxCont_Picked=idxCont_Picked,
        idxCont_Opp=idxCont_Opp,
        bDeleteInput=bDeleteInput, bEcho=bEcho, bDebug=bDebug
    )

    if gC_Res is None: return
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()