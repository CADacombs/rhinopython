#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
Entry point for the spb_EndBulge tool suite.
Routes to the curve or surface solver based on user selection.
See Crv or Srf scripts for more descriptions and notes.

This script was partially developed by Google Gemini 3.1 Pro.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
260720: Created.
260721: Updated noted.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import spb_EndBulge_Kernel as ebk
import spb_EndBulge_Crv as eb_crv
import spb_EndBulge_Srf as eb_srf


def main():
    go = ri.Custom.GetObject()


    def master_filter(rdObj, geom, compIdx):
        # Surface logic
        if isinstance(geom, rg.BrepTrim):
            if bEdge_forCrv_notSrf:
                return False
            rgB = geom.Brep
            if rgB.Faces.Count == 1 and rgB.Faces[0].IsSurface:
                return True

        # Curve logic
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


    go.SetCustomGeometryFilter(master_filter)
    go.DisablePreSelect()
    go.AcceptNumber(True, acceptZero=True)
    
    bEdge_forCrv_notSrf = False
    riOpt_Edge_forCrv_notSrf = ri.Custom.OptionToggle(bEdge_forCrv_notSrf, "Srf", "Crv")
    if "bEdge_forCrv_notSrf" in sc.sticky:
        bEdge_forCrv_notSrf = riOpt_Edge_forCrv_notSrf.CurrentValue = sc.sticky["bEdge_forCrv_notSrf"]
        
    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = ebk.Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opt.clear()

        if bEdge_forCrv_notSrf:
            go.SetCommandPrompt("Select curve to adjust")
            go.GeometryFilter = rd.ObjectType.Curve
            go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.AcceptAllAttributes
        else:
            go.SetCommandPrompt("Select curve or surface edge to adjust")
            go.GeometryFilter = rd.ObjectType.Curve | rd.ObjectType.EdgeFilter
            go.GeometryAttributeFilter = (
                ri.Custom.GeometryAttributeFilter.SurfaceBoundaryEdge |
                ri.Custom.GeometryAttributeFilter.SeamEdge
                )

        
        addOption('bDialog')
        idx_Edge_forCrv_notSrf = go.AddOptionToggle("Edge", riOpt_Edge_forCrv_notSrf)[0]

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
            objref_In = go.Object(0)
            go.Dispose()
            break
            
        if res == ri.GetResult.Number:
            key = 'fScale_Picked'
            ebk.Opts.riOpts[key].CurrentValue = go.Number()
            ebk.Opts.setValue(key)
            continue
            
        if res == ri.GetResult.Option:
            idx = go.Option().Index
            if idx == idx_Edge_forCrv_notSrf:
                bEdge_forCrv_notSrf = sc.sticky["bEdge_forCrv_notSrf"] = riOpt_Edge_forCrv_notSrf.CurrentValue
            else:
                for key in idxs_Opt:
                    if idx == idxs_Opt[key]:
                        ebk.Opts.setValue(key, go.Option().CurrentListOptionIndex)
                        break

    if objref_In.Trim() and not bEdge_forCrv_notSrf:
        eb_srf.main(objref_In)
    else:
        eb_crv.main(objref_In)


if __name__ == '__main__': main()
