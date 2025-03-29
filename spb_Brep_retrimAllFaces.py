"""
This script uses Brep.CutUpSurface to attempt to retrim all faces of a brep
in order to simplify edges and BrepTrim curves and to identify problem areas.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250327-28: Created

TODO:
    Try tighter tolerances on trim fails.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid


def createBrepOfUntrimmedFace(rgFace_In):
    rgB_Out = rgFace_In.DuplicateSurface().ToBrep()
    rgB_Out.Faces[0].PerFaceColor = rgFace_In.PerFaceColor
    return rgB_Out


def cutUpSurface(surface=None, curves=None, useEdgeCurves=None, tolerance=None, iFs_TrimFail=None, iFs_InvalidRes=None):
    """
    Parameters:
    Returns:
    """
    
    breps_Res = rg.Brep.CutUpSurface(
        surface=surface,
        curves=curves,
        useEdgeCurves=useEdgeCurves,
        tolerance=tolerance)

    if len(breps_Res) == 0:
        #raise Exception("No result of CutUpSurface for duplicated face.")
        iFs_TrimFail.append(iF)
        return

    if len(breps_Res) > 1:
        [sc.doc.Objects.AddBrep(b) for b in breps_Res]
        print("Multiple breps from CutUpSurface. They were added to model.")
        #sc.doc.Views.Redraw()
        iFs_TrimFail.append(iF)
        return

    rgB_1F_Out = breps_Res[0]

    if not rgB_1F_Out.IsValid:
        rgB_1F_Out.Dispose()
        iFs_InvalidRes.append(iF)
        sLogs.append(sLog)
        return

    return rgB_1F_Out


def ncptct(crv):
    nc = crv.ToNurbsCurve()
    ct = nc.Points.Count
    nc.Dispose()
    return ct


def processBrepObject(rdBrep_In, tolerance=None):
    if rdBrep_In.BrepGeometry.Faces.Count == 1:
        rdBs_1F = [rdBrep_In]
    else:
        rdBs_1F = rdBrep_In.GetSubObjects()
        if not rdBs_1F:
            print("No GetSubObjects result for {}.".format(rdBrep_In.Id))
            return
    
    iFs_TrimFail = []
    iFs_InvalidRes = []
    sLogs = []
    iFs_Success = []
    rgBs_1F_Res = []
    nCPs_Es_Pre = nCPs_Es_Post = nCPs_Ts_Pre = nCPs_Ts_Post = 0
    
    rgBs_1F_DoNotJoin = []
    rgBs_toManuallyTrim = []
    
    for iF, rdB_1F in enumerate(rdBs_1F):
        Rhino.RhinoApp.Wait()
        Rhino.RhinoApp.SetCommandPromptMessage(
            "Processing face {} of {}...".format(iF+1, len(rdBs_1F)))

        rgB_1F_In = rdB_1F.Geometry
        rgF_In = rgB_1F_In.Faces[0]
        rgS_In = rgF_In.UnderlyingSurface()
        
        area_In = rgB_1F_In.GetArea()
        if not area_In:
            rgBs_toManuallyTrim.append(createBrepOfUntrimmedFace(rgF_In))
            rdB_1F.Dispose()
            continue
        
        rgCs_In = [rgB_1F_In.Edges[iE].DuplicateCurve() for iE in rgF_In.AdjacentEdges()]
        #[sc.doc.Objects.AddCurve(c) for c in rgCs_In]

        rv = cutUpSurface(
            surface=rgS_In,
            curves=rgCs_In,
            useEdgeCurves=False,
            tolerance=tolerance,
            iFs_TrimFail=iFs_TrimFail,
            iFs_InvalidRes=iFs_InvalidRes,
            )

        if rv is None:
            continue

        # Success, but check area first.
        rgB_1F_Res_A = rv
        
        area_Phase1 = rgB_1F_Res_A.GetArea()
        if not area_Phase1:
            rgBs_toManuallyTrim.append(createBrepOfUntrimmedFace(rgF_In))
            rgB_1F_Res_A.Dispose()
            rdB_1F.Dispose()
            continue

        relDiff = (abs(area_In-area_Phase1)/max((area_In, area_Phase1)))
        if relDiff > 0.01:
            # Shrink surface and try again.
            bShrinkFace = rgF_In.ShrinkFace(
                disableSide=rg.BrepFace.ShrinkDisableSide.ShrinkAllSides)

            if not bShrinkFace:
                rgBs_toManuallyTrim.append(createBrepOfUntrimmedFace(rgF_In))
                rgB_1F_Res_A.Dispose()
                rdB_1F.Dispose()
                continue

            rv = cutUpSurface(
                surface=rgS_In,
                curves=rgCs_In,
                useEdgeCurves=False,
                tolerance=tolerance,
                iFs_TrimFail=iFs_TrimFail,
                iFs_InvalidRes=iFs_InvalidRes,
                )
            
            if rv is None:
                continue
            
            # Success, but check area first.
            rgB_1F_Res_A = rv
            
            area_Phase1 = rgB_1F_Res_A.GetArea()
            if not area_Phase1:
                rgBs_toManuallyTrim.append(createBrepOfUntrimmedFace(rgF_In))
                rgB_1F_Res_A.Dispose()
                rdB_1F.Dispose()
                continue

            relDiff = (abs(area_In-area_Phase1)/max((area_In, area_Phase1)))
            if relDiff > 0.01:
                rgBs_toManuallyTrim.append(createBrepOfUntrimmedFace(rgF_In))
                rgB_1F_Res_A.Dispose()
                rdB_1F.Dispose()
                continue

        # Success, so repeat on result to simplify BrepTrim curves.
        rgF_Res = rgB_1F_Res_A.Faces[0]
        rgCs_Res = [rgB_1F_Res_A.Edges[iE].DuplicateCurve() for iE in rgF_Res.AdjacentEdges()]
        
        rv = cutUpSurface(
            surface=rgS_In,
            curves=rgCs_Res,
            useEdgeCurves=False,
            tolerance=tolerance,
            iFs_TrimFail=iFs_TrimFail,
            iFs_InvalidRes=iFs_InvalidRes,)

        if rv is None:
            continue

        rgB_1F_Res = rv

        if rgF_In.PerFaceColor is not None:
            rgB_1F_Res.Faces[0].PerFaceColor = rgF_In.PerFaceColor

        rgBs_1F_Res.append(rgB_1F_Res)
        iFs_Success.append(iF)

        nCPs_Es_Pre += sum([ncptct(rgE) for rgE in rgB_1F_In.Edges])
        nCPs_Es_Post += sum([ncptct(rgE) for rgE in rgB_1F_Res.Edges])
        nCPs_Ts_Pre += sum([ncptct(rgT) for rgT in rgB_1F_In.Trims])
        nCPs_Ts_Post += sum([ncptct(rgT) for rgT in rgB_1F_Res.Trims])

    if len(rgBs_1F_Res) == 0:
        pass
    elif len(rgBs_1F_Res) == 1:
        sc.doc.Objects.AddBrep(rgBs_1F_Res[0])
    else:
        breps_Joined = rg.Brep.JoinBreps(
            brepsToJoin=rgBs_1F_Res,
            tolerance=2.0*tolerance)
        
        if len(breps_Joined) == 0:
            raise Exception("No result from JoinBreps.")
        if len(breps_Joined) > 1:
            [sc.doc.Objects.AddBrep(b) for b in breps_Joined]
            print("Multiple breps output of JoinBreps.")
            gB_Outs = [
                sc.doc.Objects.AddBrep(
                    breps_Joined[0], attributes=rdBrep_In.Attributes)
                    for brep_Joined in breps_Joined]
        else:
            sc.doc.Objects.Replace(
                objectId=rdBrep_In.Id,
                brep=breps_Joined[0])
    
    gB_toManuallyTrim = [
        sc.doc.Objects.AddBrep(
            rgB, attributes=rdBrep_In.Attributes) for rgB in rgBs_toManuallyTrim]
    if gB_toManuallyTrim:
        if not all(gB==Guid.Empty for gB in gB_toManuallyTrim):
            print("Added {} untrimmed faces. Trim this manually.".format(
                len(gB_toManuallyTrim) - gB_toManuallyTrim.count(Guid.Empty)))
        if Guid.Empty in gB_toManuallyTrim:
            print("{} surfaces could not be added. Check results.".format(
                gB_toManuallyTrim.count(Guid.Empty)))

    print("Retrimmed {} out of {} faces.  {}".format(
        len(iFs_Success), len(rdBs_1F), iFs_Success[:10]))
    print("Trim fail count: {}".format(len(iFs_TrimFail)))
    print("Invalid brep result (skipped) count: {}".format(len(iFs_InvalidRes)))
    if sLogs:
        print("Logs of invalid breps:")
        print(*sLogs)

    print("Edge CP count: {} -> {}, {:+}.".format(
        nCPs_Es_Pre,
        nCPs_Es_Post,
        nCPs_Es_Post-nCPs_Es_Pre))

    print("Trim CP count: {} -> {}, {:+}.".format(
        nCPs_Ts_Pre,
        nCPs_Ts_Post,
        nCPs_Ts_Post-nCPs_Ts_Pre))


def main():
    
    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select breps",
        acceptNothing=False,
        filter=rd.ObjectType.Brep)
    if res != Rhino.Commands.Result.Success: return
    
    tolerance = sc.doc.ModelAbsoluteTolerance
    
    for objref in objrefs:
        rdBrep_In = objref.Object()
        processBrepObject(
            rdBrep_In,
            tolerance=tolerance)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()