"""
This script will select grips interior of those already selected in NURBS curves and surfaces.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250129: Created.

Apparently, overloads for GripObject.Select for syncHighlight and persistentSelect are not needed
for both to be true after script completes. Was there a modification to the method since V5?
"""

import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc


def _getOwnersOfSelectedGrips():
    gOwners = []
    rdOwners = []
    rdSelected = list(sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=True))
    for rdObj in rdSelected:
        if rdObj.ObjectType != rd.ObjectType.Grip:
            continue
        gOwner = rdObj.OwnerId
        if gOwner in gOwners:
            continue
        rdOwner = sc.doc.Objects.FindId(gOwner)
        if rdOwner.ObjectType not in (rd.ObjectType.Brep, rd.ObjectType.Curve):
            continue
        rgOwner = rdOwner.Geometry
        if not isinstance(rgOwner, (rg.Brep, rg.NurbsCurve)):
            continue
        if isinstance(rgOwner, rg.Brep):
            if rgOwner.Faces.Count != 1:
                continue
            if not isinstance(rgOwner.Faces[0].UnderlyingSurface(), rg.NurbsSurface):
                continue
        gOwners.append(gOwner)
        rdOwners.append(rdOwner)

    return rdOwners


def _addInteriorGripsToSelected_NurbsCrv(rdCurve):

    idxs_CPs_All = []
    rdGrips = rdCurve.GetGrips()
    for rdGrip in rdGrips:
        if rdGrip.IsSelected(checkSubObjects=False):
            iCt_Cps, idxs_CPs_Now = rdGrip.GetCurveCVIndices()
            for i in idxs_CPs_Now:
                if i not in idxs_CPs_All:
                    idxs_CPs_All.append(i)

    if len(idxs_CPs_All) < 2:
        return []

    nc = rdCurve.CurveGeometry

    rdGrips_Added = []

    for i in range(min(idxs_CPs_All)+1, max(idxs_CPs_All)):
        rdGrip = rdGrips[i]
        if not rdGrip.IsSelected(checkSubObjects=False):
            if rdGrip.Select(True):
                rdGrips_Added.append(rdGrip)

    return rdGrips_Added


def _addInteriorGripsToSelection_Brep_1Face(rdBrep):

    iVs_perU = {}
    iUs_perV = {}
    rdGrips = rdBrep.GetGrips()
    for rdGrip in rdGrips:
        if rdGrip.IsSelected(checkSubObjects=False):
            iCt_Cps, idxs_UVs = rdGrip.GetSurfaceCVIndices()
            for idxs_UV in idxs_UVs:
                iU, iV = idxs_UV
                if iU not in iVs_perU:
                    iVs_perU[iU] = [iV]
                else:
                    iVs_perU[iU].append(iV)
                if iV not in iUs_perV:
                    iUs_perV[iV] = [iU]
                else:
                    iUs_perV[iV].append(iU)

    if not iVs_perU and not iUs_perV:
        return []

    ns = rdBrep.BrepGeometry.Faces[0].UnderlyingSurface()

    rdGrips_Added = []

    for iU in iVs_perU:
        if len(iVs_perU[iU]) == 1:
            continue
        for iV in range(min(iVs_perU[iU])+1, max(iVs_perU[iU])):
            rdGrip = rdGrips[iU*ns.Points.CountV + iV]
            if not rdGrip.IsSelected(checkSubObjects=False):
                if rdGrip.Select(True):
                    rdGrips_Added.append(rdGrip)

    for iV in iUs_perV:
        if len(iUs_perV[iV]) == 1:
            continue
        for iU in range(min(iUs_perV[iV])+1, max(iUs_perV[iV])):
            rdGrip = rdGrips[iU*ns.Points.CountV + iV]
            if not rdGrip.IsSelected(checkSubObjects=False):
                if rdGrip.Select(True):
                    rdGrips_Added.append(rdGrip)

    return rdGrips_Added


def _addInteriorGripsToSelected(rdOwner):
    if isinstance(rdOwner, rd.BrepObject):
        return _addInteriorGripsToSelection_Brep_1Face(rdBrep=rdOwner)
    elif isinstance(rdOwner, rd.CurveObject):
        return _addInteriorGripsToSelected_NurbsCrv(rdCurve=rdOwner)
    else:
        raise Exception("{} not suppored input for _addInteriorGripsToSelected.".format(rdOwner))


def main():

    rdOwners = _getOwnersOfSelectedGrips()
    if not rdOwners: return
    
    rdGrips_Added_All = []

    for rdOwner in rdOwners:
        rdGrips_Added_All.extend(_addInteriorGripsToSelected(rdOwner))

    if not rdGrips_Added_All:
        print("No grips were selected.")
        return

    print("{} surface grips are selected.".format(len(rdGrips_Added_All)))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()