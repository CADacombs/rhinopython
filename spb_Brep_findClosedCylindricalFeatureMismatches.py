"""
Use this script to check for possible mismatches in closed cylindrical face sets
between breps to find modeling errors between overlapping holes or holes with
their respective fasteners.

_ExplodeBlock all instances that should be checked. Remember to _UndoMultiple
the _ExplodeBlock later. It should respect the layer state, not process invisible objects.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250203-05: Created.

TODO:
    Add support for block instances so that _ExplodeBlock does not need to be
    used before this script.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color


def _getAllNormal_breps_and_instances():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep | rd.ObjectType.InstanceReference
    return list(sc.doc.Objects.GetObjectList(oes))


def collect_cylinders_per_face_index(rgBrep, tolCyl):
    """
    Parameters:
    Returns:
    """
    cyls_perF = {}

    for rgF in rgBrep.Faces:
        #if rg.Surface.IsCylinder(
        #    rgF,
        #    tolerance=tolCyl)
        bSuccess, cyl = rg.Surface.TryGetCylinder(
            rgF,
            tolerance=tolCyl)
        if bSuccess:
            cyls_perF[rgF.FaceIndex] = cyl
            continue

        if isinstance(rgF.UnderlyingSurface(), rg.NurbsSurface):
            continue

        # Due to bug in TryGetCylinder, try one more time with the NURBS-equivalent surface.
        ns = rgF.ToNurbsSurface()
        bSuccess, cyl = rg.Surface.TryGetCylinder(
            ns,
            tolerance=tolCyl)
        if bSuccess:
            cyls_perF[rgF.FaceIndex] = cyl

    return cyls_perF


def _areCoCylindric(cylA, cylB, tolCyl=Rhino.RhinoMath.ZeroTolerance):
    if abs(cylA.Radius - cylB.Radius) > tolCyl:
        return False

    if not cylA.Axis.EpsilonEquals(cylB.Axis, epsilon=tolCyl):
        if not cylA.Axis.EpsilonEquals(-cylB.Axis, epsilon=tolCyl):
            return False

    isParallelTo = cylA.Axis.IsParallelTo(
        cylB.Axis,
        angleTolerance=0.1*sc.doc.ModelAngleToleranceRadians)

    if isParallelTo == 0:
        sEval = "cylA.Axis"; print(sEval,'=',eval(sEval))
        sEval = "cylB.Axis"; print(sEval,'=',eval(sEval))
        angle = rg.Vector3d.VectorAngle(cylA.Axis, cylB.Axis)
        #sEval = "angle"; print(sEval,'=',eval(sEval))
        #sEval = "Rhino.RhinoMath.ToDegrees(angle)"; print(sEval,'=',eval(sEval))
        #sEval = "tolCyl"; print(sEval,'=',eval(sEval))
        sc.doc.Objects.AddCircle(cylA.CircleAt(linearParameter=0))
        sc.doc.Objects.AddCircle(cylB.CircleAt(linearParameter=0))
        raise Exception("isParallelTo of {} should not be so because EpsilonEquals of previous block passed.".format(
            isParallelTo))
        return False

    cirA = cylA.CircleAt(0.0)
    planeA = cirA.Plane
    cirB = cylB.CircleAt(0.0)
    xf = rg.Transform.PlanarProjection(planeA)
    ncB = cirB.ToNurbsCurve() # So that curve can deform as needed.
    ncB.Transform(xf)
    #sc.doc.Objects.AddCurve(ncB)
    bSuccess, cirB = ncB.TryGetCircle(tolCyl)
    if not bSuccess:
        return False

    if not cirA.Center.EpsilonEquals(cirB.Center, tolCyl):
        return False

    if not Rhino.RhinoMath.EpsilonEquals(cirA.Radius, cirB.Radius, tolCyl):
        return False

    return True


    ac_A = rg.ArcCurve(cirA)
    ac_B = rg.ArcCurve(cirB)
    
    if not cirA.EpsilonEquals(cirB, tolCyl):
        rvs = rg.Curve.GetDistancesBetweenCurves(ac_A, ac_B, 0.1*tolCyl)
        if not rvs[0]:
            return False

        if rvs[1] > tolCyl:
            return False
        
        #sc.doc.Objects.AddCircle(rgCirA)
        #sc.doc.Objects.AddCircle(rgCirB)

    return True


def findContiguousCocylindricFacesOf1Face(rgBrep, cyls_perF, idx_F_Start, idxs_Fs_Cocyl_In, idxs_Fs_Fails_In, tolCyl):

    idxs_Fs_Cocyl_Out = idxs_Fs_Cocyl_In[:]
    idxs_Fs_Fails_Out = idxs_Fs_Fails_In[:]

    idxs_FAs = sorted(rgBrep.Faces[idx_F_Start].AdjacentFaces())
    #sEval = "idxs_FAs"; print(sEval,'=',eval(sEval))

    if len(idxs_FAs) == 0:
        return idxs_Fs_Cocyl_Out, idxs_Fs_Fails_Out

    cyl_This = cyls_perF[idx_F_Start]

    for iF_Adj in idxs_FAs:
        if iF_Adj not in cyls_perF:
            continue

        if iF_Adj in idxs_Fs_Cocyl_Out:
            continue

        if iF_Adj in idxs_Fs_Fails_In:
            continue

        #sEval = "iF, iF_Adj"; print(sEval,'=',eval(sEval))

        #sEval = "cyls_perF[iF].IsFinite"; print(sEval,'=',eval(sEval))

        cyl_Other = cyls_perF[iF_Adj]

        if _areCoCylindric(cyl_This, cyl_Other, tolCyl=tolCyl):
            idxs_Fs_Cocyl_Out.append(iF_Adj)
        else:
            idxs_Fs_Fails_Out.append(iF_Adj)

        #if not cyls_perF[iF].IsFinite and not cyls_perF[iF_Adj].IsFinite:
        #    if cyls_perF[iF].Axis.EpsilonEquals(
        #        cyls_perF[iF_Adj].Axis, epsilon=tolCyl
        #    ):
        #        if cyls_perF[iF].Center.EpsilonEquals(
        #            cyls_perF[iF_Adj].Center, epsilon=tolCyl
        #        ):
        #            cyls_Out.append(cyls_perF[iF])
        #            continue

    return idxs_Fs_Cocyl_Out, idxs_Fs_Fails_Out


def get_lists_of_face_indices_for_cocylindric_faces(rgBrep, cyls_perF, tolCyl):

    idxs_Fs_perCyl_Out = [] # List of lists of int indices.

    idxs_Fs_Processed = [] # Flat list.

    for i, iF in enumerate(sorted(cyls_perF)):
        #print("{} of {}".format(i+1, len(cyls_perF)))
        if i == 4:
            pass

        if iF in idxs_Fs_Processed:
            continue

        idxs_Fs_Cocyls = [iF]
        idxs_Fs_Fails = idxs_Fs_Processed

        i = 0

        while True:
            sc.escape_test()

            iF_Start = idxs_Fs_Cocyls[i]

            idxs_Fs_Cocyls_Res, idxs_Fs_Fails_Res = findContiguousCocylindricFacesOf1Face(
                rgBrep,
                cyls_perF,
                iF_Start,
                idxs_Fs_Cocyls,
                idxs_Fs_Fails,
                tolCyl)

            idxs_Fs_Cocyls = idxs_Fs_Cocyls_Res
            idxs_Fs_Fails = idxs_Fs_Fails_Res

            i += 1

            if i > (len(idxs_Fs_Cocyls) - 1):
                break

        idxs_Fs_perCyl_Out.append(idxs_Fs_Cocyls)
        idxs_Fs_Processed.extend(idxs_Fs_Cocyls)

    return idxs_Fs_perCyl_Out


def find_closed_cylinders_per_face_sets(rgBrep, cyls_perF, idxFs_perCyl, tolCyl):
    rgB_Dup = rgBrep.DuplicateBrep()
    rgB_Dup.Faces.ShrinkFaces()

    cyls_Out = []

    for lst_idxF in idxFs_perCyl:
        if len(lst_idxF) == 1:
            iF = lst_idxF[0]
            if rgB_Dup.Faces[iF].IsClosed(0) or rgB_Dup.Faces[iF].IsClosed(1):
                cyls_Out.append(cyls_perF[iF])
            continue

        rgBs_DupFaces = [rgB_Dup.Faces[iF].DuplicateFace(duplicateMeshes=False) for iF in lst_idxF]

        rgBs_Joined = rg.Brep.JoinBreps(
            brepsToJoin=rgBs_DupFaces,
            tolerance=2.0*sc.doc.ModelAbsoluteTolerance)

        for _ in rgBs_DupFaces: _.Dispose()

        if len(rgBs_Joined) != 1:
            raise Exception("Number of breps: {}".format(len(rgBs_Joined)))

        rgB_Joined = rgBs_Joined[0]
        #sc.doc.Objects.AddBrep(rgB_Joined)

        cyl = cyls_perF[lst_idxF[0]]
        #sEval = "cyl.IsFinite"; print(sEval,'=',eval(sEval))
        cyl.Height1 = 1000.0
        #sEval = "cyl.IsFinite"; print(sEval,'=',eval(sEval))
        cyl.Height2 = -1000.0
        #sEval = "cyl.IsFinite"; print(sEval,'=',eval(sEval))
        srf_cyl = cyl.ToRevSurface()
        #sc.doc.Objects.AddSurface(srf_cyl)

        rgBs_CutUpSurface = rg.Brep.CutUpSurface(
            surface=srf_cyl,
            curves=[_ for _ in rgB_Joined.Edges if _.Valence == rg.EdgeAdjacency.Naked],
            useEdgeCurves=True,
            tolerance=sc.doc.ModelAbsoluteTolerance)

        rgB_Joined.Dispose()

        if len(rgBs_CutUpSurface) != 1:
            raise Exception("Number of breps: {}".format(len(rgBs_CutUpSurface)))

        rgB_CutUpSurface = rgBs_CutUpSurface[0]

        rgB_CutUpSurface.Faces.ShrinkFaces()

        #sc.doc.Objects.AddBrep(rgB_CutUpSurface)

        if rgB_CutUpSurface.Faces[0].IsClosed(0) or rgB_CutUpSurface.Faces[0].IsClosed(1):
            cyls_Out.append(cyl)

        rgB_CutUpSurface.Dispose()

    rgB_Dup.Dispose()

    return cyls_Out


def markPossibleMismatches(cyls_Closed_perB, tol_Min, tol_Max):
    tol_Deg = 5.0
    tol_Rad = Rhino.RhinoMath.ToRadians(tol_Deg)

    attrib_Red = rd.ObjectAttributes()
    attrib_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attrib_Red.ObjectColor = Color.Red

    gPts_Out = []
    gCrvs_Out = []

    for iThisBrep in range(len(cyls_Closed_perB)-1):
        for iOtherBrep in range(1, len(cyls_Closed_perB)):
            sc.escape_test()
            #print(iThisBrep, iOtherBrep)
            cyls_This = cyls_Closed_perB[iThisBrep]
            cyls_Other = cyls_Closed_perB[iOtherBrep]
            for iCyl_This in range(len(cyls_This)):
                cyl_This = cyls_This[iCyl_This]
                center_This = cyl_This.Center
                line_This = rg.Line(center_This, cyl_This.Axis)
                for iCyl_Other in range(len(cyls_Other)):
                    cyl_Other = cyls_Other[iCyl_Other]

                    # WIP:

                    isParallelTo = cyl_This.Axis.IsParallelTo(
                        cyl_Other.Axis,
                        angleTolerance=tol_Rad)
                    #sEval = "isParallelTo"; print(sEval,'=',eval(sEval))

                    if isParallelTo == 0:
                        continue

                    #angle_Rad = rg.Vector3d.VectorAngle(cyl_This.Axis, cyl_Other.Axis) # Values range from 0 to 180 degrees.
                    #angle_Deg = Rhino.RhinoMath.ToDegrees(angle_Rad)
                    #sEval = "angle_Deg"; print(sEval,'=',eval(sEval))
                    #if angle_Rad > tol_Rad:
                    #    continue

                    center_Other = cyl_Other.Center
                    line_Other = rg.Line(cyl_Other.Center, cyl_Other.Axis)
                    bSuccess, tA, tB = rg.Intersect.Intersection.LineLine(line_This, line_Other)
                    if not bSuccess:
                        continue
                    ptA = line_This.PointAt(tA)
                    ptB = line_Other.PointAt(tB)
                    dist = ptA.DistanceTo(ptB)
                    if tol_Min < dist <= tol_Max:
                        gPt = sc.doc.Objects.AddPoint(ptA, attrib_Red)
                        if gPt == Guid.Empty: raise Exception("Point could not be added.")
                        gPts_Out.append(gPt)
                        gPt = sc.doc.Objects.AddPoint(ptB, attrib_Red)
                        if gPt == Guid.Empty: raise Exception("Point could not be added.")
                        gPts_Out.append(gPt)
                        if ptA.DistanceTo(center_This) > 100.0 * tol_Min:
                            gCrv = sc.doc.Objects.AddLine(rg.Line(center_This, ptA), attrib_Red)
                        else:
                            gCrv = sc.doc.Objects.AddLine(line_This, attrib_Red)
                        if gCrv == Guid.Empty: raise Exception("Curve could not be added.")
                        gCrvs_Out.append(gCrv)
                        if ptB.DistanceTo(center_Other) > 100.0 * tol_Min:
                            gCrv = sc.doc.Objects.AddLine(rg.Line(center_Other, ptB), attrib_Red)
                        else:
                            gCrv = sc.doc.Objects.AddLine(line_Other, attrib_Red)
                        if gCrv == Guid.Empty: raise Exception("Curve could not be added.")
                        gCrvs_Out.append(gCrv)


                    #sEval = "cyl_This.Center"; print(sEval,'=',eval(sEval))
                    #sEval = "cyl_Other.Center"; print(sEval,'=',eval(sEval))
                    #dist = center_This.DistanceTo(center_Other)
                    #sEval = "dist"; print(sEval,'=',eval(sEval))

    return gPts_Out, gCrvs_Out


def main():

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select breps <All normal>",
        acceptNothing=True,
        filter=rd.ObjectType.Brep)
    if res != Rhino.Commands.Result.Success: return


    if objrefs:
        rdBreps = [objref.Object() for objref in objrefs]
    else:
        rdBreps = _getAllNormal_breps_and_instances()
        sEval = "len(rdBreps)"; print(sEval,'=',eval(sEval))

    if len(rdBreps) == 0:
        print("No breps.")
        return

    if len(rdBreps) == 1:
        print("Only 1 brep, so no further processing.")
        return

    for rdB in rdBreps:
        try:
            rdB.BrepGeometry
        except:
            ss = [
                "Block instances found",
                "_ExplodeBlock the instances and rerun this script.",
                "(Remember to _UndoMultiple the _ExplodeBlock later.)",
                ]
            raise Exception("\n".join(ss))

    tolCyl_mm = 1e-6
    sEval = "tolCyl_mm"; print(sEval,'=',eval(sEval))

    tolCyl = tolCyl_mm * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    sEval = "tolCyl"; print(sEval,'=',eval(sEval))

    tol_Min = tolCyl
    sEval = "tol_Min"; print(sEval,'=',eval(sEval))

    tol_Max = 0.25 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Inches,
        sc.doc.ModelUnitSystem)
    sEval = "tol_Max"; print(sEval,'=',eval(sEval))

    rdBs_withCyls = []
    dicts_face_cyl_perB = []
    Rhino.RhinoApp.SetCommandPromptMessage("Collecting cylinders in each brep...")
    for rdB in rdBreps:
        dict_face_cyl = collect_cylinders_per_face_index(
            rdB.BrepGeometry,
            tolCyl)
        if dict_face_cyl:
            rdBs_withCyls.append(rdB)
            dicts_face_cyl_perB.append(dict_face_cyl)

            #for idx_F in dict_face_cyl:
                #compIdx = rg.ComponentIndex(
                #        rg.ComponentIndexType.BrepFace,
                #        index=idx_F)
                #rdB.SelectSubObject(
                #        componentIndex=compIdx,
                #        select=True, syncHighlight=True, persistentSelect=True)


    idxFs_perCyl_perB = []
    Rhino.RhinoApp.SetCommandPromptMessage("Finding cocylindric faces in each brep...")
    for iB, rdB in enumerate(rdBs_withCyls):
        #Rhino.RhinoApp.SetCommandPromptMessage(
        #    "Processing cylinders from brep {} against brep {}, both out of {}".format(
        #        iThisBrep+1, iOtherBrep+1, len(cyls_Closed_perB)))

        dict_face_cyl = dicts_face_cyl_perB[iB]
        idxFs_perCyl = get_lists_of_face_indices_for_cocylindric_faces(
            rdB.BrepGeometry,
            dict_face_cyl,
            tolCyl)
        #sEval = "idxFs_perCyl"; print(sEval,'=',eval(sEval))
        idxFs_perCyl_perB.append(idxFs_perCyl)

    cyls_Closed_perB = []
    #sCPM = 
    #Rhino.RhinoApp.SetCommandPromptMessage("Finding closed, cocylindric faces in each brep...")

    for iB, rdB in enumerate(rdBs_withCyls):
        sc.escape_test()
        Rhino.RhinoApp.SetCommandPromptMessage(
            "Finding closed, cocylindric faces in brep {} of {}".format(
                iB+1, len(rdBs_withCyls)))

        dict_face_cyl = dicts_face_cyl_perB[iB]
        idxFs_perCyl = idxFs_perCyl_perB[iB]
        cyls_Closed = find_closed_cylinders_per_face_sets(
            rdB.BrepGeometry,
            dict_face_cyl,
            idxFs_perCyl,
            tolCyl)
        cyls_Closed_perB.append(cyls_Closed)

    #sEval = "len(rdBs_withCyls)"; print(sEval,'=',eval(sEval))
    #sEval = "len(cyls_Closed_perB)"; print(sEval,'=',eval(sEval))

    Rhino.RhinoApp.SetCommandPromptMessage("Marking possible mismatches...")
    gPts_Res, gCrvs_Res = markPossibleMismatches(cyls_Closed_perB, tol_Min, tol_Max)

    if not (gPts_Res or gCrvs_Res):
        print("No possible mismatches found.")
        return

    print("Possible mismatches found. {} points and {} curves added.".format(len(gPts_Res), len(gCrvs_Res)))

    sc.doc.Views.Redraw()

    #import pprint
    #pprint.pprint(rdBs_withCyls)
    #pprint.pprint(dicts_face_cyl_perB)


if __name__ == '__main__': main()