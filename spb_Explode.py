"""
This script is to be used as an alternative to _Explode and _ExplodeBlocks.

Breps
_Explode doesn't retain face colors.
This script retains both face and object colors.

Curves (wires)
_Explode splits NURBS curves at knots with Gsmooth_continuous ("Aesthetic discontinuity")
https://developer.rhino3d.com/api/rhinocommon/rhino.geometry.continuity
This script explodes NURBS curves at knots that are G2 discontinuous
per Rhino's _GCon and non-overload Curve.GetNextDiscontinuity.
G2 discontinuous is where at least one of the following two items is true:
    1. abs(k0-k1) / max(k0,k1) > 0.05, where k0 and k1 are magnitudes of the curvature
    vectors below and above the knot.
    2. The difference in curvature vectors is > 2.0 degrees.

This script explode polycurves to their actual segments.

Blocks
Same function as _Explode, except:
    Also report the resultant object types. (_Explode only reports the number of objects.)
    Display Color, Linetype, Print Color, Print Width, or Section Style of block objects that are
    set to ByParent are set to the instances setting after this script execution.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
230906-07: Created.
231009: Bug fix in working with PolylineCurves.
231024: Improved printed feedback.
231113: Now, polyline curve segments in polycurves are exploded immediately.
250113: Curves are now exploded similar to that by _Explode, except
        NurbsCurves are exploded to G2 (per Rhino (see above)) instead of Gsmooth_continuous.
250811: Updated notes above concerning block objects that are ByParent before this script execution.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import spb_Brep_explode
import spb_ObjectDescriptions


def get_G2_discontinuities(nc):
    """
    Parameters:
    Returns: list of float: Parameters of discontinuities
    """

    t0 = nc.Domain.Min
    t1 = nc.Domain.Max

    if nc.IsClosed:
        continuityType = rg.Continuity.G2_locus_continuous
    else:
        continuityType = rg.Continuity.G2_continuous

    ts = []

    while True:
        sc.escape_test()
        bSuccess, t = nc.GetNextDiscontinuity(continuityType, t0, t1)
        if not bSuccess:
            return ts

        ts.append(t)
        t0 = t


def explode_curve(rgC_In):
    """
    Parameters:
    Returns:
    """

    sc.escape_test()

    if isinstance(rgC_In, (rg.ArcCurve, rg.LineCurve)):
        return

    if isinstance(rgC_In, rg.PolylineCurve):
        return rgC_In.DuplicateSegments()

    if isinstance(rgC_In, rg.PolyCurve):
        #pc_In = rgC_In.Duplicate()
        #bNestingWasRemoved = pc_In.RemoveNesting()
        segs = rgC_In.DuplicateSegments()
        segs_Out = []
        for seg in segs:
            rv = explode_curve(seg)
            if rv is None:
                segs_Out.append(seg)
            else:
                segs_Out.extend(rv)
        return segs_Out

    if not isinstance(rgC_In, rg.NurbsCurve):
        raise Exception("Why is {} here?".format(rgC_In.GetType().Name))

    # Process NurbsCurve.

    ts_discont = get_G2_discontinuities(rgC_In)
    if not ts_discont:
        return
    if len(ts_discont) == 1 and rgC_In.IsClosed:
        return

    return rgC_In.Split(ts_discont)


def doesInstDefContainInvalidGeometry(rdInstDef):
    for rdO in rdInstDef.GetObjects():
        if not rdO.Geometry.IsValid:
            return True
    return False


def main():

    rdObjs_Pre = list(sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False))

    rdObjs_Post = []

    if not rdObjs_Pre:
        res, objrefs = ri.RhinoGet.GetMultipleObjects(
            "Select objects to explode",
            acceptNothing=False,
            filter=rd.ObjectType.Curve | rd.ObjectType.Brep | rd.ObjectType.InstanceReference)
        if res != Rhino.Commands.Result.Success: return

        for objref in objrefs:
            rdObjs_Post.append(objref.Object())


    iCt_Exploded = 0
    iCt_ExplodeResults = 0

    gOuts = []
    sLogs = []

    for rdObj in rdObjs_Pre + rdObjs_Post:

        if rdObj.ObjectType == rd.ObjectType.Brep:
            if not rdObj.Geometry.IsValid:
                sLogs.append("Brep with invalid geometry skipped.  Repair and try again.")
                continue
            gBs_Res, sLogs = spb_Brep_explode.explodeBrepObject(rdObj, bEcho=False, bDebug=False)
            if gBs_Res:
                gOuts.extend(gBs_Res)
                iCt_ExplodeResults += len(gBs_Res)
                iCt_Exploded += 1

        elif rdObj.ObjectType == rd.ObjectType.Curve:
            crvs_ret = explode_curve(rdObj.CurveGeometry)
            if not crvs_ret: continue
            
            bFail_to_add_object = False

            for crv in crvs_ret:
                gOut = sc.doc.Objects.AddCurve(crv, rdObj.Attributes)
                if gOut == gOut.Empty:
                    bFail_to_add_object = True
                else:
                    gOuts.append(gOut)
                    iCt_ExplodeResults += 1

            if not bFail_to_add_object:
                sc.doc.Objects.Delete(rdObj)

            iCt_Exploded += 1

        elif rdObj.ObjectType == rd.ObjectType.InstanceReference:
            rdInstRef = rdObj
            rdInstDef = rdInstRef.InstanceDefinition
            #print(rdInstDef.ObjectCount)
            if doesInstDefContainInvalidGeometry(rdInstDef):
                print("Block instance with invalid geometry skipped.  Repair first or use _Explode.")
                continue
            gObjs_from_Inst = sc.doc.Objects.AddExplodedInstancePieces(
                rdInstRef,
                explodeNestedInstances=False,
                deleteInstance=True)
            if gObjs_from_Inst is None:
                if iCt_Exploded:
                    print("Error exploding one block instance, but {} other objects were exploded into {} objects.".format(
                        iCt_Exploded, len(gOuts)))
                raise Exception(
                    "Block instance could not explode." \
                    "  It possibly has invalid b-reps." \
                    "  _Undo to retrieve possbly lost data and use _Explode instead.")
            for g in gObjs_from_Inst:
                if g == g.Empty:
                    raise Exception("Object in instance could not be added to document.  Undo and investigate.")
            gOuts.extend(gObjs_from_Inst)
            iCt_Exploded += 1

    if sLogs:
        for sLog in sorted(set(sLogs)):
            print("{} of {}".format(sLogs.count(sLog), sLog)) 

    if not iCt_Exploded:
        print("Nothing was exploded.")
        return

    print("{} objects were exploded into {} objects.".format(
        iCt_Exploded, len(gOuts)))

    sc.doc.Objects.UnselectAll()

    sc.doc.Objects.Select(objectIds=gOuts)

    spb_ObjectDescriptions.main()

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()