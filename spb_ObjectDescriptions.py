"""
This script prints a summary of selected objects, such as types, whether the objects
are open, etc., to the command history. This can be easier to review than some of the
output of _What or _List.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
150417: Created.
...
230906-07: Recreated.  NurbsCurves are now listed with their matching shapes.
230915: Now when the only non-curve output is one InstanceReference, its definition name is also reported.
250113: Added 'closed' and 'open' to curve report. Refactored.

TODO:
    Add closed, single surface count for breps.
    Add single surface and polysuface counts for breps.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc



def createCrvReport(myCrvCts, myNC_shape_cts):

    if not any(value > 0 for value in myCrvCts.values()): return

    iCt_Crvs = sum((value for value in myCrvCts.values()))

    s = "{} curve{} (".format(iCt_Crvs, 's' if iCt_Crvs > 1 else '')

    s_List_more = []

    if myCrvCts['lc']:
        s_List_more.append("{} line".format(myCrvCts['lc']))

    iCt_plcs = myCrvCts['plc_Closed'] + myCrvCts['plc_Open']
    if iCt_plcs:
        if iCt_plcs == 1:
            s_List_more.append("{} polyline ({})".format(
                iCt_plcs,
                "c" if myCrvCts['plc_Closed'] else "o"))
        else:
            ss = []
            if myCrvCts['plc_Closed']:
                if myCrvCts['plc_Closed'] == iCt_plcs:
                    ss.append("c")
                else:
                    ss.append("{}c".format(myCrvCts['plc_Closed']))
            if myCrvCts['plc_Open']:
                if myCrvCts['plc_Open'] == iCt_plcs:
                    ss.append("o")
                else:
                    ss.append("{}o".format(myCrvCts['plc_Open']))
            s_List_more.append("{} polyline ({})".format(
                iCt_plcs,
                ','.join(ss)))

    iCt_acs = myCrvCts['ac_Closed'] + myCrvCts['ac_Open']
    if iCt_acs:
        if iCt_acs == 1:
            s_List_more.append("{} arc ({})".format(
                iCt_acs,
                "c" if myCrvCts['ac_Closed'] else "o"))
        else:
            ss = []
            if myCrvCts['ac_Closed']:
                if myCrvCts['ac_Closed'] == iCt_acs:
                    ss.append("c")
                else:
                    ss.append("{}c".format(myCrvCts['ac_Closed']))
            if myCrvCts['ac_Open']:
                if myCrvCts['ac_Open'] == iCt_acs:
                    ss.append("o")
                else:
                    ss.append("{}o".format(myCrvCts['ac_Open']))
            s_List_more.append("{} arc ({})".format(
                iCt_acs,
                ','.join(ss)))

    iCt_pcs = myCrvCts['pc_Closed'] + myCrvCts['pc_Open']
    if iCt_pcs:
        if iCt_pcs == 1:
            s_List_more.append("{} polycrv ({})".format(
                iCt_pcs,
                "c" if myCrvCts['pc_Closed'] else "o"))
        else:
            ss = []
            if myCrvCts['pc_Closed']:
                if myCrvCts['pc_Closed'] == iCt_pcs:
                    ss.append("c")
                else:
                    ss.append("{}c".format(myCrvCts['pc_Closed']))
            if myCrvCts['pc_Open']:
                if myCrvCts['pc_Open'] == iCt_pcs:
                    ss.append("o")
                else:
                    ss.append("{}o".format(myCrvCts['pc_Open']))
            s_List_more.append("{} polycrv ({})".format(
                iCt_pcs,
                ','.join(ss)))

    iCt_ncs = myCrvCts['nc_Closed'] + myCrvCts['nc_Open']
    if iCt_ncs:
        if iCt_ncs == 1:
            s_List_more.append("{} NURBS ({})".format(
                iCt_ncs,
                "C" if myCrvCts['nc_Closed'] else "O"))
        else:
            ss = []
            if myCrvCts['nc_Closed']:
                if myCrvCts['nc_Closed'] == iCt_ncs:
                    ss.append("c")
                else:
                    ss.append("{}c".format(myCrvCts['nc_Closed']))
            if myCrvCts['nc_Open']:
                if myCrvCts['nc_Open'] == iCt_ncs:
                    ss.append("o")
                else:
                    ss.append("{}o".format(myCrvCts['nc_Open']))

            s_List_more.append("{} NURBS ({})".format(
                iCt_ncs,
                ','.join(ss)))

        if any(value > 0 for value in myNC_shape_cts.values()):
            ss = []
            if myNC_shape_cts['nc_Linear']:
                ss.append("{} linear".format(myNC_shape_cts['nc_Linear']))
            if myNC_shape_cts['nc_Circle']:
                ss.append("{} circular".format(myNC_shape_cts['nc_Circle']))
            if myNC_shape_cts['nc_Arc']:
                ss.append("{} arc-shaped".format(myNC_shape_cts['nc_Arc']))
            if myNC_shape_cts['nc_Ellipse']:
                ss.append("{} elliptical".format(myNC_shape_cts['nc_Ellipse']))
            if myNC_shape_cts['nc_EllipticalArc']:
                ss.append("{} elliptical arc-shaped".format(myNC_shape_cts['nc_EllipticalArc']))

            s_List_more.append("({})".format(', '.join(ss)))


    s += ", ".join(s_List_more)

    s += ")"

    return s


def main():

    myCrvCts = {
        'lc': 0,
        'plc_Closed': 0,
        'plc_Open': 0,
        'ac_Closed': 0,
        'ac_Open': 0,
        'pc_Closed': 0,
        'pc_Open': 0,
        'nc_Closed': 0,
        'nc_Open': 0,
        }


    myNC_shape_cts = {
        'nc_Linear': 0,
        'nc_Circle': 0,
        'nc_Arc': 0,
        'nc_Ellipse': 0,
        'nc_EllipticalArc': 0,
        }


    def createReport_ExceptCrvs():

        sDescrsWithCts = []

        otype = rd.ObjectType.Brep
        if otype in otypes:
            sDescrsWithCts.append(
                "{} {} ({} closed, {} open)".format(
                    otypes.count(otype),
                    otype,
                    iCt_ClosedBs,
                    otypes.count(otype) - iCt_ClosedBs))

        for otype in sorted(set(otypes)):
            if otype in (rd.ObjectType.Curve, rd.ObjectType.Brep):
                continue
            sDescrsWithCts.append("{} {}".format(otypes.count(otype), otype))

        if len(otypes) == 1 and otype == rd.ObjectType.InstanceReference:
            return "{}[{}]".format(sDescrsWithCts[0], rdO.InstanceDefinition.Name)

        return ", ".join(sDescrsWithCts)



    rdObjs_Pre = list(sc.doc.Objects.GetSelectedObjects(includeLights=True, includeGrips=False))

    rdObjs_Post = []

    if not rdObjs_Pre:
        res, objrefs = ri.RhinoGet.GetMultipleObjects(
            "Select objects",
            acceptNothing=False,
            filter=rd.ObjectType.AnyObject)
        if res != Rhino.Commands.Result.Success: return

        for objref in objrefs:
            rdObjs_Post.append(objref.Object())


    iCt_ClosedBs = 0


    otypes = []

    for rdO in rdObjs_Pre + rdObjs_Post:
        objType = rdO.ObjectType
        otypes.append(objType)

        #print(objType)

        if objType == rd.ObjectType.Curve:
            crv = rdO.Geometry

            if isinstance(crv, rg.ArcCurve):
                if crv.IsClosed:
                    myCrvCts['ac_Closed'] += 1
                else:
                    myCrvCts['ac_Open'] += 1
                continue

            if isinstance(crv, rg.LineCurve):
                myCrvCts['lc'] += 1
                continue

            if isinstance(crv, rg.NurbsCurve):
                if crv.IsClosed:
                    myCrvCts['nc_Closed'] += 1
                else:
                    myCrvCts['nc_Open'] += 1

                if crv.IsLinear():
                    myNC_shape_cts['nc_Linear'] += 1
                elif crv.IsArc():
                    if crv.IsCircle():
                        myNC_shape_cts['nc_Circle'] += 1
                    else:
                        myNC_shape_cts['nc_Arc'] += 1
                elif crv.IsEllipse():
                    if crv.IsClosed:
                        myNC_shape_cts['nc_Ellipse'] += 1
                    else:
                        myNC_shape_cts['nc_EllipticalArc'] += 1

                continue

            if isinstance(crv, rg.PolyCurve):
                if crv.IsClosed:
                    myCrvCts['pc_Closed'] += 1
                else:
                    myCrvCts['pc_Open'] += 1
                continue

            if isinstance(crv, rg.PolylineCurve):
                if crv.IsClosed:
                    myCrvCts['plc_Closed'] += 1
                else:
                    myCrvCts['plc_Open'] += 1
                continue

            raise Exception("Why are we here?")

        if objType == rd.ObjectType.Brep:
            brep = rdO.Geometry
            if brep.IsSolid: iCt_ClosedBs += 1


    sReport_Crvs = createCrvReport(myCrvCts, myNC_shape_cts)
    sReport_NoCrvs = createReport_ExceptCrvs()

    if sReport_Crvs: print(sReport_Crvs)
    if sReport_NoCrvs: print(sReport_NoCrvs)


if __name__ == '__main__': main()