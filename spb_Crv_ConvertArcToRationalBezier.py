#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script converts
1. ArcCurves
2. NurbsCurves within ZeroTolerance (~2.328e-10) deviation tolerance of an arc
3. Optionally, NurbsCurves within an input deviation tolerance of an arc
to rational, single-spanned (Bezier) NurbsCurves.

Options:
    AllowNurbs: Yes or No whether to convert NurbsCurves at user-input tolerance (IsArcTol)
    IsArcTol: User-input tolerance for NurbsCurve that pass as arcs
    TargetDegree: 2, 3, or 5
    AngleLimitAction: RaiseDegree (keeps single span) or SplitSpans (keeps degree)
    MaxAngleDeg2: Maximum angle a single degree-2 span can represent (90-150)
    MaxAngleDeg3: Maximum angle a single degree-3 span can represent (180-210)
    DeleteInput: Yes to replace the existing curve or No to add a new curve

This script was created with the use of Google Gemini 3.1 Pro.

Send any questions, comments, or script development service needs to @spb on the McNeel Forums: https://discourse.mcneel.com/
"""

"""
260614: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.DocObjects as rd
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import math


def create_single_span_deg2(arc):
    """Generates an exact degree-2 single-span rational Bezier arc."""
    r = arc.Radius
    theta = arc.Angle
    w1 = math.cos(theta / 2.0)

    p0 = arc.StartPoint
    p2 = arc.EndPoint

    mid_pt = arc.MidPoint
    mid_dir = mid_pt - arc.Center
    mid_dir.Unitize()

    dist = r / w1
    p1 = arc.Center + mid_dir * dist

    nc = rg.NurbsCurve(3, True, 3, 3)
    nc.Points.SetPoint(0, p0, 1.0)
    nc.Points.SetPoint(1, p1, w1)
    nc.Points.SetPoint(2, p2, 1.0)

    nc.Knots[0] = 0.0
    nc.Knots[1] = 0.0
    nc.Knots[2] = 1.0
    nc.Knots[3] = 1.0
    return nc


def create_single_span_deg3(arc):
    """Generates an exact degree-3 single-span rational Bezier arc."""
    r = arc.Radius
    theta = arc.Angle

    # The exact interior weight for a cubic arc
    w1 = (1.0 + 2.0 * math.cos(theta / 2.0)) / 3.0

    p0 = arc.StartPoint
    p3 = arc.EndPoint

    # THE FIX: The numerator coefficient must be 2.0, not 3.0
    L = r * (2.0 * math.sin(theta / 2.0)) / (1.0 + 2.0 * math.cos(theta / 2.0))

    arc_crv = rg.ArcCurve(arc)
    t0_vec = arc_crv.TangentAt(arc_crv.Domain.Min)
    t3_vec = arc_crv.TangentAt(arc_crv.Domain.Max)

    # P2 goes backwards from the end, hence the subtraction
    p1 = p0 + t0_vec * L
    p2 = p3 - t3_vec * L

    nc = rg.NurbsCurve(3, True, 4, 4)
    nc.Points.SetPoint(0, p0, 1.0)
    nc.Points.SetPoint(1, p1, w1)
    nc.Points.SetPoint(2, p2, w1)
    nc.Points.SetPoint(3, p3, 1.0)

    nc.Knots[0] = 0.0
    nc.Knots[1] = 0.0
    nc.Knots[2] = 0.0
    nc.Knots[3] = 1.0
    nc.Knots[4] = 1.0
    nc.Knots[5] = 1.0
    return nc


def create_quintic_arc(arc):
    """Generates an exact degree-5 single-span rational Bezier arc."""
    r = arc.Radius
    center = arc.Center
    pt_start = arc.StartPoint

    x_axis = pt_start - center
    x_axis.Unitize()
    z_axis = arc.Plane.ZAxis
    y_axis = rg.Vector3d.CrossProduct(a=z_axis, b=x_axis)

    start_plane = rg.Plane(
        origin=center,
        xDirection=x_axis,
        yDirection=y_axis)

    pts_2d = [
        rg.Point3d(r, 0, 0),
        rg.Point3d(r, 4 * r, 0),
        rg.Point3d(-3 * r, 2 * r, 0),
        rg.Point3d(-3 * r, -2 * r, 0),
        rg.Point3d(r, -4 * r, 0),
        rg.Point3d(r, 0, 0)
    ]
    weights = [5.0, 1.0, 1.0, 1.0, 1.0, 5.0]

    full_circ = rg.NurbsCurve(
        dimension=3,
        rational=True,
        order=6,
        pointCount=6)

    for i in range(6):
        pt_3d = start_plane.PointAt(u=pts_2d[i].X, v=pts_2d[i].Y)
        full_circ.Points.SetPoint(index=i, point=pt_3d, weight=weights[i])

    full_circ.Knots.CreateUniformKnots(knotSpacing=1.0)

    if arc.IsCircle:
        return full_circ

    pt_end = arc.EndPoint
    rc, t1 = full_circ.ClosestPoint(testPoint=pt_end)

    trimmed_crv = full_circ.Trim(t0=0.0, t1=t1)
    if not trimmed_crv:
        return None

    return trimmed_crv.ToNurbsCurve()


def get_sub_arc(arc, i, num_spans):
    """Safely extracts a mathematically perfect sub-arc for multi-span curves."""
    arc_crv = rg.ArcCurve(arc)
    dom = arc_crv.Domain
    step = dom.Length / num_spans
    t0 = dom.Min + i * step
    t1 = dom.Min + (i + 1) * step

    trimmed = arc_crv.Trim(rg.Interval(t0, t1))
    bSuccess, sub_arc = trimmed.TryGetArc()
    return sub_arc


def create_multispan_rational_arc(arc, degree, num_spans):
    """Constructs a continuous NURBS curve with exact C0 knot junctions."""
    if num_spans == 1:
        if degree == 2: crv = create_single_span_deg2(arc)
        elif degree == 3: crv = create_single_span_deg3(arc)
        elif degree == 5: crv = create_quintic_arc(arc)
        
        # Force exact double-precision snap to original arc endpoints
        last_idx = crv.Points.Count - 1
        crv.Points.SetPoint(0, arc.StartPoint, crv.Points[0].Weight)
        crv.Points.SetPoint(last_idx, arc.EndPoint, crv.Points[last_idx].Weight)
        return crv

    cvs = []
    for i in range(num_spans):
        sub_arc = get_sub_arc(arc, i, num_spans)

        if degree == 2:
            crv = create_single_span_deg2(sub_arc)
        elif degree == 3:
            crv = create_single_span_deg3(sub_arc)
        elif degree == 5:
            crv = create_quintic_arc(sub_arc)

        crv_pts = crv.Points
        start_idx = 1 if i > 0 else 0
        for j in range(start_idx, crv_pts.Count):
            pt = crv_pts[j].Location
            w = crv_pts[j].Weight
            cvs.append((pt, w))

    cv_count = len(cvs)
    nc = rg.NurbsCurve(3, True, degree + 1, cv_count)

    for i, (pt, w) in enumerate(cvs):
        nc.Points.SetPoint(i, pt, w)

    knot_idx = 0
    for i in range(num_spans + 1):
        for _ in range(degree):
            nc.Knots[knot_idx] = float(i)
            knot_idx += 1

    # Force exact double-precision snap to original overall arc endpoints
    last_idx = nc.Points.Count - 1
    nc.Points.SetPoint(0, arc.StartPoint, nc.Points[0].Weight)
    nc.Points.SetPoint(last_idx, arc.EndPoint, nc.Points[last_idx].Weight)

    return nc


class RationalBezierGetter(ri.Custom.GetObject):
    def __init__(self):
        super(RationalBezierGetter, self).__init__()

        init_nurbs = sc.sticky.get("ConvertArcToRB_NurbsCrvs", True)
        self.opt_nurbs = ri.Custom.OptionToggle(init_nurbs, "No", "Yes")

        # Retrieve sticky values or use defaults
        init_tol = sc.sticky.get("ConvertArcToRB_IsArcTol", Rhino.RhinoMath.ZeroTolerance)
        self.tol_val = init_tol
        self.opt_tol = ri.Custom.OptionDouble(self.tol_val)

        # Using OptionInteger natively avoids Rhino list parsing bugs for numeric strings
        self.deg_val = sc.sticky.get("ConvertArcToRB_TargetDegree", 2)
        self.opt_deg = ri.Custom.OptionInteger(self.deg_val, 2, 5)

        # Using OptionToggle neatly wraps the internal boolean logic you requested
        self.raise_degree_val = sc.sticky.get("ConvertArcToRB_RaiseDegree", True)
        self.opt_action = ri.Custom.OptionToggle(self.raise_degree_val, "SplitSpans", "RaiseDegree")

        self.deg2_val = sc.sticky.get("ConvertArcToRB_MaxDeg2", 120.0)
        self.opt_deg2_limit = ri.Custom.OptionDouble(self.deg2_val, 90.0, 150.0)

        self.deg3_val = sc.sticky.get("ConvertArcToRB_MaxDeg3", 180.0)
        self.opt_deg3_limit = ri.Custom.OptionDouble(self.deg3_val, 180.0, 210.0)

        init_delete = sc.sticky.get("ConvertArcToRB_DeleteInput", True)
        self.opt_delete = ri.Custom.OptionToggle(init_delete, "No", "Yes")

        self.SetCommandPrompt("Select arcs to convert")
        self.GeometryFilter = rd.ObjectType.Curve

        self.AcceptNumber(True, True)
        self.BuildCustomOptions()

    def BuildCustomOptions(self):
        self.ClearCommandOptions()

        self.opt_idx_action = -1
        self.opt_idx_deg2 = -1
        self.opt_idx_deg3 = -1
        self.opt_idx_tol = -1

        self.AddOptionToggle(
            englishName="AllowNurbs", 
            toggleValue=self.opt_nurbs
        )

        if self.opt_nurbs.CurrentValue:
            self.opt_idx_tol = self.AddOptionDouble(
                englishName="IsArcTol", 
                numberValue=self.opt_tol
            )[0]

        self.opt_idx_degree = self.AddOptionInteger(
            englishName="TargetDegree", 
            intValue=self.opt_deg
        )[0]

        if self.opt_deg.CurrentValue in [2, 3]:
            self.opt_idx_action = self.AddOptionToggle(
                englishName="AngleLimitAction", 
                toggleValue=self.opt_action
            )[0]

        if self.opt_deg.CurrentValue == 2:
            self.opt_idx_deg2 = self.AddOptionDouble(
                englishName="MaxAngleDeg2", 
                numberValue=self.opt_deg2_limit
            )[0]

        if self.opt_deg.CurrentValue == 3 or (self.opt_deg.CurrentValue == 2 and self.opt_action.CurrentValue):
            self.opt_idx_deg3 = self.AddOptionDouble(
                englishName="MaxAngleDeg3", 
                numberValue=self.opt_deg3_limit
            )[0]

        self.AddOptionToggle(
            englishName="DeleteInput", 
            toggleValue=self.opt_delete
        )[0]

    def CustomGeometryFilter(self, rh_object, geometry, component_index):
        if not self.opt_nurbs.CurrentValue:
            return geometry.IsArc(tolerance=Rhino.RhinoMath.ZeroTolerance)
        return True


def save_sticky_options(go):
    """Saves the validated command options to the scriptcontext sticky dictionary."""
    sc.sticky["ConvertArcToRB_TargetDegree"] = go.deg_val
    sc.sticky["ConvertArcToRB_RaiseDegree"] = go.raise_degree_val
    sc.sticky["ConvertArcToRB_MaxDeg2"] = go.deg2_val
    sc.sticky["ConvertArcToRB_MaxDeg3"] = go.deg3_val
    sc.sticky["ConvertArcToRB_NurbsCrvs"] = go.opt_nurbs.CurrentValue
    sc.sticky["ConvertArcToRB_DeleteInput"] = go.opt_delete.CurrentValue
    sc.sticky["ConvertArcToRB_IsArcTol"] = go.tol_val


def main():
    go = RationalBezierGetter()

    while True:
        go.BuildCustomOptions()
        res = go.GetMultiple(1, 0)

        if res == Rhino.Input.GetResult.Option:
            op = go.Option()
            #sEval = "op.Index"; print(sEval, ':', eval(sEval))
            #sEval = "go.opt_idx_degree"; print(sEval, ':', eval(sEval))
            if op.Index == go.opt_idx_degree:
                val = go.opt_deg.CurrentValue
                # Snap 4 to 5 because degree-4 is unoptimized for exact arcs
                if val == 4:
                    val = go.opt_deg.CurrentValue = 5
                go.deg_val = val
                go.opt_deg = ri.Custom.OptionInteger(val, 2, 5)
            elif op.Index == go.opt_idx_action:
                go.raise_degree_val = go.opt_action.CurrentValue
            elif op.Index == go.opt_idx_deg2:
                go.deg2_val = go.opt_deg2_limit.CurrentValue
            elif op.Index == go.opt_idx_deg3:
                go.deg3_val = go.opt_deg3_limit.CurrentValue
            elif op.Index == go.opt_idx_tol:
                val = go.opt_tol.CurrentValue
                if val < 0.0:
                    val = Rhino.RhinoMath.ZeroTolerance
                    go.opt_tol = ri.Custom.OptionDouble(val)
                go.tol_val = val
            continue

        elif res == Rhino.Input.GetResult.Number:
            val = go.Number()
            if val < 0.0:
                val = Rhino.RhinoMath.ZeroTolerance
            go.tol_val = val
            go.opt_tol = ri.Custom.OptionDouble(val)
            continue

        elif res != Rhino.Input.GetResult.Object:
            save_sticky_options(go=go)
            return

        break

    save_sticky_options(go=go)

    tol = go.tol_val if go.opt_nurbs.CurrentValue else Rhino.RhinoMath.ZeroTolerance

    target_deg = go.deg_val
    raise_degree = go.raise_degree_val
    max_d2 = go.deg2_val
    max_d3 = go.deg3_val

    gs_Replaced = []
    gs_Added = []
    gs_AlreadyConverted = []

    for objref in go.Objects():
        crv = objref.Curve()
        if not crv:
            continue

        if isinstance(crv, rg.ArcCurve):
            arc = crv.Arc
        else:
            bSuccess, arc = crv.TryGetArc(tolerance=tol)
            if not bSuccess:
                continue

        # Mathematical logic cascade for degree limits
        arc_angle_deg = math.degrees(arc.Angle)
        final_deg = target_deg
        splits = 1

        if target_deg == 2:
            if arc_angle_deg <= max_d2:
                pass 
            elif raise_degree:
                if arc_angle_deg <= max_d3:
                    final_deg = 3
                else:
                    final_deg = 5
            else:
                splits = int(math.ceil(arc_angle_deg / max_d2))

        elif target_deg == 3:
            if arc_angle_deg <= max_d3:
                pass 
            elif raise_degree:
                final_deg = 5
            else:
                splits = int(math.ceil(arc_angle_deg / max_d3))

        # Check if geometry already exactly matches target specifications
        if isinstance(crv, rg.NurbsCurve):
            if crv.Degree == final_deg and crv.SpanCount == splits and crv.IsRational:
                gs_AlreadyConverted.append(objref.ObjectId)
                continue

        delete_input = False if isinstance(crv, rg.BrepEdge) else go.opt_delete.CurrentValue

        new_crv = create_multispan_rational_arc(
            arc=arc,
            degree=final_deg,
            num_spans=splits)
        if not new_crv:
            continue

        if delete_input:
            if sc.doc.Objects.Replace(objectId=objref.ObjectId, curve=new_crv):
                gs_Replaced.append(objref.ObjectId)
        else:
            attr = rd.ObjectAttributes()
            attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
            g_Added = sc.doc.Objects.AddCurve(curve=new_crv, attributes=attr)
            if g_Added != Guid.Empty:
                gs_Added.append(g_Added)

    if not go.ObjectsWerePreselected:
        for objref in go.Objects():
            sc.doc.Objects.Select(objref=objref, select=False)
        for objectId in gs_Replaced:
            sc.doc.Objects.Select(objectId=objectId, select=True)
        for objectId in gs_Added:
            sc.doc.Objects.Select(objectId=objectId)

    sc.doc.Views.Redraw()

    if gs_AlreadyConverted:
        print("{} curves already matched the target criteria.".format(len(gs_AlreadyConverted)))

    if gs_Replaced:
        print("Replaced {} curves with exact rational Bezier NURBS geometry.".format(len(gs_Replaced)))

    if gs_Added:
        print("Added {} exact rational Bezier NURBS curves.".format(len(gs_Added)))

    if not gs_Replaced and not gs_Added and not gs_AlreadyConverted:
        print("No curves were converted.")


if __name__ == '__main__': main()