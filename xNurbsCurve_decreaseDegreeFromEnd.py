"""
201124-30: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bConservePickedEnd'; keys.append(key)
    values[key] = False
    names[key] = 'ConserveContinuityEnd'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'OppPicked', 'Picked')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = False
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
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
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
    def setValues(cls):
        for key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
            if key == 'fTol':
                if cls.riOpts[key].CurrentValue < 0.0:
                    cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick curve near an end")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve


    def geomFilter_Curve(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if isinstance(geom, rg.BrepEdge):
            # DuplicateCurve gets the edge as a curve, which may be a subset of the EdgeCurve.
            rgC = geom.DuplicateCurve()
        elif isinstance(geom, rg.Curve):
            rgC = geom
        else:
            return False

        if isinstance(rgC, rg.NurbsCurve):
            nc = rgC
        elif isinstance(rgC, rg.PolyCurve):
            nc = rgC.ToNurbsCurve()
        else:
            return False

        if nc.IsPeriodic:
            print "Periodic curves are not supported."
            return False

        if nc.Degree == 1:
            print "Ignorded edge because NURBS curve is degree 1."
            return False

        return True


    go.SetCustomGeometryFilter(geomFilter_Curve)

    go.DisablePreSelect()

    idxs_Opt = {}

    while True:
        key = 'bConservePickedEnd'; idxs_Opt[key] = Opts.addOption(go, key)
        key = 'bReplace'; idxs_Opt[key] = Opts.addOption(go, key)
        key = 'bEcho'; idxs_Opt[key] = Opts.addOption(go, key)
        key = 'bDebug'; idxs_Opt[key] = Opts.addOption(go, key)

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            rgC = objref.Curve()
            iDegree_Start = rgC.Degree

            if iDegree_Start == 2:
                bDegree_Target = 1
            else:
                rc, bDegree_Target = ri.RhinoGet.GetInteger(
                    "Change degree from {} to".format(iDegree_Start),
                    acceptNothing=True,
                    outputNumber=iDegree_Start-1 if iDegree_Start < 5 else 3,
                    lowerLimit=1,
                    upperLimit=iDegree_Start-1)

                if rc == Rhino.Commands.Result.Cancel:
                    return

            return tuple([objref, bDegree_Target] + [Opts.values[key] for key in Opts.keys])

        # An option was selected.
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getScaleFactorForDegreeChange(nc_From, bConserveT1End, iDeg_To):
    """
    """

    scale = 1.0


    for i in range(nc_From.Degree-1, iDeg_To-1, -1):
        scale *= (float(i)**2 - 1.0) / (float(i)**2)


    if nc_From.SpanCount > 1:
        if bConserveT1End:
            iKnot_AdjToEnd = nc_From.Knots.Count - nc_From.Degree - 1
        else:
            iKnot_AdjToEnd = nc_From.Degree

        if nc_From.Knots.KnotMultiplicity(iKnot_AdjToEnd) == 1:
            # Adjust scale factor for knot vector.
            if bConserveT1End:
                domSpan = nc_From.SpanDomain(nc_From.SpanCount-1)
                domSpan_Adj = nc_From.SpanDomain(nc_From.SpanCount-2)
            else:
                domSpan = nc_From.SpanDomain(0)
                domSpan_Adj = nc_From.SpanDomain(1)

            m = (domSpan.Length + domSpan_Adj.Length) / domSpan.Length

            scale *= m


    if nc_From.IsRational:
        # Adjust scale factor for rationality.
        if bConserveT1End:
            w0 = nc_From.Points.GetWeight(nc_From.Points.Count-1)
            w1 = nc_From.Points.GetWeight(nc_From.Points.Count-2)
            w2 = nc_From.Points.GetWeight(nc_From.Points.Count-3)
        else:
            w0 = nc_From.Points.GetWeight(0)
            w1 = nc_From.Points.GetWeight(1)
            w2 = nc_From.Points.GetWeight(2)

        m = w1**2 / (w0 * w2)

        scale *= m


    return scale


def createCurve(nc_From, iDeg_To, bConserveT1End, bDebug=False):
    """
    Parameters:
        nc_From: rg.NurbsCurve
        iDeg_To: int of target degree
        bConserveT1End: True for T1, False for T0
        bDebug: bool
    Returns on success: rg.NurbsCurve
    Returns on fail: None
    """

    if not isinstance(nc_From, rg.NurbsCurve): return


    nc_In = nc_From


    # Periodic curves are not supported.
    if nc_In.IsPeriodic: return


    pts_New = []

    scale = getScaleFactorForDegreeChange(nc_In, bConserveT1End, iDeg_To) if iDeg_To > 1 else 1.0

    if not bConserveT1End:
        # Don't include points on T1 end.
        pts_New.append(nc_In.Points[0].Location)
        xform = rg.Transform.Scale(pts_New[0], scale)
        for i in range(1, iDeg_To+1):
            pt_New = rg.Point3d(nc_In.Points[i].Location)
            pt_New.Transform(xform)
            pts_New.append(pt_New)
    else:
        # Don't include points on T0 end.
        pt_Anchor = nc_In.Points[nc_In.Points.Count-1].Location
        xform = rg.Transform.Scale(pt_Anchor, scale)
        for i in range(nc_In.Points.Count-iDeg_To-1, nc_In.Points.Count-1):
            pt_New = rg.Point3d(nc_In.Points[i].Location)
            pt_New.Transform(rg.Transform.Scale(pt_Anchor, scale))
            pts_New.append(pt_New)
        pts_New.append(pt_Anchor)

    if bDebug:
        for i in range(len(pts_New)):
            pt = pts_New[i]
            rgDot = rg.TextDot("{}".format(i), pt)
            rgDot.FontHeight = 11
            sc.doc.Objects.AddTextDot(rgDot)
        sc.doc.Views.Redraw()

    nc_Out = rg.NurbsCurve.Create(
        periodic=False,
        degree=iDeg_To,
        points=pts_New)

    return nc_Out


def processCurveObject(objref_In, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iDegree = getOpt('iDegree')
    bConservePickedEnd = getOpt('bConservePickedEnd')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rgC_In = objref_In.Curve()

    if isinstance(rgC_In, rg.BrepEdge):
        nc_In = rgC_In.ToNurbsCurve()
    elif isinstance(rgC_In, rg.NurbsCurve):
        nc_In = rgC_In
    elif isinstance(rgC_In, rg.PolyCurve):
        nc_In = rgC_In.ToNurbsCurve()
    else:
        return


    bSuccess, t_AtPicked = nc_In.ClosestPoint(objref_In.SelectionPoint())
    if not bSuccess:
        return

    bPickedEndIsT1 = t_AtPicked > nc_In.Domain.Mid

    rgT_Sel = objref_In.Trim()

    bConserveT1End = bPickedEndIsT1 == bConservePickedEnd



    nc_Res = createCurve(
        nc_From=nc_In,
        bConserveT1End=bConserveT1End,
        iDeg_To=iDegree,
        bDebug=bDebug,
        )

    if nc_Res is None:
        print "Curve could not be created."
        return


    if not bReplace or objref_In.Edge():
        gC_Out = sc.doc.Objects.AddCurve(nc_Res)
        if gC_Out == gC_Out.Empty:
            print "Could not add curve."
        else:
            print "Curve was added."
    else:
        if sc.doc.Objects.Replace(objref_In.ObjectId, nc_Res):
            gC_Out = objref_In.ObjectId
            print "Replaced curve."
        else:
            print "Could not replace curve."

    return gC_Out


def main():
    
    rc = getInput()
    if rc is None: return

    (
        objref_In,
        iDegree,
        bConservePickedEnd,
        bReplace,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObject(
        objref_In=objref_In,
        iDegree=iDegree,
        bConservePickedEnd=bConservePickedEnd,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if gC_Res is None:
        return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()