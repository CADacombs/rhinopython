"""
201209: Created.
210129: Bug fix in option input.
231223: Modified option input for ease of use.
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


    key = 'bSameDegree'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDegree'; keys.append(key)
    values[key] = 3
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=1)
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
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):
        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("Why is key, {}, here?  Value was not set or sticky-saved.".format(key))
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick curve near its end for start of extension")

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

        if rgC.IsPeriodic:
            print "Periodic curves are not supported."
            return False

        return True


    go.SetCustomGeometryFilter(geomFilter_Curve)

    go.DisablePreSelect()

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bSameDegree')
        if not Opts.values['bSameDegree']:
            addOption('iDegree')
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
            key = 'iDegree'
            if go.Number() == Opts.riOpts[key].CurrentValue:
                continue
            Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 1:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                go.ClearCommandOptions()
                break


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


def createCurve(nc_Ref, bExtFromT1, iDeg_To=0, bDebug=False):
    """
    Parameters:
        nc_Ref: rg.NurbsCurve
        bExtFromT1: True for T1, False for T0
        iDeg_To: int(target degree)  0 is to use same degree as nc_Ref.
        bDebug: bool
    Returns on success: rg.NurbsCurve
    Returns on fail: None
    """

    if not isinstance(nc_Ref, rg.NurbsCurve): return



    # Periodic curves are not supported.
    if nc_Ref.IsPeriodic: return



    # Create curve at existing degree.
    if bExtFromT1:
        t0 = nc_Ref.Domain.T1
        t1 = t0 + nc_Ref.Domain.Length
    else:
        t1 = nc_Ref.Domain.T0
        t0 = t1 - nc_Ref.Domain.Length

    intrvl = rg.Interval(t0, t1)

    nc_Out = nc_Ref.Extend(intrvl)

    nc_Out = nc_Out.Trim(intrvl)

    if iDeg_To == 0 or iDeg_To == nc_Ref.Degree:
        return nc_Out
    elif iDeg_To > nc_Ref.Degree:
        nc_Out.IncreaseDegree(iDeg_To)
        return nc_Out


    # Decrease degree.

    bConserveT1End = not bExtFromT1

    nc_WIP = nc_Out


    pts_New = []

    scale = getScaleFactorForDegreeChange(nc_WIP, bConserveT1End, iDeg_To) if iDeg_To > 1 else 1.0

    if not bConserveT1End:
        # Don't include points on T1 end.
        pts_New.append(nc_WIP.Points[0].Location)
        xform = rg.Transform.Scale(pts_New[0], scale)
        for i in range(1, iDeg_To+1):
            pt_New = rg.Point3d(nc_WIP.Points[i].Location)
            pt_New.Transform(xform)
            pts_New.append(pt_New)
    else:
        # Don't include points on T0 end.
        pt_Anchor = nc_WIP.Points[nc_WIP.Points.Count-1].Location
        xform = rg.Transform.Scale(pt_Anchor, scale)
        for i in range(nc_WIP.Points.Count-iDeg_To-1, nc_WIP.Points.Count-1):
            pt_New = rg.Point3d(nc_WIP.Points[i].Location)
            pt_New.Transform(rg.Transform.Scale(pt_Anchor, scale))
            pts_New.append(pt_New)
        pts_New.append(pt_Anchor)

    nc_WIP.Dispose()

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


def processCurveObject(objref_In, iDeg_To=0, bEcho=True, bDebug=False):
    """
    """

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

    bExtFromT1 = t_AtPicked > nc_In.Domain.Mid

    rgT_Sel = objref_In.Trim()


    nc_Res = createCurve(
        nc_Ref=nc_In,
        bExtFromT1=bExtFromT1,
        iDeg_To=iDeg_To,
        bDebug=bDebug,
        )

    if nc_Res is None:
        print "Curve could not be created."
        return


    gC_Out = sc.doc.Objects.AddCurve(nc_Res)
    if gC_Out == gC_Out.Empty:
        print "Could not add curve."
    else:
        print "Curve was added."

    return gC_Out


def main():

    objref_In = getInput()
    if objref_In is None: return

    bSameDegree = Opts.values['bSameDegree']
    iDegree = Opts.values['iDegree']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObject(
        objref_In=objref_In,
        iDeg_To=0 if bSameDegree else iDegree,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if gC_Res is None:
        return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()