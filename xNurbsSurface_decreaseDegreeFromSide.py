"""
201119-30: Created.
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


    key = 'bConservePickedSide'; keys.append(key)
    values[key] = False
    names[key] = 'ConserveContinuitySide'
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

    go.SetCommandPrompt("Select edge of untrimmed NURBS surface")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.EdgeFilter


    def geomFilter_TrimOfFull4TrimNS(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx.ComponentIndexType, compIdx.Index

        if not isinstance(geom, rg.BrepTrim):
            return False

        rgT = geom
        rgB = rgT.Brep
        #if rgB.Faces.Count > 1:
        #    return False

        rgF = rgT.Face
        if not rgF.IsSurface:
            return False

        if not rgF.Loops.Count == 1:
            return False

        rgL = rgF.Loops[0]

        rgS = rgF.UnderlyingSurface()

        if isinstance(rgS, rg.NurbsSurface):
            ns = rgS
        elif isinstance(rgS, rg.SumSurface):
            ns = rgS.ToNurbsSurface()
        else:
            return False

        if ns.IsPeriodic(0) or ns.IsPeriodic(1):
            print "Periodic surfaces are not supported."
            return False

        if ns.IsRational:
            print "Rational NURBS surfaces are not supported."
            return False

        if not rgL.Trims.Count == 4:
            return False

        # Reject any NS with singular trims.
        for iT in range(4):
            if rgB.Trims[iT].TrimType == rg.BrepTrimType.Singular:
                print "NurbsSurface with singular trim is not supported."
                return False

        if rgT.IsoStatus == rg.IsoStatus.East:
            iDegree_Start = ns.Degree(0)
            iSpanCt = ns.SpanCount(0)
        elif rgT.IsoStatus == rg.IsoStatus.West:
            iDegree_Start = ns.Degree(0)
            iSpanCt = ns.SpanCount(0)
        else:
            iDegree_Start = ns.Degree(1)
            iSpanCt = ns.SpanCount(1)

        if iDegree_Start == 1:
            print "Ignorded edge because its surface is degree 1 from edge."
            return False

        #if iSpanCt > 1:
        #    print "Ignored edge because surface has {} spans from edge.".format(iSpanCt)
        #    return False


        return True


    go.SetCustomGeometryFilter(geomFilter_TrimOfFull4TrimNS)

    idxs_Opt = {}

    while True:
        key = 'bConservePickedSide'; idxs_Opt[key] = Opts.addOption(go, key)
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
            rdBrep_In = objref.Object()
            rgSrf = objref.Surface()
            rgTrim = objref.Trim()
            if rgTrim.IsoStatus == rg.IsoStatus.East:
                iDegree_Start = rgSrf.Degree(0)
            elif rgTrim.IsoStatus == rg.IsoStatus.West:
                iDegree_Start = rgSrf.Degree(0)
            else:
                iDegree_Start = rgSrf.Degree(1)
            
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


def getScaleFactorForDegreeChange(ns_From, side_ToConserve, iDeg_To):
    """
    """

    iDir = int(side_ToConserve in (rg.IsoStatus.South, rg.IsoStatus.North))

    degree = ns_From.Degree(iDir)

    scale = 1.0


    for i in range(degree-1, iDeg_To-1, -1):
        scale *= (float(i)**2 - 1.0) / (float(i)**2)


    if ns_From.SpanCount(iDir) > 1:
        knots = ns_From.KnotsV if iDir else ns_From.KnotsU
        if side_ToConserve in (rg.IsoStatus.East, rg.IsoStatus.North):
            iKnot_AdjToEnd = knots.Count - degree - 1
        else:
            iKnot_AdjToEnd = degree

        multy_AdjKnot = knots.KnotMultiplicity(iKnot_AdjToEnd)

        if multy_AdjKnot == 1:
            # Adjust scale factor for knot vector.
            iCt_Span = ns_From.SpanCount(iDir)
            t_Knots_Span = ns_From.GetSpanVector(iDir)

            if side_ToConserve in (rg.IsoStatus.East, rg.IsoStatus.North):
                fLen_spanEnd = t_Knots_Span[iCt_Span] - t_Knots_Span[iCt_Span-1]
                fLen_spanAdj = t_Knots_Span[iCt_Span-1] - t_Knots_Span[iCt_Span-2]
            else:
                fLen_spanEnd = t_Knots_Span[1] - t_Knots_Span[0]
                fLen_spanAdj = t_Knots_Span[2] - t_Knots_Span[1]

            m = (fLen_spanEnd + fLen_spanAdj) / fLen_spanEnd

            scale *= m


    return scale


def createSurface(ns_From, iDeg_To, side_ToConserve, bDebug=False):
    """
    Parameters:
        ns_From: rg.NurbsSurface
        iDeg_To: int of target degree
        side_ToConserve: rg.IsoStatus
        bDebug: bool
    Returns on success: rg.NurbsSurface
    Returns on fail: None
    """

    if not isinstance(ns_From, rg.NurbsSurface): return

    ns_In = ns_From

    if ns_In.IsRational: return

    if ns_In.IsPeriodic(0) or ns_In.IsPeriodic(1): return


    pts_New = [] # Order for NurbsSurface.CreateFromPoints is u0,v0, u0,v1, ...

    scale = getScaleFactorForDegreeChange(ns_In, side_ToConserve, iDeg_To) if iDeg_To > 1 else 1.0

    if side_ToConserve in (rg.IsoStatus.South, rg.IsoStatus.North):

        if side_ToConserve == rg.IsoStatus.South:
            # Don't include northernmost rows.
            for iU in range(ns_In.Points.CountU):
                pt_Anchor = ns_In.Points.GetControlPoint(iU, 0).Location
                pts_New.append(pt_Anchor)
                xform = rg.Transform.Scale(pt_Anchor, scale)
                for iV in range(1, iDeg_To+1):
                    pt_New = rg.Point3d(ns_In.Points.GetControlPoint(iU, iV).Location)
                    pt_New.Transform(rg.Transform.Scale(pt_Anchor, scale))
                    pts_New.append(pt_New)
        else:
            # Don't include southernmost rows.
            for iU in range(ns_In.Points.CountU):
                pt_Anchor = ns_In.Points.GetControlPoint(iU, ns_In.Points.CountV-1).Location
                xform = rg.Transform.Scale(pt_Anchor, scale)
                for iV in range(ns_In.Points.CountV-iDeg_To-1, ns_In.Points.CountV-1):
                    pt_New = rg.Point3d(ns_In.Points.GetControlPoint(iU, iV).Location)
                    pt_New.Transform(rg.Transform.Scale(pt_Anchor, scale))
                    pts_New.append(pt_New)
                pts_New.append(pt_Anchor)

        uCount_New = ns_In.Points.CountU
        vCount_New = iDeg_To + 1
        uDegree_New = ns_In.Degree(0)
        vDegree_New = iDeg_To

    elif side_ToConserve in (rg.IsoStatus.West, rg.IsoStatus.East):

        if side_ToConserve == rg.IsoStatus.West:
            # Don't include easternmost columns.
            for iV in range(ns_In.Points.CountV):
                pts_New.append(ns_In.Points.GetControlPoint(0, iV).Location)

            for iU in range(1, iDeg_To+1):
                for iV in range(ns_In.Points.CountV):
                    pt_New = rg.Point3d(ns_In.Points.GetControlPoint(iU, iV).Location)
                    pt_New.Transform(rg.Transform.Scale(pts_New[iV], scale))
                    pts_New.append(pt_New)
        else:
            # Don't include westernmost columns.
            for iU in range(ns_In.Points.CountU-iDeg_To-1, ns_In.Points.CountU-1):
                for iV in range(ns_In.Points.CountV):
                    pt_Anchor = ns_In.Points.GetControlPoint(ns_In.Points.CountU-1, iV).Location
                    pt_New = rg.Point3d(ns_In.Points.GetControlPoint(iU, iV).Location)
                    pt_New.Transform(rg.Transform.Scale(pt_Anchor, scale))
                    pts_New.append(pt_New)

            for iV in range(ns_In.Points.CountV):
                pts_New.append(ns_In.Points.GetControlPoint(ns_In.Points.CountU-1, iV).Location)

        uCount_New = iDeg_To + 1
        vCount_New = ns_In.Points.CountV
        uDegree_New = iDeg_To
        vDegree_New = ns_In.Degree(1)

    else:
        ns_Out = None

    if bDebug:
        for i in range(len(pts_New)):
            pt = pts_New[i]
            rgDot = rg.TextDot("{}".format(i), pt)
            rgDot.FontHeight = 11
            sc.doc.Objects.AddTextDot(rgDot)
        sc.doc.Views.Redraw()

    ns_Out = rg.NurbsSurface.CreateFromPoints(
        points=pts_New,
        uCount=uCount_New,
        vCount=vCount_New,
        uDegree=uDegree_New,
        vDegree=vDegree_New)

    ns_In.Dispose()

    return ns_Out


def processBrepObject(objref_In, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    iDegree = getOpt('iDegree')
    bConservePickedSide = getOpt('bConservePickedSide')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rdBrep_In = objref_In.Object()
    #rgBrep_In = rdBrep_In.Geometry.Duplicate() # So it is not DocumentControlled.
    rgF_In = objref_In.Surface()
    rgS_In = rgF_In.UnderlyingSurface()
    if isinstance(rgS_In, rg.NurbsSurface):
        ns_In = rgS_In.Duplicate()
    elif isinstance(rgS_In, rg.SumSurface):
        ns_In = rgS_In.ToNurbsSurface()
    else:
        return


    rgT_Sel = objref_In.Trim()

    if bConservePickedSide:
        side_Conserve = rgT_Sel.IsoStatus
    else:
        side_Conserve = rg.IsoStatus.ToObject(
            rg.IsoStatus, ((rgT_Sel.IsoStatus.value__ - 1) % 4) + 3)


    ns_Res = createSurface(
        ns_From=ns_In,
        iDeg_To=iDegree,
        side_ToConserve=side_Conserve,
        bDebug=bDebug,
        )

    if ns_Res is None:
        print "Surface could not be created."
        return

    if not bReplace:
        gB_Out = sc.doc.Objects.AddSurface(ns_Res)
        if gB_Out == gB_Out.Empty:
            print "Could not add modified surface."
        else:
            print "Surface was added."
    else:
        rdB_In = objref_In.Object()
        rgB_In = objref_In.Brep()
        rgB_WIP = rgB_In.DuplicateBrep()
        rgB_WIP.Faces.RemoveAt(objref_In.Face().FaceIndex)
        rgBs_Out = rg.Brep.CreateBooleanUnion([rgB_WIP], sc.doc.ModelAbsoluteTolerance)
        rgB_WIP.Dispose()

        if rgB_In.Faces.Count == 1:
            rgB_Out = ns_Res.ToBrep()
            if sc.doc.Objects.Replace(objref_In.ObjectId, rgB_Out):
                gB_Out = objref_In.ObjectId
                print "Replaced monoface brep with new surface."
            else:
                print "Could not replace monoface brep with new surface."
        else:
            attr = rdB_In.Attributes
            gB_Out = sc.doc.Objects.AddSurface(ns_Res, attr)
            if gB_Out == gB_Out.Empty:
                print "Could not add modified surface."
            else:
                if len(rgBs_Out) == 1:
                    if sc.doc.Objects.Replace(objref_In.ObjectId, rgBs_Out[0]):
                        print "Added new surface and deleted face of brep."
                    else:
                        print "Added new surface but could not delete face of brep."
                else:
                    gBs_Out = []
                    for rgB in rgBs_Out:
                        gBs_Out.append(sc.doc.Objects.AddBrep(rgB, attr))

                    if gBs_Out[0].Empty in gBs_Out:
                        for gB_Out in gBs_Out:
                            if gB_Out != gB_Out.Empty:
                                sc.doc.Objects.Delete(objectId=gB_Out, quiet=False)
                        print "Added new surface but could not delete face of brep."
                    else:
                        sc.doc.Objects.Delete(rdB_In)
                        print "Added new surface and deleted face of brep." \
                            "  Remainder of brep is now {} breps.".format(len(gBs_Out))

    return gB_Out


def main():
    
    rc = getInput()
    if rc is None: return

    (
        objref_In,
        iDegree,
        bConservePickedSide,
        bReplace,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gB_Res = processBrepObject(
        objref_In=objref_In,
        iDegree=iDegree,
        bConservePickedSide=bConservePickedSide,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if gB_Res is None:
        return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()