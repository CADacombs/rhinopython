#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script will translate the G1 [1] (1st from the end) [0] from the end [0] control point
by a scale (1 is no change).
For G2 setting of MaintainPicked or MaintainOpp, the respective G2 [2] control points will also be translated.
"""

"""
210303, 0307: Created.
260420: WIP Adding preview and number slider to dynamically change preview
        before accepting value.
        Refactoring.
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


    key = 'bBothEnds'; keys.append(key)
    values[key] = True
    names[key] = 'Adjust'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'PickedEnd', 'BothEnds')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fScale'; keys.append(key)
    values[key] = 0.5
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCont_Picked'; keys.append(key)
    listValues[key] = 'G0', 'G1', 'G2', 'None' # All items must be strings.
    values[key] = 2
    names[key] = 'MaintainPicked'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCont_Opp'; keys.append(key)
    listValues[key] = 'G0', 'G1', 'G2', 'None' # All items must be strings.
    values[key] = 2
    names[key] = 'MaintainOpp'
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

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fScale':
            if cls.riOpts[key].CurrentValue <= 1e-9:
                print("Invalid input for scale value.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            if cls.riOpts[key].CurrentValue == 1:
                print("A scale of 1 will not modify the curve. Scale was not modified.")
                cls.riOpts[key].CurrentValue = cls.values[key]
                return
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            print("What happened?")
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve with picked end and optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick curve near an end")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve


    def geomFilter_Curve(rdObj, geom, compIdx):
        #print(rdObj, geom, compIdx.ComponentIndexType, compIdx.Index)

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
            print("Ignorded degree 1 NURBS curve.")
            return False

        #if nc.Degree + 1 < nc.Points.Count:
        #    print("Ignored non-Bezier curve."
        #    return


        return True


    go.SetCustomGeometryFilter(geomFilter_Curve)

    go.DisablePreSelect()

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bBothEnds')
        addOption('fScale')
        addOption('iCont_Picked')
        addOption('iCont_Opp')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            return (
                objref,
                Opts.values['bBothEnds'],
                Opts.values['fScale'],
                -1 if Opts.values['iCont_Picked'] == 3 else Opts.values['iCont_Picked'],
                -1 if Opts.values['iCont_Opp'] == 3 else Opts.values['iCont_Opp'],
                Opts.values['bReplace'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fScale'
            goNumber = go.Number()
            if goNumber == 1:
                print("A scale of 1 will not modify the curve. Scale was not modified.")
                continue
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createCurve(nc_In, fScale, iEndToScale=2, iG_T0=2, iG_T1=2, bDebug=False):
    """
    Parameters:
        nc_In: rg.NurbsCurve
        fScale: float of scale factor
        iTEndsToScale: int(0 for T0, 1 for T1, or 3 for both ends nc_In)
        iG_T0: int(-1 for No continuity, 0 for G0, 1 for G1, or 2 for G2)
        iG_T1: int(-1 for No continuity, 0 for G0, 1 for G1, or 2 for G2)
        bDebug: bool
    Returns on success: rg.NurbsCurve
    Returns on fail: None
    """

    if iG_T0 is None and iG_T1 is None: return

    if nc_In.IsPeriodic: return
    if not isinstance(nc_In, rg.NurbsCurve): return
    #if nc_In.Degree + 1 < nc_In.Points.Count: return # because not Bezier.

    if nc_In.Points.Count < (iG_T0 + 1) + (iG_T1 + 1):
        print("Curve needs {} more points to maintain continuities.".format(
            (iG_T0 + 1) + (iG_T1 + 1) - nc_In.Points.Count))
        return

    pts_Prime = []


    if iEndToScale in (0,2):
        # Scale T0 end.


        if iG_T0 in (0,1,2):
            p0 = nc_In.Points[0].Location
            pts_Prime.append(p0)

            if iG_T0 in (1,2):
                p1 = nc_In.Points[1].Location

                xform = rg.Transform.Scale(p0, fScale)

                p1p = rg.Point3d(p1)
                p1p.Transform(xform)
                pts_Prime.append(p1p)

                if iG_T0 == 2:
                    p2 = nc_In.Points[2].Location

                    m = ((p1p - p0).Length/(p1 - p0).Length)**2.0

                    p2p = 2.0*p1p + -p0 + m*(-2.0*p1 + p2 + p0)

                    pts_Prime.append(p2p)



    # Duplicate the points that are not needed by either end scale.
    iCt_PtsNotNeededByT1Scale = nc_In.Points.Count - (iG_T1 + 1)
    for i in range(len(pts_Prime), iCt_PtsNotNeededByT1Scale):
        pts_Prime.append(nc_In.Points[i].Location)


    if iEndToScale in (1,2):
        # Scale T1 end.


        if iG_T1 in (0,1,2):
            pts_New_T1 = [] # To be reversed and extended on pts_Prime.

            p0 = nc_In.Points[nc_In.Points.Count-1].Location
            pts_New_T1.append(p0)

            if iG_T1 in (1,2):
                p1 = nc_In.Points[nc_In.Points.Count-2].Location

                xform = rg.Transform.Scale(p0, fScale)

                p1p = rg.Point3d(p1) # G1 CP location prime.
                p1p.Transform(xform)
                pts_New_T1.append(p1p)

                if iG_T1 == 2:
                    p2 = nc_In.Points[nc_In.Points.Count-3].Location

                    m = ((p1p - p0).Length/(p1 - p0).Length)**2.0

                    p2p = 2.0*p1p + -p0 + m*(-2.0*p1 + p2 + p0)

                    pts_New_T1.append(p2p)

            pts_New_T1.reverse()

            pts_Prime.extend(pts_New_T1)


    #for pt in pts_Prime:
    #    sc.doc.Objects.AddPoint(pt)
    #sc.doc.Views.Redraw(); return


    nc_Out = nc_In.Duplicate()

    for i in range(nc_Out.Points.Count):
        nc_Out.Points.SetPoint(
            index=i,
            point=pts_Prime[i],
            weight=nc_In.Points.GetWeight(i))


    if bDebug:
        for i in range(len(pts_Prime)):
            pt = pts_Prime[i]
            rgDot = rg.TextDot("{}".format(i), pt)
            rgDot.FontHeight = 11
            sc.doc.Objects.AddTextDot(rgDot)
        sc.doc.Views.Redraw()



    return nc_Out


def processCurveObject(objref_In, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bBothEnds = getOpt('bBothEnds')
    fScale = getOpt('fScale')
    iCont_Picked = getOpt('iCont_Picked')
    iCont_Opp = getOpt('iCont_Opp')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if iCont_Picked == iCont_Opp == 0:
        print("Both continuities are set to 0, so curve cannot be modified. Script canceled.")
        return

    if fScale == 1:
        print("Scale set to 1, so curve cannot be modified. Script canceled.")
        return

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


    if t_AtPicked > nc_In.Domain.Mid:
        iEndToScale = 2 if bBothEnds else 1
        iG_T0, iG_T1 = iCont_Opp, iCont_Picked
    else:
        iEndToScale = 2 if bBothEnds else 0
        iG_T0, iG_T1 = iCont_Picked, iCont_Opp


    nc_Res = createCurve(
        nc_In=nc_In,
        fScale=fScale,
        iEndToScale=iEndToScale,
        iG_T0=iG_T0,
        iG_T1=iG_T1,
        bDebug=bDebug,
        )

    if nc_Res is None:
        print("Curve could not be created.")
        return


    if not bReplace or objref_In.Edge():
        gC_Out = sc.doc.Objects.AddCurve(nc_Res)
        if gC_Out == gC_Out.Empty:
            print("Could not add curve.")
        else:
            print("Curve was added.")
    else:
        if sc.doc.Objects.Replace(objref_In.ObjectId, nc_Res):
            gC_Out = objref_In.ObjectId
            print("Replaced curve.")
        else:
            print("Could not replace curve.")

    return gC_Out


def main():
    
    rv = getInput()
    if rv is None: return
    #print(rv; return

    (
        objref_In,
        bBothEnds,
        fScale,
        iCont_Picked,
        iCont_Opp,
        bReplace,
        bEcho,
        bDebug,
        ) = rv

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gC_Res = processCurveObject(
        objref_In=objref_In,
        bBothEnds=bBothEnds,
        fScale=fScale,
        iCont_Picked=iCont_Picked,
        iCont_Opp=iCont_Opp,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if gC_Res is None:
        return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()