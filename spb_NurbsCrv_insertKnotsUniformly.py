"""
Insert knots uniformly within existing spans.
Therefore, when the starting surface is uniform, uniformity will be maintained.
"""
"""
211105: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Drawing


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iQtyToAddInEachSpan'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=1)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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

        if not idxOpt: print "Add option for {} failed.".format(key)

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fScale':
            if cls.riOpts[key].CurrentValue <= 1e-9:
                print "Invalid input for scale value."
                cls.riOpts[key].CurrentValue = cls.values[key]

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    go.AcceptNumber(True, acceptZero=False)

    res = go.Get()

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    objref = go.Object(0)
    go.Dispose()

    crv = objref.Curve()

    if not isinstance(crv, rg.NurbsCurve):
        s = "Underlying surface is a {}".format(crv.GetType().Name)
        s += " and will be converted to a NurbsCurve."
        print s

        nc = crv.ToNurbsSurface()
    else:
        nc = crv

    print "Degree:{}, SpanCt:{}".format(
        nc.Degree, nc.SpanCount)

    go = ri.Custom.GetInteger()

    go.AcceptNothing(True)

    go.SetCommandPrompt("Quantity to add between each span")

    idxs_Opt = {}

    while True:
        go.SetCommandPromptDefault(str(Opts.values['iQtyToAddInEachSpan']))

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            break

        if res == ri.GetResult.Number:
            key = 'iQtyToAddInEachSpan'
            if int(go.Number()) < 1:
                continue
            Opts.riOpts[key].CurrentValue = int(go.Number())
            Opts.setValue(key)
            break

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break # out of for look back to while loop.


    return (
        objref,
        Opts.values['iQtyToAddInEachSpan'],
        Opts.values['bReplace'],
        Opts.values['bEcho'],
        Opts.values['bDebug'],
        )


def createCurve(nc_In, iQtyToAddInEachSpan, bDebug=False):
    """
    Parameters:
        nc_In: rg.NurbsSurface,
        iQtyToAddInEachSpan: int
        bDebug: bool

    Returns: rg.NurbsSurface
    """

    if not isinstance(nc_In, rg.NurbsCurve):
        return

    if iQtyToAddInEachSpan < 1:
        return

    nc_Out = nc_In.Duplicate()

    # Add knots from Domain end to beginning.

    for iK in range(nc_In.Knots.Count-1, 0, -1):

        k_R = nc_In.Knots[iK]
        k_L = nc_In.Knots[iK-1]

        for i in range(iQtyToAddInEachSpan):
            fraction_from_R = float(i+1) / float(iQtyToAddInEachSpan+1)
            k_M = fraction_from_R*k_R + (1.0-fraction_from_R)*k_L
            nc_Out.Knots.InsertKnot(k_M)

    return nc_Out


def processCurveObject(objref, iQtyToAddInEachSpan=1, bReplace=True, bEcho=True, bDebug=False):
    """
    objref: rd.ObjRef of wire rg.Curve.
    """

    rgC_In = rs.coercecurve(objref, raise_if_missing=False)
    if rgC_In is None:
        print "{} skipped.".format(rs.coercerhinoobject(objref).ObjectType)
        return

    nc_Res = createCurve(
        rgC_In.ToNurbsCurve(),
        iQtyToAddInEachSpan,
        bDebug)

    if nc_Res is None:
        return


    if bReplace:
        gC_In = rs.coerceguid(objref)

        if sc.doc.Objects.Replace(
            objectId=gC_In,
            curve=nc_Res
        ):
            if bEcho:
                if not isinstance(rgC_In, rg.NurbsCurve):
                    print "Replaced {} with a NurbsCurve with {} knots.".format(
                        rgC_In.GetType().Name, nc_Res.Knots.Count)
                else:
                    print "Changed knot count of NurbsCurve from {} to {}.".format(
                        rgC_In.Knots.Count, nc_Res.Knots.Count)
            return gC_In


    # Add new curve.
    gC_Out = sc.doc.Objects.AddCurve(nc_Res)
    if gC_Out != gC_Out.Empty:
        if bEcho:
            print "Added a NurbsCurve with {} knots.".format(
                nc_Res.Knots.Count)
        return gC_Out
    else:
        print "Could not add curve."


def main():

    rc = getInput()
    if rc is None: return

    (
        objref,
        iQtyToAddInEachSpan,
        bReplace,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBs_Res = processCurveObject(
        objref,
        iQtyToAddInEachSpan,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()