#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
"""

"""
260526-27: Created.
260529: Now, the PointObjects and Geometry.Points are tested with IsValidWithLog.
"""

import Rhino.Geometry as rg
import Rhino.DocObjects as rd
import Rhino.Input as ri
import scriptcontext as sc


def getPointObjects_noDecimalPlaceOption():

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select 2 point objects")
    go.GeometryFilter = rd.ObjectType.Point
    go.AcceptNothing(False)
    go.AcceptNumber(False, acceptZero=True)

    res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Object:
        rdPts = go.Object(0).Object(), go.Object(1).Object()
        go.Dispose()
        return rdPts

    go.Dispose()


def getPointObjects_withDecimalPlaceOption():

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select 2 point objects")
    go.GeometryFilter = rd.ObjectType.Point
    go.AcceptNothing(False)
    go.AcceptNumber(True, acceptZero=True)

    key = 'iDecPlaces'
    stickyKey = '{}({})'.format(key, __file__)
    iDecPlaces_Start = sc.sticky[stickyKey] if stickyKey in sc.sticky else 17

    riOpt = ri.Custom.OptionInteger(initialValue=iDecPlaces_Start)

    while True:
        idxOpt = go.AddOptionInteger('DecPlaces', riOpt)[0]
        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            rdPts = go.Object(0).Object(), go.Object(1).Object()
            go.Dispose()
            return rdPts, riOpt.CurrentValue

        if res == ri.GetResult.Number:
            riOpt.CurrentValue = go.Number()
        elif idxOpt != 1:
            print("What happened?")

        sc.sticky[stickyKey] = riOpt.CurrentValue

        go.ClearCommandOptions()


def sFloat_OLD(fDistance, iDecPlaces):
    """float to str with decimal places."""
    if fDistance is None:
        return "(No deviation provided)"
    elif 0.0 < abs(fDistance) < 0.001: #10.0**(-(iDecPlaces)):
        return "{:.{}e}".format(fDistance, iDecPlaces)
    else:
        return "{:.{}g}".format(fDistance, iDecPlaces)


def sFloat(fFloat):
    if fFloat is None:
        return "(No deviation provided)"

    if fFloat == 0.0:
        if str(fFloat) == "-0.0":
            return "-0.0"
        return "0.0"

    if abs(fFloat) >= 0.001:
        return "{:.17g}".format(fFloat)

    # Exponential formatting.
    raw_e = "{:.16e}".format(fFloat)

    mantissa, exp = raw_e.split('e')

    if '.' in mantissa:
        mantissa = mantissa.rstrip('0').rstrip('.')

    exp_sign = exp[0]
    exp_num = exp[1:].lstrip('0')
    if not exp_num:
        exp_num = '0'
    exp_clean = exp_sign + exp_num

    return "{}e{}".format(mantissa, exp_clean)


def main():

    rv = getPointObjects_noDecimalPlaceOption()
    if rv is None: return
    rdPt1, rdPt2 = rv

    bIsValid, sLog = rdPt1.IsValidWithLog()
    if not bIsValid:
        print("PointObject 1 is not valid:")
        print(sLog)
    bIsValid, sLog = rdPt1.PointGeometry.IsValidWithLog()
    if not bIsValid:
        print("Geometry.Point 1 is not valid:")
        print(sLog)

    bIsValid, sLog = rdPt2.IsValidWithLog()
    if not bIsValid:
        print("PointObject 2 is not valid:")
        print(sLog)
    bIsValid, sLog = rdPt2.PointGeometry.IsValidWithLog()
    if not bIsValid:
        print("Geometry.Point 2 is not valid:")
        print(sLog)

    #    rv = getPointObjects_withDecimalPlaceOption()
    #    if rv is None: return
    #    (rdPt1, rdPt2), iDecPlaces = rv

    pt1 = rdPt1.PointGeometry.Location
    pt2 = rdPt2.PointGeometry.Location

    print("dx = {}  dy = {}  dz = {}".format(
        sFloat(pt2.X - pt1.X),
        sFloat(pt2.Y - pt1.Y),
        sFloat(pt2.Z - pt1.Z),
        ))
    print("Distance = {} {}".format(
        sFloat(pt1.DistanceTo(pt2)),
        sc.doc.GetUnitSystemName(
            modelUnits=True,
            capitalize=False,
            singular=False,
            abbreviate=False)))


if __name__ == '__main__': main()