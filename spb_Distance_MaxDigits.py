#! python 2
from __future__ import absolute_import, division, print_function, unicode_literals

"""
This script inputs 2 sets of 2 points each for 2 distances then reports some
relationship between the distances.
"""

"""
260526-27: Created.
"""

import Rhino.Input as ri
import scriptcontext as sc


def getPoints_noDecimalPlaceOption():

    gp = ri.Custom.GetPoint()
    gp.SetCommandPrompt("First point for distance")

    gp.AcceptNothing(False)
    gp.AcceptNumber(True, acceptZero=True)

    res = gp.Get()

    if res == ri.GetResult.Cancel:
        gp.Dispose()
        return

    if res == ri.GetResult.Point:
        pt1 = gp.Point()

    gp.Dispose()

    gp = ri.Custom.GetPoint()
    gp.SetCommandPrompt("Second point for distance")

    gp.AcceptNothing(False)
    gp.AcceptNumber(True, acceptZero=True)

    gp.DrawLineFromPoint(startPoint=pt1, showDistanceInStatusBar=True)
    gp.EnableDrawLineFromPoint(True)

    res = gp.Get()

    if res == ri.GetResult.Cancel:
        gp.Dispose()
        return

    if res == ri.GetResult.Point:
        pt2 = gp.Point()
        gp.Dispose()
        return pt1, pt2

    if res == ri.GetResult.Number:
        riOpt.CurrentValue = gp.Number()


def getPoints_withDecimalPlaceOption():

    gp = ri.Custom.GetPoint()
    gp.SetCommandPrompt("First point for distance")

    gp.AcceptNothing(False)
    gp.AcceptNumber(True, acceptZero=True)

    key = 'iDecPlaces'
    stickyKey = '{}({})'.format(key, __file__)
    iDecPlaces_Start = sc.sticky[stickyKey] if stickyKey in sc.sticky else 17

    riOpt = ri.Custom.OptionInteger(initialValue=iDecPlaces_Start)

    while True:
        idxOpt = gp.AddOptionInteger('DecPlaces', riOpt)[0]
        res = gp.Get()

        if res == ri.GetResult.Cancel:
            gp.Dispose()
            return

        if res == ri.GetResult.Point:
            pt1 = gp.Point()
            break

        if res == ri.GetResult.Number:
            riOpt.CurrentValue = gp.Number()
        elif idxOpt != 1:
            print("What happened?")

        sc.sticky[stickyKey] = riOpt.CurrentValue

        gp.ClearCommandOptions()


    gp.Dispose()

    gp = ri.Custom.GetPoint()
    gp.SetCommandPrompt("Second point for distance")

    gp.AcceptNothing(False)
    gp.AcceptNumber(True, acceptZero=True)

    gp.DrawLineFromPoint(startPoint=pt1, showDistanceInStatusBar=True)
    gp.EnableDrawLineFromPoint(True)

    while True:
        idxOpt = gp.AddOptionInteger('DecPlaces', riOpt)[0]
        res = gp.Get()

        if res == ri.GetResult.Cancel:
            gp.Dispose()
            return

        if res == ri.GetResult.Point:
            pt2 = gp.Point()
            gp.Dispose()
            return (pt1, pt2), riOpt.CurrentValue

        if res == ri.GetResult.Number:
            riOpt.CurrentValue = gp.Number()
        elif idxOpt != 1:
            print("What happened?")

        sc.sticky[stickyKey] = riOpt.CurrentValue

        gp.ClearCommandOptions()


def sFloat_OLD(fFloat, iDecPlaces):
    """float to str with decimal places."""
    if 0.0 < abs(fFloat) < 0.001: #10.0**(-(iDecPlaces)):
        return "{:.{}e}".format(fFloat, iDecPlaces)
    else:
        return "{:.{}g}".format(fFloat, iDecPlaces)


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
    
    #    rv = getPoints_withDecimalPlaceOption()
    #    if rv is None: return
    #    (pt1, pt2), iDecPlaces = rv

    rv = getPoints_noDecimalPlaceOption()
    if rv is None: return
    pt1, pt2 = rv

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