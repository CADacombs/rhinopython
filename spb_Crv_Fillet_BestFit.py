"""
This script will find the best-fit fillet given a reference curve to replace and
the 2 curves to fillet.
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
230902-06: Created.
231004: Output points are selected to immediately be used for trimming the 2 curves to fillet.
231024: Bug fix in determining default resolution.  Improved printed output.
250307: Bug fix in a reported failure.

TODO:
    Add option whether to trim the curves to the fillet.
    Add option whether to add end points of fillet.
    Handle bad input such as a linear curve for curve to replace.
    Replace brute force finder with more robust one.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fResolution'; keys.append(key)
    values[key] = 10**(-int(abs(math.log10(sc.doc.ModelAbsoluteTolerance)))) # Tol. of 0.01 -> Res. of 0.1; 0.0004 -> 0.001; etc.
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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
    def setValue(cls, key, idxList=None):

        if key == 'fResolution':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 0.0

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_Crv_to_replace():
    """
    Get curve with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve to approximate with fillet")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    #go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    def customGeometryFilter(rdObj, rgObj, compIdx):
        return not rg.Curve.IsLinear(rgObj)

    go.SetCustomGeometryFilter(customGeometryFilter)

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('fResolution')
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
            key = 'fResolution'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Crvs_to_fillet():
    """
    Get curve with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves to fillet")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    #go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    #def customGeometryFilter(rdObj, rgObj, compIdx):
    #    return rdObj.Id != gCrv_1st

    #go.SetCustomGeometryFilter(customGeometryFilter)


    go.OneByOnePostSelect = True
    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _formatDistance(fDistance, iPrecision=15):
    if fDistance is None: return "(No deviation provided)"
    
    if fDistance < 1e-4:
        return "{:.{}e}".format(fDistance, iPrecision)
    
    #if fDistance < 0.1:
    #    return "{:.{}g}".format(fDistance, iPrecision)
    
    return "{:.{}g}".format(fDistance, iPrecision)


def _createArcCurveOfFillet(curve0, curve1, radius, t0Base, t1Base):
    arc = rg.Curve.CreateFillet(curve0, curve1, radius, t0Base, t1Base)
    if not arc.IsValid:
        raise Exception("CreateFillet failed at fillet R{}.  _Connect curves to each other and try again.".format(radius))
    #sEval = "arc.Length"; print("{}: {}".format(sEval, eval(sEval)))
    #sEval = "arc.IsValid"; print("{}: {}".format(sEval, eval(sEval)))
    return rg.ArcCurve(arc)


def findFillet(crv_Ref, crv_A_ToFillet, crv_B_ToFillet, tA, tB, fResolution, bEcho=True, bDebug=False):
    """
    """


    def _getDistancesBetweenCurves(curveA, curveB):
        rc = rg.Curve.GetDistancesBetweenCurves(curveA, curveB, tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
        if not rc[0]:
            sc.doc.Objects.AddCurve(curveB)
            raise Exception("GetDistancesBetweenCurves failed at fillet R{}.  _Connect curves to each other and try again.".format(fRadius_Start))
        return rc[1]

    fRadius_Raw_Start = 1.0/rg.Curve.CurvatureAt(crv_Ref, crv_Ref.Domain.Mid).Length

    if bDebug:
        print("Raw radius start: {}".format(fRadius_Raw_Start))
    m = round(1.0/fResolution, 0)
    fRadius_Start = round(m * fRadius_Raw_Start, 0) / m

    cA = crv_A_ToFillet
    cB = crv_B_ToFillet

    arccrv_Start = _createArcCurveOfFillet(cA, cB, fRadius_Start, tA, tB)
    dev_Start = _getDistancesBetweenCurves(crv_Ref, arccrv_Start)

    if bDebug:
        sEval = "fRadius_Start"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dev_Start"; print("{}: {}".format(sEval, eval(sEval)))

    # Determine whether the target radius is >, <, or equal to the start radius.
    fRadius_Plus = fRadius_Start + fResolution
    
    arccrv_Plus = _createArcCurveOfFillet(cA, cB, fRadius_Plus, tA, tB)
    dev_Plus = _getDistancesBetweenCurves(crv_Ref, arccrv_Plus)

    fRadius_Minus = fRadius_Start - fResolution
    
    arccrv_Minus = _createArcCurveOfFillet(cA, cB, fRadius_Minus, tA, tB)
    dev_Minus = _getDistancesBetweenCurves(crv_Ref, arccrv_Minus)


    if dev_Start < dev_Plus and dev_Start < dev_Minus:
        if bEcho: print("Start radius is the best radius.")
        return fRadius_Start, arccrv_Start, dev_Start

    if dev_Plus < dev_Start and dev_Plus < dev_Minus:
        if bEcho: print("Best radius is greater than the start radius.")
        fResolution_Signed = +fResolution
        fRadius_ToTry = fRadius_Plus
        arccrv_ToTry = arccrv_Plus
        dev_ToTry = dev_Plus
    elif dev_Minus < dev_Start and dev_Minus < dev_Plus:
        if bEcho: print("Best radius is less than the start radius.")
        fResolution_Signed = -fResolution
        fRadius_ToTry = fRadius_Minus
        arccrv_ToTry = arccrv_Minus
        dev_ToTry = dev_Minus
    else:
        sEval = "dev_Minus"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dev_Start"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "dev_Plus"; print("{}: {}".format(sEval, eval(sEval)))
        raise Exception("What happened?")

    # This is just to accomodate entering the while loop.
    arccrv_Best = arccrv_Start
    dev_Best = dev_Start

    iCt_While = 0

    while dev_ToTry < dev_Best:
        iCt_While += 1
        arccrv_Best.Dispose()
        arccrv_Best = arccrv_ToTry
        fRadius_Best = fRadius_ToTry
        dev_Best = dev_ToTry
        fRadius_ToTry += fResolution_Signed
        arccrv_ToTry = _createArcCurveOfFillet(cA, cB, fRadius_ToTry, tA, tB)
        dev_ToTry = _getDistancesBetweenCurves(crv_Ref, arccrv_ToTry)
        if bDebug:
            sEval = "fRadius_ToTry"; print("{}: {}".format(sEval, eval(sEval)))
            sEval = "dev_ToTry"; print("{}: {}".format(sEval, eval(sEval)))

    print("After {} iterations of while loop.".format(iCt_While))

    return fRadius_Best, arccrv_Best, dev_Best


def main():

    objref_Ref = getInput_Crv_to_replace()
    if objref_Ref is None: return

    objrefs_ToFillet = getInput_Crvs_to_fillet()
    if objrefs_ToFillet is None: return

    fResolution = Opts.values['fResolution']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    crv_Ref = objref_Ref.Curve()

#    crv_A_ToFillet = objrefs_ToFillet[0].Curve()
#    crv_B_ToFillet = objrefs_ToFillet[1].Curve()

    crv_A_ToFillet, tA = rd.ObjRef.CurveParameter(objrefs_ToFillet[0])
    crv_B_ToFillet, tB = rd.ObjRef.CurveParameter(objrefs_ToFillet[1])

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    rc = findFillet(
        crv_Ref=crv_Ref,
        crv_A_ToFillet=crv_A_ToFillet,
        crv_B_ToFillet=crv_B_ToFillet,
        tA=tA,
        tB=tB,
        fResolution=fResolution,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if rc is not None:
        fRadius_Winner, fArcCrve_Winner, fDev_Winner = rc

    print("Best fillet at {} resolution is R{}".format(fResolution, fRadius_Winner))

    print("Deviation from reference curve: {}".format(
        _formatDistance(fDev_Winner, iPrecision=max((6, sc.doc.ModelDistanceDisplayPrecision)))))

    sc.doc.Objects.AddCurve(fArcCrve_Winner)
    sc.doc.Objects.Select(sc.doc.Objects.AddPoint(fArcCrve_Winner.PointAtStart))
    sc.doc.Objects.Select(sc.doc.Objects.AddPoint(fArcCrve_Winner.PointAtEnd))

    print("{} points added and are selected.".format(
        len(list(sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)))))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()