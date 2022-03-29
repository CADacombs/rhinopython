"""
"""

"""
190909-13, 18: Created.
200113-16: Import-related updates.  Bug fix.
200119: Removed an option.
200120: Import-related updates resulting in less code in this script.
211029: Import-related update.  Added fAngleTol_Deg.  Disabled bC1 option since G1 is currently only possible.
220328: Import-related update.

TODO:
After scripts are completed for them, remove bMaximizeMinRadius.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

import xCurve_inflections
import spb_Crv_match
import xCurve_radiusMinima


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'bC1'; keys.append(key)
    values[key] = False
    names[key] = 'Continuity'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'G1', 'C1')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fAngleTol_Deg'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSkipIfAlreadyAtContinuity'; keys.append(key)
    values[key] = True
    names[key] = 'IfAlreadyAtContinuity'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Process', 'Skip')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitCrvDev'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMaximizeMinRadius'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get PolyCurve and option values.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select polycurve")


    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter.WireCurve
    
    # Custom geometry filter to only select NurbsCurve wire curves.
    def geomFilter_PolyCurve(rdObj, geom, compIdx):
        if rdObj.ObjectType == rd.ObjectType.InstanceReference: return False
        return isinstance(geom, rg.PolyCurve)
    go.SetCustomGeometryFilter(geomFilter_PolyCurve)    


    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    while True:
#        key = 'bC1'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'fAngleTol_Deg'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bSkipIfAlreadyAtContinuity'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bLimitCrvDev'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'fDevTol'; idxs_Opts[key] = (Opts.riAddOpts[key](go)
                                           if Opts.values['bLimitCrvDev'] else None)
        key = 'bMaximizeMinRadius'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bReplace'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])

        # An option was selected or a number was entered.

        key = 'fAngleTol_Deg'
        if Opts.riOpts[key].CurrentValue <= 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        key = 'fDevTol'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getMaximumDeviation(rgCrvA, rgCrvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            rgCrvA,
            rgCrvB,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
    if rc[0]:
        return rc[1]


def formatDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def processPolyCurve(rgPolyCrv0, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bC1 = getOpt('bC1')
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bSkipIfAlreadyAtContinuity = getOpt('bSkipIfAlreadyAtContinuity')
    fDevTol = getOpt('fDevTol')
    bMaximizeMinRadius = getOpt('bMaximizeMinRadius')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    c0 = rgPolyCrv0


    if not isinstance(c0, rg.PolyCurve):
        if bEcho: print "Curve is not a PolyCurve."
        c0.Dispose()
        return


    def isCurveOkToMod(seg):
        sType = seg.GetType().Name
        if sType == 'LineCurve':
            return False
        elif sType == 'ArcCurve':
            return False
        elif seg.IsLinear(1e-9):
            return False
        elif sType == 'NurbsCurve' and seg.IsLinear(1e-9):
            return False

        nc = seg.ToNurbsCurve()

        # Check whether the NurbsCurve contains any internal fully-multiple knots.
        # Don't bother with code for a periodic curve
        # since only PolyCurve segments are checked.
        if sType == 'NurbsCurve' and seg.SpanCount > 1:
            knots = seg.Knots
            iK = seg.Degree
            while True:
                sc.escape_test()
                m = knots.KnotMultiplicity(iK)
                if m == seg.Degree:
                    print "NurbsCurve segment containing some full multiplicity knots" \
                        " will not be processed."
                    nc.Dispose()
                    return False
                iK += m
                if iK == knots.Count - seg.Degree:
                    break

        nc.Dispose()

        return True



    segs0 = c0.DuplicateSegments()

    segs1 = []

    segA = segs0[0]
    
    bSomeSegsReplaced = False

    iSegs_Skipped = []

    for iSegB in xrange(1, len(segs0) + (1 if c0.IsClosed else 0)):
        
        if bDebug: sEval='iSegB'; print '#'*20+sEval+': ',eval(sEval)
        iSegB_toWatch = 18
        if iSegB == iSegB_toWatch:
            pass
        
        Rhino.RhinoApp.CommandPrompt = "Processing joint {} of {}".format(
                iSegB,
                len(segs0) - 0 if c0.IsClosed else 1
        )

        if iSegB == len(segs0):
            # Only occurs with closed curves.
            # Start with the modified segment.
            # TODO: Deviation may exceed limit.  Fix this.
            segB = segs1[0]
        else:
            segB = segs0[iSegB]

        bOkToModifyA = isCurveOkToMod(segA)
        bOkToModifyB = isCurveOkToMod(segB)

        if bDebug:
            sEval='bOkToModifyA'; print sEval+': ',eval(sEval)
            sEval='bOkToModifyB'; print sEval+': ',eval(sEval)

        if not bOkToModifyA and not bOkToModifyB:
            if bDebug:
                s  = "Neither of segments {} or {}".format(
                    segA.GetType().Name,
                    segB.GetType().Name)
                s += " can be modified within tolerance,"
                s + " and will be skipped."
                print s
            segs1.append(segA)
            segA = segB
            iSegs_Skipped.extend([iSegB-1, iSegB])
            continue

        rc = spb_Curve_match.createNurbsCurves(
                rgCurveA=segA,
                rgCurveB=segB,
                bT1WorkEnd_A=True,
                bT1WorkEnd_B=False,
                bModifyA=bOkToModifyA,
                bModifyB=bOkToModifyB,
                sContinuity='C1' if bC1 else 'G1',
                fAngleTol_Deg=fAngleTol_Deg,
                bSkipIfAlreadyAtContinuity=bSkipIfAlreadyAtContinuity,
                iPreserveOtherEnd=2,
                fDevTol=fDevTol,
                bMaximizeMinRadius=bMaximizeMinRadius,
                bDeformLines=False,
                bRebuildRationals=False,
                bEcho=False,
                bDebug=bDebug,
        )
        if rc is None:
            segs1.append(segA)
            segA = segB
            continue

        rgNurbsCrv1_A, rgNurbsCrv1_B = rc

        bSomeSegsReplaced = True
        
        #if iSegB == iSegB_toWatch:
        #    if rgNurbsCrv1_A: sc.doc.Objects.AddCurve(rgNurbsCrv1_A)
        #    if rgNurbsCrv1_B: sc.doc.Objects.AddCurve(rgNurbsCrv1_B)
        #    sc.doc.Views.Redraw(); return
        
        segs1.append(segA if rgNurbsCrv1_A is None else rgNurbsCrv1_A)

        segA = segB if rgNurbsCrv1_B is None else rgNurbsCrv1_B

    if c0.IsClosed:
        segs1[0] = segA
    else:
        segs1.append(segA) # Last curve.

    if not bSomeSegsReplaced:
        print "No segments were replaced."
        for c in segs0: c.Dispose()
        for c in segs1: c.Dispose()
        return
    else:
        s  = "Segments skipped at least once per joint: "
        s += str(sorted(iSegs_Skipped))
        print s

    print "Segment count: {} -> {}".format(len(segs0), len(segs1))

    rc = rg.Curve.JoinCurves(segs1)
    for c in segs0: c.Dispose()
    for c in segs1: c.Dispose()

    if bDebug: sEval='rc'; print sEval+': ',eval(sEval)

    if len(rc) > 1:
        print "{} curves returned from JoinCurve.".format(len(rc))
        for c in rc:
            sc.doc.Objects.AddCurve(c)
        sc.doc.Views.Redraw()
        return

    c1 = rc[0]
    
    if isinstance(c1, rg.PolyCurve):
        print "Segments: {}".format(c1.SegmentCount)
    else:
        print "New curve is not a PolyCurve.  Modification will be skipped."
        return

    return c1


def processCurveObject(objref, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bC1 = getOpt('bC1')
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bSkipIfAlreadyAtContinuity = getOpt('bSkipIfAlreadyAtContinuity')
    fDevTol = getOpt('fDevTol')
    bMaximizeMinRadius = getOpt('bMaximizeMinRadius')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    g0 = objref.ObjectId
    c0 = objref.Curve()

    if not isinstance(c0, rg.PolyCurve):
        if bEcho: print "Curve is not a PolyCurve."
        c0.Dispose()
        return


    c1 = processPolyCurve(
            rgPolyCrv0=c0,
            bC1=bC1,
            fAngleTol_Deg=fAngleTol_Deg,
            fDevTol=fDevTol,
            bMaximizeMinRadius=bMaximizeMinRadius,
            bReplace=bReplace,
            bEcho=bEcho,
            bDebug=bDebug,
    )


    def processResults(gCrv0_ToMod, rgCrv0_ToMod, rgNurbsCrv1):

        if bReplace:
            if sc.doc.Objects.Replace(gCrv0_ToMod, rgNurbsCrv1):
                if bEcho: print "Curve was replaced."
                g1 = gCrv0_ToMod
            else:
                g1 = None
        else:
            g1 = sc.doc.Objects.AddCurve(rgNurbsCrv1)
            if g1 == Guid.Empty: g1 = None 
            if not g1:
                if bEcho: print "Curve was added."
        sc.doc.Views.Redraw()
    
        rc = xCurve_inflections.getInflectionParameters(rgCrv0_ToMod)
        i_Infl_Ct_Pre = len(rc) if rc else None
        rc = xCurve_inflections.getInflectionParameters(rgNurbsCrv1)
        i_Infl_Ct_Post = len(rc) if rc else None
        s = "InflectionCount:{}->{}".format(i_Infl_Ct_Pre, i_Infl_Ct_Post)
    
        rc = xCurve_radiusMinima.getMinimumRadius(rgCrv0_ToMod)
        fRadius_Min_Original = rc if rc is not None else None
        rc = xCurve_radiusMinima.getMinimumRadius(rgNurbsCrv1)
        fRadius_Min_New = rc if rc is not None else None
        sRadius_Min_Original = '{:.{}f}'.format(fRadius_Min_Original, sc.doc.ModelDistanceDisplayPrecision)
        sRadius_Min_New = '{:.{}f}'.format(fRadius_Min_New, sc.doc.ModelDistanceDisplayPrecision)
        s += "  MinimumRadius:{}->{}".format(sRadius_Min_Original, sRadius_Min_New)
    
        fDev = getMaximumDeviation(rgCrv0_ToMod, rgNurbsCrv1)
        s += "  Deviation:{}".format(formatDistance(fDev))
        print s


    if c1:
        processResults(g0, c0, c1)
        c1.Dispose()

    c0.Dispose()


def main():
    
    rc = getInput()
    if rc is None: return
    objref = rc[0]
    #(
    #        objref,
    #        iContinuity,
    #        iPreserveOtherEnd,
    #        bC1,
    #        bMaximizeMinRadius,
    #        fRadius_MinTol,
    #        bReplace,
    #        bEcho,
    #        bDebug,
    #) = rc

    if Opts.values['bDebug']:
        pass


    processCurveObject(
            objref=objref,
            fDevTol=Opts.values['fDevTol'] if Opts.values['bLimitCrvDev'] else None,
    )


if __name__ == '__main__': main()