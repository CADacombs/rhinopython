"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
170902: Created.
...
210730: Now, printing of all knot values only happens in Debug mode.
211002: Refactored main routines.  Bug fix.  Now, count of points added is printed.
211126: Now, interior knot intersection of surfaces are included, not just along the edges.
211226: Added dot output.  Now adds curves along surface isocurves at target knot.
        Integers of min. and max. continuities to mark replaced bool input.
        Removed an option for surface input.
211127: Minor bug fixes.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iGCont_min'; keys.append(key)
    values[key] = 0
    names[key] = 'MinG'
    riOpts[key] = ri.Custom.OptionInteger(values[key], lowerLimit=0, upperLimit=10)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iGCont_max'; keys.append(key)
    values[key] = 1
    names[key] = 'MaxG'
    riOpts[key] = ri.Custom.OptionInteger(values[key], lowerLimit=0, upperLimit=10)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDot'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHt'; keys.append(key)
    values[key] = 12
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
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

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue <= 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get Curves or BrepFaces with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curve or faces")

    go.GeometryFilter = rd.ObjectType.Curve | rd.ObjectType.Surface
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    bPreselectedObjsChecked = False

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iGCont_min')
        addOption('iGCont_max')
        addOption('bDot')
        if Opts.values['bDot']:
            addOption('iDotHt')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        #        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
        #            bPreselectedObjsChecked = True
        #            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        #            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])

        if res == ri.GetResult.Number:
            if abs(int(go.Number())) > 10:
                print("Numeric input is invalid.  "
                      "Geometric Continuity checking is between 0 and 10.")
                continue
            Opts.setValue('iGCont_max')
            if Opts.riOpts['iGCont_max'].CurrentValue < Opts.riOpts['iGCont_min'].CurrentValue:
                Opts.riOpts['iGCont_min'].CurrentValue = Opts.riOpts['iGCont_max'].CurrentValue
                Opts.setValue('iGCont_min')
            continue

        # An option was selected.
        if go.Option().Index == idxs_Opt['iGCont_min']:
            Opts.setValue('iGCont_min')
            if Opts.riOpts['iGCont_min'].CurrentValue > Opts.riOpts['iGCont_max'].CurrentValue:
                Opts.riOpts['iGCont_max'].CurrentValue = Opts.riOpts['iGCont_min'].CurrentValue
                Opts.setValue('iGCont_max')
            continue

        if go.Option().Index == idxs_Opt['iGCont_max']:
            Opts.setValue('iGCont_max')
            if Opts.riOpts['iGCont_max'].CurrentValue < Opts.riOpts['iGCont_min'].CurrentValue:
                Opts.riOpts['iGCont_min'].CurrentValue = Opts.riOpts['iGCont_max'].CurrentValue
                Opts.setValue('iGCont_min')
            continue

        if go.Option().Index in idxs_Opt:
            Opts.setValue(key, go.Option().CurrentListOptionIndex)


def addPointsAtNurbsCrvKnots(nc, iGCont_min=0, iGCont_max=1, bDot=False, iDotHt=11, bEcho=True, bDebug=False):
    """
    """

    if nc.SpanCount == 1:
        nc.Dispose()
        return

    iK = 0 if nc.IsClosed else nc.Degree

    iCt_Pt = 0

    while iK < (nc.Knots.Count - nc.Degree):
        sc.escape_test()

        tK = nc.Knots[iK]
        mK = nc.Knots.KnotMultiplicity(iK)

        if (nc.Degree - mK) < iGCont_min or iGCont_max < (nc.Degree - mK):
            iK += mK
            continue

        if bDebug: print(tK)

        pt = nc.PointAt(tK)

        if bDot:
            dot = rg.TextDot("G{}\nknot".format(nc.Degree-mK), pt)
            dot.FontHeight = iDotHt
            sc.doc.Objects.AddTextDot(dot)
        else:
            sc.doc.Objects.AddPoint(pt)
        iCt_Pt += 1

        iK += mK

    return iCt_Pt


def addPointsAtNurbsSrfKnots(ns, iGCont_min=0, iGCont_max=1, bDot=False, iDotHt=12, bEcho=True, bDebug=False):
    """
    """

    def epsEquals(a, b, epsilon=2**(-52)):
        return abs(a-b) <= epsilon


    attr_Red = rd.ObjectAttributes()
    attr_Red.LayerIndex = sc.doc.ActiveDoc.Layers.CurrentLayerIndex
    attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr_Red.ObjectColor = Color.Red

    attr_Green = attr_Red.Duplicate()
    attr_Green.ObjectColor = Color.Lime


    gPts_Out = []

    tKs = [], []
    mKs = [], []

    # Accumulate data by traversing the knots once along U then once along V,
    # not through all intersections since that would be redundant.
    for iDir in 0,1:

        if iDir == 0:
            knots = ns.KnotsU
            degree = ns.Degree(0)
            domain = ns.Domain(0)
            pointAt = lambda a, b: ns.PointAt(a, b)
        else:
            knots = ns.KnotsV
            degree = ns.Degree(1)
            domain = ns.Domain(1)
            pointAt = lambda b, a: ns.PointAt(a, b)

        iK = 0 if ns.IsClosed(iDir) else degree

        while iK < (knots.Count - degree):
            sc.escape_test()

            tK = knots[iK]

            if tK < ns.Domain(iDir).T0:
                print("Skipped {}".format(tK))
                iK += 1
                continue

            if ns.IsClosed(iDir) and epsEquals(tK, ns.Domain(iDir).T1):
                break

            if tK > ns.Domain(iDir).T1:
                print("Skipped {}".format(tK))
                break

            mK = knots.KnotMultiplicity(iK)

            if (degree - mK) < iGCont_min or iGCont_max < (degree - mK):
                iK += 1
                continue

            tKs[iDir].append(tK)
            mKs[iDir].append(mK)

            if bDot:
                if iDir == 0:
                    pt = ns.PointAt(tK, ns.Domain(1).T0)
                else:
                    pt = ns.PointAt(ns.Domain(0).T0, tK)

                dot = rg.TextDot("G{}\nknot".format(degree-mK), pt)
                dot.FontHeight = iDotHt
                gOut = sc.doc.Objects.AddTextDot(
                    dot,
                    attributes=attr_Red if iDir == 1 else attr_Green)
                if gOut != gOut.Empty: gPts_Out.append(gOut)

            if iDir == 1:
                gOut = sc.doc.Objects.AddCurve(ns.IsoCurve(0, tK), attributes=attr_Red)
            else:
                gOut = sc.doc.Objects.AddCurve(ns.IsoCurve(1, tK), attributes=attr_Green)
            if gOut != gOut.Empty: gPts_Out.append(gOut)


            iK += mK


    def isCornerPoint(ns, tU, tV):
        return (
            (epsEquals(tU, ns.Domain(0).T0) or epsEquals(tU, ns.Domain(0).T1)) and
            (epsEquals(tV, ns.Domain(1).T0) or epsEquals(tV, ns.Domain(1).T1)) )


    # Traverse all the knot isocurve intersections per accumulated data.
    #for iU, (tU, mU) in enumerate(zip(tKs[0], mKs[0])):
    #    for iV, (tV, mV) in enumerate(zip(tKs[1], mKs[0])):

    #        if isCornerPoint(ns, tU, tV): continue

    #        pt = ns.PointAt(tU, tV)

    #        gOut = sc.doc.Objects.AddPoint(pt)

    #        if gOut != gOut.Empty: gPts_Out.append(gOut)

    return len(gPts_Out)


def main():

    rc = getInput()
    if rc is None: return
    (
        objrefs,
        iGCont_min,
        iGCont_max,
        bDot,
        iDotHt,
        bEcho,
        bDebug,
        ) = rc
    
    if iGCont_min > iGCont_max:
        print("Min. G continuity target is greater than the max. G continuity target."
              "  Script canceled.")
        return
    
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    iCt_Pt_All = 0
    
    for objref in objrefs:
        rgObj = objref.Curve()

        if rgObj is None:
            rgObj = objref.Surface()
            if rgObj is None: contiue
        
        if isinstance(rgObj, rg.Curve):
            rgNurbsCrv = rgObj.Duplicate() if isinstance(rgObj, rg.NurbsCurve) else rgObj.ToNurbsCurve()
            
            iCt_Pts = addPointsAtNurbsCrvKnots(
                rgNurbsCrv,
                iGCont_min=iGCont_min,
                iGCont_max=iGCont_max,
                bDot=bDot,
                iDotHt=iDotHt,
                bEcho=bEcho,
                bDebug=bDebug)

            rgNurbsCrv.Dispose()

        elif isinstance(rgObj, rg.Surface):
            rgNurbsSrf = rgObj if isinstance(rgObj, rg.NurbsSurface) else rgObj.ToNurbsSurface()

            iCt_Pts = addPointsAtNurbsSrfKnots(
                rgNurbsSrf,
                iGCont_min=iGCont_min,
                iGCont_max=iGCont_max,
                bDot=bDot,
                iDotHt=iDotHt,
                bEcho=bEcho,
                bDebug=bDebug)

            rgNurbsSrf.Dispose()

        if iCt_Pts is None: continue

        iCt_Pt_All += iCt_Pts


    if bEcho:
        print("{} objects added.".format(iCt_Pt_All))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()