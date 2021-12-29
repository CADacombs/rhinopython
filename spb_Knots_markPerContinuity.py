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
211127: Added preview.  Minor bug fixes.
211229: Bug fixes in getInput.
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


    key = 'iGCont_max'; keys.append(key)
    values[key] = 1
    names[key] = 'MaxGContToMark'
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

    key = 'bAddObjs'; keys.append(key)
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


def getInput(gas):
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

    go.AcceptNothing(True)

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:

        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('iGCont_max')
        addOption('bDot')
        if Opts.values['bDot']:
            addOption('iDotHt')
        addOption('bAddObjs')
        addOption('bEcho')
        addOption('bDebug')

        if Opts.values['bAddObjs']:
            go.SetCommandPromptDefault("Add objs")
        else:
            go.SetCommandPromptDefault("Cancel")

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

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return False

        if res == ri.GetResult.Object:
            geoms_In = []
            for objref in go.Objects():
                geom = objref.Curve()
                if geom is None:
                    geom = objref.Surface()
                if geom: geoms_In.append(geom)
            go.Dispose()
            return geoms_In

        if res == ri.GetResult.Number:
            iNum = abs(int(go.Number()))
            if iNum > 10:
                print("Numeric input is invalid.  "
                        "Geometric continuity checking is between 0 and 10.")
                continue
            Opts.riOpts['iGCont_max'].CurrentValue = iNum
            Opts.setValue('iGCont_max')
            if gas:
                return True
            continue

        # An option was selected.
        idx = go.Option().Index
        if idx in idxs_Opt.values():
            key = idxs_Opt.keys()[idxs_Opt.values().index(idx)]
            Opts.setValue(key, go.Option().CurrentListOptionIndex)
            if key in ('iDotHt', 'bAddObjs', 'bEcho', 'bDebug'): continue
            if gas: return True


def addPointsAtNurbsCrvKnots(nc, iGCont_max=1, bDot=False, iDotHt=11, bEcho=True, bDebug=False):
    """
    """

    if nc.SpanCount == 1:
        nc.Dispose()
        return

    iK = 0 if nc.IsClosed else nc.Degree

    gas = [] # Each is tuple(rg, attr)

    while iK < (nc.Knots.Count - nc.Degree):
        sc.escape_test()

        tK = nc.Knots[iK]
        mK = nc.Knots.KnotMultiplicity(iK)

        if iGCont_max < (nc.Degree - mK):
            iK += mK
            continue

        if bDebug: print(tK)

        pt = nc.PointAt(tK)

        if bDot:
            dot = rg.TextDot("G{}\nknot".format(nc.Degree-mK), pt)
            dot.FontHeight = iDotHt
            gas.append((dot, None))
            #sc.doc.Objects.AddTextDot(dot)
        else:
            gas.append((pt, None))
            #sc.doc.Objects.AddPoint(pt)

        iK += mK

    return gas


def addPointsAtNurbsSrfKnots(ns, iGCont_max=1, bDot=False, iDotHt=12, bEcho=True, bDebug=False):
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


    gas = [] # Each is tuple(rg, attr)

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

            if iGCont_max < (degree - mK):
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
                gas.append((dot, attr_Red if iDir == 1 else attr_Green))
                #gOut = sc.doc.Objects.AddTextDot(
                #    dot,
                #    attributes=attr_Red if iDir == 1 else attr_Green)
                #if gOut != gOut.Empty: gPts_Out.append(gOut)

            if iDir == 1:
                gas.append((ns.IsoCurve(0, tK), attr_Red))
                #gOut = sc.doc.Objects.AddCurve(ns.IsoCurve(0, tK), attributes=attr_Red)
            else:
                gas.append((ns.IsoCurve(1, tK), attr_Green))
                #gOut = sc.doc.Objects.AddCurve(ns.IsoCurve(1, tK), attributes=attr_Green)
            #if gOut != gOut.Empty: gPts_Out.append(gOut)


            iK += mK


    #def isCornerPoint(ns, tU, tV):
    #    return (
    #        (epsEquals(tU, ns.Domain(0).T0) or epsEquals(tU, ns.Domain(0).T1)) and
    #        (epsEquals(tV, ns.Domain(1).T0) or epsEquals(tV, ns.Domain(1).T1)) )


    # Traverse all the knot isocurve intersections per accumulated data.
    #for iU, (tU, mU) in enumerate(zip(tKs[0], mKs[0])):
    #    for iV, (tV, mV) in enumerate(zip(tKs[1], mKs[0])):

    #        if isCornerPoint(ns, tU, tV): continue

    #        pt = ns.PointAt(tU, tV)

    #        gOut = sc.doc.Objects.AddPoint(pt)

    #        if gOut != gOut.Empty: gPts_Out.append(gOut)

    return gas


def createGeoms(geoms_In, iGCont_max=1, bDot=False, iDotHt=12, bEcho=True, bDebug=False):
    """
    """

    outs = []

    for rgObj in geoms_In:

        if isinstance(rgObj, rg.Curve):
            rgNurbsCrv = rgObj.Duplicate() if isinstance(rgObj, rg.NurbsCurve) else rgObj.ToNurbsCurve()
            
            rc = addPointsAtNurbsCrvKnots(
                rgNurbsCrv,
                iGCont_max=iGCont_max,
                bDot=bDot,
                iDotHt=iDotHt,
                bEcho=bEcho,
                bDebug=bDebug)

            rgNurbsCrv.Dispose()

            if rc: outs.extend(rc)

        elif isinstance(rgObj, rg.Surface):
            rgNurbsSrf = rgObj if isinstance(rgObj, rg.NurbsSurface) else rgObj.ToNurbsSurface()

            rc = addPointsAtNurbsSrfKnots(
                rgNurbsSrf,
                iGCont_max=iGCont_max,
                bDot=bDot,
                iDotHt=iDotHt,
                bEcho=bEcho,
                bDebug=bDebug)

            rgNurbsSrf.Dispose()

            if rc: outs.extend(rc)

    return outs


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.gas = None
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 1

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if not self.gas: return

        for ga in self.gas:
            geom, attr = ga
            self.bbox = geom.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):
        if not self.gas: return

        for ga in self.gas:
            geom, attr = ga

            if attr is None:
                objColor = sc.doc.Layers.CurrentLayer.Color
            else:
                objColor = attr.ObjectColor

            if isinstance(geom, rg.TextDot):
                dot = geom
                textColor = Color.Black if objColor != Color.Black else Color.White
                drawEventArgs.Display.DrawDot(
                    worldPosition=dot.Point,
                    text=dot.Text,
                    dotColor=objColor,
                    textColor=textColor)
            elif isinstance(geom, rg.Curve):
                crv = geom
                drawEventArgs.Display.DrawCurve(
                    curve=crv,
                    color=objColor,
                    thickness=self.crv_thk)
            elif isinstance(geom, rg.Point3d):
                pt = geom
                drawEventArgs.Display.DrawPoint(
                    point=pt)


def main():

    sk_conduit = 'conduit({})'.format(__file__)
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit

    conduit.Enabled = False
    sc.doc.Views.Redraw()


    gas = None # list(tuple(rg, attr))


    while True:
        rc = getInput(gas)
        if rc is None:
            conduit.Enabled = False
            return

        iGCont_max = Opts.values['iGCont_max']
        bDot = Opts.values['bDot']
        iDotHt = Opts.values['iDotHt']
        bAddObjs = Opts.values['bAddObjs']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        if not rc:
            conduit.Enabled = False
            break


        try:
            iter(rc)
            geoms_In = rc
        except:
            pass

        rc = createGeoms(
            geoms_In=geoms_In,
            iGCont_max=iGCont_max,
            bDot=bDot,
            iDotHt=iDotHt,
            bEcho=bEcho,
            bDebug=bDebug)
        if rc is None:
            conduit.Enabled = False
            conduit = None
            return

        gas = rc

        sc.doc.Objects.UnselectAll()

        conduit.gas = gas
        conduit.Enabled = True
        sc.doc.Views.Redraw()


    if bAddObjs:
        gOuts = []
        for geom, attr in gas:
            if isinstance(geom, rg.Point3d):
                gOut = sc.doc.Objects.AddPoint(geom)
            else:
                gOut = sc.doc.Objects.Add(geom, attr)
            if gOut != gOut.Empty:
                gOuts.append(gOut)
        if bEcho:
            print("{} objects added.".format(len(gOuts)))

        sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()