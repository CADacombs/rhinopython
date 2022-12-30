"""
Alternative to _Fair:
  1. When LimitDev is True, results closer to the distance tolerance are achieved.
  2. The clamp order is not limited to 0, a la Rhino 7, or 2, a la Rhino 8.
  3. NURBS curve structures for degrees other than 3 are preserved.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
170814-19: Created.
...
211113: Script now requires Rhino 7 or higher.  Replaced an import with a local function.
221229-30: Modified main function routine.  Refactored.
        Due to a discovery today, clamp orders are no longer limited to 2.

TODO: Investigate results of non-3-degree input.
"""

import Rhino
import Rhino.DocObjects as rd
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


    key = 'bLimitDev'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fAngleTol_Deg'; keys.append(key)
    values[key] = 50.0 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTolForKinks'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEqualEndClamps'; keys.append(key)
    values[key] = True
    names[key] = 'EqualEndClampCurvatureOrder'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClampDerivativeOrder'; keys.append(key)
    values[key] = 1
    names[key] = 'Order'
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClampStartDerivativeOrder'; keys.append(key)
    values[key] = 1
    names[key] = 'StartOrder'
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iClampEndDerivativeOrder'; keys.append(key)
    values[key] = 1
    names[key] = 'EndOrder'
    riOpts[key] = ri.Custom.OptionInteger(initialValue=values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace_NotAdd'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='Add', onValue='Replace')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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

        if key == 'fDevTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key == 'iClampDerivativeOrder':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = 0

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


def getInput():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    #go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False
    
    print("For ClampOrder(s), value represents maximum index"
          " (from relative end) of control points whose locations are not modified."
          "\nFor example, 2 means that first 3 controls point locations are preserved.")

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bLimitDev')
        if Opts.values['bLimitDev']:
            addOption('fDevTol')
            addOption('fAngleTol_Deg')
        addOption('bEqualEndClamps')
        if Opts.values['bEqualEndClamps']:
            addOption('iClampDerivativeOrder')
        else:
            addOption('iClampStartDerivativeOrder')
            addOption('iClampEndDerivativeOrder')
        addOption('bReplace_NotAdd')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            if Opts.values['bLimitDev']:
                key = 'fDevTol'
            elif Opts.values['bEqualEndClamps']:
                key = 'iClampDerivativeOrder'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getPreselectedCurves():
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    if gObjs_Preselected:
        gCrvs_Preselected = []
        iter = rd.ObjectEnumeratorSettings()
        iter.NormalObjects = True
        iter.LockedObjects = False
        iter.IncludeLights = False
        iter.IncludeGrips = False
        for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
            if rdRhinoObject.Id in gObjs_Preselected:
                if rdRhinoObject.ObjectType == rd.ObjectType.Curve:
                    gCrvs_Preselected.append(rdRhinoObject.Id)
        if gCrvs_Preselected:
            if Opts.values['bEcho']:
                s  = "({} curves".format(len(gCrvs_Preselected))
                s += " were preselected.)"
                print(s)
            return tuple(gCrvs_Preselected)


def getDistancesBetweenCurves(crvA, crvB):
    rc = rg.Curve.GetDistancesBetweenCurves(
            crvA, crvB, 0.1*sc.doc.ModelAbsoluteTolerance)

    if not rc[0]:
        raise Exception("GetDistancesBetweenCurves returned None.")
        return None

    return rc[1]


def formatDistance(fDistance, iPrecision=None):
    if iPrecision is None: iPrecision = sc.doc.ModelDistanceDisplayPrecision

    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, iPrecision)


def fairCurve_NoDevLimit(rgCrv_In, iClampStartDerivativeOrder, iClampEndDerivativeOrder, bDebug=False):
    """
    """

    nc_Start = rgCrv_In.ToNurbsCurve()

    if bDebug:
        sEval = "nc_Start.Points.Count"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "nc_Start.Degree"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "iClampStartDerivativeOrder"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "iClampEndDerivativeOrder"; print("{}: {}".format(sEval, eval(sEval)))


    if nc_Start.Points.Count < (iClampStartDerivativeOrder + iClampEndDerivativeOrder + 2):
        nc_Start.Dispose()
        return None, "Clamping continuity is too high for curve's point count."

    #if nc_Start.Points.Count < (iClampStartDerivativeOrder + iClampEndDerivativeOrder) + 2*nc_Start.Degree:
    #    nc_Start.Dispose()
    #    return None, "Clamping continuity is too high for curve's point count."



    bb = rgCrv_In.GetBoundingBox(accurate=True)
    distanceTolerance = bb.Diagonal.Length
    angleTolerance = Rhino.RhinoMath.ToRadians(90.0)

    if nc_Start.Degree == 3:
        nc_Prev = nc_Start.DuplicateCurve()
    else:
        nc_Prev = rg.NurbsCurve.Create(
            periodic=False,
            degree=3,
            points=[cp.Location for cp in nc_Start.Points])
        #sc.doc.Objects.AddCurve(nc_Prev); sc.doc.Views.Redraw(); 1/0


    iWhile = 0

    while True:
        sc.escape_test()

        nc_Res = nc_Prev.Fair(
            distanceTolerance=1000*distanceTolerance,
            angleTolerance=angleTolerance,
            clampStart=iClampStartDerivativeOrder,
            clampEnd=iClampEndDerivativeOrder,
            iterations=1000)

        if nc_Res is None:
            if iWhile == 0:
                return None, "Clamping continuity may be too high for curve's point count."
            nc_Res = nc_Prev
            break

        if nc_Res.EpsilonEquals(nc_Prev, epsilon=Rhino.RhinoMath.ZeroTolerance):
            if bDebug:
                sEval = "iWhile"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dev"; print("{}: {}".format(sEval, eval(sEval)))
            nc_Prev.Dispose()
            break

        nc_Prev.Dispose()
        nc_Prev = nc_Res

        iWhile += 1

    if nc_Start.Degree == 3:
        nc_Start.Dispose()
        return nc_Res, None
    else:
        nc_Out = nc_Start
        for iCp in range(nc_Out.Points.Count):
            nc_Out.Points.SetPoint(
                iCp,
                point=nc_Res.Points[iCp].Location,
                weight=nc_Res.Points[iCp].Weight)
        return nc_Out, None


def fairCurve_DevLimit(rgCrv_In, fDevTol, fAngleTol_Deg, iClampStartDerivativeOrder, iClampEndDerivativeOrder, bDebug=False):
    """
    """

    nc_Start = rgCrv_In.ToNurbsCurve()

    if nc_Start.Points.Count < (iClampStartDerivativeOrder + iClampEndDerivativeOrder + 2):
        nc_Start.Dispose()
        return None, "Clamping continuity is too high for curve's point count."

    #if nc_Start.Points.Count < (iClampStartDerivativeOrder + iClampEndDerivativeOrder) + 2*nc_Start.Degree:
    #    nc_Start.Dispose()
    #    return None, "Clamping continuity is too high for curve's point count."

    nc_Prev = rgCrv_In.ToNurbsCurve()

    distanceTolerance = fDevTol
    angleTolerance = Rhino.RhinoMath.ToRadians(fAngleTol_Deg)



    nc_Res = nc_Start.Fair(
        distanceTolerance=2.0*fDevTol,
        angleTolerance=angleTolerance,
        clampStart=iClampStartDerivativeOrder,
        clampEnd=iClampEndDerivativeOrder,
        iterations=1)

    if nc_Res is None:
        return None, "Clamping continuity may be too high for curve's point count."

    dev = getDistancesBetweenCurves(nc_Start, nc_Res)

    if dev == fDevTol:
        return nc_Res, None

    if dev > fDevTol:
        iterations = 1

    else:
        # Try higher iterations to get dev closer to fDevTol.

        iWhile = 0

        iter_L = 1
        iter_M = None
        iter_H = 2**16


        dev = None

        while True:
            sc.escape_test()

            if bDebug: sEval = "iWhile"; print("{}: {}".format(sEval, eval(sEval)))

            iter_M = (iter_L + iter_H) // 2

            if iter_M in (iter_L, iter_H):
                if bDebug:
                    sEval = "iter_L"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "iter_M"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "iter_H"; print("{}: {}".format(sEval, eval(sEval)))
                iterations = iter_L
                break

            dev_Last = dev

            nc_Res = nc_Start.Fair(
                distanceTolerance=distanceTolerance,
                angleTolerance=angleTolerance,
                clampStart=iClampStartDerivativeOrder,
                clampEnd=iClampEndDerivativeOrder,
                iterations=iter_M)

            if nc_Res is None:
                raise Exception("Fair returned None.")

            dev = getDistancesBetweenCurves(nc_Start, nc_Res)
            if bDebug: sEval = "dev"; print("{}: {}".format(sEval, eval(sEval)))

            if dev_Last and abs(dev - dev_Last) < 1e-6:
                # No improvement.
                iterations = iter_H
                break


            if dev <= fDevTol:
                iter_L = iter_M
            else:
                iter_H = iter_M

            nc_Res.Dispose()

            iWhile += 1
        else:
            iterations = iter_H

        if bDebug: sEval = "iterations"; print("{}: {}".format(sEval, eval(sEval)))









    dist_L = fDevTol
    dist_M = None
    dist_H = 2.0*fDevTol



    nc_Res = nc_Start.Fair(
        distanceTolerance=dist_H,
        angleTolerance=angleTolerance,
        clampStart=iClampStartDerivativeOrder,
        clampEnd=iClampEndDerivativeOrder,
        iterations=1)


    if nc_Res is None:
        raise Exception("Fair returned None.")

    dev = getDistancesBetweenCurves(nc_Start, nc_Res)

    if dev <= fDevTol:
        return nc_Res, None



    while True:
        sc.escape_test()

        dist_M = (dist_L + dist_H) / 2.0

        if dist_M in (dist_L, dist_H):
            if bDebug:
                sEval = "dist_L"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_M"; print("{}: {}".format(sEval, eval(sEval)))
                sEval = "dist_H"; print("{}: {}".format(sEval, eval(sEval)))
            break


        #sEval = "dist_M"; print("{}: {}".format(sEval, eval(sEval)))

        nc_Res = nc_Start.Fair(
            distanceTolerance=dist_M,
            angleTolerance=angleTolerance,
            clampStart=iClampStartDerivativeOrder,
            clampEnd=iClampEndDerivativeOrder,
            iterations=1)

        if nc_Res is None:
            raise Exception("Fair returned None.")

        dev = getDistancesBetweenCurves(nc_Start, nc_Res)

        if dev > fDevTol:
            dist_H = dist_M
        else:
            dist_L = dist_M

        nc_Res.Dispose()



    distanceTolerance = dist_L


    nc_Res = nc_Start.Fair(
        distanceTolerance=distanceTolerance,
        angleTolerance=angleTolerance,
        clampStart=iClampStartDerivativeOrder,
        clampEnd=iClampEndDerivativeOrder,
        iterations=iterations)

    return nc_Res, None



def isKnotVectorUniform(knots):
    return (
        (knots.KnotStyle == rg.KnotStyle.Uniform) or
        (knots.KnotStyle == rg.KnotStyle.QuasiUniform) or
        (
            (knots.KnotStyle == rg.KnotStyle.PiecewiseBezier) and
            knots.Count == knots.KnotMultiplicity(0) * 2)
        )


def printOutput(rgCrv, nc_End):

    print("Original curve is a {}.".format(rgCrv.GetType().Name))

    nc_Start = rgCrv.ToNurbsCurve()

    if nc_End:
        s = "Prop:I-O"
        s += "  {}:{}-{}".format("Deg", nc_Start.Degree, nc_End.Degree)
        s += "  {}:{}-{}".format("PtCt", nc_Start.Points.Count, nc_End.Points.Count)
        s += "  {}:{}-{}".format("IsUniform",
                str(isKnotVectorUniform(nc_Start.Knots))[0],
                str(isKnotVectorUniform(nc_End.Knots))[0])
        s += "  {}:{}-{}".format("IsRational",
                str(nc_Start.IsRational)[0],
                str(nc_End.IsRational)[0])
        s += "  {}:{}-{}".format("IsClosed",
                str(nc_Start.IsClosed)[0],
                str(nc_End.IsClosed)[0])
        if nc_Start.IsClosed or nc_End.IsClosed:
            s += "  {}:{}-{}".format("IsPeriodic",
                    str(nc_Start.IsPeriodic)[0],
                    str(nc_End.IsPeriodic)[0])

        crv_dev = getDistancesBetweenCurves(rgCrv, nc_End)

        if crv_dev:
            s += "  Deviation: {0:.{1}f}".format(
                    crv_dev, sc.doc.ModelDistanceDisplayPrecision)
        else:
            s += "  Curve deviation cannot be calculated!"
        
        print(s)
    
    else:
        sEval = "rgNurbsCrv1"; print("{}: {}".format(sEval, eval(sEval)))
        if nc_End is None:
            print("Invalid input for fairCurve?")
        else: # False
            print("Curve could not be found for entered parameters.")
            print("Input curve CpCt:{},  IsUniform:{},  IsPeriodic:{}".format(
                    nc_Start.Points.Count,
                    str(isKnotVectorUniform(nc_Start.Knots))[0],
                    nc_Start.IsPeriodic))

    nc_Start.Dispose()


def main(bEcho=False, bDebug=False):

    gCs_Preselected = getPreselectedCurves()

    objrefs = getInput()
    if objrefs is None: return

    fDevTol = Opts.values['fDevTol'] if Opts.values['bLimitDev'] else None
    fAngleTol_Deg = Opts.values['fAngleTol_Deg']
    if Opts.values['bEqualEndClamps']:
        iClampStartDerivativeOrder = Opts.values['iClampDerivativeOrder']
        iClampEndDerivativeOrder = Opts.values['iClampDerivativeOrder']
    else:
        iClampStartDerivativeOrder = Opts.values['iClampStartDerivativeOrder']
        iClampEndDerivativeOrder = Opts.values['iClampEndDerivativeOrder']
    bReplace_NotAdd = Opts.values['bReplace_NotAdd']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    gCs_Replaced = []
    gCs_Added = []
    sLogs = []
    
    for objref in objrefs:
        rdC_In = objref.Object()

        rgC_In = rdC_In.Geometry

        if isinstance(rgC_In, rg.LineCurve):
            print("LineCurve skipped.")
            continue
        if rgC_In.IsLinear():
            print("Linear curve skipped.")
            continue

        nc_Start = rgC_In.ToNurbsCurve()


        crv_dev = None
        
        if fDevTol is None:
            nc_Res, sLog = fairCurve_NoDevLimit(
                rgCrv_In=nc_Start,
                iClampStartDerivativeOrder=iClampStartDerivativeOrder,
                iClampEndDerivativeOrder=iClampEndDerivativeOrder,
                bDebug=bDebug)
        else:
            nc_Res, sLog = fairCurve_DevLimit(
                rgCrv_In=nc_Start,
                fDevTol=fDevTol,
                fAngleTol_Deg=fAngleTol_Deg,
                iClampStartDerivativeOrder=iClampStartDerivativeOrder,
                iClampEndDerivativeOrder=iClampEndDerivativeOrder,
                bDebug=bDebug)

        if sLog is not None:
            sLogs.append(sLog)

        if nc_Res is None:
            pass
        elif rg.NurbsCurve.EpsilonEquals(nc_Res, nc_Start, epsilon=Rhino.RhinoMath.ZeroTolerance):
            print("No change in curve geometry within {}.".format(
                formatDistance(Rhino.RhinoMath.ZeroTolerance)))
            continue
        else:
            if bReplace_NotAdd:
                if sc.doc.Objects.Replace(rdC_In.Id, nc_Res):
                    gCs_Replaced.append(rdC_In.Id)
            else:
                gC_Added = sc.doc.Objects.AddCurve(nc_Res)
                if gC_Added == gC_Added.Empty:
                    print("Could not add new curve.")
                else:
                    gCs_Added.append(gC_Added)

            crv_dev = getDistancesBetweenCurves(nc_Start, nc_Res)

            printOutput(rgC_In, nc_Res)


    for sLog in set(sLogs):
        print("[{}]: {}".format(sLogs.count(sLog), sLog))

    if gCs_Preselected:
        [sc.doc.Objects.Select(objectId=_) for _ in gCs_Preselected]
    else:
        if gCs_Replaced:
            [sc.doc.Objects.Select(objectId=_) for _ in gCs_Replaced]
        if gCs_Added:
            [sc.doc.Objects.Select(objectId=_) for _ in gCs_Added]
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main(bEcho=0, bDebug=0)