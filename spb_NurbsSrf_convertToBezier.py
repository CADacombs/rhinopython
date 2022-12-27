"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
190204-05: Created.
...
190817: Changed face iteration status from count to percent.
190907: Added bLimitDeviation.  Now checks areas of surface against brep before attempting trim.
191031, 1109, 1206, 0202, 0224: Import-related update.
200515: Added options including changing surface degrees.
200518: Minor bug fix.
200527-28: Bug fix.  Refactored.  Import-related update.
        Now checks for and processes differently surfaces with Bezier patterns in one direction.
200616: Modified printed output.
200619: Import-related update.
210209: Refactored.  createSurface split into 4 functions.
210220: Bug fix.
210301: Bug fix.
210909: Bug fix.  Refactored same degree routine.
211113: Replaced an import with a local function.
220701: Refactored.
221227: WIP: Refactored.

TODO:
    If bShrinkFirst and the Rebuild fails, try Rebuild before shrinking.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xBrep_getDistancesBetween2
import xBrepFace
import xBrepObject


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bLimitDeviation'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMatchInSrfDegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg1'; keys.append(key)
    values[key] = True
    names[key] = 'TryDeg1'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg2'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg3'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeg5'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSkipConical'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bShrinkFirst'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get brep faces without optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps and/or faces")

    go.GeometryFilter = rd.ObjectType.Brep | rd.ObjectType.Surface

    go.AcceptNumber(True, acceptZero=True)

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bLimitDeviation')
        if Opts.values['bLimitDeviation']:
            addOption('fDevTol')
        addOption('bMatchInSrfDegs')
        if not Opts.values['bMatchInSrfDegs']:
            addOption('bDeg1')
            addOption('bDeg2')
            addOption('bDeg3')
            addOption('bDeg5')
        addOption('bSkipConical')
        addOption('bShrinkFirst')
        addOption('bReplace')
        if Opts.values['bReplace']:
            addOption('bExtract')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()

            if Opts.values['bLimitDeviation']:
                fDevTol = Opts.values['fDevTol']
            else:
                fDevTol = None

            if Opts.values['bMatchInSrfDegs']:
                iDegsToTry = None
            else:
                iDegsToTry = []
                if Opts.values['bDeg1']:
                    iDegsToTry.append(1)
                if Opts.values['bDeg2']:
                    iDegsToTry.append(2)
                if Opts.values['bDeg3']:
                    iDegsToTry.append(3)
                if Opts.values['bDeg5']:
                    iDegsToTry.append(5)
                iDegsToTry = tuple(iDegsToTry)

            return (
                objrefs,
                fDevTol,
                iDegsToTry,
                Opts.values['bSkipConical'],
                Opts.values['bShrinkFirst'],
                Opts.values['bReplace'],
                Opts.values['bExtract'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            if Opts.values['bLimitDeviation']:
                key = 'fDevTol'
                Opts.riOpts[key].CurrentValue = go.Number()
                Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def doesSurfaceContainConicSection(rgSrf):
    if isinstance(rgSrf, rg.RevSurface):
        return "RevSurface."
    if isinstance(rgSrf, rg.NurbsSurface):
        if rgSrf.IsRational:
            #return "rational NurbsSurface."

            c_UDom_Start = rgSrf.IsoCurve(1, rgSrf.Domain(0).T0)
            if isinstance(c_UDom_Start, rg.NurbsCurve):
                if c_UDom_Start.IsRational:
                    if c_UDom_Start.IsArc(1e-9):
                        return "NurbsSurface with circular arc NurbsCurve isocurve."
                    if c_UDom_Start.IsEllipse(1e-9):
                        return "NurbsSurface with elliptical NurbsCurve isocurve."
            c_VDom_Start = rgSrf.IsoCurve(0, rgSrf.Domain(1).T0)
            if isinstance(c_VDom_Start, rg.NurbsCurve):
                if c_VDom_Start.IsRational:
                    if c_VDom_Start.IsArc(1e-9):
                        return "NurbsSurface with circular arc NurbsCurve isocurve."
                    if c_VDom_Start.IsEllipse(1e-9):
                        return "NurbsSurface with elliptical NurbsCurve isocurve."
    if isinstance(rgSrf, rg.SumSurface):
        #ns = rgSrf.ToNurbsSurface()
        #if ns.IsRational:
        #    ns.Dispose()
        #    return "SumSurface with rational isocurve."

        c_UDom_Start = rgSrf.IsoCurve(1, rgSrf.Domain(0).T0)
        if isinstance(c_UDom_Start, rg.ArcCurve):
            return "SumSurface with conic section isocurve."
            return "SumSurface with ArcCurve isocurve at U T0."
        if isinstance(c_UDom_Start, rg.NurbsCurve):
            if c_UDom_Start.IsRational:
                if c_UDom_Start.IsEllipse(1e-9):
                    return "SumSurface with elliptical NurbsCurve isocurve."
                return "SumSurface with rational NurbsCurve isocurve at U T0."
        c_VDom_Start = rgSrf.IsoCurve(0, rgSrf.Domain(1).T0)
        if isinstance(c_VDom_Start, rg.ArcCurve):
            return "SumSurface with conic section isocurve."
            return "SumSurface with ArcCurve isocurve at V T0."
        if isinstance(c_VDom_Start, rg.NurbsCurve):
            if c_VDom_Start.IsRational:
                if c_VDom_Start.IsEllipse(1e-9):
                    return "SumSurface with elliptical NurbsCurve isocurve."
                return "SumSurface with rational NurbsCurve isocurve at V T0."
    return False


def isSrfValidForConversion(rgSrf, bSkipConical=True):
    """
    return tuple of bool and string (log) or None)
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return
    
    if isinstance(rgSrf, rg.PlaneSurface):
        return False, "Skipped Planesurface."

    if bSkipConical:
        rc = doesSurfaceContainConicSection(rgSrf)
        if rc:
            return False, "Skipped " + rc

    ns = rgSrf.ToNurbsSurface()
    
    # Skip if no interior knots.
    for iDir in 0, 1:
        spanCt = ns.SpanCount(iDir)
        if ns.IsPeriodic(iDir):
            iCt_Pt = ns.Points.CountU if iDir == 0 else ns.Points.CountV
            if spanCt > iCt_Pt - ns.Degree(iDir):
                ns.Dispose()
                return True, None
        else:
            if spanCt > 1:
                ns.Dispose()
                return True, None
    else:
        ns.Dispose()
        return False, "Skipped {} with single Bezier patch.".format(rgSrf.GetType().Name)


def createSurface_SameDegs(rgSrf_In, fDevTol=None, bDebug=False):
    """
    Create surface at same degrees as input.
    """

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return

    ns_fromIn = rgSrf_In.ToNurbsSurface()


    bSkipRebuildU = ns_fromIn.SpanCount(0)
    bSkipRebuildV = ns_fromIn.SpanCount(1)

    if bSkipRebuildU and bSkipRebuildV:
        return None, "Surface is already a Bezier patch."

    nss_Res = []
    devs = []
    sLogs = []


    def getDistancesBetweenBreps(nsA, nsB, bDebug=False):
        rc = xBrep_getDistancesBetween2.getDistancesBetweenBreps(
                nsB,
                nsA,
                bIncludeEdges=True,
                rgMeshParams=None,
                bFineMesh=False,
                bCalcBrepIntersection=False,
                bDebug=bDebug
                )

        return rc[1] if rc[0] else None


    def createSurface_RemoveInteriorKnots(ns_In, bDebug=False):
        """
        """

        # Duplicating because RemoveKnots modifies the NurbsSurface of the provided knot vector.
        ns_Out = ns_In.Duplicate()
    
        degrees = ns_Out.OrderU-1, ns_Out.OrderV-1
    
        knots = ns_Out.KnotsU, ns_Out.KnotsV
    
        # Remove all interior knots.
        # Notice that index 1 is one more than the maximum index to remove.
    
        for iDir in 0,1:
            if knots[iDir].Count == 2 * degrees[iDir]:
                # There are no interior knots in iDir direction.
                continue
        
            if not knots[iDir].RemoveKnots(
                    index0=degrees[iDir],
                    index1=knots[iDir].Count-degrees[iDir]
            ):
                ns_Out.Dispose()
                if bDebug: print("RemoveKnots failed.")
                return

        return ns_Out

    ns_Res = createSurface_RemoveInteriorKnots(ns_fromIn, bDebug)

    if ns_Res is not None:
        nss_Res.append(ns_Res)
        dev = getDistancesBetweenBreps(ns_Res, ns_fromIn, bDebug)
        devs.append(dev)
        sLogs.append(None)


    if (
        (bSkipRebuildU and ns_fromIn.Degree(1) == 3)
        or
        (bSkipRebuildV and ns_fromIn.Degree(0) == 3)
    ):
        def createSurface_RebuildOneDirection(ns_In, direction, bDebug=False):
            return ns_fromIn.RebuildOneDirection(
                direction=direction,
                pointCount=4,
                loftType=rg.LoftType.Normal,
                refitTolerance=0.5*sc.doc.ModelAbsoluteTolerance)

        ns_Res = createSurface_RebuildOneDirection(
            ns_fromIn,
            direction=1 if bSkipRebuildU else 0,
            bDebug=bDebug)

        if ns_Res is not None:
            nss_Res.append(ns_Res)
            dev = getDistancesBetweenBreps(ns_Res, ns_fromIn, bDebug)
            devs.append(dev)
            sLogs.append(None)


    def createSurface_Rebuild(ns_In):
        ns_Out = ns_In.Rebuild(
            uDegree=ns_In.Degree(0),
            vDegree=ns_In.Degree(1),
            uPointCount=ns_In.Degree(0)+1,
            vPointCount=ns_In.Degree(1)+1)
        if ns_Out is None:
            return None, "Rebuild returned None.  Check for internal knots near surface boundary."
        return ns_Out, None

    ns_Res, sLog = createSurface_Rebuild(ns_fromIn)

    if ns_Res is not None:
        nss_Res.append(ns_Res)
        dev = getDistancesBetweenBreps(ns_Res, ns_fromIn, bDebug)
        devs.append(dev)
        sLogs.append(sLog)
    elif bDebug:
        print(sLog)


    ns_fromIn.Dispose()


    if len(nss_Res) == 0:
        return None, sLogs


    #for ns in nss_Res:
    #    sc.doc.Objects.AddSurface(ns)
    #sc.doc.Views.Redraw()
    #return


    dev_Min = min(devs)
    idx_MinDev = devs.index(dev_Min)

    ns_MinDev = nss_Res[idx_MinDev]


    for i, ns in enumerate(nss_Res):
        if i != idx_MinDev:
            ns.Dispose()


    if fDevTol is None:
        return (ns_MinDev, dev_Min), None


    if dev_Min <= fDevTol:
        return (ns_MinDev, dev_Min), None

    # Fail.
    ns_MinDev.Dispose()

    return (None, dev_Min), None


def createSurface_NoDevLimit_NewDegs(rgSrf_In, iDeg, bDebug=False):
    """
    """

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return

    ns_fromIn = rgSrf_In.ToNurbsSurface()

    bSkipRebuildU = (iDeg == ns_fromIn.Degree(direction=0)) and ns_fromIn.SpanCount(0)
    bSkipRebuildV = (iDeg == ns_fromIn.Degree(direction=1)) and ns_fromIn.SpanCount(1)

    if bSkipRebuildU and bSkipRebuildV:
        return None, "Surface is already a Bezier patch."

    if iDeg == 3 and bSkipRebuildU or bSkipRebuildV:
        rgNurbsSrf_Res = ns_fromIn.RebuildOneDirection(
            direction=1 if bSkipRebuildU else 0,
            pointCount=4,
            loftType=rg.LoftType.Normal,
            refitTolerance=0.1*fDevTol)
        if rgNurbsSrf_Res is None:
            ns_fromIn.Dispose()
            return None, "RebuildOneDirection returned None."
    else:
        rgNurbsSrf_Res = ns_fromIn.Rebuild(
            uDegree=iDeg,
            vDegree=iDeg,
            uPointCount=iDeg+1,
            vPointCount=iDeg+1)

        if rgNurbsSrf_Res is None:
            ns_fromIn.Dispose()
            return None, "Rebuild returned None.  Check for internal knots near surface boundary."

    rc = xBrep_getDistancesBetween2.getDistancesBetweenBreps(
            ns_fromIn,
            rgNurbsSrf_Res,
            bIncludeEdges=True,
            rgMeshParams=None,
            bFineMesh=False,
            bCalcBrepIntersection=False,
            bDebug=bDebug
            )
        
    srf_dev = rc[1] if rc[0] else None
        
    ns_fromIn.Dispose()
    return (rgNurbsSrf_Res, srf_dev), None


def createSurface_DevLimit_TryDegs(rgSrf_In, fDevTol, iDegs_toTry=(1,2,3,5), bDebug=False):
    """
    Parameters:
        rgSrf_In
        fDevTol=,
        iDegsToTry=,
        bDebug=,
    Returns when passing surface found:
        tuple(
            tuple(rg.NurbsSurface, float(Maximum deviation from rgSrf_In)),
            None)
    Returns when passing surface NOT found:
        tuple(
            tuple(None, float(Maximum deviation from closest solution to rgSrf_In))
            None)
    Returns on other fails:
        tuple(None, str(Description of failure))
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return


    ns_fromIn = rgSrf_In.ToNurbsSurface()

    iDegs_toTry = sorted(iDegs_toTry)


    bSkipRebuildU = ns_fromIn.Degree(direction=0) in iDegs_toTry and ns_fromIn.SpanCount(0)
    bSkipRebuildV = ns_fromIn.Degree(direction=1) in iDegs_toTry and ns_fromIn.SpanCount(1)

    if bSkipRebuildU and bSkipRebuildV:
        return None, "Surface is already a Bezier patch."

    srf_devs = []


    if 1 not in iDegs_toTry:
        iDegs_toTry_U = iDegs_toTry[:]
        iDegs_toTry_V = iDegs_toTry[:]
    else:
        # Quick check for not linear.

        c_VDom_Start = ns_fromIn.IsoCurve(0, ns_fromIn.Domain(1).T0)
        c_VDom_End = ns_fromIn.IsoCurve(0, ns_fromIn.Domain(1).T1)
        #sc.doc.Objects.AddCurve(c); sc.doc.Views.Redraw(); return
    
        if c_VDom_Start.IsLinear(fDevTol) and c_VDom_End.IsLinear(fDevTol):
            iDegs_toTry_U = iDegs_toTry[:]
        else:
            if len(iDegs_toTry) > 1:
                # Next degree after 1.
                iDegs_toTry_U = iDegs_toTry[1:]
            else:
                return None, "Surface is not linear along U."

        c_UDom_Start = ns_fromIn.IsoCurve(1, ns_fromIn.Domain(0).T0)
        c_UDom_End = ns_fromIn.IsoCurve(1, ns_fromIn.Domain(0).T1)
        #map(sc.doc.Objects.AddCurve, (c_UDom_Start, c_UDom_End)); sc.doc.Views.Redraw(); return

        if c_UDom_Start.IsLinear(fDevTol) and c_UDom_End.IsLinear(fDevTol):
            iDegs_toTry_V = iDegs_toTry[:]
        else:
            if len(iDegs_toTry) > 1:
                # Next degree after 1.
                iDegs_toTry_V = iDegs_toTry[1:]
            else:
                return None, "Surface is not linear along V."


    # Try rebuilding at various allowed degrees.
    for iDeg_U in iDegs_toTry_U:
        for iDeg_V in iDegs_toTry_V:
            
            #print(iDegreeU, iDegreeV

            if bSkipRebuildU or bSkipRebuildV:
                # Surface.RebuildOneDirection only seems to output degree 3
                # in rebuilt direction.  Therefore, use pointCount of 4 for Bezier.
                pointCount = ns_fromIn.Degree(1 if bSkipRebuildU else 0)+1
                # For loftType, Developable, Loose, Normal, Straight, and Tight work.
                # Uniform doesn't.
                rgNurbsSrf_Res = ns_fromIn.RebuildOneDirection(
                    direction=1 if bSkipRebuildU else 0,
                    pointCount=4,
                    loftType=rg.LoftType.Normal,
                    refitTolerance=0.1*fDevTol)
                if rgNurbsSrf_Res is None:
                    continue # to next degree V.
            else:
                rgNurbsSrf_Res = ns_fromIn.Rebuild(
                    uDegree=iDeg_U,
                    vDegree=iDeg_V,
                    uPointCount=iDeg_U+1,
                    vPointCount=iDeg_V+1)

                if rgNurbsSrf_Res is None:
                    ns_fromIn.Dispose()
                    return None, "Rebuild returned None.  Check for internal knots near surface boundary."

            rc = xBrep_getDistancesBetween2.getDistancesBetweenBreps(
                    ns_fromIn,
                    rgNurbsSrf_Res,
                    bIncludeEdges=True,
                    rgMeshParams=None,
                    bFineMesh=False,
                    bCalcBrepIntersection=False,
                    bDebug=bDebug
                    )
        
            srf_dev = rc[1] if rc[0] else None
        
            if fDevTol is None:
                ns_fromIn.Dispose()
                return (rgNurbsSrf_Res, srf_dev), None
        
            if rc[0] and srf_dev <= fDevTol:
                # Success.
                ns_fromIn.Dispose()
                return (rgNurbsSrf_Res, srf_dev), None
        
            # Surface isn't within deviation tolerance.
            rgNurbsSrf_Res.Dispose()
            srf_devs.append(srf_dev)
            
            if bSkipRebuildU or bSkipRebuildV:
                # Don't bother looping because RebuildOneDirection only rebuilds
                # to degree 3.
                return(None, min(srf_devs) if srf_devs else None), None

    return (None, min(srf_devs) if srf_devs else None), None


def processFace(rgFace_In, fDevTol=None, iDegsToTry=None, bDebug=False):
    """
    Parameters:
        BrepFace
        fDevTol
        iDegsToTry
        bDebug
    Returns:
        (rg.Brep (1-face), float(deviation)), None
        (None, float(deviation needed)), sLog
        None, None
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return


    rgSrf_In = rgFace_In.UnderlyingSurface()


    rgNurbsSrf0 = rgSrf_In

    if iDegsToTry is None:
        rc = createSurface_SameDegs(
            rgNurbsSrf0,
            fDevTol=fDevTol,
            bDebug=bDebug
            )
    else:
        # New degrees.
        if fDevTol is None:
            rc = createSurface_NoDevLimit_NewDegs(
                rgNurbsSrf0,
                iDeg=min(iDegsToTry),
                bDebug=bDebug
                )
        else:
            rc = createSurface_DevLimit_TryDegs(
                rgNurbsSrf0,
                fDevTol=fDevTol,
                iDegs_toTry=iDegsToTry,
                bDebug=bDebug
                )
    
    if rc[0] is None or rc[0][0] is None: return rc
    
    rgNurbsSrf_Converted, srf_dev = rc[0]
    
    # Success in creating NurbsSurface.  Now, create correctly trimmed brep.
    
    rgBrep0_1Face = rgFace_In.DuplicateFace(duplicateMeshes=False)
    
    rgBrep0_1Face.Faces[0].RebuildEdges(0.1*sc.doc.ModelAbsoluteTolerance, True, True)
    #rgBrep0_1Face.Faces.ShrinkFaces()
    
    if rgBrep0_1Face.IsSurface:
        rgBrep0_1Face.Dispose()
        rgBrep1_1Face = rgNurbsSrf_Converted.ToBrep()
        rgNurbsSrf_Converted.Dispose()
        if not rgBrep1_1Face.IsValid:
            rgBrep1_1Face.Dispose()
            return None, "Invalid brep geometry after ToBrep."
        return (rgBrep1_1Face, srf_dev), None
    
    # Test areas before trimming.
    fArea_Trimmed = rgBrep0_1Face.GetArea()
    if fArea_Trimmed:
        rgBrep_Untrimmed = rgBrep0_1Face.Faces[0].UnderlyingSurface().ToBrep()
        fArea_Untrimmed = rgBrep_Untrimmed.GetArea()
        rgBrep_Untrimmed.Dispose()
        if fArea_Untrimmed:
            if abs(fArea_Trimmed - fArea_Untrimmed) <= sc.doc.ModelAbsoluteTolerance:
                rgBrep1_1Face = rgNurbsSrf_Converted.ToBrep()
                rgNurbsSrf_Converted.Dispose()
                if not rgBrep1_1Face.IsValid:
                    rgBrep1_1Face.Dispose()
                    return None, "Invalid brep geometry after ToBrep."
                return (rgBrep1_1Face, srf_dev), None

    rgBrep1_1Face = xBrepFace.retrimFace(
            rgBrep0_1Face.Faces[0],
            rgSrf_Replacement=rgNurbsSrf_Converted,
            fSplitTol=1.0*sc.doc.ModelAbsoluteTolerance if fDevTol is None else fDevTol,
            bDebug=bDebug
    )
    rgNurbsSrf_Converted.Dispose()
    rgBrep0_1Face.Dispose()

    if rgBrep1_1Face is None:
        return None, "xBrepFace.createMonofaceBrep returned None."

    if not rgBrep1_1Face.IsValid:
        rgBrep1_1Face.Dispose()
        return None, "An invalid brep was skipped."

    return (rgBrep1_1Face, srf_dev), None


def processBrep(rgBrep_In, idxs_rgFaces, fDevTol=None, iDegsToTry=None, **kwargs):
    """
    Returns:
        (rgBreps_1F_Mod, idxs_rgFaces_Rebuilt, srf_devs), sLogs
        None
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bSkipConical = getOpt('bSkipConical')
    bShrinkFirst = getOpt('bShrinkFirst')
    bDebug = getOpt('bDebug')


    rgBreps_1F_Mod = []
    idxs_rgFaces_Rebuilt = []
    srf_devs = []
    srf_devs_Needed = []
    
    rgB_WIP = rgBrep_In.DuplicateBrep()
    
    if bShrinkFirst:
        rgB_WIP.Faces.ShrinkFaces()
    
    sLogs = []
    
    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt
    
    idxs_AtTenths = [int(round(0.1*i*len(idxs_rgFaces),0)) for i in range(10)]
    
    for iF, idx_rgFace in enumerate(idxs_rgFaces):
        if sc.escape_test(False):
            print("Searching interrupted by user.")
            return
        
        if iF in idxs_AtTenths:
            s = sCmdPrompt0 + ", {:d}% of {} faces in current brep ...".format(
                int(100.0 * (iF+1) / len(idxs_rgFaces)), len(idxs_rgFaces))
            
            if bDebug:
                print(s)
            else:
                Rhino.RhinoApp.SetCommandPrompt(s)
        
        rgFace_In = rgB_WIP.Faces[idx_rgFace]


        if bSkipConical:
            rc = doesSurfaceContainConicSection(rgFace_In.UnderlyingSurface())
            if rc:
                sLogs.append("Skipped " + rc)
                continue

            #if isinstance(rgFace_In.UnderlyingSurface(), rg.RevSurface):
            #    sLogs.append("Skipped RevSurface.")
            #    continue
            #elif isinstance(rgFace_In.UnderlyingSurface(), rg.NurbsSurface):
            #    if rgFace_In.UnderlyingSurface().IsRational:
            #        sLogs.append("Skipped rational NurbsSurface.")
            #        continue


        rc = processFace(
            rgFace_In,
            fDevTol=fDevTol,
            iDegsToTry=iDegsToTry,
            bDebug=bDebug
            )


        if rc[0] is None:
            sLogs.append(rc[1])
            continue
        
        rgBrep_1F_Converted, srf_dev = rc[0]
        if rgBrep_1F_Converted is None:
            if srf_dev is not None:
                srf_devs_Needed.append(srf_dev)
            continue

        rgBreps_1F_Mod.append(rgBrep_1F_Converted)
        idxs_rgFaces_Rebuilt.append(idx_rgFace)
        srf_devs.append(srf_dev)
    
    rgB_WIP.Dispose()
    
    if srf_devs_Needed:
        s  = "Need tolerances of {:.2e} through {:.2e}".format(
                min(srf_devs_Needed), max(srf_devs_Needed))
        s += " to convert remaining convertible surfaces."
        sLogs.append(s)
    
    return (rgBreps_1F_Mod, idxs_rgFaces_Rebuilt, srf_devs), sLogs


def processBrepObjects(rhBreps, idx_Faces=None, fDevTol=None, iDegsToTry=None, **kwargs):
    """
    Parameters:
        rhBreps: Objrefs of brep with face components, GUIDs, rd.Breps.
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bSkipConical = getOpt('bSkipConical')
    bShrinkFirst = getOpt('bShrinkFirst')
    bReplace = getOpt('bReplace')
    bExtract = getOpt('bExtract')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getRhinoObject(rhObj):
        """
        'Deleted objects cannot be found by id.'
        (https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_DocObjects_Tables_ObjectTable_FindId.htm)
        """
        rdObj = None
        if isinstance(rhObj, rd.RhinoObject):
            rdObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        return rdObj


    def getBrepObject(rhObj):
        rdObj = getRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def getSortedBrepIdsAndFaces(objrefs):
        """
        Parameters:
            list(objrefs)
        Returns:
            list(Brep GUIDs)
            list(lists(integers of Face indices) per brep)
        """
        
        gBreps0 = []
        idxs_Faces_perBrep = []
    
        for o in objrefs:
            gBrep0 = o.ObjectId
            rdBrep_In = o.Object()
            rgBrep_In = o.Brep()
        
            if not rgBrep_In.IsValid:
                print("Brep {} is invalid.  Fix first.".format(gBrep0))
                rgBrep_In.Dispose()
                continue
        
            idx_CompIdx = o.GeometryComponentIndex.Index
            if idx_CompIdx == -1:
                if gBrep0 in gBreps0:
                    idxs_Faces_perBrep[gBreps0.index(gBrep0)] = range(rgBrep_In.Faces.Count)
                else:
                    gBreps0.append(gBrep0)
                    idxs_Faces_perBrep.append(range(rgBrep_In.Faces.Count))
            else:
                rgFace_Brep0 = o.Face()
                if gBrep0 in gBreps0:
                    if rgFace_Brep0 in idxs_Faces_perBrep[gBreps0.index(gBrep0)]:
                        continue
                    else:
                        idxs_Faces_perBrep[gBreps0.index(gBrep0)].append(rgFace_Brep0.FaceIndex)
                else:
                    gBreps0.append(gBrep0)
                    idxs_Faces_perBrep.append([rgFace_Brep0.FaceIndex])

        return gBreps0, idxs_Faces_perBrep


    def isKnotVectorUniform(knots):
        return (
            (knots.KnotStyle == rg.KnotStyle.Uniform) or
            (knots.KnotStyle == rg.KnotStyle.QuasiUniform) or
            (
                (knots.KnotStyle == rg.KnotStyle.PiecewiseBezier) and
                knots.Count == knots.KnotMultiplicity(0) * 2)
            )


    def getNurbsSurfaceChangeDescription(rgNurbsSrf1, rgNurbsSrf2):
        s  = "  Prop:I->O"
        s += "  {}:({}->{})x({}->{})".format(
                "Deg",
                rgNurbsSrf1.OrderU-1,
                rgNurbsSrf2.OrderU-1,
                rgNurbsSrf1.OrderV-1,
                rgNurbsSrf2.OrderV-1,
        )
        s += "  {}:({}->{})x({}->{})".format(
                "PtCt",
                rgNurbsSrf1.Points.CountU,
                rgNurbsSrf2.Points.CountU,
                rgNurbsSrf1.Points.CountV,
                rgNurbsSrf2.Points.CountV,
        )
        if Rhino.RhinoApp.ExeVersion >= 7:
            s += "  {}:({}->{})x({}->{})".format(
                    "IsUniform",
                    str(isKnotVectorUniform(rgNurbsSrf1.KnotsU))[0],
                    str(isKnotVectorUniform(rgNurbsSrf2.KnotsU))[0],
                    str(isKnotVectorUniform(rgNurbsSrf1.KnotsV))[0],
                    str(isKnotVectorUniform(rgNurbsSrf2.KnotsV))[0],
            )
        s += "  {}:{}->{}".format(
                "IsRational",
                str(rgNurbsSrf1.IsRational)[0],
                str(rgNurbsSrf2.IsRational)[0],
        )
        s += "  {}:({}->{})x({}->{})".format(
                "IsClosed",
                str(rgNurbsSrf1.IsClosed(0))[0],
                str(rgNurbsSrf2.IsClosed(0))[0],
                str(rgNurbsSrf1.IsClosed(1))[0],
                str(rgNurbsSrf2.IsClosed(1))[0],
        )
        if (    rgNurbsSrf1.IsClosed(0) or rgNurbsSrf1.IsClosed(1) or
                rgNurbsSrf2.IsClosed(0) or rgNurbsSrf2.IsClosed(1)
        ):
            s += "  {}:{}->{}".format(
                    "IsPeriodic",
                    str(rgNurbsSrf1.IsPeriodic)[0],
                    str(rgNurbsSrf2.IsPeriodic)[0],
            )
        return s


    gBreps0, idxs_rgFace_PerBrep = getSortedBrepIdsAndFaces(rhBreps)
    if not gBreps0: return

    gBs1_perB0 = []
    srf_devs_All = []
    sLogs_All = []
    
    len_gBreps0 = len(gBreps0)
    idxs_AtTenths = [int(round(0.1*i*len_gBreps0,0)) for i in range(10)]
    
    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt
    
    if len(rhBreps) == 1:
        s = sCmdPrompt0 + "Brep"
        Rhino.RhinoApp.SetCommandPrompt(s)
    
    for iB, (gBrep0, idxFaces) in enumerate(zip(gBreps0, idxs_rgFace_PerBrep)):
        rdBrep_In = getBrepObject(gBrep0)
        rgBrep_In = rdBrep_In.Geometry

        if len(rhBreps) > 1:
            if iB in idxs_AtTenths:
                s = sCmdPrompt0 + "  At {:d}% of {} breps".format(
                    int(100.0 * (iB+1) / len_gBreps0), len_gBreps0)
        
            if bDebug:
                print(s)
            else:
                Rhino.RhinoApp.SetCommandPrompt(s)
        
        rc = processBrep(
                rgBrep_In,
                idxFaces,
                fDevTol=fDevTol,
                iDegsToTry=iDegsToTry,
                bSkipConical=bSkipConical,
                bShrinkFirst=bShrinkFirst,
                bDebug=bDebug)
        if rc is None: return
        sLogs_All.extend(rc[1])
        if not rc[0]: continue
        rgBreps_1F_Mod_thisBrep, idxsFs_Mod, srf_devs = rc[0]
        
        if not rgBreps_1F_Mod_thisBrep:
            rgBrep_In.Dispose()
            continue
        
        srf_devs_All.extend(srf_devs)
        
        if not bReplace:
            gBrep1_thisBrep = []
            for rgBrep_1F_New, idxF in zip(rgBreps_1F_Mod_thisBrep, idxsFs_Mod):
                gBrep1 = sc.doc.Objects.AddBrep(rgBrep_1F_New)
                if gBrep1 != Guid.Empty:
                    gBrep1_thisBrep.append(gBrep1)
            gBs1_perB0.append(gBrep1_thisBrep)
        else:
            # bReplace==True.
            if bExtract:
                rc = xBrepObject.replaceFaces(
                        rdBrep_In,
                        idxsFs_Mod,
                        rgBreps_1F_Mod_thisBrep,
                        bExtract=True,
                        fTolerance_Join=max(srf_devs))
                if rc:
                    gBreps1_NewFaces_thisBrep, gBreps_RemainingBrep = rc
                    gBs1_perB0.append(gBreps1_NewFaces_thisBrep)
            else:
                fTols_Edges = [edge.Tolerance for edge in rgBrep_In.Edges]
                fTolerance_Join = max((
                        2.0*fDevTol if fDevTol is not None else 0.0,
                        1.1*max(fTols_Edges),
                        sc.doc.ModelAbsoluteTolerance))
                rc = xBrepObject.replaceFaces(
                        rdBrep_In,
                        idxsFs_Mod,
                        rgBreps_1F_Mod_thisBrep,
                        bExtract=False,
                        fTolerance_Join=fTolerance_Join)
                if rc:
                    gBreps_withReplacedFaces_thisBrep = rc
                    gBs1_perB0.append(gBreps_withReplacedFaces_thisBrep)

        if bDebug or bEcho and len(gBreps0) == 1 and len(idxsFs_Mod)==1:
            rgBrep_1F_Mod = rgBreps_1F_Mod_thisBrep[0]
            s  = getNurbsSurfaceChangeDescription(
                    rgBrep_In.Faces[idxsFs_Mod[0]].UnderlyingSurface(),
                    rgBrep_1F_Mod.Surfaces[0])
            s += "  Deviation: {0:.2e}".format(srf_devs[0])
            print(s)

        for brep in rgBreps_1F_Mod_thisBrep: brep.Dispose()
        rgBrep_In.Dispose()

        if (
                any(g for gs in gBs1_perB0 for g in gs) and
                len(gBreps0)==1
                and (bDebug or bEcho)
        ):
            sc.doc.Objects.UnselectAll()
            if bReplace:
                if bExtract:
                    print("{} face(s) extracted from brep and replaced.".format(
                        len(gBreps1_NewFaces_thisBrep)))
                else:
                    print("{} face(s) replaced in brep.".format(len(idxsFs_Mod)))
            else:
                print("{} monoface brep(s) added.".format(
                    sum(len(bs) for bs in rgBreps_1F_Mod_thisBrep)))
    
    return (gBs1_perB0, srf_devs_All), sLogs_All


def main():
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script doesn't work in Rhino versions previous to 6.")
        return
    
    rc = getInput()
    if rc is None: return

    (
        objrefs0,
        fDevTol,
        iDegsToTry,
        bSkipConical,
        bShrinkFirst,
        bReplace,
        bExtract,
        bEcho,
        bDebug,
        ) = rc


    if not bDebug:
        sc.doc.Views.RedrawEnabled = False
    
    Rhino.RhinoApp.SetCommandPrompt("Working ...")
    
    rc = processBrepObjects(
        rhBreps=objrefs0,
        idxs_Faces=None,
        fDevTol=fDevTol,
        iDegsToTry=iDegsToTry,
        bSkipConical=bSkipConical,
        bShrinkFirst=bShrinkFirst,
        bReplace=bReplace,
        bExtract=bExtract,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    (gBreps1, srf_devs_All), sLogs = rc

    iCt_NeedLargerTol = 0
    
    for sLog in set(sLogs):
        if not bDebug and len(sLogs) > 1 and "Need tolerances" in sLog:
            iCt_NeedLargerTol += 1
        else:
            print("[{}] {}".format(sLogs.count(sLog), sLog))
    
    if iCt_NeedLargerTol:
        print("[{}] Need larger tolerance to convert.".format(iCt_NeedLargerTol))

    if not gBreps1:
        print("Nothing was converted.")
    elif gBreps1 and len(gBreps1)>1 and (bDebug or bEcho):
        print("Maximum deviations: [{:.2e},{:.2e}]".format(
            min(srf_devs_All), max(srf_devs_All)))
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()