"""
"""

"""
181120: Created.
190528: Added Opts and support for edge input.
190529: Added bFitResult.
190712-14: Refactored.  Polished the UX.
190724: Bug fixes in UI.
190924-200121, 220328: Import-related updates.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Enum
from System import Guid

import spb_Crv_fitRebuild


sOpts = (
        'bRebuildPath',
        'fDraft_Deg',
        'bCPlane',
        'fDistance',
        'bExplodePolyCrvs',
        'bSplitPathsAtKnots',
        'iCornerType',
        'fDistTol',
        'fAngleTol_Deg',
        'bFitResult',
        'bEcho',
        'bDebug',
)


class Opts():
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    
    for key in sOpts:
        keys.append(key)
        names[key] = key[1:] # Overwrite as wanted in the following.
    
    key = 'bRebuildPath'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fDraft_Deg'
    values[key] = 45.0
    names[key] = 'DraftAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCPlane'
    values[key] = True
    names[key] = 'PosDirPerZAxisOf'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fDistance'
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bExplodePolyCrvs'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSplitPathsAtKnots'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iCornerType'
    values[key] = 2
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fDistTol'
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'fAngleTol_Deg'
    values[key] = sc.doc.ModelAngleToleranceDegrees
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bFitResult'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    
    @classmethod
    def setValues(cls):
        for key in sOpts:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                pass
    
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(bFirstGetObjects=False):
    """
    Get curves or edges with optional input.
    
    Returns
        None to cancel.
        False to end create/GetObject cycle.
        tuple of ObjRefs and options to continue create/GetObject cycle.
    """
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select curve or edges")
    
    go.GeometryFilter = rd.ObjectType.Curve
    
    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    go.AddOptionToggle(Opts.names['bRebuildPath'], Opts.riOpts['bRebuildPath'])
    go.AddOptionDouble(Opts.names['fDraft_Deg'], Opts.riOpts['fDraft_Deg'])
    idxs_Opts['FlipAngle'] = go.AddOption('FlipAngle')
    go.AddOptionToggle(Opts.names['bCPlane'], Opts.riOpts['bCPlane'])
    go.AddOptionDouble(Opts.names['fDistance'], Opts.riOpts['fDistance'])
    idxs_Opts['FlipDir'] = go.AddOption('FlipDir')
    go.AddOptionToggle(Opts.names['bExplodePolyCrvs'], Opts.riOpts['bExplodePolyCrvs'])
    go.AddOptionToggle(Opts.names['bSplitPathsAtKnots'], Opts.riOpts['bSplitPathsAtKnots'])
    idxs_Opts['iCornerType'] = go.AddOptionList(
            englishOptionName=Opts.names['iCornerType'],
            listValues=Enum.GetNames(rg.ExtrudeCornerType),
            listCurrentIndex=Opts.values['iCornerType'])
    go.AddOptionDouble(Opts.names['fDistTol'], Opts.riOpts['fDistTol'])
    go.AddOptionDouble(Opts.names['fAngleTol_Deg'], Opts.riOpts['fAngleTol_Deg'])
    go.AddOptionToggle(Opts.names['bFitResult'], Opts.riOpts['bFitResult'])
    go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
    go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])

    go.AlreadySelectedObjectSelect = True # So objects can be reselected after being unselected in same go.GetMultiple.
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

    if go.ObjectsWerePreselected:
        if bFirstGetObjects:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in sOpts])

        iCt_Crvs_PreSelctd = go.ObjectCount
        objrefs = go.Objects()
        gCrvs_PreSelctd = []
        gBreps_PreSelctd = []
        idxs_Edges_PerBrep = []
        for objref in objrefs:
            rgCrv = objref.Curve()
            if isinstance(rgCrv, rg.BrepEdge):
                gBreps_PreSelctd.append(objref.ObjectId)
                idxs_Edges_PerBrep.append(rgCrv.EdgeIndex)
            else:
                # Curves other than BrepEdges.
                gCrvs_PreSelctd.append(objref.ObjectId)
            rgCrv.Dispose()

        go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
    else:
        objrefs = None
        iCt_Crvs_PreSelctd = 0
    
    if res == ri.GetResult.Object:
        if iCt_Crvs_PreSelctd == go.ObjectCount:
            def wasThereChangeInObjectsSelected(gCrvs_PreSelctd, gBreps_PreSelctd, idxs_Edges_PerBrep):
                objrefs = go.Objects()
                for objref in objrefs:
                    rgCrv = objref.Curve()
                    if isinstance(rgCrv, rg.BrepEdge):
                        rgCrv.Dispose()
                        if not objref.ObjectId in gBreps_PreSelctd:
                            return True
                        zipped = zip(gBreps_PreSelctd, idxs_Edges_PerBrep)
                        for gBrep_PreSelctd, idx_Edge_PerBrep in zipped:
                            if gBrep_PreSelctd == objref.ObjectId:
                                if idx_Edge_PerBrep == idx_Edge_PerBrep:
                                    break
                        else:
                            # Edge not found.
                            return True
                    else:
                        # Curves other than BrepEdges.
                        rgCrv.Dispose()
                        if not objref.ObjectId in gCrvs_PreSelctd:
                            return True
                return False
            if not wasThereChangeInObjectsSelected(
                    gCrvs_PreSelctd=gCrvs_PreSelctd,
                    gBreps_PreSelctd=gBreps_PreSelctd,
                    idxs_Edges_PerBrep=idxs_Edges_PerBrep
            ):
                go.Dispose()
                return tuple([False] + [Opts.values[key] for key in sOpts])
        objrefs = go.Objects()
        go.Dispose()
        return tuple([objrefs] + [Opts.values[key] for key in sOpts])
    elif res == ri.GetResult.Cancel:
        go.Dispose()
        return
    elif res == ri.GetResult.Number:
        Opts.riOpts['fDraft_Deg'].CurrentValue = go.Number()
    elif go.OptionIndex() == idxs_Opts['FlipAngle']:
        Opts.riOpts['fDraft_Deg'].CurrentValue = -Opts.riOpts['fDraft_Deg'].CurrentValue
    elif go.OptionIndex() == idxs_Opts['FlipDir']:
        Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
    elif go.OptionIndex() == idxs_Opts['iCornerType']:
        Opts.values['iCornerType'] = go.Option().CurrentListOptionIndex
    elif Opts.riOpts['fDistTol'].CurrentValue < 0.0:
        Opts.riOpts['fDistTol'].CurrentValue = Opts.riOpts['fDistTol'].InitialValue
    elif Opts.riOpts['fAngleTol_Deg'].CurrentValue < 0.0:
        Opts.riOpts['fAngleTol_Deg'].CurrentValue = Opts.riOpts['fAngleTol_Deg'].InitialValue

    Opts.setValues()
    Opts.saveSticky()
    go.ClearCommandOptions()
    
    return tuple([objrefs] + [Opts.values[key] for key in sOpts])


def prepareCurves_OLD(rgCrvs0, bExplodePolyCrvs=True, bSplitPathsAtKnots=False):
    """
    Prepare curves, including splitting per options.

    Returns: list of new curves
    """
    
    rc = rg.Curve.JoinCurves(
            rgCrvs0,
            joinTolerance=sc.doc.ModelAbsoluteTolerance)
    if not rc: return
    rgCrvs_Joined = rc

    rgCrvs_Final = []
    for rgCrv_Joined in rgCrvs_Joined:
        if isinstance(rgCrv_Joined, rg.PolyCurve):
            rgCrv_Joined.RemoveNesting()

        if bExplodePolyCrvs:
            if isinstance(rgCrv_Joined, rg.PolyCurve):
                rc = rgCrv_Joined.Explode()
                if rc: rgCrvs_SplitPoly = rc
            else:
                rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()]
        else:
            rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()]

        if bSplitPathsAtKnots:
            for c in rgCrvs_SplitPoly:
                ts_SpanBoundaries = [c.Domain.T0]
                for iSpan in xrange(c.SpanCount):
                    ts_SpanBoundaries.append(c.SpanDomain(iSpan).T1)
                rc = c.Split(ts_SpanBoundaries)
                if rc: rgCrvs_Final.extend(rc)
                c.Dispose()
        else:
            rgCrvs_Final.extend(rgCrvs_SplitPoly)
        rgCrv_Joined.Dispose()
    
    return rgCrvs_Final


def rebuildCurves(rgCrvs):
    """
    Returns bool for if any rebuilds occur.  If so, the input curves are updated.
    """
    cs1 = []
    bRebuilt = False

    for c0 in rgCrvs:
        rc = spb_Crv_fitRebuild.rebuildCurve(
                c0,
                0.1*sc.doc.ModelAbsoluteTolerance,
                iDegree=3,
                bPreserveEndTans=True,
                bFurtherTranslateCps=False,
                iMinCpCt=None,
                iMaxCpCt=40,
                bDebug=False,
        )
        if rc[0]:
            cs1.append(rc[0])
            bRebuilt = True
        else:
            cs1.append(c0.Duplicate())

    if bRebuilt:
        for c in rgCrvs: c.Dispose()
        rgCrvs[:] = cs1[:]
    else:
        for c in cs1: c.Dispose()

    return bRebuilt


def joinCurves(rgCrvs):
    """
    Returns bool for if any rebuilds occur.  If so, the input curves are updated.
    """

    bJoined = False

    cs1 = rg.Curve.JoinCurves(
            rgCrvs,
            joinTolerance=sc.doc.ModelAbsoluteTolerance)
    if cs1: bJoined = len(cs1) < len(rgCrvs)

    if bJoined:
        for c in rgCrvs: c.Dispose()
        rgCrvs[:] = cs1[:]
    else:
        for c in cs1: c.Dispose()

    return bJoined


def explodePolyCurves(rgCrvs):
    """
    Returns bool for if any explosions occur.  If so, the input curves are updated.
    """
    cs1 = []
    bExplosion = False

    for c0 in rgCrvs:
        if isinstance(c0, rg.PolyCurve):
            c0.RemoveNesting()
            rc = c0.Explode()
            if rc:
                cs1.extend(rc)
                bExplosion = True
                continue
        cs1.append(c0.Duplicate())

    if bExplosion:
        for c in rgCrvs: c.Dispose()
        rgCrvs[:] = cs1[:]
    else:
        for c in cs1: c.Dispose()

    return bExplosion


def splitCurvesAtSpans(rgCrvs):
    """
    Returns bool for if any splits occur.  If so, the input curves are updated.
    """
    cs1 = []
    bSplit = False

    for c0 in rgCrvs:
        ts_SpanBoundaries = []
        for iSpan in xrange(c0.SpanCount if c0.IsPeriodic else c0.SpanCount-1):
            ts_SpanBoundaries.append(c0.SpanDomain(iSpan).T1)
        print len(ts_SpanBoundaries)
        rc = c0.Split(ts_SpanBoundaries)
        if rc:
            cs1.extend(rc)
            bSplit = True
        else:
            cs1.append(c0.Duplicate())

    if bSplit:
        for c in rgCrvs: c.Dispose()
        rgCrvs[:] = cs1[:]
    else:
        for c in cs1: c.Dispose()

    return bSplit


def createBreps(rgCrvs_ToExtrude, **opts):
    """
    """

    bRebuildPath = opts['bRebuildPath'] if 'bRebuildPath' in opts else Opts.values['bRebuildPath']
    fDraft_Deg = opts['fDraft_Deg'] if 'fDraft_Deg' in opts else Opts.values['fDraft_Deg']
    bCPlane = opts['bCPlane'] if 'bCPlane' in opts else Opts.values['bCPlane']
    fDistance = opts['fDistance'] if 'fDistance' in opts else Opts.values['fDistance']
    bExplodePolyCrvs = opts['bExplodePolyCrvs'] if 'bExplodePolyCrvs' in opts else Opts.values['bExplodePolyCrvs']
    bSplitPathsAtKnots = opts['bSplitPathsAtKnots'] if 'bSplitPathsAtKnots' in opts else Opts.values['bSplitPathsAtKnots']
    iCornerType = opts['iCornerType'] if 'iCornerType' in opts else Opts.values['iCornerType']
    fDistTol = opts['fDistTol'] if 'fDistTol' in opts else Opts.values['fDistTol']
    bFitResult = opts['bFitResult'] if 'bFitResult' in opts else Opts.values['bFitResult']
    fAngleTol_Deg = opts['fAngleTol_Deg'] if 'fAngleTol_Deg' in opts else Opts.values['fAngleTol_Deg']
    bEcho = opts['bEcho'] if 'bEcho' in opts else Opts.values['bEcho']
    bDebug = opts['bDebug'] if 'bDebug' in opts else Opts.values['bDebug']


    if bCPlane:
        view = sc.doc.Views.ActiveView
        plane = view.ActiveViewport.ConstructionPlane()
        direction = plane.Normal
    else:
        direction = rg.Vector3d.ZAxis

    basePoint = rg.Point3d.Origin


    def prepareCurves_NestedFunc(rgCrvs0):
        """
        Returns list of new curves.
        """

        cs0 = [c.Duplicate() for c in rgCrvs0]

        for c in cs0:
            if isinstance(c, rg.PolyCurve): c.RemoveNesting()

        if bRebuildPath: rebuildCurves(cs0)

        if joinCurves(cs0) and bRebuildPath: rebuildCurves(cs0)

        if bExplodePolyCrvs:
            if explodePolyCurves(cs0) and bRebuildPath: rebuildCurves(cs0)

        if bSplitPathsAtKnots:
            if splitCurvesAtSpans(cs0) and bRebuildPath: rebuildCurves(cs0)

        return cs0
    rc = prepareCurves_NestedFunc(rgCrvs_ToExtrude)

    #rc = prepareCurves_OLD(
    #        rgCrvs0=rgCrvs_ToExtrude,
    #        bExplodePolyCrvs=bExplodePolyCrvs,
    #        bSplitPathsAtKnots=bSplitPathsAtKnots)
    #if not rc: return
    for c in rgCrvs_ToExtrude: c.Dispose()
    rgCrvs_ToExtrude = rc

    rgBreps2_Extrds_ToJoin = []
    for rgCrv_ToExtrude in rgCrvs_ToExtrude:
        rgBreps1_Extrds_Raw = rg.Brep.CreateFromTaperedExtrude(
                curveToExtrude=rgCrv_ToExtrude,
                distance=fDistance,
                direction=direction,
                basePoint=basePoint,
                draftAngleRadians=Rhino.RhinoMath.ToRadians(fDraft_Deg),
                cornerType=Enum.ToObject(rg.ExtrudeCornerType, iCornerType),
                tolerance=fDistTol,
                angleToleranceRadians=Rhino.RhinoMath.ToRadians(fAngleTol_Deg))
        if not rgBreps1_Extrds_Raw:
            print "An Extrude could not be created."
            continue

        if not bFitResult:
            rgBreps2_Extrds_ToJoin.extend(rgBreps1_Extrds_Raw)
        else:
            for rgBrep1_Extrd_Raw in rgBreps1_Extrds_Raw:
                for rgSrf1 in rgBrep1_Extrd_Raw.Surfaces:
                    rgNurbsSrf2 = rgSrf1.ToNurbsSurface()
                    rgNurbsSrf3 = rgSrf1.Fit(
                            uDegree=rgNurbsSrf2.Degree(0),
                            vDegree=rgNurbsSrf2.Degree(1),
                            fitTolerance=fDistTol)
                    if rgNurbsSrf3 is None:
                        print "A surface could not be Fit within {}.".format(fDistTol)
                        rgBreps2_Extrds_ToJoin.append(rgSrf1.ToBrep())
                    elif rgNurbsSrf3.EpsilonEquals(rgNurbsSrf2, epsilon=1e-9):
                        rgBreps2_Extrds_ToJoin.append(rgSrf1.ToBrep())
                        rgNurbsSrf3.Dispose()
                        print "A surface Fit didn't produce a new surface."
                    else:
                        rgBreps2_Extrds_ToJoin.append(rgNurbsSrf3.ToBrep())
                        rgNurbsSrf3.Dispose()
                    rgSrf1.Dispose()
                    rgNurbsSrf2.Dispose()
                rgBrep1_Extrd_Raw.Dispose()
    
    rgBreps_Extrds_Joined = rg.Brep.JoinBreps(
            rgBreps2_Extrds_ToJoin,
            tolerance=0.1*sc.doc.ModelAbsoluteTolerance) # A tighter tolerance can help identify problem areas.


    return rgBreps_Extrds_Joined


def createBrepObjects(rhObjects_CurveOrEdge, **opts):
    """
    """

    for key in Opts.keys:
        if not key in opts:
            opts[key] = Opts.values[key]
        if Opts.values['bDebug']: print key, opts[key]


    rgCrvs_ToExtrude = []
    for o in rhObjects_CurveOrEdge:
        if isinstance(o, rd.ObjRef):
            rgCrv_ToExtrude = o.Curve()
        else:
            if isinstance(o, Guid):
                rdCrv_ToExtrude = sc.doc.Objects.FindId(o) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(o)
            else:
                rdCrv_ToExtrude = o
            rgCrv_ToExtrude = rdCrv_ToExtrude.CurveGeometry

        if rgCrv_ToExtrude:
            if rgCrv_ToExtrude.IsValid:
                rgCrvs_ToExtrude.append(rgCrv_ToExtrude)
            else:
                print "Geometry of {} not valid!".format(o)
                continue
        else:
            print "Geometry of {} not found!".format(o)
            continue
    
    rc = createBreps(
            rgCrvs_ToExtrude=rgCrvs_ToExtrude,
            bRebuildPath=opts['bRebuildPath'],
            fDistance=opts['fDistance'],
            bCPlane=opts['bCPlane'],
            fDraft_Deg=opts['fDraft_Deg'],
            iCornerType=opts['iCornerType'],
            fDistTol=opts['fDistTol'],
            fAngleTol_Deg=opts['fAngleTol_Deg'],
            bEcho=opts['bEcho'],
            bDebug=opts['bDebug'])
    if rc is None:
        print "Breps not created!"
        for c in rgCrvs_ToExtrude: c.Dispose()
        sc.doc.Views.Redraw()
        return []

    return rc


def main():
    """
    """

    rc = getInput(bFirstGetObjects=True)
    if rc is None: return
    rhObjects_CurveOrEdge = rc[0]
    for key, value in zip(sOpts, rc[1:]):
        exec("{} = {}".format(key, value))

    gBreps1 = []
    
    while True:
        if rhObjects_CurveOrEdge is not None:
            rc = createBrepObjects(rhObjects_CurveOrEdge=rhObjects_CurveOrEdge)
            if rc is None: break
            elif rc:
                rgBreps1 = rc

                gBreps1 = []
                iCt_Faces_All = 0
                for rgBrep1 in rgBreps1:
                    gBrep1 = sc.doc.Objects.AddBrep(rgBrep1)
                    if gBrep1 != Guid.Empty:
                        gBreps1.append(gBrep1)
                        iCt_Faces_All += rgBrep1.Faces.Count
                    rgBrep1.Dispose()
                sc.doc.Views.Redraw()
            
                print "{} brep(s) with {} total face(s) were created.".format(
                        len(gBreps1), iCt_Faces_All)

        rc = getInput(bFirstGetObjects=False)
        if rc is None:
            for gBrep1 in gBreps1:
                sc.doc.Objects.Delete(gBrep1, True)
            return
        elif rc[0] is False:
            # Means there was no change in curve(s) selected.
            break
        else:
            rhObjects_CurveOrEdge = rc[0]
        
        for key, value in zip(sOpts, rc[1:]):
            exec("{} = {}".format(key, value))

        for gBrep1 in gBreps1:
            sc.doc.Objects.Delete(gBrep1, True)


if __name__ == '__main__': main()