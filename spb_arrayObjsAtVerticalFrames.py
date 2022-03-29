"""
"""

"""
190530-0531: Created.
190602: Replaced Sweep1 creation with Loft.
190604: Moved some code into its own script.  Simplified options.
190712: Fixed frame orientation when CPlane is not parallel to CPlane World Top.
190713-14: Refactored.  Polished the UX.
190727, 0903: Bug fix for closed curves.
190903: Added bRebuildPath.
190903-04: WIP - Start of creating smooth taper for variable taper angles.
191020: Import-related update.
191126: Bug fix.
200115-200121, 220328: Import-related updates.

TODO: Finish creating smooth taper for variable taper angles.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum
from System import Guid
from System.Drawing import Color

import spb_Crv_fitRebuild


sOpts = (
        'bRebuildPath',
        'fTaper_Start_Deg',
        'bVariableTaper',
        'fTaper_End_Deg',
        'bTaperChangePerCrvParam',
        'bCPlane',
        'pt_Base',
        'bAtGrevilles',
        'bAtKnots',
        'bAtEqualDivisions',
        'iDivisionCt',
        'bSplitPolyCrvToSegs',
        'bSplitPathsAtKnots',
        'bAddArrayedObjects',
        'bLoftCrvs',
        'iLoftType',
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
    
    key = 'fTaper_Start_Deg'
    values[key] = 0.0
    names[key] = 'TaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bVariableTaper'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTaper_End_Deg'
    values[key] = 45.0
    names[key] = 'EndTaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bTaperChangePerCrvParam'
    values[key] = False
    names[key] = 'TaperChangePerCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Length', 'Param')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCPlane'
    values[key] = False
    names[key] = 'PlanView'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'pt_Base'
    names[key] = 'BasePoint'
    values[key] = rg.Point3d(0.0,0.0,0.0)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAtGrevilles'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAtKnots'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAtEqualDivisions'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iDivisionCt'
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(
            initialValue=values[key],
            setLowerLimit=True,
            limit=1)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSplitPolyCrvToSegs'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSplitPathsAtKnots'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddArrayedObjects'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bLoftCrvs'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iLoftType'
    values[key] = 0
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
    
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput(sCmdPromp, rdObjectType):
    """
    Get objects with optional input.
    
    Returns
        None to cancel.
        False to end create/GetObject cycle.
        tuple of ObjRefs and options to continue create/GetObject cycle.
    """
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt(sCmdPromp)

    go.GeometryFilter = rdObjectType

    go.AcceptNumber(True, acceptZero=False)
    
    idxs_Opts = {}

    go.AddOptionToggle(Opts.names['bRebuildPath'], Opts.riOpts['bRebuildPath'])
    go.AddOptionDouble(Opts.names['fTaper_Start_Deg'], Opts.riOpts['fTaper_Start_Deg'])
    go.AddOptionToggle(Opts.names['bVariableTaper'], Opts.riOpts['bVariableTaper'])
    if Opts.values['bVariableTaper']:
        go.AddOptionDouble(Opts.names['fTaper_End_Deg'], Opts.riOpts['fTaper_End_Deg'])
        idxs_Opts['SwapAngles'] = go.AddOption('SwapAngles')
        go.AddOptionToggle(Opts.names['bTaperChangePerCrvParam'], Opts.riOpts['bTaperChangePerCrvParam'])
    idxs_Opts['FlipAngle'] = go.AddOption('FlipAngle')
    idxs_Opts['FlipDir'] = go.AddOption('FlipDir')
    idxs_Opts['pt_Base'] = go.AddOption(Opts.names['pt_Base'])
    go.AddOptionToggle(Opts.names['bCPlane'], Opts.riOpts['bCPlane'])
    go.AddOptionToggle(Opts.names['bAtGrevilles'], Opts.riOpts['bAtGrevilles'])
    go.AddOptionToggle(Opts.names['bAtKnots'], Opts.riOpts['bAtKnots'])
    go.AddOptionToggle(Opts.names['bAtEqualDivisions'], Opts.riOpts['bAtEqualDivisions'])
    if Opts.values['bAtEqualDivisions']:
        go.AddOptionInteger(Opts.names['iDivisionCt'], Opts.riOpts['iDivisionCt'])
    go.AddOptionToggle(Opts.names['bSplitPolyCrvToSegs'], Opts.riOpts['bSplitPolyCrvToSegs'])
    go.AddOptionToggle(Opts.names['bSplitPathsAtKnots'], Opts.riOpts['bSplitPathsAtKnots'])
    go.AddOptionToggle(Opts.names['bLoftCrvs'], Opts.riOpts['bLoftCrvs'])
    if Opts.values['bLoftCrvs']:
        idxs_Opts['iLoftType'] = go.AddOptionList(
                englishOptionName=Opts.names['iLoftType'],
                listValues=Enum.GetNames(rg.LoftType),
                listCurrentIndex=Opts.values['iLoftType'])
    if Opts.values['bLoftCrvs']:
        go.AddOptionToggle(Opts.names['bAddArrayedObjects'], Opts.riOpts['bAddArrayedObjects'])
    go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
    go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
    go.AlreadySelectedObjectSelect = True # So objects can be reselected after being unselected in same go.GetMultiple.
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
    if go.ObjectsWerePreselected:
        objrefs = None
        iCt_Crvs_PreSelctd = 0
    else:
        if rdObjectType == rd.ObjectType.AnyObject:
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
                return False
        objrefs = go.Objects()
        go.Dispose()
        return tuple([objrefs] + [Opts.values[key] for key in sOpts])
    elif res == ri.GetResult.Cancel:
        go.Dispose()
        return
    elif res == ri.GetResult.Number:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue = go.Number()
    elif Opts.values['bVariableTaper'] and go.OptionIndex() == idxs_Opts['SwapAngles']:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue, Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                Opts.riOpts['fTaper_End_Deg'].CurrentValue, Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
    elif go.OptionIndex() == idxs_Opts['FlipAngle']:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
        Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
    elif go.OptionIndex() == idxs_Opts['FlipDir']:
        if Opts.riOpts['fTaper_Start_Deg'].CurrentValue == 0.0:
            Opts.riOpts['fTaper_Start_Deg'].CurrentValue = 180.0
        elif Opts.riOpts['fTaper_Start_Deg'].CurrentValue < 0.0:
            Opts.riOpts['fTaper_Start_Deg'].CurrentValue = (
                    -180.0 - Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
        else:
            Opts.riOpts['fTaper_Start_Deg'].CurrentValue = (
                    180.0 - Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
            
        if Opts.riOpts['fTaper_End_Deg'].CurrentValue == 0.0:
            Opts.riOpts['fTaper_End_Deg'].CurrentValue = 180.0
        elif Opts.riOpts['fTaper_End_Deg'].CurrentValue < 0.0:
            Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                    -180.0 - Opts.riOpts['fTaper_End_Deg'].CurrentValue)
        else:
            Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                    180.0 - Opts.riOpts['fTaper_End_Deg'].CurrentValue)
    elif go.OptionIndex() == idxs_Opts['pt_Base']:
        cmdres, pt = ri.RhinoGet.GetPoint(
            prompt="Base point for objects to array", acceptNothing=False)
        if cmdres == Rhino.Commands.Result.Success:
            Opts.values['pt_Base'] = pt

    Opts.setValues()
    Opts.saveSticky()
    
    return tuple([objrefs] + [Opts.values[key] for key in sOpts])


def prepareCurves(rgCrvs0, bRebuild=False, bSplitPolyCrvToSegs=True, bSplitPathsAtKnots=False):
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

        if not bSplitPolyCrvToSegs:
            if not bRebuild:
                rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()]
            else:
                rc = spb_Crv_fitRebuild.rebuildCurve(
                    rgCrv_Joined,
                    0.1*sc.doc.ModelAbsoluteTolerance,
                    iDegree=3,
                    bPreserveEndTans=True,
                    bFurtherTranslateCps=False,
                    iMinCpCt=None,
                    iMaxCpCt=50,
                    bDebug=False,
                    )
                rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()] if rc[0] is None else [rc[0]]
        else:
            if isinstance(rgCrv_Joined, rg.PolyCurve):
                rc = rgCrv_Joined.Explode()
                if rc:
                    rgCrvs_Exploded = rc
                    if not bRebuild:
                        rgCrvs_SplitPoly = rc
                    else:
                        rgCrvs_SplitPoly = []
                        for c in rgCrvs_Exploded:
                            rc = spb_Crv_fitRebuild.rebuildCurve(
                                c,
                                0.1*sc.doc.ModelAbsoluteTolerance,
                                iDegree=3,
                                bPreserveEndTans=True,
                                bFurtherTranslateCps=False,
                                iMinCpCt=None,
                                iMaxCpCt=50,
                                bDebug=False,
                                )
                            rgCrvs_SplitPoly.append(c if rc[0] is None else rc[0])
            else:
                # not PolyCurve.
                if not bRebuild:
                    rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()]
                else:
                    rc = spb_Crv_fitRebuild.rebuildCurve(
                        rgCrv_Joined,
                        0.1*sc.doc.ModelAbsoluteTolerance,
                        iDegree=3,
                        bPreserveEndTans=True,
                        bFurtherTranslateCps=False,
                        iMinCpCt=None,
                        iMaxCpCt=50,
                        bDebug=False,
                        )
                    rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()] if rc[0] is None else [rc[0]]

        if not bSplitPathsAtKnots:
            if not bRebuild:
                rgCrvs_Final.extend(rgCrvs_SplitPoly)
            else:
                rgCrvs_Final = []
                for rgCrv_SplitPoly in rgCrvs_SplitPoly:
                    rc = spb_Crv_fitRebuild.rebuildCurve(
                        rgCrv_SplitPoly,
                        0.1*sc.doc.ModelAbsoluteTolerance,
                        iDegree=3,
                        bPreserveEndTans=True,
                        bFurtherTranslateCps=False,
                        iMinCpCt=None,
                        iMaxCpCt=50,
                        bDebug=False,
                        )
                    rgCrvs_Final.append(rgCrv_SplitPoly.Duplicate() if rc[0] is None else rc[0])
        else:
            for rgCrv_SplitPoly in rgCrvs_SplitPoly:
                ts_SpanBoundaries = [rgCrv_SplitPoly.Domain.T0]
                for iSpan in xrange(rgCrv_SplitPoly.SpanCount):
                    ts_SpanBoundaries.append(rgCrv_SplitPoly.SpanDomain(iSpan).T1)
                rc = rgCrv_SplitPoly.Split(ts_SpanBoundaries)
                if rc:
                    rgCrvs_SplitAtKnots = rc
                    if not bRebuild:
                        rgCrvs_Final.extend(rgCrvs_SplitAtKnots)
                    else:
                        rgCrvs_Final = []
                        for rgCrv_SplitAtKnots in rgCrvs_SplitAtKnots:
                            rc = spb_Crv_fitRebuild.rebuildCurve(
                                rgCrv_SplitAtKnots,
                                0.1*sc.doc.ModelAbsoluteTolerance,
                                iDegree=3,
                                bPreserveEndTans=True,
                                bFurtherTranslateCps=False,
                                iMinCpCt=None,
                                iMaxCpCt=50,
                                bDebug=False
                                )
                            rgCrvs_Final.append(rgCrv_SplitAtKnots if rc[0] is None else rc[0])

                rgCrv_SplitPoly.Dispose()
        rgCrv_Joined.Dispose()
    
    return rgCrvs_Final


def createArrayedGeometry(rgCrv_Path, rgObjs_ToArray, plane_Proj, fTaper_Start_Deg, fTaper_End_Deg, bTaperChangePerCrvParam, bAtGrevilles, bAtKnots, bAtEqualDivisions, iDivisionCt=None, bDebug=False):
    """
    """

    xform_Proj = rg.Transform.PlanarProjection(plane_Proj)

    rgObjs1_Arrayed_PerSection = []

    nc2_Path = rgCrv_Path.ToNurbsCurve()
    if nc2_Path is None:
        print "NurbsCurve could not be calculated from curve."
        return
    #sc.doc.Objects.AddCurve(nc2_Path)
        
    ts = []

    if bAtGrevilles:
        rc = nc2_Path.GrevilleParameters()
        if rc: ts.extend(rc)
        if nc2_Path.IsClosed:
            ts.pop()

    if bAtEqualDivisions:
        rc = nc2_Path.DivideByCount(
                segmentCount=iDivisionCt,
                includeEnds=True)
        if rc: ts.extend(rc)
        if nc2_Path.IsClosed:
            ts.append(nc2_Path.Domain.T1)
    
    if bAtKnots:
        rc = set(nc2_Path.Knots)
        if rc: ts.extend(rc)
        
    if ts is None:
        print "No parameters were obtained."
        return

    ts = sorted(set(ts)) # Remove duplicates and sort.

    # Remove overlaps for closed (including periodic) curves.
    if nc2_Path.IsClosed:
        if bDebug:
            print nc2_Path.Domain
            print ts
        ts_WIP = []
        for t in ts:
            if t >= nc2_Path.Domain.T0 and t < nc2_Path.Domain.T1:
                ts_WIP.append(t)
        ts = ts_WIP
        if bDebug: print ts

    if fTaper_Start_Deg == fTaper_End_Deg:
        angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
    elif not bTaperChangePerCrvParam:
        length_Full = rgCrv_Path.GetLength()

    rgCrv_Path_Flattened = nc2_Path.Duplicate()
    rgCrv_Path_Flattened.Transform(xform_Proj)

    #sc.doc.Objects.AddCurve(rgCrv_Path_Flattened); return
    rgObjs1_Arrayed_PerSection = []
    for iT, t in enumerate(ts):
        bSuccess, frame = rgCrv_Path_Flattened.PerpendicularFrameAt(t=t)
        if not bSuccess:
            print "Perpendicular frame could not be calculated."
            continue
                
        angle_StraightenFrame_Rad = rg.Vector3d.VectorAngle(frame.YAxis, plane_Proj.ZAxis, frame)

        if fTaper_Start_Deg != fTaper_End_Deg:
            if t == rgCrv_Path.Domain.T0:
                angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
            elif t == rgCrv_Path.Domain.T1:
                angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_End_Deg)
            #elif len(ts) >= 4 and iT == 1:
            #    # Same angle so that the brep starts at a tangent.
            #    angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_Start_Deg)
            #elif len(ts) >= 4 and iT == len(ts) - 2:
            #    # Same angle so that the brep end at a tangent.
            #    angle_Rad = Rhino.RhinoMath.ToRadians(fTaper_End_Deg)

            else:
                if bTaperChangePerCrvParam:
                    t_Normalized = rgCrv_Path.Domain.NormalizedParameterAt(t)
                    angle_Rad = Rhino.RhinoMath.ToRadians(
                            fTaper_Start_Deg * (1.0 - t_Normalized) +
                            fTaper_End_Deg * t_Normalized)
                else:
                    length_to_t = rgCrv_Path.GetLength(
                            subdomain=rg.Interval(rgCrv_Path.Domain.T0, t))
                    angle_Rad = Rhino.RhinoMath.ToRadians(
                            fTaper_Start_Deg * (1.0 - length_to_t/length_Full) +
                            fTaper_End_Deg * length_to_t/length_Full)


            if bDebug: print angle_Rad
        frame.Rotate(angle=angle_StraightenFrame_Rad+angle_Rad, axis=frame.ZAxis)
                
        # Debug frame orientation.
        #                sc.doc.Objects.AddPoint(frame.PointAt(0.0,0.0,0.0))
        #                attr_Red = rd.ObjectAttributes()
        #                attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
        #                attr_Red.ObjectColor = Color.Red
        #                sc.doc.Objects.AddPoint(frame.PointAt(1.0,0.0,0.0), attr_Red)
        #                attr_Green = rd.ObjectAttributes()
        #                attr_Green.ColorSource = rd.ObjectColorSource.ColorFromObject
        #                attr_Green.ObjectColor = Color.Lime
        #                sc.doc.Objects.AddPoint(frame.PointAt(0.0,1.0,0.0), attr_Green)
                
        xform1 = rg.Transform.PlaneToPlane(plane_Proj, frame)
        xform2 = rg.Transform.Translation(
                nc2_Path.PointAt(t)-rgCrv_Path_Flattened.PointAt(t))
        xform3 = xform2 * xform1
        rgObjs1_Arrayed_PerSection.append([])

        for rgObj0_ToArray in rgObjs_ToArray:
            rgObj1_Arrayed_WIP = rgObj0_ToArray.Duplicate()
            rgObj1_Arrayed_WIP.Transform(xform3)
            rgObjs1_Arrayed_PerSection[-1].append(rgObj1_Arrayed_WIP)
        if not rgObjs1_Arrayed_PerSection[-1]:
            del rgObjs1_Arrayed_PerSection[-1]

    nc2_Path.Dispose()
    rgCrv_Path_Flattened.Dispose()

    return rgObjs1_Arrayed_PerSection


def createLofts(rgObjs_ToLoft_PerSect, iLoftType, bClosedLoft):

    rgBreps_Lofts = []

    for i in xrange(len(rgObjs_ToLoft_PerSect[0])):
        if not isinstance(rgObjs_ToLoft_PerSect[0][i], rg.Curve):
            continue

        rgCrvsToLoft = []
        for rgObjs_ToLoft_1Sect in rgObjs_ToLoft_PerSect:
            rgCrvsToLoft.append(rgObjs_ToLoft_1Sect[i])

        rc = rg.Brep.CreateFromLoft(
                curves=rgCrvsToLoft,
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=Enum.ToObject(rg.LoftType, iLoftType),
                closed=bClosedLoft)
        if rc:
            rgBreps_Lofts.extend(rc)

    return rgBreps_Lofts


def main():
    
    while True:
        rc = getInput("Select objects to array", rd.ObjectType.AnyObject)
        if rc is None: return
        if rc[0] is None: continue
        objrefs_ToArray = rc[0]
        for key, value in zip(sOpts, rc[1:]):
            exec("{} = {}".format(key, value))
        break

    sc.doc.Objects.UnselectAll()

    rc = getInput("Select path curves", rd.ObjectType.Curve)
    if rc is None: return
    objrefs_Paths = rc[0]
    for key, value in zip(sOpts, rc[1:]):
        exec("{} = {}".format(key, value))

    rgObjs0_ToArray = []
    bCrvInToArray = False
    for gObj0_ToArray in objrefs_ToArray:
        rgObj0_ToArray = gObj0_ToArray.Geometry()
        if rgObj0_ToArray:
            rgObjs0_ToArray.append(rgObj0_ToArray)
            if not bCrvInToArray and isinstance(rgObj0_ToArray, rg.Curve):
                bCrvInToArray = True

    rgCrvs0_Path = []

    while True:

        while not (bAtGrevilles or bAtKnots or bAtEqualDivisions):
            print "No path point sampling is enabled."
            sc.doc.Views.Redraw()
            rc = getInput("Select path curves", rd.ObjectType.Curve)
            if rc is None or rc is False:
                for o in rgObjs0_ToArray: o.Dispose()
                for c in rgCrvs0_Path: c.Dispose()
                return
            objrefs_Paths = rc[0]
            for key, value in zip(sOpts, rc[1:]):
                exec("{} = {}".format(key, value))
        
        rgCrvs0_Path = []
        for objref_Path in objrefs_Paths:
            c = objref_Path.Curve()
            rgCrvs0_Path.append(c)

        Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")
        rc = prepareCurves(
                rgCrvs0=rgCrvs0_Path,
                bRebuild=bRebuildPath,
                bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
                bSplitPathsAtKnots=bSplitPathsAtKnots)
        if not rc: return
        rgCrvs1_Path = rc

        if bCPlane:
            view_Active = sc.doc.Views.ActiveView
            plane_Proj = view_Active.ActiveViewport.ConstructionPlane()
            #plane_Proj.Origin = pt_Base
            #sc.doc.Objects.AddPoint(plane_Proj.Origin)
        else:
            plane_Proj = rg.Plane.WorldXY

        gObjs1_Arrayed_All = []
        rgBreps1_Lofts_PerSect_PerSeg = []

        Rhino.RhinoApp.SetCommandPrompt("Creating geometry ...")

        for rgCrv1_Path in rgCrvs1_Path:
            rc = createArrayedGeometry(
                    rgCrv_Path=rgCrv1_Path,
                    rgObjs_ToArray=rgObjs0_ToArray,
                    plane_Proj=plane_Proj,
                    fTaper_Start_Deg=fTaper_Start_Deg,
                    fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                    bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                    bAtGrevilles=bAtGrevilles,
                    bAtKnots=bAtKnots,
                    bAtEqualDivisions=bAtEqualDivisions,
                    iDivisionCt=iDivisionCt,
                    bDebug=bDebug)
            for rgO in rgObjs0_ToArray: rgO.Dispose()
            if rc is None: continue
            rgObjs_Arrayed_PerSect_1PathSeg = rc
    
            if bAddArrayedObjects or not bLoftCrvs:
                for rgObjs_Arrayed_1Sect in rgObjs_Arrayed_PerSect_1PathSeg:
                    for i, rgObj_Arrayed_1Sect in enumerate(rgObjs_Arrayed_1Sect):
                        gObj_Arrayed = sc.doc.Objects.Add(rgObj_Arrayed_1Sect)
                        if gObj_Arrayed != Guid.Empty:
                            gObjs1_Arrayed_All.append(gObj_Arrayed)
                            #rgDot_ = rg.TextDot(str(i), rgObj1_Arrayed_1Seg.PointAtStart)
                            #sc.doc.Objects.AddTextDot(rgDot_)

            if bLoftCrvs and bCrvInToArray:
                rc = createLofts(
                        rgObjs_ToLoft_PerSect=rgObjs_Arrayed_PerSect_1PathSeg,
                        iLoftType=iLoftType,
                        bClosedLoft=rgCrv1_Path.IsClosed)
                if rc:
                    rgBreps1_Lofts_PerSect_PerSeg.extend(rc)
            for rgOs in rgObjs_Arrayed_PerSect_1PathSeg:
                for rgO in rgOs:
                    rgO.Dispose()

        s = ""
        if gObjs1_Arrayed_All and bEcho:
            s += "Added {} arrayed objects.".format(len(gObjs1_Arrayed_All))

        gBreps1 = []
        iCt_Faces_All = 0
        if rgBreps1_Lofts_PerSect_PerSeg:
            rc = rg.Brep.JoinBreps(
                    rgBreps1_Lofts_PerSect_PerSeg,
                    tolerance=0.5*sc.doc.ModelAbsoluteTolerance)
            if rc:
                for rgBrep1 in rc:
                    gBrep1 = sc.doc.Objects.AddBrep(rgBrep1)
                    if gBrep1 != Guid.Empty:
                        gBreps1.append(gBrep1)
                        iCt_Faces_All += rgBrep1.Faces.Count

        if gBreps1 and bEcho:
            if s: s += "  "
            s += "Added {} brep(s) with {} total face(s).".format(
                    len(gBreps1), iCt_Faces_All)
            print s

        sc.doc.Views.Redraw()

        rc = getInput("Select path curves", rd.ObjectType.Curve)
        if rc is None:
            for gO in gObjs1_Arrayed_All:
                sc.doc.Objects.Delete(gO, True)
            for gBrep1 in gBreps1:
                sc.doc.Objects.Delete(gBrep1, True)
            break
        elif rc is False:
            break
        objrefs_Paths = rc[0]
        for key, value in zip(sOpts, rc[1:]):
            exec("{} = {}".format(key, value))

        for gO in gObjs1_Arrayed_All:
            sc.doc.Objects.Delete(gO, True)
        for gBrep1 in gBreps1:
            sc.doc.Objects.Delete(gBrep1, True)

    for c in rgCrvs0_Path: c.Dispose()
    for c in rgCrvs1_Path: c.Dispose()
    for _s in rgObjs_Arrayed_PerSect_1PathSeg:
        for _ in _s:
            _.Dispose()



if __name__ == '__main__': main()