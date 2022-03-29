"""
"""

"""
200402-04: Created, starting with other scripts.
220328: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum
from System import Guid

import spb_Crv_fitRebuild


sBrepMethods = 'Loft2Crvs', 'LoftSectionLines', 'Sweep2A', 'Sweep2B', 'Network'

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


    key = 'bRebuildRails'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCPlane'; keys.append(key)
    values[key] = False
    names[key] = 'PlanView'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAtGrevilles'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAtEqualDivisions'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iDivisionCt'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(
            initialValue=values[key],
            setLowerLimit=True,
            limit=2)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSplitPolyCrvToSegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bSplitPathsAtKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddArrayedLines'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddBrep'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iBrepMethod'; keys.append(key)
    values[key] = 0
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iLoftType'; keys.append(key)
    values[key] = 0
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fBrepTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
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


def get2Curves():
    """
    Get 2 curves.  No options.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves")

    go.GeometryFilter = rd.ObjectType.Curve

    while True:
        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs


def prepareCurves(rgCs_In, bRebuild=False, bSplitPolyCrvToSegs=True, bSplitPathsAtKnots=False):
    """
    Prepare curves, including splitting per options.

    Returns: list of new curves
    """
    
    rc = rg.Curve.JoinCurves(
            rgCs_In,
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
                    0.25*sc.doc.ModelAbsoluteTolerance,
                    iDegree=3,
                    bPreserveEndTans=True,
                    bFurtherTranslateCps=False,
                    iMinCpCt=None,
                    iMaxCpCt=50,
                    bDebug=False,
                    )
                rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()] if rc[0] is None else [rc[0]]
        else:
            # bSplitPolyCrvToSegs == True.
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
                                0.25*sc.doc.ModelAbsoluteTolerance,
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
                        0.25*sc.doc.ModelAbsoluteTolerance,
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
                        0.25*sc.doc.ModelAbsoluteTolerance,
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
                                0.25*sc.doc.ModelAbsoluteTolerance,
                                iDegree=3,
                                bPreserveEndTans=True,
                                bFurtherTranslateCps=False,
                                iMinCpCt=None,
                                iMaxCpCt=50,
                                bDebug=False
                                )
                            rgCrvs_Final.append(rgCrv_SplitAtKnots if rc[0] is None else rc[0])

                rgCrv_SplitPoly.Dispose()
    
    return rgCrvs_Final


def getParameters(nc, bAtGrevilles, iDivisionCt=None, fDivisionLength=None):


    ts = []

    #if bAtGrevilles:
    #    rc = nc.GrevilleParameters()
    #    if rc: ts.extend(rc)
    #    if nc.IsClosed:
    #        ts.pop()


    #if iDivisionCt:
    #    rc = nc.DivideByCount(
    #            segmentCount=iDivisionCt,
    #            includeEnds=True)
    #    if rc: ts.extend(rc)
    #    if nc.IsClosed:
    #        ts.append(nc.Domain.T1)
    
    fDivisionLength = 10.0*sc.doc.ModelAbsoluteTolerance
    if fDivisionLength:
        rc = nc.DivideByLength(
            segmentLength=fDivisionLength,
            includeEnds=True)
        print rc[-1]
        if rc: ts.extend(rc)
        if not nc.IsClosed:
            # DivideByLength doesn't add the T1 segment
            # even though includeEnds == True.
            # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
            # shows the parameter labeled as 'includeStart'.
            ts.append(nc.Domain.T1)

    if ts is None:
        print "No parameters were obtained."
        return

    ts = sorted(set(ts)) # Remove duplicates and sort.

    # Remove overlaps for closed (including periodic) curves.
    if nc.IsClosed:
        ts_WIP = []
        for t in ts:
            if t >= nc.Domain.T0 and t < nc.Domain.T1:
                ts_WIP.append(t)
        ts = ts_WIP

    return ts


def createCrossSectionLines_NoProjection(ncs_A, ncs_B, ts_perA, bDebug=False):
    """
    """

    lines_Out_perA = []
    ncs_B_Out_perA = []

    for iA in range(len(ncs_A)):
        ncA = ncs_A[iA]
        tsA = ts_perA[iA]
        
        lines_Out_perA.append([])
        ncs_B_Out_perA.append([None])

        t_B_StartForA = None
        nc_B_StartForA = None
        t_B_EndForA = None
        nc_B_EndForA = None
            
        for iT, tA in enumerate(tsA):

            if tA > 26.1:
                pass
                #bDebug = True

            bSuccess, frame = ncA.PerpendicularFrameAt(t=tA)
            if not bSuccess:
                print "Perpendicular frame could not be calculated."
                continue

            pt_A = ncA.PointAt(tA)

            for iB in range(len(ncs_B)):
                
                nc_B = ncs_B[iB]

                rc = rg.Intersect.Intersection.CurvePlane(
                    nc_B,
                    frame,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if not rc:
                    # No intersection, so don't bother proceeding with this nc_B.
                    continue # to next ncs_B.

                res, tB = nc_B.ClosestPoint(pt_A)
                if not res:
                    raise ValueError("ClosestPoint returned {}.".format(res))

                pt_B_ClosestPt = nc_B.PointAt(tB)
                if bDebug: sc.doc.Objects.AddPoint(pt_B_ClosestPt)

                for intersection in rc:
                    ptB_PerpToA = intersection.PointA
                    if bDebug: sc.doc.Objects.AddPoint(ptB_PerpToA)

                    dist = ptB_PerpToA.DistanceTo(pt_B_ClosestPt)
                    if bDebug: print dist

                    if dist <= 1.0*sc.doc.ModelAbsoluteTolerance:
                        if bDebug: sc.doc.Objects.AddLine(rg.Line(pt_A, ptB_PerpToA))
                        line = rg.Line(ncA.PointAt(tA), to=nc_B.PointAt(tB))
                        if bDebug: sc.doc.Objects.AddLine(line)
                        lines_Out_perA[-1].append(line)
                        if t_B_StartForA is None:
                            nc_B_StartForA = nc_B
                            t_B_StartForA = tB
                        if bDebug: sc.doc.Views.Redraw(); return
                        break # to next ncB.

    return lines_Out_perA


def createCrossSectionLines_PerProjection(ncs_A, ncs_B, ts_perA, plane_Proj, bDebug=False):
    """
    """

    xform_Proj = rg.Transform.PlanarProjection(plane_Proj)

    ncs_B_Flattened = []

    for nc_B in ncs_B:
        ncB_Flattened = nc_B.Duplicate()
        ncB_Flattened.Transform(xform_Proj)
        #sc.doc.Objects.AddCurve(ncB_Flattened)
        ncs_B_Flattened.append(ncB_Flattened)

    lines_Out_perA = []

    for iA in range(len(ncs_A)):
        ncA = ncs_A[iA]
        tsA = ts_perA[iA]
        
        ncA_Flattened = ncA.Duplicate()
        ncA_Flattened.Transform(xform_Proj)
        #sc.doc.Objects.AddCurve(ncA_Flattened)

        lines_Out_perA.append([])

        for iT, tA in enumerate(tsA):

            if tA > 26.1:
                pass
                #bDebug = True

            bSuccess, frame = ncA_Flattened.PerpendicularFrameAt(t=tA)
            if not bSuccess:
                print "Perpendicular frame could not be calculated."
                continue

            ptA_Flat = ncA_Flattened.PointAt(tA)
            
            for iB in range(len(ncs_B_Flattened)):
                
                ncB_Flat = ncs_B_Flattened[iB]
                if bDebug: sc.doc.Objects.AddCurve(ncB_Flat)
                nc_B = ncs_B[iB]

                rc = rg.Intersect.Intersection.CurvePlane(
                    ncB_Flat,
                    frame,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
                if not rc:
                    # No intersection, so don't bother getting ClosestPoint.
                    continue

                res, tB = ncB_Flat.ClosestPoint(ptA_Flat)
                if not res: continue

                ptB_Flat_ClosestPt = ncB_Flat.PointAt(tB)
                if bDebug: sc.doc.Objects.AddPoint(ptB_Flat_ClosestPt)

                for intersection in rc:
                    ptB_Flat_PerpToA = intersection.PointA
                    if bDebug: sc.doc.Objects.AddPoint(ptB_Flat_PerpToA)

                    dist = ptB_Flat_PerpToA.DistanceTo(ptB_Flat_ClosestPt)
                    if bDebug: print dist

                    if dist <= 1.0*sc.doc.ModelAbsoluteTolerance:
                        if bDebug: sc.doc.Objects.AddLine(rg.Line(ptA_Flat, ptB_Flat_PerpToA))
                        line = rg.Line(ncA.PointAt(tA), to=nc_B.PointAt(tB))
                        if bDebug: sc.doc.Objects.AddLine(line)
                        lines_Out_perA[-1].append(line)
                        if bDebug: sc.doc.Views.Redraw(); return
                        break # to next ncB_Flat.

        ncA_Flattened.Dispose()


    for nc in ncs_B_Flattened:
        nc.Dispose()

    return lines_Out_perA


def createBrep(iBrepMethod, iLoftType, fBrepTol, ncs_A, ncs_B, lines_perA):
    """
    """

    rgBreps1 = []

    for iA in range(len(ncs_A)):
        nc_A = ncs_A[iA]
        lines_A = lines_perA[iA]
    
        if sBrepMethods[Opts.values['iBrepMethod']] == 'Loft2Crvs':
            rgBreps1 = rg.Brep.CreateFromLoft(
                    curves=[nc_A, rgNurbsCrv_TaperEnd_1Seg],
                    start=rg.Point3d.Unset,
                    end=rg.Point3d.Unset,
                    loftType=rg.LoftType.Straight,
                    closed=False)
        elif sBrepMethods[Opts.values['iBrepMethod']] == 'LoftSectionLines':
            for L in lines_A:
                print L.PointAtStart, L.PointAtEnd
            print lines_A[0].PointAtEnd.EpsilonEquals(lines_A[-1].PointAtEnd, epsilon=1e-12)
            rgBreps1 = rg.Brep.CreateFromLoft(
                    curves=lines_A,
                    start=rg.Point3d.Unset,
                    end=rg.Point3d.Unset,
                    loftType=Enum.ToObject(rg.LoftType, iLoftType),
                    closed=rgNurbsCrv1_PathSeg.IsClosed)
            print rgBreps1
        elif sBrepMethods[Opts.values['iBrepMethod']] == 'Sweep2A':
            rgBreps1 = rg.Brep.CreateFromSweep(
                    rail1=rgNurbsCrv1_PathSeg,
                    rail2=rgNurbsCrv_TaperEnd_1Seg,
                    shapes=lines_A,
                    closed=rgNurbsCrv1_PathSeg.IsClosed,
                    tolerance=fBrepTol)
        elif sBrepMethods[Opts.values['iBrepMethod']] == 'Sweep2B':
            rgSweep2 = rg.SweepTwoRail()
            #rgSweep2.AngleToleranceRadians
            rgSweep2.ClosedSweep = rgNurbsCrv1_PathSeg.IsClosed
            rgSweep2.MaintainHeight = False
            rgSweep2.SweepTolerance = fBrepTol
            rgBreps1 = rgSweep2.PerformSweep(
                    rail1=rgNurbsCrv1_PathSeg,
                    rail2=rgNurbsCrv_TaperEnd_1Seg,
                    crossSections=lines_A)
        elif sBrepMethods[Opts.values['iBrepMethod']] == 'Network':
            rgNurbsSrf, iError = rg.NurbsSurface.CreateNetworkSurface(
                    curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg]+lines_A,
                    continuity=1,
                    edgeTolerance=fBrepTol,
                    interiorTolerance=fBrepTol,
                    angleTolerance=0.1*sc.doc.ModelAngleToleranceDegrees)
            if iError:
                print "CreateNetworkSurface error code: {}".format(iError)
            else:
                rgBreps1 = [rgNurbsSrf.ToBrep()]
                rgNurbsSrf.Dispose()

    return rgBreps1


def createGeometry(rgCs_A_In, rgCs_B_In, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bRebuildRails = getOpt('bRebuildRails')
    bSplitPolyCrvToSegs = getOpt('bSplitPolyCrvToSegs')
    bSplitPathsAtKnots = getOpt('bSplitPathsAtKnots')
    bAtGrevilles = getOpt('bAtGrevilles')
    iDivisionCt = getOpt('iDivisionCt')
    bCPlane = getOpt('bCPlane')
    iBrepMethod = getOpt('iBrepMethod')
    iLoftType = getOpt('iLoftType')
    fBrepTol = getOpt('fBrepTol')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')



    Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")

    rgCs_A_Prepped = prepareCurves(
        rgCs_A_In,
        bRebuild=bRebuildRails,
        bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
        bSplitPathsAtKnots=bSplitPathsAtKnots)

    rgCs_B_Prepped = prepareCurves(
        rgCs_B_In,
        bRebuild=bRebuildRails,
        bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
        bSplitPathsAtKnots=bSplitPathsAtKnots)


    ncs_A = [c.ToNurbsCurve() for c in rgCs_A_Prepped]
    ncs_B = [c.ToNurbsCurve() for c in rgCs_B_Prepped]


    ts_perA = []
    for nc in ncs_A:
        ts = getParameters(
            nc,
            bAtGrevilles=bAtGrevilles,
            iDivisionCt=iDivisionCt,
            fDivisionLength=None)
        if ts is None:
            print "Parameters could not be obtained for curve."
            return
        ts_perA.append(ts)


    if bCPlane:
        view_Active = sc.doc.Views.ActiveView
        plane_Proj = view_Active.ActiveViewport.ConstructionPlane()
        
        xform1 = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, plane_Proj)
        rgLine_ToArray.Transform(xform1)
    else:
        plane_Proj = rg.Plane.WorldXY


    rc = createCrossSectionLines_NoProjection(
        ncs_A,
        ncs_B,
        ts_perA,
        bDebug=bDebug)

    #rc = createCrossSectionLines_PerProjection(
    #    ncs_A,
    #    ncs_B,
    #    ts_perA,
    #    plane_Proj,
    #    bDebug=bDebug)

    lines_perA = rc

    #brep = createBrep(
    #        iBrepMethod=iBrepMethod,
    #        iLoftType=iLoftType,
    #        fBrepTol=fBrepTol,
    #        ncs_A=ncs_A,
    #        ncs_B=ncs_B,
    #        lines_perA=lines_perA)

    return lines_perA

    return lines_perA, brep


class DrawLinesConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.color = sc.doc.Layers.CurrentLayer.Color
        self.lines = None

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if self.lines:
            self.bbox = rg.BoundingBox(points=[line.From for line in self.lines])
            calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):
        if self.lines:
            drawEventArgs.Display.DrawLines(
                lines=self.lines,
                color=self.color,
                thickness=1)


def getOptions(rgCs_A_In, rgCs_B_In):
    """
    Get options.
    
    Returns
        None to cancel.
        False to indicate to create objects with current options.
        True to indicate to regenerate geometry and return to this function.
    """
    
    valuesBefore = {}

    go = ri.Custom.GetOption()
    go.SetCommandPrompt("Set options")

    go.AcceptNothing(True)

    rc = createGeometry(
        rgCs_A_In,
        rgCs_B_In,
        iDivisionCt=Opts.values['iDivisionCt'] if Opts.values['bAtEqualDivisions'] else None,
        )

    conduit = DrawLinesConduit()

    if rc:
        conduit.lines = [line for lines in rc for line in lines]
        conduit.Enabled = True
        sc.doc.Views.Redraw()


    idxs_Opts = {}

    while True:
        sc.escape_test()

        for key in Opts.keys:
            valuesBefore[key] = Opts.values[key]

        go.AddOptionToggle(Opts.names['bRebuildRails'], Opts.riOpts['bRebuildRails'])
        go.AddOptionToggle(Opts.names['bCPlane'], Opts.riOpts['bCPlane'])
        go.AddOptionToggle(Opts.names['bAtGrevilles'], Opts.riOpts['bAtGrevilles'])
        go.AddOptionToggle(Opts.names['bAtEqualDivisions'],
                            Opts.riOpts['bAtEqualDivisions'])
        if Opts.values['bAtEqualDivisions']:
            go.AddOptionInteger(Opts.names['iDivisionCt'],
                                Opts.riOpts['iDivisionCt'])
        go.AddOptionToggle(Opts.names['bSplitPolyCrvToSegs'], Opts.riOpts['bSplitPolyCrvToSegs'])
        go.AddOptionToggle(Opts.names['bSplitPathsAtKnots'], Opts.riOpts['bSplitPathsAtKnots'])
        go.AddOptionToggle(Opts.names['bAddArrayedLines'], Opts.riOpts['bAddArrayedLines'])
        go.AddOptionToggle(Opts.names['bAddBrep'], Opts.riOpts['bAddBrep'])
        if Opts.values['bAddBrep']:
            idxs_Opts['iBrepMethod'] = go.AddOptionList(
                    englishOptionName=Opts.names['iBrepMethod'],
                    listValues=sBrepMethods,
                    listCurrentIndex=Opts.values['iBrepMethod'])
            if sBrepMethods[Opts.values['iBrepMethod']] == 'LoftSectionLines':
                idxs_Opts['iLoftType'] = go.AddOptionList(
                        englishOptionName=Opts.names['iLoftType'],
                        listValues=Enum.GetNames(rg.LoftType),
                        listCurrentIndex=Opts.values['iLoftType'])
            go.AddOptionDouble(Opts.names['fBrepTol'], Opts.riOpts['fBrepTol'])
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        
        res = go.Get()

        conduit.Enabled = False
        sc.doc.Views.Redraw()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Nothing:
            # Accept current result.
            return rc
        elif Opts.values['bAddBrep'] and go.OptionIndex() == idxs_Opts['iBrepMethod']:
            Opts.values['iBrepMethod'] = go.Option().CurrentListOptionIndex
        elif Opts.values['bAddBrep'] and Opts.values['iBrepMethod'] == 1 and go.OptionIndex() == idxs_Opts['iLoftType']:
            Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
        elif Opts.riOpts['fBrepTol'].CurrentValue < 0.0:
            Opts.riOpts['fBrepTol'].CurrentValue = Opts.riOpts['fBrepTol'].InitialValue

        Opts.setValues()
        Opts.saveSticky()

        for key in Opts.keys:
            if valuesBefore[key] != Opts.values[key]:
                rc = createGeometry(
                    rgCs_A_In,
                    rgCs_B_In,
                    iDivisionCt=Opts.values['iDivisionCt'] if Opts.values['bAtEqualDivisions'] else None,
                    )
                if rc:
                    conduit.lines = [line for lines in rc for line in lines]
                    conduit.Enabled = True
                    sc.doc.Views.Redraw()
                break
        else:
            for key in Opts.keys:
                print valuesBefore[key], Opts.values[key]
            print "No options were changed."


def main():

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select first rail curves",
        acceptNothing=False,
        filter=rd.ObjectType.Curve)
    if res == Rhino.Commands.Result.Cancel: return
    
    rgCs_A_In = [o.Curve() for o in objrefs]

    sc.doc.Objects.UnselectAll()

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select second rail curves",
        acceptNothing=False,
        filter=rd.ObjectType.Curve)
    if res == Rhino.Commands.Result.Cancel: return

    rgCs_B_In = [o.Curve() for o in objrefs]

    sc.doc.Objects.UnselectAll()

    #rc = get2Curves()
    #if rc is None: return

    #rgC_A_In = rc[0].Curve()
    #rgC_B_In = rc[1].Curve()


    bRebuildRails = Opts.values['bRebuildRails']
    bCPlane = Opts.values['bRebuildRails']
    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
    bSplitPathsAtKnots = Opts.values['bSplitPathsAtKnots']
    bAddArrayedLines = Opts.values['bAddArrayedLines']
    bAddBrep = Opts.values['bAddBrep']
    iBrepMethod = Opts.values['iBrepMethod']
    iLoftType = Opts.values['iLoftType']
    fBrepTol = Opts.values['fBrepTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rc = getOptions(rgCs_A_In, rgCs_B_In)

    if rc is None:
        # Cancel.
        return
    if not rc:
        # Create objects.
        return

    for lines in rc:
        for line in lines:
            sc.doc.Objects.AddCurve(rg.LineCurve(line))

    sc.doc.Views.Redraw()

    return

    bRebuildRails = Opts.values['bRebuildRails']
    bCPlane = Opts.values['bRebuildRails']
    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
    bSplitPathsAtKnots = Opts.values['bSplitPathsAtKnots']
    bAddArrayedLines = Opts.values['bAddArrayedLines']
    bAddBrep = Opts.values['bAddBrep']
    iBrepMethod = Opts.values['iBrepMethod']
    iLoftType = Opts.values['iLoftType']
    fBrepTol = Opts.values['fBrepTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    rgLine_ToArray = None
    rgCrvs0_Path = []
    rgCrvs1_Path = []

    while not (bAtGrevilles or bAtEqualDivisions):
        print "No path point sampling is enabled."
        sc.doc.Views.Redraw()
        rc = getOptions(bFirstGetObjects=False)
        if rc is None or rc is False:
            if rgLine_ToArray: rgLine_ToArray.Dispose()
            for c in rgCrvs0_Path: c.Dispose()
            for c in rgCrvs1_Path: c.Dispose()
            return
        objrefs_Paths = rc[0]

    rgCrvs0_Path = []
    for objref_Path in objrefs_Paths:
        c = objref_Path.Curve()
        rgCrvs0_Path.append(c)


    ncs_A = []
    ncs_B = []

    nc2_Path = ncA.ToNurbsCurve()
    if nc2_Path is None:
        print "NurbsCurve could not be calculated from curve."
        return
    #sc.doc.Objects.AddCurve(nc2_Path)




    Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")
    rc = prepareCurves(
            rgCrvs0=rgCrvs0_Path,
            bRebuild=bRebuildRails,
            bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
            bSplitPathsAtKnots=bSplitPathsAtKnots)
    if not rc: return
    rgCrvs1_Path = rc

    rgLine_ToArray = rg.Line(
            rg.Point3d(0.0, 0.0, 0.0),
            rg.Point3d(0.0, fDistance, 0.0))
    
    rgLineCrv_ToArray = rg.LineCurve(rgLine_ToArray)

    gLineCrvs1_Arrayed = []
    gCrvs_TaperEnds_All = []
    rgBreps1 = []
    gBreps1 = []

    rgCrvs_PathSegs = None
        
    Rhino.RhinoApp.SetCommandPrompt("Creating geometry ...")

    for rgCrv1_Path_1Seg in rgCrvs1_Path:
        rc = createCrossSectionLines_PerProjection(
            ncA=rgCrv1_Path_1Seg,
            rgObjs_ToArray=[rgLineCrv_ToArray],
            plane_Proj=plane_Proj,
            fTaper_Start_Deg=fTaper_Start_Deg,
            fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
            bTaperChangePerCrvParam=bTaperChangePerCrvParam,
            bAtGrevilles=bAtGrevilles,
            bAtEqualDivisions=bAtEqualDivisions,
            iDivisionCt=iDivisionCt,
            bDebug=bDebug)
        if rc is None: return
        # Flatten list.
        rgLineCrvs_Arrayed_1PathSeg = [rgLs[0] for rgLs in rc]
    
        if bAddArrayedLines:
            for i, rgL in enumerate(rgLineCrvs_Arrayed_1PathSeg):
                gL = sc.doc.Objects.AddCurve(rgL)
                if gL != Guid.Empty:
                    gLineCrvs1_Arrayed.append(gL)
                    #rgDot_ = rg.TextDot(str(i), rgObj1_Arrayed_1Seg.PointAtStart)
                    #sc.doc.Objects.AddTextDot(rgDot_)
        else:
            # bAddArrayedLines == False.
            if bAtEqualDivisions:
                if bAtEqualDivisions and iBrepMethod == 0:
                    s  = "AtEqualDivisions"
                    s += " only affects added arrayed lines"
                    s += " when BrepMethod == {},".format(sBrepMethods[iBrepMethod])
                    s += " but AddArrayed option is disabled."
                    print s

        if bAddBrep:

            if not bAtGrevilles or bAtEqualDivisions:
                rc = createCrossSectionLines_PerProjection(
                    ncA=rgCrv1_Path_1Seg,
                    rgObjs_ToArray=[rgLineCrv_ToArray],
                    plane_Proj=plane_Proj,
                    fTaper_Start_Deg=fTaper_Start_Deg,
                    fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                    bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                    bAtGrevilles=True,
                    bAtEqualDivisions=False,
                    iDivisionCt=0,
                    bDebug=bDebug)
                if rc is None: continue
                # Flatten list.
                rgLineCrvs_Arrayed_1PathSeg_GrevsOnly = [L[0] for L in rc]
            else:
                rgLineCrvs_Arrayed_1PathSeg_GrevsOnly = [L.Duplicate() for L in rgLineCrvs_Arrayed_1PathSeg]
                
            pts_EndOf_LineCrvs_Arrayed = []
            for rgLineCrv in rgLineCrvs_Arrayed_1PathSeg_GrevsOnly:
                pts_EndOf_LineCrvs_Arrayed.append(rgLineCrv.PointAtEnd)

                rgNurbsCrv_TaperEnd = rgCrv1_Path_1Seg.ToNurbsCurve()
                rgNurbsCrv_TaperEnd.SetGrevillePoints(pts_EndOf_LineCrvs_Arrayed)

            if bAddBrep:
                rc = createBrep(
                        iBrepMethod=iBrepMethod,
                        iLoftType=iLoftType,
                        fBrepTol=fBrepTol,
                        rgCrv_Path=rgCrv1_Path_1Seg,
                        rgNurbsCrv_TaperEnd_1Seg=rgNurbsCrv_TaperEnd,
                        rgLineCrvs_Arrayed=rgLineCrvs_Arrayed_1PathSeg_GrevsOnly)
                if rc is None:
                    print "Cannot create brep(s).  Check input."
                else:
                    rgBreps1.extend(rc)
                rgCrv1_Path_1Seg.Dispose()
                
            for c in rgLineCrvs_Arrayed_1PathSeg_GrevsOnly:
                c.Dispose()


    if rgBreps1:
        rgBreps1_Joined = rg.Brep.JoinBreps(
                rgBreps1,
                tolerance=0.5*sc.doc.ModelAbsoluteTolerance)
        for b in rgBreps1: b.Dispose()
            
        if rgBreps1_Joined:
            gBreps1 = []
            for rgBrep1_Joined in rgBreps1_Joined:
                gBrep1 = sc.doc.Objects.AddBrep(rgBrep1_Joined)
                rgBrep1_Joined.Dispose()
                if gBrep1 != Guid.Empty:
                    gBreps1.append(gBrep1)
            if bEcho:
                print "{} brep(s) with {} face(s) created.".format(
                len(gBreps1), len(rgBreps1))

    sc.doc.Views.Redraw()

    rc = getOptions(bFirstGetObjects=False)
    if rc is None:
        for g_ in gLineCrvs1_Arrayed:
            sc.doc.Objects.Delete(g_, True)
        for g_ in gCrvs_TaperEnds_All:
            sc.doc.Objects.Delete(g_, True)
        for g_ in gBreps1:
            sc.doc.Objects.Delete(g_, True)
    objrefs_Paths = rc[0]

    for g_ in gLineCrvs1_Arrayed:
        sc.doc.Objects.Delete(g_, True)
    for g_ in gCrvs_TaperEnds_All:
        sc.doc.Objects.Delete(g_, True)
    for g_ in gBreps1:
        sc.doc.Objects.Delete(g_, True)

    for c in rgCrvs0_Path: c.Dispose()
    for c in rgCrvs1_Path: c.Dispose()
    rgLineCrv_ToArray.Dispose()
    for c in rgLineCrvs_Arrayed_1PathSeg:
        c.Dispose()


if __name__ == '__main__': main()