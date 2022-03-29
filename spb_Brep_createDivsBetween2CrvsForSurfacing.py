"""
"""

"""
200402-04: Created, starting with other scripts.
200517: Further development.
220328: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from clr import StrongBox
from System import Array
from System import Enum
from System import Guid


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


    key = 'fDivisionLength'; keys.append(key)
    values[key] = 1000.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
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
    
    key = 'bAddCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
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


def prepareCurves(rgCs_In, bSplitPolyCrvToSegs=True, bSplitPathsAtKnots=False):
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
            rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()]
        else:
            # bSplitPolyCrvToSegs == True.
            if isinstance(rgCrv_Joined, rg.PolyCurve):
                rc = rgCrv_Joined.Explode()
                if rc:
                    rgCrvs_Exploded = rc
                    rgCrvs_SplitPoly = rc
            else:
                # not PolyCurve.
                rgCrvs_SplitPoly = [rgCrv_Joined.Duplicate()]

        if not bSplitPathsAtKnots:
            rgCrvs_Final.extend(rgCrvs_SplitPoly)
        else:
            for rgCrv_SplitPoly in rgCrvs_SplitPoly:
                ts_SpanBoundaries = [rgCrv_SplitPoly.Domain.T0]
                for iSpan in xrange(rgCrv_SplitPoly.SpanCount):
                    ts_SpanBoundaries.append(rgCrv_SplitPoly.SpanDomain(iSpan).T1)
                rc = rgCrv_SplitPoly.Split(ts_SpanBoundaries)
                if rc:
                    rgCrvs_SplitAtKnots = rc
                    rgCrvs_Final.extend(rgCrvs_SplitAtKnots)

                rgCrv_SplitPoly.Dispose()
    
    return rgCrvs_Final


def getParameters(nc, fDivisionLength=None):


    ts = []

    if fDivisionLength is None:
        fDivisionLength = 1000.0*sc.doc.ModelAbsoluteTolerance

    rc = nc.DivideByLength(
        segmentLength=fDivisionLength,
        includeEnds=True)
    if rc: ts.extend(rc)
    if not nc.IsClosed:
        # DivideByLength doesn't add the T1 segment
        # even though includeEnds == True.
        # https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_DivideByLength.htm
        # shows the parameter labeled as 'includeStart'.
        ts.append(nc.Domain.T1)

        lastDivLength = nc.GetLength(subdomain=rg.Interval(ts[-2], ts[-1]))
        if lastDivLength < 0.5 * fDivisionLength:
            del ts[-2]

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


def getPointPairsBetweenCrvs_ClosestPts(rgCrvA, rgCrvB, bDebug=False):
    """
    """

    def getPoints(nc_From, nc_To):

        pts_From = []
        ts_From = []
        pts_To = []
        ts_To = []

        ts_Grevilles = nc_From.GrevilleParameters()
        pts_Grevilles = nc_From.GrevillePoints()

        ts_Knots = set([nc_From.Knots[i] for i in range(nc_From.Knots.Count)])
        pts_Knots = [nc_From.PointAt(t) for t in ts_Knots]

        ts_From_IntvEnds = list(set(list(ts_Grevilles) + list(ts_Knots)))

        ts_From_IntvEnds.sort()

        for i in range(1, len(ts_From_IntvEnds)):
            interval = rg.Interval(ts_From_IntvEnds[i-1], ts_From_IntvEnds[i])

            subC_From = nc_From.Trim(interval)

            rc = subC_From.ClosestPoints(nc_To)
            if not rc[0]: continue

            pt_From = rc[1]
            b, t_From = nc_From.ClosestPoint(pt_From)
            if not b: continue

            # Do not accept points at each interval's T0 and T1.
            # TODO: Is this always correct to do?
            for pt_From_IntvEnd in list(pts_Grevilles) + pts_Knots:
                dist = pt_From_IntvEnd.DistanceTo(pt_From)
                if dist <= sc.doc.ModelAbsoluteTolerance:
                    break # to next interval.

            #for t_From_IntvEnd in ts_From_IntvEnds:
            #    param_dist = abs(t_From_IntvEnd-t_From)
            #    if param_dist <= 1.0/(2.0**32):
            #        break # to next interval.
            else:
                pt_To = rc[2]
                if (
                    pt_To.EpsilonEquals(nc_To.PointAtStart, sc.doc.ModelAbsoluteTolerance) or
                    pt_To.EpsilonEquals(nc_To.PointAtEnd, sc.doc.ModelAbsoluteTolerance)
                ):
                    continue # to next interval.

                b, t_To = nc_To.ClosestPoint(pt_To)
                if not b: continue

                pts_From.append(pt_From)
                ts_From.append(t_From)
                pts_To.append(pt_To)
                ts_To.append(t_To)
                #sc.doc.Objects.AddPoint(pt_From)
                #sc.doc.Objects.AddPoint(pt_To)

        return pts_From, ts_From, pts_To, ts_To


    def areEpsilonEqual(a, b, epsilon):
        # This is a relative comparison.
        delta = abs(a - b)
        fRelComp = delta / max(abs(a), abs(b))
        bRhEpsEquals = Rhino.RhinoMath.EpsilonEquals(a, b, epsilon)
        return fRelComp < epsilon


    ts_Out_A = []
    ts_Out_B = []

    nc_A = rgCrvA.ToNurbsCurve()
    nc_B = rgCrvB.ToNurbsCurve()


    pts_Out_A, ts_Out_A, pts_Out_B, ts_Out_B = getPoints(nc_A, nc_B)

    pts_B_2ndRound, ts_B_2ndRound, pts_A_2ndRound, ts_A_2ndRound = getPoints(nc_B, nc_A)


    # Add parameters from 2nd round to return data if not redundant.
    pts_A_toAdd = []
    ts_A_toAdd = []
    for i in range(len(pts_A_2ndRound)):
        pt_toCheck = pts_A_2ndRound[i]
        for pt in pts_Out_A:
            if pt.EpsilonEquals(pt_toCheck, sc.doc.ModelAbsoluteTolerance):
                break
        else:
            pts_A_toAdd.append(pt_toCheck)
            ts_A_toAdd.append(ts_A_2ndRound[i])

    pts_Out_A.extend(pts_A_toAdd) 
    ts_Out_A.extend(ts_A_toAdd) 

    pts_B_toAdd = []
    ts_B_toAdd = []
    for i in range(len(pts_B_2ndRound)):
        pt_toCheck = pts_B_2ndRound[i]
        for pt in pts_Out_B:
            if pt.EpsilonEquals(pt_toCheck, sc.doc.ModelAbsoluteTolerance):
                break
        else:
            pts_B_toAdd.append(pt_toCheck)
            ts_B_toAdd.append(ts_B_2ndRound[i])

    pts_Out_B.extend(pts_B_toAdd) 
    ts_Out_B.extend(ts_B_toAdd) 

    return pts_Out_A, ts_Out_A, pts_Out_B, ts_Out_B


def createCrossSectionLines_PerpTo(ncs_A, ncs_B, ts_perA, bDebug=False):
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

            pt_A = ncA.PointAt(tA)

            for iB in range(len(ncs_B)):
                
                nc_B = ncs_B[iB]

                res, tB_forSeed = nc_B.ClosestPoint(pt_A)
                if not res:
                    raise ValueError("ClosestPoint returned {}.".format(res))

                bSuccess, tB = nc_B.GetLocalPerpPoint(
                    testPoint=pt_A,
                    seedParmameter=tB_forSeed)

                if not bSuccess:
                    continue

                line = rg.Line(ncA.PointAt(tA), to=nc_B.PointAt(tB))
                if bDebug: sc.doc.Objects.AddLine(line)

                lines_Out_perA[-1].append(line)

    return lines_Out_perA


def createCrossSectionLines_EqualDivs(rgCrv_DivByLength, rgCrv_DivByCount, fDivisionLength=None, bDebug=False):
    """
    """

    lines_Out = []

    if fDivisionLength is None:
        fDivisionLength = 1000.0*sc.doc.ModelAbsoluteTolerance

    tAs = getParameters(rgCrv_DivByLength, fDivisionLength=fDivisionLength)

    strongBox_points = StrongBox[Array[rg.Point3d]]()

    tBs = rgCrv_DivByCount.DivideByCount(
        segmentCount=len(tAs)-1,
        includeEnds=True,
        points=strongBox_points)

    pts_B = list(strongBox_points.Value)

    for iT in range(len(tAs)):
        tA = tAs[iT]

        pt_A = rgCrv_DivByLength.PointAt(tA)
        pt_B = pts_B[iT]

        line = rg.Line(pt_A, to=pt_B)
        if bDebug: sc.doc.Objects.AddLine(line)

        lines_Out.append(line)

    return lines_Out


def createConnectingCurves(rgCs_A_In, rgCs_B_In, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bSplitPolyCrvToSegs = getOpt('bSplitPolyCrvToSegs')
    bSplitPathsAtKnots = getOpt('bSplitPathsAtKnots')
    fDivisionLength = getOpt('fDivisionLength')
    bAtGrevilles = getOpt('bAtGrevilles')
    iDivisionCt = getOpt('iDivisionCt')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')



    Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")

    rgCs_A_Prepped = prepareCurves(
        rgCs_A_In,
        bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
        bSplitPathsAtKnots=bSplitPathsAtKnots)

    rgCs_B_Prepped = prepareCurves(
        rgCs_B_In,
        bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
        bSplitPathsAtKnots=bSplitPathsAtKnots)


    joinedCrvs_A = rg.Curve.JoinCurves(rgCs_A_Prepped)
    if len(joinedCrvs_A) != 1:
        return

    joinedCrvs_B = rg.Curve.JoinCurves(rgCs_B_Prepped)
    if len(joinedCrvs_B) != 1:
        return

    joined_A = joinedCrvs_A[0]
    joined_B = joinedCrvs_B[0]

    if not rg.Curve.DoDirectionsMatch(joined_A, joined_B):
        joined_B.Reverse()

    rc = getPointPairsBetweenCrvs_ClosestPts(
        joined_A,
        joined_B,
        bDebug=bDebug)

    pts_A, ts_A, pts_B, ts_B = rc

    for pt in pts_A:
        sc.doc.Objects.AddPoint(pt)

    for pt in pts_B:
        sc.doc.Objects.AddPoint(pt)


    if len(ts_A) != len(ts_B):
        raise ValueError("Parameter count mismatch: {} and {}".format(
            len(ts_A), len(ts_B)))


    # Split curves at parameters.
    splitCrvs_A = joined_A.Split(ts_A)
    #for c in splitCrvs_A: sc.doc.Objects.AddCurve(c)

    splitCrvs_B = joined_B.Split(ts_B)

    lines_Out = []

    for i in range(len(splitCrvs_A)):
        splitA = splitCrvs_A[i]
        splitB = splitCrvs_B[i]
        
        lengthA = splitA.GetLength()
        lengthB = splitB.GetLength()

        if lengthA > lengthB:
            # Divide A.  Divisions of B will be less than fDivisionLength.
            crv_DivByLength = splitA
            ts_DivByLength = ts_A
            crv_DivByCount = splitB
            bReverseCrvs = False
        else:
            crv_DivByLength = splitB
            ts_DivByLength = ts_B
            crv_DivByCount = splitA
            bReverseCrvs = True

        lines = createCrossSectionLines_EqualDivs(
            crv_DivByLength,
            crv_DivByCount,
            fDivisionLength,
            bDebug=bDebug)

        if bReverseCrvs:
            for line in lines:
               line.Flip()

        if i == 0:
            lines_Out.append(lines)
        else:
            lines_Out.append(lines[1:])

    return lines_Out


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


def getInput_Options(rgCs_A_In, rgCs_B_In):
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

    rgCrvs_Connecting = createConnectingCurves(
        rgCs_A_In,
        rgCs_B_In,
        fDivisionLength=Opts.values['fDivisionLength'],
        )

    conduit = DrawLinesConduit()

    if rgCrvs_Connecting:
        conduit.lines = [line for lines in rgCrvs_Connecting for line in lines]
        conduit.Enabled = True
        sc.doc.Views.Redraw()


    idxs_Opts = {}

    while True:
        sc.escape_test()

        for key in Opts.keys:
            valuesBefore[key] = Opts.values[key]

        go.AddOptionDouble(Opts.names['fDivisionLength'], Opts.riOpts['fDivisionLength'])
        go.AddOptionToggle(Opts.names['bAtGrevilles'], Opts.riOpts['bAtGrevilles'])
        #go.AddOptionToggle(Opts.names['bAtEqualDivisions'],
        #                    Opts.riOpts['bAtEqualDivisions'])
        #if Opts.values['bAtEqualDivisions']:
        #    go.AddOptionInteger(Opts.names['iDivisionCt'],
        #                        Opts.riOpts['iDivisionCt'])
        go.AddOptionToggle(Opts.names['bSplitPolyCrvToSegs'], Opts.riOpts['bSplitPolyCrvToSegs'])
        go.AddOptionToggle(Opts.names['bSplitPathsAtKnots'], Opts.riOpts['bSplitPathsAtKnots'])
        go.AddOptionToggle(Opts.names['bAddCrvs'], Opts.riOpts['bAddCrvs'])
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
            return rgCrvs_Connecting

        Opts.setValues()
        Opts.saveSticky()

        for key in Opts.keys:
            if valuesBefore[key] != Opts.values[key]:
                rgCrvs_Connecting = createConnectingCurves(
                    rgCs_A_In,
                    rgCs_B_In,
                    fDivisionLength=Opts.values['fDivisionLength'],
                    )
                if rgCrvs_Connecting:
                    conduit.lines = [line for lines in rgCrvs_Connecting for line in lines]
                    conduit.Enabled = True
                    sc.doc.Views.Redraw()
                break
        else:
            for key in Opts.keys:
                print valuesBefore[key], Opts.values[key]
            print "No options were changed."


def main():

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select first set of rail curves",
        acceptNothing=False,
        filter=rd.ObjectType.Curve)
    if res == Rhino.Commands.Result.Cancel: return
    
    rgCs_A_In = [o.Curve() for o in objrefs]

    sc.doc.Objects.UnselectAll()

    res, objrefs = ri.RhinoGet.GetMultipleObjects(
        "Select second set of rail curves",
        acceptNothing=False,
        filter=rd.ObjectType.Curve)
    if res == Rhino.Commands.Result.Cancel: return

    rgCs_B_In = [o.Curve() for o in objrefs]

    sc.doc.Objects.UnselectAll()

    bAtGrevilles = Opts.values['bAtGrevilles']
    bAtEqualDivisions = Opts.values['bAtEqualDivisions']
    iDivisionCt = Opts.values['iDivisionCt']
    bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
    bSplitPathsAtKnots = Opts.values['bSplitPathsAtKnots']
    bAddCrvs = Opts.values['bAddCrvs']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    rc = getInput_Options(rgCs_A_In, rgCs_B_In)

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


if __name__ == '__main__': main()