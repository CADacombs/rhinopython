"""

Notes on UI's _Match:
A unitized tangent vector, as obtained with Curve.CurvatureAt, seems to be used for
Tangent and Curvature matches (at least with Average disabled)
For Tangent matches, NurbsCurve.SetEndCondition results are identical.
For Curvature matches, identical results were not found but were closer when using
the unitized tangent vector than when using the 1st derivative vector of Curve B directly.
"""

"""
190708-25: Created.
190824: Now, a unit vector is used when target continuity is G1.  Otherwise, SetEndCondition may move the control point quite far from the closest possible position.
        Added curve deviation feedback.
190831: Now tests whether joint is already at desired continuity before calculating new curve(s).
190907: Refactored.  Added bLimitCrvDev and fDevTol.
190908-12: More refactoring and optimization of curves.
190913: Added sContinuity in place of iContinuity in most functions.  Various debugging and refactoring.
190918: createForCurveA: Now curves with dev of None will not pass.
        createForBothCurves is functional for position and tangency continuities.
190919: Bug fix.
190924, 1010, 1020: Import-related update.
191208: Printed output bug fix.
200112-16: Updated Opts to latest code format.  Import-related updates.
200120: Refactored, nesting some functions for simpler outside code calling this module.
200505: Import-related updates.
200630: Now correctly compares control point positions between curves of different degrees for C1 matching.
        Improved printed feedback when curve is not modified.
201122, 210908: Import-related updates.
211029: Simplified curve to modify option from 3 to 2 choices.  Added fAngleTol_Deg.
220328, 0425: Import-related update.

TODO:
Try to adjust for G1 to same G0G1 length as original curve.  Check its deviation against SetEndCondition result.
Determine whether to keep C1 continuity (and add C2).
    If so, one scenario is to set all curve domains to [0.0,1.0] before matching.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum
from System import Guid
from System.Drawing import Color

import spb_Crv_continuityBetween2
import xCurve_inflections
import xCurve_radiusMinima
import spb_Crv_fitRebuild
import xNurbsCurve_maximizeMinimumRadius


sOpts_Continuity = ['G0', 'G1', 'C1', 'G2']


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


    key = 'iContinuity'; keys.append(key)
    values[key] = 1
    riAddOpts[key] = addOptionList(key, names, sOpts_Continuity, values)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fAngleTol_Deg'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSkipIfAlreadyAtContinuity'; keys.append(key)
    values[key] = True
    names[key] = 'IfAlreadyAtContinuity'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Process', 'Skip')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iPreserveOtherEnd'; keys.append(key)
    values[key] = 2
    names[key] = 'PreserveOtherEnd'
    riAddOpts[key] = addOptionList(key, names,
        Enum.GetNames(rg.NurbsCurve.NurbsCurveEndConditionType),
        values)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLimitCrvDev'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMaximizeMinRadius'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bModOnly1stPicked'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDeformLines'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRebuildRationals'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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


def addGeoms(geoms, bRedraw=True):
    """ For debugging. """

    if not hasattr(geoms, '__iter__'):
        geoms = [geoms]

    for geom in geoms:
        if isinstance(geom, tuple):
            gOut = sc.doc.Objects.AddSurface(geom[1])
        elif isinstance(geom, rg.Curve):
            gOut = sc.doc.Objects.AddCurve(geom)
        elif isinstance(geom, rg.Surface):
            gOut = sc.doc.Objects.AddSurface(geom)
        elif isinstance(geom, rg.Point3d):
            gOut = sc.doc.Objects.AddPoint(geom)
        elif isinstance(geom, rg.Plane):
            intrvl = rg.Interval(-sc.doc.ModelAbsoluteTolerance*1000.0, sc.doc.ModelAbsoluteTolerance*1000.0)
            psrf = rg.PlaneSurface(geom, intrvl, intrvl)
            gOut = sc.doc.Objects.AddSurface(psrf)
        else:
            raise ValueError("Method to add {} missing from addGeoms.".format(geom))
    if bRedraw: sc.doc.Views.Redraw()
    return gOut


def getInput():
    """
    Get curve and parameter with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves near their ends to match")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.OneByOnePostSelect = True
    go.DisablePreSelect()

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    while True:
        idxs_Opts['iContinuity'] = Opts.riAddOpts['iContinuity'](go)
        idxs_Opts['fAngleTol_Deg'] = Opts.riAddOpts['fAngleTol_Deg'](go)
        Opts.riAddOpts['bSkipIfAlreadyAtContinuity'](go)
        idxs_Opts['iPreserveOtherEnd'] = Opts.riAddOpts['iPreserveOtherEnd'](go)
        Opts.riAddOpts['bLimitCrvDev'](go)
        if Opts.values['bLimitCrvDev']: Opts.riAddOpts['fDevTol'](go)
        if sOpts_Continuity[Opts.values['iContinuity']] in ('G0', 'G1'):
            Opts.riAddOpts['bMaximizeMinRadius'](go)
        Opts.riAddOpts['bModOnly1stPicked'](go)
        Opts.riAddOpts['bDeformLines'](go)
        Opts.riAddOpts['bRebuildRationals'](go)
        Opts.riAddOpts['bReplace'](go)
        Opts.riAddOpts['bEcho'](go)
        Opts.riAddOpts['bDebug'](go)

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])

        # An option was selected or a number was entered.

        key = 'fAngleTol_Deg'
        if Opts.riOpts[key].CurrentValue <= 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        if res == ri.GetResult.Number:
            Opts.riOpts['fDevTol'].CurrentValue = go.Number()
        elif res == ri.GetResult.Option:
            if go.Option().Index == idxs_Opts['iContinuity']:
                Opts.values['iContinuity'] = (
                    go.Option().CurrentListOptionIndex)
            elif go.Option().Index == idxs_Opts['iPreserveOtherEnd']:
                Opts.values['iPreserveOtherEnd'] = (
                    go.Option().CurrentListOptionIndex)

        if Opts.riOpts['bLimitCrvDev'].CurrentValue:
            key = 'fDevTol'
            if Opts.riOpts[key].CurrentValue <= 0.0:
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
    if fDistance is None:
        return "(None)"
    if fDistance == Rhino.RhinoMath.UnsetValue:
        return "(Infinite)"
    if fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-3)):
        return "{:.2e}".format(fDistance)
    return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def createNurbsCurves(rgCurveA, rgCurveB, bT1WorkEnd_A, bT1WorkEnd_B, bModifyA, bModifyB, **kwargs):
    """
    Parameters:
        rgCurveA
        rgCurveB
        bT1WorkEnd_A
        bT1WorkEnd_B
        bModifyA: If both bModifyA and bModifyB are False, will return the curve with less deviation result.
        bModifyB: (See bModifyA.)

    Returns:
        Success (One of the following):
            rg.NurbsCurve, rg.NurbsCurve
            rg.NurbsCurve, None
            None, rg.NurbsCurve
        Fail: None
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if 'sContinuity' in kwargs:
        sContinuity = kwargs['sContinuity']
    else:
        sContinuity = sOpts_Continuity[Opts.values['iContinuity']]
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bSkipIfAlreadyAtContinuity = getOpt('bSkipIfAlreadyAtContinuity')
    iPreserveOtherEnd = getOpt('iPreserveOtherEnd')
    fDevTol = getOpt('fDevTol')
    bMaximizeMinRadius = getOpt('bMaximizeMinRadius')
    bDeformLines = getOpt('bDeformLines')
    bRebuildRationals = getOpt('bRebuildRationals')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if not bModifyA and not bModifyB:
        print "Neither curve is supposed to be modified."
        return


    def areCurvesAlreadyAtDesiredContinuity():
        """
        """

        nc0_A = rgCurveA.ToNurbsCurve()
        if bT1WorkEnd_A:
            t_WorkEnd_A = nc0_A.Domain.T1
            crvEvalSide_A = rg.CurveEvaluationSide.Below
            idxCp_Pos_A = nc0_A.Points.Count - 1
            idxCp_Tan_A = nc0_A.Points.Count - 2
        else:
            t_WorkEnd_A = nc0_A.Domain.T0
            crvEvalSide_A = rg.CurveEvaluationSide.Above
            idxCp_Pos_A = 0
            idxCp_Tan_A = 1
    
        nc0_B = rgCurveB.ToNurbsCurve()
        if bT1WorkEnd_B:
            t_WorkEnd_B = nc0_B.Domain.T1
            crvEvalSide_B = rg.CurveEvaluationSide.Below
            idxCp_Pos_B = nc0_B.Points.Count - 1
            idxCp_Tan_B = nc0_B.Points.Count - 2
        else:
            t_WorkEnd_B = nc0_B.Domain.T0
            crvEvalSide_B = rg.CurveEvaluationSide.Above
            idxCp_Pos_B = 0
            idxCp_Tan_B = 1

        rc = spb_Crv_continuityBetween2.processCurves(
                nc0_A,
                nc0_B,
                bEvalT1End_A=bT1WorkEnd_A,
                bEvalT1End_B=bT1WorkEnd_B,
                fDistTol=Rhino.RhinoMath.ZeroTolerance,
                fAngleTol_Vector_Deg=fAngleTol_Deg,
                bAlignCrvDirs=True,
                bDebug=bDebug,
        )
        if not rc:
            print "Error in obtaining data from spb_Crv_continuityBetween2.processCurves."
            return
        iContinuity0_G, iContinuity0_C, sContinuityDescr = rc
        if bDebug:
            print sContinuityDescr
            sEval='sContinuity'; print sEval+': ',eval(sEval)
        s  = "Starting continuity: G{}/C{}  ".format(iContinuity0_G, iContinuity0_C)
        if sContinuity == 'G0':
            if iContinuity0_G:
                if bEcho:
                    s += "  Curves' endpoints already meet at G{} continuity.".format(
                            iContinuity0_G)
                    print s
                nc0_A.Dispose()
                nc0_B.Dispose()
                return True
        elif sContinuity == 'G1':
            if iContinuity0_G >= 1:
                if bEcho:
                    s += "  Curves' endpoints already meet at G{} continuity.".format(
                            iContinuity0_G)
                    print s
                nc0_A.Dispose()
                nc0_B.Dispose()
                return True
        elif sContinuity == 'C1':
            if iContinuity0_C >= 1:
                if bEcho:
                    s += "  Curves' endpoints already meet at G{} continuity.".format(
                            iContinuity0_G)
                    print s
                nc0_A.Dispose()
                nc0_B.Dispose()
                return True
        elif sContinuity == 'G2':
            if iContinuity0_G >= 2:
                if bEcho:
                    s += " Curves' endpoints already meet at G{} continuity.".format(
                            iContinuity0_G)
                    print s
                nc0_A.Dispose()
                nc0_B.Dispose()
                return True
        if bEcho: print s

        return False


    def rebuildRationalCurve(crv):
        if bRebuildRationals and crv.IsRational:
            degree = 3 if crv.Degree < 3 else crv.Degree

            if not fDevTol:
                # Create Bezier.
                rc = crv.Rebuild(
                    pointCount=degree+1,
                    degree=degree,
                    preserveTangents=True)
                if rc:
                    return rc
            else:
                rc = spb_Crv_fitRebuild.rebuildCurve(
                        crv,
                        fDevTol=0.5*fDevTol,
                        iDegree=degree,
                        bPreserveEndTans=True,
                        bFurtherTranslateCps=True,
                        iMinCpCt=None,
                        iMaxCpCt=40,
                        bDebug=False,
                        )
                if rc[0] is not None:
                    return rc[0]


    def doCrvsMatchWithinDistDev(cA, cB, fDevTol=fDevTol):
        if fDevTol:
            fDev = getMaximumDeviation(cA, cB)
            if fDev is None:
                print "Curve is not acceptable because its deviation cannot be determined."
                return False
            elif fDev > fDevTol:
                print "Curve is not acceptable because deviation would be {}.".format(
                        formatDistance(fDev))
                return False
        return True


    def createForBothCurves(rgCurveA, rgCurveB, bT1WorkEnd_A, bT1WorkEnd_B):
        """
        Returns on success: (rg.NurbsCurve, rg.NurbsCurve)
        Returns on fail: None
        """

        c0_A = rgCurveA
        c0_B = rgCurveB
        c0s = [rgCurveA, rgCurveB]
        bT1WorkEnds = bT1WorkEnd_A, bT1WorkEnd_B

        ncs_PreMatch = [c0_A.ToNurbsCurve(), c0_B.ToNurbsCurve()]

        if bRebuildRationals:
            for i in 0,1:
                rc = rebuildRationalCurve(ncs_PreMatch[i])
                if rc:
                    ncs_PreMatch[i].Dispose()
                    ncs_PreMatch[i] = rc

        nc_A_PreMatch = ncs_PreMatch[0]
        if bT1WorkEnd_A:
            t_WorkEnd_A = nc_A_PreMatch.Domain.T1
            crvEvalSide_A = rg.CurveEvaluationSide.Below
            idxCp_Pos_A = nc_A_PreMatch.Points.Count - 1
            idxCp_Tan_A = nc_A_PreMatch.Points.Count - 2
            idxCp_Crv_A = nc_A_PreMatch.Points.Count - 3
        else:
            t_WorkEnd_A = nc_A_PreMatch.Domain.T0
            crvEvalSide_A = rg.CurveEvaluationSide.Above
            idxCp_Pos_A = 0
            idxCp_Tan_A = 1
            idxCp_Crv_A = 2

        nc_B_PreMatch = ncs_PreMatch[1]
        if bT1WorkEnd_B:
            t_WorkEnd_B = nc_B_PreMatch.Domain.T1
            crvEvalSide_B = rg.CurveEvaluationSide.Below
            idxCp_Pos_B = nc_B_PreMatch.Points.Count - 1
            idxCp_Tan_B = nc_B_PreMatch.Points.Count - 2
            idxCp_Crv_B = nc_B_PreMatch.Points.Count - 3
        else:
            t_WorkEnd_B = nc_B_PreMatch.Domain.T0
            crvEvalSide_B = rg.CurveEvaluationSide.Above
            idxCp_Pos_B = 0
            idxCp_Tan_B = 1
            idxCp_Crv_B = 2


        #vectsDerivatives_A = nc0_A.DerivativeAt(
        #        t_WorkEnd_A,
        #        derivativeCount=2,
        #        side=crvEvalSide_A)
        #print vectsDerivatives_A[2], nc0_A.CurvatureAt(t_WorkEnd_A)

        #vectsDerivatives_B = nc0_B.DerivativeAt(
        #        t_WorkEnd_B,
        #        derivativeCount=2,
        #        side=crvEvalSide_B)
        #print vectsDerivatives_B[2], nc0_B.CurvatureAt(t_WorkEnd_B)


        # (PC to maintain condition type at opposite end of curve) + (PC to modify working side of curve), where PC is "required min. point count".
        if sContinuity == 'G0':
            iCt_Cp_MinNeeded = iPreserveOtherEnd + 1
        elif sContinuity in ('G1', 'C1'):
            iCt_Cp_MinNeeded = iPreserveOtherEnd + 2
        elif sContinuity == 'G2':
            iCt_Cp_MinNeeded = iPreserveOtherEnd + 3

        nc_WIP_A = None
        nc_WIP_B = None


        if bDeformLines:
            for i in 0,1:
                if ncs_PreMatch[i].Points.Count < iCt_Cp_MinNeeded:
                    if ncs_PreMatch[i].Points.Count == 2:
                        ts = ncs_PreMatch[i].DivideByCount(
                            segmentCount=iCt_Cp_MinNeeded-1,
                            includeEnds=True)
                        rc = rg.NurbsCurve.Create(
                            periodic=False,
                            degree=3,
                            points=[ncs_PreMatch[i].PointAt(t) for t in ts])
                        ncs_PreMatch[i].Dispose()
                        ncs_PreMatch[i] = rc
            nc_A_PreMatch = ncs_PreMatch[0]
            nc_B_PreMatch = ncs_PreMatch[1]


        if (
                nc_A_PreMatch.Points.Count < iCt_Cp_MinNeeded
                and
                nc_B_PreMatch.Points.Count < iCt_Cp_MinNeeded
        ):
            print "Neither curve has enough control points for this continuity modification."
            for i in 0,1:
                if ncs_PreMatch[i] is not None: ncs_PreMatch[i].Dispose()
            return
        elif nc_A_PreMatch.Points.Count < iCt_Cp_MinNeeded:
            print "Curve A doesn't have enough control points for this continuity modification."
            for i in 0,1:
                if ncs_PreMatch[i] is not None: ncs_PreMatch[i].Dispose()
            return
        elif nc_B_PreMatch.Points.Count < iCt_Cp_MinNeeded:
            print "Curve B doesn't have enough control points for this continuity modification."
            for i in 0,1:
                if ncs_PreMatch[i] is not None: ncs_PreMatch[i].Dispose()
            return

        nc_WIP_A = nc_A_PreMatch.Duplicate()
        nc_WIP_B = nc_B_PreMatch.Duplicate()


        # First, modify position.
        pt_Average = (0.5 * nc_WIP_A.Points[idxCp_Pos_A].Location +
                      0.5 * nc_WIP_B.Points[idxCp_Pos_B].Location)
        if nc_WIP_A.Points[idxCp_Pos_A] != pt_Average:
            nc_WIP_A.Points[idxCp_Pos_A] = pt_Average
            nc_WIP_B.Points[idxCp_Pos_B] = pt_Average


        nc_WIPs = [nc_WIP_A, nc_WIP_B]

        if not (
            doCrvsMatchWithinDistDev(nc_WIPs[0], ncs_PreMatch[0]) and
            doCrvsMatchWithinDistDev(nc_WIPs[1], ncs_PreMatch[1])
        ):
            ncs_PreMatch[0].Dispose()
            ncs_PreMatch[1].Dispose()
            nc_WIPs[0].Dispose()
            nc_WIPs[1].Dispose()
            return


        # Determine tangent vector(s).
        vectTan_A1 = nc_A_PreMatch.TangentAt(t_WorkEnd_A)
        vectTan_B1 = nc_B_PreMatch.TangentAt(t_WorkEnd_B)
        if bT1WorkEnd_A == bT1WorkEnd_B:
            vectTan_Target_B = (0.5 * -vectTan_A1) + (0.5 * vectTan_B1)
            vectTan_Target_A = -vectTan_Target_B
        else:
            vectTan_Target_A = (0.5 * vectTan_A1) + (0.5 * vectTan_B1)
            vectTan_Target_B = vectTan_Target_A


        # Modify tangency.
        if sContinuity in ('G1', 'C1', 'G2'):
            # for G1 or (C1 and the adjacent curve is linear or circular-arc-shaped).
            bEndConditionSet_A = nc_WIP_A.SetEndCondition(
                    bSetEnd=bT1WorkEnd_A,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
                    point=pt_Average,
                    tangent=vectTan_Target_A)
            bEndConditionSet_B = nc_WIP_B.SetEndCondition(
                    bSetEnd=bT1WorkEnd_B,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
                    point=pt_Average,
                    tangent=vectTan_Target_B)

            # Deviation check.
            if not (
                doCrvsMatchWithinDistDev(nc_WIPs[0], ncs_PreMatch[0]) and
                doCrvsMatchWithinDistDev(nc_WIPs[1], ncs_PreMatch[1])
            ):
                #addGeoms(nc_WIPs); 1/0
                ncs_PreMatch[0].Dispose()
                ncs_PreMatch[1].Dispose()
                nc_WIPs[0].Dispose()
                nc_WIPs[1].Dispose()
                return


            if sContinuity == 'C1' and not nc_A_PreMatch.CurvatureAt(t_WorkEnd_A).IsTiny() and not nc_B_PreMatch.CurvatureAt(t_WorkEnd_B).IsTiny():
                # TODO: Make average.
                # In the following, some of the paretheses are included to obtain Point3d instead of Vector3d.
            
                pt_Tan1_A = 0.5 * nc_A_PreMatch.Points[idxCp_Tan_A].Location + 0.5 * (
                        nc_B_PreMatch.Points[idxCp_Pos_B].Location -
                        float(nc_B_PreMatch.Degree / nc_WIP_A.Degree) *
                        (nc_B_PreMatch.Points[idxCp_Tan_B].Location -
                        nc_B_PreMatch.Points[idxCp_Pos_B].Location))
                #sc.doc.Objects.AddPoint(pt_Tan1_A)
            
                pt_Tan1_B = 0.5 * nc_B_PreMatch.Points[idxCp_Tan_B].Location + 0.5 * (
                        nc_A_PreMatch.Points[idxCp_Pos_A].Location -
                        float(nc_A_PreMatch.Degree / nc_WIP_B.Degree) *
                        (nc_A_PreMatch.Points[idxCp_Tan_A].Location -
                        nc_A_PreMatch.Points[idxCp_Pos_A].Location))
                #sc.doc.Objects.AddPoint(pt_Tan1_B)
            
                if nc_WIP_A.Points[idxCp_Tan_A] != pt_Tan1_A:
                    nc_WIP_A.Points[idxCp_Tan_A] = pt_Tan1_A
                if nc_WIP_B.Points[idxCp_Tan_B] != pt_Tan1_B:
                    nc_WIP_B.Points[idxCp_Tan_B] = pt_Tan1_B

                # Deviation check.
                if not (
                    doCrvsMatchWithinDistDev(nc_WIPs[0], ncs_PreMatch[0]) and
                    doCrvsMatchWithinDistDev(nc_WIPs[1], ncs_PreMatch[1])
                ):
                    ncs_PreMatch[0].Dispose()
                    ncs_PreMatch[1].Dispose()
                    nc_WIPs[0].Dispose()
                    nc_WIPs[1].Dispose()
                    return


            if bMaximizeMinRadius and sContinuity not in ('C1', 'C2', 'G2'):
                s = "Adjusting tangent control point spread ..."
                if bEcho: print s
                nc_WIPs = [nc_WIP_A.ToNurbsCurve(), nc_WIP_B.ToNurbsCurve()]
                for i in 0,1:
                    rc = xNurbsCurve_maximizeMinimumRadius.adjustTanCpSpread_OneEndOnly(
                            nc_WIPs[i],
                            bT1WorkEnds[i],
                            fDevTol,
                            ncs_PreMatch[i],
                            bDebug=bDebug,
                    )
                    if rc:
                        nc_WIPs[i].Dispose()
                        nc_WIPs[i] = rc[0]
                nc_WIP_A, nc_WIP_B = nc_WIPs

                # Deviation has already been checked in adjustTanCpSpread.


        if sContinuity == 'G2':
            # Determine curvature vector(s).
            vectCrv_A1 = nc_A_PreMatch.CurvatureAt(t_WorkEnd_A)
            vectCrv_B1 = nc_B_PreMatch.CurvatureAt(t_WorkEnd_B)
            if bT1WorkEnd_A == bT1WorkEnd_B:
                vectCrv_Target_B = (0.5 * -vectCrv_A1) + (0.5 * vectCrv_B1)
                vectCrv_Target_A = -vectCrv_Target_B
            else:
                vectCrv_Target_A = vectCrv_Target_B = (0.5 * vectCrv_A1) + (0.5 * vectCrv_B1)


            bEndConditionSet_A = nc_WIP_A.SetEndCondition(
                    bSetEnd=bT1WorkEnd_A,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                    point=pt_Average,
                    tangent=vectTan_Target_A,
                    curvature=vectCrv_Target_A)
            bEndConditionSet_B = nc_WIP_B.SetEndCondition(
                    bSetEnd=bT1WorkEnd_B,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                    point=pt_Average,
                    tangent=vectTan_Target_B,
                    curvature=vectCrv_Target_B)

        # Deviation check.
        if not (
            doCrvsMatchWithinDistDev(nc_WIPs[0], ncs_PreMatch[0]) and
            doCrvsMatchWithinDistDev(nc_WIPs[1], ncs_PreMatch[1])
        ):
            ncs_PreMatch[0].Dispose()
            ncs_PreMatch[1].Dispose()
            nc_WIPs[0].Dispose()
            nc_WIPs[1].Dispose()
            return

        # Curve is within tolerance.

        for i in 0,1:
            if nc_WIPs[i].EpsilonEquals(ncs_PreMatch[i], 1e-12):
                print "Curve was not modified."
                ncs_PreMatch[0].Dispose()
                ncs_PreMatch[1].Dispose()
                nc_WIPs[0].Dispose()
                nc_WIPs[1].Dispose()
                return

        ncs_PreMatch[0].Dispose()
        ncs_PreMatch[1].Dispose()

        return tuple(nc_WIPs)


    def createForOneCurve(rgCurve_One, rgCurve_Ref, bT1WorkEnd_One, bT1WorkEnd_Ref):
        """
        Parameters:
            rgCurve_One: For which to create a new curve.
            rgCurve_Ref: Reference.
        Returns:
            Single NurbsCurve on success.
            None on fail.
        """
        
        kwargs_reviewDuringDebug = kwargs

        ncOne_PreMatch = rgCurve_One.ToNurbsCurve()


        # (PC to maintain condition type at opposite end of curve) + (PC to modify working side of curve), where PC is "required min. point count".
        if sContinuity == 'G0':
            iCt_Cp_MinNeeded = iPreserveOtherEnd + 1
        elif sContinuity in ('G1', 'C1'):
            iCt_Cp_MinNeeded = iPreserveOtherEnd + 2
        elif sContinuity == 'G2':
            iCt_Cp_MinNeeded = iPreserveOtherEnd + 3


        if bDeformLines:
            if ncOne_PreMatch.Points.Count < iCt_Cp_MinNeeded:
                if ncOne_PreMatch.Points.Count == 2:
                    ts = ncOne_PreMatch.DivideByCount(
                        segmentCount=iCt_Cp_MinNeeded-1,
                        includeEnds=True)
                    rc = rg.NurbsCurve.Create(
                        periodic=False,
                        degree=3,
                        points=[ncOne_PreMatch.PointAt(t) for t in ts])
                    ncOne_PreMatch.Dispose()
                    ncOne_PreMatch = rc


        if bRebuildRationals and ncOne_PreMatch.IsRational:
            rc = rebuildRationalCurve(ncOne_PreMatch)
            if rc:
                ncOne_PreMatch.Dispose()
                ncOne_PreMatch = rc
            else:
                print "Rational curve was not rebuit."
                return


        if bT1WorkEnd_One:
            t_WorkEnd_One = ncOne_PreMatch.Domain.T1
            crvEvalSide_One = rg.CurveEvaluationSide.Below
            idxCp_Pos_One = ncOne_PreMatch.Points.Count - 1
            idxCp_Tan_One = ncOne_PreMatch.Points.Count - 2
            idxCp_Crv_One = ncOne_PreMatch.Points.Count - 3
        else:
            t_WorkEnd_One = ncOne_PreMatch.Domain.T0
            crvEvalSide_One = rg.CurveEvaluationSide.Above
            idxCp_Pos_One = 0
            idxCp_Tan_One = 1
            idxCp_Crv_One = 2
    
        nc_Ref = rgCurve_Ref.ToNurbsCurve()
        if bT1WorkEnd_Ref:
            t_WorkEnd_Ref = nc_Ref.Domain.T1
            crvEvalSide_Ref = rg.CurveEvaluationSide.Below
            idxCp_Pos_Ref = nc_Ref.Points.Count - 1
            idxCp_Tan_Ref = nc_Ref.Points.Count - 2
            idxCp_Crv_Ref = nc_Ref.Points.Count - 3
        else:
            t_WorkEnd_Ref = nc_Ref.Domain.T0
            crvEvalSide_Ref = rg.CurveEvaluationSide.Above
            idxCp_Pos_Ref = 0
            idxCp_Tan_Ref = 1
            idxCp_Crv_Ref = 2


        #vectsDerivatives_A = nc0_A.DerivativeAt(
        #        t_WorkEnd_A,
        #        derivativeCount=2,
        #        side=crvEvalSide_A)
        #print vectsDerivatives_A[2], nc0_A.CurvatureAt(t_WorkEnd_A)

        vectsDerivatives_B = nc_Ref.DerivativeAt(
                t_WorkEnd_Ref,
                derivativeCount=2,
                side=crvEvalSide_Ref)
        #print vectsDerivatives_B[2], nc0_B.CurvatureAt(t_WorkEnd_B)


        if ncOne_PreMatch.Points.Count < iCt_Cp_MinNeeded:
            if nc_Ref.Points.Count < iCt_Cp_MinNeeded:
                print "Neither curve has enough control points" \
                    " for this continuity modification."
                ncOne_PreMatch.Dispose()
                nc_Ref.Dispose()
                return
            if bEcho:
                s = "Curve A has {} control points, but needs {}.".format(
                        ncOne_PreMatch.Points.Count,
                        iCt_Cp_MinNeeded)
                if nc_Ref.Points.Count < iCt_Cp_MinNeeded:
                    s += "  Curve B also doesn't have enough."
                    if bRebuildRationals and nc_Ref.IsRational:
                        s += ", but unrationalization was not performed on rational B yet."
                    else:
                        s += "."
                else:
                    s += ", but Curve B does."
                print s
            ncOne_PreMatch.Dispose()
            nc_Ref.Dispose()
            return


        def createForCurvature():
            rgNurbsCrv_WIP_A = ncOne_PreMatch.Duplicate()
            vectTanForA = nc_Ref.TangentAt(t_WorkEnd_Ref)
        
            if bT1WorkEnd_One == bT1WorkEnd_Ref:
                bReversedTanForA = vectTanForA.Reverse()
                if bDebug: sEval='bReversedTanForA'; print sEval+': ',eval(sEval)

            bEndConditionSet_A = rgNurbsCrv_WIP_A.SetEndCondition(
                    bSetEnd=bT1WorkEnd_One,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Curvature,
                    point=nc_Ref.PointAtEnd if bT1WorkEnd_Ref else nc_Ref.PointAtStart,
                    tangent=vectTanForA,
                    curvature=nc_Ref.CurvatureAt(t_WorkEnd_Ref))
            return [rgNurbsCrv_WIP_A]


        def createForC1():
        
            if nc_Ref.CurvatureAt(t_WorkEnd_Ref).IsTiny(): # or not vectsDerivatives_B[2].IsTiny(): # Does it matter which one is used?
                return
        
            nc1_A = ncOne_PreMatch.Duplicate()

            # First, modify the end control point's position.
            if nc1_A.Points[idxCp_Pos_One] != nc_Ref.Points[idxCp_Pos_Ref]:
                nc1_A.Points[idxCp_Pos_One] = nc_Ref.Points[idxCp_Pos_Ref]

            # Now check control point spans to avoid distorting curve.
            dist_Cp0Cp2_A = nc1_A.Points[idxCp_Pos_One].Location.DistanceTo(
                    nc1_A.Points[idxCp_Crv_One].Location)
            dist_Cp0Cp1_R = nc_Ref.Points[idxCp_Pos_Ref].Location.DistanceTo(
                    nc_Ref.Points[idxCp_Tan_Ref].Location)
            degA = nc1_A.Degree
            degR = nc_Ref.Degree
            dist_Cp0Cp2_A_Scaled = float(degR)*dist_Cp0Cp2_A
            dist_Cp0Cp1_R_Scaled = float(degA)*dist_Cp0Cp1_R
            if dist_Cp0Cp2_A_Scaled < dist_Cp0Cp1_R_Scaled:
                if bDebug:
                    print "Avoided moving Tan Cp beyond Crv Cp."
                return

            # Now, modify the tangent (second from end) control point's position.
            # In the following, some of the paretheses are included
            # to obtain Point3d instead of Vector3d.
        
            pt_Target = (
                nc_Ref.Points[idxCp_Pos_Ref].Location -
                (
                    (float(nc_Ref.Degree) / float(nc1_A.Degree)) *
                    (nc_Ref.Points[idxCp_Tan_Ref].Location -
                    nc_Ref.Points[idxCp_Pos_Ref].Location)
                )
            )
            if nc1_A.Points[idxCp_Tan_One] != pt_Target:
                nc1_A.Points[idxCp_Tan_One] = pt_Target
            #sc.doc.Objects.AddCurve(nc1_A); sc.doc.Views.Redraw()
        
            return nc1_A


        def createForTangency():
            ncs = []
        
            # Version: C1 continuity.  Also used as an optional result for G1.
            rc = createForC1()
            if rc: ncs.append(rc)
            #

            if sContinuity == 'C1' and not nc_Ref.CurvatureAt(t_WorkEnd_Ref).IsTiny():
                # If IsTiny, control point span can be anything for C1,
                # so test other tangent solutions.
                return ncs

            #
            # Version: Move tangency control point to closest point inline with
            # tangency line of Curve B.
            ncs.append(ncOne_PreMatch.Duplicate())

            # First, modify the end control point's position.
            if ncs[-1].Points[idxCp_Pos_One] != nc_Ref.Points[idxCp_Pos_Ref]:
                ncs[-1].Points[idxCp_Pos_One] = nc_Ref.Points[idxCp_Pos_Ref]

            # Now, modify the tangent (second from end) control point's position
            # to be inline with the 2 end ones of Curve B.
            # Cannot modify the ControlPoint Location directly.  Use SetPoint instead.
            line = rg.Line(
                    nc_Ref.Points[idxCp_Pos_Ref].Location,
                    nc_Ref.Points[idxCp_Tan_Ref].Location)
            pt_Target = line.ClosestPoint(
                    ncs[-1].Points[idxCp_Tan_One].Location,
                    limitToFiniteSegment=False)
            ncs[-1].Points.SetPoint(idxCp_Tan_One, pt_Target)
            #sc.doc.Objects.AddCurve(ncs[-1]); sc.doc.Views.Redraw(); return
            #

            #
            # Version: Using NurbsCurve.SetEndCondition.
            ncs.append(ncOne_PreMatch.Duplicate())

            vectTanForA = nc_Ref.TangentAt(t_WorkEnd_Ref)
            if bT1WorkEnd_One == bT1WorkEnd_Ref:
                bReversedTanForA = vectTanForA.Reverse()
                if bDebug: sEval='bReversedTanForA'; print sEval+': ',eval(sEval)
            
            bEndConditionSet_A = ncs[-1].SetEndCondition(
                    bSetEnd=bT1WorkEnd_One,
                    continuity=rg.NurbsCurve.NurbsCurveEndConditionType.Tangency,
                    point=nc_Ref.PointAtEnd if bT1WorkEnd_Ref else nc_Ref.PointAtStart,
                    tangent=vectTanForA)
            #sc.doc.Objects.AddCurve(ncs[-1]); sc.doc.Views.Redraw(); return
            #

            return ncs


        def createForPosition():

            ncs = []

            #if nc0_A.Points[idxCp_Pos_A] == nc0_B.Points[idxCp_Pos_B]:
            #    print "Curves' endpoints already match."
            #    return

            # Version: Move only end control point.
            ncs.append(ncOne_PreMatch.Duplicate())
            ncs[-1].Points[idxCp_Pos_One] = nc_Ref.Points[idxCp_Pos_Ref]

            fRadius_Min_MoveEndCp = xCurve_radiusMinima.getMinimumRadius(ncs[-1])

            # Version: Using Curve.SetStartPoint / SetEndPoint method.
            ncs.append(ncOne_PreMatch.Duplicate())
            #attr_Red = rd.ObjectAttributes()
            #attr_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
            #attr_Red.ObjectColor = Color.Red
            if bT1WorkEnd_One:
                ncs[-1].SetEndPoint(nc_Ref.Points[idxCp_Pos_Ref].Location)
            else:
                ncs[-1].SetStartPoint(nc_Ref.Points[idxCp_Pos_Ref].Location)

            fRadius_Min_RcSet = xCurve_radiusMinima.getMinimumRadius(ncs[-1])

            if fRadius_Min_MoveEndCp >= fRadius_Min_RcSet:
                print "Moving just the end control point produces a result with" \
                    " minimum radius >= that from Curve.SetStartPoint/SetEndPoint."

            return ncs


        if sContinuity == 'G2':
            ncs = createForCurvature()
        elif sContinuity in ('G1', 'C1'):
            ncs = createForTangency()
        elif sContinuity == 'G0':
            ncs = createForPosition()

        if not ncs:
            ncOne_PreMatch.Dispose()
            nc_Ref.Dispose()
            return

        #for nc in ncs: sc.doc.Objects.AddCurve(nc)
        #sc.doc.Views.Redraw()

        if bEcho: print "{} curves initially generated.".format(len(ncs))

        if bMaximizeMinRadius and sContinuity not in ('C1', 'C2'):
            s = "Adjusting tangent control point spread ..."
            if bEcho: print s #Rhino.RhinoApp.CommandPrompt = s
            nc_MaxMinRads = []
            for nc in ncs:
                rc = xNurbsCurve_maximizeMinimumRadius.adjustTanCpSpread_OneEndOnly(
                        nc,
                        bT1WorkEnd_One,
                        fDevTol,
                        ncOne_PreMatch,
                        bDebug=bDebug,
                )
                #sc.doc.Objects.AddCurve(nc); sc.doc.Views.Redraw()
                if rc[0]: nc_MaxMinRads.append(rc[0])

            ncs.extend(nc_MaxMinRads)
            if bDebug: sEval='len(ncs)'; print sEval+': ',eval(sEval)


        devs = [getMaximumDeviation(ncOne_PreMatch, nc) for nc in ncs]
        if bDebug: sEval='devs'; print sEval+': ',eval(sEval)


        minRadii = []
        for nc in ncs:
            rad = xCurve_radiusMinima.getMinimumRadius(nc)
            minRadii.append(rad)
        if bDebug: sEval='minRadii'; print sEval+': ',eval(sEval)


        ncs_PassingTols = []
        devs_PassingTols = []
        minRadii_PassTols = []

        # Eliminate those curve which do not fall within the deviation tolerance, if any.
        for i in range(len(ncs)):
            if devs[i] is None:
                # Bad curve.
                ncs[i].Dispose()
                continue
            elif fDevTol is None:
                # All deviations except None will pass.
                pass
            elif devs[i] > fDevTol:
                ncs[i].Dispose()
                continue
            else:
                # Deviation is within tolerance.
                pass

            ncs_PassingTols.append(ncs[i])
            devs_PassingTols.append(devs[i])
            minRadii_PassTols.append(minRadii[i])

        if not ncs_PassingTols:
            min_dev = min(devs)
            idx_minDev = devs.index(min_dev)
            max_minRad = max(minRadii)
            idx_max_minRad = minRadii.index(max_minRad)

            s  = "No curves pass tolerances.  In set of curves,"
            s += "\nMinimumDev:{}".format(formatDistance(min_dev))
            s += " with MinRadius:{}".format(formatDistance(minRadii[idx_minDev]))
            s += "\nMaxMinRadius:{}".format(formatDistance(max_minRad))
            s += " with Dev:{}".format(formatDistance(devs[idx_max_minRad]))
            print s

            ncOne_PreMatch.Dispose()
            nc_Ref.Dispose()
            return

        if len(ncs_PassingTols) == 1:
                ncOne_PreMatch.Dispose()
                nc_Ref.Dispose()
                return ncs_PassingTols[0]


        if bMaximizeMinRadius:
            # TODO: Change this check from absolute max() to all values
            # within a practical range, e.g., 10.0*sc.doc.ModelAbsoluteTolerance.
            # This way, a curve similar in min. radius but with less deviation has
            # a better chance of being the winner.
            idx_Winner = minRadii_PassTols.index(max(minRadii_PassTols))
        else:
            # TODO: Change this check from absolute max() to all values
            # within a practical range, e.g., 0.1*sc.doc.ModelAbsoluteTolerance.
            # This way, a curve similar in deviation but with larger min. radius has
            # a better chance of being the winner.
            idx_Winner = devs_PassingTols.index(min(devs_PassingTols))

        for idx in range(len(ncs_PassingTols)):
            if idx == idx_Winner:
                nc = ncs_PassingTols[idx]
            else:
                ncs_PassingTols[idx].Dispose()


        if nc.EpsilonEquals(ncOne_PreMatch, 1e-12):
            if bDebug:
                print "New curve is within 1e-12 of old." \
                    "  The curve will not be modified."
            ncOne_PreMatch.Dispose()
            nc_Ref.Dispose()
            nc.Dispose()
            return

        ncOne_PreMatch.Dispose()
        nc_Ref.Dispose()

        return nc


    def createNcWithLessDevResult(rgCurveA, rgCurveB, bT1WorkEnd_A, bT1WorkEnd_B):
        """
        Returns on success:
            (NurbsCurve, None), (NurbsCurve, None), or (None, NurbsCurve)
        Returns on fail: None
        """

        c0_A = rgCurveA
        c0_B = rgCurveB


        rc = createForBothCurves(
            rgCurveA=rgCurveA,
            rgCurveB=rgCurveB,
            bT1WorkEnd_A=bT1WorkEnd_A,
            bT1WorkEnd_B=bT1WorkEnd_B,
            )
        ncs_Both = rc


        ncA_Alone = createForOneCurve(
                rgCurve_One=c0_A,
                rgCurve_Ref=c0_B,
                bT1WorkEnd_One=bT1WorkEnd_A,
                bT1WorkEnd_Ref=bT1WorkEnd_B,
        )

        ncB_Alone = createForOneCurve(
                rgCurve_One=c0_B,
                rgCurve_Ref=c0_A,
                bT1WorkEnd_One=bT1WorkEnd_B,
                bT1WorkEnd_Ref=bT1WorkEnd_A,
        )



        if not any((ncA_Alone, ncB_Alone, ncs_Both)):
            if bDebug: print "No curve could be created (within deviation?)."
            return

        if ncs_Both and (not ncA_Alone) and (not ncB_Alone):
            return ncs_Both

        if ncA_Alone and (not ncB_Alone) and (not ncs_Both):
            if bDebug:
                devA = getMaximumDeviation(c0_A, ncA_Alone)
                print "Curve could only be created for A (Deviation:{}).".format(
                        formatDistance(devA))
            return ncA_Alone, None

        if ncB_Alone and (not ncA_Alone) and (not ncs_Both):
            if bDebug:
                devB = getMaximumDeviation(c0_B, ncB_Alone)
                print "Curve could only be created for B (Deviation:{}).".format(
                        formatDistance(devB))
            return None, ncB_Alone


        # There are 2 or 3 sets of results.

        devA_Alone = getMaximumDeviation(c0_A, ncA_Alone)
        devB_Alone = getMaximumDeviation(c0_B, ncB_Alone)
        devA_Both = getMaximumDeviation(c0_A, ncs_Both[0])
        devB_Both = getMaximumDeviation(c0_B, ncs_Both[1])
        devs_Both = (devA_Both, devB_Both) if devA_Both and devB_Both else None

        if ncA_Alone and (devA_Alone < devB_Alone) and (devA_Alone < max(devs_Both)):
            if bDebug:
                s  = "Curve created for A has less deviation ({})".format(
                        formatDistance(devA_Alone))
                s += " than that for B ({}).".format(
                        formatDistance(devB_Alone))
                print s
            ncB_Alone.Dispose()
            return ncA_Alone, None

        if ncB_Alone and (devB_Alone < devA_Alone) and (devB_Alone < max(devs_Both)):
            if bDebug:
                s  = "Curve created for B has less deviation ({})".format(
                        formatDistance(devB_Alone))
                s += " than that for A ({}).".format(
                        formatDistance(devA_Alone))
                print s
            ncA_Alone.Dispose()
            return None, ncB_Alone

        if ncs_Both and (max(devs_Both) < devA_Alone) and (max(devs_Both) < devB_Alone):
            if bDebug:
                s  = "Curve created for both has less deviation ({})".format(
                        formatDistance(max(devs_Both)))
                s += " than that for either alone ({} and {}).".format(
                        formatDistance(devA_Alone), formatDistance(devB_Alone))
                print s
            if ncA_Alone: ncA_Alone.Dispose()
            if ncB_Alone: ncB_Alone.Dispose()
            return ncs_Both


            #            s  = "Created curves have the same deviation ({}).".format(
            #                    formatDistance(devA_Alone))
            #            s += "  Will modify A."
            #            ncB_Alone.Dispose()
            #            return ncA_Alone, None



    if bSkipIfAlreadyAtContinuity:
        if areCurvesAlreadyAtDesiredContinuity():
            return


    if bModifyA and not bModifyB:
        if bDebug: print "Modifying A and not B ..."
        rc = createForOneCurve(
            rgCurve_One=rgCurveA,
            rgCurve_Ref=rgCurveB,
            bT1WorkEnd_One=bT1WorkEnd_A,
            bT1WorkEnd_Ref=bT1WorkEnd_B,
            )
        if rc is None:
            return
        return rc, None
    elif not bModifyA and bModifyB:
        if bDebug: print "Modifying B and not A ..."
        rc = createForOneCurve(
            rgCurve_One=rgCurveB,
            rgCurve_Ref=rgCurveA,
            bT1WorkEnd_One=bT1WorkEnd_B,
            bT1WorkEnd_Ref=bT1WorkEnd_A,
            )
        if rc is None:
            return
        return None, rc
    else:
        # not bModifyA and not bModifyB
        if bDebug: print "Will look for result with less deviation ..."
        rc = createNcWithLessDevResult(
            rgCurveA=rgCurveA,
            rgCurveB=rgCurveB,
            bT1WorkEnd_A=bT1WorkEnd_A,
            bT1WorkEnd_B=bT1WorkEnd_B,
            )
        if rc is None:
            if bDebug:
                print "None from createNcWithLessDevResult."
            return
        return rc


def processCurveObjects(objrefs, **kwargs):

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    if 'sContinuity' in kwargs:
        sContinuity = kwargs['sContinuity']
    else:
        sContinuity = sOpts_Continuity[Opts.values['iContinuity']]
    fAngleTol_Deg = getOpt('fAngleTol_Deg')
    bSkipIfAlreadyAtContinuity = getOpt('bSkipIfAlreadyAtContinuity')
    iPreserveOtherEnd = getOpt('iPreserveOtherEnd')
    fDevTol = getOpt('fDevTol')
    bMaximizeMinRadius = getOpt('bMaximizeMinRadius')
    bModOnly1stPicked = getOpt('bModOnly1stPicked')
    bDeformLines = getOpt('bDeformLines')
    bRebuildRationals = getOpt('bRebuildRationals')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')
    
    if bDebug:
        for key in kwargs: print key, kwargs[key]

    gCrv0_A = objrefs[0].ObjectId
    c0_A = objrefs[0].Curve()
    
    if not bRebuildRationals:
        if isinstance(c0_A, rg.ArcCurve):
            print "ArcCurve skipped."
            c0_A.Dispose()
            return
        elif isinstance(c0_A, rg.NurbsCurve) and c0_A.IsRational and c0_A.IsArc(1e-9):
            print "Arc-shaped (within 1e-9), Rational NurbsCurve skipped."
            c0_A.Dispose()
            return


    bSuccess, tA = c0_A.ClosestPoint(objrefs[0].SelectionPoint())
    if not bSuccess:
        c0_A.Dispose()
        return
    

    gCrv0_B = objrefs[1].ObjectId
    c0_B = objrefs[1].Curve()
    bSuccess, tB = c0_B.ClosestPoint(objrefs[1].SelectionPoint())
    if not bSuccess:
        c0_A.Dispose()
        c0_B.Dispose()
        return
    
    bT1WorkEnd_A = tA > c0_A.Domain.Mid
    bT1WorkEnd_B = tB > c0_B.Domain.Mid

    if bModOnly1stPicked:
        bModifyA = True
        bModifyB = False
    else:
        bModifyA = True
        bModifyB = True


    rc = createNurbsCurves(
        rgCurveA=c0_A,
        rgCurveB=c0_B,
        bT1WorkEnd_A=bT1WorkEnd_A,
        bT1WorkEnd_B=bT1WorkEnd_B,
        bModifyA = bModifyA,
        bModifyB = bModifyB,
        sContinuity=sContinuity,
        fAngleTol_Deg=fAngleTol_Deg,
        bSkipIfAlreadyAtContinuity=bSkipIfAlreadyAtContinuity,
        iPreserveOtherEnd=iPreserveOtherEnd,
        fDevTol=fDevTol,
        bMaximizeMinRadius=bMaximizeMinRadius,
        bDeformLines=bDeformLines,
        bRebuildRationals=bRebuildRationals,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if rc is None:
        c0_A.Dispose()
        c0_B.Dispose()
        return
    rgNurbsCrv1_A, rgNurbsCrv1_B = rc


    if rgNurbsCrv1_A is None and rgNurbsCrv1_B is None:
        print "Curve(s) within tolerance not created, so they were not modified."

    def processResults(gCrv0_ToMod, rgCrv0_ToMod, rgNurbsCrv1):

        g1 = None

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

        s = ""

        rc = xCurve_inflections.getInflectionParameters(rgCrv0_ToMod)
        i_Infl_Ct_Pre = len(rc) if rc else None
        rc = xCurve_inflections.getInflectionParameters(rgNurbsCrv1)
        i_Infl_Ct_Post = len(rc) if rc else None
        s = "InflectionCount:{}->{}".format(i_Infl_Ct_Pre, i_Infl_Ct_Post)
    
        rc = xCurve_radiusMinima.getMinimumRadius(rgCrv0_ToMod)
        fRadius_Min_Original = rc if rc is not None else None
        rc = xCurve_radiusMinima.getMinimumRadius(rgNurbsCrv1)
        fRadius_Min_New = rc if rc is not None else None
        sRadius_Min_Original = formatDistance(fRadius_Min_Original)
        if bDebug: sEval='fRadius_Min_New'; print sEval+': ',eval(sEval)
        sRadius_Min_New = formatDistance(fRadius_Min_New)
        #sRadius_Min_New = '{:.{}f}'.format(fRadius_Min_New, sc.doc.ModelDistanceDisplayPrecision)
        if s: s += "  "
        s += "  MinimumRadius:{}->{}".format(sRadius_Min_Original, sRadius_Min_New)
    
        fDev = getMaximumDeviation(rgCrv0_ToMod, rgNurbsCrv1)
        s += "  Deviation:{}".format(formatDistance(fDev))
        print s

        return g1


    gs1 = []
    if rgNurbsCrv1_A:
        g1 = processResults(gCrv0_A, c0_A, rgNurbsCrv1_A)
        if g1: gs1.append(g1)
        rgNurbsCrv1_A.Dispose()
    if rgNurbsCrv1_B:
        g1 = processResults(gCrv0_B, c0_B, rgNurbsCrv1_B)
        if g1: gs1.append(g1)
        rgNurbsCrv1_B.Dispose()


    c0_A.Dispose()
    c0_B.Dispose()
    if rgNurbsCrv1_A: rgNurbsCrv1_A.Dispose()
    if rgNurbsCrv1_B: rgNurbsCrv1_B.Dispose()

    return gs1


def main():
    
    rc = getInput()
    if rc is None: return
    objrefs = rc[0]

    if Opts.values['bDebug']:
        print "Running with Debug mode on."
        import sys
        for sModule in list(sys.modules):
            if sModule[0] == 'x':
                try:
                    reload(sys.modules[sModule])
                    print "{} reloaded.".format(sModule)
                except:
                    s  = "{} NOT reloaded.".format(sModule)
                    s += "  Does the module contain a bug,"
                    s += " or was its name changed?"
                    print s
    else:
        sc.doc.Views.RedrawEnabled = False

    rc = processCurveObjects(
            objrefs,
            fDevTol=Opts.values['fDevTol'] if Opts.values['bLimitCrvDev'] else None,
    )

    if not rc: print "No curves were modified."

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()