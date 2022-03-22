"""
This script creates lofts or loft-like surfaces by choice of 3 methods:
    1. Standard loft.  This was added to the script to utilize the script's UX.
    2. 
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
200404: Created.
200515: Now aligns duplicates of the input curves.
211127: Added 2 more lofting methods from other scripts as options.
220322: Bug fix in resetting options.  UX improvements.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAlignCrvDirs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftMethod'; keys.append(key)
    values[key] = 1
    listValues[key] = 'Standard', 'Iterative', 'SrfPtGrid'
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftType'; keys.append(key)
    values[key] = 0
    listValues[key] = Enum.GetNames(rg.LoftType)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDegree'; keys.append(key)
    values[key] = 5
    riOpts[key] = ri.Custom.OptionInteger(values[key], lowerLimit=1, upperLimit=11)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bClosed'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPeriodic'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iGContinuity'; keys.append(key)
    values[key] = -1
    names[key] = 'MaxGContAdjustmentAtLoftEnds'
    riOpts[key] = ri.Custom.OptionInteger(values[key])
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

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'iDegree':
            if cls.riOpts[key].CurrentValue <= 0:
                cls.riOpts[key].CurrentValue = cls.values[key]
                return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_Curves():
    """
    """

    sk_PrevSel = 'UsePrevSelection({})({})'.format(__file__, sc.doc.Name)

    go = ri.Custom.GetObject()

    go.GeometryFilter = rd.ObjectType.Curve


    while True:
        go.ClearCommandOptions()

        idxs_Opt = {}

        go.SetCommandPrompt("Select curves in lofting order")

        if (sk_PrevSel in sc.sticky) and sc.sticky[sk_PrevSel]:
            key = 'UsePrevSelection'; idxs_Opt[key] = go.AddOption(key)

        key = 'ResetOptions'; idxs_Opt[key] = go.AddOption(key)

        res = go.GetMultiple(minimumNumber=2, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            sc.sticky[sk_PrevSel] = objrefs
            go.Dispose()
            return objrefs

        if go.Option().Index == idxs_Opt['UsePrevSelection']:
            objrefs = sc.sticky[sk_PrevSel]
            go.Dispose()
            return objrefs


        if go.Option().Index == idxs_Opt['ResetOptions']:
            for key in Opts.keys:
                if key in Opts.riOpts:
                    Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)


def _alignCrvDirs(crvs_In):
    """
    Returns list of new curves.
    """
    crvs_Out = [c.Duplicate() for c in crvs_In]
    for i in range(1, len(crvs_Out)):
        if not rg.Curve.DoDirectionsMatch(crvs_Out[0], crvs_Out[i]):
            crvs_Out[i].Reverse()
    return crvs_Out


def createLoft_CreateFromLoft(crvs, **kwargs):
    """
    Parameters:
        crvs
        bAlignCrvDirs
        iLoftType
        bClosed
        bPeriodic
        iGContinuity
        bEcho
        bDebug
    Returns on success: Rhino.Geometry.NurbsSurface
    Returns on fail: None
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAlignCrvDirs = getOpt('bAlignCrvDirs')
    iLoftType = getOpt('iLoftType')
    bClosed = getOpt('bClosed')
    bPeriodic = getOpt('bPeriodic')
    iGContinuity = getOpt('iGContinuity')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    loftType = Enum.ToObject(rg.LoftType, iLoftType)


    if bAlignCrvDirs:
        crvs = _alignCrvDirs(crvs)

    for i in range(1, len(crvs)):
        if not rg.Curve.DoDirectionsMatch(crvs[0], crvs[i]):
            print(crvs[i].Reverse())

    rgBreps_fromLoft = rg.Brep.CreateFromLoft(
        curves=crvs,
        start=rg.Point3d.Unset,
        end=rg.Point3d.Unset,
        loftType=loftType,
        closed=bClosed if len(crvs) > 2 else False)

    if not rgBreps_fromLoft:
        return

    if len(rgBreps_fromLoft) == 1:
        return rgBreps_fromLoft[0]

    rgBs_Joined = rg.Brep.JoinBreps(rgBreps_fromLoft)
    if len(rgBs_Joined) > 1:
        raise ValueError("More than one Brep resulted from Join.")
    return rgBs_Joined[0]


def createLoft_MatchGrevilles(ncs, **kwargs):
    """
    Parameters:
        ncs
        bAlignCrvDirs
        iDegree
        fDistTol
        bClosed
        bPeriodic
        iGContinuity
        bEcho
        bDebug
    Returns on success: Rhino.Geometry.NurbsSurface
    Returns on fail: None
    """


    if len(ncs) < 2:
        return

    for nc in ncs:
        if not isinstance(nc, Rhino.Geometry.NurbsCurve):
            print("Not all input curves are NURBS curves.")
            return

    for i in range(1, len(ncs)):
        if ncs[i].Degree != ncs[0].Degree:
            print("Various degrees in input curves.")
            return
        if ncs[i].Points.Count != ncs[0].Points.Count:
            print("Various point counts in input curves.")
            return

    for i in range(1, len(ncs)):
        if (
            ncs[i].SpanCount > 1 and
            ncs[i].Knots.KnotStyle != rg.KnotStyle.QuasiUniform
        ):
            print("Non-uniform curve in input.  Check results.")
            break


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAlignCrvDirs = getOpt('bAlignCrvDirs')
    iDegree = getOpt('iDegree')
    fDistTol = getOpt('fDistTol')
    bClosed = getOpt('bClosed')
    bPeriodic = getOpt('bPeriodic')
    iGContinuity = getOpt('iGContinuity')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if bAlignCrvDirs:
        ncs = _alignCrvDirs(ncs)


    if (iDegree == 2) and (len(ncs) == 2):
        breps = rg.Brep.CreateFromLoft(ncs, start=rg.Point3d.Unset, end=rg.Point3d.Unset, loftType=rg.LoftType.Straight)
        if len(breps) > 1:
            return
        brep = breps[0]
        if brep.Surfaces.Count > 1:
            return
        return brep.Surfaces[0]


    pts_InCols = []

    for nc in ncs:
        # TODO: Add support to maintain periodic knots in V direction.
        if nc.IsPeriodic: nc.Knots.ClampEnd(rg.CurveEnd.Both)
        pts_Col = [pt for pt in nc.GrevillePoints(all=nc.IsClosed)]
        pts_InCols.append(pts_Col)

    pts_Flat = [pt for pts_Col in pts_InCols for pt in pts_Col]


    uDegree = iDegree
    if (len(ncs) - 1) < uDegree:
        uDegree = len(ncs) - 1
        print("Degree along loft is made {} due to number of input curves.".format(uDegree))


    # Using CreateThroughPoints to create starting surface due to resultant
    # proximity of Grevilles points between surface and curves for later matching.
    ns = Rhino.Geometry.NurbsSurface.CreateFromPoints(
        points=pts_Flat,
        uCount=len(pts_InCols),
        vCount=len(pts_InCols[0]),
        uDegree=uDegree,
        vDegree=ncs[0].Degree)
        #uClosed=False,
        #vClosed=ncs[0].IsClosed)

    ns.KnotsU.CreateUniformKnots(knotSpacing=1.0)

    if ns.SpanCount(1) > 1:
        if ncs[0].Knots.KnotStyle == rg.KnotStyle.QuasiUniform:
            ns.KnotsV.CreateUniformKnots(knotSpacing=1.0)
        else:
            for iK in range(ncs[0].Knots.Count):
                ns.KnotsV[iK] = ncs[0].Knots[iK]

    #for i, pt in enumerate(pts_Flat):
    #    dot = rg.TextDot("[{}]".format(i), pt)
    #    sc.doc.Objects.AddTextDot(dot)
    #return



    Rhino.RhinoApp.SetCommandPrompt("Matching surface to curves ...")


    # Match first and last column CP locations to those of first and last curves, respectively.
    for iV in range(ncs[0].Points.Count):
        ns.Points.SetControlPoint(0,iV,ncs[0].Points[iV])
        ns.Points.SetControlPoint(ns.Points.CountU-1,iV,ncs[-1].Points[iV])


    iteration = 0

    zipped = zip(range(ns.Points.CountU), range(ns.Points.CountV))
    
    while True:
        if sc.escape_test(throw_exception=False):
            Rhino.RhinoApp.Wait()
            return
        
        bMovedPt = False
        
        iPt = 0

        for u in range(ns.Points.CountU):
            for v in range(ns.Points.CountV):
                pt_Target = pts_Flat[iPt]
                pt2d_Gr = ns.Points.GetGrevillePoint(u,v)
                pt3d_Gr = ns.PointAt(pt2d_Gr.X, pt2d_Gr.Y)
                dist = pt3d_Gr.DistanceTo(pt_Target)
                #print(pt3d_Gr, pt_Target, dist)
                if dist <= fDistTol:
                    iPt += 1
                    continue
                vect = pt_Target - pt3d_Gr
                pt_cp = ns.Points.GetControlPoint(u,v).Location
                ns.Points.SetControlPoint(u,v,rg.ControlPoint(pt_cp+vect))
                bMovedPt = True
                iPt += 1

        if not bMovedPt:
            print("Solution found after {} iterations.".format(iteration))
            break # out of while loop.

        iteration += 1


    return ns


def createLoft_CreateThroughPoints(ncs, **kwargs):
    """
    Parameters:
        ncs
        iDegree
        bClosed
        bPeriodic
        iGContinuity
        bEcho
        bDebug
    Returns on success: Rhino.Geometry.NurbsSurface
    Returns on fail: None
    """


    if len(ncs) < 2:
        return

    for nc in ncs:
        if not isinstance(nc, Rhino.Geometry.NurbsCurve):
            print("Not all input curves are NURBS curves.")
            return

    for i in range(1, len(ncs)):
        if ncs[i].Degree != ncs[0].Degree:
            print("Various degrees in input curves.")
            return
        if ncs[i].Points.Count != ncs[0].Points.Count:
            print("Various point counts in input curves.")
            return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bAlignCrvDirs = getOpt('bAlignCrvDirs')
    iDegree = getOpt('iDegree')
    bClosed = getOpt('bClosed')
    bPeriodic = getOpt('bPeriodic')
    iGContinuity = getOpt('iGContinuity')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if bAlignCrvDirs:
        ncs = _alignCrvDirs(ncs)


    if bClosed:
        ncs.append(ncs[0])
    else:
        bPeriodic = False


    uDegree = iDegree

    if bClosed or iGContinuity < 0:
        iAdditionalCps = 0
        if (len(ncs) - 1) < uDegree:
            uDegree = len(ncs) - 1
            print("Degree along loft is made {} due to number of input curves and" \
                " required end continuity.".format(uDegree))
    else:
        iAdditionalCps = uDegree//2 + iGContinuity + uDegree%2

        if len(ncs) + iAdditionalCps <= uDegree:
            while True:
                sc.escape_test()
                uDegree -= 1
                iAdditionalCps = uDegree//2 + iGContinuity + uDegree%2
                if len(ncs) + iAdditionalCps > uDegree:
                    break
            print("Degree along loft is made {} due to number of input curves and" \
                " required end continuity.".format(uDegree))

    if bDebug: sEval = "iAdditionalCps"; print("{}: {}".format(sEval, eval(sEval)))


    pts_InCols = []

    for nc in ncs:
        pts_Col = [pt for pt in nc.GrevillePoints(all=False)]
        pts_InCols.append(pts_Col)


    def insertAdditionalCps():
        cols_to_insert_NearT0 = []
        cols_to_insert_NearT1 = []

        numerator = 0.0
        denominator = float(iAdditionalCps**2 + 3*iAdditionalCps + 2) / 2.0
        if bDebug: sEval = "denominator"; print("{}: {}".format(sEval, eval(sEval)))

        for iDiv in range(1, iAdditionalCps+1):
            col_to_insert_NearT0 = []
            col_to_insert_NearT1 = []
            numerator += float(iDiv)
            ratioFromEnd = numerator/denominator

            for i in range(len(pts_InCols[0])):
                pt = (1.0 - ratioFromEnd) * pts_InCols[0][i] + (ratioFromEnd) * pts_InCols[1][i]
                col_to_insert_NearT0.append(pt)

                pt = (1.0 - ratioFromEnd) * pts_InCols[-1][i] + (ratioFromEnd) * pts_InCols[-2][i]
                col_to_insert_NearT1.append(pt)

            cols_to_insert_NearT0.append(col_to_insert_NearT0)
            cols_to_insert_NearT1.append(col_to_insert_NearT1)

        for i, col_to_insert_NearT0 in enumerate(cols_to_insert_NearT0):
            pts_InCols.insert(1+i, col_to_insert_NearT0)

        for i, col_to_insert_NearT1 in enumerate(cols_to_insert_NearT1):
            pts_InCols.insert(-1-i, col_to_insert_NearT1)

    if iAdditionalCps:
        insertAdditionalCps()


    #for pts_Col in pts_InCols:
    #    for pt in pts_Col:
    #        sc.doc.Objects.AddPoint(pt)

    if bDebug:
        sEval = "len(ncs)"; print("{}: {}".format(sEval, eval(sEval)))
        sEval = "len(pts_Col[0])"; print("{}: {}".format(sEval, eval(sEval)))

    #print(pts_InCols)



    def idx_flat_to_uv(iCountV, idx):
        return idx // iCountV, idx % iCountV


    def idx_uv_to_flat(iCountV, idxU, idxV):
        return idxU * iCountV + idxV



    pts_Flat = [pt for pts_Col in pts_InCols for pt in pts_Col]

    if bDebug:
        for pt in pts_Flat:
            sc.doc.Objects.AddPoint(pt)


    return Rhino.Geometry.NurbsSurface.CreateThroughPoints(
        points=pts_Flat,
        uCount=len(pts_InCols),
        vCount=len(pts_InCols[0]),
        uDegree=uDegree,
        vDegree=ncs[0].Degree,
        uClosed=bPeriodic,
        vClosed=ncs[0].IsClosed)


class DrawBrepConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.color = sc.doc.Layers.CurrentLayer.Color
        self.brep = None

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if self.brep:
            self.bbox = self.brep.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):
        if self.brep:
            drawEventArgs.Display.DrawBrepWires(
                brep=self.brep,
                color=self.color,
                wireDensity=1)
            drawEventArgs.Display.DrawBrepShaded(
                brep=self.brep,
                material=Rhino.Display.DisplayMaterial(diffuse=self.color))


def _createLoft_interactively(ncs_In):
    """
    Returns
        None to cancel.
        False to indicate to create objects with current options.
        True to indicate to regenerate geometry and return to this function.
    """
    
    valuesBefore = {}

    go = ri.Custom.GetOption()

    go.SetCommandPromptDefault("Accept") 

    go.AcceptNothing(True)


    sk_conduit = 'conduit({})'.format(__file__)
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
    else:
        conduit = DrawBrepConduit()
        sc.sticky[sk_conduit] = conduit


    idxs_Opt = {}

    bRecalc = True

    while True:

        sCurrentMethod = Opts.listValues['iLoftMethod'][Opts.values['iLoftMethod']]

        rgB_Out = None

        if not bRecalc:
            conduit.Enabled = True
            sc.doc.Views.Redraw()
        else:
            if sCurrentMethod == 'Standard':
                rgB_Out = createLoft_CreateFromLoft(ncs_In)
            elif sCurrentMethod == 'Iterative':
                ns_Res = createLoft_MatchGrevilles(ncs_In)
                if ns_Res:
                    rgB_Out = ns_Res.ToBrep()
            elif sCurrentMethod == 'SrfPtGrid':
                ns_Res = createLoft_CreateThroughPoints(ncs_In)
                if ns_Res:
                    rgB_Out = ns_Res.ToBrep()
            else:
                raise ValueError("LoftMethod value error.")

            if rgB_Out:
                conduit.brep = rgB_Out
                conduit.Enabled = True
                sc.doc.Views.Redraw()
            else:
                print("Loft was not created at current settings.")


        go.ClearCommandOptions()

        idxs_Opt.clear()

        go.SetCommandPrompt("{}".format(sCurrentMethod))

        for i, key in enumerate(Opts.listValues['iLoftMethod']):
            #if i == Opts.values['iLoftMethod']:
            #    idxs_Opt[key] = None
            #    continue
            idxs_Opt[key] = go.AddOption(key)

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bAlignCrvDirs')
        if sCurrentMethod == 'Standard':
            go.AcceptNumber(False, acceptZero=True)

            addOption('iLoftType')
            if len(ncs_In) > 2:
                addOption('bClosed')
        elif sCurrentMethod == 'Iterative':
            go.AcceptNumber(True, acceptZero=True)

            addOption('iDegree')
            addOption('fDistTol')
        elif sCurrentMethod == 'SrfPtGrid':
            go.AcceptNumber(True, acceptZero=True)

            addOption('iDegree')
            if len(ncs_In) > 2:
                addOption('bClosed')
                if Opts.values['bClosed']:
                    addOption('bPeriodic')
                else:
                    addOption('iGContinuity')
            else:
                addOption('iGContinuity')
        else:
            raise ValueError("LoftMethod value error.")

        addOption('bEcho')
        addOption('bDebug')


        for key in Opts.keys:
            valuesBefore[key] = Opts.values[key]

        res = go.Get()

        conduit.Enabled = False
        sc.doc.Views.Redraw()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            # Accept current result.
            go.Dispose()
            return rgB_Out

        if res == ri.GetResult.Number:
            key = 'iDegree'
            iNewDeg = abs(int(go.Number()))
            if not (0 < iNewDeg < 12):
                continue
            Opts.riOpts[key].CurrentValue = iNewDeg
            Opts.setValue(key)
            continue

        # An option was selected.

        def getLoftMethodChange():
            for key in Opts.listValues['iLoftMethod']:
                if go.Option().Index == idxs_Opt[key]:
                    return key

        rc = getLoftMethodChange()
        if rc:
            key = rc
            if key == Opts.listValues['iLoftMethod'][Opts.values['iLoftMethod']]:
                print("No change.")
                bRecalc = False
            Opts.setValue(
                'iLoftMethod',
                Opts.listValues['iLoftMethod'].index(key))
            bRecalc = True
            continue


        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                bRecalc = True
                break


def main():

    objrefs = getInput_Curves()
    if not objrefs: return
    
    ncs_FromIn = [o.Curve().ToNurbsCurve() for o in objrefs]

    rc = _createLoft_interactively(ncs_FromIn)

    if rc is None:
        # Cancel.
        return
    if not rc:
        # Create objects.
        return

    sc.doc.Objects.AddBrep(rc)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()