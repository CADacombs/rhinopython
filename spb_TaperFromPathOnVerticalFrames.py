"""
"""

"""
190604: Created from a split from another script.  Added bAddLoftEndCrvs.
190624: Now, command waits for Enter or options when curves are preselected.
190625: Added missing print statement.
190627: Fixed a bug in getInput.
190712,14: Polished the UX.  Added various options.  Corrected line orientation when CPlane is not parallel to CPlane World Top.
190727: Bug fixes.  Modified for improved debugging.
190903-04: Added bRebuildPath.  Modified for improved debugging.
210302: Modified an option default value.
220328: Import-related update.
220824: Added an option.
220828: Added an option.  Refactored.

TODO: Correctly create curve opposite path and brep when path curve is closed.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

import math

from System import Enum
from System import Guid

import spb_arrayObjsAtVerticalFrames


sBrepMethods = 'Loft2Crvs', 'LoftSectionLines', 'Sweep2A', 'Sweep2B', 'Network'

sOpts = (
    'bRebuildPath',
    'bEnterNum_TrueForDist_FalseForAngle',
    'fDistance',
    'bProjDist',
    'fTaper_Start_Deg',
    'bVariableTaper',
    'fTaper_End_Deg',
    'bTaperChangePerCrvParam',
    'bCPlane',
    'bAtGrevilles',
    'bAtKnots',
    'bAtEqualDivisions',
    'iDivisionCt',
    'bSplitPolyCrvToSegs',
    'bSplitPathsAtKnots',
    'bAddArrayedLines',
    'bAddLoftEndCrvs',
    'bAddBrep',
    'iBrepMethod',
    'iLoftType',
    'fBrepTol',
    'bEcho',
    'bDebug',
    )


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bRebuildPath'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEnterNum_TrueForDist_FalseForAngle'; keys.append(key)
    values[key] = True
    names[key] = 'NumberEntry'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'DraftAngle', 'Dist')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDistance'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bProjDist'; keys.append(key)
    values[key] = True
    names[key] = 'DistType'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'True', 'Projected')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTaper_Start_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'TaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bVariableTaper'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTaper_End_Deg'; keys.append(key)
    values[key] = 45.0
    names[key] = 'EndTaperAngle'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTaperChangePerCrvParam'; keys.append(key)
    values[key] = False
    names[key] = 'TaperChangePerCrv'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Length', 'Param')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCPlane'; keys.append(key)
    values[key] = True
    names[key] = 'PlanView'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtGrevilles'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAtEqualDivisions'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDivisionCt'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(
            initialValue=values[key],
            setLowerLimit=True,
            limit=2)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitPolyCrvToSegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSplitPathsAtKnots'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddArrayedLines'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddLoftEndCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddBrep'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iBrepMethod'; keys.append(key)
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iLoftType'; keys.append(key)
    values[key] = 0
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fBrepTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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


def getInput(bFirstGetObjects=False):
    """
    Get objects with optional input.
    
    Returns
        None to cancel.
        tuple of (list(ObjRefs), bool(Generate new geometry), bool(Accept results):
    """


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


    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select path curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    # Due to the object collection below, do not use AcceptNothing.
    go.AcceptNumber(True, acceptZero=False)

    idxs_Opts = {}

    go.AddOptionToggle(Opts.names['bRebuildPath'], Opts.riOpts['bRebuildPath'])
    go.AddOptionToggle(Opts.names['bEnterNum_TrueForDist_FalseForAngle'], Opts.riOpts['bEnterNum_TrueForDist_FalseForAngle'])
    go.AddOptionDouble(Opts.names['fDistance'], Opts.riOpts['fDistance'])
    if not Opts.values['bVariableTaper']:
        go.AddOptionToggle(Opts.names['bProjDist'], Opts.riOpts['bProjDist'])
    go.AddOptionDouble(Opts.names['fTaper_Start_Deg'], Opts.riOpts['fTaper_Start_Deg'])
    go.AddOptionToggle(Opts.names['bVariableTaper'], Opts.riOpts['bVariableTaper'])
    if Opts.values['bVariableTaper']:
        go.AddOptionDouble(Opts.names['fTaper_End_Deg'], Opts.riOpts['fTaper_End_Deg'])
        idxs_Opts['SwapAngles'] = go.AddOption('SwapAngles')
        go.AddOptionToggle(Opts.names['bTaperChangePerCrvParam'], Opts.riOpts['bTaperChangePerCrvParam'])
    idxs_Opts['FlipAngle'] = go.AddOption('FlipAngle')
    idxs_Opts['FlipDir'] = go.AddOption('FlipDir')
    go.AddOptionToggle(Opts.names['bCPlane'], Opts.riOpts['bCPlane'])
    go.AddOptionToggle(Opts.names['bAtGrevilles'], Opts.riOpts['bAtGrevilles'])
    go.AddOptionToggle(Opts.names['bAtKnots'], Opts.riOpts['bAtKnots'])
    go.AddOptionToggle(Opts.names['bAtEqualDivisions'],
                        Opts.riOpts['bAtEqualDivisions'])
    if Opts.values['bAtEqualDivisions']:
        go.AddOptionInteger(Opts.names['iDivisionCt'],
                            Opts.riOpts['iDivisionCt'])
    go.AddOptionToggle(Opts.names['bSplitPolyCrvToSegs'], Opts.riOpts['bSplitPolyCrvToSegs'])
    go.AddOptionToggle(Opts.names['bSplitPathsAtKnots'], Opts.riOpts['bSplitPathsAtKnots'])
    go.AddOptionToggle(Opts.names['bAddArrayedLines'], Opts.riOpts['bAddArrayedLines'])
    go.AddOptionToggle(Opts.names['bAddLoftEndCrvs'],
                        Opts.riOpts['bAddLoftEndCrvs'])
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
    go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
    go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
    go.AlreadySelectedObjectSelect = True # So objects can be reselected after being unselected in same go.GetMultiple.
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

    if not go.ObjectsWerePreselected:
        objrefs = None
        iCt_Crvs_PreSelctd = 0
    else:
        if bFirstGetObjects:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs, True, False

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

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Object:
        objrefs = go.Objects()
        if iCt_Crvs_PreSelctd == go.ObjectCount:
            if not wasThereChangeInObjectsSelected(
                    gCrvs_PreSelctd=gCrvs_PreSelctd,
                    gBreps_PreSelctd=gBreps_PreSelctd,
                    idxs_Edges_PerBrep=idxs_Edges_PerBrep
            ):
                go.Dispose()
                return objrefs, False, True
        go.Dispose()
        return objrefs, True, False
    elif res == ri.GetResult.Number:
        if Opts.values['bEnterNum_TrueForDist_FalseForAngle']:
            Opts.riOpts['fDistance'].CurrentValue = go.Number()
        else:
            Opts.riOpts['fTaper_Start_Deg'].CurrentValue = go.Number()
    elif Opts.values['bVariableTaper'] and go.OptionIndex() == idxs_Opts['SwapAngles']:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue, Opts.riOpts['fTaper_End_Deg'].CurrentValue = (
                Opts.riOpts['fTaper_End_Deg'].CurrentValue, Opts.riOpts['fTaper_Start_Deg'].CurrentValue)
    elif go.OptionIndex() == idxs_Opts['FlipAngle']:
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
        Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
    elif go.OptionIndex() == idxs_Opts['FlipDir']:
        Opts.riOpts['fDistance'].CurrentValue = -Opts.riOpts['fDistance'].CurrentValue
        Opts.riOpts['fTaper_Start_Deg'].CurrentValue = -Opts.riOpts['fTaper_Start_Deg'].CurrentValue
        Opts.riOpts['fTaper_End_Deg'].CurrentValue = -Opts.riOpts['fTaper_End_Deg'].CurrentValue
    elif Opts.values['bAddBrep'] and go.OptionIndex() == idxs_Opts['iBrepMethod']:
        Opts.values['iBrepMethod'] = go.Option().CurrentListOptionIndex
    elif Opts.values['bAddBrep'] and Opts.values['iBrepMethod'] == 1 and go.OptionIndex() == idxs_Opts['iLoftType']:
        Opts.values['iLoftType'] = go.Option().CurrentListOptionIndex
    elif Opts.riOpts['fBrepTol'].CurrentValue < 0.0:
        Opts.riOpts['fBrepTol'].CurrentValue = Opts.riOpts['fBrepTol'].InitialValue

    Opts.setValues()
    Opts.saveSticky()

    return objrefs, True, False


def createBrep(iBrepMethod, iLoftType, fBrepTol, rgCrv_Path, rgNurbsCrv_TaperEnd_1Seg, rgLineCrvs_Arrayed):
    """
    """

    rgNurbsCrv1_PathSeg = rgCrv_Path.ToNurbsCurve()

    rgBreps1 = []
    
    if sBrepMethods[Opts.values['iBrepMethod']] == 'Loft2Crvs':
        rgBreps1 = rg.Brep.CreateFromLoft(
                curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg],
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=rg.LoftType.Straight,
                closed=False)
    elif sBrepMethods[Opts.values['iBrepMethod']] == 'LoftSectionLines':
        for L in rgLineCrvs_Arrayed:
            print L.PointAtStart, L.PointAtEnd
        print rgLineCrvs_Arrayed[0].PointAtEnd.EpsilonEquals(rgLineCrvs_Arrayed[-1].PointAtEnd, epsilon=1e-12)
        rgBreps1 = rg.Brep.CreateFromLoft(
                curves=rgLineCrvs_Arrayed,
                start=rg.Point3d.Unset,
                end=rg.Point3d.Unset,
                loftType=Enum.ToObject(rg.LoftType, iLoftType),
                closed=rgNurbsCrv1_PathSeg.IsClosed)
        print rgBreps1
    elif sBrepMethods[Opts.values['iBrepMethod']] == 'Sweep2A':
        rgBreps1 = rg.Brep.CreateFromSweep(
                rail1=rgNurbsCrv1_PathSeg,
                rail2=rgNurbsCrv_TaperEnd_1Seg,
                shapes=rgLineCrvs_Arrayed,
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
                crossSections=rgLineCrvs_Arrayed)
    elif sBrepMethods[Opts.values['iBrepMethod']] == 'Network':
        rgNurbsSrf, iError = rg.NurbsSurface.CreateNetworkSurface(
                curves=[rgNurbsCrv1_PathSeg, rgNurbsCrv_TaperEnd_1Seg]+rgLineCrvs_Arrayed,
                continuity=1,
                edgeTolerance=fBrepTol,
                interiorTolerance=fBrepTol,
                angleTolerance=0.1*sc.doc.ModelAngleToleranceDegrees)
        if iError:
            print "CreateNetworkSurface error code: {}".format(iError)
        else:
            rgBreps1 = [rgNurbsSrf.ToBrep()]
            rgNurbsSrf.Dispose()

    rgNurbsCrv1_PathSeg.Dispose()

    return rgBreps1


def main():

    while True:
        sc.escape_test()

        rc = getInput(bFirstGetObjects=True)
        if rc is None: return

        objrefs_Paths, bGenerateNew, bAcceptResults = rc
        if not objrefs_Paths: continue

        bRebuildPath = Opts.values['bRebuildPath']
        fDistance = Opts.values['fDistance']
        bProjDist = Opts.values['bProjDist']
        fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
        bVariableTaper = Opts.values['bVariableTaper']
        fTaper_End_Deg = Opts.values['fTaper_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        bCPlane = Opts.values['bCPlane']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
        bSplitPathsAtKnots = Opts.values['bSplitPathsAtKnots']
        bAddArrayedLines = Opts.values['bAddArrayedLines']
        bAddLoftEndCrvs = Opts.values['bAddLoftEndCrvs']
        bAddBrep = Opts.values['bAddBrep']
        iBrepMethod = Opts.values['iBrepMethod']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']

        break
    
    if bDebug:
        reload(spb_arrayObjsAtVerticalFrames)
    
    rgLine_ToArray = None
    rgCrvs0_Path = []
    rgCrvs1_Path = []

    while True:

        if bDebug:
            reload(spb_arrayObjsAtVerticalFrames)

        while not (bAtGrevilles or bAtKnots or bAtEqualDivisions):
            print "No path point sampling is enabled."
            sc.doc.Views.Redraw()
            rc = getInput(bFirstGetObjects=False)
            if rc is None:
                if rgLine_ToArray: rgLine_ToArray.Dispose()
                for c in rgCrvs0_Path: c.Dispose()
                for c in rgCrvs1_Path: c.Dispose()
                return

            objrefs_Paths, bGenerateNew, bAcceptResults = rc

            if bAcceptResults:
                if rgLine_ToArray: rgLine_ToArray.Dispose()
                for c in rgCrvs0_Path: c.Dispose()
                for c in rgCrvs1_Path: c.Dispose()
                return

            bRebuildPath = Opts.values['bRebuildPath']
            fDistance = Opts.values['fDistance']
            bProjDist = Opts.values['bProjDist']
            fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
            bVariableTaper = Opts.values['bVariableTaper']
            fTaper_End_Deg = Opts.values['fTaper_End_Deg']
            bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
            bCPlane = Opts.values['bCPlane']
            bAtGrevilles = Opts.values['bAtGrevilles']
            bAtKnots = Opts.values['bAtKnots']
            bAtEqualDivisions = Opts.values['bAtEqualDivisions']
            iDivisionCt = Opts.values['iDivisionCt']
            bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
            bSplitPathsAtKnots = Opts.values['bSplitPathsAtKnots']
            bAddArrayedLines = Opts.values['bAddArrayedLines']
            bAddLoftEndCrvs = Opts.values['bAddLoftEndCrvs']
            bAddBrep = Opts.values['bAddBrep']
            iBrepMethod = Opts.values['iBrepMethod']
            iLoftType = Opts.values['iLoftType']
            fBrepTol = Opts.values['fBrepTol']
            bEcho = Opts.values['bEcho']
            bDebug = Opts.values['bDebug']


        rgCrvs0_Path = []
        for objref_Path in objrefs_Paths:
            c = objref_Path.Curve()
            rgCrvs0_Path.append(c)

        Rhino.RhinoApp.SetCommandPrompt("Preparing path curves ...")
        rc = spb_arrayObjsAtVerticalFrames.prepareCurves(
                rgCrvs0=rgCrvs0_Path,
                bRebuild=bRebuildPath,
                bSplitPolyCrvToSegs=bSplitPolyCrvToSegs,
                bSplitPathsAtKnots=bSplitPathsAtKnots)
        if not rc: return
        rgCrvs1_Path = rc

        if bProjDist and not bVariableTaper:
            rgLine_ToArray = rg.Line(
                rg.Point3d(0.0, 0.0, 0.0),
                rg.Point3d(
                    0.0,
                    fDistance/math.cos(math.radians(fTaper_Start_Deg)),
                    0.0))
        else:
            rgLine_ToArray = rg.Line(
                rg.Point3d(0.0, 0.0, 0.0),
                rg.Point3d(0.0, fDistance, 0.0))
    
        if bCPlane:
            view_Active = sc.doc.Views.ActiveView
            plane_Proj = view_Active.ActiveViewport.ConstructionPlane()
        
            xform1 = rg.Transform.PlaneToPlane(rg.Plane.WorldXY, plane_Proj)
            rgLine_ToArray.Transform(xform1)
        else:
            plane_Proj = rg.Plane.WorldXY

        rgLineCrv_ToArray = rg.LineCurve(rgLine_ToArray)

        gLineCrvs1_Arrayed = []
        gCrvs_TaperEnds_All = []
        rgBreps1 = []
        gBreps1 = []

        rgCrvs_PathSegs = None
        
        Rhino.RhinoApp.SetCommandPrompt("Creating geometry ...")

        for rgCrv1_Path_1Seg in rgCrvs1_Path:
            rc = spb_arrayObjsAtVerticalFrames.\
                    createArrayedGeometry(
                        rgCrv_Path=rgCrv1_Path_1Seg,
                        rgObjs_ToArray=[rgLineCrv_ToArray],
                        plane_Proj=plane_Proj,
                        fTaper_Start_Deg=fTaper_Start_Deg,
                        fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                        bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                        bAtGrevilles=bAtGrevilles,
                        bAtKnots=bAtKnots,
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
                if bAtKnots or bAtEqualDivisions:
                    if bAtKnots and iBrepMethod == 0:
                        s  = "AtKnots"
                        s += " only affects added arrayed lines"
                        s += " when BrepMethod == {},".format(sBrepMethods[iBrepMethod])
                        s += " but AddArrayed option is disabled."
                        print s
                    if bAtEqualDivisions and iBrepMethod == 0:
                        s  = "AtEqualDivisions"
                        s += " only affects added arrayed lines"
                        s += " when BrepMethod == {},".format(sBrepMethods[iBrepMethod])
                        s += " but AddArrayed option is disabled."
                        print s

            if bAddLoftEndCrvs or bAddBrep:

                if not bAtGrevilles or bAtKnots or bAtEqualDivisions:
                    rc = spb_arrayObjsAtVerticalFrames.\
                            createArrayedGeometry(
                                rgCrv_Path=rgCrv1_Path_1Seg,
                                rgObjs_ToArray=[rgLineCrv_ToArray],
                                plane_Proj=plane_Proj,
                                fTaper_Start_Deg=fTaper_Start_Deg,
                                fTaper_End_Deg=fTaper_End_Deg if bVariableTaper else fTaper_Start_Deg,
                                bTaperChangePerCrvParam=bTaperChangePerCrvParam,
                                bAtGrevilles=True,
                                bAtKnots=False,
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

                if bAddLoftEndCrvs:
                    gCrv_TaperEnd_1Seg = sc.doc.Objects.AddCurve(rgNurbsCrv_TaperEnd)
                    if gCrv_TaperEnd_1Seg != Guid.Empty:
                        gCrvs_TaperEnds_All.append(gCrv_TaperEnd_1Seg)
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

        rc = getInput(bFirstGetObjects=False)
        if rc is None:
            for g_ in gLineCrvs1_Arrayed:
                sc.doc.Objects.Delete(g_, True)
            for g_ in gCrvs_TaperEnds_All:
                sc.doc.Objects.Delete(g_, True)
            for g_ in gBreps1:
                sc.doc.Objects.Delete(g_, True)
            break

        objrefs_Paths, bGenerateNew, bAcceptResults = rc

        if bAcceptResults:
            break

        if not bGenerateNew:
            continue

        bRebuildPath = Opts.values['bRebuildPath']
        fDistance = Opts.values['fDistance']
        bProjDist = Opts.values['bProjDist']
        fTaper_Start_Deg = Opts.values['fTaper_Start_Deg']
        bVariableTaper = Opts.values['bVariableTaper']
        fTaper_End_Deg = Opts.values['fTaper_End_Deg']
        bTaperChangePerCrvParam = Opts.values['bTaperChangePerCrvParam']
        bCPlane = Opts.values['bCPlane']
        bAtGrevilles = Opts.values['bAtGrevilles']
        bAtKnots = Opts.values['bAtKnots']
        bAtEqualDivisions = Opts.values['bAtEqualDivisions']
        iDivisionCt = Opts.values['iDivisionCt']
        bSplitPolyCrvToSegs = Opts.values['bSplitPolyCrvToSegs']
        bSplitPathsAtKnots = Opts.values['bSplitPathsAtKnots']
        bAddArrayedLines = Opts.values['bAddArrayedLines']
        bAddLoftEndCrvs = Opts.values['bAddLoftEndCrvs']
        bAddBrep = Opts.values['bAddBrep']
        iBrepMethod = Opts.values['iBrepMethod']
        iLoftType = Opts.values['iLoftType']
        fBrepTol = Opts.values['fBrepTol']
        bEcho = Opts.values['bEcho']
        bDebug = Opts.values['bDebug']


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