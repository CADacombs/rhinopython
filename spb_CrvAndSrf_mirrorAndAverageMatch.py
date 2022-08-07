"""
Behavior notes:
    Symmetry planes are per World, not CPlane.
    For curves, only wires are modified.
    For surfaces, only untrimmed ones are modified.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220807: Created by combining 2 other scripts.
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


    key = 'fDistTol'; keys.append(key)
    values[key] = 1.0
    names[key] = 'MaxDist'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iSymPlane'; keys.append(key)
    values[key] = 0
    listValues[key] = 'X', 'Y', 'Z'
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

        if key == 'fDistTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            if cls.riOpts[key].CurrentValue <= Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = cls.values[key] = Rhino.RhinoMath.ZeroTolerance
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
        else:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            elif key in cls.listValues:
                cls.values[key] = idxList
            else:
                return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def _getInput_WireCrvsAndUntrimmedSrfs():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select wire curves and/or untrimmed surfaces")

    go.GeometryFilter = (
        rd.ObjectType.Curve
        |
        rd.ObjectType.Surface
        )
    go.GeometryAttributeFilter = (
        ri.Custom.GeometryAttributeFilter.WireCurve
        |
        ri.Custom.GeometryAttributeFilter.UntrimmedSurface
        )

    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fDistTol')
        #addOption('bEcho')
        #addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fDistTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _getInput_SymmetryPlane():
    """
    Get option.
    """

    idxOpt_Default = Opts.values['iSymPlane']

    go = ri.Custom.GetOption()

    go.SetCommandPrompt("Pick mirror plane")
    go.SetCommandPromptDefault(Opts.listValues['iSymPlane'][idxOpt_Default])
    go.AcceptNothing(True)

    idxs_Opt = {}

    for sOpt in Opts.listValues['iSymPlane']:
        idxs_Opt[sOpt] = go.AddOption(sOpt)

    res = go.Get()


    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    if res == ri.GetResult.Nothing:
        idxOpt = idxOpt_Default
    elif res == ri.GetResult.Option:
        idxOpt = go.Option().Index - 1 # '- 1' because go.Option().Index is base 1.
        Opts.setValue('iSymPlane', idxOpt)

    go.Dispose()

    return Opts.listValues['iSymPlane'][idxOpt]


def _getPerpCompFromPlane(sPlane):
    if sPlane == 'XY': return 'Z'
    if sPlane == 'YZ': return 'X'
    if sPlane == 'ZX': return 'Y'


def _getLocComp(cp, sLocComp):
    if sLocComp == 'X':
        return cp.X
    if sLocComp == 'Y':
        return cp.Y
    if sLocComp == 'Z':
        return cp.Z


def getSide(ns_In, sPlane, fDistTol):

    sCompZero = _getPerpCompFromPlane(sPlane)

    for iu in range(ns_In.Points.CountU):
        cp = ns_In.Points.GetControlPoint(iu, 0)
        if abs(_getLocComp(cp, sCompZero)) > fDistTol:
            break
    else:
        return rg.IsoStatus.South

    for iv in range(ns_In.Points.CountV):
        cp = ns_In.Points.GetControlPoint(ns_In.Points.CountU-1, iv)
        if abs(_getLocComp(cp, sCompZero)) > fDistTol:
            break
    else:
        return rg.IsoStatus.East

    for iu in range(ns_In.Points.CountU):
        cp = ns_In.Points.GetControlPoint(iu, ns_In.Points.CountV-1)
        if abs(_getLocComp(cp, sCompZero)) > fDistTol:
            break
    else:
        return rg.IsoStatus.North

    for iv in range(ns_In.Points.CountV):
        cp = ns_In.Points.GetControlPoint(0, iv)
        if abs(_getLocComp(cp, sCompZero)) > fDistTol:
            break
    else:
        return rg.IsoStatus.West


def createSymmetricPair_Brep(ns_In, side, sPlane_Sym):
    
    ns_Out_A = ns_In.Duplicate()
    
    sComp_Perp = _getPerpCompFromPlane(sPlane_Sym)
    
    ius_P0 = []
    ivs_P0 = []
    ius_P1 = []
    ivs_P1 = []
    
    if side == rg.IsoStatus.South:
        for iu in range(ns_Out_A.Points.CountU):
            ius_P0.append(iu)
            ivs_P0.append(0)
            ius_P1.append(iu)
            ivs_P1.append(1)
    elif side == rg.IsoStatus.East:
        for iv in range(ns_Out_A.Points.CountV):
            ius_P0.append(ns_Out_A.Points.CountU-1)
            ivs_P0.append(iv)
            ius_P1.append(1)
            ivs_P1.append(ns_Out_A.Points.CountU-2)
    elif side == rg.IsoStatus.North:
        for iu in range(ns_Out_A.Points.CountU):
            ius_P0.append(iu)
            ivs_P0.append(ns_Out_A.Points.CountV-1)
            ius_P1.append(iu)
            ivs_P1.append(ns_Out_A.Points.CountV-2)
    elif side == rg.IsoStatus.West:
        for iv in range(ns_Out_A.Points.CountV):
            ius_P0.append(0)
            ivs_P0.append(iv)
            ius_P1.append(1)
            ivs_P1.append(iv)
    
    
    def setLocComp(cp, sLocComp, fComp=0.0):
        if sLocComp == 'X':
            cp.X = fComp
            return cp.X == fComp
        if sLocComp == 'Y':
            cp.Y = fComp
            return cp.Y == fComp
        if sLocComp == 'Z':
            cp.Z = fComp
            return cp.Z == fComp
    
    
    for iu_P0, iv_P0, iu_P1, iv_P1 in zip(ius_P0, ivs_P0, ius_P1, ivs_P1):
        cpA = ns_Out_A.Points.GetControlPoint(iu_P0, iv_P0)
        cpB = ns_Out_A.Points.GetControlPoint(iu_P1, iv_P1)
        
        delta_Comp = -_getLocComp(cpA, sComp_Perp)
        
        if delta_Comp != 0.0:
            setLocComp(cpA, sComp_Perp, 0.0)
            ns_Out_A.Points.SetControlPoint(iu_P0, iv_P0, cpA)
            
            setLocComp(
                cpB,
                sComp_Perp,
                delta_Comp + _getLocComp(cpB, sComp_Perp)
                )
            
            ns_Out_A.Points.SetControlPoint(iu_P1, iv_P1, cpB)
        
        for sComp in sPlane_Sym:
            compA = _getLocComp(cpA, sComp)
            if compA != _getLocComp(cpB, sComp):
                setLocComp(cpB, sComp, compA)
                
                ns_Out_A.Points.SetControlPoint(iu_P1, iv_P1, cpB)
    
    
    ns_Out_B = ns_Out_A.Duplicate()
    
    if sPlane_Sym == 'XY':
        plane_sym = rg.Plane.WorldXY
    elif sPlane_Sym == 'YZ':
        plane_sym = rg.Plane.WorldYZ
    elif sPlane_Sym == 'ZX':
        plane_sym = rg.Plane.WorldZX
    else:
        raise Exception("What happened?") 
    
    xform = rg.Transform.Mirror(plane_sym)
    
    if not ns_Out_B.Transform(xform):
        print("Transformation failed.")
        ns_Out_A.Dispose()
        ns_Out_B.Dispose()
        return
    
    return ns_Out_A, ns_Out_B


def createSymmetricPair_Curve(nc, plane_sym, bT1):
    
    nc_Out_A = nc.DuplicateCurve()
    nc_Out_B = nc.DuplicateCurve()
    
    xform = rg.Transform.Mirror(plane_sym)
    
    if not nc_Out_B.Transform(xform):
        print("Transformation failed.")
        nc_Out_A.Dispose()
        nc_Out_B.Dispose()
        return
    
    return rg.Curve.CreateMatchCurve(
        curve0=nc_Out_A,
        reverse0=bT1,
        continuity=rg.BlendContinuity.Tangency,
        curve1=nc_Out_B,
        reverse1=bT1,
        preserve=rg.PreserveEnd.Tangency,
        average=True)


def processBrepObject(objref_In, sMirrorOpt, fDistTol):
    rgF_In = objref_In.Face()
    rgS_In = rgF_In.UnderlyingSurface()
    if not isinstance(rgS_In, rg.NurbsSurface): return
    
    ns_In = rgS_In
    
    if sMirrorOpt == 'X':
        sPlane_Sym = 'ZX'
    elif sMirrorOpt == 'Y':
        sPlane_Sym = 'YZ'
    elif sMirrorOpt == 'Z':
        sPlane_Sym = 'XY'
    else:
        raise ValueError("{} not a value symmetry plane description.".format(
            sMirrorOpt))

    side = getSide(ns_In, sPlane_Sym, fDistTol)
    if side is None: return
    
    nss_WIP = createSymmetricPair_Brep(ns_In, side, sPlane_Sym)
    if nss_WIP is None: return

    if not sc.doc.Objects.Replace(objref_In, nss_WIP[0]):
        return

    gMirror = sc.doc.Objects.AddSurface(nss_WIP[1], objref_In.Object().Attributes)
    if gMirror == gMirror.Empty: return

    return objref_In.ObjectId, gMirror


def processCurveObject(objref_In, sMirrorOpt, fDistTol):
    rgC_In = objref_In.Curve()
    if not isinstance(rgC_In, rg.NurbsCurve): return
    
    nc_In = rgC_In
    
    if sMirrorOpt == 'X':
        plane_sym = rg.Plane.WorldZX
        fStartValue = nc_In.PointAtStart.Y
        fEndValue = nc_In.PointAtEnd.Y
    elif sMirrorOpt == 'Y':
        plane_sym = rg.Plane.WorldYZ
        fStartValue = nc_In.PointAtStart.X
        fEndValue = nc_In.PointAtEnd.X
    elif sMirrorOpt == 'Z':
        plane_sym = rg.Plane.WorldXY
        fStartValue = nc_In.PointAtStart.Z
        fEndValue = nc_In.PointAtEnd.Z
    else:
        raise ValueError("{} not a value symmetry plane description.".format(
            sMirrorOpt))

    if abs(fStartValue) <= fDistTol:
        ncs_WIP = createSymmetricPair_Curve(nc_In, plane_sym, False)
        if not ncs_WIP: return
        if abs(fEndValue) <= fDistTol:
            ncs_WIP[1].Dispose()
            ncs_WIP = createSymmetricPair_Curve(ncs_WIP[0], plane_sym, True)
            if not ncs_WIP: return
    elif abs(fEndValue) <= fDistTol:
        ncs_WIP = createSymmetricPair_Curve(nc_In, plane_sym, True)
    else:
        return

    if not sc.doc.Objects.Replace(objref_In, ncs_WIP[0]):
        return

    gMirror = sc.doc.Objects.AddCurve(ncs_WIP[1], objref_In.Object().Attributes)
    if gMirror == gMirror.Empty: return

    return objref_In.ObjectId, gMirror


def main():
    
    objrefs_In = _getInput_WireCrvsAndUntrimmedSrfs()
    if objrefs_In is None: return

    sMirrorOpt = _getInput_SymmetryPlane()
    if sMirrorOpt is None: return
    if sMirrorOpt not in ('X', 'Y', 'Z'):
        raise ValueError("{} not a value symmetry plane description.".format(
                sMirrorOpt))
    
    fDistTol = Opts.values['fDistTol']
    
    iCt_Results_Bs = iCt_Results_Cs = 0
    gBs_In = []
    gCs_In = []
    
    for objref_In in objrefs_In:
        if isinstance(objref_In.Object(), rd.BrepObject):
            gBs_In.append(objref_In.ObjectId)
            rc = processBrepObject(objref_In, sMirrorOpt, fDistTol)
            if rc is None: continue
            iCt_Results_Bs += 1
        elif isinstance(objref_In.Object(), rd.CurveObject):
            gCs_In.append(objref_In.ObjectId)
            rc = processCurveObject(objref_In, sMirrorOpt, fDistTol)
            if rc is None: continue
            iCt_Results_Cs += 1


    sc.doc.Views.Redraw()


    if gBs_In and gCs_In:
        print("{} out of {} curves and {} out of {} surfaces made into symmetric pairs.".format(
            iCt_Results_Cs,
            len(gCs_In),
            iCt_Results_Bs,
            len(gBs_In),
            ))
    elif gBs_In:
        print("{} out of {} surfaces made into symmetric pairs.".format(
            iCt_Results_Bs, len(gBs_In)))
    else:
        print("{} out of {} curves made into symmetric pairs.".format(
            iCt_Results_Cs, len(gCs_In)))


if __name__ == '__main__': main()