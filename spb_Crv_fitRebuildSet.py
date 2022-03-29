"""

Input curves and options:
    Curve deviation tolerance
    Whether degree is fixed.  If so:
        Degree for rebuild
    Whether preserve tangents
    Whether replace curve
If degree is fixed:
    Iterate from lowest controls points until tolerance is achieved.
Else:
    If curve is not closed: Try to simplify starting with a line (degree 1 with 2 CP's).
    At degrees 2, 3, & 5, iterate from lowest controls points until tolerance is achieved
    Choose the curve with the least amount of control points
"""

"""
190417: This script started as a split from another.
191019-200121: Import-related updates.
200210: Bug fix.  Import-related update.
210316: Bug fix for BrepEdge support.
220328: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

import spb_Crv_fitRebuild


sOpts = (
        'fDevTol',
        'iDegree',
        'bPreserveTan',
        'iMaxCpCt',
        'bReplace',
        'bEcho',
        'bDebug',
)

class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    idxOpt = {}
    optIndices = {}
    stickyKeys = {}
    
    @classmethod
    def init(cls):
        
        key = 'fDevTol'
        cls.keys.append(key)
        cls.values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
        cls.names[key] = 'DevTol'
        cls.riOpts[key] = ri.Custom.OptionDouble(initialValue=cls.values[key], setLowerLimit=True, limit=0.0)
        cls.stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
        
        key = 'iDegree'
        cls.keys.append(key)
        cls.values[key] = 3
        cls.names[key] = key[1:]
        cls.riOpts[key] = ri.Custom.OptionInteger(initialValue=cls.values[key])
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bPreserveTan'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = 'PreserveEndTangents'
        cls.riOpts[key] = ri.Custom.OptionToggle(initialValue=cls.values[key], offValue='No', onValue='Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'iMaxCpCt'
        cls.keys.append(key)
        cls.values[key] = 20
        cls.names[key] = 'MaxCpCt'
        cls.riOpts[key] = ri.Custom.OptionInteger(initialValue=cls.values[key], setLowerLimit=True, limit=2)
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bOnlyReplaceUniformWithLessCps'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = 'OnlyReplaceUniformCrvsWithLessCps'
        cls.riOpts[key] = ri.Custom.OptionToggle(initialValue=cls.values[key], offValue='No', onValue='Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bReplace'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = 'Action'
        cls.riOpts[key] = ri.Custom.OptionToggle(initialValue=cls.values[key], offValue='Add', onValue='Replace')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bEcho'
        cls.keys.append(key)
        cls.values[key] = True
        cls.names[key] = 'Echo'
        cls.riOpts[key] = ri.Custom.OptionToggle(initialValue=cls.values[key], offValue='No', onValue='Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        key = 'bDebug'
        cls.keys.append(key)
        cls.values[key] = False
        cls.names[key] = 'Debug'
        cls.riOpts[key] = ri.Custom.OptionToggle(initialValue=cls.values[key], offValue='No', onValue='Yes')
        cls.stickyKeys[key] = '{}({})'.format(key, __file__)
        
        # Load sticky.
        for key in cls.stickyKeys:
            if sc.sticky.has_key(cls.stickyKeys[key]):
                if key in cls.riOpts:
                    cls.values[key] = cls.riOpts[key].CurrentValue = sc.sticky[cls.stickyKeys[key]]
                else:
                    # For iDegs.
                    cls.values[key] = sc.sticky[cls.stickyKeys[key]]
                    for i in 1,2,3,5:
                        cls.riOpts['bDeg{}'.format(i)].CurrentValue = (i in cls.values[key])
    
    @classmethod
    def addOptions(cls, go, sInclude=None):
        
        go.ClearCommandOptions()
        sInclude = None
        if sInclude is None:
            sInclude = cls.keys[:]
        if not sInclude: return
        pass
        
        def addOptionDouble(key):
            if key in sInclude:
                go.AddOptionDouble(cls.names[key], cls.riOpts[key])
        
        def addOptionInteger(key):
            if key in sInclude:
                go.AddOptionInteger(cls.names[key], cls.riOpts[key])
        
        def addOptionToggle(key):
            if key in sInclude:
                go.AddOptionToggle(cls.names[key], cls.riOpts[key])
        
        cls.setValues()
        
        addOptionDouble('fDevTol')
        addOptionToggle('bPreserveTan')
        addOptionInteger('iDegree')
        addOptionInteger('iMaxCpCt')
        addOptionToggle('bReplace')
        addOptionToggle('bEcho')
        addOptionToggle('bDebug')
    
    @classmethod
    def setValues(cls):
        for key in cls.stickyKeys:
            if cls.riOpts[key]:
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
                key
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
    
    @classmethod
    def processInput(cls, go):
        res = go.Result()
        
        key = 'fDevTol'
        if res == ri.GetResult.Number:
            cls.riOpts[key].CurrentValue = go.Number()
        if cls.riOpts[key].CurrentValue < 0.0:
            cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
        
        cls.setValues()
        cls.saveSticky()
        
        # Clear and add options regardless if a number was entered or options were modified in another way.
        Opts.addOptions(go)
Opts.init()


def getInput():
    """Get curves with optional input."""
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select curves")
    
    go.GeometryFilter = rd.ObjectType.Curve
    
    go.AcceptNumber(True, True)
    
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    Opts.addOptions(go)
    
    s  = "For MaxCpCt, \"0\" will use the curve's current control point count."
    s += "  Other values are the absolute maximum to allow."
    print s
    
    while True:
        if sc.escape_test(False): break
        
        res = go.GetMultiple(minimumNumber=2, maximumNumber=0)
        
        if res == ri.GetResult.Cancel: return # Esc key was pressed.  (Don't use Rhino.Commands.Result.Cancel)
        
        if res == ri.GetResult.Object:
            break
        elif res == ri.GetResult.Cancel:
            return
        else:
            # An option was selected or a number was entered.
            Opts.processInput(go)
    
    objrefs = go.Objects()
    
    go.Dispose()
    
    return (
            objrefs,
            Opts.values['bEcho'],
            Opts.values['bDebug'],
    )


def processCurves(curvesAndEdges0, fDevTol=None, iDegree=None, bPreserveTan=None, iMaxCpCt=None, bReplace=None, bEcho=None, bDebug=None):
    """
    curvesAndEdges0 = [ObjRefs or GUIDs_of_CurveObjects]
    """
    
    if fDevTol is None: fDevTol=Opts.values['fDevTol']
    if iDegree is None: iDegree=Opts.values['iDegree']
    if bPreserveTan is None: bPreserveTan=Opts.values['bPreserveTan']
    if iMaxCpCt is None: iMaxCpCt=Opts.values['iMaxCpCt']
    if bReplace is None: bReplace=Opts.values['bReplace']
    if bEcho is None: bEcho=Opts.values['bEcho']
    if bDebug is None: bDebug=Opts.values['bDebug']
    
    s  = "Maximum CP counts to allow: "
    s += "(Same as input curve)" if iMaxCpCt == 0 else str(iMaxCpCt)
    print s
    
    gObjs_In = []
    gBreps_In = []
    bIs_wire = []
    ncs_Start = []
    ncs1 = []
    devs1 = []
    
    for curveOrEdge0 in curvesAndEdges0:
        rdObj_In = rs.coercerhinoobject(curveOrEdge0)
        if isinstance(rdObj_In, rd.CurveObject):
            bIs_wire.append(True)
        elif isinstance(rdObj_In, rd.BrepObject):
            bIs_wire.append(False)
        else:
            raise ValueError("Invalid input: {}".format(rdObj_In))
        
        gObjs_In.append(rdObj_In.Id)
        
        ncs_Start.append(rs.coercecurve(curveOrEdge0).ToNurbsCurve())
    
    for i, nc0 in enumerate(ncs_Start):
        nc1, dev1, sLog = spb_Crv_fitRebuild.rebuildCurve(
                rgCurve0=nc0,
                fDevTol=fDevTol,
                iDegree=iDegree,
                bPreserveEndTans=bPreserveTan,
                bFurtherTranslateCps=False,
                iMinCpCt=None,
                iMaxCpCt=iMaxCpCt,
                bDebug=bDebug,
        )
        if nc1 is None:
            if bEcho:
                print "Dev:{}  Log:{} for {}".format(
                    dev1,
                    sLog,
                    "Edge" if gObjs_In[i] is None else gObjs_In[i])
            return
        
        ncs1.append(nc1)
        devs1.append(dev1)
    
    if iMaxCpCt == 0:
        iMaxCpCt = max([nc0.Points.Count for nc0 in ncs_Start]) # Greatest control point count of all curves.
    
    iMax_cpCt_ncs1 = max([nc1.Points.Count for nc1 in ncs1])
    if iMax_cpCt_ncs1 > iMaxCpCt:
        s  = "Required control point count to create all curves within tolerance:"
        s += " {}".format(iMax_cpCt_ncs1)
        for c in ncs_Start: c.Dispose()
        for c in ncs1: c.Dispose()
        return

    if min([nc1.Points.Count for nc1 in ncs1]) == iMax_cpCt_ncs1:
        # NurbsCurve are already at same control point count.
        ncs2 = [nc1.DuplicateCurve() for nc1 in ncs1]
        devs2 = devs1[:]
    else:
        pointCount = iMax_cpCt_ncs1
        while pointCount <= iMaxCpCt:
            sc.escape_test()
            ncs2 = []
            devs2 = []
            for nc0 in ncs_Start:
                nc2 = nc0.Rebuild(
                        pointCount=pointCount,
                        degree=iDegree,
                        preserveTangents=bPreserveTan)

                if nc2 is None:
                    raise ValueError("NurbsCurve.Rebuild returned None.")
            
                ncs2.append(nc2)

                rc = rg.Curve.GetDistancesBetweenCurves(nc0, nc2, 0.1*fDevTol)
                if rc[0]:
                    dev2 = rc[1]
                    if dev2 > fDevTol:
                        for c in ncs2: c.Dispose()
                        pointCount += 1
                        break # for loop back to while to try to Rebuild all curves to new pointCount.
                    devs2.append(dev2)
            else:
                # All curves were rebuilt within fDevTol.
                break


    s = ""
    
    ct_Replaced = 0
    ct_Added = 0
    
    for gObj_In, nc2, bIsWire in zip(gObjs_In, ncs2, bIs_wire):
        if bReplace and bIsWire:
            if sc.doc.Objects.Replace(gObj_In, nc2):
                ct_Replaced += 1
        else:
            gCrv1 = sc.doc.Objects.AddCurve(nc2)
            if gCrv1 != Guid.Empty:
                ct_Added += 1

    if bReplace and any(bIs_wire):
        if ct_Replaced == 0:
            s += "No wires were replaced.".format(ct_Replaced)
        elif ct_Replaced == len(gObjs_In):
            s += "All {} wires were replaced.".format(ct_Replaced)
        else:
            s +=  "Only {} out of {} wires were replaced.".format(ct_Replaced, len(gObjs_In))

    if ct_Added == len(gObjs_In):
        s += "All {} curves were added.".format(ct_Added)
    else:
        s += "Only {} out of {} curves were added.".format(ct_Added, len(gObjs_In))

    s += "  Degree: {}".format(nc2.Degree)
    s += "  CpCt: {}".format(nc2.Points.Count)
    s += "  Maximum deviation: {0:.{1}f}".format(
            max(devs2), sc.doc.ModelDistanceDisplayPrecision)
    
    print s


def main():
    
    rc = getInput()
    if rc is None: return
    
    (
            objrefs,
            bEcho,
            bDebug,
    ) = rc
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")
    
    if bDebug:
        reload(spb_Crv_fitRebuild)
    else:
        sc.doc.Views.RedrawEnabled = False
    
    processCurves(curvesAndEdges0=objrefs)
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()