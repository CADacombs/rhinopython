"""
200521: Created, starting with another script.
200729, 220317: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xCurve


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


    key = 'bKeepExterior'; keys.append(key)
    values[key] = False
    names[key] = "Keep"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Interior', 'Exterior')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bReplace'; keys.append(key)
    values[key] = True
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


def getInput_Curves():
    """
    Get wires with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select wire curves")
    go.SetCommandPromptDefault("Enter for all normal wires")

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve
    
    go.AcceptNothing(True)
    
    go.AcceptNumber(enable=True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False)
    go.EnableUnselectObjectsOnExit(False)
    
    idxs_Opts = {}

    while True:
        key = 'bKeepExterior'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bReplace'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            return
        elif res == ri.GetResult.Nothing:
            rdCrvs = []
            rgEdges = []
            settings = rd.ObjectEnumeratorSettings()
            settings.NormalObjects = True
            settings.LockedObjects = False
            for rdObj in sc.doc.Objects.GetObjectList(settings):
                if rdObj.ObjectType == rd.ObjectType.Curve:
                    rdCrvs.append(rdObj)
            return tuple([rdCrvs + rgEdges] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fTolerance'].CurrentValue = go.Number()
        
        if Opts.riOpts['fTolerance'].CurrentValue < 0.0:
            Opts.riOpts['fTolerance'].CurrentValue = Opts.riOpts['fTolerance'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInput_Face():
    """Get Brepface with optional input."""
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(enable=True, acceptZero=True)
    
    idxs_Opts = {}

    while True:
        key = 'bKeepExterior'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bReplace'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            return
        
        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            Opts.riOpts['fTolerance'].CurrentValue = go.Number()
        
        if Opts.riOpts['fTolerance'].CurrentValue < 0.0:
            Opts.riOpts['fTolerance'].CurrentValue = Opts.riOpts['fTolerance'].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def processCurveObjects(rhCrvObjs, rhSrf, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bKeepExterior = getOpt('bKeepExterior')
    fTolerance = getOpt('fTolerance')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getCurve(rhObj):
        if isinstance(rhObj, rd.CurveObject):
            return rhObj, rhObj.CurveGeometry

        if isinstance(rhObj, rg.Curve):
            return None, rhObj

        if isinstance(rhObj, rg.GeometryBase):
            rdObj = None
            rgObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
            rgObj = rhObj.Geometry()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
            rgObj = rdObj.Geometry
        else:
            return

        if isinstance(rgObj, rg.Curve):
            return rdObj, rgObj


    def getSurface(rhSrf):
        """
        Parameters:
            rhFace
        Returns:
            rg.BrepFace or rg.Surface
        """
        
        if isinstance(rhSrf, rg.BrepFace):
            return rhSrf

        if isinstance(rhSrf, rg.Surface):
            return rhSrf

        if isinstance(rhSrf, rg.Brep):
            rgBrep = rhSrf
            if rgBrep.Faces.Count == 1:
                return rgBrep.Faces[0]
            else:
                return

        if isinstance(rhSrf, rd.ObjRef):
            objref = rhSrf
            compIdxType = objref.GeometryComponentIndex.ComponentIndexType

            if compIdxType == rg.ComponentIndexType.BrepFace:
                return objref.Face()

            if compIdxType == rg.ComponentIndexType.InvalidType:
                rgBrep = objref.Brep()
                if rgBrep.Faces.Count == 1:
                    return rgBrep.Faces[0]
                else:
                    return

        if isinstance(rhSrf, Guid):
            rdObj = sc.doc.Objects.FindId(rhSrf) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhSrf)
            if isinstance(rdObj, rd.BrepObject):
                brep = rdObj.Geometry
                if rhSrf.Faces.Count == 1:
                    return rhSrf.Faces[0]
                else:
                    return

        srf = None
        if isinstance(geom, rg.BrepFace):
            srf = geom.UnderlyingSurface()
        elif isinstance(geom, rg.Surface):
            srf = geom

        return srf


    rdCrvs_In = []
    rgCrvs_In = []
    for o in rhCrvObjs:
        rdCrv_In, rgCrv_In = getCurve(o)
        if rdCrv_In:
            rdCrvs_In.append(rdCrv_In)
            rgCrvs_In.append(rgCrv_In)

    rgSrf_In = getSurface(rhSrf)

    gCrvs_Out = []
    
    iCt_Split_AllCrv_In = 0
    iCt_CrvsSplit_NoError = 0
    
    for i, rgCrv in enumerate(rgCrvs_In):
        
        attrs = rdCrvs_In[i].Attributes if bReplace else None

        rc = xCurve.splitCurvesWithSurfaceEdges(
            [rgCrv],
            rgSrf_In=rgSrf_In,
            fTolerance=fTolerance)

        if not rc: continue
        
        (
            cs_Interior,
            cs_Boundary,
            cs_Exterior
            ) = rc

        for c in cs_Boundary: c.Dispose()

        iCt_Added_ThisCrv_In = 0
        
        if bKeepExterior:
            crvs_toAdd = cs_Exterior
            for c in cs_Interior: c.Dispose()
        else:
            crvs_toAdd = cs_Interior
            for c in cs_Exterior: c.Dispose()
        
        if not crvs_toAdd: continue
        
        for c in crvs_toAdd:
            gC = sc.doc.Objects.AddCurve(c, attributes=attrs)
            if gC != Guid.Empty:
                gCrvs_Out.append(gC)
                iCt_Added_ThisCrv_In += 1
        
        if iCt_Added_ThisCrv_In == len(crvs_toAdd):
            if bReplace:
                sc.doc.Objects.Delete(objectId=rdCrvs_In[i].Id, quiet=False)
            iCt_Split_AllCrv_In += iCt_Added_ThisCrv_In
            iCt_CrvsSplit_NoError += 1
        else:
            print "Added {} out of {} curves.  Original curve was not removed.".format(
                iCt_Added_ThisCrv_In, len(crvs_toAdd))

        for c in crvs_toAdd: c.Dispose()

    if iCt_CrvsSplit_NoError:
        print "Trimmed {} curves into {} curves.".format(
            iCt_CrvsSplit_NoError, iCt_Split_AllCrv_In)
    else:
        print "No curves were split."

    return gCrvs_Out


def main():


    rc = getInput_Curves()
    if rc is None: return
    rhObjects_Wires = rc[0]

    rc = getInput_Face()
    if rc is None: return
    objrefs_Face = rc[0]

    #sc.doc.Objects.UnselectAll()

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gCrvs_Ret = processCurveObjects(
        rhCrvObjs=rhObjects_Wires,
        rhSrf=objrefs_Face,
        )

    if gCrvs_Ret:
        sc.doc.Objects.UnselectAll()
        for gC in gCrvs_Ret:
            sc.doc.Objects.Select(objectId=gC)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()