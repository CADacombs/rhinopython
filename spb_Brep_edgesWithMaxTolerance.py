"""
An alternative to _testMarkOTEdges, this script has more options.
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
191102-04: Created.
...
191107: Bug fix.  Added printed feedback.
...
231110: Now dots the edges with maximum tolerance per all input breps regardless of tolerance.
240802: Modified an option default. Modified some printed output.
241124: Added and modified some options.
241125: Refactored.
241128: Modified some printed output.
250223: Now, only the absolute maximum is the target, whether it be the single maximum of all or each of the input breps.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bPerBrep'; keys.append(key)
    values[key] = False
    names[key] = 'MaxTol'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'OfAllBreps', 'PerBrep')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOnlyOutOfSearchTol'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSearchTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDupCrv'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDot'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHeight'; keys.append(key)
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
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
        _debug = sc.sticky
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
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fSearchTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def getAllNormalBreps():
    oes = rd.ObjectEnumeratorSettings()
    oes.NormalObjects = True
    oes.LockedObjects = False # Default is True.
    oes.IncludeLights = False
    oes.IncludeGrips = False
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))


def getInput():
    """
    Get Breps with optional input
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False

    go.AcceptNothing(True)

    #sc.doc.Views.Redraw()

    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('bPerBrep')
        addOption('bOnlyOutOfSearchTol')
        if Opts.values['bOnlyOutOfSearchTol']:
            addOption('fSearchTol')
        addOption('bDupCrv')
        addOption('bDot')
        if Opts.values['bDot']:
            addOption('iDotHeight')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return getAllNormalBreps()

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            if 'fSearchTol' not in idxs_Opts:
                print("Numerical input was not applied.")
                continue
            key = 'fSearchTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.Brep):
        return rhObj
    elif isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.BrepObject):
        return rhObj.BrepGeometry
    elif isinstance(rhObj, rd.ObjRef):
        #print(rhObj.GeometryComponentIndex.ComponentIndexType
        geom = rhObj.Geometry()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    if isinstance(geom, rg.Brep):
        return geom


def formatDistance(fDistance):
    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 10**(-sc.doc.ModelDistanceDisplayPrecision):
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def _addDot(edge, iDotHeight, attr_Out):
    # Dot edge at midpoint.
    ts = edge.DivideByCount(2, includeEnds=False)
    pt = edge.PointAtStart if ts is None else edge.PointAt(ts[0])

    rgDot = rg.TextDot(
        text='e[{}]\nTol={}'.format(
            edge.EdgeIndex, formatDistance(edge.Tolerance)),
        location=pt)
    rgDot.FontHeight = iDotHeight
    return sc.doc.Objects.AddTextDot(rgDot, attr_Out)


def processBreps(rhBreps, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bPerBrep = getOpt('bPerBrep')
    bOnlyOutOfSearchTol = getOpt('bOnlyOutOfSearchTol')
    fSearchTol = getOpt('fSearchTol')
    bDupCrv = getOpt('bDupCrv')
    bDot = getOpt('bDot')
    iDotHeight = getOpt('iDotHeight')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if bDupCrv or bDot:
        attr_Out = rd.ObjectAttributes()
        attr_Out.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        attr_Out.ColorSource = rd.ObjectColorSource.ColorFromObject
        attr_Out.ObjectColor = Color.FromArgb(255,0,0)

    rgB_maxTol_AllBs = []
    fTols_Max_AllBs = []
    idxEs_Maxima_AllBs = []

    iCt_Found = 0

    sLogs = []

    gCrvs_Out = []
    gDots_Out = []

    for rhBrep in rhBreps:
        rgB = coerceBrep(rhBrep)

        if not rgB.IsValid:
            sLogs.append("Warning: Brep is invalid but will still be processed. ")

        tols = [edge.Tolerance for edge in rgB.Edges]
        maxTol_ThisB = max(tols)

        if (
            bOnlyOutOfSearchTol and
            (round(maxTol_ThisB, sc.doc.ModelDistanceDisplayPrecision) <=
                round(fSearchTol, sc.doc.ModelDistanceDisplayPrecision))
        ):
            if len(rhBreps) == 1:
                sLogs.append("Maximum Edge.Tolerance is {} ({} less than search tol.).".format(
                        formatDistance(maxTol_ThisB),
                        formatDistance(fSearchTol-maxTol_ThisB)))
            else:
                sLogs.append(
                    "Maximum Edge.Tolerance of the brep is <= {} (search tol.).".format(
                        fSearchTol))
            continue # to next brep.

        fTols_Max_AllBs.append(maxTol_ThisB)

        rgB_maxTol_AllBs.append(rgB)

        idxEs_Maxima = []
        for edge in rgB.Edges:
            if round(edge.Tolerance, sc.doc.ModelDistanceDisplayPrecision) == round(maxTol_ThisB, sc.doc.ModelDistanceDisplayPrecision):
                idxEs_Maxima.append(edge.EdgeIndex)

        idxEs_Maxima_AllBs.append(idxEs_Maxima)

        bOutOfTolFound = maxTol_ThisB > fSearchTol

        if not bOutOfTolFound:
            if len(rhBreps) == 1:
                sLogs.append("Maximum Edge.Tolerance is {} ({} less than search tol.).".format(
                        formatDistance(maxTol_ThisB),
                        formatDistance(fSearchTol-maxTol_ThisB)))
            else:
                sLogs.append(
                    "Maximum Edge.Tolerance is <= {} (search tol.).".format(
                        fSearchTol))

            continue # to next brep.

    if sLogs:
        s = "\n".join("[{}] {}".format(sLogs.count(sLog), sLog) for sLog in set(sLogs))
    else:
        s = ""

    if bPerBrep:
        iCt_Found = 0
        for rgB, idxs_E in zip(rgB_maxTol_AllBs, idxEs_Maxima_AllBs):
            for idx_E in idxs_E:
                edge = rgB.Edges[idx_E]
                iCt_Found += 1
                if bDupCrv:
                    gCrv_Out = sc.doc.Objects.AddCurve(edge, attr_Out)
                    if gCrv_Out != gCrv_Out.Empty:
                        gCrvs_Out.append(gCrv_Out)
                if bDot:
                    gDot_Out = _addDot(edge, iDotHeight, attr_Out)
                    if gDot_Out != gDot_Out.Empty:
                        gDots_Out.append(gDot_Out)
        s += " {} edges within {} of max. edge tol. per Brep.".format(
            iCt_Found, 10**-sc.doc.ModelDistanceDisplayPrecision)
    else:
        fTol_Max_All = max(fTols_Max_AllBs)
        iCt_Found = 0
        for rgB, idxs_E in zip(rgB_maxTol_AllBs, idxEs_Maxima_AllBs):
            for idx_E in idxs_E:
                edge = rgB.Edges[idx_E]
                if round(edge.Tolerance, sc.doc.ModelDistanceDisplayPrecision) >= round(fTol_Max_All, sc.doc.ModelDistanceDisplayPrecision):
                    iCt_Found += 1
                    if bDupCrv:
                        gCrv_Out = sc.doc.Objects.AddCurve(edge, attr_Out)
                        if gCrv_Out != gCrv_Out.Empty:
                            gCrvs_Out.append(gCrv_Out)
                    if bDot:
                        gDot_Out = _addDot(edge, iDotHeight, attr_Out)
                        if gDot_Out != gDot_Out.Empty:
                            gDots_Out.append(gDot_Out)

        if iCt_Found:
            s += " {} edges found within {} of {}.".format(
                iCt_Found, 10.0**-sc.doc.ModelDistanceDisplayPrecision, formatDistance(fTol_Max_All))
        else:
            s += " Max. tol. of all edges is {}.".format(
                formatDistance(max(fTols_Max_AllBs)))

    if bDupCrv:
        s += " {} curves added.".format(len(gCrvs_Out))
    if bDot:
        s += " {} dots added.".format(len(gDots_Out))

    print(s)


def main():

    rhObjs = getInput()
    if not rhObjs: return

    bPerBrep = Opts.values['bPerBrep']
    bOnlyOutOfSearchTol = Opts.values['bOnlyOutOfSearchTol']
    fSearchTol = Opts.values['fSearchTol']
    bDupCrv = Opts.values['bDupCrv']
    bDot = Opts.values['bDot']
    iDotHeight = Opts.values['iDotHeight']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    processBreps(
        rhObjs,
        bPerBrep=bPerBrep,
        bOnlyOutOfSearchTol=bOnlyOutOfSearchTol,
        fSearchTol=fSearchTol,
        bDupCrv=bDupCrv,
        bDot=bDot,
        iDotHeight=iDotHeight,
        bEcho=bEcho,
        bDebug=bDebug,
        )


    sc.doc.Views.Redraw()


if __name__ == '__main__': main()