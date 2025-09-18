"""
This script selects curves (wires or all curves (including edges)) and optionally adds curves of naked edges that lie
on selected faced.
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
200520-23: Created, starting with another script.
...
220317: Moved main function to a main library module.
230106: Bug fix.
230108: Fixed selection routine.
250917: Added option to filter edges.

TODO:
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xBrepFace
import xCurve


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bUnderlyingSrf'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPartiallyOnSrf'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fSamplingResolution'; keys.append(key)
    values[key] = 100.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bIncludeEdges'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddCrvsForEdges'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_Faces():
    """
    Get Brepfaces with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select faces")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(True, acceptZero=True)
    
    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bUnderlyingSrf')
        addOption('fTolerance')
        addOption('bPartiallyOnSrf')
        addOption('fSamplingResolution')
        addOption('bIncludeEdges')
        if Opts.values['bIncludeEdges']:
            addOption('bAddCrvsForEdges')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fTolerance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_Curves(objrefs_Face):
    """
    Get curves or edges with optional input.
    """


    def getBrep(rhBrep):
        if isinstance(rhBrep, rg.Brep):
            return None, rhBrep
        elif isinstance(rhBrep, rg.GeometryBase):
            rdObj = None
            rgObj = rhBrep
        elif isinstance(rhBrep, rd.ObjRef):
            rdObj = rhBrep.Object()
            rgObj = rhBrep.Geometry()
        elif isinstance(rhBrep, Guid):
            rdObj = sc.doc.Objects.FindId(rhBrep) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhBrep)
            rgObj = rdObj.Geometry
        else:
            return

        if isinstance(rgObj, (rg.Brep, rg.BrepFace)):
            return rdObj, rgObj


    idxs_EdgesOfFaceToSplit = []
    for objref_Face in objrefs_Face:
        rdBrep_withFaceToCheck, rgBrep_withFaceToCheck = getBrep(objref_Face)
        gBrep_withFaceToCheck = rdBrep_withFaceToCheck.Id
        idxs_EdgesOfFaceToSplit.extend(objref_Face.Face().AdjacentEdges())


    go = ri.Custom.GetObject()

    go.GeometryFilter = rd.ObjectType.Curve

    def notEdgeOfFaceToCheck(rdObj, geom, compIdx):
        #print rdObj, geom, compIdx
        if isinstance(rdObj, rd.BrepObject) and rdObj.Id == gBrep_withFaceToCheck:
            if geom.EdgeIndex in idxs_EdgesOfFaceToSplit:
                print("An edge of a face to split was picked and will not be used.")
                return False
        return True
    go.SetCustomGeometryFilter(notEdgeOfFaceToCheck)


    go.AcceptNothing(True)
    
    go.AcceptNumber(True, acceptZero=True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False)
    go.EnableUnselectObjectsOnExit(False)
    
    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    #print(go.GeometryAttributeFilter)

    while True:
        bIncludeEdges = Opts.values['bIncludeEdges']
        if bIncludeEdges:
            go.SetCommandPrompt("Select wires or edges")
            go.SetCommandPromptDefault("Enter for all normal wires and brep naked edges")
            go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.AcceptAllAttributes
        else:
            go.SetCommandPrompt("Select wires")
            go.SetCommandPromptDefault("Enter for all normal wires")
            go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bUnderlyingSrf')
        addOption('fTolerance')
        addOption('bPartiallyOnSrf')
        addOption('fSamplingResolution')
        if Opts.values['bIncludeEdges']:
            addOption('bAddCrvsForEdges')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return []

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fTolerance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _sortFacesByBrep(objrefs):
    """
    Parameters:
        list(objrefs)
    Returns:
        list(Brep GUIDs)
        list(lists(integers of Face indices) per brep)
    """

    gBs = []
    rdBs = []
    idxs_Fs_perB = []

    for o in objrefs:
        gB = o.ObjectId
        rdB = o.Object()
        rgB = o.Brep()
        
        if not rgB.IsValid:
            print("Brep {} is invalid.  Fix first.".format(gB))
            rgB.Dispose()
            continue
        
        idx_CompIdx = o.GeometryComponentIndex.Index
        if idx_CompIdx == -1:
            if gB in gBs:
                idxs_Fs_perB[gBs.index(gB)] = range(rgB.Faces.Count)
                continue

            gBs.append(gB)
            rdBs.append(rdB)
            idxs_Fs_perB.append(range(rgB.Faces.Count))
        else:
            rgF = o.Face()
            if gB in gBs:
                if rgF.FaceIndex in idxs_Fs_perB[gBs.index(gB)]:
                    continue
                idxs_Fs_perB[gBs.index(gB)].append(rgF.FaceIndex)
                continue

            gBs.append(gB)
            rdBs.append(rdB)
            idxs_Fs_perB.append([rgF.FaceIndex])

    return rdBs, idxs_Fs_perB


def _getDataForAllNormalCrvs(bIncludeEdges=False):
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False

    rdOs_Out = []
    rgCs_Out = []

    if bIncludeEdges:
        oes.ObjectTypeFilter = rd.ObjectType.Brep | rd.ObjectType.Curve
        for rdO in sc.doc.Objects.GetObjectList(oes):
            if isinstance(rdO, rd.BrepObject):
                rgB = rdO.BrepGeometry
                for rgE in rgB.Edges:
                    if rgE.Valence == rg.EdgeAdjacency.Naked:
                        rdOs_Out.append(rdO)
                        rgCs_Out.append(rgE)
            else:
                rdOs_Out.append(rdO)
                rgCs_Out.append(rdO.CurveGeometry)
    else:
        oes.ObjectTypeFilter = rd.ObjectType.Curve

        for rdO in sc.doc.Objects.GetObjectList(oes):
            rdOs_Out.append(rdO)
            rgCs_Out.append(rdO.CurveGeometry)

    return rdOs_Out, rgCs_Out


def main():

    objrefs_Face = getInput_Faces()
    if objrefs_Face is None: return

    sc.doc.Objects.UnselectAll()

    objrefs_CurvesOrEdges = getInput_Curves(objrefs_Face)
    if objrefs_CurvesOrEdges is None: return

    sc.doc.Objects.UnselectAll()

    bUnderlyingSrf = Opts.values['bUnderlyingSrf']
    fTolerance = Opts.values['fTolerance']
    bPartiallyOnSrf = Opts.values['bPartiallyOnSrf']
    fSamplingResolution = Opts.values['fSamplingResolution']
    bIncludeEdges = Opts.values['bIncludeEdges']
    bAddCrvsForEdges = Opts.values['bAddCrvsForEdges']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    Rhino.RhinoApp.SetCommandPrompt("Working ...")


    #rdBs, idxs_Fs_perB = _sortFacesByBrep(objrefs_Face)
    #if not rdBs: return

    rgSs = []
    for objref_F in objrefs_Face:
        rgF = objref_F.Face()
        rgS = rgF.UnderlyingSurface() if bUnderlyingSrf else rgF
        rgSs.append(rgS)


    if len(objrefs_CurvesOrEdges) == 0:
        rdCsOrBs, rgCs = _getDataForAllNormalCrvs(bIncludeEdges=bIncludeEdges)
        if len(rdCsOrBs) == 0:
            print("No objects with curves.")
            return
    else:
        rdCsOrBs = []
        rgCs = []
        for objref in objrefs_CurvesOrEdges:
            rdC, rgC = objref.Object(), objref.Curve()
            if rdC:
                rdCsOrBs.append(rdC)
                rgCs.append(rgC)


    iCBs_Found = []

    for iS, rgS in enumerate(rgSs):

        for iCB, (rdC_or_EB, rgC) in enumerate(zip(rdCsOrBs, rgCs)):

            if rdC_or_EB.Id in iCBs_Found:
                continue

            rc = xCurve.filterCurvesOnSurface(
                rgC,
                rgS,
                fSamplingResolution=fSamplingResolution,
                fTolerance=fTolerance,
                bDebug=bDebug)

            if not rc or (len(rc[0]) + len(rc[1]) == 0): continue

            cs_AllOnSrf, cs_PartiallyOnSrf = rc

            if cs_AllOnSrf:
                iCBs_Found.append(iCB)
                continue

            if bPartiallyOnSrf and cs_PartiallyOnSrf:
                iCBs_Found.append(iCB)
                continue


    if len(iCBs_Found) == 0:
        print("No curves found.")
        return


    sOuts = []

    sOuts.append("{} curves found.".format(len(iCBs_Found)))


    gCs_Out = []

    for iCBs in iCBs_Found:
        rdC_or_B = rdCsOrBs[iCBs]
        if isinstance(rdC_or_B, rd.CurveObject):
            rdC_or_B.Select(rdC_or_B)
        elif bAddCrvsForEdges:
            gC_Out = sc.doc.Objects.AddCurve(rgCs[iCBs])
            if gC_Out != gC_Out.Empty:
                gCs_Out.append(gC_Out)
                sc.doc.Objects.Select(objectId=gC_Out)
        else:
            rgE = rgCs[iCBs]
            compIdx = rg.ComponentIndex(
                type=rg.ComponentIndexType.BrepEdge,
                index=rgE.EdgeIndex)
            rdC_or_B.SelectSubObject(
                componentIndex=compIdx,
                select=True,
                syncHighlight=True,
                persistentSelect=True)

    if gCs_Out:
        sOuts.append("{} curves added.".format(len(gCs_Out)))

    print("  ".join(sOuts))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()