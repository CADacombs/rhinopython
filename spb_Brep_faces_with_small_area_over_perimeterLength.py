"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250514: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import xBrepObject


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fMaxTol'; keys.append(key)
    #values[key] = max((sc.doc.ModelAbsoluteTolerance)**2.0, 1e-6)
    #values[key] = (1.5 * sc.doc.ModelAbsoluteTolerance)**2.0
    #values[key] = 0.5 * sc.doc.ModelAbsoluteTolerance**2.0 # Triangle.
    values[key] = 0.5 * sc.doc.ModelAbsoluteTolerance
    #names[key] = 'MaxTol'
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

        if key == 'fMaxTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.values[key] = cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            if cls.riOpts[key].CurrentValue <= 1e-12:
                cls.values[key] = cls.riOpts[key].CurrentValue = 1e-12
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


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
    Get breps with optional input
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False
    
    go.AcceptNothing(True)

    go.AcceptNumber(True, acceptZero=True)
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    bPreselectedObjsChecked = False

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        addOption('fMaxTol')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return getAllNormalBreps()

        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fMaxTol'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getFaces(rgBrep0, fMaxTol=None, bEcho=None, bDebug=None):
    """
    """
    
    if fMaxTol is None: fMaxTol = Opts.values['fMaxTol']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    idx_rgFaces_Found = []
    fArea_Min_Brep = fArea_Max_Brep = None

    sCP_Base = Rhino.RhinoApp.CommandPrompt

    # Search all faces of brep for slivers.
    for iF, rgFace in enumerate(rgBrep0.Faces):
        if sc.escape_test(throw_exception=False):
            raise Exception("Search canceled by user.")

        Rhino.RhinoApp.CommandPrompt = "{}  Face {} of {}".format(
            sCP_Base,
            iF + 1,
            rgBrep0.Faces.Count)

        rgBrepX_1Face = rgFace.DuplicateFace(duplicateMeshes=True)
        
        area = rgBrepX_1Face.GetArea()
        
        if area is None:
            raise Exception("Area is None.")
        
        # GetArea returns 0 when the area cannot be calculated.
        if area == 0:
            rgBrepX_1Face.Dispose()
            continue
        
        fTotalEdgeLength = sum([edge.GetLength() for edge in rgBrepX_1Face.Edges])
        
        if bDebug:
            sEval = "fMaxTol"; print(sEval, '=', eval(sEval))
            sEval = "area"; print(sEval, '=', eval(sEval))
            sEval = "fTotalEdgeLength"; print(sEval, '=', eval(sEval))
            sEval = "area/fTotalEdgeLength"; print(sEval, '=', eval(sEval))
        
        if (area/fTotalEdgeLength) > fMaxTol:
            continue
        
        idx_rgFaces_Found.append(rgFace.FaceIndex)
        
        if (
                fArea_Min_Brep is None or
                area < fArea_Min_Brep
        ):
            fArea_Min_Brep = area
        if (
                fArea_Max_Brep is None or
                area > fArea_Max_Brep
        ):
            fArea_Max_Brep = area
        
        rgBrepX_1Face.Dispose()

    Rhino.RhinoApp.CommandPrompt = sCP_Base

    return sorted(idx_rgFaces_Found), fArea_Min_Brep, fArea_Max_Brep


def formatDistance(fDistance):
    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def processBrepObjects(rhBreps_In, fMaxTol=None, bEcho=None, bDebug=None):
    """
    """
    
    if fMaxTol is None: fMaxTol = Opts.values['fMaxTol']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']
    
    fArea_Min_All = fArea_Max_All = None
    iCt_Found = 0
    gBreps1_Found_All = [] # Accumulation of duplicated faces (breps)
    
    nBreps_In = len(rhBreps_In)

    idxs_AtTenths = [int(round(0.1*i*nBreps_In,0)) for i in range(10)]

    for iB, rhBrep0 in enumerate(rhBreps_In):
        if sc.escape_test(False):
            print("Search canceled by user.")
            return

        if nBreps_In == 1:
            Rhino.RhinoApp.CommandPrompt = "Analyzing 1 brep"
        elif iB in idxs_AtTenths:
            Rhino.RhinoApp.CommandPrompt = "Analysis at {:d}% of {} breps".format(
                int(100.0 * (iB+1) / nBreps_In), nBreps_In)
        
        # Obtain GUID, RhinoObject, and geometry.
        if isinstance(rhBrep0, Guid):
            gBrep0 = rhBrep0
            rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gBrep0)
            if rdBrep0 is None: continue
        elif isinstance(rhBrep0, rd.ObjRef):
            gBrep0 = rhBrep0.ObjectId
            rdBrep0 = rhBrep0.Object()
        elif isinstance(rhBrep0, rd.RhinoObject):
            rdBrep0 = rhBrep0
            gBrep0 = rhBrep0.Id
        else:
            print("GUID could not be obtained.")
            continue
        rgBrep0 = rdBrep0.Geometry
        if rgBrep0 is None:
            print("Brep geometry could not be obtained for {}.".format(gBrep0))
            continue
        if not rgBrep0.IsValid and bEcho:
            print("Brep {} is invalid and will be skipped because rg.Curve.GetDistancesBetweenCurves may hang.".format(
                gBrep0))
            continue

        rc = getFaces(rgBrep0, fMaxTol=fMaxTol)

        if rc is None: continue
        (
            idx_rgFaces_Found,
            fArea_Min_Brep,
            fArea_Max_Brep,
            ) = rc
    
        if not idx_rgFaces_Found: continue
        
        iCt_Found += len(idx_rgFaces_Found)
    
        # Extract faces from polyface brep.
        rc = xBrepObject.extractFaces(
                gBrep0,
                idx_rgFaces_Found,
                bEcho=bEcho)
        if rc is None: continue

        gBreps1_Found_All.extend(rc[0]) # To be selected later.
        
        if fArea_Min_All is None or fArea_Min_Brep < fArea_Min_All:
            fArea_Min_All = fArea_Min_Brep
        if fArea_Max_All is None or fArea_Max_Brep > fArea_Max_All:
            fArea_Max_All = fArea_Max_Brep
    
    if not iCt_Found:
        print("No faces with area less than {} square {} found.".format(
                formatDistance(fMaxTol),
                sc.doc.GetUnitSystemName(
                        modelUnits=True,
                        capitalize=False,
                        singular=False,
                        abbreviate=False
                    )
                ))
        return
    
    sc.doc.Objects.UnselectAll()
    ct_Selected = rs.SelectObjects(object_ids=gBreps1_Found_All)
    if len(gBreps1_Found_All) == ct_Selected:
        s = "{} monoface brep(s) selected.".format(ct_Selected)
    else:
        s  = "{} faces found,".format(
                len(gBreps1_Found_All))
        s += " but only {} selected.".format(ct_Selected)

    if iCt_Found == 1:
        s += "\nArea: {} square {}".format(
            formatDistance(fArea_Min_All),
            sc.doc.GetUnitSystemName(
                modelUnits=True,
                capitalize=False,
                singular=False,
                abbreviate=False
                )
            )
    else:
        s += "\nRange of areas found:"
        s += " [{}, {}] square {}".format(
            formatDistance(fArea_Min_All),
            formatDistance(fArea_Max_All),
            sc.doc.GetUnitSystemName(
                modelUnits=True,
                capitalize=False,
                singular=False,
                abbreviate=False
                )
            )
    print(s)


def main():

    rhObjs = getInput()
    if rhObjs is None: return

    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    processBrepObjects(rhBreps_In=rhObjs)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()