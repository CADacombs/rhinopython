"""
This script retrims faces. This is often helpful in brep repair.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160629-30: Created.
...
200518-28, 0613,15: Refactored.  Simplified main function and removed option to create oppositely trimmed faces.
200618: Added PerFaceColor support for V7.
200619, 0729: Import-related update.
210503: Import-related update.
210630, 250324: Modified an option default value.


When joined curves are used for splitting, edges traversing the joins will be 
merged where possible.

Duplicate naked edges.
*Join all curves.  (Also used to order curves in loops.)
if not closed: Feedback and return fail.
if not allow merging of edges:
    Explode (Curve.DuplicateSegments).
    *Only join short curves with contiguous ones.
Split.
if not all edges are at least min. allowed length:
    Duplicate naked edges.
    *Join only curves >= limit (Curve.JoinCurves(crvs))
    if not closed: Feedback and return fail.
    if not allow merging of edges:
        Explode (Curve.DuplicateSegments).
    Split.
    if not all edges are at least min. allowed length: Feedback and return fail.

TODO: Finish refactoring and replace some local code with imported modules.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xBrepFace
import xBrepObject


class Opts:
    
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


    key = 'bTrimToSegs'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fEdgeLen_Min'; keys.append(key)
    values[key] = 1.8 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'EdgeLengthMinTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fSplitTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bOutputSplitOnFail'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOutputTrimmingCrvs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDotFailures'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotFontHt'; keys.append(key)
    values[key] = 33 if Rhino.RhinoApp.ExeVersion >= 6 else 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], True, 3)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtractNotAdd'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Extract')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
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
        for key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get single face breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select faces to retrim")
    go.SetCommandPromptDefault("All normal breps when none are selected")

    go.GeometryFilter = rd.ObjectType.Surface
    #go.GeometryAttributeFilter = (
    #        ri.Custom.GeometryAttributeFilter.OpenSurface)

    go.AcceptNothing(True)

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.

    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    go.AcceptNumber(enable=True, acceptZero=True)
    
    bPreselectedObjsChecked = False

    idxs_Opts = {}
    
    while True:
        for key in Opts.keys: idxs_Opts[key] = None
        key = 'bTrimToSegs'; idxs_Opts[key] = Opts.riAddOpts[key](go)
        key = 'fEdgeLen_Min'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'fSplitTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bOutputSplitOnFail'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bOutputTrimmingCrvs'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDotFailures'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        if Opts.values['bDotFailures']:
            key = 'iDotFontHt'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bExtractNotAdd'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Nothing:
            oes = rd.ObjectEnumeratorSettings()
            oes.LockedObjects = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep
            rdBs = list(sc.doc.Objects.GetObjectList(oes))
            go.Dispose()
            if len(rdBs) == 0: return

            return rdBs

        # An option was selected or a number was entered.

        key = 'fEdgeLen_Min'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()

        for key in 'fEdgeLen_Min', 'fSplitTol':
            if Opts.riOpts[key].CurrentValue <= 0.0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def _getSortedBrepIdsAndFaces(rhObjs):
    """
    Parameters:
        list(objrefs and/or rd.BrepObjects)
    Returns:
        list(Brep GUIDs)
        list(lists(integers of Face indices) per brep)
    """

    gBreps_In = []
    idxs_Faces_perBrep = []

    for rhObj in rhObjs:

        if isinstance(rhObj, rd.ObjRef):
            objref = rhObj
            gB_In = objref.ObjectId
            rdBrep_In = objref.Object()
            rgBrep_In = rdBrep_In.BrepGeometry
            idx_CompIdx = objref.GeometryComponentIndex.Index
        elif isinstance(rhObj, rd.BrepObject):
            rdBrep_In = rhObj
            gB_In = rdBrep_In.Id
            rgBrep_In = rdBrep_In.BrepGeometry
            idx_CompIdx = None
        else:
            raise Exception("{} passed to _getSortedBrepIdsAndFaces. Needs to be Objref or BrepObject.".format(
                rhObj.GetType().Name))

        #            if not rgBrep_In.IsValid:
        #                print("Brep {} is invalid.  Fix first.".format(gB_In)
        #                rgBrep_In.Dispose()
        #                continue

        if idx_CompIdx in (None, -1):
            if gB_In in gBreps_In:
                idxs_Faces_perBrep[gBreps_In.index(gB_In)] = range(rgBrep_In.Faces.Count)
            else:
                gBreps_In.append(gB_In)
                idxs_Faces_perBrep.append(range(rgBrep_In.Faces.Count))
        else:
            rgFace_Brep0 = objref.Face()
            if gB_In in gBreps_In:
                if rgFace_Brep0 in idxs_Faces_perBrep[gBreps_In.index(gB_In)]:
                    continue
                else:
                    idxs_Faces_perBrep[gBreps_In.index(gB_In)].append(rgFace_Brep0.FaceIndex)
            else:
                gBreps_In.append(gB_In)
                idxs_Faces_perBrep.append([rgFace_Brep0.FaceIndex])

    return gBreps_In, idxs_Faces_perBrep


def processBrepObjects(rhFaces, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bTrimToSegs = getOpt('bTrimToSegs')
    fEdgeLen_Min = getOpt('fEdgeLen_Min')
    fSplitTol = getOpt('fSplitTol')
    bOutputSplitOnFail = getOpt('bOutputSplitOnFail')
    bOutputTrimmingCrvs = getOpt('bOutputTrimmingCrvs')
    bDotFailures = getOpt('bDotFailures')
    iDotFontHt = getOpt('iDotFontHt')
    bExtractNotAdd = getOpt('bExtractNotAdd')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getRhinoObject(rhObj):
        """
        'Deleted objects cannot be found by id.'
        (https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_DocObjects_Tables_ObjectTable_FindId.htm)
        """
        if isinstance(rhObj, rd.RhinoObject):
            return rhObj
        elif isinstance(rhObj, rd.ObjRef):
            return rhObj.Object()
        elif isinstance(rhObj, Guid):
            return sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)


    def getBrepObject(rhObj):
        rdObj = getRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def dotSrf(rgSrf, text='', iDotHeight=14, rgb=(255, 255, 255)):
    
        rgDot = dotGeometryAtSurfaceCentroid(
                rgSrf=rgSrf,
                text=text,
                iDotHeight=iDotHeight)
        if rgDot is None: return
    
        attr = rd.ObjectAttributes()
        attr.ColorSource = rd.ObjectColorSource.ColorFromObject
        attr.ObjectColor = Color.FromArgb(*rgb)
    
        gDot = sc.doc.Objects.AddTextDot(rgDot, attr)
    
        rgDot.Dispose()
    
        if gDot != Guid.Empty: return rgDot


    gBs_In, idxs_F_perB = _getSortedBrepIdsAndFaces(rhFaces)
    if not gBs_In: return

    idxs_AtTenths = [int(round(0.1*i*len(gBs_In),0)) for i in range(10)]

    gBs_1F_Success_allBreps = []

    for iB, gB_In in enumerate(gBs_In):
        if bDebug:
            print('-'*80 + '\n' + '-'*80)
            print("Phase 1")
        
        if sc.escape_test(False):
            print("Searching interrupted by user.")
            return
        
        if iB in idxs_AtTenths:
            Rhino.RhinoApp.CommandPrompt = "Processed {:d}% of {} breps ...".format(
                int(100.0 * (iB+1) / len(gBs_In)), len(gBs_In))
        
        rdBrep_In = getBrepObject(gB_In)
        rgBrep_In = rdBrep_In.Geometry

        idxFaces_TrimSuccess = []
        rgBs_1F_Pos_thisBrep = []

        for iF, idxFace in enumerate(idxs_F_perB[iB]):

            rgFace_In = rgBrep_In.Faces[idxFace]
        
            rgB_1F_Retrimmed = xBrepFace.retrimFace(
                rgFace_In,
                fSplitTol=fSplitTol,
                bTrimToSegs=bTrimToSegs,
                fEdgeLen_Min=fEdgeLen_Min,
                bOutputSplitOnFail=bOutputSplitOnFail,
                bOutputTrimmingCrvs=bOutputTrimmingCrvs,
                bEcho=bEcho,
                bDebug=bDebug,
                )
            if rgB_1F_Retrimmed is None:
                if bDotFailures: dotSrf(rgFace_In, '', (255, 0, 0))
                continue

            rgBs_1F_Pos_thisBrep.append(rgB_1F_Retrimmed)

            idxFaces_TrimSuccess.append(idxFace)

        if not rgBs_1F_Pos_thisBrep:
            continue # to next brep.

        if bEcho:
            if len(idxFaces_TrimSuccess) != len(idxs_F_perB[iB]):
                print("{} faces of {} could not be calculated.".format(
                    len(idxs_F_perB[iB]) - len(idxFaces_TrimSuccess),
                    gB_In))

        if bExtractNotAdd:
            rc = xBrepObject.replaceFaces(
                gB_In,
                idxs_rgFaces=idxFaces_TrimSuccess,
                rgBreps_NewGeom=rgBs_1F_Pos_thisBrep,
                bExtract=True)
            gB_Replaced = rc[0]
            if gB_Replaced:
                gBs_1F_Success_allBreps.extend(gB_Replaced)
            else:
                if bDotFailures:
                    dotSrf(rgFace_In, '', (255, 0, 0))
                if bEcho:
                    print("Attempt to replace {} has failed.".format(gB_In))
        else:
            gB = sc.doc.Objects.AddBrep(rgB_1F_Retrimmed)
            if gB == Guid.Empty:
                bScrapAddFail = True
                if bDotFailures:
                    dotSrf(rgFace_In, '', (255, 0, 0))
                if bEcho:
                    print("Attempt to add a new monoface brep has failed.")
            else:
                gBs_1F_Success_allBreps.append(gB)


        for rgB in rgBs_1F_Pos_thisBrep: rgB.Dispose()

    if bEcho:
        if bExtractNotAdd:
            print("Replaced {} faces with retrimmed monoface breps.".format(
                len(set(gBs_1F_Success_allBreps))))
        else:
            print("Added {} retrimmed monoface breps.".format(
                len(gBs_1F_Success_allBreps)))

    return gBs_1F_Success_allBreps


def main():
    
    rhFaces = getInput()
    if rhFaces is None: return

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    rc = processBrepObjects(rhFaces)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()