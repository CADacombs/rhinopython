"""
This script is an alternative to _ToNURBS, which has problems with some RevSurfaces.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180206-07: Created, starting with fitSurface.py
...
190519: Updated an import name.
200207: Now accepts BrepFaces.  Refactored.  Removed reporting routine.
200414: Added bAcceptParamMismatch and now uses an external module to trim the face.
200422: Now, doesn't attempt to trim closed surfaces from closed monoface breps.
200428: Bug fix.
200527, 0619: Import-related update.
231108: Improved some debug reporting.
250326: Bug fix for processing individual faces of a polysrf. Refactored.
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


    key = 'fTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'Tolerance'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAcceptParamMismatch'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='Add', onValue='Replace')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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


def getInput():
    """
    Get surfaces and/or breps and option values.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps and/or faces")

    go.GeometryFilter = rd.ObjectType.Brep

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opts = {}


    while True:
        key = 'fTol'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bAcceptParamMismatch'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bReplace'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])

        # An option was selected or a number was entered.
        
        key = 'fTol'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def processBrepObjects(rhBreps, idxFaces_perBrep, **kwargs):
    """
    rhBreps: Can be rd.BrepObjects or GUIDs,
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTol = getOpt('fTol')
    bAcceptParamMismatch = getOpt('bAcceptParamMismatch')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


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


    def areAllBrepTrimIsostatusesSENW(rgFace0):
        if rgFace0.IsSurface: return True
        iCt_Trims = 0
        for rgTrim_ in rgFace0.Brep.Trims:
            if rgTrim_.Face.FaceIndex == rgFace0.FaceIndex:
                if rgTrim_.IsoStatus == rg.IsoStatus.None:
                    return False
                elif rgTrim_.IsoStatus == rg.IsoStatus.X:
                    return False
                elif rgTrim_.IsoStatus == rg.IsoStatus.Y:
                    return False
                iCt_Trims += 1
        if iCt_Trims != 4:
            return False
        return True


    iFs_perB = idxFaces_perBrep


    bNon1AccuracyProduced = False

    gBs_Out = []
    iCt_Fs_Success = 0
    sLogs = []

    for iB in range(len(rhBreps)):
        rhB = rhBreps[iB]
        rc = getBrep(rhB)
        if not rc: continue
        rdB_In, rgB_In = rc
        if iFs_perB is None or iFs_perB[iB] is None:
            iFs = range(rgB_In.Faces.Count)
        else:
            iFs = iFs_perB[iB]

        iFs_Processed = []
        rgBs_1Fs_NewGeom = []

        for iF in iFs:
            rgSrf_In = rgB_In.Faces[iF].UnderlyingSurface()
            if isinstance(rgSrf_In, rg.NurbsSurface):
                continue

            rgNSrf, iAccuracy = rgSrf_In.ToNurbsSurface(tolerance=fTol)
            if rgNSrf is None:
                continue

            # accuracy result per https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Surface_ToNurbsSurface_1.htm
            if bDebug:
                print("accuracy set (out parameter) by ToNurbsSurface: {}".format(iAccuracy))
                if iAccuracy == 0:
                    print("0 = unable to create NURBS representation with desired accuracy.")
                elif iAccuracy == 1:
                    print("1 = success - returned NURBS parameterization matches the surface's to the desired accuracy.")
                elif iAccuracy == 2:
                    print("2 = success - returned NURBS point locus matches the surface's to the desired accuracy and the domain of the NURBS surface is correct. However, this surface's parameterization and the NURBS surface parameterization may not match to the desired accuracy. This situation happens when getting NURBS representations of surfaces that have a transendental parameterization like spheres, cylinders, and cones.")
                else:
                    raise ValueError("What happened?")

            if iAccuracy == 1:
                if rgB_In.Faces.Count == 1 and rgB_In.IsSolid:
                    rgB_1F_Out = rgNSrf.ToBrep()
                    rgBs_1Fs_NewGeom.append(rgB_1F_Out)
                    iFs_Processed.append(iF)
                elif rgB_In.Faces.Count == 1 and rgB_In.Faces[iF].IsSurface:
                    rgB_1F_Out = rgNSrf.ToBrep()
                    rgBs_1Fs_NewGeom.append(rgB_1F_Out)
                    iFs_Processed.append(iF)
                elif rgB_In.Faces.Count == 1 and areAllBrepTrimIsostatusesSENW(rgB_In.Faces[iF]):
                    rgB_1F_Out = rgNSrf.ToBrep()
                    rgBs_1Fs_NewGeom.append(rgB_1F_Out)
                    iFs_Processed.append(iF)
                else:
                    rgB_1F_WIP = rgB_In.Faces[iF].DuplicateFace(duplicateMeshes=False)
                    
                    rgB_1F_WIP.AddSurface(rgNSrf)
    
                    rgNSrf.Dispose()
                    if not rgB_1F_WIP.Faces[0].ChangeSurface(1):
                        print("Brep.ChangeSurface failed.")
                    else:
                        rgB_1F_WIP.Compact()
    
                        rgBs_1Fs_NewGeom.append(rgB_1F_WIP)
                        iFs_Processed.append(iF)
            elif iAccuracy == 2 and bAcceptParamMismatch:
                
                if rgB_In.Faces.Count == 1 and rgB_In.IsSolid:
                    rgB_1F_Out = rgNSrf.ToBrep()
                    rgBs_1Fs_NewGeom.append(rgB_1F_Out)
                    iFs_Processed.append(iF)
                elif rgB_In.Faces.Count == 1 and rgB_In.Faces[iF].IsSurface:
                    rgB_1F_Out = rgNSrf.ToBrep()
                    rgBs_1Fs_NewGeom.append(rgB_1F_Out)
                    iFs_Processed.append(iF)
                elif rgB_In.Faces.Count == 1 and areAllBrepTrimIsostatusesSENW(rgB_In.Faces[iF]):
                    rgB_1F_Out = rgNSrf.ToBrep()
                    rgBs_1Fs_NewGeom.append(rgB_1F_Out)
                    iFs_Processed.append(iF)
                else:
                    rgB_1F_WIP_Retrimmed = xBrepFace.retrimFace(
                        rgB_In.Faces[iF],
                        rgSrf_Replacement=rgNSrf,
                        bDebug=bDebug
                        )

                    rgNSrf.Dispose()

                    if not rgB_1F_WIP_Retrimmed:
                        rgB_In.Dispose()
                        continue

                    rgBs_1Fs_NewGeom.append(rgB_1F_WIP_Retrimmed)
                    iFs_Processed.append(iF)
            else:
                if iAccuracy == 0:
                    sLogs.append(
                        "unable to create NURBS representation with desired accuracy.")
                    #sc.doc.Objects.AddSurface(rgSrf_In)
                elif iAccuracy == 2:
                    sLogs.append(
                        "this surface's parameterization and the NURBS surface parameterization may not match to the desired accuracy." \
                            "  This situation happens when getting NURBS representations of surfaces that have a transendental parameterization like spheres, cylinders, and cones.")
                    #sc.doc.Objects.AddSurface(rgSrf_In)

                if not bNon1AccuracyProduced:
                    s  = "For Accuracy level, see "
                    s += "http://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Surface_ToNurbsSurface_1.htm"
                    print(s)
                    bNon1AccuracyProduced = True


        if not iFs_Processed:
            rgB_In.Dispose()
            continue


        if bReplace:
            rc = xBrepObject.replaceFaces(
                rdB_In, iFs_Processed, rgBs_1Fs_NewGeom,
                bExtract=False,
                )
            if not rc: continue
            gBs_Out.extend(rc)
            iCt_Fs_Success += len(iFs_Processed)
        else:
            for iB in range(len(rgBs_1Fs_NewGeom)):
                gB_Out = sc.doc.Objects.AddBrep(rgBs_1Fs_NewGeom[iB])
                if gB_Out == Guid.Empty: continue
                gBs_Out.append(rc)
                iCt_Fs_Success += 1


    if sLogs:
        for sLog in set(sLogs):
            print("[{}] {}".format(sLogs.count(sLog), sLog))


    print("{} faces {}.".format(
        iCt_Fs_Success,
        "replaced" if bReplace else "added"))

    return gBs_Out


def processObjRefs(objrefs, **kwargs):
    """
    objrefs
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTol = getOpt('fTol')
    bAcceptParamMismatch = getOpt('bAcceptParamMismatch')
    bReplace = getOpt('bReplace')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if not all([isinstance(objref, rd.ObjRef) for objref in objrefs]):
        print("Only input ObjRefs. Script canceled.")
        return


    def sortBrepObjects_and_face_indices(objrefs):
        gBs = []; iFs_perB = []
    
        for objref in objrefs:
            if objref.Object().ObjectType != rd.ObjectType.Brep: continue

            if objref.GeometryComponentIndex.Index == -1:
                if objref.ObjectId in gBs:
                    iB = gBs.index(objref.ObjectId)
                    iFs_perB[iB] = None
                else:
                    gBs.append(objref.ObjectId)
                    iFs_perB.append(None)
            elif (
                objref.GeometryComponentIndex.ComponentIndexType ==
                rg.ComponentIndexType.BrepFace
            ):
                if objref.ObjectId in gBs:
                    iB = gBs.index(objref.ObjectId)
                    if iFs_perB[iB] is None:
                        pass
                    else:
                        iFs_perB[iB].append(objref.GeometryComponentIndex.Index)
                else:
                    gBs.append(objref.ObjectId)
                    iFs_perB.append([objref.GeometryComponentIndex.Index])


        iFs_perB_Sorted = []
        for iFs in iFs_perB:
            if iFs is None:
                iFs_perB_Sorted.append(None)
            else:
                iFs_perB_Sorted.append(sorted(iFs))

        return gBs, iFs_perB_Sorted


    rdBreps, idxFaces_perBrep = sortBrepObjects_and_face_indices(objrefs)

    return processBrepObjects(
        rhBreps=rdBreps, idxFaces_perBrep=idxFaces_perBrep)


def main(bDebug=False):
    
    rc = getInput()
    if rc is None: return

    objrefs = rc[0]
    bDebug = Opts.values['bDebug']

    if bDebug:
        pass
    else:
        sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    rc = processObjRefs(
        objrefs=objrefs)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main(bDebug=bool(0))