"""
To reduce execution time when testing for both cases,
No-area faces and slow-area faces are tested in same loop in this script.

170712-13: Created.
170715: Now uses Brep.GetArea() and also checks for unsuccessful computes.
...
180604: Changed default for bRetainAttr from True to False.
...
190405: Replaced print with prompt.
190513: Added sc.escape_test in brep loop.
190605: Import-related update.
190625: To speed up execution, feedback now reports at 10% intervals.
190711: Now will select all normal breps with Enter and no objects are selected.  Added bExtract.
190716: Made printed output more readable.
190719: Bug fixes.
191023-24: Refactored.
191101,07: Import-related update.
191121-22: Modified some options and refactored dotting routine.
200701: Import-related update.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List
from System.Drawing import Color

import time

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


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'fTime_MaxTol'; keys.append(key)
    values[key] = 1.0
    names[key] = 'MinSecondsForSlowArea'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRetainAttr'; keys.append(key)
    values[key] = True
    names[key] = 'RetainAttributes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtractNoArea'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtractSlowArea'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHeight'; keys.append(key)
    values[key] = 11
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    riAddOpts[key] = addOptionInteger(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]


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


def getInput():
    """
    Get breps with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none selected")
    
    go.GeometryFilter = rd.ObjectType.Brep
    
    go.AcceptNothing(True)
    go.AcceptNumber(True, True)
    
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.
    
    bPreselectedObjsChecked = False

    print "For MaxSecondsAllowed, enter 0 or less to skip extracting or selecting faces over a time limit."
    
    while True:
        Opts.riAddOpts['fTime_MaxTol'](go)
        Opts.riAddOpts['bExtract'](go)
        if Opts.values['bExtract']:
            Opts.riAddOpts['bRetainAttr'](go)
        Opts.riAddOpts['bExtractNoArea'](go)
        Opts.riAddOpts['bAddDot'](go)
        if Opts.values['bAddDot']:
            Opts.riAddOpts['iDotHeight'](go)
        Opts.riAddOpts['bEcho'](go)
        Opts.riAddOpts['bDebug'](go)

        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])
        if res == ri.GetResult.Nothing:
            iter = rd.ObjectEnumeratorSettings()
            iter.NormalObjects = True
            iter.LockedObjects = False
            iter.IncludeLights = False
            iter.IncludeGrips = False
            rdRhinoObjects = []
            for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
                if rdRhinoObject.ObjectType == rd.ObjectType.Brep:
                    rdRhinoObjects.append(rdRhinoObject)
            go.Dispose()
            return tuple([rdRhinoObjects] + [Opts.values[key] for key in Opts.keys])
        else:
            # An option was selected or a number was entered.
            key = 'fTime_MaxTol'
            if res == ri.GetResult.Number:
                Opts.riOpts[key].CurrentValue = go.Number()

            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()


def getFacesAreaTimes(rgBrep, fTime_MaxTol=Opts.riOpts['fTime_MaxTol'].InitialValue):
    """
    Parameters:
        rgBrep: rg.Brep
        fTime_MaxTol: float for seconds
        bEcho: bool
    Returns tuple:
        fTimes (float or None),
        ptCentroids) of all faces
    """

    # Check whether brep is invalid.
    #if not rgBrep.IsValid:
    #    print "Brep is invalid and will be skipped.  " \
    #            "_ExtractBadSrf and repair the bad faces first."
    #    return
    
    fTimes = []
    
    iCt_Fs = rgBrep.Faces.Count
    
    idxs_AtTenths = [int(round(0.1*i*iCt_Fs,0)) for i in range(10)]

    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt
    
    for idxF in range(rgBrep.Faces.Count):

        if iCt_Fs > 10:
            if idxF in idxs_AtTenths:
                Rhino.RhinoApp.SetCommandPrompt(
                    sCmdPrompt0 +
                    "  Analysis at {:d}% of {} faces ...".format(
                        int(100.0 * (idxF+1) / iCt_Fs), iCt_Fs))
        elif iCt_Fs > 1:
            Rhino.RhinoApp.SetCommandPrompt(
                sCmdPrompt0 +
                "  Analyzing {} of {} faces".format(
                    idxF+1, iCt_Fs))

        rgBrep_1Face = rgBrep.Faces[idxF].DuplicateFace(True)

        time0 = time.time()
        area = rgBrep_1Face.GetArea() # Value is 0.0, not None, for no area compute.
        timeElapsed = time.time() - time0

        if not area:
            # Value is 0.0, not None, for no area compute.
            fTimes.append(None)
        else:
            fTimes.append(timeElapsed)

        rgBrep_1Face.Dispose()

    Rhino.RhinoApp.SetCommandPrompt(sCmdPrompt0)

    return fTimes


def processBrepObjects(rhBreps0, **kwargs):
    """
    Parameters:
        rhBreps0: rd.Breps or ObjRefs of Brep.
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTime_MaxTol = getOpt('fTime_MaxTol')
    bExtract = getOpt('bExtract')
    bRetainAttr = getOpt('bRetainAttr')
    bExtractNoArea = getOpt('bExtractNoArea')
    bAddDot = getOpt('bAddDot')
    iDotHeight = getOpt('iDotHeight')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if bDebug: bEcho = True


    def coerceRhinoObject(rhObj):
        rdObj = None
        if isinstance(rhObj, rd.RhinoObject):
            rdObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        return rdObj


    def coerceBrepObject(rhObj):
        rdObj = coerceRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    rgDots = []
    rdBreps0 = []
    gBreps_Extracted_From1Brep = [] # Accumulation of test positive, extracted single-face breps.
    iCt_NoAreaFaces_All = 0
    iCt_SlowAreaFaces_All = 0
    time0_Loop = time.time()
    time_MaxBrep = (0.0, None, None)
    iCt_Selected = 0

    len_rhBs0 = len(rhBreps0)
    
    idxs_AtTenths = [int(round(0.1*i*len(rhBreps0),0)) for i in range(10)]
    
    for iB, rhBrep in enumerate(rhBreps0):
        if sc.escape_test(False):
            print "Searching interrupted by user."
            return
        
        if len_rhBs0 > 10:
            if iB in idxs_AtTenths:
                Rhino.RhinoApp.SetCommandPrompt("Analysis at {:d}% of {} breps ...".format(
                    int(100.0 * (iB+1) / len_rhBs0), len_rhBs0))
        elif len_rhBs0 > 1:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing {} of {} breps".format(
                iB+1, len_rhBs0))
        else:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing brep ...")
        
        rdBrep0 = coerceBrepObject(rhBrep)
        if not rdBrep0:
            if bEcho and len_rhBs0 == 1:
                print "Non-brep DocObject passed to processBrepObjects."
            continue
        
        rgBrep0 = rdBrep0.BrepGeometry

        if not rgBrep0.IsValid and bExtract:
            if bEcho and len_rhBs0 == 1:
                s  = "Brep {} is invalid and will be skipped.".format(rdBrep0.Id)
                s += "  _ExtractBadSrf and repair the bad faces first."
                print s
            rgBrep0.Dispose()
            continue
        
        rdBreps0.append(rdBrep0)

        iCt_Faces_B0 = rgBrep0.Faces.Count
        
        # Get indices and times of slow Brep.GetArea brep faces.
        time0_Brep = time.time()


        rc = getFacesAreaTimes(
                rgBrep0,
                fTime_MaxTol,
        )
        if rc is None:
            rgBrep0.Dispose()
            continue

        idx_rgFaces_NoArea = []
        idx_rgFaces_SlowArea = []
        fTimes = []
        for idxF, fTime in enumerate(rc):
            if fTime is None:
                idx_rgFaces_NoArea.append(idxF)
            elif  0.0 < fTime_MaxTol <= fTime:
                idx_rgFaces_SlowArea.append(idxF)
                fTimes.append(fTime)

        time_Entire_B0 = time.time() - time0_Brep
        if time_Entire_B0 > time_MaxBrep[0]:
            time_MaxBrep = (time_Entire_B0, rdBrep0.Id, iCt_Faces_B0)
        
        if not (idx_rgFaces_NoArea + idx_rgFaces_SlowArea): # No matches in this brep
            rgBrep0.Dispose()
            continue


        if bAddDot and (idx_rgFaces_NoArea + idx_rgFaces_SlowArea):
            rgBrep_ShrunkFs = rgBrep0.DuplicateBrep()
            rgBrep_ShrunkFs.Faces.ShrinkFaces()
        else:
            rgBrep_ShrunkFs = None

        for idxF in idx_rgFaces_NoArea:
            rgFace = rgBrep0.Faces[idxF]
                
            if bAddDot:
                rgSrf = rgBrep_ShrunkFs.Faces[idxF].UnderlyingSurface()
                ptMidDomain = rgSrf.PointAt(rgSrf.Domain(0).Mid, rgSrf.Domain(1).Mid)

                if ptMidDomain is None:
                    s  = "Dot could not be added for Face"
                    if rgBrep0.Faces.Count == 1:
                        s += " of {}.".format(rdBrep0.Id)
                    else:
                        s += "[{}] of {}.".format(idxF, rdBrep0.Id)
                    print s
                else:
                    rgDot = rg.TextDot(
                        text="NoArea", location=ptMidDomain)
                    rgDot.FontHeight = iDotHeight
                    rgDots.append(rgDot)

            if bEcho and len_rhBs0 == 1:
                print "No area compute for " \
                        "Face {} of {} of {}.".format(
                        idxF, rgBrep0.Faces.Count, rdBrep0.Id)

        for idxF, fTime in zip(idx_rgFaces_SlowArea, fTimes):
            rgFace = rgBrep0.Faces[idxF]
                
            if bEcho and len_rhBs0 == 1:
                if 0 < len(idx_rgFaces_SlowArea) <= 10:
                    s  = "    Face[{}]: ".format(idxF)
                    s += "{:.3f} seconds".format(fTime)
                    print s

            if bAddDot:
                rgSrf = rgBrep_ShrunkFs.Faces[idxF].UnderlyingSurface()
                ptMidDomain = rgSrf.PointAt(rgSrf.Domain(0).Mid, rgSrf.Domain(1).Mid)

                if ptMidDomain is None:
                    if bEcho:
                        s  = "Dot could not be added for Face"
                        if rgBrep0.Faces.Count == 1:
                            s += " of {}.".format(rdBrep0.Id)
                        else:
                            s += "[{}] of {}.".format(idxF, rdBrep0.Id)
                        print s
                else:
                    rgDot = rg.TextDot(
                        text="{:.3f}s".format(fTime), location=ptMidDomain)
                    rgDot.FontHeight = iDotHeight
                    rgDots.append(rgDot)

        if rgBrep_ShrunkFs: rgBrep_ShrunkFs.Dispose()

        iCt_NoAreaFaces_All += len(idx_rgFaces_NoArea)
        iCt_SlowAreaFaces_All += len(idx_rgFaces_SlowArea)

        rgBrep0.Dispose()

        idx_Fs_toExtractOrSelect = (
            (idx_rgFaces_NoArea if bExtractNoArea else []) +
            (idx_rgFaces_SlowArea if fTime_MaxTol > 0.0 else []))

        if not idx_Fs_toExtractOrSelect: continue

        if not bExtract and iCt_Selected == 0: sc.doc.Objects.UnselectAll()

        if bExtract:
            rc = xBrepObject.extractFaces(
                    rdBrep0,
                    idx_Fs_toExtractOrSelect,
                    bAddOnlyMonofaces=True,
                    bRetainLayer=bRetainAttr,
                    bRetainColor=bRetainAttr,
                    bEcho=bEcho,
                    bDebug=bDebug)
            if rc is None:
                if bEcho and len_rhBs0 == 1: print "Extract faces failed."
                continue
            gBreps_Extracted_From1Brep.extend(rc[0])
        else:
            if iCt_Faces_B0 == 1:
                sc.doc.Objects.Select(objectId=rdBrep0.Id)
                iCt_Selected += 1
            else:
                if Rhino.RhinoApp.ExeVersion >= 6:
                    for idx_rgFace in idx_Fs_toExtractOrSelect:
                        compIdx = rg.ComponentIndex(
                                rg.ComponentIndexType.BrepFace,
                                index=idx_rgFace)
                        rdBrep0.SelectSubObject(
                                componentIndex=compIdx,
                                select=True, syncHighlight=True, persistentSelect=True)
                        iCt_Selected += 1
        
        # End of rhBreps0 loop.
    
    if time_MaxBrep[1] is None: return
    
    if bDebug:
        print "Maximum time: {} seconds for findAllBadAreaFaces for " \
                "{} with {} faces.".format(*time_MaxBrep)
    
    if bEcho:
        s  = "{:.3f} seconds to search for and report on".format(
                time.time()-time0_Loop)
        s += " {} breps.".format(len(rdBreps0))
    
    if iCt_NoAreaFaces_All + iCt_SlowAreaFaces_All == 0:
        if bEcho:
            s += "\n  No faces found whose Brep.GetArea cannot be computed"
            if fTime_MaxTol > 0.0:
                s += " or took longer than {:.1f} seconds.".format(fTime_MaxTol)
            else:
                s += "."
            print s
        return
    
    if bEcho:
        s += "\n  {} faces found whose area cannot be computed.".format(
                iCt_NoAreaFaces_All)
    
        if fTime_MaxTol > 0.0:
            s += "\n  {} faces found whose area calculation took longer " \
                    "than {:.3f} seconds.".format(
                        iCt_SlowAreaFaces_All,
                        fTime_MaxTol)


    if gBreps_Extracted_From1Brep:
        sc.doc.Objects.UnselectAll()

        iCt_Selected = sc.doc.Objects.Select(List[Guid](gBreps_Extracted_From1Brep))
        if iCt_Selected == 0:
            if bEcho:
                s += "Error none of the faces are selected!"
                print s
            return
    
        s += "\n  {} monoface breps are selected.".format(iCt_Selected)

    if not bExtract and iCt_Selected:
        s += "\n  {} monoface breps / faces are selected.".format(iCt_Selected)

    if bEcho: print s

    if rgDots:
        attrDot = rd.ObjectAttributes()
        attrDot.ColorSource = rd.ObjectColorSource.ColorFromObject
        attrDot.ObjectColor = Color.FromArgb(255,0,0)
        gDots = []
        for rgDot in rgDots:
            gDots.append(sc.doc.Objects.AddTextDot(rgDot, attrDot))
            rgDot.Dispose()

        if gDots: sc.doc.Objects.Select(List[Guid](gDots))

    sc.doc.Views.Redraw()


def main():
    
    rc = getInput()
    if rc is None: return
    rhBreps0 = rc[0]
    
    processBrepObjects(rhBreps0)


if __name__ == '__main__': main()