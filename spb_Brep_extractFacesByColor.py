"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
210324: Created.
210511: Now breps are processed whose faces are all the same color.
230915: Modified an option default value.
250514: Added command prompts during processing.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List
from System.Drawing import Color

import xBrepObject


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'bSingleColorByPick'; keys.append(key)
    values[key] = True
    names[key] = 'Extract'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'AllColors', 'PickedFaceColor')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTransferFaceColorToObject'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
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

        #if key == 'fPairMatchingDistTol':
        #    if cls.riOpts[key].CurrentValue < 0:
        #        cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
        #    elif cls.riOpts[key].CurrentValue < 1e-9:
        #        cls.riOpts[key].CurrentValue = 1e-9

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select reference breps")

    go.GeometryFilter = rd.ObjectType.Brep

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    idxs_Opt = {}

    colorPicked = None

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bSingleColorByPick')
        addOption('bTransferFaceColorToObject')
        addOption('bEcho')
        addOption('bDebug')

        if Opts.values['bSingleColorByPick']:
            go.DisablePreSelect()
            go.GeometryFilter = rd.ObjectType.Surface
            go.SubObjectSelect = True
            go.SetCommandPrompt("Select face")

            res = go.Get()
            sc.escape_test()

            if res == ri.GetResult.Object:
                objref = go.Object(0)
                rgFace = objref.Face()
                if rgFace:
                    colorPicked = rgFace.PerFaceColor
                    rdBrep = objref.Object()
                    go.Dispose()
                    return (
                        [rdBrep.Id],
                        colorPicked,
                        Opts.values['bTransferFaceColorToObject'],
                        Opts.values['bEcho'],
                        Opts.values['bDebug'],
                        )

                else:
                    print("Face not picked!")
                    print(sc.doc.Objects.UnselectAll())
        else:
            go.GeometryFilter = rd.ObjectType.Brep
            go.GeometryAttributeFilter = (
                ri.Custom.GeometryAttributeFilter.ClosedPolysrf |
                ri.Custom.GeometryAttributeFilter.OpenPolysrf
                )
            go.SubObjectSelect = False
            go.SetCommandPrompt("Select breps")

            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

            if res == ri.GetResult.Object:
                objrefs = go.Objects()
                go.Dispose()

                gBreps = [o.ObjectId for o in objrefs]

                return (
                    gBreps,
                    None,
                    Opts.values['bTransferFaceColorToObject'],
                    Opts.values['bEcho'],
                    Opts.values['bDebug'],
                    )

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return


        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def setObjectColor(gObj, color):
    rdObj = sc.doc.Objects.FindId(gObj)
    rdObj.Attributes.ColorSource = rd.ObjectColorSource.ColorFromObject
    rdObj.Attributes.ObjectColor = color
    return rdObj.CommitChanges()


def removeAllFaceColors(gBrep):
    rdB = sc.doc.Objects.FindId(gBrep)
    rgB_Out = rdB.Geometry.Duplicate()
    for f in rgB_Out.Faces:
        f.PerFaceColor = Color.Empty
    return sc.doc.Objects.Replace(objectId=gBrep, brep=rgB_Out)


def collectFaceColors(rgBrep):
    colors_InBrep = []
    idxFs_perColor = []

    sCP_Base = Rhino.RhinoApp.CommandPrompt

    for iF, f in enumerate(rgBrep.Faces):
        if sc.escape_test(throw_exception=False):
            raise Exception("Search canceled by user.")

        Rhino.RhinoApp.CommandPrompt = "{}  Face {} of {}".format(
            sCP_Base,
            iF + 1,
            rgBrep.Faces.Count)

        if f.PerFaceColor not in colors_InBrep:
            colors_InBrep.append(f.PerFaceColor)
            idxFs_perColor.append([f.FaceIndex])
        else:
            idxFs_perColor[colors_InBrep.index(f.PerFaceColor)].append(f.FaceIndex)

    Rhino.RhinoApp.CommandPrompt = sCP_Base

    return colors_InBrep, idxFs_perColor


def processBrepObjects(gBreps_In, colors_toExtract=None, bXferFaceColorToObj=False, bEcho=True, bDebug=False):
    """
    """

    gBs_Out_AllBs = []

    for iB, gB_In in enumerate(gBreps_In):

        if len(gBreps_In) == 1:
            sCP_Base = "Analyzing 1 brep"
        else:
            sCP_Base = "Analysis at {:d}% of {} breps".format(
                int(100.0 * (iB+1) / len(gBreps_In)), len(gBreps_In))

        if colors_toExtract is not None and len(colors_toExtract) == 1:
            colorFind = colors_toExtract[0]
            sCP_Base += " for ARGB[{},{},{},{}]".format(
                colorFind.A,
                colorFind.R,
                colorFind.G,
                colorFind.G,
                )

        Rhino.RhinoApp.CommandPrompt = sCP_Base

        rdB_In = sc.doc.Objects.FindId(gB_In)
        rgB_In = rdB_In.Geometry

        colors_InBrep, idxFs_perColor = collectFaceColors(rgB_In)


        if len(colors_InBrep) == 1:
            if colors_InBrep[0] == Color.Empty:
                print("No PerFaceColor in brep.")
                continue

            if bDebug:
                print("All faces are the same PerFaceColor in brep.")

            if bXferFaceColorToObj:
                if setObjectColor(gB_In, colors_InBrep[0]):
                    removeAllFaceColors(gB_In)
                    gBs_Out_AllBs.append(gB_In)

            continue


        # More than 1 PerFaceColor in brep.

        Rhino.RhinoApp.CommandPrompt = "{}  {} face colors found".format(
            sCP_Base,
            len(colors_InBrep))


        gBs_Out_ThisB = []

        if colors_toExtract is None:
            # Extract all colors.
            for color, idxFs in zip(colors_InBrep, idxFs_perColor):
                rc = xBrepObject.addFromSubsetOfFaces(
                    gB_In,
                    idxFs,
                    bAddOnlyMonofaces=False,
                    bRetainLayer=True,
                    bRetainColor=True,
                    bDebug=False)
                if rc is not None:
                    gBs_Out_ThisB.extend(rc)
                    if bXferFaceColorToObj and color != Color.Empty:
                        for gB_Out in rc:
                            if setObjectColor(gB_Out, color):
                                removeAllFaceColors(gB_Out)

            if rgB_In.Faces.Count != sum(sc.doc.Objects.FindId(g).Geometry.Faces.Count for g in gBs_Out_ThisB):
                print("Before and after faces counts do not match.  Original brep was not deleted.")
            else:
                if sc.doc.Objects.Delete(objectId=gB_In, quiet=False):
                    print("Original brep was replaced by {} breps of {} colors.".format(
                        len(gBs_Out_ThisB), len(colors_InBrep)))
        else:
            # Extract only colors in colors_InBrep.
            idxs_Fs_NotTargetColor = []
            for i, color in enumerate(colors_InBrep):
                if color in colors_toExtract:
                    idxFs = idxFs_perColor[i]
                    rc = xBrepObject.addFromSubsetOfFaces(
                        rhBrep=gB_In,
                        idxFaces=idxFs,
                        bAddOnlyMonofaces=False,
                        bRetainLayer=True,
                        bRetainColor=True,
                        bDebug=False)
                    if rc is not None:
                        gBs_Out_ThisB.extend(rc)
                        if bXferFaceColorToObj:
                            for gB_Out in rc:
                                if setObjectColor(gB_Out, color):
                                    removeAllFaceColors(gB_Out)
                else:
                    idxs_Fs_NotTargetColor.extend(idxFs_perColor[i])

            rc = xBrepObject.addFromSubsetOfFaces(
                rhBrep=gB_In,
                idxFaces=idxs_Fs_NotTargetColor,
                bAddOnlyMonofaces=False,
                bRetainLayer=True,
                bRetainColor=True,
                bDebug=False)
            if rc is not None:
                gBs_Out_ThisB.extend(rc)

            if rgB_In.Faces.Count != sum(sc.doc.Objects.FindId(g).Geometry.Faces.Count for g in gBs_Out_ThisB):
                print("Before and after faces counts do not match.  Original brep was not deleted.")
            else:
                if sc.doc.Objects.Delete(objectId=gB_In, quiet=False):
                    print("Original brep was replaced by {} breps of {} colors.".format(
                        len(gBs_Out_ThisB), len(colors_InBrep)))


        gBs_Out_AllBs.extend(gBs_Out_ThisB)

    return gBs_Out_AllBs


def main():

    if Rhino.RhinoApp.ExeVersion < 7:
        print("This script requires Rhino Version 7 or higher.")
        return

    rc = getInput()
    if rc is None: return
    (
        gBreps,
        color_Picked,
        bTransferFaceColorToObject,
        bEcho,
        bDebug,
        ) = rc
    
    if not bDebug: sc.doc.Views.RedrawEnabled = False
    
    processBrepObjects(
        gBreps_In=gBreps,
        colors_toExtract=[color_Picked],
        bXferFaceColorToObj=bTransferFaceColorToObject,
        bEcho=bEcho,
        bDebug=bDebug
        )
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()