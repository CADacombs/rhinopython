"""
For a STEP export to include block instances, not their exploded results,
the block instances need to be at full scale and not mirrored.
Note: For block definitions that are represented by a single block instance,
those instances will also be exploded on STEP export.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180709: Created, replacing EvaluateBlockInstanceScale.rvb.
180725: Now supports multiple selected objects.
181117: Refactored.  Added scale routine.
190208: Various bug fixes.
230718: Split from another script.
251027: Refactored many places and removed reliance on rhinoscriptsyntax.
        Added options for modifying block definitions and instances for STEP export.
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


    key = 'bFixScalingForStepExport'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bScaleNonDocUnitBlocks'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fScale'; keys.append(key)
    values[key] = 1.0/25.4
    riOpts[key] = ri.Custom.OptionDouble(values[key], True, Rhino.RhinoMath.SqrtEpsilon)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bInverselyScaleModDefsInsts'; keys.append(key)
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
            if riOpts[key]:
                values[key] = riOpts[key].CurrentValue = sc.sticky[stickyKeys[key]]
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

        if key == 'fScale':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key == 'fTol_IsEllipse':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < cls.riOpts[key].InitialValue:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def _staticBlockDefs():
    """
    Static == Embedded only
    """

    ret = []
    for rdDef in sc.doc.InstanceDefinitions.GetList(ignoreDeleted=True):
        if rdDef.UpdateType == rd.InstanceDefinitionUpdateType.Static:
            ret.append(rdDef)
    return ret


def getInput():
    """
    Get options.
    """

    go = ri.Custom.GetOption()
    go.SetCommandPrompt("Set options")

    go.AcceptNothing(True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        addOption('bFixScalingForStepExport')
        if Opts.values['bFixScalingForStepExport']:
            go.AcceptNumber(False, acceptZero=False)
        else:
            addOption('bScaleNonDocUnitBlocks')
            if Opts.values['bScaleNonDocUnitBlocks']:
                go.AcceptNumber(False, acceptZero=False)
            else:
                addOption('fScale')
                go.AcceptNumber(True, acceptZero=False)
            addOption('bInverselyScaleModDefsInsts')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return True

        if res == ri.GetResult.Number:
            key = 'fScale'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def scaleContentsOfBlockDefinition(rdDef, fScale, bEcho=True, bDebug=False):
    """
    """

    if fScale <= 0.0:
        print("fScale must be > 0.0. {} was provided".format(fScale))
        return

    rdObjs_InBlock = rdDef.GetObjects()

    if not rdObjs_InBlock:
        if bEcho:
            print("Block definition, {}, contains no objects, so it will not be processed.".format(
                rdDef.Name))
        return

    scale = rg.Transform.Scale(
        anchor=rg.Point3d.Origin,
        scaleFactor=fScale)

    geoms = []
    attrs = []

    for rdObj in rdObjs_InBlock:
        geom = rdObj.Geometry
        attr = rdObj.Attributes
        if isinstance(rdObj, rd.InstanceObject):
            pt = rg.Point3d.Origin
            pt.Transform(rdObj.InstanceXform)
            translation = rg.Transform.Translation(fScale*pt - pt)
            geom.Transform(translation)
        else:
            geom.Transform(scale)

        geoms.append(rdObj.Geometry)
        attrs.append(rdObj.Attributes)

    return sc.doc.InstanceDefinitions.ModifyGeometry(rdDef.Index, geoms, attrs)


def _getInstanceScaleComponents(rdInst):
    xform = rdInst.InstanceXform
    return tuple(
        [(xform[0,c]**2.0 + xform[1,c]**2.0 + xform[2,c]**2.0)**0.5 for c in (0,1,2)]
        )


def scaleRootLevelInstances(rdDef, fScale, bEcho=True, bDebug=False):
    """
    """

    if fScale <= 0.0:
        print("fScale must be > 0.0. {} was provided".format(fScale))
        return

    rdInsts = rdDef.GetReferences(wheretoLook=0)
    if not rdInsts:
        print("No instances at root level for rdDef.Name.")
        return

    #sEval = "fScale"; print(sEval, '=', eval(sEval))

    gInsts_Scaled = []
    gFails = []

    for rdInst in rdInsts:
        pt = rg.Point3d.Origin
        pt.Transform(rdInst.InstanceXform)
        #sEval = "pt"; print(sEval, '=', eval(sEval))

        xform_Scale = rg.Transform.Scale(
            anchor=pt,
            scaleFactor=fScale)

        #id = scriptcontext.doc.Objects.Transform(old_id, xform, not copy)
        rv = sc.doc.Objects.Transform(
            obj=rdInst,
            xform=xform_Scale,
            deleteOriginal=True)

        #if not rdInst.Geometry.Transform(xform_Scale):
        #    print("Failed to scale the instance geometry of {}.".format(rdDef.Name))
        #    return

        #rv = rdInst.CommitChanges() # CommitChanges doesn't work.

        if bDebug: print(rv)
        if rv != rdInst.Id:
            print("Scaling {} failed.".format(rdInst.Id))
            gFails.append(rdInst.Id)
        else:
            gInsts_Scaled.append(rdInst.Id)
            rdInst_Out = sc.doc.Objects.FindId(rv)
            svs = _getInstanceScaleComponents(rdInst_Out)
            for sv in svs:
                if not Rhino.RhinoMath.EpsilonEquals(sv, 1.0, epsilon=Rhino.RhinoMath.Epsilon):
                    print("Warning, instance {} of block {} is not at full scale, but instead {}, {}, {}".format(
                        rdInst_Out.Id, rdDef.Name, *svs))
                    break

    return gInsts_Scaled

    #print("Scaled {} root-level instance(s) by {:.20g} about their insertion points.".format(
    #        n_Insts_processed, 1.0/fScale))


def main():

    rdDefs_Static = _staticBlockDefs()

    if not rdDefs_Static:
        print("Document has no static block definitions.")
        return

    if sc.doc.Modified:
        showMessageResult = Rhino.UI.Dialogs.ShowMessage(
            message="This document has been modified since the last save." \
                "\n\nIn the case that this script produces erroneous results," \
                " it is recommended to first press the Cancel button" \
                " then _Save or _SaveACopy before proceeding." \
                "\n\nPress OK to continue fixing the scaling.",
            title="Document Not Saved",
            buttons=Rhino.UI.ShowMessageButton.OKCancel,
            icon=Rhino.UI.ShowMessageIcon.Warning)
        if showMessageResult == Rhino.UI.ShowMessageResult.Cancel:
            return

    rv = getInput()
    if rv is None: return

    bFixScalingForStepExport = Opts.values['bFixScalingForStepExport']
    if bFixScalingForStepExport:
        bScaleNonDocUnitBlocks = True
        fScale = None
        bInverselyScaleModDefsInsts = True
    else:
        bScaleNonDocUnitBlocks = Opts.values['bScaleNonDocUnitBlocks']
        fScale = Opts.values['fScale']
        bInverselyScaleModDefsInsts = Opts.values['bInverselyScaleModDefsInsts']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = True

    if bScaleNonDocUnitBlocks:
        rdDefs_toProcess = [rdDef for rdDef in rdDefs_Static if rdDef.UnitSystem != sc.doc.ModelUnitSystem]
        if not rdDefs_toProcess:
            print("No block definitions exist whose units are not {}, the document units.".format(
                sc.doc.ModelUnitSystem))
            return
    else:
        rdDefs_toProcess = rdDefs_Static


    #res, bProceed = ri.RhinoGet.GetBool(
    #    prompt="{} block definitions to have their contents scaled. Proceed?".format(
    #            len(rdDefs_toProcess)),
    #    acceptNothing=True,
    #    offPrompt="No", onPrompt="Yes", boolValue=False)
    #if res != Rhino.Commands.Result.Success or not bProceed:
    #    print("Nothing was modified.")
    #    return


    rdDefs_Scaled = []
    rdInsts_Scaled = []

    if bScaleNonDocUnitBlocks:
        if bDebug: print("Doc unit: {}".format(sc.doc.ModelUnitSystem))

        for rdDef in rdDefs_toProcess:
            if bDebug: print("{} unit: {}".format(rdDef.Name, rdDef.UnitSystem))
            if rdDef.UnitSystem == sc.doc.ModelUnitSystem:
                if bDebug: print("Units match, so this block definition will not be modified.")
                continue
            fScale_ThisBlock = Rhino.RhinoMath.UnitScale(
                rdDef.UnitSystem,
                sc.doc.ModelUnitSystem)
            if bDebug: print("Will scale block objects by {}.".format(fScale_ThisBlock))

            rv = scaleContentsOfBlockDefinition(
                rdDef,
                fScale=fScale_ThisBlock,
                bEcho=bEcho,
                bDebug=bDebug)
            if rv:
                rdDefs_Scaled.append(rdDef)
            else:
                print("Failed to scale contents of {}.".format(rdDef.Name))
                continue

            rdDef.UnitSystem = sc.doc.ModelUnitSystem

            if bInverselyScaleModDefsInsts:
                rv = scaleRootLevelInstances(
                    rdDef,
                    fScale=1.0/fScale_ThisBlock,
                    bEcho=bEcho,
                    bDebug=bDebug)
                if bDebug: print(rv)
                if rv is None: continue
                rdInsts_Scaled.extend(rv)
    else:
        for rdDef in rdDefs_toProcess:

            rv = scaleContentsOfBlockDefinition(
                rdDef,
                fScale=fScale,
                bEcho=bEcho,
                bDebug=bDebug)
            if rv:
                rdDefs_Scaled.append(rdDef)
            else:
                print("Failed to scale contents of {}.".format(rdDef.Name))
                continue

            if bInverselyScaleModDefsInsts:
                rv = scaleRootLevelInstances(
                    rdDef,
                    fScale=1.0/fScale,
                    bEcho=bEcho,
                    bDebug=bDebug)
                if bDebug: print(rv)
                if rv is None: continue
                rdInsts_Scaled.extend(rv)


    print("Scaled objects in {} definitions.".format(len(rdDefs_Scaled)))
    print("Scaled {} instances.".format(len(rdInsts_Scaled)))


if __name__ == '__main__': main()