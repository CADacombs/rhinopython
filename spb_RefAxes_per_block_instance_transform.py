"""
This script adds 3 lines, and optionally, 3 planes, and 4 points, representing
the X, Y, and Z axes per the transformation of the selected target block instance.
"""

#!  python 2 Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
241007-08: Created, starting with another script.
250219: Replaced most functions with those of an import.
250910: Import-related bug fixes.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import spb_RefAxes


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAddToDef'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddAxesBlock'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fLineLength'; keys.append(key)
    values[key] = Rhino.RhinoMath.UnitScale(Rhino.UnitSystem.Inches, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPlanes'; keys.append(key)
    values[key] = True
    names[key] = 'IncludePlanes'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bOriginPt'; keys.append(key)
    values[key] = True
    names[key] = 'IncludeOriginPt'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAxesPts'; keys.append(key)
    values[key] = False
    names[key] = 'IncludeAxesPts'
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
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def addOption(cls, go, key):

        idxOpt = None

        if key in cls.riOpts:
            if key[0] == 'b':
                idxOpt = go.AddOptionToggle(cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                idxOpt = go.AddOptionDouble(cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                idxOpt = go.AddOptionInteger(englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        elif key in cls.listValues:
            idxOpt = go.AddOptionList(englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])
        else:
            print("{} is not a valid key in Opts.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick block instances")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.InstanceReference

    def addOption(key):
        Opts.addOption(go, key)
        keys_in_order_added.append(key)

    while True:
        go.ClearCommandOptions()
        keys_in_order_added = [None] # None is just a placeholder since option indices are base 1.

        addOption('bAddToDef')
        addOption('bAddAxesBlock')
        addOption('fLineLength')
        addOption('bPlanes')
        addOption('bOriginPt')
        addOption('bAxesPts')
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

        Opts.setValue(keys_in_order_added[go.Option().Index])


def main():

    objrefs = getInput()
    if objrefs is None: return

    bAddToDef = Opts.values['bAddToDef']
    bAddAxesBlock = Opts.values['bAddAxesBlock']
    fLineLength = Opts.values['fLineLength']
    bPlanes = Opts.values['bPlanes']
    bOriginPt = Opts.values['bOriginPt']
    bAxesPts = Opts.values['bAxesPts']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if bAddAxesBlock:
        bDefineBlock = True
        sAxesBlockName = "RefAxes"
        if sc.doc.InstanceDefinitions.Find(sAxesBlockName) is not None:
            iExistingBlockChoice = rs.MessageBox(
                "Press 'Yes' to redefine it or 'No' to use existing block",
                3,
                sAxesBlockName + " already exists.")
            if iExistingBlockChoice == 2:
                return
            elif iExistingBlockChoice == 7:
                bDefineBlock = False

    sc.doc.Views.RedrawEnabled = False

    if bAddAxesBlock and bDefineBlock:
        strLayer = "Default"
        if not rs.IsLayer(strLayer):
            rs.AddLayer(strLayer)
        gGeomsForAxesBlock = spb_RefAxes.addLineCurves(fLineLength)
        if bPlanes:
            gGeomsForAxesBlock.extend(spb_RefAxes.addSurfaceObjects_Planes(fLineLength))
        if bOriginPt:
            gGeomsForAxesBlock.append(spb_RefAxes.addPointObject_Origin())
        if bAxesPts:
            gGeomsForAxesBlock.extend(spb_RefAxes.addPointObjects_Axes(fLineLength))
        rs.ObjectLayer(gGeomsForAxesBlock , strLayer)
        rs.AddBlock(gGeomsForAxesBlock, rg.Point3d.Origin, name=sAxesBlockName, delete_input=True)


    if bAddToDef:
        gDefs = list(set(_.Geometry().ParentIdefId for _ in objrefs))
        for gDef in gDefs:
            idef = sc.doc.InstanceDefinitions.FindId(gDef)
            geoms = []
            attrs = []
            for rdObj in idef.GetObjects():
                geoms.append(rdObj.Geometry)
                attrs.append(rdObj.Attributes)

            if bAddAxesBlock:
                idef_Axes = sc.doc.InstanceDefinitions.Find(sAxesBlockName)
                gAxes = sc.doc.Objects.AddInstanceObject(
                    idef_Axes.Index,
                    instanceXform=rg.Transform.Identity,
                    attributes=sc.doc.CreateDefaultAttributes())
                rdAxes = sc.doc.Objects.FindId(gAxes)
                geoms.append(rdAxes.Geometry)
                attrs.append(sc.doc.CreateDefaultAttributes())
                sc.doc.Objects.Delete(gAxes, quiet=False)
            else:
                for Line, attr in zip(*spb_RefAxes.createLines_and_attributes(fLineLength)):
                    geoms.append(Line)
                    attrs.append(attr)
                if bPlanes:
                    for plane, attr in zip(*spb_RefAxes.createPlaneSurfaces_and_attributes(fLineLength)):
                        geoms.append(plane)
                        attrs.append(attr)
                if bOriginPt:
                    pt, attr = spb_RefAxes.createPoint3d_and_attributes_OriginPt()
                    geoms.append(pt)
                    attrs.append(attr)
                if bAxesPts:
                    for pt, attr in zip(*spb_RefAxes.createPoint3ds_and_attributes_Axes(fLineLength)):
                        geoms.append(pt)
                        attrs.append(attr)

            if not sc.doc.InstanceDefinitions.ModifyGeometry(idef.Index, geoms, attrs):
                if bEcho: print("Could not modify block instance, {}.".format(idef.Name))
    else:
        xforms = [_.Geometry().Xform for _ in objrefs]

        for xform in xforms:
            if bAddAxesBlock:
                rs.InsertBlock2(sAxesBlockName, xform)
            else:
                gGeoms_Out = spb_RefAxes.addLineCurves(fLineLength)
                if bPlanes:
                    gGeoms_Out.extend(spb_RefAxes.addSurfaceObjects_Planes(fLineLength))
                if bOriginPt:
                    gGeoms_Out.append(spb_RefAxes.addPointObject_Origin())
                if bAxesPts:
                    gGeoms_Out.extend(spb_RefAxes.addPointObjects_Axes(fLineLength))
                if xform != rg.Transform.Identity:
                    rs.TransformObjects(gGeoms_Out, xform)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()