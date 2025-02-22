"""
This script adds 3 lines, and optionally, 3 planes, and 4 points, representing
the X, Y, and Z axes aligned to the current CPlane or WorldXY plane.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
150909: Created.
170529: Changed default of Block option to No.
241007-08: Streamlined input. Refactored. Added options.
241017: Fixed typo.
241027: Split an option into two.
250219: Replaced most functions with those of an import.
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


    key = 'bAddAxesBlock'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCPlane_Not_World'; keys.append(key)
    values[key] = True
    names[key] = 'Plane'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'ViewCPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fLineLength'; keys.append(key)
    values[key] = Rhino.RhinoMath.UnitScale(
            Rhino.UnitSystem.Inches, sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPlanes'; keys.append(key)
    values[key] = False
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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def setOptions():
    """
    Returns:
        bool indicating whether to proceed
    """

    go = ri.Custom.GetOption()

    go.SetCommandPrompt("Set options <Enter to create axes>")

    go.AcceptNothing(True)

    def addOption(key):
        Opts.addOption(go, key)
        keys_in_order_added.append(key)

    while True:
        go.ClearCommandOptions()
        keys_in_order_added = [None] # None is just a placeholder since option indices are base 1.

        addOption('bAddAxesBlock')
        addOption('bCPlane_Not_World')
        addOption('fLineLength')
        addOption('bPlanes')
        addOption('bOriginPt')
        addOption('bAxesPts')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return False

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return True

        Opts.setValue(keys_in_order_added[go.Option().Index])


def main():

    if not setOptions():
        return

    bAddAxesBlock = Opts.values['bAddAxesBlock']
    bCPlane_Not_World = Opts.values['bCPlane_Not_World']
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


    # To avoid unnecessary transforms.
    if bCPlane_Not_World:
        plane_CPlane = rs.ViewCPlane()
        bCPlane_Not_World = plane_CPlane != rg.Plane.WorldXY

    if bAddAxesBlock:
        if bCPlane_Not_World:
            xform = rs.XformChangeBasis(rs.ViewCPlane(), rs.WorldXYPlane())
        else:
            xform = rs.XformIdentity()
        rs.InsertBlock2(sAxesBlockName, xform)
    else:
        gGeoms_Out = spb_RefAxes.addLineCurves(fLineLength)
        if bPlanes:
            gGeoms_Out.extend(spb_RefAxes.addSurfaceObjects_Planes(fLineLength))
        if bOriginPt:
            gGeoms_Out.append(spb_RefAxes.addPointObject_Origin())
        if bAxesPts:
            gGeoms_Out.extend(spb_RefAxes.addPointObjects_Axes(fLineLength))
        if bCPlane_Not_World:
            xform = rs.XformChangeBasis(rs.ViewCPlane(), rs.WorldXYPlane())
            rs.TransformObjects(gGeoms_Out, xform)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()