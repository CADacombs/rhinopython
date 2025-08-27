"""
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
200701: Created as a union of other scripts.
210825: Added functionality to get the target color from a face.
231111: Modified some printed output.
250826: Modified a selectionfilter.
"""

import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


sMatchTypes = (
    'OtherLayer',
    'ItsFirstObject',
    'ItsObject',
    'Object',
    'Face',
    )


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    idxOpt = {}
    stickyKeys = {}


    key = 'iDefaultOpt'; keys.append(key)
    values[key] = 2
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

        cls.idxOpt[key] = None

        if key in cls.riOpts:
            if key[0] == 'b':
                cls.idxOpt[key] = go.AddOptionToggle(
                        cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'f':
                cls.idxOpt[key] = go.AddOptionDouble(
                    cls.names[key], cls.riOpts[key])[0]
            elif key[0] == 'i':
                cls.idxOpt[key] = go.AddOptionInteger(
                    englishName=cls.names[key], intValue=cls.riOpts[key])[0]
        else:
            cls.idxOpt[key] = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return cls.idxOpt[key]


    @classmethod
    def setValues(cls):
        for key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
            if key == 'fTol':
                if cls.riOpts[key].CurrentValue < 0.0:
                    cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_MatchType():
    """
    """

    go = ri.Custom.GetOption()

    go.SetCommandPrompt("Match layer color to")
    go.SetCommandPromptDefault(sMatchTypes[Opts.values['iDefaultOpt']])

    go.AcceptNothing(True)

    for sMatchType in sMatchTypes:
        go.AddOption(sMatchType)

    res = go.Get()

    if res == ri.GetResult.Option:
        idx_Base0 = go.Option().Index - 1
        Opts.values['iDefaultOpt'] = idx_Base0
        Opts.saveSticky()
        go.Dispose()
        return idx_Base0
    
    if res == ri.GetResult.Nothing:
        go.Dispose()
        return Opts.values['iDefaultOpt']


    go.Dispose()
    return


def matchToOtherLayer():

    sLayer_Source = rs.GetLayer("Select color source layer")
    if sLayer_Source is None: return

    color = rs.LayerColor(sLayer_Source)

    sLayers_toMod = rs.GetLayers("Select layers to edit")
    if sLayers_toMod is None: return

    for sLayer in sLayers_toMod:
        print("{} color changed from {} to {}.".format(
            sLayer,
            rs.LayerColor(sLayer, color=color),
            rs.LayerColor(sLayer)
            ))


def blockObjectOnLayer(sLayer):
    for sBlock in rs.BlockNames(sort=True):
        for g in rs.BlockObjects(sBlock):
            if rs.ObjectLayer(g) == sLayer:
                return g


def matchTo1stObj():
    
    bs = rs.GetBoolean("Match Layer Colors to Their 1st Object",
    ("AllLayers", "No", "Yes"),
    (False))
    if bs is None: return
    
    bAllLayers = bs[0]
    
    if bAllLayers:
        sLayers = rs.LayerNames()
    else:
        sLayers = rs.GetLayers("Select layer names to set the color of each to that of their 1st object")
        if sLayers is None: return
    
    sc.doc.Views.RedrawEnabled = False
    
    for sLayer in sLayers:
        if sLayer == "Default":
            print("Skipped color modification of Default.")
            continue
        gObjs_on_layer = rs.ObjectsByLayer(sLayer)
        if not gObjs_on_layer:
            g = blockObjectOnLayer(sLayer)
            if g is None:
                print("No objects found on {}.".format(sLayer))
                continue # to next layer.
        else:
            g = gObjs_on_layer[0]
        color = rs.ObjectColor(g)
        rs.LayerColor(sLayer, color)
        print(str(color) + " assigned to layer " + str(sLayer))
    
    rs.EnableRedraw(True)


def stringOfColor(color):
    knownColor = color.ToKnownColor()
    if str(knownColor) == '0':
        return "({} {} {})".format(color.R, color.G, color.B)
    else:
        return "{} ({} {} {})".format(knownColor, color.R, color.G, color.B)


def matchToItsObj():
    """
    Object is post-picked.
    """

    # Why was instance not included?
    #gObjs = rs.GetObjects("Select objects to obtain their colors",
    #        filter= 4294967295 - rs.filter.instance,
    #        preselect=True) # 4294967295 is value of Rhino.DocObjects.ObjectType.AnyObject.

    gObjs = rs.GetObjects("Select objects to obtain their colors",
            filter=0,
            preselect=True)
    if gObjs is None: return
    
    sc.doc.Views.RedrawEnabled = False
    
    for gObj in gObjs:
        colorObj = rs.ObjectColor(gObj)
        sLayer = rs.ObjectLayer(gObj)
        if rs.LayerColor(sLayer) != colorObj:
            colorLayer0 = rs.LayerColor(sLayer, colorObj)
            colorLayer1 = rs.LayerColor(sLayer)
            print("Changed color of {}".format(sLayer))
            print("  from " + stringOfColor(colorLayer0))
            print("  to   " + stringOfColor(colorLayer1))
    
    sc.doc.Views.RedrawEnabled = True


def matchToObj():
    
    gObj = rs.GetObject("Select object to obtain its color",
            filter= 4294967295 - rs.filter.instance,
            preselect=True) # 4294967295 is value of Rhino.DocObjects.ObjectType.AnyObject.
    if gObj is None: return
    
    sLayers = rs.GetLayers("Select layer(s) for color change")
    if sLayers is None: return False
    
    colorObj = rs.ObjectColor(gObj)

    sc.doc.Views.RedrawEnabled = False

    sLayers_AlreadyHaveTargetColor = []
    sLayers_AssignedTargetColor = []
    sLayers_FailedToBeAssignedTargetColor = []

    for sLayer in sLayers:
        if rs.LayerColor(sLayer) == colorObj:
            sLayers_AlreadyHaveTargetColor.append(sLayer)
        else:
            colorLayer0 = rs.LayerColor(sLayer, colorObj)
            colorLayer1 = rs.LayerColor(sLayer)
            if colorLayer1 == colorObj:
                sLayers_AssignedTargetColor.append(sLayer)
            else:
                sLayers_FailedToBeAssignedTargetColor.append(sLayer)

    if sLayers_AlreadyHaveTargetColor:
        print("Layers skipped because they already have color {}:".format(stringOfColor(colorObj)))
        for sLayer in sLayers_AlreadyHaveTargetColor:
            print(sLayer)

    if sLayers_AssignedTargetColor:
        print("Color {} assigned to layers:".format(stringOfColor(colorObj)))
        for sLayer in sLayers_AssignedTargetColor:
            print(sLayer)

    if sLayers_FailedToBeAssignedTargetColor:
        print("Color {} could not be assigned to layers:".format(stringOfColor(colorObj)))
        for sLayer in sLayers_FailedToBeAssignedTargetColor:
            print(sLayer)

    sc.doc.Views.RedrawEnabled = True


def matchToFace():
    """
    """

    objref = rs.GetObject(
        "Select face to obtain its color",
        filter=rs.filter.surface,
        preselect=True,
        subobjects=True
        )
    if objref is None: return


    sLayers = rs.GetLayers("Select layer(s) for color change")
    if sLayers is None: return False


    rgF = objref.Surface()
    color_of_face = rgF.PerFaceColor

    sc.doc.Views.RedrawEnabled = False

    sLayers_AlreadyHaveTargetColor = []
    sLayers_AssignedTargetColor = []
    sLayers_FailedToBeAssignedTargetColor = []

    for sLayer in sLayers:
        if rs.LayerColor(sLayer) == color_of_face:
            sLayers_AlreadyHaveTargetColor.append(sLayer)
        else:
            colorLayer0 = rs.LayerColor(sLayer, color_of_face)
            colorLayer1 = rs.LayerColor(sLayer)
            if colorLayer1 == color_of_face:
                sLayers_AssignedTargetColor.append(sLayer)
            else:
                sLayers_FailedToBeAssignedTargetColor.append(sLayer)

    if sLayers_AlreadyHaveTargetColor:
        print("Layers skipped because they already have color {}:".format(stringOfColor(color_of_face)))
        for sLayer in sLayers_AlreadyHaveTargetColor:
            print(sLayer)

    if sLayers_AssignedTargetColor:
        print("Color {} assigned to layers:".format(stringOfColor(color_of_face)))
        for sLayer in sLayers_AssignedTargetColor:
            print(sLayer)

    if sLayers_FailedToBeAssignedTargetColor:
        print("Color {} could not be assigned to layers:".format(stringOfColor(color_of_face)))
        for sLayer in sLayers_FailedToBeAssignedTargetColor:
            print(sLayer)

    sc.doc.Views.RedrawEnabled = True


def main():

    idx = getInput_MatchType()
    if idx is None: return


    funcsPerType = (
        matchToOtherLayer,
        matchTo1stObj,
        matchToItsObj,
        matchToObj,
        matchToFace,
        )

    funcsPerType[idx]()


if __name__ == '__main__': main()