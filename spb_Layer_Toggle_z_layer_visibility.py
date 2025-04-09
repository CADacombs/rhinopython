"""
This script will toggle the visibility of all the layers whose names begin with
'z' and a digit to a common state opposite of the current state of the majority.

For example, if 10 'z' layers are on and 8 'z' layers are off, all are set to off.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250409: Created.

GetPersistentVisibility is missing from rhinoscriptsyntax, at least through 8.17.
In the context of this script, its result is correct versus IsVisible.

IsVisibleInUserInterface is only for the display of the layer name in the Layer panel,
not the visibilty state of the layer.
"""

import Rhino
import scriptcontext as sc


def _get_layers(layers_to_check):
    list_z_layers_On = []
    list_z_layers_Off = []
    list_z_layers_Visible = []
    list_z_layers_NotVisible = []

    for layer in layers_to_check:
        if layer.IsReference: continue
        if layer.IsDeleted: continue
        sPath = layer.FullPath
        sLayerName = sPath.split('::')[-1]
        if sLayerName.lower()[0] != 'z': continue
        if not sLayerName[1].isdigit(): continue

        #if layer.FullPath == "#":
        #    sEval = "layer.FullPath"; print(sEval,'=',eval(sEval))
        #    sEval = "layer.IsVisible"; print(sEval,'=',eval(sEval))
        #    sEval = "layer.GetPersistentVisibility()"; print(sEval,'=',eval(sEval))

        if layer.GetPersistentVisibility():
            list_z_layers_On.append(layer)
        else:
            list_z_layers_Off.append(layer)
        if layer.IsVisible:
            list_z_layers_Visible.append(layer)
        else:
            list_z_layers_NotVisible.append(layer)
    
    return (
        list_z_layers_On,
        list_z_layers_Off,
        list_z_layers_Visible,
        list_z_layers_NotVisible,
        )


def _create_report(list_z_layers_On, list_z_layers_Off, list_z_layers_Visible, list_z_layers_NotVisible):
    nOn = len(list_z_layers_On)
    nOff = len(list_z_layers_Off)
    
    if (
        (nOn == len(list_z_layers_Visible)) and
        (nOff == len(list_z_layers_NotVisible))
    ):
        return "{} 'z'n layers are on, {} are off.".format(nOn, nOff)
    else:
        return("{} 'z'n layers are on, {} are off, although only {} are visible"
            " due to parent visibility state.".format(
                nOn, nOff, len(list_z_layers_Visible)))


def main():

    bDocWasAlreadyModified = sc.doc.Modified

    rvs = _get_layers(sc.doc.Layers)
    if rvs is None: return
    (
        list_z_layers_On_Start,
        list_z_layers_Off_Start,
        list_z_layers_Visible_Start,
        list_z_layers_NotVisible_Start,
        ) = rvs

    if not (list_z_layers_On_Start or list_z_layers_Off_Start):
        print("No 'z' layers in document.")
        return
    
    sReport = _create_report(
        list_z_layers_On_Start,
        list_z_layers_Off_Start,
        list_z_layers_Visible_Start,
        list_z_layers_NotVisible_Start,
        )
    if sReport: print("Start:", sReport)

    nOn = len(list_z_layers_On_Start)
    nOff = len(list_z_layers_Off_Start)

    sc.doc.Views.RedrawEnabled = False

    if nOn < nOff: # Notice that ties also result in setting visiblity to False.
        #print("Will turn all on.")
        for layer in list_z_layers_Off_Start:
            layer.IsVisible = True
            #sc.doc.Layers.ForceLayerVisible(layer.Id)
    else:
        #print("Will turn all off.")
        for layer in list_z_layers_On_Start:
            layer.IsVisible = False
            #if Rhino.RhinoApp.ExeVersion == 8: # and layer.PersistentVisibility:
            #    sEval = "layer.FullPath"; print(sEval,'=',eval(sEval))
            #    sEval = "layer.PersistentVisibility"; print(sEval,'=',eval(sEval))
            layer.SetPersistentVisibility(False)

    sc.doc.Views.RedrawEnabled = True

    rvs = _get_layers(list_z_layers_On_Start + list_z_layers_Off_Start)
    if rvs is None: return
    (
        list_z_layers_On_End,
        list_z_layers_Off_End,
        list_z_layers_Visible_End,
        list_z_layers_NotVisible_End,
        ) = rvs

    sReport = _create_report(
        list_z_layers_On_End,
        list_z_layers_Off_End,
        list_z_layers_Visible_End,
        list_z_layers_NotVisible_End,
        )
    if sReport: print("End:", sReport)

    if not bDocWasAlreadyModified: sc.doc.Modified = False


if __name__ == '__main__': main()