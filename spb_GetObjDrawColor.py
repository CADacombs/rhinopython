"""
Per https://developer.rhino3d.com/api/rhinocommon/rhino.docobjects.objectattributes/drawcolor
"object's draw color, which is based on the object's color source"

This is not the ObjectColor of the RhinoObject when ColorSource is one of
ColorFromLayer, ColorFromMaterial, or ColorFromParent
"""

#! python 2  Must be on a line less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
150830: Created, starting with SetObjectColor.py.
170926: Changed file name from getObjectColor.py to listObjectColor.py.
171003: Added ColorSource to output.
250922: Now, prints name of color in NamedColorList instead of from System.Drawing.Color.
"""

import Rhino
import rhinoscriptsyntax as rs


def getStringFromColor(color):

    sKnownColor = color.ToKnownColor().ToString()

    # '0' is returned from ToKnownColor when color is not named.
    if sKnownColor == '0':
        return "({} {} {})".format(color.R, color.G, color.B)

    return sKnownColor


def _tryFindRhinoColorName(color):
    for rh_UI_NamedColor in Rhino.UI.NamedColorList.Default:
        if rh_UI_NamedColor.Color == color:
            return "{} ({} {} {})".format(
                rh_UI_NamedColor.Name,
                color.R,
                color.G,
                color.B)
    return "({} {} {})".format(
        color.R,
        color.G,
        color.B)


def main():

    gObjs = rs.GetObjects("Select objects to list their colors", preselect=True, select=True)
    if gObjs is None: return

    colorsFound = []

    for gObj in gObjs:
        rdObj = rs.coercerhinoobject(gObj)

        colorDisplayed = rs.ObjectColor(gObj)
        sColorDisplayed = _tryFindRhinoColorName(colorDisplayed)
        #sColorDisplayed = getStringFromColor(colorDisplayed)

        colorSource = rdObj.Attributes.ColorSource.ToString()

        if (sColorDisplayed, colorSource) not in colorsFound:
            colorsFound.append((sColorDisplayed, colorSource))

    for a,b in colorsFound:
        print(a, b)


if __name__ == '__main__': main()