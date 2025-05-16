"""
This script is to be used more as a module by other scripts.
Running this script directly provides a command for debugging this script.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
200526: Created.
250515: Import-related updates.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid

import spb_BoundingBox
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


    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    riAddOpts[key] = addOptionDouble(key, names, riOpts)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bDelOtherFacesNotExtract'; keys.append(key)
    values[key] = False
    names[key] = 'Action'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Extract', 'DeleteOtherFaces')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
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


def getInput_GeomForBBoxToMatch():
    """
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select objects to create bounding box to match")

    #go.GeometryFilter = rd.ObjectType

    go.GroupSelect = True

    go.SubObjectSelect = False

    idxs_Opts = {}

    while True:
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDelOtherFacesNotExtract'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return [objrefs] + [Opts.values[key] for key in Opts.keys]
        else:
            # An option was selected or a number was entered.
            key = 'fTolerance'
            if res == ri.GetResult.Number:
                Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 1e-12:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInput_BrepWithFacesToSearch():
    """
    Get Brep Surface with optional input
    """

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps to search")

    go.GeometryFilter = rd.ObjectType.Brep

    idxs_Opts = {}

    while True:
        key = 'fTolerance'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDelOtherFacesNotExtract'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bEcho'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]
        key = 'bDebug'; idxs_Opts[key] = Opts.riAddOpts[key](go)[0]

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return [objref] + [Opts.values[key] for key in Opts.keys]
        else:
            # An option was selected or a number was entered.
            key = 'fTolerance'
            if res == ri.GetResult.Number:
                Opts.riOpts[key].CurrentValue = go.Number()
            if Opts.riOpts[key].CurrentValue < 1e-12:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def findFace(rgBrep_toSearch, rgObjs_forBBoxToMatch, fTolerance=None, bDebug=False):
    """
    Returns:
        Index of brep in breps.

    This function is used by other scripts, so leave it at the root level.
    """

    if fTolerance is None:
        fTolerance = 2.0 * sc.doc.ModelAbsoluteTolerance

    bbox_toMatch = spb_BoundingBox.createBoundingBox(rgObjs_forBBoxToMatch)

    #sc.doc.Objects.AddBox(rg.Box(bbox=bbox_toMatch)); sc.doc.Views.Redraw()

    bboxes_Faces_perBrep = []
    rgB = rgBrep_toSearch

    for idxF, rgFace in enumerate(rgB.Faces):
        rgEs_forBorder = []

        # This only uses the OutLoop.
        #for rgT in rgFace.OuterLoop.Trims:
        #    if rgT != rg.BrepTrimType.Seam:
        #        rgEs_forBorder.append(rgT.Edge)

        # Accumulate the edge curves of all non-seam trims.
        for idxE in rgFace.AdjacentEdges():
            rgE = rgB.Edges[idxE]
            idxsTs = rgE.TrimIndices()
            if rgB.Trims[idxsTs[0]].TrimType != rg.BrepTrimType.Seam:
                rgEs_forBorder.append(rgE)

        bbox_ofFaceBorder = spb_BoundingBox.createBoundingBox(rgEs_forBorder)
        bboxes_Faces_perBrep.append(bbox_ofFaceBorder)

        #sc.doc.Objects.AddBox(rg.Box(bbox=bbox_ofFaceBorder)); sc.doc.Views.Redraw()

    idxF_Match = spb_BoundingBox.findMatchingBoundingBox(
        bboxes_Faces_perBrep,
        bbox_toMatch,
        fTolerance=fTolerance,
        bDebug=bDebug)

    return idxF_Match


def processBrepObject(rhBrep_toSearch, rhObjs_forBBoxToMatch, **kwargs):
    """
    Parameters:
        fTolerance
        bDelOtherFacesNotExtract
        bEcho
        bDebug
    Returns on success:
        GUID of monoface brep
    Returns on fail:
        None
    """

    if not rhBrep_toSearch or rhObjs_forBBoxToMatch is None: return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    bDelOtherFacesNotExtract = getOpt('bDelOtherFacesNotExtract')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    try: rhBrep_toSearch = list(rhBrep_toSearch)
    except: rhBrep_toSearch = [rhBrep_toSearch]

    try: rhObjs_forBBoxToMatch = list(rhObjs_forBBoxToMatch)
    except: rhObjs_forBBoxToMatch = [rhObjs_forBBoxToMatch]


    geom_forBBox_toMatch = [rs.coercegeometry(o) for o in rhObjs_forBBoxToMatch]

    gB_toSearch = rs.coerceguid(rhBrep_toSearch)
    rgB = rs.coercebrep(rhBrep_toSearch)

    idxFace_Match = findFace(rgB, geom_forBBox_toMatch)

    if idxFace_Match is None:
        print("Matching face was not found.")
        return

    if bDelOtherFacesNotExtract:
        rdBrep = rs.coercerhinoobject(rhBrep_toSearch)
        rgB_1F = rgB.Faces[idxFace_Match].DuplicateFace(duplicateMeshes=True)
        gB_1F = sc.doc.Objects.AddBrep(
            rgB_1F,
            attributes=rdBrep.Attributes)
        rgB_1F.Dispose()
        if gB_1F == Guid.Empty:
            print("Brep could not be added.")
            return

        if sc.doc.Objects.Delete(objectId=gB_toSearch, quiet=False):
            print("Monoface brep replaced polyface brep.")
        return gB_1F
    else:
        rc = xBrepObject.extractFaces(gB_toSearch, [idxFace_Match])
        if rc is None: return
        gExtracted, gRemaining = rc
        return gExtracted[0]


def main():

    rc = getInput_GeomForBBoxToMatch()
    if rc is None: return
    objrefs_forBBox_toMatch = rc[0]

    rc = getInput_BrepWithFacesToSearch()
    if rc is None: return
    objref_Brep_toSearch = rc[0]

    if not Opts.values['bDebug']:
        sc.doc.Views.RedrawEnabled = False

    rc = processBrepObject(
            rhBrep_toSearch=objref_Brep_toSearch,
            rhObjs_forBBoxToMatch=objrefs_forBBox_toMatch)

    if not rc: return

    #if isinstance(rc, int):
    #    idxF = rc
    #    sc.doc.Objects.UnselectAll()
    #    if Rhino.RhinoApp.ExeVersion >= 6:
    #        rdBrep = rs.coercerhinoobject(objref_Brep_toSearch)

    #        compIdx = rg.ComponentIndex(
    #            rg.ComponentIndexType.BrepFace,
    #            index=idxF)
    #        rdBrep.SelectSubObject(
    #            componentIndex=compIdx,
    #            select=True,
    #            syncHighlight=True,
    #            persistentSelect=True)
    
    if isinstance(rc, Guid):
        gB_Extracted = rc
        sc.doc.Objects.Select(objectId=gB_Extracted)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()