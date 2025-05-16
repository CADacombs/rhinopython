"""
This script is to be used both as a module and a stand-along script.
Unlike _BoundingBox, the script offers an Accurate option.
Unlike the script, _BoundingBox has an Output object type option.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
...
190505: Removed a function.
200525-26: Added functionality so this script can be run on its own or continued being used as a module.
250515: Added support to create bounding box per CPlane. Added command options.

TODO:
    Merge spb_BoundingBox2 with this script. It includes options for center point, etc.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bAccurate'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCPlane_notWorld'; keys.append(key)
    values[key] = True
    names[key] = 'CoordinateSystem'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'World', 'CPlane')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bCumulative'; keys.append(key)
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

        #if key == 'fTol_Edge_and_trim':
        #    if cls.riOpts[key].CurrentValue < 0.0:
        #        cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
        #    elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
        #        cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        #    sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
        #    return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList
            return

        print("Invalid key?")


def getInput():
    """
    Get objects.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select objects to frame with a box")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Brep

    #go.AcceptNumber(True, acceptZero=False)
    go.EnableClearObjectsOnEntry(False) # If not set to False, faces will be unselected when result == ri.GetResult.Object 

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()
        idxs_Opts.clear()

        addOption('bAccurate')
        addOption('bCPlane_notWorld')
        addOption('bCumulative')
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

        if res == ri.GetResult.Number:
            key = 'fTol_Edge_and_trim'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createInverseTransform(xform_In):
    bSuccess, xform_Out = xform_In.TryGetInverse()
    return xform_Out if bSuccess else None


def createBoundingBox(rgObjs, bAccurate=True, xform=rg.Plane.WorldXY):
    """
    Returns:
        rg.BoundingBox as transformed per xform.
    """

    try: rgObjs = list(rgObjs)
    except: rgObjs = [rgObjs]

    bbox = rg.BoundingBox.Empty

    if xform is None or xform == rg.Plane.WorldXY:
        for geom in rgObjs:
            bbox.Union(geom.GetBoundingBox(accurate=bAccurate))
    else:
        for geom in rgObjs:
            dup = geom.Duplicate()
            dup.Transform(xform)
            #sc.doc.Objects.Add(dup); #1/0
            bbox.Union(dup.GetBoundingBox(accurate=bAccurate))
            #sEval = "geom.IsDocumentControlled"; print(sEval, '=', eval(sEval))
            dup.Dispose()
        sEval = "bbox"; print(sEval, '=', eval(sEval))
        #sc.doc.Objects.AddPoints(bbox.GetCorners())
        xform_Inverse = createInverseTransform(xform)
        if xform_Inverse is None:
            raise Exception("Could not get inverse transform.")
        #bbox.Transform(xform_Inverse)
        #sc.doc.Objects.AddPoints(bbox.GetCorners())
        sEval = "bbox"; print(sEval, '=', eval(sEval))

    if bbox is None:
        raise ValueError("Bounding box not created.")

    return bbox


def epsilonEquals(bBoxA, bBoxB, fTol=sc.doc.ModelAbsoluteTolerance):
    """
    Created because BoundingBox structure does not contain an EpsilonEquals method.
    """
    boxA = Rhino.Geometry.Box(bBoxA)
    if not boxA: return
    boxB = Rhino.Geometry.Box(bBoxB)
    if not boxB: return
    return boxA.EpsilonEquals(boxB, fTol)


def findMatchingBoundingBox(bboxes_toSearch, bbox_toMatch, fTolerance=None, bDebug=False):
    """
    Parameters:
        bboxes_toSearch: Iterable of BoundingBoxes from which to search.
        bbox_toMatch: BoundingBox to match.
        fTolerance
        bDebug
    Return on success:
        int(Index of single matching bounding box)
    Return on fail:
        None
    """

    if bDebug: print('findMatchingBoundingBox()...')

    if fTolerance is None:
        fTolerance = sc.doc.ModelAbsoluteTolerance

    idx_Found = None

    for i, bbox_toCheck in enumerate(bboxes_toSearch):
        if (
            bbox_toMatch.Min.EpsilonEquals(bboxes_toSearch[i].Min, epsilon=fTolerance) and
            bbox_toMatch.Max.EpsilonEquals(bboxes_toSearch[i].Max, epsilon=fTolerance)
        ):
            if idx_Found is not None:
                print("More than one match found. Check results.")
                return
            idx_Found = i
        
    if idx_Found is None:
        if bDebug:
            print("Matching bounding box NOT found.")
            #sc.doc.Objects.AddBrep(bbox_toMatch.ToBrep())
            #for bb in bboxes_toSearch:
                #sc.doc.Objects.AddBrep(bb.ToBrep())

    return idx_Found


def addDocObject_of_boundingBox(bbox, xform_for_bbbox, bEcho=True, bDebug=False):
    iCt_MatchingCoords = 0
    bCoordMatches = [False, False, False]
    if abs(bbox.Min.X - bbox.Max.X) < 1e-9:
        iCt_MatchingCoords += 1
        bCoordMatches[0] = True
    if abs(bbox.Min.Y - bbox.Max.Y) < 1e-9:
        iCt_MatchingCoords += 1
        bCoordMatches[1] = True
    if abs(bbox.Min.Z - bbox.Max.Z) < 1e-9:
        iCt_MatchingCoords += 1
        bCoordMatches[2] = True

    iCt_MatchingCoords = bCoordMatches.count(True)

    if bDebug: print(bbox)

    if iCt_MatchingCoords == 0:
        box = rg.Box(bbox=bbox)
        if xform_for_bbbox is not None:
            box.Transform(xform_for_bbbox)
        gBox = sc.doc.Objects.AddBox(box)
        if gBox == Guid.Empty:
            if bEcho: print("Box could not be added.")
            return
        else:
            if bEcho: print("Box (extrusion) was added.")
            return gBox

    if iCt_MatchingCoords == 1:
        if bDebug: print("Bounding box is planar.")
        if bCoordMatches[0]:
            plane = rg.Plane.WorldYZ
            plane.Translate(rg.Vector3d(bbox.Min.X, 0.0, 0.0))
        elif bCoordMatches[1]:
            plane = rg.Plane.WorldZX
            plane.Translate(rg.Vector3d(0.0, bbox.Min.Y, 0.0))
        elif bCoordMatches[2]:
            plane = rg.Plane.WorldXY
            plane.Translate(rg.Vector3d(0.0, 0.0, bbox.Min.Z))
        rect = rg.Rectangle3d(
            plane=plane,
            cornerA=bbox.Min,
            cornerB=bbox.Max)
        if xform_for_bbbox is not None:
            rect.Transform(xform_for_bbbox)
        gPLCrv = sc.doc.Objects.AddRectangle(rect)
        if gPLCrv == Guid.Empty:
            if bEcho: print("Polyline could not be added.")
            return
        else:
            if bEcho: print("Polyline was added.")
            return gPLCrv

    if iCt_MatchingCoords == 2:
        if bEcho: print("Bounding box is linear.")
        line = rg.Line(bbox.Min, to=bbox.Max)
        if xform_for_bbbox is not None:
            line.Transform(xform_for_bbbox)
        gLineCrv = sc.doc.Objects.AddLine(line)
        if gLineCrv == Guid.Empty:
            if bEcho: print("Line could not be added.")
            return
        else:
            sc.doc.Views.Redraw()
            if bEcho: print("Line was added.")
            return gLineCrv

    if iCt_MatchingCoords == 3:
        if bEcho: print("Bounding box is a point. Point was not added.")
        return

    raise Exception("What happened?")


def main():

    #res, objrefs = ri.RhinoGet.GetMultipleObjects(
    #    "Select objects to create bounding box",
    #    acceptNothing=False,
    #    filter=rd.ObjectType.AnyObject)
    #if res == Rhino.Commands.Result.Cancel: return

    objrefs = getInput()
    if objrefs is None: return

    bAccurate = Opts.values['bAccurate']
    bCPlane_notWorld = Opts.values['bCPlane_notWorld']
    bCumulative = Opts.values['bCumulative']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    xform_WtoC = xform_CtoW = None
    if bCPlane_notWorld:
        plane = sc.doc.Views.ActiveView.ActiveViewport.ConstructionPlane()
        if plane != rg.Plane.WorldXY:
            xform_WtoC = rg.Transform.ChangeBasis(rg.Plane.WorldXY, plane)
            xform_CtoW = createInverseTransform(xform_WtoC)
            if xform_CtoW is None:
                raise Exception("Could not get inverse of Cplane-to-World transform.")

    sc.doc.Objects.UnselectAll()

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gResults = []

    if bCumulative:
        rgObjs = [o.Geometry() for o in objrefs]
        bbox = createBoundingBox(rgObjs, bAccurate=bAccurate, xform=xform_WtoC)
        if bbox is None:
            if bEcho: print("Could not create bounding box.")
            return
        rv = addDocObject_of_boundingBox(
            bbox=bbox,
            xform_for_bbbox=xform_CtoW,
            bEcho=bEcho,
            bDebug=bDebug)
        if rv:
            gResults.append(rv)
    else:
        for objref in objrefs:
            rgObjs = [objref.Geometry()]
            bbox = createBoundingBox(rgObjs, bAccurate=bAccurate, xform=xform_WtoC)
            if bbox is None:
                if bEcho: print("Could not create bounding box.")
                continue
            rv = addDocObject_of_boundingBox(
                bbox=bbox,
                xform_for_bbbox=xform_CtoW,
                bEcho=bEcho,
                bDebug=bDebug)
            if rv:
                gResults.append(rv)

        #box = rg.Box(bbox=bbox)
        #box.Transform(xform_CtoW)
        #sc.doc.Objects.AddBox(box)
        #return

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()