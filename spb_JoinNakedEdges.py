"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160303: Created.  'JoinNakedEdges' apparently is an undocumented command.
160304: Changed minimum tolerance allowed to be entered to 0.  Changed formatting of printed decimals.
160728: Modularized and added more printed output.
190625: Added Opts.  Refactored getInput.
190822: Modified an option default value.
200630: Import-related update and modified an optino default value.
230928: Bug fixes.  Added a command option.  Refactored.
231001: Refactored.  Added text (indices of bad brep components) to dots.
250215: '+' now is included in reporting of changes in edge counts.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System.Drawing import Color

#import re


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}

    key = 'fTolerance'; keys.append(key)
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMarkBadGeom'; keys.append(key)
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
        if sc.sticky.has_key(stickyKeys[key]):
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue < 0:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        if key in cls.stickyKeys:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getAllNormalBreps():
    oes = rd.ObjectEnumeratorSettings()
    oes.NormalObjects = True
    oes.LockedObjects = False # Default is True.
    oes.IncludeLights = False
    oes.IncludeGrips = False
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))


def getInput():
    """
    Get brep face.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select polysrf breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.PolysrfFilter
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.OpenPolysrf
    go.SubObjectSelect = False

    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.EnableClearObjectsOnEntry(False) # Do not clear objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False) # Do not unselect object when an option selected, a number is entered, etc.

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fTolerance')
        addOption('bMarkBadGeom')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return getAllNormalBreps()


        if res == ri.GetResult.Number:
            key = 'fTolerance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def countsOfEdgeValences(rgBrep0):
    numEdgeValence = [0]*4
    for e in rgBrep0.Edges:
        numEdgeValence[e.Valence.value__] += 1
        e.Dispose()
    return numEdgeValence 


def getMaxDeviationToJoin(rgBrep0, fMaxDevToFlag):
    rgCrvs = rgBrep0.DuplicateNakedEdgeCurves(True, True)
    rgBrep0.Dispose()
    
    rs.Prompt(message="Checking deviations of naked edges ...")
    
    # Check all combinations of curves for maximum deviation.
    fMaxDev = 0.
    for i in range(len(rgCrvs)):
        for j in range(i + 1, len(rgCrvs)):
            rc = rg.Curve.GetDistancesBetweenCurves(
                    rgCrvs[i], rgCrvs[j], sc.doc.ModelAbsoluteTolerance)
            if rc[0]:
                fDev = rc[1]
                if fDev <= fMaxDevToFlag:
                    if fDev > fMaxDev: fMaxDev = fDev
        rgCrvs[i].Dispose()
    
    if fMaxDev:
        print("Largest deviation <= {:f}: {:f}".format(fMaxDevToFlag, fMaxDev))
        fTolerance = fMaxDev
    else:
        print("No deviations <= {:f} exist.".format(fMaxDevToFlag))
        fTolerance = 2. * fMaxDevToFlag
    
    return rs.GetReal("Enter maximum deviation to join", fTolerance, 0.)


def _dotEdge(rgBrep, idxEdge):
    rgEdge = rgBrep.Edges[idxEdge]
    pt = rgEdge.PointAtNormalizedLength(0.5)
    if pt.X == Rhino.RhinoMath.UnsetValue:
        pt = rgEdge.PointAtStart

    rgDot = rg.TextDot(text="E{}".format(idxEdge), location=pt)
    rgDot.FontHeight = 11 if Rhino.RhinoApp.ExeVersion >= 6 else 14

    attr = rd.ObjectAttributes()
    attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex # If not done, layer of index 0 will be used.
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr.ObjectColor = Color.Red

    gDot = sc.doc.Objects.AddTextDot(rgDot, attr)
    
    if gDot != gDot.Empty: return rgDot

    print("AddTextDot failed for E[{}}.".format(idxEdge))


def _dotFace(rgBrep, idxFace):

    rgFace = rgBrep.Faces[idxFace]

    amp = rg.AreaMassProperties.Compute(rgFace)
    if amp is None:
        print("Surface not dotted because its AreaMassProperties,"
            "and thus its centroid, cannot be calculated.")
        return

    ptCentroid = amp.Centroid
    getrc, u, v = rgFace.ClosestPoint(ptCentroid)
    if rg.BrepFace.IsPointOnFace(u, v):
        location = rgFace.PointAt(u, v)
    else:
        pts = []
        for idxE in rgFace.AdjacentEdges():
            rgE = rgBrep.Edges[idxE]
            bSuccess, t = rgE.ClosestPoint(ptCentroid)
            if bSuccess: pts.append(rgE.PointAt(t))
        pt3dlist = Rhino.Collections.Point3dList(pts)
        location = pt3dlist.ClosestPointInList(ptCentroid)

    rgDot = rg.TextDot(text="F{}".format(idxFace), location=location)
    rgDot.FontHeight = 11 if Rhino.RhinoApp.ExeVersion >= 6 else 14

    attr = rd.ObjectAttributes()
    attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex # If not done, layer of index 0 will be used.
    attr.ColorSource = rd.ObjectColorSource.ColorFromObject
    attr.ObjectColor = Color.Red

    gDot = sc.doc.Objects.AddTextDot(rgDot, attr)

    if gDot != gDot.Empty: return rgDot

    print("AddTextDot failed for F[{}].".format(idxFace))


def processBrepObject(gBrep, fTolerance=None, bMarkBadGeom=True, bEcho=None, bDebug=None):
    """
    """

    if fTolerance is None: fTolerance = Opts.values['fTolerance']
    if bMarkBadGeom is None: bMarkBadGeom = Opts.values['bMarkBadGeom']
    if bEcho is None: bEcho = Opts.values['bEcho']
    if bDebug is None: bDebug = Opts.values['bDebug']

    rgBrepX = rs.coercebrep(gBrep)
 
    # Report face count and edge valence counts.
    numF_B0 = rgBrepX.Faces.Count
    numEValence_B0 = countsOfEdgeValences(rgBrepX)
    
    print("Counts: {} faces, {} none edges, {} naked edges,"
          "{} interior edges, {} non-manifold edges".format(
              numF_B0,
              numEValence_B0[0],
              numEValence_B0[1],
              numEValence_B0[2],
              numEValence_B0[3],
              ))
    #    print "Faces: {}".format(numF_B0)
    #    print "None edges: {}".format(numEValence_B0[0])
    #    print "Naked edges: {}".format(numEValence_B0[1])
    #    print "Interior edges: {}".format(numEValence_B0[2])
    #    print "Non-manifold edges: {}".format(numEValence_B0[3])
    
    #    # Get upper limit of naked edge deviation to flag.
    #    fTol = sc.doc.ModelAbsoluteTolerance
    #    fMaxDevToFlag = rs.GetReal(
    #            "Enter upper limit of naked edge deviation to flag", 10.*fTol, 0.)
    #    if fMaxDevToFlag is None:
    #        oRef.Dispose()
    #        return
    #    
    #    # Get maximum deviation of naked edges to join.
    #    fTolerance = getMaxDeviationToJoin(rgBrepX, fMaxDevToFlag)
    #    if fTolerance is None:
    #        rgBrepX.Dispose()
    #        oRef.Dispose()
    #        return
    
    Rhino.RhinoApp.SetCommandPrompt(prompt="Joining naked edges ...")
    
    # Join naked edges and replace polysurface with brep.
    if rgBrepX.JoinNakedEdges(fTolerance) > 0:
        if not rgBrepX.IsValid:
            if bEcho:
                if bMarkBadGeom:
                    print("New brep geometry is invalid.  Fixed marked area(s).")
                else:
                    print("New brep geometry is invalid.")
            sLog = rgBrepX.IsValidWithLog()[1]
            print(sLog)
            if bMarkBadGeom:
                ## Alternative
                #m = re.search(r"\[([A-Za-z0-9_]+)\]", sLog)
                #print(m.group(1))
                if 'ON_Brep.m_F[' in sLog:
                    idxF = int(
                            sLog.split('ON_Brep.m_F[', 1)[1].\
                            split('] is invalid.')[0])
                    _dotFace(rgBrepX, idxF)
                if 'ON_Brep.m_E[' in sLog:
                    idxE = int(
                            sLog.split('ON_Brep.m_E[', 1)[1].\
                            split('] is invalid.')[0])
                    _dotEdge(rgBrepX, idxE)
            return

        if not sc.doc.Objects.Replace(gBrep, rgBrepX):
            print("Brep was not replaced.")
            return False

        sc.doc.Views.Redraw()
        numF_B1 = rgBrepX.Faces.Count
        numEValence_B1 = countsOfEdgeValences(rgBrepX)
        
        print("Counts: Before -> After = Change")
        print("Faces: {} - {} = {:+}".format(
                numF_B0,
                numF_B1,
                numF_B1-numF_B0,
        ))
        print("None edges: {} -> {} = {:+}".format(
                numEValence_B0[0],
                numEValence_B1[0],
                numEValence_B1[0]-numEValence_B0[0],
        ))
        print("Naked edges: {} -> {} = {:+}".format(
                numEValence_B0[1],
                numEValence_B1[1],
                numEValence_B1[1]-numEValence_B0[1],
        ))
        print("Interior edges: {} -> {} = {:+}".format(
                numEValence_B0[2],
                numEValence_B1[2],
                numEValence_B1[2]-numEValence_B0[2],
        ))
        print("Non-manifold edges: {} -> {} = {:+}".format(
                numEValence_B0[3],
                numEValence_B1[3],
                numEValence_B1[3]-numEValence_B0[3],
        ))
        
        return rgBrepX


    print("No edges were joined within a deviation of {:f}.".format(fTolerance))
    return False

    return rgBrepX


def main():
    rhObjs_In = getInput()
    if rhObjs_In is None: return
    
    fTolerance = Opts.values['fTolerance']
    bMarkBadGeom = Opts.values['bMarkBadGeom']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt(prompt="Working ...")

    gBs_Res = []

    for rhObj in rhObjs_In:
        gBrep = rs.coerceguid(rhObj)

        rc = processBrepObject(
            gBrep,
            fTolerance=fTolerance,
            bMarkBadGeom=bMarkBadGeom,
            bEcho=bEcho,
            bDebug=bDebug)

        if rc is None:
            print("Error in new brep geometry.")
        elif not rc:
            print("No error in brep geometry, but no edges were joined.")
        else:
            print("Edges were joined.")
            gBs_Res.append(rc)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()