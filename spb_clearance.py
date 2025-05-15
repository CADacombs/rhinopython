"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
160718-22: Created.
...
230927: Added option to disable reflection loop.  Previously, reflection count had to be set to 1.
231013-14: Refactored.  Relabeled command options.
        Replaced V5-era library functions to obtain geometry within blocks with RhinoCommon 6+ methods.
        Reversed the order of using picked objects in main routine to be more intuitive.
231030: Added a option to add a point at intersection for trimming curves.
250427: Modified some option default values.

TODO: Fix functionality of post pick point for non-Brep/BrepFaces.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bFaceForA'; keys.append(key)
    values[key] = False
    names[key] = 'FaceFilter'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPostPickPtOnA'; keys.append(key)
    values[key] = False
    names[key] = 'StartPt'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'ReusePicked', 'PickNew')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bFaceForB'; keys.append(key)
    values[key] = False
    names[key] = 'FaceFilter'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bPostPickPtOnB'; keys.append(key)
    values[key] = False
    names[key] = 'StartPt'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'ReusePicked', 'PickNew')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bLoopToMin'; keys.append(key)
    values[key] = True
    names[key] = 'LoopToLocalMin'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iIters_Max'; keys.append(key)
    values[key] = 1000
    names[key] = 'MaxLoopIters'
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fMinPtDelta'; keys.append(key)
    values[key] = max((2.0*sc.doc.ModelAbsoluteTolerance, 1e-6))
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAddMarks'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddLine'; keys.append(key)
    values[key] = True
    names[key] = 'Line'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddPt'; keys.append(key)
    values[key] = False
    names[key] = 'Pt'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = True
    names[key] = 'Dot'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotDecPlaces'; keys.append(key)
    values[key] = 3 if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches else sc.doc.DistanceDisplayPrecision - 2
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'iDotFontHt'; keys.append(key)
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
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

        if key == 'iIters_Max':
            if cls.riOpts[key].CurrentValue <= 0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def areLayerAndAllAncestorsVisible(idx_rdLayer_Child):
    
    rdLayer_Child = sc.doc.Layers[idx_rdLayer_Child]
    if not rdLayer_Child.IsVisible: return False
    
    idLayer_Parent = rdLayer_Child.ParentLayerId
    if idLayer_Parent == idLayer_Parent.Empty: return True
    else:
        idx_rdLayer_Parent = sc.doc.Layers.Find(idLayer_Parent, True)
        return areLayerAndAllAncestorsVisible(idx_rdLayer_Parent)


def _addTextDot(pt, text='', iDotFontHt=11):
    rgDot = rg.TextDot(text, pt)
    rgDot.FontHeight = iDotFontHt
    sc.doc.Objects.AddTextDot(rgDot)


def _xformedRgObjOfBlock_WireframePick(rdInstRef, pt, xform=rg.Transform.Identity):
    
    xform *= rdInstRef.InstanceXform
    rdInstDef = rdInstRef.InstanceDefinition
    rdObjs = rdInstDef.GetObjects()
    rgObj_BreporEdge = None
    
    for rdObj in rdObjs:
        if not areLayerAndAllAncestorsVisible(rdObj.Attributes.LayerIndex):
            continue
        if rdObj.ObjectType == rd.ObjectType.InstanceReference:
            rgObj = _xformedRgObjOfBlock_WireframePick(rdObj, pt, xform)
            if rgObj is not None: return rgObj
        else: # Object is not a block instance
            rgObj = rdObj.Geometry
            rgObj.Transform(xform)
            
            # Points and curves have precedence over breps and extrusions.
            if rgObj.ObjectType == rd.ObjectType.Point:
                ptOn = rgObj.Location
                if pt.EpsilonEquals(ptOn, .01*sc.doc.ModelAbsoluteTolerance):
                    return rgObj
                continue
            if rgObj.ObjectType == rd.ObjectType.Curve:
                bPt, t = rgObj.ClosestPoint(pt)
                ptOn = rgObj.PointAt(t)
                if pt.EpsilonEquals(ptOn, .01*sc.doc.ModelAbsoluteTolerance):
                    return rgObj
                continue
            
            # Not point or curve
            if rgObj.ObjectType == rd.ObjectType.Brep:
                ptOn = rgObj.ClosestPoint(pt)
            elif rgObj.ObjectType == rd.ObjectType.Extrusion:
                bPt, u, v = rgObj.ClosestPoint(pt)
                ptOn = rgObj.PointAt(u, v)
            else: # Unsupported object
                rgObj.Dispose()
                continue
            
            # Check for point on brep or extrusion.
            if pt.EpsilonEquals(ptOn, .001*sc.doc.ModelAbsoluteTolerance):
                rgObj_BreporEdge = rgObj
    
    if rgObj_BreporEdge is not None: return rgObj_BreporEdge


def getInput():

    # Set object type filters for face vs. other object selection.
    #objTypeForFaceFilter = (rd.ObjectType.InstanceReference |
    #                        rd.ObjectType.Surface)
    objTypeForFaceFilter = rd.ObjectType.Surface
    objTypeForNotFaceFilter = ( rd.ObjectType.InstanceReference |
                                rd.ObjectType.Point |
                                rd.ObjectType.Curve |
                                rd.ObjectType.Brep |
                                rd.ObjectType.Extrusion)
    objTypeForNotFaceFilter = ( rd.ObjectType.Point |
                                rd.ObjectType.Curve |
                                rd.ObjectType.Brep |
                                rd.ObjectType.Extrusion)


    def getObj(obj_A, bAcceptPreselection=False):
        # Get object with optional input.

        go = ri.Custom.GetObject()

        go.AcceptNumber(True, acceptZero=True)
        go.EnableHighlight(False)

        idxs_Opt = {} # To be reused.
        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        while True:
            go.ClearCommandOptions()

            idxs_Opt.clear()

            if obj_A is None:
                addOption('bFaceForA')
                addOption('bPostPickPtOnA')
                if Opts.values['bFaceForA']:
                    go.SetCommandPrompt("Select face of 1st object")
                    go.GeometryFilter = objTypeForFaceFilter
                else:
                    go.SetCommandPrompt("Select 1st object")
                    go.GeometryFilter = objTypeForNotFaceFilter
            else:
                addOption('bFaceForB')
                addOption('bPostPickPtOnB')
                if Opts.values['bFaceForB']:
                    go.SetCommandPrompt("Select face of 2nd object")
                    go.GeometryFilter = objTypeForFaceFilter
                else:
                    go.SetCommandPrompt("Select 2nd object")
                    go.GeometryFilter = objTypeForNotFaceFilter

            addOption('bLoopToMin')
            if Opts.values['bLoopToMin']:
                addOption('iIters_Max')
                addOption('fMinPtDelta')
            addOption('bAddMarks')
            if Opts.values['bAddMarks']:
                addOption('bAddLine')
                addOption('bAddPt')
                addOption('bAddDot')
                if Opts.values['bAddDot']:
                    addOption('iDotDecPlaces')
                    addOption('iDotFontHt')
            addOption('bEcho')
            addOption('bDebug')

            res = go.Get()

            if res == ri.GetResult.Cancel:
                go.Dispose()
                return

            if res == ri.GetResult.Object:
                objref = go.Object(0)
                if Opts.values['bDebug']:
                    sEval = "objref.Object()"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "objref.Geometry()"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "objref.InstanceDefinitionPart()"; print("{}: {}".format(sEval, eval(sEval)))
                    # Avoid using objref.Brep(), because once it is, the result of objref.Face() will become None.
                    rgF = objref.Face()
                    if rgF is None:
                        sEval = "objref.Brep()"; print("{}: {}".format(sEval, eval(sEval)))
                    else:
                        sEval = "rgF.Brep.Faces.Count"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "objref.Face()"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "objref.Surface()"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "objref.Curve()"; print("{}: {}".format(sEval, eval(sEval)))
                    sEval = "objref.Point()"; print("{}: {}".format(sEval, eval(sEval)))
                go.Dispose()
                return objref

                # Prepare for repeat of go.Get().
                #sc.doc.Objects.UnselectAll()
                #sc.doc.Views.Redraw()
                #continue

            # An option was selected or a number was entered.
            if res == ri.GetResult.Number:
                if not Opts.values['bLoopToMin']:
                    continue
                key = 'iIters_Max'
                Opts.riOpts[key].CurrentValue = int(go.Number())
                Opts.setValue(key)
                continue
                #optI_iIters_Max.CurrentValue = go.Number()

            for key in idxs_Opt:
                if go.Option().Index == idxs_Opt[key]:
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)
                    break

    objref_A = getObj(obj_A=None, bAcceptPreselection=True)
    if not objref_A: return
    rgObjA = objref_A.Geometry()
    if isinstance(rgObjA, rg.BrepFace):
        rgObjA = rgObjA.DuplicateFace(duplicateMeshes=False)

    if not Opts.values['bPostPickPtOnA']:
        sc.doc.Views.Redraw()
        ptA = objref_A.SelectionPoint()
    else:
        #if not isinstance(objref_A.Geometry(), rg.BrepFace):
        #    print(
        gp = ri.Custom.GetPoint()
        gp.SetCommandPrompt("Select analysis start point on 1st object")
        gp.Constrain(brep=rgObjA, wireDensity=-1, faceIndex=-1, allowPickingPointOffObject=False)
        gp.Get()
        if gp.CommandResult() != Rhino.Commands.Result.Success:
            return
        else:
            ptA = gp.Point()

    sc.doc.Objects.UnselectAll(); sc.doc.Views.Redraw()


    objref_B = getObj(obj_A=rgObjA, bAcceptPreselection=False)
    if not objref_B: return
    rgObjB = objref_B.Geometry()
    if isinstance(rgObjB, rg.BrepFace):
        rgObjB = rgObjB.DuplicateFace(duplicateMeshes=False)

    if not Opts.values['bPostPickPtOnB']:
        sc.doc.Views.Redraw()
        ptB = objref_B.SelectionPoint()
    else:
        gp = ri.Custom.GetPoint()
        gp.SetCommandPrompt("Select analysis start point on 2nd object")
        gp.Constrain(brep=rgObjB, wireDensity=-1, faceIndex=-1, allowPickingPointOffObject=False)
        gp.Get()
        if gp.CommandResult() != Rhino.Commands.Result.Success:
            return
        else:
            ptB = gp.Point()

    return (
        rgObjA,
        ptA,
        rgObjB,
        ptB,
        )


def _closestPoint(rgObj, testpoint):
    if isinstance(rgObj, rg.Point):
        return rgObj.Location

    if isinstance(rgObj, rg.Curve):
        b, t = rgObj.ClosestPoint(testpoint)
        if not b:
            raise Exception("Closest point not found for curve.")
        return rgObj.PointAt(t)

    if isinstance(rgObj, rg.Brep):
        return rgObj.ClosestPoint(testpoint)

    raise Exception("{} passed to _closestPoint.".format(rgObj.GetType().Name))


def main():
    
    rc = getInput()
    if rc is None: return
    (
        rgObjA,
        ptA,
        rgObjB,
        ptB,
    ) = rc


    bLoopToMin = Opts.values['bLoopToMin']
    iIters_Max = Opts.values['iIters_Max']
    fMinPtDelta = Opts.values['fMinPtDelta']
    bAddMarks = Opts.values['bAddMarks']
    bAddLine = Opts.values['bAddLine']
    bAddPt = Opts.values['bAddPt']
    bAddDot = Opts.values['bAddDot']
    iDotDecPlaces = Opts.values['iDotDecPlaces']
    iDotFontHt = Opts.values['iDotFontHt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    ptB_Prev = ptB
    ptB = _closestPoint(rgObjB, ptA)

    if bLoopToMin:
        # Loop until difference between old and new B points are within a specific distance.

        if (    rgObjA.ObjectType == rd.ObjectType.Brep or 
                rgObjB.ObjectType == rd.ObjectType.Brep):
            tol = fMinPtDelta
        else:
            tol = Rhino.RhinoMath.ZeroTolerance

        if bDebug: sEval = "tol"; print("{}: {}".format(sEval, eval(sEval)))

        distB_Delta = ptB.DistanceTo(ptB_Prev)

        iIters = 1

        while True:
            sc.escape_test()

            ptA_Prev = ptA
            ptA = _closestPoint(rgObjA, ptB)
            if bDebug: sc.doc.Objects.AddPoint(ptA)

            iIters += 1

            if iIters == iIters_Max:
                break

            distA_Delta = ptA.DistanceTo(ptA_Prev)
            if distA_Delta < tol:
                if distB_Delta < tol:
                    break

            ptB_Prev = ptB
            ptB = _closestPoint(rgObjB, ptA)
            if bDebug: sc.doc.Objects.AddPoint(ptB)

            iIters += 1

            if iIters == iIters_Max:
                break

            distB_Delta = ptB.DistanceTo(ptB_Prev)
            if distB_Delta < tol:
                if distA_Delta < tol:
                    break

    dist = ptB.DistanceTo(ptA)
    if dist <= sc.doc.ModelAbsoluteTolerance:
        sPrint = ["Objects intersect."]
        
        if bAddMarks:
            if bAddDot:
                _addTextDot(
                    (ptA+ptB)/2.0,
                    "X",
                    iDotFontHt=iDotFontHt)
            if bAddPt:
                gPt = sc.doc.Objects.AddPoint((ptA+ptB)/2.0)
                if gPt != gPt.Empty:
                    if sc.doc.Objects.Select(gPt):
                        sPrint.append("Point at intersection was created and is selected.")
        
        print("  ".join(sPrint))
        
        sc.doc.Views.Redraw()
        return

    s = "Clearance: {:.{}f} {}".format(
        dist,
        sc.doc.ModelDistanceDisplayPrecision,
        str(sc.doc.ModelUnitSystem).lower())

    if sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
        s += " [{:.{}f} millimeters]".format(
            dist*25.4,
            sc.doc.ModelDistanceDisplayPrecision-2)
    elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
        s += " [{:.{}f} inches]".format(
            dist/25.4,
            sc.doc.ModelDistanceDisplayPrecision+2)

    if bLoopToMin:
        s += " found in {} iterations.".format(iIters)

    print(s)

    if bAddMarks:
        if bAddLine:
            sc.doc.Objects.AddLine(ptA, ptB)
        if bAddDot:
            _addTextDot(
                (ptA+ptB)/2.0,
                '{0:.{1}f}'.format(dist, iDotDecPlaces),
                iDotFontHt=iDotFontHt)

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()