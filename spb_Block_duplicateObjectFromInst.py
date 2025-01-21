"""
This script duplicates object of block instance by a single pick.
Non-uniform blocks are supported.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160416-27: Created.
...
170208: ArcCurve, PolyCurve, and Brep geometries now are made deformable when the
        transform's SimilarityType is NotSimilarity, e.g., non-uniform scale.
180309: In response to V6 allowing objects to be added to reference layers,
        option to retain the duplicated object's layer was added.
180926: Now, geometry from linked (not embedded) blocks will only be created on current layer.
...
230718: Now, block instances from references are placed on current layer.
230817: Modified some command option text.
231009: Bug fix for when BlockInstance is the target.
240629: Annotations are now skipped during search of lowest level picked geometry.
250120: Bug fix in command options.

TODO:
    Point and TextDots require very accurate picking.
            Possibly allow a distance tolerance for these objects.
    Support adding BlockInstance of reference object.

Coding notes:
    ArcCurve is not supported by MakeDeformable.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


sObjTypes = (
        'bBlockInstance',
        'bBrep',
        'bCurve',
        'bExtrusion',
        'bPoint',
        'bTextDot',
)


sOpts = sObjTypes + (
        'bFace',
        'bRetainLayer_NotCurrent',
        'bEcho',
        'bDebug',
)


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bBlockInstance'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bBrep'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bCurve'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bExtrusion'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bPoint'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bTextDot'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bFace'; keys.append(key)
    values[key] = True
    names[key] = 'BrepComponent'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'EntireBrep', 'Face')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bRetainLayer_NotCurrent'; keys.append(key)
    values[key] = False
    names[key] = 'OutputLayer'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Current', 'TryToRetain') # TryToRetain because current layer is used when an object from a reference block is picked.
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

        #if key == 'fRadius':
        #    if cls.riOpts[key].CurrentValue < 2.0*sc.doc.ModelAbsoluteTolerance:
        #        cls.riOpts[key].CurrentValue = 0.0

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def findPickedObjOfBlockInst(rdInstObj, objTypes_ToAccept, pt_Picked, xform=Rhino.Geometry.Transform.Identity):
    xform *= rdInstObj.InstanceXform
    rdInstDef = rdInstObj.InstanceDefinition
    rdObjs = rdInstDef.GetObjects()
    xformSimilarityType = xform.SimilarityType

    for rdObj in rdObjs:

        # No ClosestPoint for Annotations so skip them.
        if rdObj.ObjectType == rd.ObjectType.Annotation:
            continue

        if rdObj.ObjectType == rd.ObjectType.InstanceReference:
            
            # Dig into block instances until other object type is found.
            
            rc = findPickedObjOfBlockInst(rdObj, objTypes_ToAccept, pt_Picked, xform)
            if rc:
                return rc
            continue # on to next rdObj.

        # Object is not a block instance.
        
        rgObj = rdObj.DuplicateGeometry()
        if not rgObj.IsValid:
            s  = "{} has invalid geometry.".format(rdObj.Id)
            s += "  Nothing will be duplicated."
            s += "  Repair invalid objects in block definitions."
            print(s)
            return
        
        if xformSimilarityType == rg.TransformSimilarityType.NotSimilarity:
            # Convert ArcCurves, PolyCurves (in case the contain arc segments), &
            # Breps (in case they contain faces of non-NURBS surfaces) into
            # deformable geometry.
            
            typeGeom = rgObj.GetType()

            # ArcCurve is not supported by MakeDeformable.
            if typeGeom == rg.ArcCurve:
                #                    # This doesn't work.
                #                    print(typeGeom)
                #                    if not rgObj.MakeDeformable():
                #                        print("Error in converting object into NURBS representation.  Object will not be duplicated.")
                #                        continue
                
                rgObjNurbsCrv = rgObj.ToNurbsCurve()
                if rgObjNurbsCrv is not None: rgObj = rgObjNurbsCrv
                else:
                    print("Error in converting arc curve into NURBS curve!")
                    continue
            elif typeGeom == rg.PolyCurve or typeGeom == rg.Brep:
                if not rgObj.MakeDeformable():
                    print("Error in making {} deformable!".format(typeGeom))
                    continue
                if not rgObj.IsValid:
                    print("{} is not valid after MakeDeformable!  Not using MakeDeformable...".format(rgObj.ObjectType))
                    b, sLog = rgObj.IsValidWithLog()
                    print(sLog)
                    rgObj = rdObj.DuplicateGeometry()
                    #sc.doc.Objects.Add(rgObj) # For debugging.  Add object without transformation.
        
        rgObj.Transform(xform)
        
        if not rgObj.IsValid:
            print("{} is not valid after transformation!  Skipping object...".format(rgObj.ObjectType))
            rgObj.Dispose()
            continue
        rdObjType_ThisObj = rgObj.ObjectType
        if objTypes_ToAccept & rd.ObjectType.InstanceReference:
            if rdObjType_ThisObj == rd.ObjectType.Brep:
                pt_ClosestOnObj = rgObj.ClosestPoint(pt_Picked)
            elif rdObjType_ThisObj == rd.ObjectType.Curve:
                bPt, u = rgObj.ClosestPoint(pt_Picked)
                pt_ClosestOnObj = rgObj.PointAt(u)
            elif rdObjType_ThisObj == rd.ObjectType.Extrusion:
                bPt, u, v = rgObj.ClosestPoint(pt_Picked)
                pt_ClosestOnObj = rgObj.PointAt(u, v)
                pt_ClosestOnObj = rgObj.ClosestPoint(pt_Picked)
            elif rdObjType_ThisObj == rd.ObjectType.Point:
                pt_ClosestOnObj = rgObj.Location
            elif rdObjType_ThisObj == rd.ObjectType.TextDot:
                pt_ClosestOnObj = rgObj.Point
            else:
                rgObj.Dispose()
                continue

            rgObj.Dispose()
            
            if pt_Picked.DistanceTo(pt_ClosestOnObj) <= sc.doc.ModelAbsoluteTolerance:
                return rdInstObj, xform
        else:
            # Target is not InstanceReference.
            if objTypes_ToAccept & rd.ObjectType.Brep & rdObjType_ThisObj:
                pt_ClosestOnObj = rgObj.ClosestPoint(pt_Picked)
            elif objTypes_ToAccept & rd.ObjectType.Curve & rdObjType_ThisObj:
                bPt, u = rgObj.ClosestPoint(pt_Picked)
                pt_ClosestOnObj = rgObj.PointAt(u)
            elif objTypes_ToAccept & rd.ObjectType.Extrusion & rdObjType_ThisObj:
                bPt, u, v = rgObj.ClosestPoint(pt_Picked)
                pt_ClosestOnObj = rgObj.PointAt(u, v)
            elif objTypes_ToAccept & rd.ObjectType.Point & rdObjType_ThisObj:
                pt_ClosestOnObj = rgObj.Location
            elif objTypes_ToAccept & rd.ObjectType.TextDot & rdObjType_ThisObj:
                pt_ClosestOnObj = rgObj.Point
            else:
                rgObj.Dispose()
                continue
            
            rgObj.Dispose()
            
            if pt_Picked.DistanceTo(pt_ClosestOnObj) <= sc.doc.ModelAbsoluteTolerance:
                return rdObj, xform


def getInput():
    """
    Get block instance with optional input
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick object of block instance")

    go.GeometryFilter = rd.ObjectType.InstanceReference

    go.DisablePreSelect()

    go.OneByOnePostSelect = True

    idxs_Opts = {}

    def addOption(ric, key): idxs_Opts[key] = Opts.addOption(ric, key)


    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption(go, 'bBlockInstance')
        if not Opts.values['bBlockInstance']:
            idxs_Opts['CrvsOnly'] = go.AddOption('CrvsOnly')
            idxs_Opts['BrepsOnly'] = go.AddOption('BrepsOnly')
            idxs_Opts['ObjectFilter'] = go.AddOption('ObjectFilter')
            if Opts.values['bBrep']:
                addOption(go, 'bFace')
        addOption(go, 'bRetainLayer_NotCurrent')
        addOption(go, 'bEcho')
        addOption(go, 'bDebug')

        sTrue = []
        sFalse = []
        for key in sObjTypes:
            if Opts.values[key]:
                sTrue.append(Opts.names[key])
            else:
                sFalse.append(Opts.names[key])

        if not Opts.values['bBlockInstance']:
            s  = "Objects that will ONLY be accepted:"
            s += " {}".format(", ".join(sTrue) if sTrue else 'None')
            s += "\t\tNOT accepted:"
            s += " {}".format(", ".join(sFalse) if sFalse else 'None')
            print(s)

        objTypes_ToAccept = rd.ObjectType.None
        if Opts.values['bBlockInstance']:
            objTypes_ToAccept |= rd.ObjectType.InstanceReference
            objTypes_ToAccept |= rd.ObjectType.Brep
            objTypes_ToAccept |= rd.ObjectType.Curve
            objTypes_ToAccept |= rd.ObjectType.Extrusion
            objTypes_ToAccept |= rd.ObjectType.Point
            objTypes_ToAccept |= rd.ObjectType.TextDot
        else:
            if Opts.values['bBrep']: objTypes_ToAccept |= rd.ObjectType.Brep
            if Opts.values['bCurve']: objTypes_ToAccept |= rd.ObjectType.Curve
            if Opts.values['bExtrusion']: objTypes_ToAccept |= rd.ObjectType.Extrusion
            if Opts.values['bPoint']: objTypes_ToAccept |= rd.ObjectType.Point
            if Opts.values['bTextDot']: objTypes_ToAccept |= rd.ObjectType.TextDot

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            
            rdInstObj = objref.Object()
            
            rc = findPickedObjOfBlockInst(rdInstObj, objTypes_ToAccept, objref.SelectionPoint())
            if not rc:
                print("Picked object not found.  Pick again.")
                continue

            if objTypes_ToAccept & rd.ObjectType.InstanceReference:
                if rc[0].InstanceDefinition.Index == rdInstObj.InstanceDefinition.Index:
                    print("Nested block instance not found.")
                    go.Dispose()
                    return

            # Success.

            rdObj, xform = rc

            go.Dispose()

            return (
                rdObj,
                xform,
                Opts.values['bRetainLayer_NotCurrent'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        # An option was selected.

        if 'CrvsOnly' in idxs_Opts and go.Option().Index == idxs_Opts['CrvsOnly']:
            for key in sObjTypes:
                if key == 'bCurve':
                    Opts.riOpts[key].CurrentValue = True
                else:
                    Opts.riOpts[key].CurrentValue = False
                Opts.setValue(key)
            continue

        if 'BrepsOnly' in idxs_Opts and go.Option().Index == idxs_Opts['BrepsOnly']:
            for key in sObjTypes:
                if key == 'bBrep':
                    Opts.riOpts[key].CurrentValue = True
                else:
                    Opts.riOpts[key].CurrentValue = False
                Opts.setValue(key)
            continue

        if 'ObjectFilter' in idxs_Opts and go.Option().Index == idxs_Opts['ObjectFilter']:

            go_ObjType = ri.Custom.GetOption()
            go_ObjType.SetCommandPrompt("Object type(s) to filter")

            while True:
                idxs_Opts.clear()

                addOption(go_ObjType, 'bBlockInstance')
                addOption(go_ObjType, 'bBrep')
                addOption(go_ObjType, 'bCurve')
                addOption(go_ObjType, 'bExtrusion')
                addOption(go_ObjType, 'bPoint')
                addOption(go_ObjType, 'bTextDot')
                idxs_Opts['YesToAll'] = go_ObjType.AddOption('YesToAll')
                idxs_Opts['NoToAll'] = go_ObjType.AddOption('NoToAll')

                res = go_ObjType.Get()
                if res != ri.GetResult.Option:
                    break
                    
                if go_ObjType.OptionIndex() == idxs_Opts['YesToAll']:
                    for key in sObjTypes:
                        Opts.riOpts[key].CurrentValue = True
                elif go_ObjType.OptionIndex() == idxs_Opts['NoToAll']:
                    for key in sObjTypes:
                        Opts.riOpts[key].CurrentValue = False

                for key in idxs_Opts:
                    if go_ObjType.Option().Index == idxs_Opts[key]:
                        Opts.setValue(key, go_ObjType.Option().CurrentListOptionIndex)
                        break

                go_ObjType.ClearCommandOptions()

            go_ObjType.Dispose()
            continue

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break

    go.Dispose()

    return rdObj, rgObj, bRetainLayer_NotCurrent


def processInput(rdObj, xform, bRetainLayer_NotCurrent=True, bEcho=True, bDebug=False):
    """
    """
    
    attr = rdObj.Attributes.Duplicate()
    
    if Opts.values['bRetainLayer_NotCurrent']:
        if rdObj.IsReference:
            s  = "{} is not embedded.".format(rdObj.ObjectType)
            s += "  Geometry will be created on the current layer."
            print(s)
            attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    else:
        attr.LayerIndex = sc.doc.Layers.CurrentLayerIndex
    
    if rdObj.ObjectType == rd.ObjectType.InstanceReference:
        idxInstDef = rdObj.InstanceDefinition.Index
        if rdObj.IsReference:
            # The following doesn't work in a script because per
            # https://discourse.mcneel.com/t/history-in-script/171964/20
            # "Because history is command based, you might consider writing a plug-in" ...
            # gObj1 = sc.doc.Objects.AddInstanceObject(
            #     instanceDefinitionIndex=idxInstDef,
            #     instanceXform=xform,
            #     attributes=attr,
            #     history=rd.HistoryRecord(
            #     reference=False)

            gObj1 = sc.doc.Objects.AddInstanceObject(idxInstDef, xform, attr)

            # Attempt to embed a reference block.
            #            gObj_eraseme = sc.doc.Objects.AddInstanceObject(idxInstDef, xform)
            #            sc.doc.Objects.Delete(objectId=gObj_eraseme, quiet=True)
            rdObj1 = sc.doc.Objects.FindId(gObj1) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(gObj1)
            # Attempt to embed a reference block.
            #            rdObj1.Attributes = rd.ObjectAttributes()
            #            print(rdObj1.CommitChanges())
            if gObj1 == gObj1.Empty:
                print("Error!  {} not added.".format(rdObj.ObjectType))
                return
            elif sc.doc.Objects.Select(gObj1):
                s  = "Duplicated {} is selected.".format(rdObj.ObjectType)
                s += "  It is a reference, so modify a duplicate of it."
                print(s)
                return gObj1
        else:
            gObj1 = sc.doc.Objects.AddInstanceObject(idxInstDef, xform, attr)
    else:
        xformSimilarityType = xform.SimilarityType
        rgObj_Xformed = rdObj.Geometry.Duplicate()
        if xformSimilarityType == rg.TransformSimilarityType.NotSimilarity:
            # ArcCurve is not supported by MakeDeformable.
            if isinstance(rgObj_Xformed, rg.ArcCurve):
                rgObj_Xformed = rgObj_Xformed.ToNurbsCurve()
                if rgObj_Xformed is None:
                    print("Error in converting arc curve into NURBS curve!")
                    return

            if isinstance(rgObj_Xformed, (rg.PolyCurve, rg.Brep)):
                if not rgObj_Xformed.MakeDeformable():
                    print("Error in making {} deformable!".format(rgObj_Xformed.ObjectType))
                    return
                if not rgObj_Xformed.IsValid:
                    print("{} is not valid after MakeDeformable!  Not using MakeDeformable...".format(rgObj_Xformed.ObjectType))
                    b, sLog = rgObj_Xformed.IsValidWithLog()
                    print(sLog)
                    rgObj_Xformed.Dispose()
                    rgObj_Xformed = rdObj.Geometry.Duplicate()

        rgObj_Xformed.Transform(xform)
        gObj1 = sc.doc.Objects.Add(rgObj_Xformed, attr)
        rgObj_Xformed.Dispose()
    
    if gObj1 == gObj1.Empty:
        print("Error!  {} not added.".format(rdObj.ObjectType))
        return
    elif sc.doc.Objects.Select(gObj1):
        print("Duplicated {} is selected.".format(rdObj.ObjectType))
    
    return gObj1


def main():

    rc = getInput()
    if rc is None: return
    (
        rdObj,
        xform,
        bRetainLayer_NotCurrent,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    processInput(
        rdObj=rdObj,
        xform=xform,
        bRetainLayer_NotCurrent=bRetainLayer_NotCurrent,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()