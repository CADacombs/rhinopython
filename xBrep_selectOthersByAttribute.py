"""
201030: Created, starting with xBreps_similar.py.
201102: Added support for Extrusions.  Now considers Attributes.ColorSource.  Added feedback for when analyzing large data sets.
201105: Improved efficency of volume checking and open breps are now skipped for volume matching.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


sAttrs = [
    'Colour',
    'Name',
    'Layer',
    'Volume',
    ]


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}
    
    key = 'fFloatTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=sc.doc.ModelAbsoluteTolerance)
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


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
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """
    
    # Get Brep/Extrusion with optional input.
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select source brep/extrusion")
    go.GeometryFilter = rd.ObjectType.Brep # Brep will also filter Extrusions ( https://discourse.mcneel.com/t/restriction-of-objecttype-in-rhinocommon/73603/5 )
    go.SubObjectSelect = False

    idxs_Opt = {}

    while True:
        key = 'fFloatTol'; idxs_Opt[key] = Opts.addOption(go, key)
        key = 'bEcho'; idxs_Opt[key] = Opts.addOption(go, key)
        key = 'bDebug'; idxs_Opt[key] = Opts.addOption(go, key)

        res = go.Get()

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            go = ri.Custom.GetOption()
            go.SetCommandPrompt("Select matching attribute")


            idxs_Opt = {}

            while True:
                for key in sAttrs:
                    idxs_Opt[key] = go.AddOption(key)
                key = 'fFloatTol'; idxs_Opt[key] = Opts.addOption(go, key)
                key = 'bEcho'; idxs_Opt[key] = Opts.addOption(go, key)
                key = 'bDebug'; idxs_Opt[key] = Opts.addOption(go, key)

                res = go.Get()

                if res == ri.GetResult.Object:
                    objref = go.Object(0)
                    go.Dispose()
                    return tuple([objref] + [Opts.values[key] for key in Opts.keys])
                elif res == ri.GetResult.Cancel:
                    go.Dispose()
                    return
                else:
                    # An option was selected or a number was entered.
                    for key in sAttrs:
                        if go.OptionIndex() == idxs_Opt[key]:
                            return tuple([objref] + [key] + [Opts.values[key] for key in Opts.keys])

                Opts.setValues()
                Opts.saveSticky()
                go.ClearCommandOptions()

        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return
        else:
            # An option was selected.
            pass

        Opts.saveSticky()


def getMatches(objref_toMatch, **kwargs):
    """
    """


    def setOpt(key, value=None):
        if key in kwargs:
            return kwargs[key]
        elif key in Opts.riOpts:
            return Opts.riOpts[key].InitialValue
        else:
            return value


    sAttr = setOpt('sAttr')
    fFloatTol = setOpt('fFloatTol')
    bEcho = setOpt('bEcho')
    bDebug = setOpt('bDebug')


    objref0 = objref_toMatch

    bToMatch_IsExtrusion = bool(objref0.Surface())

    if sAttr == "Colour":
        rdRhObj0 = objref0.Object()

        if rdRhObj0.Attributes.ColorSource == rd.ObjectColorSource.ColorFromMaterial:
            print "Colour from material not supported yet."
            return
        elif rdRhObj0.Attributes.ColorSource == rd.ObjectColorSource.ColorFromObject:
            colorRhObj0 = rdRhObj0.Attributes.ObjectColor
        else:
            if rdRhObj0.Attributes.ColorSource == rd.ObjectColorSource.ColorFromLayer:
                pass
            else:
                # Color from parent.
                # Use layer color unless object is within a block definition.
                pass

            li = rdRhObj0.Attributes.LayerIndex
            layer = sc.doc.Layers.FindIndex(li)
            colorRhObj0 = layer.Color


    elif sAttr == 'Volume':
        rgBrep0 = objref0.Brep() # Brep is also returned when objref0 contains an Extrusion.

        if not rgBrep0.IsValid:
            print "Reference {} {} is invalid.  Repair it then rerun this script.".format(
                "Extrusion" if bToMatch_IsExtrusion else "Brep",
                objref0.ObjectId)
            rgBrep0.Dispose()
            return

        if not rgBrep0.IsSolid:
            print "Reference {} {} is open.  Its 'Volume' will not be matched.".format(
                "Extrusion" if bToMatch_IsExtrusion else "Brep",
                objref0.ObjectId)
            rgBrep0.Dispose()
            return

        fVol0 = rgBrep0.GetVolume()
        if not fVol0:
            print "Volume can not be calculated.  Repair {} {}.".format(
                "Extrusion" if bToMatch_IsExtrusion else "Brep",
                objref0.ObjectId)
            rgBrep0.Dispose()
            return

        iToMatch_FaceCt = rgBrep0.Faces.Count

        rgBrep0.Dispose()
    
    gBreps_ToSearch = []
    rgBreps_ToSearch = []
    iCts_Faces_ToSearch = []
    fEdgeLens_ToSearch = []

    iter = rd.ObjectEnumeratorSettings()
    iter.NormalObjects = True
    iter.LockedObjects = False
    iter.IncludeLights = False
    iter.IncludeGrips = False
    iter.ObjectTypeFilter = rd.ObjectType.Brep | rd.ObjectType.Extrusion

    gBreps_MatchesFound = [] # Also can include Extrusions.

    iBrepExtrCt = 0
    for rdObj in sc.doc.Objects.GetObjectList(iter):
        iBrepExtrCt += 1

    iOpenBrepCt = 0

    idxs_AtTenths = [int(round(0.1*i*iBrepExtrCt,0)) for i in range(10)]

    for i, rdRhObjX in enumerate(sc.doc.Objects.GetObjectList(iter)):
        if sc.escape_test(throw_exception=False):
            print "*** Analysis interrupted by user." \
                "  Selected breps/extrusions are of partial results."
            return gBreps_MatchesFound

        if rdRhObjX.Id == objref0.ObjectId: continue

        if iBrepExtrCt > 10:
            if i in idxs_AtTenths:
                Rhino.RhinoApp.SetCommandPrompt("Analysis at {:d}% of {} breps/extrusions ...".format(
                    int(100.0 * (i+1) / iBrepExtrCt), iBrepExtrCt))
        elif iBrepExtrCt > 1:
            Rhino.RhinoApp.SetCommandPrompt(
                "Analysis at {} of {} breps/extrusions".format(
                    i+1, iBrepExtrCt))
        else:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing other brep/extrusion ...")


        if sAttr == 'Colour':
            if rdRhObjX.Attributes.ColorSource == rd.ObjectColorSource.ColorFromObject:
                colorRhObjX = rdRhObjX.Attributes.ObjectColor
            elif rdRhObjX.Attributes.ColorSource == rd.ObjectColorSource.ColorFromMaterial:
                print "Colour from material not supported yet."
                return
            else:
                if rdRhObjX.Attributes.ColorSource == rd.ObjectColorSource.ColorFromLayer:
                    pass
                else:
                    # Color from parent.
                    # Use layer color unless object is within a block definition.
                    pass

                li = rdRhObjX.Attributes.LayerIndex
                layer = sc.doc.Layers.FindIndex(li)
                colorRhObjX = layer.Color

            if colorRhObjX == colorRhObj0:
                gBreps_MatchesFound.append(rdRhObjX.Id)
        elif sAttr == 'Name':
            if rdRhObjX.Attributes.Name == objref0.Object().Attributes.Name:
                gBreps_MatchesFound.append(rdRhObjX.Id)
        elif sAttr == 'Layer':
            if rdRhObjX.Attributes.LayerIndex == objref0.Object().Attributes.LayerIndex:
                gBreps_MatchesFound.append(rdRhObjX.Id)
        elif sAttr == 'Volume':

            bToCheck_IsExtrusion = isinstance(rdRhObjX, rd.ExtrusionObject)

            rgGeomX = rdRhObjX.Geometry

            if not rgGeomX.IsValid:
                print "{} {} is invalid.  Fix first.".format(
                    rdRhObjX.GetType().Name,
                    rdRhObjX.Id)
                rgGeomX.Dispose()
                continue

            if bToCheck_IsExtrusion:
                rgBrepX = rgGeomX.ToBrep(splitKinkyFaces=True)
                rgGeomX.Dispose()
            else:
                rgBrepX = rgGeomX

            if not rgBrepX.IsSolid:
                iOpenBrepCt += 1
                rgBrepX.Dispose()
                continue

            ## This significantly speeds up the analysis.
            #if rdRhObjX.ObjectType == rd.ObjectType.Brep:
            #    if rgBrepX.Faces.Count != iToMatch_FaceCt:
            #        rgBrepX.Dispose()
            #        continue

            fVol = rgBrepX.GetVolume() # GetVolume may be faster than VolumeMassProperties.Compute.
            if not fVol:
                print "Volume can not be calculated.  Repair {} {}.".format(
                    rdRhObjX.GetType().Name,
                    rdRhObjX.Id)
                rgBrepX.Dispose()
                continue

            rgBrepX.Dispose()

            fVolDiff = abs(fVol - fVol0)
            if bDebug: print "Volume:", fVol0, fVolDiff
            if fVolDiff > fFloatTol:
                continue

            gBreps_MatchesFound.append(rdRhObjX.Id)


    if sAttr == 'Volume':
        if iOpenBrepCt:
            print "{} open breps skipped for volume matching.".format(iOpenBrepCt)
        else:
            print "No open breps are present."


    return gBreps_MatchesFound


def main():
    """
    """
    
    rc = getInput()
    if rc is None: return
    (
            objref_toMatch,
            sAttr,
            fFloatTol,
            bEcho,
            bDebug,
    ) = rc
    
    if bDebug:
        pass
        #reload()
    else:
        sc.doc.Views.RedrawEnabled = False
    
    #print list(sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False))

    sc.doc.Objects.UnselectAll()
    
    rc = getMatches(
            objref_toMatch=objref_toMatch,
            sAttr=sAttr,
            fFloatTol=fFloatTol,
            bEcho=bEcho,
            bDebug=bDebug,
            )
    if not rc:
        print "No matches found.  No objects are selected."
        return
    
    gBreps_MatchesFound = rc
    
    sc.doc.Objects.Select(objectIds=[objref_toMatch.ObjectId] + gBreps_MatchesFound)

    rdObjs_Sel = list(sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False))

    if rdObjs_Sel:
        print "{} [{}] breps / [{}] extrusions are selected.".format(
            len(rdObjs_Sel),
            len([rdObj for rdObj in rdObjs_Sel if rdObj.ObjectType == rd.ObjectType.Brep]),
            len([rdObj for rdObj in rdObjs_Sel if rdObj.ObjectType == rd.ObjectType.Extrusion]))
    else:
        print "No objects are selected."

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()