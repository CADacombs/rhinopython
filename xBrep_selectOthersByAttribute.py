"""
201030: Created, starting with xBreps_similar.py.
"""

import Rhino
import Rhino.DocObjects as rd
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
    
    # Get brep with optional input.
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select source brep")
    go.GeometryFilter = rd.ObjectType.Brep
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

    rgBrep0 = objref0.Brep()
    if not rgBrep0.IsValid:
        print "Brep {} is invalid.  Repair it then rerun this script.".format(
            objref0.ObjectId)
        rgBrep0.Dispose()
        return
        
    if sAttr == 'Volume':
        fVol0 = rgBrep0.GetVolume()
        if not fVol0:
            print "Volume can not be calculated.  Repair brep {}.".format(objref0.ObjectId)
            rgBrep0.Dispose()
        
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
    rdBrepObjects = []
    
    Rhino.RhinoApp.CommandPrompt = "Scanning objects ..."
    
    gBreps_MatchesFound = []
    
    for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
        sc.escape_test()
        
        if rdRhinoObject.ObjectType != rd.ObjectType.Brep: continue
        if rdRhinoObject.Id == objref0.ObjectId: continue
        
        rgBrepX = rdRhinoObject.BrepGeometry
        
        if not rgBrepX.IsValid:
            print "Brep {} is invalid.  Fix first.".format(rdRhinoObject.Id)
            rgBrepX.Dispose()
            continue
        
        if sAttr == 'Colour':
            if rdRhinoObject.Attributes.ObjectColor == objref0.Object().Attributes.ObjectColor:
                gBreps_MatchesFound.append(rdRhinoObject.Id)
        if sAttr == 'Name':
            if rdRhinoObject.Attributes.Name == objref0.Object().Attributes.Name:
                gBreps_MatchesFound.append(rdRhinoObject.Id)
        if sAttr == 'Layer':
            if rdRhinoObject.Attributes.LayerIndex == objref0.Object().Attributes.LayerIndex:
                gBreps_MatchesFound.append(rdRhinoObject.Id)
        elif sAttr == 'Volume':
            fVol = rgBrepX.GetVolume()
            if not fVol:
                print "Volume can not be calculated.  Repair brep {}.".format(rdRhinoObject.Id)
                rgBrepX.Dispose()
                continue

            fVolDiff = abs(fVol - fVol0)
            if bDebug: print "Volume:", fVol0, fVolDiff
            if fVolDiff > fFloatTol:
                continue

            gBreps_MatchesFound.append(rdRhinoObject.Id)

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
    
    sc.doc.Objects.UnselectAll()
    
    nSelected_All = 0
    
    rc = getMatches(
            objref_toMatch=objref_toMatch,
            sAttr=sAttr,
            fFloatTol=fFloatTol,
            bEcho=bEcho,
            bDebug=bDebug,
            )
    if not rc: return
    
    gBreps_MatchesFound = rc
    
    sc.doc.Objects.Select(objectIds=[objref_toMatch.ObjectId] + gBreps_MatchesFound)
    
    rdObjs_Sel = sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)
    
    print "{} breps are selected.".format(len(list(rdObjs_Sel)))
    
    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()