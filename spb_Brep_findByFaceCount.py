"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
180608: Created.
220707,09: Added option to use selected breps for face count instead of entering the value.
250604: Bug fixes.
"""

import Rhino.DocObjects as rd
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iTargetCt'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iCtTol'; keys.append(key)
    values[key] = 0
    riOpts[key] = ri.Custom.OptionInteger(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

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

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'iTargetCt':
            if cls.riOpts[key].CurrentValue < 1:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
        elif key == 'iCtTol':
            if cls.riOpts[key].CurrentValue < 0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
        else:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            elif key in cls.listValues:
                cls.values[key] = idxList
            else:
                return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_FaceCounts():
    """
    Get open breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Optionally select breps to use their face counts")

    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    go.AcceptNothing(True)
    go.AcceptNumber(True, acceptZero=True)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        go.SetCommandPromptDefault("{}".format(Opts.values['iTargetCt']))

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iTargetCt')
        addOption('iCtTol')
        addOption('bEcho')
        addOption('bDebug')
    
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
    
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return [Opts.values['iTargetCt']]

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            iCts_NoTol = sorted(set([o.Brep().Faces.Count for o in objrefs]))
            if len(iCts_NoTol) == 1:
                Opts.setValue('iTargetCt')
            iCtTol = Opts.values['iCtTol']
            if iCtTol == 0:
                return iCts_NoTol
            iCts_WithTol = []
            for iCt in iCts_NoTol:
                for iTol in range(-iCtTol, iCtTol+1):
                    iCts_WithTol.append(iCt + iTol)
            return sorted(set(iCts_WithTol))

        if res == ri.GetResult.Number:
            if Opts.values['bLimitDev']:
                key = 'fDevTol'
                Opts.riOpts[key].CurrentValue = go.Number()
            else:
                key = 'iMinTargetGContinuity'
                Opts.riOpts[key].CurrentValue = int(go.Number())

            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def main():
    
    for gObj in rs.NormalObjects():
        if rs.ObjectType(gObj) in (8,16): break
    else:
        print("No breps are selectable.")
        return

    iFaceCts_Target = getInput_FaceCounts()
    if iFaceCts_Target is None: return

    gBreps = []
    iCts_Faces = []

    for gObj in rs.NormalObjects():
        if rs.ObjectType(gObj) not in (8,16): continue
        gBreps.append(gObj)
        brep = rs.coercebrep(gObj)
        iCts_Faces.append(brep.Faces.Count)
    
    if len(gBreps) == 0:
        print("No breps are selectable.")
        return

    iCtTol = Opts.values['iCtTol']

    idxs_Found = []
    
    for idxB, iCt in enumerate(iCts_Faces):
        for iFaceCt_Target in iFaceCts_Target:
            if (iFaceCt_Target - iCtTol) <= iCt <= (iFaceCt_Target + iCtTol):
                idxs_Found.append(idxB)
                break

    #idxs_Found = [i for i, v in enumerate(iCts_Faces) if v in iFaceCts_Target]
    
    sc.doc.Objects.UnselectAll()
    
    iCts_Faces_Found = []
    for i in idxs_Found:
        sc.doc.Objects.Select(gBreps[i])
        iCts_Faces_Found.append(iCts_Faces[i])
    
    nSelected = len(list(
            sc.doc.Objects.GetSelectedObjects(
            includeLights=False, includeGrips=False)))
    print("{} breps with {} faces are selected.".format(nSelected, iCts_Faces_Found))
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main()