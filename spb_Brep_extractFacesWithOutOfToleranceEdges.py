"""
This script extact faces with edges that are outside of a specified tolerance.
The extract faces (monofaced breps) can then have their edges rebuilt, e.g., _RebuildEdges.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250602: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xBrepObject


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fMaxTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    lowerLimit_Tol = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    riOpts[key] = ri.Custom.OptionDouble(
        values[key],
        setLowerLimit=True,
        limit=lowerLimit_Tol)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

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

        if key == 'fMaxTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            if cls.riOpts[key].CurrentValue < cls.lowerLimit_Tol:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.lowerLimit_Tol
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none selected")
    
    go.GeometryFilter = rd.ObjectType.Brep

    go.AcceptNothing(True)
    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fMaxTol')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            sc.doc.Objects.UnselectAll()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            sc.doc.Objects.UnselectAll()
            oes = rd.ObjectEnumeratorSettings()
            oes.LockedObjects = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep
            return list(sc.doc.Objects.GetObjectList(oes))

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fMaxTol'
            if key not in idxs_Opts:
                print("Numerical input was not applied.")
                continue
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def findFacesWithOTEdges(rgBrep, fMaxTol, bDebug=False):
    iFs_Found = []

    for rgE in rgBrep.Edges:
        if rgE.Tolerance > fMaxTol:
            iFs_Found.extend(rgE.AdjacentFaces())

    return sorted(set(iFs_Found))


def processBrepObjects(rhBreps, fMaxTol, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBreps0: rd.Breps or ObjRefs of Brep.
    """

    if bDebug: bEcho = True

    def coerceRhinoObject(rhObj):
        rdObj = None
        if isinstance(rhObj, rd.RhinoObject):
            rdObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj)
        return rdObj


    def coerceBrepObject(rhObj):
        rdObj = coerceRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    gBs_Out = [] # Accumulation of test positive, extracted single-face breps.
    iCt_FacesFound_All = 0
    iCt_1FBs_Selected = 0

    len_rhBs0 = len(rhBreps)

    idxs_AtTenths = [int(round(0.1*i*len(rhBreps),0)) for i in range(10)]

    for iB, rhBrep in enumerate(rhBreps):
        if sc.escape_test(False):
            print("Searching interrupted by user.")
            return

        if len_rhBs0 > 10:
            if iB in idxs_AtTenths:
                Rhino.RhinoApp.SetCommandPrompt("Analysis at {:d}% of {} breps ...".format(
                    int(100.0 * (iB+1) / len_rhBs0), len_rhBs0))
        elif len_rhBs0 > 1:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing {} of {} breps".format(
                iB+1, len_rhBs0))
        else:
            Rhino.RhinoApp.SetCommandPrompt("Analyzing brep ...")
        
        rdBrep_In = coerceBrepObject(rhBrep)
        if not rdBrep_In:
            if bEcho and len_rhBs0 == 1:
                print("Non-brep DocObject passed to processBrepObjects.")
            continue
        
        rgB = rdBrep_In.BrepGeometry

        idxs_Found_ThisBrep = findFacesWithOTEdges(rgB, fMaxTol=fMaxTol, bDebug=bDebug)
        if idxs_Found_ThisBrep is None: continue
        if len(idxs_Found_ThisBrep) == 0:
            if bDebug: print("No faces found in brep.")
            continue

        iCt_FacesFound_All += len(idxs_Found_ThisBrep)

        if rgB.Faces.Count == 1:
            if sc.doc.Objects.Select(objectId=rdBrep_In.Id):
                iCt_1FBs_Selected += 1
                gBs_Out.append(rdBrep_In.Id)
            else:
                print("Could not select {}.".format(rdBrep_In.Id))

            continue

        # Faces were extracted from polyfaced brep.

        rvs = xBrepObject.extractFaces(
            rdBrep_In,
            idxFaces=idxs_Found_ThisBrep,
            bAddOnlyMonofaces=True,
            bRetainLayer=True,
            bRetainColor=True,
            bEcho=False,
            bDebug=bDebug)
        if rvs is None:
            continue
        gBreps_Extracted, gBreps_Remaining = rvs

        for gBrep_Extracted in gBreps_Extracted:
            if sc.doc.Objects.Select(objectId=gBrep_Extracted):
                iCt_1FBs_Selected += 1
                gBs_Out.append(rdBrep_In.Id)
            else:
                print("Could not select {}.".format(gBrep_Extracted))


        # End of rhBreps0 loop.

    if not iCt_FacesFound_All:
        if bEcho:
            print("No faces with out-of-tolerance edges found.")
        return gBs_Out

    if iCt_1FBs_Selected:
        print("\n{} monofaced breps are selected.".format(iCt_1FBs_Selected))

    return gBs_Out


def main():

    rhBreps_In = getInput()
    if rhBreps_In is None: return

    fMaxTol = Opts.values['fMaxTol']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBreps_Res = processBrepObjects(
        rhBreps_In,
        fMaxTol=fMaxTol,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()