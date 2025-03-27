"""
This script reports and selects faces of breps with either invalid surfaces
or the surfaces convert to invalid breps.

After running the script, convert any bad, non-NURBS surfaces found to NurbsSurfaces.

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250325-26: Created.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List
from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    idxOpt = {}
    stickyKeys = {}


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
    def addOptions(cls, go):
        go.ClearCommandOptions()

        keys_Current = []
        keys_Current.append('bEcho')
        keys_Current.append('bDebug')

        for key in keys_Current:
            cls.idxOpt[key] = None

            if key in cls.riOpts:
                if key[0] == 'b':
                    cls.idxOpt[key] = go.AddOptionToggle(
                            cls.names[key], cls.riOpts[key])[0]
                elif key[0] == 'f':
                    cls.idxOpt[key] = go.AddOptionDouble(
                        cls.names[key], cls.riOpts[key])[0]
                elif key[0] == 'i':
                    cls.idxOpt[key] = go.AddOptionInteger(
                        englishName=cls.names[key], intValue=cls.riOpts[key])[0]
            else:
                cls.idxOpt[key] = go.AddOptionList(
                    englishOptionName=cls.names[key],
                    listValues=cls.listValues[key],
                    listCurrentIndex=cls.values[key])


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


    @classmethod
    def processInput(cls, go):
        res = go.Result()

        cls.setValues()
        cls.saveSticky()

        # Refresh options regardless if a number was entered or options were modified in another way.
        Opts.addOptions(go)


def getInput():
    """
    Get breps with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none selected")
    
    go.GeometryFilter = rd.ObjectType.Brep

    go.AcceptNothing(True)

    Opts.addOptions(go)

    while True:

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            sc.doc.Objects.UnselectAll()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Nothing:
            go.Dispose()
            sc.doc.Objects.UnselectAll()
            oes = rd.ObjectEnumeratorSettings()
            oes.LockedObjects = False
            oes.ObjectTypeFilter = rd.ObjectType.Brep
            rdBs = list(sc.doc.Objects.GetObjectList(oes))
            return rdBs

        # An option was selected.
        Opts.processInput(go)


def findFacesWithBadSrfsForBreps(rgBrep, bDebug=False):
    dict_idxFs = {
        'InvalidNurbsSrfs': [],
        'InvalidPlaneSrfs': [],
        'InvalidRevSrfs': [],
        'InvalidSumSrfs': [],
        'NurbsSrfsMakeInvalidBrep': [],
        'PlaneSrfsMakeInvalidBrep': [],
        'RevSrfsMakeInvalidBrep': [],
        'SumSrfsMakeInvalidBreps': [],
        }

    for rgF in rgBrep.Faces:
        rgS = rgF.UnderlyingSurface()
        if not rgS.IsValid:
            if isinstance(rgS, rg.NurbsSurface):
                dict_idxFs['InvalidNurbsSrfs'].append(rgF.FaceIndex)
                continue
            if isinstance(rgS, rg.PlaneSurface):
                dict_idxFs['InvalidPlaneSrfs'].append(rgF.FaceIndex)
                continue
            if isinstance(rgS, rg.RevSurface):
                dict_idxFs['InvalidRevSrfs'].append(rgF.FaceIndex)
                continue
            if isinstance(rgS, rg.SumSurface):
                dict_idxFs['InvalidSumSrfs'].append(rgF.FaceIndex)
                continue
        _rgB = rgS.ToBrep()
        if not _rgB.IsValid:
            _rgB.Dispose()
            if isinstance(rgS, rg.NurbsSurface):
                dict_idxFs['NurbsSrfsMakeInvalidBrep'].append(rgF.FaceIndex)
                continue
            if isinstance(rgS, rg.PlaneSurface):
                dict_idxFs['PlaneSrfsMakeInvalidBrep'].append(rgF.FaceIndex)
                continue
            if isinstance(rgS, rg.RevSurface):
                dict_idxFs['RevSrfsMakeInvalidBrep'].append(rgF.FaceIndex)
                continue
            if isinstance(rgS, rg.SumSurface):
                dict_idxFs['SumSrfsMakeInvalidBreps'].append(rgF.FaceIndex)
                continue
        _rgB.Dispose()

    return dict_idxFs


def processBrepObjects(rhBreps, bEcho=True, bDebug=False):
    """
    Parameters:
        rhBreps0: rd.Breps or ObjRefs of Brep.
    """

    if bDebug: bEcho = True

    sc.doc.Objects.UnselectAll()

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


    keys = [
        'InvalidNurbsSrfs',
        'InvalidPlaneSrfs',
        'InvalidRevSrfs',
        'InvalidSumSrfs',
        'NurbsSrfsMakeInvalidBrep',
        'PlaneSrfsMakeInvalidBrep',
        'RevSrfsMakeInvalidBrep',
        'SumSrfsMakeInvalidBreps',
        ]

    dict_counts = {}
    for key in keys:
        dict_counts[key] = 0


    gBs_Out = [] # Accumulation of test positive, extracted single-face breps.
    iCt_FacesFound_All = 0
    iCt_Bs_Selected = 0
    iCt_Fs_Selected = 0

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

        rv = findFacesWithBadSrfsForBreps(rgB)
        if rv is None: continue
        if len(rv) == 0:
            print("Empty dict.")
            continue

        dict_Res = rv
        idxs_Found_ThisBrep = []

        if bEcho:
            ss = []
            ss.append("Counts found:")
        for key in keys:
            if dict_Res[key]:
                idxs_Found_ThisBrep.extend(dict_Res[key])
                dict_counts[key] += len(dict_Res[key])
                if bEcho: ss.append("  {}: {}".format(key, len(dict_Res[key])))
        set_Found = set(idxs_Found_ThisBrep)
        if len(set_Found) < len(idxs_Found_ThisBrep):
            print("Duplicate faces found. How did that happen?")

        if len(idxs_Found_ThisBrep) == 0:
            continue # to next brep.

        if bEcho:
            print("\n".join(ss))

        iCt_FacesFound_All += len(idxs_Found_ThisBrep)

        # What was the purpose of this?
        #if iCt_Bs_Selected + iCt_Fs_Selected == 0:
        #    sc.doc.Objects.UnselectAll()

        if rgB.Faces.Count == 1:
            sc.doc.Objects.Select(objectId=rdBrep_In.Id)
            iCt_Bs_Selected += 1
        else:
            for idxF in idxs_Found_ThisBrep:
                compIdx = Rhino.Geometry.ComponentIndex(
                        Rhino.Geometry.ComponentIndexType.BrepFace,
                        index=idxF)
                rdBrep_In.SelectSubObject(
                        componentIndex=compIdx,
                        select=True, syncHighlight=True, persistentSelect=True)
                iCt_Fs_Selected += 1
        gBs_Out.append(rdBrep_In.Id)
        
        # End of rhBreps0 loop.

    if not iCt_FacesFound_All:
        if bEcho:
            print("No faces with bad surfaces found.")
        return []
    
    if iCt_Bs_Selected and iCt_Fs_Selected:
        print("\n{} monoface breps and {} faces in polysrfs are selected.".format(iCt_Bs_Selected, iCt_Fs_Selected))
    elif iCt_Bs_Selected:
        print("\n{} monoface breps are selected.".format(iCt_Bs_Selected))
    elif iCt_Fs_Selected:
        print("\n{} faces in polysrfs are selected.".format(iCt_Fs_Selected))

    if bEcho and (iCt_FacesFound_All < (iCt_Bs_Selected + iCt_Fs_Selected)):
        print("But {} faces with bad surfaces found. Why?".format(iCt_FacesFound_All))

    return gBs_Out


def main():

    rhBreps_In = getInput()
    if rhBreps_In is None: return

    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    gBreps_Res = processBrepObjects(
        rhBreps_In,
        bEcho=bEcho,
        bDebug=bDebug)

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()