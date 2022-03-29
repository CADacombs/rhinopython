"""
Matches the end condition (G0, G1, or G2) of one curve to another.

NurbsCurve.SetEndCondition works as so:
    For tangency, the G1 CP([1]) is rotated about the G0 CP.
    For curvature, the G2 CP projection onto the G0-G1 position vector is at
        G0-G1 CP distance from G1 for Bezier curves
        2 * (G0-G1) distance from G1 for non-Bezier curves
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220328: Created.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Enum


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iContinuity'; keys.append(key)
    values[key] = 1
    listValues[key] = Enum.GetNames(rg.NurbsCurve.NurbsCurveEndConditionType)[1:]
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

        if not idxOpt: print("Add option for {} failed.".format(key))

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fDevTol':
            if cls.riOpts[key].CurrentValue == 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = Rhino.RhinoMath.ZeroTolerance
            elif cls.riOpts[key].CurrentValue <= 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
            else:
                cls.values[key] = cls.riOpts[key].CurrentValue
        elif key == 'iMinTargetGContinuity':
            cls.values[key] = cls.riOpts[key].CurrentValue
        else:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            elif key in cls.listValues:
                cls.values[key] = idxList
            else:
                return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve and parameter with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves near their matching ends, 1st matching to 2nd")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.OneByOnePostSelect = True
    go.DisablePreSelect()

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

        addOption('iContinuity')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected or a number was entered.

        if res == ri.GetResult.Number:
            if int(go.Number()) not in (0,1,2):
                print("Numeric input is invalid.")
                continue
            Opts.setValue('iContinuity', idxList=int(go.Number()))
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def main():
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script is supported only in Rhino V6 and above.")
        return

    objrefs_In = getInput()
    if objrefs_In is None: return


    iContinuity = Opts.values['iContinuity']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    gCrv0_A = objrefs_In[0].ObjectId
    c0_A = objrefs_In[0].Curve()
    
    bSuccess, tA = c0_A.ClosestPoint(objrefs_In[0].SelectionPoint())
    if not bSuccess:
        c0_A.Dispose()
        return

    gCrv0_B = objrefs_In[1].ObjectId
    c0_B = objrefs_In[1].Curve()
    bSuccess, tB = c0_B.ClosestPoint(objrefs_In[1].SelectionPoint())
    if not bSuccess:
        c0_A.Dispose()
        c0_B.Dispose()
        return
    
    bT1WorkEnd_A = tA > c0_A.Domain.Mid
    bT1WorkEnd_B = tB > c0_B.Domain.Mid

    nc_In = c0_A.ToNurbsCurve()
    nc_WIP = c0_A.DuplicateCurve()
    nc_Ref = c0_B.ToNurbsCurve()

    continuity = Enum.ToObject(
        rg.NurbsCurve.NurbsCurveEndConditionType, iContinuity+1) # +1 because None is 0.

    if bT1WorkEnd_B:
        point = nc_Ref.PointAtEnd
        tangent = nc_Ref.TangentAtEnd
        curvature = nc_Ref.CurvatureAt(nc_Ref.Domain.T1)
    else:
        point = nc_Ref.PointAtStart
        tangent = nc_Ref.TangentAtStart
        curvature = nc_Ref.CurvatureAt(nc_Ref.Domain.T0)

    bSuccess = nc_WIP.SetEndCondition(
        bT1WorkEnd_A,
        continuity,
        point,
        tangent,
        curvature)
    if not bSuccess:
        print("SetEndCondition failed.")

    if nc_WIP.EpsilonEquals(nc_In, epsilon=Rhino.RhinoMath.ZeroTolerance):
        print("No change in NurbsCurve.")
        return

    if sc.doc.Objects.Replace(objrefs_In[0], curve=nc_WIP):
        sc.doc.Views.Redraw()
        print("Curve was replaced.")
    else:
        print("Curve could not be replaced.")


if __name__ == '__main__': main()