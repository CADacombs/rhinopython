"""
This script allows the modification of a curve's length to be set to a entered value or
to match the length of another curve.

It is an alternative script to Pascal Golay's
RelimitCurve.py (June 2024) https://discourse.mcneel.com/t/questions-about-the-extend-command-for-curves/183480/4
SetCurveLength (November 2022) https://discourse.mcneel.com/t/matching-curve-lengths-with-a-macro/151353/6

Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250727: Started with Pascal Golay's RelimitCurve.py (see above).

TODO:
    Support BrepEdges as input Curves. Their resultant Curves will be added to the document (No Brep modification).
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


    key = 'bBoth'; keys.append(key)
    values[key] = False
    names[key] = 'ModifyEnd'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'NearPick', 'Both')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fLength'; keys.append(key)
    values[key] = 1.0
    names[key] = 'fNewLength'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.ModelUnitSystem)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def _formatDistance(fDistance, iPrecision=None):
    if iPrecision is None:
        iPrecision = sc.doc.ModelDistanceDisplayPrecision

    if fDistance is None:
        return "(No deviation provided)"

    if fDistance < 10**-iPrecision:
        return "{:.{}e}".format(fDistance, 0)

    if fDistance < 0.1:
        return "{:.{}g}".format(fDistance, iPrecision)

    return "{:.{}f}".format(fDistance, iPrecision)


def getInput_CrvsToMod():
    """
    Get curves with optional input.
    """

    go = ri.Custom.GetObject()

    go.GeometryFilter = rd.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    #go.OneByOnePostSelect = True
    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bBoth')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        if Opts.values['bBoth']:
            go.SetCommandPrompt("Select curves to modify")
            res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        else:
            go.SetCommandPrompt("Select curve near its end to modify")
            res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_TargetLength(fStartLength=None):
    """
    Get float with optional input.
    """

    go = ri.Custom.GetObject()


    if fStartLength is None:
        go.SetCommandPrompt("Enter target length or select curve to match its length")
    else:
        go.SetCommandPrompt("Starting length is {}. Enter target length or select curve to match its length".format(
            _formatDistance(fStartLength)))

    go.SetCommandPromptDefault(str(Opts.values['fLength']))

    go.GeometryFilter = rd.ObjectType.Curve

    go.AcceptNumber(True, acceptZero=True)
    go.AcceptNothing(True)

    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return Opts.values['fLength']

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            fLength = objref.Curve().GetLength()
            sc.sticky[Opts.stickyKeys['fLength']] = fLength
            return fLength

        if res == ri.GetResult.Number:
            fLength = go.Number()
            if fLength < Rhino.RhinoMath.ZeroTolerance:
                print("Value must be greater than {}".format(Rhino.RhinoMath.ZeroTolerance))
                continue
            sc.sticky[Opts.stickyKeys['fLength']] = fLength
            return fLength

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _curveEnd_closestToPick(objref_Crv):
    rgCrv = objref_Crv.Curve()

    bSuccess, t = rgCrv.ClosestPoint(objref_Crv.SelectionPoint())
    if not bSuccess:
        return

    t_MidLength = rgCrv.DivideByCount(segmentCount=2, includeEnds=False)[0]

    return rg.CurveEnd.End if t >= t_MidLength else rg.CurveEnd.Start


def modifyCurveLength(crv_In, fLength_Target, end):
    """
    Returns a new rg.Curve
    """

    fLength_Start = crv_In.GetLength()

    fLength_Delta = fLength_Target - fLength_Start

    if fLength_Target < fLength_Start:
        fLength_toRemove = fLength_Start - fLength_Target
        return rg.Curve.Trim(
            crv_In,
            side=end,
            length=fLength_toRemove/2.0 if end==rg.CurveEnd.Both else fLength_toRemove)
    else:
        fLength_toAdd = fLength_Target - fLength_Start
        return rg.Curve.Extend(
            crv_In,
            side=end,
            length=fLength_toAdd/2.0 if end == rg.CurveEnd.Both else fLength_toAdd,
            style=rg.CurveExtensionStyle.Smooth)


def main():

    objrefs = getInput_CrvsToMod()
    if objrefs is None: return

    bBoth = Opts.values['bBoth']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    fLength_Target = getInput_TargetLength(objrefs[0].Curve().GetLength() if objrefs.Count==1 else None)
    if fLength_Target is None: return
    #sEval = "fLength_Target"; print(sEval,'=', eval(sEval))

    bReplace = Opts.values['bReplace']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    if bBoth:
        end = rg.CurveEnd.Both

    gAlreadyAtTarget = []
    gFailed = []
    gMissedTarget_Arc = []
    gMissedTarget_Other = []
    gAdded = []
    gModified = []

    for objref in objrefs:

        crv_In = objref.Curve()

        if not bBoth:
            end = _curveEnd_closestToPick(objrefs[0])


        fLength_In = crv_In.GetLength()

        fLength_Delta = fLength_Target - fLength_In
        if abs(fLength_Delta) <= (0.1 * sc.doc.ModelAbsoluteTolerance):
            gAlreadyAtTarget.append(objref.ObjectId)
            continue

        crv_Res = modifyCurveLength(crv_In, fLength_Target, end)
        if crv_Res is None:
            gFailed.append(objref.ObjectId)
            continue

        fLength_Res = crv_Res.GetLength()

        if abs(fLength_Res - fLength_Target) > (0.1 * sc.doc.ModelAbsoluteTolerance):
            if isinstance(crv_In, rg.ArcCurve):
                gMissedTarget_Arc.append(objref.ObjectId)
            else:
                gMissedTarget_Other.append(objref.ObjectId)
            crv_Res.Dispose()
            continue

        if bReplace:
            if sc.doc.Objects.Replace(objref, crv_Res):
                gModified.append(objref.ObjectId)
            else:
                gFailed.append(objref.ObjectId)
        else:
            gCrv_Out = sc.doc.Objects.AddCurve(crv_Res)
            if gCrv_Out != gCrv_Out.Empty:
                gAdded.append(gCrv_Out)
            else:
                gFailed.append(objref.ObjectId)

    if bEcho:
        if objrefs.Count == 1:
            if gAdded:
                print("Added curve with length of {}".format(
                    _formatDistance(fLength_Res)))
            elif gModified:
                print("Curve length is now {}".format(
                    _formatDistance(fLength_Res)))
            elif gAlreadyAtTarget:
                print("Curve is already at target length, {}.".format(
                    _formatDistance(fLength_Target)))
            elif gMissedTarget_Arc:
                print("Skipped arc curve with length of {}.".format(
                    _formatDistance(fLength_Res)))
            elif gMissedTarget_Other:
                print("Skipped curve with length of {}.".format(
                    _formatDistance(fLength_Res)))
            elif gFailed:
                print("Failed to create curve with new length, {}.".format(
                    _formatDistance(fLength_Target)))
        elif objrefs.Count > 1:
            if gAdded:
                print("Added {} curves each with length of {}".format(
                    len(gAdded),
                    _formatDistance(fLength_Target)))
            if gModified:
                print("Modified length of {} curves to {}".format(
                    len(gModified),
                    _formatDistance(fLength_Target)))
            if gAlreadyAtTarget:
                print("{} curves are already at target length and were not modified.".format(
                    len(gAlreadyAtTarget)))
            if gMissedTarget_Arc:
                print("{} arc curves could not be modified to target length.".format(
                    len(gMissedTarget_Arc)))
            if gMissedTarget_Other:
                print("{} non-arc curves could not be modified to target length.".format(
                    len(gMissedTarget_Other)))
            if gFailed:
                print("{} curves were not modified.".format(len(gFailed)))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()