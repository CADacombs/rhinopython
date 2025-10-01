"""
Send any questions, comments, or script development service needs to
@spb on the McNeel Forums, https://discourse.mcneel.com/
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
221114: Created.
221118: Bug fix: Added a View.Redraw.
221119: Added option to dot edges.
221204: Bug fix: Now dots are added to current layer.
251001: Modified an option default value.
"""

import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fMaxLength'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance #1e-4
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bShortestOnly'; keys.append(key)
    values[key] = False
    names[key] = 'Mark'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'All', 'Shortest')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDupCrv'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bIncludeIndex'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHt'; keys.append(key)
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

        if key == 'fMaxLength':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

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
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False

    go.AcceptNothing(True)

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()


        addOption('fMaxLength')
        addOption('bShortestOnly')
        addOption('bDupCrv')
        addOption('bAddDot')
        if Opts.values['bAddDot']:
            addOption('bIncludeIndex')
            addOption('iDotHt')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return getAllNormalBreps()

        if res == ri.GetResult.Object:
            rdObjs = [o.Object() for o in go.Objects()]
            go.Dispose()
            return rdObjs

        if res == ri.GetResult.Number:
            key = 'fMaxLength'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def formatDistance(fDistance, iPrecision=None):
    if iPrecision is None:
        iPrecision = sc.doc.ModelDistanceDisplayPrecision
    if fDistance is None:
        return "(No value provided)"
    if fDistance == 0.0:
        return "0"
    if fDistance < 0.01:
        return "{:.2e}".format(fDistance)
    return "{:.{}g}".format(fDistance, iPrecision)


def main():

    rdBs_In = getInput()
    if rdBs_In is None: return

    fMaxLength = Opts.values['fMaxLength']
    bShortestOnly = Opts.values['bShortestOnly']
    bDupCrv = Opts.values['bDupCrv']
    bAddDot = Opts.values['bAddDot']
    bIncludeIndex = Opts.values['bIncludeIndex']
    iDotHt = Opts.values['iDotHt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    fLengths_Short = []
    gBs_WithShort = []
    idxEs_Short = []

    fMinLength_All = float("inf")

    bReportEveryBrep = False
    if 1 <= len(rdBs_In) <= 10:
        bReportEveryBrep = True
    elif 10 < len(rdBs_In):
        s  = "More than 10 breps are selected,"
        s += " so only maximum of all edge tolerances,"
        s += " is reported."
        print(s)
    else:
        return

    gCrvs_Out = []
    gDots_Out = []
    if bAddDot:
        attrib_Dot = rd.ObjectAttributes()
        attrib_Dot.ColorSource = rd.ObjectColorSource.ColorFromObject
        attrib_Dot.ObjectColor = Color.Red
        attrib_Dot.LayerIndex = sc.doc.Layers.CurrentLayerIndex


    epsilon_Shortest = 0.5 * 10.0**-sc.doc.ModelDistanceDisplayPrecision
    print("TODO: determine epsilon for shortest.  Now is {}.".format(epsilon_Shortest))

    for iB, rdB in enumerate(rdBs_In):
        gB = rdB.Id
        rgB = rdB.BrepGeometry

        for rgE in rgB.Edges:
            fLength = rgE.GetLength()
            if fLength > fMaxLength: continue

            fLengths_Short.append(fLength)
            gBs_WithShort.append(gB)
            idxEs_Short.append([rgE.EdgeIndex])

            if fLength < fMinLength_All:
                fMinLength_All = fLength

            if bShortestOnly:
                continue # Shortest will be marked after for loops.

            if bDupCrv:
                gCrvs_Out.append(sc.doc.Objects.AddCurve(rgE))

            if bAddDot:
                text = formatDistance(fLength)
                if bIncludeIndex:
                    text = "e[{}]:".format(rgE.EdgeIndex) + text
                location = rgE.PointAt(rgE.Domain.Mid)
                rgDot = rg.TextDot(text, location)
                rgDot.FontHeight = iDotHt
                gDots_Out.append(sc.doc.Objects.AddTextDot(rgDot, attrib_Dot))



    if not fLengths_Short:
        print("No edge lengths <= {} found.".format(formatDistance(fMaxLength)))
        return

    print("Found {} edges in {} breps with lengths [{},{}].".format(
        len(idxEs_Short),
        len(set(gBs_WithShort)),
        formatDistance(min(fLengths_Short)),
        formatDistance(max(fLengths_Short))))

    if not bShortestOnly:
        sc.doc.Views.Redraw()
        return


    # Find edges closest to maximum.
    fLengths_Shortest = []
    gBs_WithShortest = []
    idxEs_Shortest = []

    for fLength, gB, idxE in zip(fLengths_Short, gBs_WithShort, idxEs_Short):
        if (fLength - fMinLength_All) < epsilon_Shortest:
            fLengths_Shortest.append(fLength)
            gBs_WithShortest.append(gBs_WithShortest)
            idxEs_Shortest.append(idxEs_Shortest)

            if bDupCrv:
                gCrvs_Out.append(sc.doc.Objects.AddCurve(rgE))

            if bAddDot:
                text = formatDistance(fLength)
                if bIncludeIndex:
                    text = "e[{}]:".format(rgE.EdgeIndex) + text
                location = rgE.PointAt(rgE.Domain.Mid)
                rgDot = rg.TextDot(text, location)
                rgDot.FontHeight = iDotHt
                gDots_Out.append(sc.doc.Objects.AddTextDot(rgDot, attrib_Dot))


    print("Found {} edges in {} breps within {} of {}.".format(
        len(idxEs_Shortest),
        len(set(gBs_WithShortest)),
        formatDistance(epsilon_Shortest),
        formatDistance(fLengths_Shortest)))

    if bDupCrv:
        print("{} curves added.".format(len(gDots_Out)))
    if bAddDot:
        print("{} dots added.".format(len(gDots_Out)))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()