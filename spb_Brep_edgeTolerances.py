"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
200131: Created.
200507: Bug fix.
221104: Added options for marking tolerances using dots.
241105: Modified printed output to include range and count of edges above mdoel tolerance.
241105: Modified printed output of when float vs. scientific notation is used.
250301: Modified formatDistance.
250602: Dots are now red.
"""

import Rhino
import Rhino.DocObjects as rd
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


    key = 'fMaxTol'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key], setLowerLimit=True, limit=1e-6)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDotOnlyAboveTol'; keys.append(key)
    values[key] = True
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

        if key == 'fMaxTol':
            if cls.riOpts[key].CurrentValue <= 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
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

    go.SetCommandPrompt("Select brep")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep
    #go.SubObjectSelect = False

    go.AcceptNothing(True)

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fMaxTol')
        addOption('bAddDot')
        if Opts.values['bAddDot']:
            addOption('bDotOnlyAboveTol')
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
            objrefs = go.Objects()
            go.Dispose()

            return objrefs

        if res == ri.GetResult.Number:
            key = 'fMaxTol'
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
    if fDistance < 0.001:
        return "{:.3e}".format(fDistance)
    return "{:.{}g}".format(fDistance, iPrecision)


def main():

    rhBs_In = getInput()
    if rhBs_In is None: return

    bAddDot = Opts.values['bAddDot']
    bDotOnlyAboveTol = Opts.values['bDotOnlyAboveTol']
    fMaxTol = Opts.values['fMaxTol']
    iDotHt = Opts.values['iDotHt']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    fMaxEdgeTol_All = 0.0
    iCt_Edges_All = 0
    iCt_Edges_Found_All = 0
    
    bReportEveryBrep = False
    if 1 <= len(rhBs_In) <= 10:
        bReportEveryBrep = True
    elif 10 < len(rhBs_In):
        s  = "More than 10 breps are selected,"
        s += " so only maximum of all edge tolerances,"
        s += " is reported."
        print(s)
    else:
        return

    if isinstance(rhBs_In[0], rd.ObjRef):
        rgBs = [_.Brep() for _ in rhBs_In]
    elif isinstance(rhBs_In[0], rd.BrepObject):
        rgBs = [_.BrepGeometry for _ in rhBs_In]

    gDots = []

    if bAddDot:
        attr_Out = rd.ObjectAttributes()
        attr_Out.LayerIndex = sc.doc.Layers.CurrentLayerIndex
        attr_Out.ColorSource = rd.ObjectColorSource.ColorFromObject
        attr_Out.ObjectColor = Color.FromArgb(255,0,0)

    fMaxEdgeTols = []
    fOutOfTols_All = []

    for rgB in rgBs:

        fEdgeTols = []
        #fMaxEdgeTol = 0.0
        fOutOfTols_1B = []

        for rgE in rgB.Edges:
            edgeTol = rgE.Tolerance
            fEdgeTols.append(edgeTol)
            if edgeTol > fMaxTol:
                fOutOfTols_1B.append(edgeTol)

            #            if edgeTol > fMaxEdgeTol:
            #                fMaxEdgeTol = edgeTol

            if not bAddDot: continue

            if bDotOnlyAboveTol and edgeTol <= fMaxTol: continue

            text = formatDistance(edgeTol)
            location = rgE.PointAt(rgE.Domain.Mid)
            rgDot = Rhino.Geometry.TextDot(text, location)
            rgDot.FontHeight = 11
            gDot = sc.doc.Objects.AddTextDot(rgDot, attributes=attr_Out)
            if gDot != gDot.Empty:
                gDots.append(gDot)

        fMaxEdgeTol = max(fEdgeTols)
        fMinEdgeTol = min(fEdgeTols)
        #        fMaxEdgeTols.append(fMaxEdgeTol)
        fOutOfTols_All.extend(fOutOfTols_1B)

        #        if fMaxEdgeTol > fMaxEdgeTol_All:
        #            fMaxEdgeTol_All = fMaxEdgeTol
        
        if bReportEveryBrep:
            if not fOutOfTols_1B:
                print("No edge tolerances above {} in brep with {} edges.".format(
                    fMaxTol,
                    rgB.Edges.Count,
                    ))
                #            print(
                #                "[{},{}] = edge tolerances found in brep with {} edges.".format(
                #                formatDistance(fMinEdgeTol),
                #                formatDistance(fMaxEdgeTol),
                #                rgB.Edges.Count,
                #                ),
            else:
                if len(fOutOfTols_1B) == 1:
                    print(
                        "{} = only edge tolerance above {} in brep with {} edges.".format(
                        formatDistance(fOutOfTols_1B[0]),
                        fMaxTol,
                        rgB.Edges.Count,
                        )
                        )
                else:
                    print(
                        "[{},{}] = range of {} edge tolerances above {} in brep with {} edges.".format(
                        formatDistance(min(fOutOfTols_1B)),
                        formatDistance(max(fOutOfTols_1B)),
                        len(fOutOfTols_1B),
                        fMaxTol,
                        rgB.Edges.Count,
                        )
                        )
            #            print("{} maximum edge tolerance found in brep with {} edges.".format(
            #                formatDistance(fMaxEdgeTol),
            #                rgB.Edges.Count,
            #                ))

        iCt_Edges_All += rgB.Edges.Count

    #    fMaxEdgeTol_All = max(fMaxEdgeTols)


    #    # Find number of edges close to maximum.
    #    for rgB in rgBs:
    #        for rgE in rgB.Edges:
    #            if (
    #                abs(fMaxEdgeTol_All - rgE.Tolerance) <
    #                (0.5 * 10.0**-sc.doc.ModelDistanceDisplayPrecision)
    #            ):
    #                iCt_Edges_Found_All += 1
    #
    #
    #    if len(rgBs) > 1:
    #        print("{} maximum edge tolerance found in {} of {} edges in {} breps.".format(
    #            formatDistance(fMaxEdgeTol_All),
    #            iCt_Edges_Found_All,
    #            iCt_Edges_All,
    #            len(objref_Bs)))

    if not bReportEveryBrep and not fOutOfTols_All:
        print('-'*80)
        print("No edge tolerances above {}.".format(
            fMaxTol))
        return

    if len(rgBs) > 1:
        if not fOutOfTols_All:
            print('-'*80)
            print("No edge tolerances above {}.".format(
                fMaxTol))
            return

        print('Total','-'*40)
        print(
            "[{},{}] = range of {} edge tolerances above {} in {} breps.".format(
            formatDistance(min(fOutOfTols_All)),
            formatDistance(max(fOutOfTols_All)),
            len(fOutOfTols_All),
            fMaxTol,
            len(rgBs),
            )
            )

    if gDots:
        print("{} dots added.".format(len(gDots)))
        sc.doc.Views.Redraw()
    elif bAddDot:
        print("No dots added.")


if __name__ == '__main__': main()