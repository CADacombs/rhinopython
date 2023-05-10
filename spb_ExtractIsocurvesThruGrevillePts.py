"""
Some of the output curves of this script can be modified (with limitations),
and a new shaped surface with the same knot structure can be created with _Sweep1
using its default options.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
230508-09: Created.
"""

import Rhino
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


    key = 'bAlongU'; keys.append(key)
    values[key] = True
    names[key] = 'U'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAlongV'; keys.append(key)
    values[key] = True
    names[key] = 'V'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bIncludeOppDirBorder'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


class DrawConduit(Rhino.Display.DisplayConduit):

    def __init__(self):
        self.crvs = None
        self.colors = None
        displayMode = Rhino.RhinoDoc.ActiveDoc.Views.ActiveView.ActiveViewport.DisplayMode
        self.crv_thk = displayMode.DisplayAttributes.CurveThickness + 2

    def CalculateBoundingBox(self, calculateBoundingBoxEventArgs):
        if not self.crvs: return

        for curve in self.crvs:
            self.bbox = curve.GetBoundingBox(accurate=False)
            calculateBoundingBoxEventArgs.IncludeBoundingBox(self.bbox)

    def PreDrawObjects(self, drawEventArgs):
        if not self.crvs: return

        for iC, curve in enumerate(self.crvs):
            color = self.colors[iC]

            drawEventArgs.Display.DrawCurve(
                curve,
                color,
                thickness=self.crv_thk)


def createCrvs(ns, iDir, bOnlyBorder=False):
    """
    """

    nc_Out = []

    if not iDir:
        if bOnlyBorder:
            u, v = ns.Points.GetGrevillePoint(0, 0)
            nc = ns.IsoCurve(0, v)
            nc_Out.append(nc)
            u, v = ns.Points.GetGrevillePoint(0, ns.Points.CountV-1)
            nc = ns.IsoCurve(0, v)
            nc_Out.append(nc)
            return nc_Out

        for iV in range(ns.Points.CountV):
            u, v = ns.Points.GetGrevillePoint(0, iV)
            nc = ns.IsoCurve(0, v)
            nc_Out.append(nc)
        return nc_Out

    if bOnlyBorder:
        u, v = ns.Points.GetGrevillePoint(0, 0)
        nc = ns.IsoCurve(1, u)
        nc_Out.append(nc)
        u, v = ns.Points.GetGrevillePoint(ns.Points.CountU-1, 0)
        nc = ns.IsoCurve(1, u)
        nc_Out.append(nc)
        return nc_Out

    for iU in range(ns.Points.CountU):
        u, v = ns.Points.GetGrevillePoint(iU, 0)
        nc = ns.IsoCurve(1, u)
        nc_Out.append(nc)

    return nc_Out


def main():

    sk_conduit = 'conduit({})'.format(__file__)
    if (sk_conduit in sc.sticky) and sc.sticky[sk_conduit]:
        conduit = sc.sticky[sk_conduit]
        conduit.Enabled = False
        sc.doc.Views.Redraw()
    else:
        conduit = DrawConduit()
        sc.sticky[sk_conduit] = conduit



    rc, objref = ri.RhinoGet.GetOneObject(
        "Select face",
        acceptNothing=False,
        filter=Rhino.DocObjects.ObjectType.Surface)
    if rc != Rhino.Commands.Result.Success:
        return


    face = objref.Surface()
    srf = face.UnderlyingSurface()

    go = ri.Custom.GetPoint()

    go.SetCommandPrompt("Pick to cycle through directions")
    go.SetCommandPromptDefault("Create curves")

    go.AcceptNothing(True)

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)


    ncs_U = []
    ncs_V = []


    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        addOption('bAlongU')
        addOption('bAlongV')
        if not (Opts.values['bAlongU'] and Opts.values['bAlongV']):
            addOption('bIncludeOppDirBorder')
        addOption('bEcho')
        addOption('bDebug')


        if ncs_U:
            for nc in ncs_U: nc.Dispose()
            ncs_U = []
        if ncs_V:
            for nc in ncs_V: nc.Dispose()
            ncs_V = []
        colors = []

        if Opts.values['bAlongU']:
            ncs_U = createCrvs(
                srf,
                iDir=0,
                bOnlyBorder=False)
            colors.extend([Color.Red] * len(ncs_U))
        elif Opts.values['bIncludeOppDirBorder']:
            ncs_U = createCrvs(
                srf,
                iDir=0,
                bOnlyBorder=True)
            colors.extend([Color.Red] * len(ncs_U))

        if Opts.values['bAlongV']:
            ncs_V = createCrvs(
                srf,
                iDir=1,
                bOnlyBorder=False)
            colors.extend([Color.Lime] * len(ncs_V))
        elif Opts.values['bIncludeOppDirBorder']:
            ncs_V = createCrvs(
                srf,
                iDir=1,
                bOnlyBorder=True)
            colors.extend([Color.Red] * len(ncs_U))

        if ncs_U or ncs_V:
            conduit.crvs = ncs_U + ncs_V
            conduit.colors = colors
            conduit.Enabled = True
            sc.doc.Views.Redraw()


        res = go.Get()

        if res == ri.GetResult.Cancel:
            conduit.Enabled = False
            sc.doc.Views.Redraw()
            go.Dispose()
            for nc in ncs_U + ncs_V: nc.Dispose()
            return

        if res == ri.GetResult.Nothing:
            conduit.Enabled = False
            sc.doc.Views.Redraw()
            go.Dispose()
            break

        if res == ri.GetResult.Point:
            if Opts.values['bAlongU'] and Opts.values['bAlongV']:
                Opts.riOpts['bAlongU'].CurrentValue = True
                Opts.setValue('bAlongU')
                Opts.riOpts['bAlongV'].CurrentValue = False
                Opts.setValue('bAlongV')
                continue
            if Opts.values['bAlongU'] and not Opts.values['bAlongV']:
                Opts.riOpts['bAlongU'].CurrentValue = False
                Opts.setValue('bAlongU')
                Opts.riOpts['bAlongV'].CurrentValue = True
                Opts.setValue('bAlongV')
                continue
            if not Opts.values['bAlongU'] and Opts.values['bAlongV']:
                Opts.riOpts['bAlongU'].CurrentValue = True
                Opts.setValue('bAlongU')
                Opts.riOpts['bAlongV'].CurrentValue = True
                Opts.setValue('bAlongV')
                continue
            # Fix when both are False.
            Opts.riOpts['bAlongU'].CurrentValue = True
            Opts.setValue('bAlongU')
            Opts.riOpts['bAlongV'].CurrentValue = True
            Opts.setValue('bAlongV')
            continue


        # An option was selected.

        if go.Option().Index == idxs_Opt['bAlongU']:
            if Opts.values['bAlongU'] and not Opts.values['bAlongV']:
                Opts.riOpts['bAlongV'].CurrentValue = True
                Opts.setValue('bAlongV')
            Opts.setValue('bAlongU')
            continue

        if go.Option().Index == idxs_Opt['bAlongV']:
            if Opts.values['bAlongV'] and not Opts.values['bAlongU']:
                Opts.riOpts['bAlongU'].CurrentValue = True
                Opts.setValue('bAlongU')
            Opts.setValue('bAlongV')
            continue

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


    gCs_Out = []

    for nc in ncs_U + ncs_V:
        gC_Out = sc.doc.Objects.AddCurve(nc)
        if gC_Out == gC_Out.Empty:
            if Opts.values['bEcho']: print("Curve could not be added.")
        else:
            gCs_Out.append(gC_Out)

    sc.doc.Views.Redraw()

    if Opts.values['bEcho']:
        print("{} curves added.".format(len(gCs_Out)))


if __name__ == '__main__': main()