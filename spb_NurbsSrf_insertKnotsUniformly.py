"""
Insert knots uniformly within existing spans.
Therefore, when the starting surface is uniform, uniformity will be maintained.
"""
"""
210802: Created.
211105: Now modifies a single face.
        Added direction option after the face is picked.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Drawing


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'iQtyToAddInEachSpan'; keys.append(key)
    values[key] = 2
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=1)
    stickyKeys[key] = '{}({})'.format(key, __file__)

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

        if not idxOpt: print "Add option for {} failed.".format(key)

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        #if key == 'fScale':
        #    if cls.riOpts[key].CurrentValue <= 1e-9:
        #        print "Invalid input for scale value."
        #        cls.riOpts[key].CurrentValue = cls.values[key]

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


class DrawUVArrowsConduit(Rhino.Display.DisplayConduit):

    def __init__(self, ns):
        self.locs_U = []
        self.dirs_U = []
        self.locs_V = []
        self.dirs_V = []
        for u in ns.Domain(0).Min, ns.Domain(0).Mid, ns.Domain(0).Max:
            for v in ns.Domain(1).Min, ns.Domain(1).Mid, ns.Domain(1).Max:
                location = ns.PointAt(u, v)
                isocurve = ns.IsoCurve(direction=0, constantParameter=v)
                direction = isocurve.TangentAt(u)
                self.locs_U.append(location)
                self.dirs_U.append(direction)
                isocurve = ns.IsoCurve(direction=1, constantParameter=u)
                direction = isocurve.TangentAt(v)
                self.locs_V.append(location)
                self.dirs_V.append(direction)

    def DrawForeground(self, e):
        rhRed = Drawing.Color.FromArgb(red=200, green=0, blue=0)
        rhGreen = Drawing.Color.FromArgb(red=0, green=127, blue=0)
        zipped = zip(self.locs_U, self.dirs_U, self.dirs_V)
        for location, direction_U, direction_V in zipped:
            e.Display.DrawDirectionArrow(
                location,
                direction=direction_U,
                color=rhRed)
            e.Display.DrawDirectionArrow(
                location,
                direction=direction_V,
                color=rhGreen)


def getInput():
    """
    Get face with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select brep face")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Surface

    res = go.Get()

    if res == ri.GetResult.Cancel:
        go.Dispose()
        return

    objref = go.Object(0)

    go.Dispose()

    srf = objref.Face().UnderlyingSurface()
            
    if not isinstance(srf, rg.NurbsSurface):
        s = "Underlying surface is a {}".format(srf.GetType().Name)
        s += " and will be converted to a NurbsSurface."
        print s

        ns = srf.ToNurbsSurface()
    else:
        ns = srf

    print "U direction: Degree:{}, SpanCt:{}".format(
        ns.Degree(0), ns.SpanCount(0))
    print "V direction: Degree:{}, SpanCt:{}".format(
        ns.Degree(1), ns.SpanCount(1))

    conduit = DrawUVArrowsConduit(ns)

    go = ri.Custom.GetOption()

    go.AcceptNothing(True)

    go.AcceptNumber(True, acceptZero=False)

    go.SetCommandPrompt("Direction")

    sDirs = 'U', 'V', 'Both'

    key = 'iDirection'
    stickyKey = '{}({})'.format(key, __file__)
    if stickyKey in sc.sticky:
        idxDir_Default = sc.sticky[stickyKey]
    else:
        idxDir_Default = 0

    go.SetCommandPromptDefault(sDirs[idxDir_Default])

    idxs_Opt = {}

    conduit.Enabled = True
    sc.doc.Views.Redraw()

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        for sDir in sDirs:
            idxs_Opt[sDir] = go.AddOption(sDir)


        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('iQtyToAddInEachSpan')
        addOption('bReplace')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        conduit.Enabled = False

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            sDir = sDirs[idxDir_Default]
            break

        if res == ri.GetResult.Number:
            key = 'iQtyToAddInEachSpan'
            if int(go.Number()) < 1:
                continue
            Opts.riOpts[key].CurrentValue = int(go.Number())
            Opts.setValue(key)
            continue

        # An option was selected.
        idx_Opt = go.Option().Index
        if idx_Opt in (1,2,3):
            idx_Dir = go.Option().Index - 1 # Notice the '-1'.
            sDir = sDirs[idx_Dir] 
            sc.sticky[stickyKey] = idx_Dir
            go.Dispose()
            break # out of while loop.

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break # out of for look back to while loop.


    return (
        objref,
        sDir,
        Opts.values['iQtyToAddInEachSpan'],
        Opts.values['bReplace'],
        Opts.values['bEcho'],
        Opts.values['bDebug'],
        )


def createSurface(ns_In, sDir, iQtyToAddInEachSpan, bDebug=False):
    """
    Parameters:
        ns_In: rg.NurbsSurface,
        sDir: str('U', 'V', or 'Both')
        iQtyToAddInEachSpan: int
        bDebug: bool

    Returns: rg.NurbsSurface
    """

    if not isinstance(ns_In, rg.NurbsSurface):
        return

    if sDir not in ('U', 'V', 'Both'):
        return

    if iQtyToAddInEachSpan < 1:
        return

    ns_Out = ns_In.Duplicate()

    # Add knots from Domain end to beginning.

    if sDir in ('U', 'Both'):
        for iK in range(ns_In.KnotsU.Count-1, 0, -1):

            k_R = ns_In.KnotsU[iK]
            k_L = ns_In.KnotsU[iK-1]

            for i in range(iQtyToAddInEachSpan):
                fraction_from_R = float(i+1) / float(iQtyToAddInEachSpan+1)
                k_M = fraction_from_R*k_R + (1.0-fraction_from_R)*k_L
                ns_Out.KnotsU.InsertKnot(k_M)

    if sDir in ('V', 'Both'):
        for iK in range(ns_In.KnotsV.Count-1, 0, -1):

            k_R = ns_In.KnotsV[iK]
            k_L = ns_In.KnotsV[iK-1]

            for i in range(iQtyToAddInEachSpan):
                fraction_from_R = float(i+1) / float(iQtyToAddInEachSpan+1)
                k_M = fraction_from_R*k_R + (1.0-fraction_from_R)*k_L
                ns_Out.KnotsV.InsertKnot(k_M)

    return ns_Out


def processBrepObject(objref, sDir, iQtyToAddInEachSpan=1, bReplace=True, bEcho=True, bDebug=False):
    """
    objref: rd.ObjRef of BrepFace.
    """

    face = rs.coercesurface(objref, raise_if_missing=False)
    if isinstance(face, rg.BrepFace):
        srf = face.UnderlyingSurface()
    else:
        print "{} skipped.".format(rs.coercerhinoobject(objref).ObjectType)
        return

    ns_Res = createSurface(
        srf.ToNurbsSurface(),
        sDir,
        iQtyToAddInEachSpan,
        bDebug)

    if ns_Res is None:
        return

    if bReplace:
        gB_In = rs.coerceguid(objref)


        if face.Brep.Faces.Count == 1:
            if sc.doc.Objects.Replace(
                objectId=gB_In,
                surface=ns_Res
            ):
                return gB_In

        # Polyfaced brep.
        face.Brep.AddSurface(ns_Res)
        face.ChangeSurface(surfaceIndex=face.Brep.Surfaces.Count-1)
        face.Brep.Compact()
        if sc.doc.Objects.Replace(
            objectId=gB_In,
            brep=face.Brep
        ):
            return gB_In


    # Add new face.
    gB_Out = sc.doc.Objects.AddSurface(ns_Res)
    if gB_Out != gB_Out.Empty:
        return gB_Out
    else:
        print "Could not add surface."


def main():

    rc = getInput()
    if rc is None: return

    (
        objref,
        sDir,
        iQtyToAddInEachSpan,
        bReplace,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gBs_Res = processBrepObject(
        objref,
        sDir,
        iQtyToAddInEachSpan,
        bReplace=bReplace,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()