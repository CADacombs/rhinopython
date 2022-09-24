"""
This script is an alternative to _ConnectSrf, except:
    The trimmed faces are not shrunk.
    Sometimes will produce results when _ConnectSrf fails.
    Sometimes produces more accurate results.

This script uses Brep.CreateFilletSurface with radius=0.0, trim=True, and extend=True
to trim 2 surfaces with each other.  extend=True allows for the surfaces to not
have to match at their sides.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
200924: Created.
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
    riAddOpts = {}
    stickyKeys = {}


    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = "Tol"
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bTryOtherTolsOnFail'; keys.append(key)
    values[key] = False
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

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fRadius':
            if cls.riOpts[key].CurrentValue < 1e-6:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key == 'fTolerance':
            if cls.riOpts[key].CurrentValue < 1e-6:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        if key in cls.stickyKeys:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_FaceA():
    """
    Select BrepFace with options.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick 1st surface at a point on area to keep")

    go.GeometryFilter = rd.ObjectType.Surface
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.TopSurface

    idxs_Opt = {}

    while True:

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fTolerance')
        addOption('bTryOtherTolsOnFail')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()

            sc.doc.Objects.UnselectAll()

            return objref

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_FaceB(objref_A):
    """
    Select BrepFaces with options.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Pick 2nd surface at a point on area to keep")

    go.GeometryFilter = rd.ObjectType.Surface
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.TopSurface

    rdObj_A = objref_A.Object()
    compIdx_A = objref_A.GeometryComponentIndex
    rdObj_A.HighlightSubObject(compIdx_A, highlight=True)
    sc.doc.Views.Redraw()


    idxs_Opt = {}

    while True:

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fTolerance')
        addOption('bTryOtherTolsOnFail')
        addOption('bEcho')
        addOption('bDebug')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            rdObj_A.HighlightSubObject(compIdx_A, highlight=False)
            sc.doc.Views.Redraw()
            return

        if res == ri.GetResult.Object:
            objref_B = go.Object(0)

            rdObj_B = objref_B.Object()
            compIdx_B = objref_B.GeometryComponentIndex

            sc.doc.Objects.UnselectAll()

            if (
                rdObj_B.Id == rdObj_A.Id and
                compIdx_B.ComponentIndexType == compIdx_A.ComponentIndexType and
                compIdx_B.Index == compIdx_A.Index
            ):
                # Face A was picked.
                rdObj_A.HighlightSubObject(compIdx_A, highlight=True)
                sc.doc.Views.Redraw()
                continue

            go.Dispose()

            rdObj_A.HighlightSubObject(compIdx_A, highlight=False)
            sc.doc.Views.Redraw()

            return objref_B

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def processBreps(rgFaceA, pt3d_A, rgFaceB, pt3d_B, fTolerance=None, bTryOtherTolsOnFail=False, bEcho=True, bDebug=False):
    """
    """

    rgBs_PerCFS = []
    rgBs_Face0_Out = []
    rgBs_Face1_Out = []

    # The following are in-sync with rgBs_PerCFS.
    fTols_Out = []
    areas = []
    boundingBoxCenters = []

    b, u, v = rgFaceA.ClosestPoint(pt3d_A)
    pt2d_uvA = rg.Point2d(u, v)
    b, u, v = rgFaceB.ClosestPoint(pt3d_B)
    pt2d_uvB = rg.Point2d(u, v)


    fTol_Nom = sc.doc.ModelAbsoluteTolerance if fTolerance is None else fTolerance

    fTols = [fTol_Nom]
    if bTryOtherTolsOnFail:
        fTols.extend([0.1*fTol_Nom, 10.0*fTol_Nom])


    for iT, fTol in enumerate(fTols):

        rgBs_Fillet, rgBs_Face0, rgBs_Face1 = rg.Brep.CreateFilletSurface(
                face0=rgFaceA,
                uv0=pt2d_uvA,
                face1=rgFaceB,
                uv1=pt2d_uvB,
                radius=0.0,
                trim=True,
                extend=True,
                tolerance=fTol
                )

        if not (rgBs_Face0 and rgBs_Face1):
            if bEcho:
                print("Faces were not trimmed at {} tolerance.".format(
                    fTol))
            continue

        return rgBs_Face0, rgBs_Face1, fTol


def processBrepObjects(objref_A, objref_B, fTolerance=None, bTryOtherTolsOnFail=False, bEcho=True, bDebug=False):
    """
    """

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    fA = objref_A.Face()
    fB = objref_B.Face()

    ptA = objref_A.SelectionPoint()
    ptB = objref_B.SelectionPoint()

    gBreps_Fillets = []


    rc = processBreps(
        fA,
        ptA,
        fB,
        ptB,
        fTolerance=fTolerance,
        bTryOtherTolsOnFail=bTryOtherTolsOnFail,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if not rc:
        if bEcho: print("Fillets were not created using provided points.")
        return

    objrefs = [objref_A, objref_B]
    rgBs_Faces_Res = [None, None]

    rgBs_Faces_Res[0], rgBs_Faces_Res[1], fTols_Ret = rc

    bReplaced = [False, False]

    for i in 0,1:
        if len(rgBs_Faces_Res[i]) == 1:
            if sc.doc.Objects.Replace(objrefs[i], rgBs_Faces_Res[i][0]):
                bReplaced[i] = True
        elif len(rgBs_Faces_Res[i]) > 1:
            print("Face {} was split into mulitple faces.".format(i+1))
            gBs_F0 = []
            for rgB in rgBs_Faces_Res[i]:
                gBs_F0.append(sc.doc.Objects.AddBrep(rgB))
            if not any(gB==gB.Empty for gB in gBs_F0):
                sc.doc.Objects.Delete(objrefs[i], quiet=False)

    if all(bReplaced):
        print("Replaced both faces.")
    elif not any(bReplaced):
        print("Neither face was replaced.")
    else:
        print("Replaced {} face.".format("1st" if bReplaced[0] else "2nd"))

    return True


def main():

    if Rhino.RhinoApp.ExeVersion == 5:
        print("This script uses Brep.CreateFilletSurface that is missing from Rhino V5's RC.")
        return

    objref_A = getInput_FaceA()
    if objref_A is None: return

    fTolerance = Opts.values['fTolerance']
    bTryOtherTolsOnFail = Opts.values['bTryOtherTolsOnFail']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    objref_B = getInput_FaceB(objref_A)
    if objref_B is None: return


    fTolerance = Opts.values['fTolerance']
    bTryOtherTolsOnFail = Opts.values['bTryOtherTolsOnFail']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    gBreps_Fillets = processBrepObjects(
        objref_A,
        objref_B,
        fTolerance=fTolerance,
        bTryOtherTolsOnFail=bTryOtherTolsOnFail,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()