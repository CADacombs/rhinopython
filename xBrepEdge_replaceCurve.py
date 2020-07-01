"""
200701: Created.
"""

import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    idxOpt = {}
    stickyKeys = {}


    key = 'fTol'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
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
            if key == 'fTol':
                if cls.riOpts[key].CurrentValue < 0.0:
                    cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_Edge():
    """
    Get BrepEdge with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edge to replace")
    
    go.GeometryFilter = rd.ObjectType.EdgeFilter

    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opts = {}
    
    while True:
        Opts.addOption(go, 'fTol')
        Opts.addOption(go, 'bEcho')
        Opts.addOption(go, 'bDebug')
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return

        # An option was selected or a number was entered.

        if res == ri.GetResult.Number:
            Opts.riOpts['fTol'].CurrentValue = go.Number()

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInput_Curve():
    """
    Get Curve with optional input.
    """

    sc.doc.Objects.UnselectAll()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select replacing curve")
    
    go.GeometryFilter = rd.ObjectType.Curve

    go.AcceptNumber(enable=True, acceptZero=True)

    idxs_Opts = {}
    
    while True:
        Opts.addOption(go, 'fTol')
        Opts.addOption(go, 'bEcho')
        Opts.addOption(go, 'bDebug')
        
        res = go.Get()
        
        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])
        elif res == ri.GetResult.Cancel:
            go.Dispose()
            return

        # An option was selected or a number was entered.

        if res == ri.GetResult.Number:
            Opts.riOpts['fTol'].CurrentValue = go.Number()

        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def replaceBrepObject(objref_Edge, objref_Curve, **kwargs):
    """
    """


    def setOpt(key, value=None):
        if key in kwargs:
            return kwargs[key]
        elif key in Opts.riOpts:
            return Opts.riOpts[key].InitialValue
        else:
            return value

    fTol = setOpt('fTol')
    bEcho = setOpt('bEcho')
    bDebug = setOpt('bDebug')


    def printValidLog():
        bValid, sLog = rgBrep_Out.IsValidWithLog()
        if bValid:
            print "Brep is valid."
        else:
            print sLog


    rdBrep_In = objref_Edge.Object()
    rgBrep_In = rdBrep_In.Geometry.Duplicate() # So it is not DocumentControlled.
    rgEdge_toMod = objref_Edge.Edge()
    rgCrv_In = objref_Curve.Curve().DuplicateCurve() # DuplicateCurve in case Curve is a BrepEdge.

    rgBrep_In.Compact() # Purge unused components to avoid duplicating them.


    rgBrep_Out = rg.Brep()


    for srf_In in rgBrep_In.Surfaces:
        srf_Out = srf_In.Duplicate()
        rgBrep_Out.AddSurface(srf_Out)


    for iEC, iE in sorted([(e.EdgeCurveIndex, i) for i, e in enumerate(rgBrep_In.Edges)]):
        edge_In = rgBrep_In.Edges[iE]
        if edge_In.EdgeIndex != rgEdge_toMod.EdgeIndex:
            edgeCrv_Out = edge_In.EdgeCurve.Duplicate()
            rgBrep_Out.AddEdgeCurve(edgeCrv_Out)
        else:
            if not rg.Curve.DoDirectionsMatch(curveA=edge_In.EdgeCurve, curveB=rgCrv_In):
                rgCrv_In.Reverse()
            rgBrep_Out.AddEdgeCurve(rgCrv_In)


    for iTC, iT in sorted([(t.TrimCurveIndex, i) for i, t in enumerate(rgBrep_In.Trims)]):
        trim_In = rgBrep_In.Trims[iT]

        if trim_In.Edge.EdgeIndex != rgEdge_toMod.EdgeIndex:
            trimCrv_Out = trim_In.TrimCurve.Duplicate()
            rgBrep_Out.AddTrimCurve(trimCrv_Out)
        else:
            face_In = trim_In.Face
            srf_In = face_In.UnderlyingSurface()
            trimCrv_Out = srf_In.Pullback(rgCrv_In, tolerance=fTol)
            if not rg.Curve.DoDirectionsMatch(curveA=trimCrv_Out, curveB=trim_In.TrimCurve):
                trimCrv_Out.Reverse()
            rgBrep_Out.AddTrimCurve(trimCrv_Out)


    #for i, c in enumerate(rgBrep_Out.Curves3D):
    #    sc.doc.Objects.AddCurve(c)
    #    sc.doc.Objects.AddTextDot(text=str(i),location=c.PointAt(c.Domain.Mid))
    #sc.doc.Views.Redraw()
    #return


    for iF, face_In in enumerate(rgBrep_In.Faces):
        rgBrep_Out.Faces.Add(face_In.SurfaceIndex)
        rgBrep_Out.Faces[rgBrep_Out.Faces.Count-1].OrientationIsReversed = (
            face_In.OrientationIsReversed)


    for iV, vertex_In in enumerate(rgBrep_In.Vertices):
        rgBrep_Out.Vertices.Add(
            point=vertex_In.Location,
            vertexTolerance=fTol)


    for iL, loop_In in enumerate(rgBrep_In.Loops):
        rgBrep_Out.Loops.Add(
            loopType=loop_In.LoopType,
            face=rgBrep_Out.Faces[loop_In.Face.FaceIndex])


    for iE, edge_In in enumerate(rgBrep_In.Edges):
        rgBrep_Out.Edges.Add(
            startVertexIndex=edge_In.StartVertex.VertexIndex,
            endVertexIndex=edge_In.EndVertex.VertexIndex,
            curve3dIndex=edge_In.EdgeCurveIndex,
            edgeTolerance=fTol)


    for iT, trim_In in enumerate(rgBrep_In.Trims):
        rgBrep_Out.Trims.Add(
            edge=rgBrep_Out.Edges[trim_In.Edge.EdgeIndex],
            rev3d=trim_In.IsReversed(),
            loop=rgBrep_Out.Loops[trim_In.Loop.LoopIndex],
            curve2dIndex=trim_In.TrimCurveIndex)
        rgBrep_Out.Trims[rgBrep_Out.Trims.Count-1].SetTolerances(0.0,0.0)
        rgBrep_Out.Trims[rgBrep_Out.Trims.Count-1].IsoStatus = (
            trim_In.Face.IsIsoparametric(trim_In.TrimCurve))


    if rgBrep_Out.IsValid:
        rgBrep_Out.SetTolerancesBoxesAndFlags(
            bLazy=False,
            bSetVertexTolerances=True,
            bSetEdgeTolerances=True,
            bSetTrimTolerances=True,
            bSetTrimIsoFlags=True,
            bSetTrimTypeFlags=True,
            bSetLoopTypeFlags=True,
            bSetTrimBoxes=True)
        if sc.doc.Objects.Replace(objectId=rdBrep_In.Id, brep=rgBrep_Out):
            if bEcho:
                print "Brep geometry was replaced."
            return rdBrep_In.Id
        else:
            if bEcho:
                print "Brep geometry was not replaced."
    else:
        if bEcho:
            printValidLog()


def main():
    
    rc = getInput_Edge()
    if rc is None: return

    (
        objref_Edge,
        fTol,
        bEcho,
        bDebug,
        ) = rc

    rc = getInput_Curve()
    if rc is None: return

    (
        objref_Curve,
        fTol,
        bEcho,
        bDebug,
        ) = rc

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()

    gBrep_Res = replaceBrepObject(
        objref_Edge=objref_Edge,
        objref_Curve=objref_Curve,
        fTol=fTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()