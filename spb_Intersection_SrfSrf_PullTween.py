"""
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
210921: Created.
210927: Bug fix: Match starting curve directions.
210929: Modified tolerance for point-point coincidence.
        Replaced CreateTweenCurvesWithMatching with CreateTweenCurvesWithSampling for better results.
210930: Now numSamples of CreateTweenCurvesWithSampling is incremented until 3 deviation results are within tolerance.
221004: Now, when start curve's endpoints cannot automatically be determined, user can pick them.
"""

import Rhino
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


    key = 'bUseUnderlyingSrfs'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTolerance'; keys.append(key)
    values[key] = sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fMinCrvLength'; keys.append(key)
    values[key] = 20.0 * sc.doc.ModelAbsoluteTolerance
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

        if key == 'fTolerance':
            if cls.riOpts[key].CurrentValue <= 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < Rhino.RhinoMath.ZeroTolerance:
                cls.riOpts[key].CurrentValue = Rhino.RhinoMath.ZeroTolerance

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Select 2 BrepFaces with options.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 faces")

    go.GeometryFilter = go.GeometryFilter.Surface

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('bUseUnderlyingSrfs')
        addOption('fTolerance')
        addOption('fMinCrvLength')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()

            return (
                objrefs,
                Opts.values['bUseUnderlyingSrfs'],
                Opts.values['fTolerance'],
                Opts.values['fMinCrvLength'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fTolerance'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getStartCrvEndPts(face_A, face_B, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    fMinCrvLength = getOpt('fMinCrvLength')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    cs_Es_A = []
    for iE in face_A.AdjacentEdges():
        cs_Es_A.append(face_A.Brep.Edges[iE].ToNurbsCurve())

    cs_Es_B = []
    for iE in face_B.AdjacentEdges():
        cs_Es_B.append(face_B.Brep.Edges[iE].ToNurbsCurve())


    # Determine tolerance for whether EdgeFace points are coincident to EdgeEdge points.
    fMaxEdgeTol_A = max([face_A.Brep.Edges[iE].Tolerance for iE in face_A.AdjacentEdges()])
    fMaxEdgeTol_B = max([face_B.Brep.Edges[iE].Tolerance for iE in face_B.AdjacentEdges()])
    fMaxEdgeTol = max((fMaxEdgeTol_A, fMaxEdgeTol_B))
    fPoint_coincident_Tol = max(fMaxEdgeTol, 10.0*sc.doc.ModelAbsoluteTolerance)


    pts_Intersects = []


    def addNewPointToIntersectList(pt):
        """
        Returns:
            True if point was added.
            False if point was not added, e.i., already in list.
            None if list already has 2 points.
        """
        if len(pts_Intersects) == 0:
            pts_Intersects.append(pt)
            return True
        for ptSaved in pts_Intersects:
            dist = pt.DistanceTo(ptSaved)
            #print(dist)
            if dist <= fPoint_coincident_Tol:
                return False
        else:
            # Match not found.
            pts_Intersects.append(pt)
            return True


    def getEdgeEdgeIntersections():
        for cA in cs_Es_A:
            for cB in cs_Es_B:
                crvinters = rg.Intersect.Intersection.CurveCurve(
                    cA,
                    cB,
                    tolerance=fTolerance,
                    overlapTolerance=0.0)
                if crvinters.Count == 0: continue # to next curve.
                for crvinter in crvinters:
                    pt = crvinter.PointA
                    #sc.doc.Objects.AddPoint(pt)
                    addNewPointToIntersectList(pt)

    getEdgeEdgeIntersections()

    #    if bDebug:
    #        for pt in pts_Intersects: sc.doc.Objects.AddPoint(pt)
    #        return

    if len(pts_Intersects) > 2:
        if bDebug:
            for pt in pts_Intersects: sc.doc.Objects.AddPoint(pt)
        return None, "More than 2 points found."


    def addEdgeFaceIntersections():
        for cs, f in zip((cs_Es_A, cs_Es_B), (face_B, face_A)):

            for c in cs:
                rc = rg.Intersect.Intersection.CurveBrepFace(
                        c,
                        f,
                        tolerance=fTolerance)

                if not rc[0]: continue

                pts = list(rc[2])
                if len(pts) == 0:
                    continue

                for pt in pts:
                    #sc.doc.Objects.AddPoint(pt)
                    addNewPointToIntersectList(pt)

    if bDebug:
        for pt in pts_Intersects:
            sc.doc.Objects.AddPoint(pt)

    if len(pts_Intersects) < 2:
        addEdgeFaceIntersections()


    if len(pts_Intersects) == 0:
        return None, "No intersection points found."
    if len(pts_Intersects) == 1:
        return None, "Only 1 intersection point found."
    if len(pts_Intersects) > 2:
        if bDebug:
            for pt in pts_Intersects: sc.doc.Objects.AddPoint(pt)
        return None, "More than 2 points found."


    pt_Start = pts_Intersects[0]
    pt_End = pts_Intersects[1]

    return (pt_Start, pt_End), None


def createCurve(face_A, face_B, pt_Start, pt_End, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    fMinCrvLength = getOpt('fMinCrvLength')
    bDebug = getOpt('bDebug')


    def doCurveStucturesMatch(a, b):
        if not isinstance(a, rg.NurbsCurve):
            return
        if not isinstance(b, rg.NurbsCurve):
            return

        if a.Degree != b.Degree:
            return False

        if a.Points.Count != b.Points.Count:
            return False

        if a.Knots.Count != b.Knots.Count:
            return False

        if a.SpanCount != b.SpanCount:
            return False

        #c = a.DuplicateCurve()
        #d = b.DuplicateCurve()

        #c.Domain = rg.Interval(0.0, 1.0)
        #d.Domain = rg.Interval(0.0, 1.0)

        iK = a.Degree

        while iK < (a.Knots.Count - a.Degree):
            sc.escape_test()

            mA = a.Knots.KnotMultiplicity(iK)

            if mA != b.Knots.KnotMultiplicity(iK):
                return False

            iK += mA

        return True


    def createTweenCurves(a, b, i, iAdditionalSamples=0):
        if doCurveStucturesMatch(a, b):
            return rg.Curve.CreateTweenCurves(
                a,
                b,
                numCurves=1,
                tolerance=0.5*fTolerance)

        # CreateTweenCurvesWithMatching refits the starting curves but doesn't
        # seem to use tolerance in doing so.
        if i == 1:
            return rg.Curve.CreateTweenCurvesWithMatching(
                a,
                b,
                numCurves=1,
                tolerance=0.5*fTolerance)

        numSamples = i
        #numSamples = max((a.SpanCount, b.SpanCount)) + iAdditionalSamples
        #numSamples = int(a.GetLength() / (100.0*fTolerance)) + 1

        return rg.Curve.CreateTweenCurvesWithSampling(
            a,
            b,
            numCurves=1,
            numSamples=numSamples,
            tolerance=0.5*fTolerance)


    def getEdgeWithPointsAtBothVertices(face):
        for iE in face.AdjacentEdges():
            e = face.Brep.Edges[iE]
            if (
                (pt_Start.DistanceTo(e.StartVertex.Location) <= fTolerance
                and
                pt_End.DistanceTo(e.EndVertex.Location) <= fTolerance)
                or
                (pt_End.DistanceTo(e.StartVertex.Location) <= fTolerance
                and
                pt_Start.DistanceTo(e.EndVertex.Location) <= fTolerance)
            ):
                return e.DuplicateCurve()


    crv_WIP = None
    iAdditionalSamples = 0

    cA_Start = getEdgeWithPointsAtBothVertices(face_A)
    if cA_Start is not None:
        cB_Start = getEdgeWithPointsAtBothVertices(face_B)
        if cB_Start is not None:
            if not rg.Curve.DoDirectionsMatch(cA_Start, cB_Start):
                # Align curves to each other and per pt_Start.
                if (
                    cA_Start.PointAtStart.DistanceTo(pt_Start) <
                    cB_Start.PointAtStart.DistanceTo(pt_Start)
                ):
                    cB_Start.Reverse()
                else:
                    cA_Start.Reverse()
            rc = createTweenCurves(cA_Start, cB_Start, i=1)
            if len(rc) == 1:
                crv_WIP = rc[0]
                i = 2 # will be the next iteration.


    if crv_WIP is None:
        crv_WIP = rg.LineCurve(pt_Start, pt_End)
        i = 1 # will be the next iteration.

    if bDebug:
        sc.doc.Objects.AddCurve(crv_WIP)#; sc.doc.Views.Redraw(); return


    crv_WIP_Prev = None
    dev = None
    dev_Prev = None

    while True:
        if sc.escape_test(False, False):
            return None, "Escape key pressed.  Last deviation was {:.{}f}.".format(
                dev, sc.doc.ModelDistanceDisplayPrecision)
            #sc.doc.Objects.AddCurve(crv_WIP)
            #sc.doc.Views.Redraw()
            return

        crv_WIP.SetStartPoint(pt_Start)
        crv_WIP.SetEndPoint(pt_End)

        rc = rg.Curve.PullToBrepFace(crv_WIP, face_A, tolerance=fTolerance)
        if len(rc) == 0:
            if bDebug:
                print("No result from pull onto face A.  Adjusting end points and retrying PullToBrepFace.")
            crv_WIP.SetStartPoint(pt_Start)
            crv_WIP.SetEndPoint(pt_End)
            rc = rg.Curve.PullToBrepFace(crv_WIP, face_A, tolerance=fTolerance)
            if len(rc) == 0:
                return None, "No result from pull onto face A."
        if len(rc) == 1:
            crv_Pull_A = rc[0]
        else:
            lengths = [c.GetLength() for c in rc]
            crv_Pull_A = rc[lengths.index(max(lengths))]
            #if bDebug:
            #    for c in rc:
            #        sc.doc.Objects.AddCurve(c)
            #return None, "More than one curve from pull onto face A."

        #if i == 4: sc.doc.Objects.AddCurve(crv_Pull_A)

        if crv_Pull_A.GetLength() < fMinCrvLength:
            return None, "Curve is too short."

        rc = rg.Curve.PullToBrepFace(crv_WIP, face_B, tolerance=fTolerance)
        if len(rc) == 0:
            if bDebug:
                print("No result from pull onto face B.  Adjusting end points and retrying PullToBrepFace.")
            crv_WIP.SetStartPoint(pt_Start)
            crv_WIP.SetEndPoint(pt_End)
            rc = rg.Curve.PullToBrepFace(crv_WIP, face_B, tolerance=fTolerance)
            if len(rc) == 0:
                return None, "No result from pull onto face B."
        if len(rc) == 1:
            crv_Pull_B = rc[0]
        else:
            lengths = [c.GetLength() for c in rc]
            crv_Pull_B = rc[lengths.index(max(lengths))]
            #if bDebug:
            #    for c in rc:
            #        sc.doc.Objects.AddCurve(c)
            #return None, "More than one curve from pull onto face B."

        #if i == 4: sc.doc.Objects.AddCurve(crv_Pull_B)

        if crv_Pull_B.GetLength() < fMinCrvLength:
            return None, "Curve is too short."

        if i == 4:
            pass

        rc = createTweenCurves(crv_Pull_A, crv_Pull_B, i, iAdditionalSamples)
        if len(rc) != 1:
            return None, "{} curves returned from CreateTweenCurves...".format(len(rc))
        crv_WIP = rc[0]
        #if bDebug:
        #    sc.doc.Objects.AddCurve(crv_WIP)
        if bDebug: sEval = "crv_WIP.SpanCount"; print(sEval+':', eval(sEval))

        if crv_WIP_Prev is None:
            crv_WIP_Prev = crv_WIP
            continue

        rc = rg.Curve.GetDistancesBetweenCurves(crv_WIP, crv_WIP_Prev, tolerance=0.1*fTolerance)
        if not rc[0]:
            raise ValueError("GetDistancesBetweenCurves could not be calculated.")

        if dev is not None:
            dev_Prev = dev

        dev = rc[1]

        if dev <= fTolerance:
            #sc.doc.Objects.AddCurve(crv_WIP_Prev)
            #sc.doc.Objects.AddCurve(crv_WIP)
            #1/0
            
            #dev_FromPrev = dev
            
            # Further checks.
            rc = rg.Curve.GetDistancesBetweenCurves(crv_WIP, crv_Pull_A, tolerance=0.1*fTolerance)
            if not rc[0]:
                raise ValueError("GetDistancesBetweenCurves could not be calculated.")
            
            dev = rc[1]
            
            if dev > 0.5*fTolerance:
                iAdditionalSamples += 1
                i += 1
                continue
            
            rc = rg.Curve.GetDistancesBetweenCurves(crv_WIP, crv_Pull_B, tolerance=0.1*fTolerance)
            if not rc[0]:
                raise ValueError("GetDistancesBetweenCurves could not be calculated.")
            
            dev = rc[1]
            
            if dev > 0.5*fTolerance:
                iAdditionalSamples += 1
                i += 1
                continue
            
            # Success.
            
            crv_WIP_Prev.Dispose()
            return crv_WIP, "Intersection found after {} Pull/TweenCurves iterations.".format(i)

        if dev_Prev is not None and abs(dev-dev_Prev) <= 1e-9:
            return (
                crv_WIP,
                "Deviation is still {:.{}f} after {} Pull/TweenCurves iterations." \
                    "  Check results.".format(
                        dev, sc.doc.ModelDistanceDisplayPrecision, i)
                )

        #print(dev)

        crv_WIP_Prev.Dispose()
        crv_WIP_Prev = crv_WIP

        i += 1


def createCurveObject(face_A, face_B, pt_Start, pt_End, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTolerance = getOpt('fTolerance')
    fMinCrvLength = getOpt('fMinCrvLength')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    rc = createCurve(
        face_A,
        face_B,
        pt_Start,
        pt_End,
        fTolerance=fTolerance,
        fMinCrvLength=fMinCrvLength,
        bDebug=bDebug,
        )

    if rc is None: return

    c_res, sLog = rc

    if bEcho:
        if c_res is None:
            print(sLog + "  No intersection was returned.")
        else:
            print(sLog)

    if c_res is None: return

    g = sc.doc.Objects.AddCurve(c_res)
    if g == g.Empty: return
    return g


def main():

    rc = getInput()
    if rc is None: return

    (
        objrefs,
        bUseUnderlyingSrfs,
        fTolerance,
        fMinCrvLength,
        bEcho,
        bDebug,
        ) = rc


    face_A = objrefs[0].Face()
    face_B = objrefs[1].Face()

    if bUseUnderlyingSrfs:
        face_A = face_A.UnderlyingSurface().ToBrep().Faces[0]
        face_B = face_B.UnderlyingSurface().ToBrep().Faces[0]


    rc = getStartCrvEndPts(
        face_A,
        face_B,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rc is None: return

    if rc[0] is None:
        print(rc[1], "Pick start curve's end points to continue.")
        res, pt_Start = ri.RhinoGet.GetPoint(
            "Pick one end point of intersection",
            acceptNothing=False)
        if res != Rhino.Commands.Result.Success: return

        res, pt_End = ri.RhinoGet.GetPoint(
            "Pick other end point of intersection",
            acceptNothing=False)
        if res != Rhino.Commands.Result.Success: return
    else:
        pt_Start, pt_End = rc[0]


    gCrv_Res = createCurveObject(
        face_A,
        face_B,
        pt_Start,
        pt_End,
        fTolerance=fTolerance,
        fMinCrvLength=fMinCrvLength,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if gCrv_Res is not None:
        sc.doc.Objects.UnselectAll()
        sc.doc.Objects.Select(objectId=gCrv_Res)
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()