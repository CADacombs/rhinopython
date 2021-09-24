"""
180922-24: Created.
190626, 1206, 200215: Modified an option default value.
210218: Refactored getInput.  Removed the single-span, degree-3 bezier option.
210224: Added an option.
210327: Bug fix.
210411: Added bTryOtherTols.  Moved some of createBreps to the new createBrepObjects.
210418: Added more printed feedback.  Modified an option default value.
210923: Disabled option for bRebuildDeg5 since its code hasn't yet been implemented.
210924: Added one-to-many selection option.

_FilletSrf and Brep.CreateFilletSurface create surfaces that are degree 2 rational in one direction.


TODO: Add code for bRebuildDeg5.
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


    key = 'fRadius'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bMultiPickB'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    # No sticky for this since it should always default to False.

    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = "Tol"
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bTryOtherTols'; keys.append(key)
    values[key] = False
    #names[key] = "IncludeUniquePlusMinusTols"
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bRebuildDeg5'; keys.append(key)
    values[key] = False
    names[key] = 'AcrossSrf'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Deg2Rat', 'Deg5NonRat')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtend'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTrimFillets'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bUseUnderlyingSrf'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bTryOtherPickPts'; keys.append(key)
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

    go.SetCommandPrompt("Select face A")

    go.GeometryFilter = rd.ObjectType.Surface

    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fRadius')
        addOption('bMultiPickB')
        addOption('fTolerance')
        addOption('bTryOtherTols')
        #addOption('bRebuildDeg5')
        addOption('bExtend')
        addOption('bTrimFillets')
        addOption('bUseUnderlyingSrf')
        addOption('bTryOtherPickPts')
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

            return (
                objref,
                Opts.values['bMultiPickB'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

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

    go.SetCommandPrompt("Select face B")

    go.GeometryFilter = rd.ObjectType.Surface

    rdObj_A = objref_A.Object()
    compIdx_A = objref_A.GeometryComponentIndex
    rdObj_A.HighlightSubObject(compIdx_A, highlight=True)
    sc.doc.Views.Redraw()


    go.AcceptNumber(True, acceptZero=True)

    idxs_Opt = {}

    while True:

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fRadius')
        addOption('fTolerance')
        addOption('bTryOtherTols')
        #addOption('bRebuildDeg5')
        addOption('bExtend')
        addOption('bTrimFillets')
        addOption('bUseUnderlyingSrf')
        addOption('bTryOtherPickPts')
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
                rdObj_A.HighlightSubObject(compIdx_A, highlight=True)
                sc.doc.Views.Redraw()
                continue

            go.Dispose()

            rdObj_A.HighlightSubObject(compIdx_A, highlight=False)
            sc.doc.Views.Redraw()

            return (
                objref_B,
                Opts.values['fRadius'],
                Opts.values['fTolerance'],
                Opts.values['bTryOtherTols'],
                Opts.values['bExtend'],
                Opts.values['bTrimFillets'],
                Opts.values['bUseUnderlyingSrf'],
                Opts.values['bTryOtherPickPts'],
                Opts.values['bEcho'],
                Opts.values['bDebug'],
                )

        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createBreps(rgFaceA, pt3d_A, rgFaceB, pt3d_B, fRadius, fTolerance=None, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bTryOtherTols = getOpt('bTryOtherTols')
    #bRebuildDeg5 = getOpt('bRebuildDeg5')
    bExtend = getOpt('bExtend')
    bTrimFillets = getOpt('bTrimFillets')
    bUseUnderlyingSrf = getOpt('bUseUnderlyingSrf')
    bTryOtherPickPts = getOpt('bTryOtherPickPts')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if bUseUnderlyingSrf:
        rgFaceA = rgFaceA.UnderlyingSurface().ToBrep().Faces[0] # UnderlyingSurface() is not needed.
        rgFaceB = rgFaceB.UnderlyingSurface().ToBrep().Faces[0]
    else:
        # No change in rgFaceA or rgFaceB.
        pass
    
    rgBs_PerCFS = []
    # The following are in-sync with rgBs_PerCFS.
    fTols_Out = []
    areas = []
    boundingBoxCenters = []

    b, u, v = rgFaceA.ClosestPoint(pt3d_A)
    pt2d_uvA = rg.Point2d(u, v)
    b, u, v = rgFaceB.ClosestPoint(pt3d_B)
    pt2d_uvB = rg.Point2d(u, v)


    fTol_Nom = sc.doc.ModelAbsoluteTolerance if fTolerance is None else fTolerance


    if bTryOtherTols:
        fTols = fTol_Nom, 10.0*fTol_Nom, 0.1*fTol_Nom
    else:
        fTols = fTol_Nom,


    for iT, fTol in enumerate(fTols):

        rgBs_Fillet = rg.Brep.CreateFilletSurface(
                face0=rgFaceA,
                uv0=pt2d_uvA,
                face1=rgFaceB,
                uv1=pt2d_uvB,
                radius=fRadius,
                extend=bExtend,
                tolerance=fTol
                )
        if not rgBs_Fillet:
            if bEcho:
                print "Fillet was not created at {} tolerance using provided points.".format(
                    fTol)
            continue


        def getTotalArea():
            total_area = 0.0
            for brep in rgBs_Fillet:
                a = brep.GetArea()
                if not a:
                    return None
                total_area += a
            return total_area

        area_ThisTol = getTotalArea()

        if area_ThisTol is None:
            print "Area could not be obtained."
            rgBs_PerCFS.append(rgBs_Fillet)
            fTols_Out.append(fTol)
            areas.append(None)
            boundingBoxCenters.append(None)
            continue # to next tolerance.


        def getBBoxCenter():
            bbox_ThisTol = rg.BoundingBox.Empty
            for rgB in rgBs_Fillet:
                bbox = rgB.GetBoundingBox(accurate=True)
                if not bbox:
                    return None
                bbox_ThisTol.Union(bbox)
            center = bbox_ThisTol.Center
            return center

        center_ThisTol = getBBoxCenter()


        if center_ThisTol is None:
            print "Bounding box could not be obtained."
            rgBs_PerCFS.append(rgBs_Fillet)
            fTols_Out.append(fTol)
            areas.append(area_ThisTol)
            boundingBoxCenters.append(None)
            continue # to next tolerance.


        if iT == 0:
            rgBs_PerCFS.append(rgBs_Fillet)
            fTols_Out.append(fTol)
            areas.append(area_ThisTol)
            boundingBoxCenters.append(center_ThisTol)
        else:
            if bDebug:
                print "Fillet with accumulative area of" \
                    " {} square units created at {} tolerance.".format(
                        area_ThisTol, fTol)


            def matchFound():
                for i in range(len(fTols_Out)):
                    if abs(area_ThisTol - areas[i]) <= fTol_Nom**2:
                        if center_ThisTol.DistanceTo(boundingBoxCenters[i]) <= fTol_Nom:
                            return True
                return False

            if not matchFound():
                rgBs_PerCFS.append(rgBs_Fillet)
                fTols_Out.append(fTol)
                areas.append(area_ThisTol)
                boundingBoxCenters.append(center_ThisTol)


        if not bTryOtherPickPts:
            continue

        # Create fillets based on using various locations of each surface.
        
        for pt2d_uvA in (
                rg.Point2d(rgFaceA.Domain(0).Mid, rgFaceA.Domain(1).T0),
                rg.Point2d(rgFaceA.Domain(0).T1, rgFaceA.Domain(1).Mid),
                rg.Point2d(rgFaceA.Domain(0).Mid, rgFaceA.Domain(1).T1),
                rg.Point2d(rgFaceA.Domain(0).T0, rgFaceA.Domain(1).Mid),
        ):
            for pt2d_uvB in (
                    rg.Point2d(rgFaceB.Domain(0).Mid, rgFaceB.Domain(1).T0),
                    rg.Point2d(rgFaceB.Domain(0).T1, rgFaceB.Domain(1).Mid),
                    rg.Point2d(rgFaceB.Domain(0).Mid, rgFaceB.Domain(1).T1),
                    rg.Point2d(rgFaceB.Domain(0).T0, rgFaceB.Domain(1).Mid),
            ):
                #print pt2d_uvA, pt2d_uvB
                rgBs_Fillet = rg.Brep.CreateFilletSurface(
                        face0=rgFaceA,
                        uv0=pt2d_uvA,
                        face1=rgFaceB,
                        uv1=pt2d_uvB,
                        radius=fRadius,
                        extend=bExtend,
                        tolerance=fTol
                        )
                if not rgBs_Fillet:
                    continue

                #sc.doc.Objects.AddPoint(rgFaceA.PointAt(pt2d_uvA.X, pt2d_uvA.Y))
                #sc.doc.Objects.AddPoint(rgFaceB.PointAt(pt2d_uvB.X, pt2d_uvB.Y))


                # Don't include any area matches.
                for rgB in rgBs_Fillet:


                    area = rgB.GetArea()
                    if not area:
                        print "Area could not be obtained."
                        continue

                    bbox = rgB.GetBoundingBox(accurate=False)
                    if not bbox:
                        print "Bounding box could not be obtained."
                        continue
                    center = bbox.Center


                    def isFilletInList():
                        for a in areas:
                            if abs(area-a) <= sc.doc.ModelAbsoluteTolerance**2:
                                # Also check the center of the bounding box before eliminating this surface.
                                for c in boundingBoxCenters:
                                    if center.DistanceTo(c) <= sc.doc.ModelAbsoluteTolerance:
                                        return True

                        return False

                    if isFilletInList():
                        rgB.Dispose()
                        continue


                    rgBs_PerCFS.append(rgBs_Fillet)
                    fTols_Out.append(fTol)
                    areas.append(area)
                    boundingBoxCenters.append(center)
        if not rgBs_PerCFS:
            if bEcho: print "No fillets were created."


    return rgBs_PerCFS, fTols_Out


def createBrepObjects(objref_A, objref_B, fRadius, fTolerance=None, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bTryOtherTols = getOpt('bTryOtherTols')
    #bRebuildDeg5 = getOpt('bRebuildDeg5')
    bExtend = getOpt('bExtend')
    bTrimFillets = getOpt('bTrimFillets')
    bUseUnderlyingSrf = getOpt('bUseUnderlyingSrf')
    bTryOtherPickPts = getOpt('bTryOtherPickPts')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    if bUseUnderlyingSrf:
        # UnderlyingSurface() may not be needed.
        fA = objref_A.Face().UnderlyingSurface().ToBrep().Faces[0] 
        fB = objref_B.Face().UnderlyingSurface().ToBrep().Faces[0]
    else:
        fA = objref_A.Face()
        fB = objref_B.Face()

    ptA = objref_A.SelectionPoint()
    ptB = objref_B.SelectionPoint()

    gBreps_Fillets = []


    rc = createBreps(
        fA,
        ptA,
        fB,
        ptB,
        fRadius=fRadius,
        fTolerance=fTolerance,
        bTryOtherTols=bTryOtherTols,
        bExtend=bExtend,
        bTrimFillets=bTrimFillets,
        bUseUnderlyingSrf=bUseUnderlyingSrf,
        bTryOtherPickPts=bTryOtherPickPts,
        bDebug=bDebug,
        )
    if not rc:
        print "Fillets were not created using provided points."
        return

    rgBs_PerCFS, fTols_Ret = rc

    for rgBs_ThisCFS, fTol in zip(rgBs_PerCFS, fTols_Ret):
        if bDebug: print rgBs_ThisCFS, fTol

        for rgB in rgBs_ThisCFS:
            if bTrimFillets:
                gBrep_Fillet = sc.doc.Objects.AddBrep(rgB)
                if gBrep_Fillet != gBrep_Fillet.Empty:
                    gBreps_Fillets.append(gBrep_Fillet)
                else:
                    if not bExtend and not rgB.IsSurface:
                        print "Non-extend fillet created that is trimmed."
            else:
                for f in rgB.Faces:
                    gBrep_Fillet = sc.doc.Objects.AddBrep(f.ToBrep())
                    if gBrep_Fillet != gBrep_Fillet.Empty:
                        gBreps_Fillets.append(gBrep_Fillet)

    if bEcho:
        s = "Added {} brep(s) to document".format(len(gBreps_Fillets))

        if bTryOtherTols:
            s += " of tolerance(s) [{}].".format(','.join(str(f) for f in fTols_Ret))
        else:
            s += "."

        print s

    return gBreps_Fillets


def main():

    if Rhino.RhinoApp.ExeVersion == 5:
        print "This script uses Brep.CreateFilletSurface that is missing from Rhino V5's RC."
        return

    rc = getInput_FaceA()
    if rc is None: return
    (
        objref_A,
        bMultiPickB,
        bEcho,
        bDebug,
        ) = rc

    while True:
        sc.escape_test()

        rc = getInput_FaceB(objref_A)
        if rc is None: break

        (
            objref_B,
            fRadius,
            fTolerance,
            bTryOtherTols,
            bExtend,
            bTrimFillets,
            bUseUnderlyingSrf,
            bTryOtherPickPts,
            bEcho,
            bDebug,
            ) = rc

        gBreps_Fillets = createBrepObjects(
            objref_A,
            objref_B,
            fRadius=fRadius,
            fTolerance=fTolerance,
            bTryOtherTols=bTryOtherTols,
            bExtend=bExtend,
            bTrimFillets=bTrimFillets,
            bUseUnderlyingSrf=bUseUnderlyingSrf,
            bTryOtherPickPts=bTryOtherPickPts,
            bEcho=bEcho,
            bDebug=bDebug,
            )

        if not gBreps_Fillets and not bMultiPickB: return

        if not bMultiPickB: break

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()