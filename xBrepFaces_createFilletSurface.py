"""
180922-24: Created.
190626, 1206, 200215: Modified an option default value.
210218: Refactored getInput.  Removed the single-span, degree-3 bezier option.
210224: Added an option.
210327: Bug fix.
210411: Added bIncludePlusMinusTols.  Moved some of createBreps to the new createBrepObjects.
210418: Added more printed feedback.  Modified an option default value.
210923: Disabled option for bRebuildDeg5 since its code hasn't yet been implemented.

_FilletSrf and Brep.CreateFilletSurface create surfaces that are degree 2 rational in one direction.


TODO: Add code for bRebuildDeg5.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

from System import Guid


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

    key = 'fTolerance'; keys.append(key)
    values[key] = 1.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = "Tol"
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bIncludePlusMinusTols'; keys.append(key)
    values[key] = False
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

    key = 'bTryMorePickPts'; keys.append(key)
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

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    """

    # Get face A with optional input.
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face A")
    
    go.GeometryFilter = rd.ObjectType.Surface
    #    goA.GeometryFilter = rd.ObjectType.Brep
    #    goA.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.SubSurface

    go.AcceptNumber(True, True)
    go.EnableHighlight(False)
    
    idxs_Opt = {}

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fRadius')
        addOption('fTolerance')
        addOption('bIncludePlusMinusTols')
        #addOption('bRebuildDeg5')
        addOption('bExtend')
        addOption('bTrimFillets')
        addOption('bUseUnderlyingSrf')
        addOption('bTryMorePickPts')
        addOption('bEcho')
        addOption('bDebug')


        res = go.Get()
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            if sc.doc.Objects.UnselectAll() > 0: sc.doc.Views.Redraw() # Allows for face highlighting.
    
            objrefA = go.Object(0)
            go.Dispose()
            gBrepA = objrefA.ObjectId
            rgBrep = objrefA.Brep()
            compIdxA = objrefA.GeometryComponentIndex.Index
            idx_rgFaceA = 0 if compIdxA == -1 else compIdxA
    
            rdBrepA = objrefA.Object()
            rdBrepA.HighlightSubObject(objrefA.GeometryComponentIndex, highlight=True)
            sc.doc.Views.Redraw()

            break


        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            key = 'fRadius'
            goNum =  go.Number()
            if goNum < 1e-6:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            else:
                Opts.riOpts[key].CurrentValue = goNum
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


    # Get face B with optional input.

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face B")

    go.GeometryFilter = rd.ObjectType.Surface

    go.AcceptNumber(True, True)
    go.EnableHighlight(False)
    
    while True:
        sc.escape_test()

        go.ClearCommandOptions()

        idxs_Opt.clear()

        def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

        addOption('fRadius')
        addOption('fTolerance')
        addOption('bIncludePlusMinusTols')
        addOption('bRebuildDeg5')
        addOption('bExtend')
        addOption('bTrimFillets')
        addOption('bUseUnderlyingSrf')
        addOption('bTryMorePickPts')
        addOption('bEcho')
        addOption('bDebug')


        if sc.escape_test(False):
            rdBrepA.HighlightSubObject(objrefA.GeometryComponentIndex, highlight=False)
            sc.doc.Views.Redraw()
            return
        
        res = go.Get()

        if res == ri.GetResult.Cancel:
            rdBrepA.HighlightSubObject(objrefA.GeometryComponentIndex, highlight=False)
            sc.doc.Views.Redraw()
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            rdBrepA.HighlightSubObject(objrefA.GeometryComponentIndex, highlight=False)
            sc.doc.Views.Redraw()

            objrefB = go.Object(0)
            go.Dispose()
            gBrepB = objrefB.ObjectId
            compIdxB = objrefB.GeometryComponentIndex.Index
            idx_rgFaceB = 0 if compIdxB == -1 else compIdxB
            rdBrepA.HighlightSubObject(objrefA.GeometryComponentIndex, highlight=False)
            sc.doc.Views.Redraw()
    
            pt_PickedA = objrefA.SelectionPoint()
            pt_PickedB = objrefB.SelectionPoint()
    
            return (
                gBrepA, idx_rgFaceA, pt_PickedA,
                gBrepB, idx_rgFaceB, pt_PickedB,
                [Opts.values[key] for key in Opts.keys]
            )


        # An option was selected or a number was entered.
        
        if res == ri.GetResult.Number:
            key = 'fRadius'
            goNum =  go.Number()
            if goNum < 1e-6:
                Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            else:
                Opts.riOpts[key].CurrentValue = goNum
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

    bIncludePlusMinusTols = getOpt('bIncludePlusMinusTols')
    bRebuildDeg5 = getOpt('bRebuildDeg5')
    bExtend = getOpt('bExtend')
    bTrimFillets = getOpt('bTrimFillets')
    bUseUnderlyingSrf = getOpt('bUseUnderlyingSrf')
    bTryMorePickPts = getOpt('bTryMorePickPts')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if bUseUnderlyingSrf:
        rgFaceA = rgFaceA.UnderlyingSurface().ToBrep().Faces[0] # UnderlyingSurface() is not needed.
        rgFaceB = rgFaceB.UnderlyingSurface().ToBrep().Faces[0]
    else:
        # No change in rgFaceA or rgFaceB.
        pass
    
    rgBs_Fillet_PerCFS_PerTol = [] # 2-level nested lists of Breps.
    fTols_Out = [] # In-syn with rgBs_Fillet_PerCFS_PerTol.
    
    areas = []
    boundingBoxCenters = []

    b, u, v = rgFaceA.ClosestPoint(pt3d_A)
    pt2d_uvA = rg.Point2d(u, v)
    b, u, v = rgFaceB.ClosestPoint(pt3d_B)
    pt2d_uvB = rg.Point2d(u, v)


    fTol_Nom = sc.doc.ModelAbsoluteTolerance if fTolerance is None else fTolerance

    fTols = (fTol_Nom, 10.0*fTol_Nom, 0.1*fTol_Nom) if bIncludePlusMinusTols else (fTol_Nom,)

    for fTol in fTols:
        rgBs_Fillet_PerCFS = [] # 1-level nested.

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

        areas_ThisTol = []

        for rgBrep in rgBs_Fillet:
            # Record area for comparison with other fillets.
            area = rgBrep.GetArea()
            if not area:
                print "Area could not be obtained."
            # Record bounding box center point for comparison with other fillets.
            bbox = rgBrep.GetBoundingBox(accurate=False)
            if not bbox:
                print "Bounding box could not be obtained."
            center = bbox.Center
            
            rgBs_Fillet_PerCFS.append(rgBs_Fillet)
            fTols_Out.append(fTol)
            areas.append(area)
            areas_ThisTol.append(area)
            boundingBoxCenters.append(center)

        rgBs_Fillet_PerCFS_PerTol.append(rgBs_Fillet_PerCFS)

        if bEcho:
            print "Fillet with accumulative area of" \
                " {} square units as created at {} tolerance.".format(
                    sum(areas_ThisTol), fTol)


        if not bTryMorePickPts:
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

                for rgBrep in rgBs_Fillet:
                    # Don't include any area matches.


                    area = rgBrep.GetArea()
                    if not area:
                        print "Area could not be obtained."
                        continue

                    bbox = rgBrep.GetBoundingBox(accurate=False)
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
                        rgBrep.Dispose()
                        continue


                    rgBs_Fillet_PerCFS.append(rgBs_Fillet)
                    areas.append(area)
                    boundingBoxCenters.append(center)
        if not rgBs_Fillet_PerCFS:
            if bEcho: print "No fillets were created."

        rgBs_Fillet_PerCFS_PerTol.append(rgBs_Fillet_PerCFS)


    return rgBs_Fillet_PerCFS_PerTol, fTols_Out


def createBrepObjects(rgFaceA, pt3d_A, rgFaceB, pt3d_B, fRadius, fTolerance=None, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bIncludePlusMinusTols = getOpt('bIncludePlusMinusTols')
    bRebuildDeg5 = getOpt('bRebuildDeg5')
    bExtend = getOpt('bExtend')
    bTrimFillets = getOpt('bTrimFillets')
    bUseUnderlyingSrf = getOpt('bUseUnderlyingSrf')
    bTryMorePickPts = getOpt('bTryMorePickPts')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    Rhino.RhinoApp.SetCommandPrompt("Working ...")
    
    if bUseUnderlyingSrf:
        rgFaceA = rgFaceA.UnderlyingSurface().ToBrep().Faces[0] # UnderlyingSurface() is not needed.
        rgFaceB = rgFaceB.UnderlyingSurface().ToBrep().Faces[0]
    else:
        # No change in rgFaceA or rgFaceB.
        pass

    rc = createBreps(
        rgFaceA, pt3d_A,
        rgFaceB, pt3d_B,
        fRadius=fRadius,
        fTolerance=fTolerance,
        bIncludePlusMinusTols=bIncludePlusMinusTols,
        bRebuildDeg5=bRebuildDeg5,
        bExtend=bExtend,
        bTrimFillets=bTrimFillets,
        bUseUnderlyingSrf=bUseUnderlyingSrf,
        bTryMorePickPts=bTryMorePickPts,
        bDebug=bDebug,
        )
    if not rc:
        print "Fillets were not created using provided points."
        return

    rgBs_Fillet_PerCFS_PerTol, fTols_Ret = rc

    gBreps_Fillets = []
    
    for rgBs_Fillet_PerCFS, fTol in zip(rgBs_Fillet_PerCFS_PerTol, fTols_Ret):
        if bDebug: print rgBs_Fillet_PerCFS, fTol

        for rgBs_Fillet in rgBs_Fillet_PerCFS:
            for rgB in rgBs_Fillet:
                if bTrimFillets:
                    gBrep_Fillet = sc.doc.Objects.AddBrep(rgB)
                    if gBrep_Fillet != Guid.Empty:
                        gBreps_Fillets.append(gBrep_Fillet)
                    else:
                        if not bExtend and not rgB.IsSurface:
                            print "Non-extend fillet created that is trimmed."
                else:
                    for f in rgB.Faces:
                        gBrep_Fillet = sc.doc.Objects.AddBrep(f.ToBrep())
                        if gBrep_Fillet != Guid.Empty:
                            gBreps_Fillets.append(gBrep_Fillet)

    if bEcho:
        print "Added {} brep(s) to document of tolerance(s) [{}].".format(
            len(gBreps_Fillets),
            ','.join(str(f) for f in fTols_Ret)
            )

    return gBreps_Fillets


def main():

    if Rhino.RhinoApp.ExeVersion == 5:
        print "This script uses Brep.CreateFilletSurface that is missing from Rhino V5's RC."
        return

    rc = getInput()
    if rc is None: return

    (
        gBrepA, idx_rgFaceA, pt3d_A,
        gBrepB, idx_rgFaceB, pt3d_B,
        (
            fRadius,
            fTolerance,
            bIncludePlusMinusTols,
            bRebuildDeg5,
            bExtend,
            bTrimFillets,
            bUseUnderlyingSrf,
            bTryMorePickPts,
            bEcho,
            bDebug,
        )
    ) = rc

    rgBrepA = rs.coercebrep(gBrepA)
    rgBrepB = rs.coercebrep(gBrepB)

    rgFaceA = rgBrepA.Faces[idx_rgFaceA]
    rgFaceB = rgBrepB.Faces[idx_rgFaceB]

    gBreps_Fillets = []

    rc = createBrepObjects(
        rgFaceA, pt3d_A,
        rgFaceB, pt3d_B,
        fRadius=fRadius,
        fTolerance=fTolerance,
        bIncludePlusMinusTols=bIncludePlusMinusTols,
        bRebuildDeg5=bRebuildDeg5,
        bExtend=bExtend,
        bTrimFillets=bTrimFillets,
        bUseUnderlyingSrf=bUseUnderlyingSrf,
        bTryMorePickPts=bTryMorePickPts,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rc:
        gBreps_Fillets.extend(rc)

    if gBreps_Fillets:
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()