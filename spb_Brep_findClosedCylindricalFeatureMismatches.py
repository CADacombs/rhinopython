"""
Use this script to check for possible mismatches in closed cylindrical face sets
between breps to find modeling errors between overlapping holes or holes with
their respective fasteners.

_ExplodeBlock all instances that should be checked. Remember to _UndoMultiple
the _ExplodeBlock later. It should respect the layer state, not process invisible objects.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
250203-06: Created.

TODO:
    Add support for block instances so that _ExplodeBlock does not need to be
    used before this script.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fTol_IsCyl_Dist'; keys.append(key)
    values[key] = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    names[key] = 'IsCylDistTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol_Cocyl_AdjFace_Angle_Deg'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'AdjFaceCocylParallelAngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fTol_Mismatch_Min'; keys.append(key)
    values[key] = 1e-6 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Millimeters,
        sc.doc.ModelUnitSystem)
    names[key] = 'MinMismatchDistTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol_Mismatch_Ignore'; keys.append(key)
    values[key] = 0.1 * Rhino.RhinoMath.UnitScale(
        Rhino.UnitSystem.Inches,
        sc.doc.ModelUnitSystem)
    names[key] = 'DistToIgnore'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'fTol_Mismatch_Deg'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAngleToleranceDegrees
    names[key] = 'MismatchParallelAngleTol'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})'.format(key, __file__)

    #key = 'bDot'; keys.append(key)
    #values[key] = True
    #riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    #stickyKeys[key] = '{}({})'.format(key, __file__)

    #key = 'iDotHeight'; keys.append(key)
    #values[key] = 11
    #riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    #stickyKeys[key] = '{}({})'.format(key, __file__)

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
        _debug = sc.sticky
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

        if key == 'fTol_Mismatch_Ignore':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            if cls.riOpts[key].CurrentValue < cls.riOpts['fTol_Mismatch_Min'].CurrentValue:
                cls.riOpts[key].CurrentValue = cls.values[key]
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def _get_all_normal_breps():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep
    return list(sc.doc.Objects.GetObjectList(oes))


def _getAllNormal_breps_and_instances():
    oes = rd.ObjectEnumeratorSettings()
    oes.LockedObjects = False # Default is True.
    oes.ObjectTypeFilter = rd.ObjectType.Brep | rd.ObjectType.InstanceReference
    return list(sc.doc.Objects.GetObjectList(oes))


def getInput():
    """
    Get Breps with optional input
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps")
    go.SetCommandPromptDefault("All normal when none are selected")

    go.GeometryFilter = rd.ObjectType.Brep
    go.SubObjectSelect = False

    go.AcceptNothing(True)

    #sc.doc.Views.Redraw()

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fTol_IsCyl_Dist')
        addOption('fTol_Cocyl_AdjFace_Angle_Deg')
        addOption('fTol_Mismatch_Min')
        addOption('fTol_Mismatch_Ignore')
        addOption('fTol_Mismatch_Deg')
        #addOption('bDot')
        #if Opts.values['bDot']:
        #    addOption('iDotHeight')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        # Not setting minimumNumber to 2 to avoid a result of GetResult.Nothing
        # for selection of 1 object.

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            return _get_all_normal_breps()

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return [objref.Object() for objref in objrefs]

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def _collect_cylinders_per_face_index(rgBrep, fTol_IsCyl_Dist):
    """
    Parameters:
    Returns:
    """
    cyls_perF = {}

    for rgF in rgBrep.Faces:
        #if rg.Surface.IsCylinder(
        #    rgF,
        #    tolerance=fTol_IsCyl_Dist)
        bSuccess, cyl = rg.Surface.TryGetCylinder(
            rgF,
            tolerance=fTol_IsCyl_Dist)
        if bSuccess:
            cyls_perF[rgF.FaceIndex] = cyl
            continue

        if isinstance(rgF.UnderlyingSurface(), rg.NurbsSurface):
            continue

        # Due to bug in TryGetCylinder, try one more time with the NURBS-equivalent surface.
        ns = rgF.ToNurbsSurface()
        bSuccess, cyl = rg.Surface.TryGetCylinder(
            ns,
            tolerance=fTol_IsCyl_Dist)
        if bSuccess:
            cyls_perF[rgF.FaceIndex] = cyl

    return cyls_perF


def _are_cylinders_colcylindrical(cylA, cylB, distanceTolerance, angleTolerance):
    if abs(cylA.Radius - cylB.Radius) > distanceTolerance:
        return False

    if not cylA.Axis.EpsilonEquals(cylB.Axis, epsilon=distanceTolerance):
        if not cylA.Axis.EpsilonEquals(-cylB.Axis, epsilon=distanceTolerance):
            return False

    isParallelTo = cylA.Axis.IsParallelTo(
        cylB.Axis,
        angleTolerance=angleTolerance)

    if isParallelTo == 0:
        sEval = "cylA.Axis"; print(sEval,'=',eval(sEval))
        sEval = "cylB.Axis"; print(sEval,'=',eval(sEval))
        angle = rg.Vector3d.VectorAngle(cylA.Axis, cylB.Axis)
        #sEval = "angle"; print(sEval,'=',eval(sEval))
        #sEval = "Rhino.RhinoMath.ToDegrees(angle)"; print(sEval,'=',eval(sEval))
        #sEval = "fTol_IsCyl_Dist"; print(sEval,'=',eval(sEval))
        sc.doc.Objects.AddCircle(cylA.CircleAt(linearParameter=0))
        sc.doc.Objects.AddCircle(cylB.CircleAt(linearParameter=0))
        raise Exception("isParallelTo of {} should not be so because EpsilonEquals of previous block passed.".format(
            isParallelTo))
        return False

    cirA = cylA.CircleAt(0.0)
    planeA = cirA.Plane
    cirB = cylB.CircleAt(0.0)
    xf = rg.Transform.PlanarProjection(planeA)
    ncB = cirB.ToNurbsCurve() # So that curve can deform as needed.
    ncB.Transform(xf)
    #sc.doc.Objects.AddCurve(ncB)
    bSuccess, cirB = ncB.TryGetCircle(distanceTolerance)
    if not bSuccess:
        return False

    if not cirA.Center.EpsilonEquals(cirB.Center, distanceTolerance):
        return False

    if not Rhino.RhinoMath.EpsilonEquals(cirA.Radius, cirB.Radius, distanceTolerance):
        return False

    return True


    ac_A = rg.ArcCurve(cirA)
    ac_B = rg.ArcCurve(cirB)
    
    if not cirA.EpsilonEquals(cirB, distanceTolerance):
        rvs = rg.Curve.GetDistancesBetweenCurves(ac_A, ac_B, 0.1*distanceTolerance)
        if not rvs[0]:
            return False

        if rvs[1] > distanceTolerance:
            return False
        
        #sc.doc.Objects.AddCircle(rgCirA)
        #sc.doc.Objects.AddCircle(rgCirB)

    return True


def _findContiguousCocylindricFacesOf1Face(rgBrep, cyls_perF, idx_F_Start, idxs_Fs_Cocyl_In, idxs_Fs_Fails_In, fTol_IsCyl_Dist, fTol_Cocyl_AdjFace_Angle_Deg):

    idxs_Fs_Cocyl_Out = idxs_Fs_Cocyl_In[:]
    idxs_Fs_Fails_Out = idxs_Fs_Fails_In[:]

    idxs_FAs = sorted(rgBrep.Faces[idx_F_Start].AdjacentFaces())
    #sEval = "idxs_FAs"; print(sEval,'=',eval(sEval))

    if len(idxs_FAs) == 0:
        return idxs_Fs_Cocyl_Out, idxs_Fs_Fails_Out

    cyl_This = cyls_perF[idx_F_Start]

    for iF_Adj in idxs_FAs:
        if iF_Adj not in cyls_perF:
            continue

        if iF_Adj in idxs_Fs_Cocyl_Out:
            continue

        if iF_Adj in idxs_Fs_Fails_In:
            continue

        #sEval = "iF, iF_Adj"; print(sEval,'=',eval(sEval))

        #sEval = "cyls_perF[iF].IsFinite"; print(sEval,'=',eval(sEval))

        cyl_Other = cyls_perF[iF_Adj]

        if _are_cylinders_colcylindrical(
            cylA=cyl_This,
            cylB=cyl_Other,
            distanceTolerance=fTol_IsCyl_Dist,
            angleTolerance=fTol_Cocyl_AdjFace_Angle_Deg
        ):
            idxs_Fs_Cocyl_Out.append(iF_Adj)
        else:
            idxs_Fs_Fails_Out.append(iF_Adj)

        #if not cyls_perF[iF].IsFinite and not cyls_perF[iF_Adj].IsFinite:
        #    if cyls_perF[iF].Axis.EpsilonEquals(
        #        cyls_perF[iF_Adj].Axis, epsilon=fTol_IsCyl_Dist
        #    ):
        #        if cyls_perF[iF].Center.EpsilonEquals(
        #            cyls_perF[iF_Adj].Center, epsilon=fTol_IsCyl_Dist
        #        ):
        #            cyls_Out.append(cyls_perF[iF])
        #            continue

    return idxs_Fs_Cocyl_Out, idxs_Fs_Fails_Out


def _get_lists_of_face_indices_for_cocylindric_faces(rgBrep, cyls_perF, fTol_IsCyl_Dist, fTol_Cocyl_AdjFace_Angle_Deg):

    idxs_Fs_perCyl_Out = [] # List of lists of int indices.

    idxs_Fs_Processed = [] # Flat list.

    for i, iF in enumerate(sorted(cyls_perF)):
        #print("{} of {}".format(i+1, len(cyls_perF)))
        if i == 4:
            pass

        if iF in idxs_Fs_Processed:
            continue

        idxs_Fs_Cocyls = [iF]
        idxs_Fs_Fails = idxs_Fs_Processed

        i = 0

        while True:
            sc.escape_test()

            iF_Start = idxs_Fs_Cocyls[i]

            idxs_Fs_Cocyls_Res, idxs_Fs_Fails_Res = _findContiguousCocylindricFacesOf1Face(
                rgBrep,
                cyls_perF,
                iF_Start,
                idxs_Fs_Cocyls,
                idxs_Fs_Fails,
                fTol_IsCyl_Dist,
                fTol_Cocyl_AdjFace_Angle_Deg=fTol_Cocyl_AdjFace_Angle_Deg)

            idxs_Fs_Cocyls = idxs_Fs_Cocyls_Res
            idxs_Fs_Fails = idxs_Fs_Fails_Res

            i += 1

            if i > (len(idxs_Fs_Cocyls) - 1):
                break

        idxs_Fs_perCyl_Out.append(idxs_Fs_Cocyls)
        idxs_Fs_Processed.extend(idxs_Fs_Cocyls)

    return idxs_Fs_perCyl_Out


def _find_closed_cylinders_per_face_sets(rgBrep, cyls_perF, idxFs_perCyl, fTol_IsCyl_Dist):
    rgB_Dup = rgBrep.DuplicateBrep()
    rgB_Dup.Faces.ShrinkFaces()

    cyls_Out = []

    for lst_idxF in idxFs_perCyl:
        if len(lst_idxF) == 1:
            iF = lst_idxF[0]
            if rgB_Dup.Faces[iF].IsClosed(0) or rgB_Dup.Faces[iF].IsClosed(1):
                cyls_Out.append(cyls_perF[iF])
            continue

        rgBs_DupFaces = [rgB_Dup.Faces[iF].DuplicateFace(duplicateMeshes=False) for iF in lst_idxF]

        rgBs_Joined = rg.Brep.JoinBreps(
            brepsToJoin=rgBs_DupFaces,
            tolerance=2.0*sc.doc.ModelAbsoluteTolerance)

        for _ in rgBs_DupFaces: _.Dispose()

        if len(rgBs_Joined) != 1:
            raise Exception("Number of breps: {}".format(len(rgBs_Joined)))

        rgB_Joined = rgBs_Joined[0]
        #sc.doc.Objects.AddBrep(rgB_Joined)

        cyl = cyls_perF[lst_idxF[0]]
        #sEval = "cyl.IsFinite"; print(sEval,'=',eval(sEval))
        cyl.Height1 = 1000.0
        #sEval = "cyl.IsFinite"; print(sEval,'=',eval(sEval))
        cyl.Height2 = -1000.0
        #sEval = "cyl.IsFinite"; print(sEval,'=',eval(sEval))
        srf_cyl = cyl.ToRevSurface()
        #sc.doc.Objects.AddSurface(srf_cyl)

        rgBs_CutUpSurface = rg.Brep.CutUpSurface(
            surface=srf_cyl,
            curves=[_ for _ in rgB_Joined.Edges if _.Valence == rg.EdgeAdjacency.Naked],
            useEdgeCurves=True,
            tolerance=sc.doc.ModelAbsoluteTolerance)

        rgB_Joined.Dispose()

        if len(rgBs_CutUpSurface) != 1:
            raise Exception("Number of breps: {}".format(len(rgBs_CutUpSurface)))

        rgB_CutUpSurface = rgBs_CutUpSurface[0]

        rgB_CutUpSurface.Faces.ShrinkFaces()

        #sc.doc.Objects.AddBrep(rgB_CutUpSurface)

        if rgB_CutUpSurface.Faces[0].IsClosed(0) or rgB_CutUpSurface.Faces[0].IsClosed(1):
            cyls_Out.append(cyl)

        rgB_CutUpSurface.Dispose()

    rgB_Dup.Dispose()

    return cyls_Out


def _areLinesParallel(lineA, lineB, distanceTolerance, angleTolerance):
    isParallelTo = lineA.Direction.IsParallelTo(
        lineB.Direction,
        angleTolerance=angleTolerance)
    return isParallelTo != 0


def _areLinesCollinear(lineA, lineB, distanceTolerance, angleTolerance):
    isParallelTo = lineA.Direction.IsParallelTo(
        lineB.Direction,
        angleTolerance=angleTolerance)
    if isParallelTo == 0:
        return False

    #bSuccess, tA, tB = rg.Intersect.Intersection.LineLine(lineA, lineB)
    ## bSuccess will be False if lines are parallel.
    #if bSuccess:
    #    return False

    ptB = lineB.From
    ptA = lineA.ClosestPoint(ptB, limitToFiniteSegment=False)
    dist = ptA.DistanceTo(ptB)
    return dist <= distanceTolerance


def _create_geom_to_mark_possible_mismatches(cyls_Closed_perB, fTol_Mismatch_Min, fTol_Mismatch_Ignore, fTol_Mismatch_Deg):
    fTol_Mismatch_Rad = Rhino.RhinoMath.ToRadians(fTol_Mismatch_Deg)

    rgPts_Out = []
    rgLines_Out = []

    for iThisBrep in range(len(cyls_Closed_perB)-1):
        for iOtherBrep in range(iThisBrep+1, len(cyls_Closed_perB)):
            sc.escape_test()
            #print(iThisBrep, iOtherBrep)
            cyls_This = cyls_Closed_perB[iThisBrep]
            cyls_Other = cyls_Closed_perB[iOtherBrep]
            for iCyl_This in range(len(cyls_This)):
                cyl_This = cyls_This[iCyl_This]
                center_This = cyl_This.Center
                line_This = rg.Line(center_This, cyl_This.Axis)
                for iCyl_Other in range(len(cyls_Other)):
                    cyl_Other = cyls_Other[iCyl_Other]

                    # WIP:

                    #isParallelTo = cyl_This.Axis.IsParallelTo(
                    #    cyl_Other.Axis,
                    #    angleTolerance=angleTolerance)
                    #sEval = "isParallelTo"; print(sEval,'=',eval(sEval))

                    #if isParallelTo == 0:
                    #    continue

                    #angle_Rad = rg.Vector3d.VectorAngle(cyl_This.Axis, cyl_Other.Axis) # Values range from 0 to 180 degrees.
                    #angle_Deg = Rhino.RhinoMath.ToDegrees(angle_Rad)
                    #sEval = "angle_Deg"; print(sEval,'=',eval(sEval))
                    #if angle_Rad > tol_Rad:
                    #    continue

                    center_Other = cyl_Other.Center
                    line_Other = rg.Line(cyl_Other.Center, cyl_Other.Axis)

                    if _areLinesParallel(
                        line_This,
                        line_Other,
                        distanceTolerance=fTol_Mismatch_Min,
                        angleTolerance=fTol_Mismatch_Rad
                    ):
                        if _areLinesCollinear(
                            line_This,
                            line_Other,
                            distanceTolerance=fTol_Mismatch_Min,
                            angleTolerance=fTol_Mismatch_Rad
                        ):
                            continue
                        ptA = line_This.From
                        ptB = line_Other.ClosestPoint(ptA, limitToFiniteSegment=False)
                    else:
                        bSuccess, tA, tB = rg.Intersect.Intersection.LineLine(line_This, line_Other)
                        if not bSuccess:
                            # Will be False if lines are parallel.
                            raise Exception("Intersection.LineLine returned False.")
                            continue
                        ptA = line_This.PointAt(tA)
                        ptB = line_Other.PointAt(tB)

                    dist = ptA.DistanceTo(ptB)

                    if fTol_Mismatch_Min < dist <= fTol_Mismatch_Ignore:
                        rgPts_Out.extend((ptA, ptB))
                        rgLines_Out.extend(())
                        if ptA.DistanceTo(center_This) > 100.0 * fTol_Mismatch_Min:
                            rgLines_Out.append(rg.Line(center_This, ptA))
                        else:
                            rgLines_Out.append(line_This)
                        if ptA.DistanceTo(center_Other) > 100.0 * fTol_Mismatch_Min:
                            rgLines_Out.append(rg.Line(center_Other, ptB))
                        else:
                            rgLines_Out.append(line_Other)

    return rgPts_Out, rgLines_Out


def processBreps(rgBreps, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTol_IsCyl_Dist = getOpt('fTol_IsCyl_Dist')
    fTol_Cocyl_AdjFace_Angle_Deg = getOpt('fTol_Cocyl_AdjFace_Angle_Deg')
    fTol_Mismatch_Min = getOpt('fTol_Mismatch_Min')
    fTol_Mismatch_Ignore = getOpt('fTol_Mismatch_Ignore')
    fTol_Mismatch_Deg = getOpt('fTol_Mismatch_Deg')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')

    if len(rgBreps) == 0:
        return

    if len(rgBreps) == 1:
        return

    rgBs_withCyls = []
    dicts_face_cyl_perB = []
    Rhino.RhinoApp.SetCommandPromptMessage("Collecting cylinders in each brep...")
    for rgB in rgBreps:
        dict_face_cyl = _collect_cylinders_per_face_index(
            rgB,
            fTol_IsCyl_Dist)
        if dict_face_cyl:
            rgBs_withCyls.append(rgB)
            dicts_face_cyl_perB.append(dict_face_cyl)

    idxFs_perCyl_perB = []
    Rhino.RhinoApp.SetCommandPromptMessage("Finding cocylindric faces in each brep...")
    for iB, rgB in enumerate(rgBs_withCyls):
        #Rhino.RhinoApp.SetCommandPromptMessage(
        #    "Processing cylinders from brep {} against brep {}, both out of {}".format(
        #        iThisBrep+1, iOtherBrep+1, len(cyls_Closed_perB)))

        dict_face_cyl = dicts_face_cyl_perB[iB]
        idxFs_perCyl = _get_lists_of_face_indices_for_cocylindric_faces(
            rgB,
            dict_face_cyl,
            fTol_IsCyl_Dist,
            fTol_Cocyl_AdjFace_Angle_Deg)
        #sEval = "idxFs_perCyl"; print(sEval,'=',eval(sEval))
        idxFs_perCyl_perB.append(idxFs_perCyl)

    cyls_Closed_perB = []
    #sCPM = 
    #Rhino.RhinoApp.SetCommandPromptMessage("Finding closed, cocylindric faces in each brep...")

    for iB, rgB in enumerate(rgBs_withCyls):
        sc.escape_test()
        Rhino.RhinoApp.SetCommandPromptMessage(
            "Finding closed, cocylindric faces in brep {} of {}".format(
                iB+1, len(rgBs_withCyls)))

        dict_face_cyl = dicts_face_cyl_perB[iB]
        idxFs_perCyl = idxFs_perCyl_perB[iB]
        cyls_Closed = _find_closed_cylinders_per_face_sets(
            rgB,
            dict_face_cyl,
            idxFs_perCyl,
            fTol_IsCyl_Dist)
        cyls_Closed_perB.append(cyls_Closed)

    #sEval = "len(rdBs_withCyls)"; print(sEval,'=',eval(sEval))
    #sEval = "len(cyls_Closed_perB)"; print(sEval,'=',eval(sEval))

    Rhino.RhinoApp.SetCommandPromptMessage("Finding possible mismatches...")

    return _create_geom_to_mark_possible_mismatches(
        cyls_Closed_perB,
        fTol_Mismatch_Min=fTol_Mismatch_Min,
        fTol_Mismatch_Ignore=fTol_Mismatch_Ignore,
        fTol_Mismatch_Deg=fTol_Mismatch_Deg)


def processBrepObjects(rdBreps, **kwargs):
    """
    """

    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    fTol_IsCyl_Dist = getOpt('fTol_IsCyl_Dist')
    fTol_Cocyl_AdjFace_Angle_Deg = getOpt('fTol_Cocyl_AdjFace_Angle_Deg')
    fTol_Mismatch_Min = getOpt('fTol_Mismatch_Min')
    fTol_Mismatch_Ignore = getOpt('fTol_Mismatch_Ignore')
    fTol_Mismatch_Deg = getOpt('fTol_Mismatch_Deg')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if len(rdBreps) == 0:
        print("No breps.")
        return

    if len(rdBreps) == 1:
        print("Only 1 brep in input, so no mismatches can be evaluated.")
        return

    rgBreps = []

    for rdB in rdBreps:
        try:
            rgBreps.append(rdB.BrepGeometry)
        except:
            ss = [
                "Block instances found",
                "_ExplodeBlock the instances and rerun this script.",
                "(Remember to _UndoMultiple the _ExplodeBlock later.)",
                ]
            raise Exception("\n".join(ss))

    return processBreps(
        rgBreps,
        fTol_IsCyl_Dist=fTol_IsCyl_Dist,
        fTol_Cocyl_AdjFace_Angle_Deg=fTol_Cocyl_AdjFace_Angle_Deg,
        fTol_Mismatch_Min=fTol_Mismatch_Min,
        fTol_Mismatch_Ignore=fTol_Mismatch_Ignore,
        fTol_Mismatch_Deg=fTol_Mismatch_Deg,
        bEcho=bEcho,
        bDebug=bDebug,
        )


def main():

    rdBreps = getInput()
    if not rdBreps: return

    fTol_IsCyl_Dist = Opts.values['fTol_IsCyl_Dist']
    fTol_Cocyl_AdjFace_Angle_Deg = Opts.values['fTol_Cocyl_AdjFace_Angle_Deg']
    fTol_Mismatch_Min = Opts.values['fTol_Mismatch_Min']
    fTol_Mismatch_Ignore = Opts.values['fTol_Mismatch_Ignore']
    fTol_Mismatch_Deg = Opts.values['fTol_Mismatch_Deg']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled=False

    #res, objrefs = ri.RhinoGet.GetMultipleObjects(
    #    "Select breps <All normal>",
    #    acceptNothing=True,
    #    filter=rd.ObjectType.Brep)
    #if res != Rhino.Commands.Result.Success: return


    rv = processBrepObjects(
        rdBreps,
        fTol_IsCyl_Dist=fTol_IsCyl_Dist,
        fTol_Cocyl_AdjFace_Angle_Deg=fTol_Cocyl_AdjFace_Angle_Deg,
        fTol_Mismatch_Min=fTol_Mismatch_Min,
        fTol_Mismatch_Ignore=fTol_Mismatch_Ignore,
        fTol_Mismatch_Deg=fTol_Mismatch_Deg,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    if not rv: return

    rgPts_Res, rgLines_Res = rv

    if not (rgPts_Res or rgLines_Res):
        print("No possible mismatches found.")
        return

    attrib_Red = rd.ObjectAttributes()
    attrib_Red.ColorSource = rd.ObjectColorSource.ColorFromObject
    attrib_Red.ObjectColor = Color.Red

    gPts_Out = []
    gCrvs_Out = []

    sc.doc.Objects.UnselectAll()

    for pt in rgPts_Res:
        gPt = sc.doc.Objects.AddPoint(pt, attrib_Red)
        if gPt == Guid.Empty: raise Exception("Point could not be added.")
        gPts_Out.append(gPt)
        sc.doc.Objects.Select(gPt)

    for line in rgLines_Res:
        gCrv = sc.doc.Objects.AddLine(line, attrib_Red)
        if gCrv == Guid.Empty: raise Exception("Curve could not be added.")
        gCrvs_Out.append(gCrv)
        sc.doc.Objects.Select(gCrv)

    if bEcho:
        print("Possible mismatches found. {} points and {} curves added.".format(len(gPts_Out), len(gCrvs_Out)))

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()