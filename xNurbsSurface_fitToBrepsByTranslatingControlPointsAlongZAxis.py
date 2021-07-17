"""
191018-19: Created.
210717: Points on breps to match are now closest to starting surface.
        Disabled options not yet implemented in main routine.

Starting surface's control points' X and Y are maintained.
Starting surface's Greville points are used for measurement.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Drawing import Color

#import itertools
import random


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    addOptions = {}
    stickyKeys = {}
    
    #key = 'bPickStartingSrf'; keys.append(key)
    #values[key] = True
    #names[key] = key[1:]
    #riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    #addOptions[key] = lambda getObj, key=key, names=names, riOpts=riOpts: (
    #        ri.Custom.GetBaseClass.AddOptionToggle(
    #            getObj, englishName=names[key], toggleValue=riOpts[key]))
    #stickyKeys[key] = '{}({})'.format(key, __file__)
    
    #key = 'fSamplingLength'; keys.append(key)
    #values[key] = 200.0 * sc.doc.ModelAbsoluteTolerance
    #names[key] = key[1:]
    #riOpts[key] = ri.Custom.OptionDouble(values[key])
    #addOptions[key] = lambda getObj, key=key, names=names, riOpts=riOpts: (
    #        ri.Custom.GetBaseClass.AddOptionDouble(
    #            getObj, englishName=names[key], numberValue=riOpts[key]))
    #stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bExtrapolateMissingPts'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    addOptions[key] = lambda getObj, key=key, names=names, riOpts=riOpts: (
            ri.Custom.GetBaseClass.AddOptionToggle(
                getObj, englishName=names[key], toggleValue=riOpts[key]))
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iExtrapolationCt'; keys.append(key)
    values[key] = 1
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=0)
    addOptions[key] = lambda getObj, key=key, names=names, riOpts=riOpts: (
            ri.Custom.GetBaseClass.AddOptionInteger(
                getObj, englishName=names[key], intValue=riOpts[key]))
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    addOptions[key] = lambda getObj, key=key, names=names, riOpts=riOpts: (
            ri.Custom.GetBaseClass.AddOptionToggle(
                getObj, englishName=names[key], toggleValue=riOpts[key]))
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    addOptions[key] = lambda getObj, key=key, names=names, riOpts=riOpts: (
            ri.Custom.GetBaseClass.AddOptionToggle(
                getObj, englishName=names[key], toggleValue=riOpts[key]))
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in cls.keys:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue


def getInput_ObjsToFit():
    """
    Get objects with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select breps to skin")

    # The following are objects to ignore instead since
    # BrepVertex cannot be added to a filter of objects to accept.
    #go.GeometryFilter &= ~ (
    #        rd.ObjectType.Annotation |
    #        rd.ObjectType.BrepLoop |
    #        rd.ObjectType.Cage |
    #        rd.ObjectType.ClipPlane |
    #        rd.ObjectType.Detail |
    #        rd.ObjectType.InstanceReference |
    #        rd.ObjectType.Light
    #)
    go.GeometryFilter = rd.ObjectType.Brep

    #go.AcceptNumber(True, acceptZero=True)

    #s  = "{}:".format(Opts.names['fSamplingLength'])
    #s += " Sampling division curve length for curves and"
    #s += " target edge length of meshes extracted from breps."
    #print s

    while True:
        #Opts.addOptions['bPickStartingSrf'](go)
        #Opts.addOptions['fSamplingLength'](go)
        Opts.addOptions['bExtrapolateMissingPts'](go)
        if Opts.values['bExtrapolateMissingPts']:
            Opts.addOptions['iExtrapolationCt'](go)
        Opts.addOptions['bEcho'](go)
        Opts.addOptions['bDebug'](go)

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in Opts.keys])

        #if res == ri.GetResult.Number:
        #    Opts.riOpts['fSamplingLength'].CurrentValue = abs(go.Number())
            
        #key = 'fSamplingLength'
        #if Opts.riOpts[key].CurrentValue <= 0.0:
        #    Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
            
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInput_srf_Starting():
    """
    Get Surface with optional input.
    """

    sc.doc.Objects.UnselectAll()
    sc.doc.Views.Redraw()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select surface to use its normals")

    # The following are objects to ignore instead since
    # BrepVertex cannot be added to a filter of objects to accept.
    go.GeometryFilter = rd.ObjectType.Surface

    #go.AcceptNumber(True, acceptZero=True)

    #s  = "{}:".format(Opts.names['fSamplingLength'])
    #s += " Sampling division curve length for curves and"
    #s += " target edge length of meshes extracted from breps."
    #print s

    while True:
        #Opts.addOptions['bPickStartingSrf'](go)
        #Opts.addOptions['fSamplingLength'](go)
        Opts.addOptions['bExtrapolateMissingPts'](go)
        if Opts.values['bExtrapolateMissingPts']:
            Opts.addOptions['iExtrapolationCt'](go)
        Opts.addOptions['bEcho'](go)
        Opts.addOptions['bDebug'](go)

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in Opts.keys])

        #if res == ri.GetResult.Number:
        #    Opts.riOpts['fSamplingLength'].CurrentValue = abs(go.Number())

        #key = 'fSamplingLength'
        #if Opts.riOpts[key].CurrentValue <= 0.0:
        #    Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def coerceBrep(rhObj):
    if isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        #print rhObj.GeometryComponentIndex.ComponentIndexType
        geom = rhObj.Geometry()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    if isinstance(geom, rg.Brep):
        return geom


def coerceSurface(rhObj):
    if isinstance(rhObj, rg.GeometryBase):
        geom = rhObj
    elif isinstance(rhObj, rd.ObjRef):
        #print rhObj.GeometryComponentIndex.ComponentIndexType
        geom = rhObj.Geometry()
    elif isinstance(rhObj, Guid):
        rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        geom = rdObj.Geometry
    else:
        return

    srf = None
    if isinstance(geom, rg.BrepFace):
        srf = geom.UnderlyingSurface()
    elif isinstance(geom, rg.Surface):
        srf = geom
    elif isinstance(geom, rg.Brep):
        if geom.Faces.Count == 1:
            srf = geom.Faces[0].UnderlyingSurface()

    return srf


def getPointsAtSurfaceParameters(rhObj_Srf):
    srf = coerceSurface(rhObj_Srf)
    if not srf: return
    
    if not srf.GetType() == rg.NurbsSurface:
        print "Not a NurbsSurface."
        return

    ns = srf
    pts_out = []

    for iU in range(ns.Points.CountU):
        pts_out.append([])
        for iV in range(ns.Points.CountV):
            uv = ns.Points.GetGrevillePoint(iU, iV)
            pt = ns.PointAt(uv[0], uv[1])
            pts_out[-1].append(pt)

    return pts_out


def getNormalsAtSurfaceParameters(rhObj_Srf, params):
    srf = coerceSurface(rhObj_Srf)
    if not srf: return
    
    if not srf.GetType() == rg.NurbsSurface:
        print "Not a NurbsSurface."
        return

    ns = srf
    vs_out = []

    for iU in range(ns.Points.CountU):
        vs_out.append([])
        for iV in range(ns.Points.CountV):
            u, v = params[iU][iV]
            normal = ns.NormalAt(u, v)
            vs_out[-1].append(normal)

    return vs_out


def createPointsProjectedToBreps(rhObjs_ProjectTo, pts_toProject):
    """
    rhObjects_Ref can include ObjRefs, DocObjects.RhinoObjects, GUIDS, or Geometry, but must be all the same type.
    """

    breps = []
    for obj in rhObjs_ProjectTo:
        brep = coerceBrep(obj)
        if brep: breps.append(brep)
    if not breps: return

    pts_Out = []
    
    for iU in range(len(pts_toProject)):
        pts_Out.append([])
        for iV in range(len(pts_toProject[0])):
            rc = rg.Intersect.Intersection.ProjectPointsToBreps(
                    breps=breps,
                    points=[pts_toProject[iU][iV]],
                    direction=rg.Vector3d.ZAxis,
                    tolerance=0.1*sc.doc.ModelAbsoluteTolerance)
            if len(rc) == 0:
                pt_to = None
            elif len(rc) == 1:
                pt_to = rc[0]
            else:
                pts = rc
                dists = []
                for pt in pts:
                    dist = pt.DistanceTo(pts_toProject[iU][iV])
                    dists.append(dist)
                pt_to = pts[dists.index(min(dists))]

            pts_Out[-1].append(pt_to)

    return pts_Out


def getFormattedDistance(fDistance):
    if fDistance is None: return "(No deviation provided)"
    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def flattenNestedPointList(pts_in):
    pts_out = []
    for iU in range(len(pts_in)):
        for iV in range(len(pts_in[iU])):
            pts_out.append(pts_in[iU][iV])
    return pts_out


def iterateFit(pts, ns_Starting):
    """
    """


    def getSeparatedPointIndices(pts, ns_Starting):
        idxs_Pts_SameAsStartingSrf = []
        idxs_pts_toFit = []
        for iU in range(len(pts)):
            for iV in range(len(pts[0])):
                if pts[iU][iV].EpsilonEquals(
                        ns_Starting.Points.GetControlPoint(iU,iV).Location,
                        epsilon=sc.doc.ModelAbsoluteTolerance
                ):
                    idxs_Pts_SameAsStartingSrf.append((iU,iV))
                else:
                    idxs_pts_toFit.append((iU,iV))
        return idxs_Pts_SameAsStartingSrf, idxs_pts_toFit

    idxs_Pts_SameAsStartingSrf, idxs_pts_toFit = getSeparatedPointIndices(pts, ns_Starting)

    #print idxs_Pts_SameAsStartingSrf
    #print idxs_pts_toFit

    pts_Projected_Flat = flattenNestedPointList(pts)

    ns1 = ns_Starting.Duplicate()

    for iU, iV in idxs_pts_toFit:
        ns1.Points.SetControlPoint(iU,iV,pts[iU][iV])
    
    for i in xrange(200):
        bTransPts = False
        for iU, iV in idxs_pts_toFit:
            uv_Gr = ns1.Points.GetGrevillePoint(iU, iV)
            pt_Gr = ns1.PointAt(uv_Gr[0], uv_Gr[1])
            dist = pt_Gr.DistanceTo(pts[iU][iV])
            vect = pts[iU][iV] - pt_Gr
            if vect.IsTiny(): continue
            #print dist, vect
            cp = ns1.Points.GetControlPoint(iU,iV)
            ns1.Points.SetControlPoint(iU,iV,cp.Location+vect)
            bTransPts = True
        if not bTransPts: break
    
    s  = "{} iterations".format(i+1)
    s += " for Grevilles to lie on target(s) within Vector3d.IsTiny()."
    print s

    return ns1


def fit_Geometry(breps_ProjectTo, srf_Starting, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bExtrapolateMissingPts = getOpt('bExtrapolateMissingPts')
    iExtrapolationCt = getOpt('iExtrapolationCt')
    #fSamplingLength = getOpt('fSamplingLength')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    ns_Starting = srf_Starting.ToNurbsSurface()

    pts_Greville = getPointsAtSurfaceParameters(ns_Starting)

    pts_Projected = createPointsProjectedToBreps(
            rhObjs_ProjectTo=breps_ProjectTo,
            pts_toProject=pts_Greville,
    )
    if not pts_Projected:
        print "Projected points were not obtained."
        return

    #for iU in range(len(pts_Nested)):
    #    for iV in range(len(pts_Nested[iU])):
    #        print iU, iV, pts_Nested[iU][iV]


    def getNeighborCount(pts, iU, iV):

        idx_MaxU = len(pts)-1
        idx_MaxV = len(pts[0])-1

        ct = 0

        if iU-1 >= 0 and pts[iU-1][iV] is not None:
            ct += 1
        if iU+1 <= idx_MaxU and pts[iU+1][iV] is not None:
            ct += 1
        if iV-1 >= 0 and pts[iU][iV-1] is not None:
            ct += 1
        if iV+1 <= idx_MaxV and pts[iU][iV+1] is not None:
            ct += 1

        return ct


    def closestPointsOfNeighborsOnNormalLines(pts_In, iU, iV, iMinNeighborCt=1, bDiag=True, bLineExts=True):
        """
        bLineExts:
            When True: If neighbor's point and its neighbor's point in the same
            direction are both available, create a line through those points
            and get the ClosestPoint of that line on the normal line.
        """
        pts_Out = []

        idx_MaxU = len(pts_In)-1
        idx_MaxV = len(pts_In[0])-1

        attr.ObjectColor = Color.FromArgb(
                red=random.randint(0, 255),
                green=random.randint(0, 255),
                blue=random.randint(0, 255))

        # West (Previous U).
        if iU-1 >= 0 and pts_In[iU-1][iV] is not None:
            pt = None
            if bLineExts and iU-2 >= 0 and pts_In[iU-2][iV] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU-2][iV], pts_In[iU-1][iV])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU-1][iV],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # East (Next U).
        if iU+1 <= idx_MaxU and pts_In[iU+1][iV] is not None:
            pt = None
            if bLineExts and iU+2 <= idx_MaxU and pts_In[iU+2][iV] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU+2][iV], pts_In[iU+1][iV])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU+1][iV],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # South (Previous V).
        if iV-1 >= 0 and pts_In[iU][iV-1] is not None:
            pt = None
            if bLineExts and iV-2 >= 0 and pts_In[iU][iV-2] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU][iV-2], pts_In[iU][iV-1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU][iV-1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # North (Next V).
        if iV+1 <= idx_MaxV and pts_In[iU][iV+1] is not None:
            pt = None
            if bLineExts and iV+2 <= idx_MaxV and pts_In[iU][iV+2] is not None:
                line_ThruNeighbors = rg.Line(pts_In[iU][iV+2], pts_In[iU][iV+1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU][iV+1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Southwest (Previous U, Previous V).
        if (
                bDiag and
                iU-1 >= 0 and
                iV-1 >= 0 and
                pts_In[iU-1][iV-1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU-2 >= 0 and
                    iV-2 >= 0 and
                    pts_In[iU-2][iV-2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU-2][iV-2], pts_In[iU-1][iV-1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU-1][iV-1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Southeast (Next U, Previous V).
        if (
                bDiag and
                iU+1 <= idx_MaxU and
                iV-1 >= 0 and
                pts_In[iU+1][iV-1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU+2 <= idx_MaxU and
                    iV-2 >= 0 and
                    pts_In[iU+2][iV-2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU+2][iV-2], pts_In[iU+1][iV-1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU+1][iV-1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Northwest (Previous U, Next V).
        if (
                bDiag and
                iU-1 >= 0 and
                iV+1 <= idx_MaxV and
                pts_In[iU-1][iV+1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU-2 >= 0 and
                    iV+2 <= idx_MaxV and
                    pts_In[iU-2][iV+2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU-2][iV+2], pts_In[iU-1][iV+1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU-1][iV+1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        # Northeast (Next U, Next V).
        if (
                bDiag and
                iU+1 <= idx_MaxU and
                iV+1 <= idx_MaxV and
                pts_In[iU+1][iV+1] is not None
        ):
            pt = None
            if (
                    bLineExts and
                    iU+2 <= idx_MaxU and
                    iV+2 <= idx_MaxV and
                    pts_In[iU+2][iV+2] is not None
            ):
                line_ThruNeighbors = rg.Line(pts_In[iU+2][iV+2], pts_In[iU+1][iV+1])
                line_ThruNeighbors.Length *= 2.0
                #sc.doc.Objects.AddLine(line_ThruNeighbors, attr)
                rc = rg.Intersect.Intersection.LineLine(
                        lineA=lines_thruStartingSrfNormals[iU][iV],
                        lineB=line_ThruNeighbors)
                if rc[0]:
                    pt = lines_thruStartingSrfNormals[iU][iV].PointAt(rc[1])
            else:
                pt = lines_thruStartingSrfNormals[iU][iV].ClosestPoint(
                        pts_In[iU+1][iV+1],
                        limitToFiniteSegment=False)
            if pt: pts_Out.append(pt)

        if len(pts_Out) < iMinNeighborCt:
            return []

        return pts_Out


    def addMissingPointsAlongBorder(pts, idxs_pt_filter=None, iMinNeighborCt=1, bDiag=False, bLineExts=True):
        """
        Modify a copy of the list so that new points do not affect
        subsequent ones in this function call.
        """

        pts0 = pts

        pts_Out = [p[:] for p in pts0]

        bModificationOccured = False

        for iU in range(len(pts0)):
            for iV in range(len(pts0[0])):
                if idxs_pt_filter and not (iU, iV) in idxs_pt_filter: continue

                if pts0[iU][iV] is not None: continue

                pts = closestPointsOfNeighborsOnNormalLines(
                        pts0,
                        iU,
                        iV,
                        iMinNeighborCt=iMinNeighborCt,
                        bDiag=bDiag,
                        bLineExts=bLineExts)
                if not pts: continue

                pt_Sum = None
                for pt in pts:
                    pt_Sum = pt if pt_Sum is None else pt_Sum + pt

                pt = pt_Sum / float(len(pts))
                #sc.doc.Objects.AddPoint(pt, attr)

                pts_Out[iU][iV] = pt

                bModificationOccured = True

        # Modify original list.
        for iU in range(len(pts0)):
            for iV in range(len(pts0[0])):
                pts0[iU][iV] = pts_Out[iU][iV]

        return bModificationOccured


    # Fill any missing points in grid by interpolation or extrapolation.

    def hasMissingPoints(pts):
        for iU in range(len(pts)):
            for iV in range(len(pts[0])):
                if pts[iU][iV] is None:
                    return True
        return False


    def missingPointIndices(pts):
        idxs_Pts_SameAsStartingSrf = []
        for iU in range(len(pts_forCTP)):
            for iV in range(len(pts_forCTP[0])):
                if pts_forCTP[iU][iV] is None:
                    idxs_Pts_SameAsStartingSrf.append((iU,iV))


    def addMissingPerStartingSrf(pts, ns_Starting):
        for iU in range(len(pts_forCTP)):
            for iV in range(len(pts_forCTP[0])):
                if pts_forCTP[iU][iV] is None:
                    pts_forCTP[iU][iV] = ns_Starting.Points.GetControlPoint(iU,iV).Location


    if not hasMissingPoints(pts_Projected):

        pts_forCTP = [p[:] for p in pts_Projected]

        pts_Projected_Flat = flattenNestedPointList(pts_forCTP)

        ns1 = iterateFit(pts_forCTP, ns_Starting)
        return ns1

    else:

        idxs_borderPts = []
        for iU in range(len(pts_Projected)):
            for iV in range(len(pts_Projected[0])):
                if (
                        pts_Projected[iU][iV] is None and
                        getNeighborCount(pts_Projected, iU, iV) > 0
                ):
                    idxs_borderPts.append((iU,iV))


        def createZAxisLinesAtGrevilles(pts_Greville_toProject):
            lines = []
            for iU in range(len(pts_Greville_toProject)):
                lines.append([])
                for iV in range(len(pts_Greville_toProject[0])):
                    line = rg.Line(
                            start=pts_Greville_toProject[iU][iV],
                            span=rg.Vector3d.ZAxis)
                    #sc.doc.Objects.AddLine(line)
                    lines[-1].append(line)
            return lines

        lines_thruStartingSrfNormals = createZAxisLinesAtGrevilles(pts_Greville)

        attr = rd.ObjectAttributes()
        attr.ColorSource = rd.ObjectColorSource.ColorFromObject


        if not bExtrapolateMissingPts:

            pts_forCTP = [p[:] for p in pts_Projected]

            # Fill any remaining missing points with the
            # starting surface control point locations.

            addMissingPerStartingSrf(pts_forCTP, ns_Starting)

            ns1 = iterateFit(pts_forCTP, ns_Starting)
            return ns1





        # Extrapolate missing points.

        #for bDiag1, bDiag2, bLineExts1, bLineExts2 in itertools.product(
        #            (False, True), (False, True), (False, True), (False, True)):
        for bLineExts1, bLineExts2 in ((True, False),):
        #for bLineExts1, bLineExts2 in itertools.product(
        #            (False, True), (False, True)):

            if not bLineExts1 and bLineExts2: continue

            pts_forCTP = [p[:] for p in pts_Projected]

            if bExtrapolateMissingPts:
                # Perform this block of code no matter the value of iExtrapolationCt.

                for iMinNeighborCt in 4,3,2,1: #(1,): #

                    idxs_pt_filter = None

                    addMissingPointsAlongBorder(
                            pts_forCTP,
                            idxs_pt_filter=idxs_borderPts,
                            iMinNeighborCt=iMinNeighborCt,
                            bDiag=False,
                            bLineExts=bLineExts1)
                    if not hasMissingPoints(pts_forCTP):
                        break


            i = 2 # For 2nd iteration (not index) to be relevant with iExtrapolationCt.

            while (
                    hasMissingPoints(pts_forCTP) and
                    ((iExtrapolationCt == 0) or (i < iExtrapolationCt))
            ):
                sc.escape_test()

                addMissingPointsAlongBorder(
                        pts_forCTP,
                        bDiag=False,
                        bLineExts=bLineExts2)
                i += 1

            # Fill any remaining missing points with the
            # starting surface control point locations.

            idxs_Pts_SameAsStartingSrf = []
            for iU in range(len(pts_forCTP)):
                for iV in range(len(pts_forCTP[0])):
                    if pts_forCTP[iU][iV] is None:
                        pts_forCTP[iU][iV] = ns_Starting.Points.GetControlPoint(iU,iV).Location
                        idxs_Pts_SameAsStartingSrf.append((iU,iV))

                #for pts_V in pts_forCTP:
                #    for pt in pts_V:
                #        if pt:
                #            sc.doc.Objects.AddPoint(pt)
                #sc.doc.Views.Redraw()
                #return


            pts_Projected_Flat = flattenNestedPointList(pts_forCTP)

            ns1 = iterateFit(pts_forCTP, ns_Starting)
            if not ns1: continue


            attr.ObjectColor = Color.FromArgb(
                    red=random.randint(0, 255),
                    green=random.randint(0, 255),
                    blue=random.randint(0, 255))

            attr.Name = (
                        'bLineExts1' + str(bLineExts1) + ',' +
                        'bLineExts2' + str(bLineExts2) + ',')
            sc.doc.Objects.AddSurface(ns1, attr)
            sc.doc.Views.Redraw()



def fit_DocObjects(objrefs_toFit, objref_srf_Starting, **kwargs):
    """
    """


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bExtrapolateMissingPts = getOpt('bExtrapolateMissingPts')
    iExtrapolationCt = getOpt('iExtrapolationCt')
    #fSamplingLength = getOpt('fSamplingLength')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    breps_ProjectTo = [coerceBrep(rhObj) for rhObj in objrefs_toFit]

    srf_Starting = coerceSurface(objref_srf_Starting)

    ns1 = fit_Geometry(
            breps_ProjectTo=breps_ProjectTo,
            srf_Starting=srf_Starting,
            )


    for brep in breps_ProjectTo: brep.Dispose()
    srf_Starting.Dispose()


    if not ns1: return

    g_ns1 = sc.doc.Objects.AddSurface(ns1)

    ns1.Dispose()

    if g_ns1 == Guid.Empty: return

    sc.doc.Views.Redraw()
    return g_ns1


def main():

    rc = getInput_ObjsToFit()
    if rc is None: return
    objrefs_toFit = rc[0]

    rc = getInput_srf_Starting()
    if rc is None: return
    objref_srf_Starting = rc[0]


    if Opts.values['bDebug']:
        pass
    else:
        pass


    fit_DocObjects(objrefs_toFit, objref_srf_Starting)


if __name__ == '__main__': main()