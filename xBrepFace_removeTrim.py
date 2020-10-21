"""
160518-20: Created.
160521: Now supports singular trims.
160523: Bug fix for joined curves.
160531-0601: Code modularizations.
160822: Added more Dispose().  Replaced rhinoscriptsyntax functions with RC.
        Name changed from UntrimWithRemoveEdge to DetachLocalTrim.
160826: Updated alias-related code.
        Name changed from DetachLocalTrim to UntrimLocal.
160831: Added removeTrim().  removeTrim() now processes only one face at a time.
171020: Renamed from untrimLocal.py to untrim.py.  Removed alias maker.
180530: Fixed bug that indexed beyond end of list.
180630: Fixed minor bug that only created the curves of one brep when multiple were selected.
190505: Moved some functions from library modules.  Updated imports.
190513-14: Refactored.  Import-related updates.
190620: Corrected printed output of curves added.
200109-10: Import-related update.  Printed feedback change.
        Added Brep.Repair as a temporary fix for some cases where invalid breps are created.
200415: Corrected GeometryAttributeFilter in getInput.
200701: Import-related update.

TODO: 
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xBrepLoop
import xBrepObject
import xBrepTrim
import xSurface

bDebugAdd2dCrvs = False
bDebugAdd3dCrvs = False


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    riAddOpts = {}
    stickyKeys = {}


    def addOptionDouble(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionDouble(
            getObj, englishName=names[key], numberValue=riOpts[key])


    def addOptionInteger(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionInteger(
            getObj, englishName=names[key], intValue=riOpts[key])


    def addOptionList(key, names, listValues, values):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionList(
            getObj,
            englishOptionName=names[key],
            listValues=listValues,
            listCurrentIndex=values[key])


    def addOptionToggle(key, names, riOpts):
        return lambda getObj: ri.Custom.GetBaseClass.AddOptionToggle(
            getObj, englishName=names[key], toggleValue=riOpts[key])


    key = 'bEcho'; keys.append(key)
    values[key] = True
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDebug'; keys.append(key)
    values[key] = False
    names[key] = key[1:]
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    riAddOpts[key] = addOptionToggle(key, names, riOpts)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                values[key] = sc.sticky[stickyKeys[key]]


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
            else:
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get brep trims.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select edge to remove")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.EdgeFilter

    go.GeometryAttributeFilter = (
        ~
        ri.Custom.GeometryAttributeFilter.SeamEdge
        )

    while True:
        Opts.riAddOpts['bEcho'](go)
        Opts.riAddOpts['bDebug'](go)
    
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        if res == ri.GetResult.Object:
            break
        else:
            Opts.setValues()
            Opts.saveSticky()
    
    gBreps0 = []
    idx_rgTrims = []
    
    for idx_Obj in range(go.ObjectCount):
        objRef = go.Object(idx_Obj)
        idx_rgTrim = objRef.GeometryComponentIndex.Index
        gBrep0 = objRef.ObjectId
    
        if not gBrep0 in gBreps0:
            gBreps0.append(gBrep0)
            idx_rgTrims.append([idx_rgTrim])
        else:
            idx_rgTrims[gBreps0.index(gBrep0)].append(idx_rgTrim)
    
    go.Dispose()
    
    return tuple(
            [gBreps0] +
            [idx_rgTrims] +
            [Opts.values[key] for key in Opts.keys]
    )


def formatDistance(fDistance):
    try:
        fDistance = float(fDistance)
    except:
        return "(No deviation provided)"

    if fDistance < 0.001:
        return "{:.2e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def removeTrim(rgBrep0, idx_rgFace, idx_rgLoops, idx_rgTrims_ToRemove, bDebug=False):
    """
    Parameters:
        rgBrep0
        idx_rgFaces (of idx_rgTrims_ToRemove)
        idx_rgLoops (of idx_rgTrims_ToRemove)
        idx_rgTrims_ToRemove
    Returns:
        Success: rgBreps1, rgCrvs1_Joined
        Fail: None
    """


    def indicesOfTrimsNotToRemove(rgLoop, idx_rgTrims_ToRemove):
        # Create list of trims of natural, not natural keep, and sublists of not natural remove.
        idx_rgTs_Alpha = [] # Natural indices in main list; not natural indices in sublists.
        for t, rgT in enumerate(rgLoop.Trims):
            if xBrepTrim.isSenw(rgT):
                idx_rgTs_Alpha.append(rgT.TrimIndex)
            else: # Place indices of contiguous not natural trims to remove in sublists.
                if t == 0:
                    idx_rgTs_Alpha.append([rgT.TrimIndex])
                elif type(idx_rgTs_Alpha[-1]) is not list:
                    idx_rgTs_Alpha.append([rgT.TrimIndex])
                elif xBrepTrim.isStartPointOnSenw(rgT):
                    idx_rgTs_Alpha.append([rgT.TrimIndex])
                else:
                    idx_rgTs_Alpha[-1].append(rgT.TrimIndex)
        
        """ If not natural trim sublist starts at end of list and continues to
        beginning, add beginning sublist to the end of the end sublist."""
        if type(idx_rgTs_Alpha[0]) is list and type(idx_rgTs_Alpha[-1]) is list:
            if not xBrepTrim.isStartPointOnSenw(
                    rgLoop.Brep.Trims[idx_rgTs_Alpha[0][0]]):
                idx_rgTs_Alpha = (idx_rgTs_Alpha[1:-1] + [idx_rgTs_Alpha[-1] +
                        idx_rgTs_Alpha[0]])
        
        """ Flatten not natural trims to keep."""
        idx_rgTs_Omega = [] # Natural and not natural but keep indices in main list; not natural, remove indices in sublists.
        numToKeep = 0
        
        for s in idx_rgTs_Alpha:
            if type(s) is int: idx_rgTs_Omega.append(s)
            else: # s is sublist of not natural trims.
                if any(idx in s for idx in idx_rgTrims_ToRemove):
                    idx_rgTs_Omega.append(s) # To remove
                else: # To keep
                    idx_rgTs_Omega.extend(s)
                    numToKeep += 1
        
        if numToKeep == 0 and not xBrepLoop.hasMultipleTrimsOnAnyNaturalEdges(rgLoop): return
        
        """
        1. Find first trim point that doesn't change (between trims to keep).
            If it doesn't exist, return None to use simple loop creator.
        2. Rotate the list so that this trim is at index 0.
            lines = lines[breakPtIdxs[0]:] + lines[:breakPtIdxs[0]]
        """
        
        if not (type(idx_rgTs_Omega[0]) is int and type(idx_rgTs_Omega[-1]) is int):
            iLen = len(idx_rgTs_Omega)
            for i in range(iLen):
                iPlus1 = (i+1) % iLen
                if (type(idx_rgTs_Omega[i]) is int
                        and type(idx_rgTs_Omega[iPlus1]) is int):
                    idx_rgTs_Omega = idx_rgTs_Omega[iPlus1:] + idx_rgTs_Omega[:iPlus1]
                    break
            else: return # Because no connecting trims to keep found.
        
        # Remove subLists.
        return [s for s in idx_rgTs_Omega if type(s) is int]


    def indexOfExistingVertexPerLocation(pt3d_ToTest, rgBrep, tol=1.e-6):
        for rgVertex in rgBrep.Vertices:
            if pt3d_ToTest.EpsilonEquals(rgVertex.Location, tol):
                return rgVertex.VertexIndex


    def tryExistingEdgeMatchingCurve(rgCrv, rgBrep):
        """
        Returns: rgEdge if vertex indices match and curves have distance deviation within a tolerance.
        """
    
        # First, get existing vertices that match or return None if either don't exist.
        idx_rgV_rgCrvStart = indexOfExistingVertexPerLocation(
                rgCrv.PointAtStart, rgBrep)
        if idx_rgV_rgCrvStart is None: return
        idx_rgV_rgCrvEnd = indexOfExistingVertexPerLocation(
                rgCrv.PointAtEnd, rgBrep)
        if idx_rgV_rgCrvEnd is None: return
    
        for rgE in rgBrep.Edges:
            # First, check whether vertex indices at both ends match.
            if not ((idx_rgV_rgCrvStart == rgE.StartVertex.VertexIndex and
                    idx_rgV_rgCrvEnd == rgE.EndVertex.VertexIndex)
                    or
                    (idx_rgV_rgCrvStart == rgE.EndVertex.VertexIndex and
                    idx_rgV_rgCrvEnd == rgE.StartVertex.VertexIndex)):
                continue
        
            # Now, check whether curves overlap, and if so, their maximum deviation.
            rc = rg.Curve.GetDistancesBetweenCurves(
                    rgCrv, rgE, sc.doc.ModelAbsoluteTolerance)
            if not rc[0] or rc[1] > Rhino.RhinoMath.ZeroTolerance:
                continue
        
            # Now, compare lengths.
            if Rhino.RhinoMath.EpsilonEquals(rgCrv.GetLength(), rgE.GetLength(), 1.0e-7):
                return rgE


    def addEdgeTrimToLoopAsDuplicateFromOtherBrep(rgL_B1, rgT_B0, bDebug=False):
    
        if bDebug: print '_\n' + 'addEdgeTrimToLoopAsDuplicateFromOtherBrep()' + '-'*10 + 'Start'
    
        rgB1 = rgL_B1.Brep
        rgTrims_L_B0 = rgT_B0.Brep.Loops[rgT_B0.Loop.LoopIndex].Trims
    
        fTol = 0.0
    
        rgEdge_B0 = rgT_B0.Edge
        rgCrv3d = rgEdge_B0.Duplicate()
    
        # Search and use existing or add a new edge.
        rgEdge_B1 = tryExistingEdgeMatchingCurve(rgCrv3d, rgB1) # Try all trims in case a previous boundary becomes a seam.
        if rgEdge_B1 is not None:
            bRevCrv3d = True # Because the previously made edge was already aligned with its Trim.
        else: # Trim's edge has not been created.
            rgEdge_B0_Tol = rgEdge_B0.Tolerance
            if rgT_B0.IsReversed(): rgCrv3d.Reverse()
            bRevCrv3d = False
        
            # If this is the 1st trim, add vertex for first trim of brep.
            if rgL_B1.Trims.Count == 0:
                rgVertex_EStart = rgB1.Vertices.Add(
                        rgCrv3d.PointAtStart, fTol)
                if bDebugAdd3dCrvs: rs.AddTextDot('V0', rgVertex_EStart.Location)
            else: # Use vertex matching the end of the last trim in loop.
                rgVertex_EStart = xBrepLoop.endVertexOfLastTrim(rgL_B1)
                if rgVertex_EStart is None:
                    idxVs_NoEdges = indicesOfVerticesWithoutEdges(rgB1)
                    if len(idxVs_NoEdges) > 0:
                        rgVertex_EStart = rgB1.Vertices[idxVs_NoEdges[-1]]
                    else:
                        print 'No existing vertex found for edge!'
                        return
        
            # Search and use existing or add a new vertex at end of curve for edge.
            idx_rgVertex_End = indexOfExistingVertexPerLocation(rgCrv3d.PointAtEnd, rgB1)
            if idx_rgVertex_End is None:
                rgVertex_EEnd = rgB1.Vertices.Add(rgCrv3d.PointAtEnd, fTol)
                idx_rgVertex_End = rgVertex_EEnd.VertexIndex
                if bDebugAdd3dCrvs:
                    rs.AddTextDot('V' + str(idx_rgVertex_End),
                            rgVertex_EEnd.Location)
            else:
                rgVertex_EEnd = rgB1.Vertices[idx_rgVertex_End]
        
            #
            # Add Crv3d to brep
            idx_Crv3d = rgB1.AddEdgeCurve(rgCrv3d)
        
            #
            # Add brep edge
            rgEdge_B1 = rgB1.Edges.Add(
                    rgVertex_EStart.VertexIndex,
                    rgVertex_EEnd.VertexIndex,
                    idx_Crv3d, rgEdge_B0_Tol)
        
            if bDebugAdd3dCrvs:
                sc.doc.Objects.AddCurve(rgEdge_B1)
                rs.AddTextDot('E' + str(rgB1.Edges[rgB1.Edges.Count - 1].EdgeIndex),
                        rgEdge_B1.PointAt(rgEdge_B1.DivideByCount(2, False)[0]))
    
        #
        # Add brep trim.
        rgCrv2d = rgT_B0.ToNurbsCurve()
    
        idx_Crv2d = rgB1.AddTrimCurve(rgCrv2d)
    
        rgTrim_B1 = rgB1.Trims.Add(rgEdge_B1, bRevCrv3d, rgL_B1, idx_Crv2d)
        rgTrim_B1.SetTolerances(fTol, fTol)
        rgTrim_B1.IsoStatus = rgTrim_B1.Face.IsIsoparametric(rgCrv2d)
    
        if bDebugAdd2dCrvs:
            sc.doc.Objects.AddCurve(rgCrv2d)
            rs.AddTextDot('T' + str(rgB1.Trims[rgB1.Trims.Count - 1].TrimIndex),
                    rgCrv2d.PointAt(rgCrv2d.DivideByCount(2, False)[0]))
        if bDebug:
            sEval='rgB1.Trims[rgB1.Trims.Count - 1].TrimIndex'; print sEval+':',eval(sEval)
            print '-'*50 + 'End'
    
        return rgB1


    def isParamPtAlphaBeforeParamPtBetaPerIsoStatusPerLoopDir(rgSrf, isoStatus, pt3d_0, pt3d_1):
        if isoStatus == rg.IsoStatus.South:
            return pt3d_0.X < pt3d_1.X
        elif isoStatus == rg.IsoStatus.East:
            return pt3d_0.Y < pt3d_1.Y
        elif isoStatus == rg.IsoStatus.North:
            return pt3d_0.X > pt3d_1.X
        elif isoStatus == rg.IsoStatus.West:
            return pt3d_0.Y > pt3d_1.Y


    def paramPt3dAtStartOfIsoStatusPerLoopDir(rgSrf, isoStatus):
        if isoStatus == rg.IsoStatus.South:
            return rg.Point3d(
                    rgSrf.Domain(0).Min, rgSrf.Domain(1).Min, 0.0)
        elif isoStatus == rg.IsoStatus.East:
            return rg.Point3d(
                    rgSrf.Domain(0).Max, rgSrf.Domain(1).Min, 0.0)
        elif isoStatus == rg.IsoStatus.North:
            return rg.Point3d(
                    rgSrf.Domain(0).Max, rgSrf.Domain(1).Max, 0.0)
        elif isoStatus == rg.IsoStatus.West:
            return rg.Point3d(
                    rgSrf.Domain(0).Min, rgSrf.Domain(1).Max, 0.0)


    def paramPt3dAtEndOfIsoStatusPerLoopDir(rgSrf, isoStatus):
        if isoStatus == rg.IsoStatus.South:
            return rg.Point3d(
                    rgSrf.Domain(0).Max, rgSrf.Domain(1).Min, 0.0)
        elif isoStatus == rg.IsoStatus.East:
            return rg.Point3d(
                    rgSrf.Domain(0).Max, rgSrf.Domain(1).Max, 0.0)
        elif isoStatus == rg.IsoStatus.North:
            return rg.Point3d(
                    rgSrf.Domain(0).Min, rgSrf.Domain(1).Max, 0.0)
        elif isoStatus == rg.IsoStatus.West:
            return rg.Point3d(
                    rgSrf.Domain(0).Min, rgSrf.Domain(1).Min, 0.0)


    def addSingularTrimToLoop(rgL_B1, isoStatus, ptT_Start=None, ptT_End=None, bDebug=False):
        if bDebug: print '_\n' + 'addSingularTrimToLoop()' + '-'*10 + 'Start'
        
        rgB1 = rgL_B1.Brep
        
        fTol = 0.0
        
        if ptT_Start is None:
            ptT_Start = paramPt3dAtStartOfIsoStatusPerLoopDir(rgL_B1.Face, isoStatus)
        
        if ptT_End is None:
            ptT_End = paramPt3dAtEndOfIsoStatusPerLoopDir(rgL_B1.Face, isoStatus)
        
        pt_V = rgL_B1.Face.PointAt(ptT_Start.X, ptT_Start.Y)
        
        idx_rgV = indexOfExistingVertexPerLocation(pt_V, rgB1)
        if idx_rgV is None:
            rgV = rgB1.Vertices.Add(pt_V, fTol)
            idx_rgV = rgV.VertexIndex
        else:
            rgV = rgB1.Vertices[idx_rgV]
        
        rgCrv2d = rg.LineCurve(ptT_Start, ptT_End)
        
        idx_Crv2d = rgB1.AddTrimCurve(rgCrv2d)
        isoStatus = rgB1.Faces[0].IsIsoparametric(rgCrv2d)
        rgB1.Trims.AddSingularTrim(
                rgV, rgL_B1, isoStatus, idx_Crv2d)
        pt_End = rgV.Location
        
        if bDebugAdd2dCrvs:
            sc.doc.Objects.AddCurve(rgCrv2d)
            rs.AddTextDot('T' + str(rgB1.Trims[rgB1.Trims.Count - 1].TrimIndex),
                    rgCrv2d.PointAt(rgCrv2d.DivideByCount(2, False)[0]))
        if bDebug:
            sEval='rgB1.Trims[rgB1.Trims.Count - 1].TrimIndex'; print sEval+':',eval(sEval)
            print '-'*50 + 'End'
        
        return rgB1


    def greaterOf2ContiguousIsoStatuses(isoStatuses):
        if isoStatuses[0] is None or isoStatuses[1] is None: return
        south = rg.IsoStatus.South
        east = rg.IsoStatus.East
        north = rg.IsoStatus.North
        west = rg.IsoStatus.West
        senw = (rg.IsoStatus.South, rg.IsoStatus.East,
                rg.IsoStatus.North, rg.IsoStatus.West)
        if (senw.index(isoStatuses[0]) - senw.index(isoStatuses[1]) == 1 or
                senw.index(isoStatuses[0]) - senw.index(isoStatuses[1]) == -3):
            return isoStatuses[0]
        return isoStatuses[1]


    def nextIsoStatus(isoStatus):
        return eval('rg.IsoStatus.' + (
                'West', 'South', 'East', 'North')[(isoStatus.value__ - 3 + 1) % 4])


    def addSenwEdgeTrimToLoop(rgL_B1, isoStatus, ptT_Start=None, ptT_End=None, bDebug=False):
        """
        Parameters:
        
        Returns:
        
        """
        if bDebug: print '_\n' + 'addSenwEdgeTrimToLoop()' + '-'*10 + 'Start'
    
        rgB1 = rgL_B1.Brep
    
        fTol = 0.0
    
        if ptT_Start is None:
            ptT_Start = paramPt3dAtStartOfIsoStatusPerLoopDir(rgL_B1.Face, isoStatus)
    
        if ptT_End is None:
            ptT_End = paramPt3dAtEndOfIsoStatusPerLoopDir(rgL_B1.Face, isoStatus)
    
        # Create rgCrv3d.
        rgCrv3d_0 = xSurface.senwCurvePerSide(rgL_B1.Face, isoStatus)
    
        if isoStatus == rg.IsoStatus.South:
            intrvl = rg.Interval(ptT_Start.X, ptT_End.X)
        elif isoStatus == rg.IsoStatus.North:
            intrvl = rg.Interval(ptT_End.X, ptT_Start.X)
        elif isoStatus == rg.IsoStatus.West:
            intrvl = rg.Interval(ptT_End.Y, ptT_Start.Y)
        elif isoStatus == rg.IsoStatus.East:
            intrvl = rg.Interval(ptT_Start.Y, ptT_End.Y)
        else: return
    
        # Trim rgCrv3d to interval.
        rgCrv3d_1 = rgCrv3d_0.Trim(intrvl)
        if rgCrv3d_1 is None:
            if bDebug:
                sEval='rgCrv3d_1'; print sEval+':',eval(sEval)
            intrvl.Swap()
            rgCrv3d_1 = rgCrv3d_0.Trim(intrvl)
            if bDebug:
                sEval='rgCrv3d_1'; print sEval+':',eval(sEval)
                sEval='formatDistance(intrvl.T1 - intrvl.T0)'; print sEval+':',eval(sEval)
            if rgCrv3d_1 is None:
                if bDebug: print "intrvl may be too short."
                return
    
        # If trim is at North or West, reverse rgCrv3d's direction to match that of loop.
        if (isoStatus == rg.IsoStatus.North or
                isoStatus == rg.IsoStatus.West):
            rgCrv3d_1.Reverse()
    
        #
        # Get / add existing edge.
        if xSurface.isIsoStatusAtSeam(rgL_B1.Face, isoStatus):
            rgEdge_B1 = tryExistingEdgeMatchingCurve(rgCrv3d_1, rgB1)
        else: rgEdge_B1 = None
        if rgEdge_B1 is not None:
            bRevCrv3d = True # Because curve for previous edge was aligned with the shared trim.
        else: # Trim's edge has not been created.
            bRevCrv3d = False # Because rgCrv3d was made as an extraction from the surface and already aligned with trim.
        
            # If this is the 1st trim, add vertex for first trim of brep.
            if rgL_B1.Trims.Count == 0:
                rgVertex_EStart = rgB1.Vertices.Add(rgCrv3d_1.PointAtStart, fTol)
                if bDebugAdd3dCrvs: rs.AddTextDot('V0', rgVertex_EStart.Location)
            else: # Get existing vertex.
                rgVertex_EStart = xBrepLoop.endVertexOfLastTrim(rgL_B1)
                if rgVertex_EStart is None:
                    idxVs_NoEdges = indicesOfVerticesWithoutEdges(rgB1)
                    if len(idxVs_NoEdges) > 0:
                        rgVertex_EStart = rgB1.Vertices[idxVs_NoEdges[-1]]
                    else:
                        print 'No existing vertex found for edge!'
                        return
        
            # Get / add end vertex.
            idx_rgVertex_End = indexOfExistingVertexPerLocation(rgCrv3d_1.PointAtEnd, rgB1)
            if idx_rgVertex_End is None:
                rgVertex_EEnd = rgB1.Vertices.Add(rgCrv3d_1.PointAtEnd, fTol)
                idx_rgVertex_End = rgVertex_EEnd.VertexIndex
                if bDebugAdd3dCrvs:
                    rs.AddTextDot('V' + str(idx_rgVertex_End),
                            rgVertex_EEnd.Location)
            else:
                rgVertex_EEnd = rgB1.Vertices[idx_rgVertex_End]
        
            #
            # Add rgCrv3d to brep.
            idx_Crv3d_B1 = rgB1.AddEdgeCurve(rgCrv3d_1)
        
            #
            # Add brep edge
            rgEdge_B1 = rgB1.Edges.Add(
                    rgVertex_EStart.VertexIndex,
                    rgVertex_EEnd.VertexIndex,
                    idx_Crv3d_B1, fTol)
        
            if bDebugAdd3dCrvs:
                sc.doc.Objects.AddCurve(rgEdge_B1)
                rs.AddTextDot('E' + str(rgB1.Edges[rgB1.Edges.Count - 1].EdgeIndex),
                        rgEdge_B1.PointAt(rgEdge_B1.DivideByCount(2, False)[0]))
    
        #
        # Add brep trim.
        rgCrv2d = rg.LineCurve(ptT_Start, ptT_End)
    
        idx_Crv2d = rgB1.AddTrimCurve(rgCrv2d)
    
        rgTrim_B1 = rgB1.Trims.Add(rgEdge_B1, bRevCrv3d,
                rgL_B1, idx_Crv2d)
        rgTrim_B1.SetTolerances(fTol, fTol)
        rgTrim_B1.IsoStatus = rgTrim_B1.Face.IsIsoparametric(rgCrv2d)
    
        if bDebugAdd2dCrvs:
            sc.doc.Objects.AddCurve(rgTrim_B1.ToNurbsCurve())
            rs.AddTextDot('T' + str(rgB1.Trims[rgB1.Trims.Count - 1].TrimIndex),
                    rgCrv2d.PointAt(rgCrv2d.DivideByCount(2, False)[0]))
        if bDebug:
            sEval='rgB1.Trims[rgB1.Trims.Count - 1].TrimIndex'; print sEval+':',eval(sEval)
            print '-'*50 + 'End'
    
        return rgB1


    def addSenwTrimToLoop(rgL_B1, isoStatus, ptT_Start=None, ptT_End=None, bDebug=False):
        if xSurface.isSrfSideSingularPerIsoStatus(rgL_B1.Face, isoStatus):
            return addSingularTrimToLoop(
                    rgL_B1,
                    isoStatus,
                    ptT_Start,
                    ptT_End,
                    bDebug=bDebug)
        else:
            return addSenwEdgeTrimToLoop(
                    rgL_B1,
                    isoStatus,
                    ptT_Start,
                    ptT_End,
                    bDebug=bDebug)


    def addOuterLoop(rgB1, rgLoop_B0, idx_rgTrims_ToRemove, bDebug=False):
        """
        Parameters:
            rgB1
            rgLoop_B0
            idx_rgTrims_ToRemove
            bDebug
        Returns:
            
        """
        
        
        idx_rgTs_B0_NotRemoved = indicesOfTrimsNotToRemove(rgLoop_B0, idx_rgTrims_ToRemove)
        if bDebug: sEval='idx_rgTs_B0_NotRemoved'; print sEval+':',eval(sEval)
        if idx_rgTs_B0_NotRemoved is None:
            rgB1.Loops.AddOuterLoop(0)
            return rgB1
        
        fTrimEndPtTol = 1.e-6
        
        #
        # Add loop to Brep1.
        rgL_B1 = rgB1.Loops.Add(
                rg.BrepLoopType.Outer, rgB1.Faces[0])
        
        #
        # Create vertices, edges, trims.
        
        #    for iT_L_B0, rgT_B0 in enumerate(rgTrims_L_B0):
        for i, idx_rgT_B0_Reduced in enumerate(idx_rgTs_B0_NotRemoved):
            if bDebug:
                print '_'*20
                print 'Trim construction stage:', i
                sEval='idx_rgT_B0_Reduced'; print sEval+':',eval(sEval)
            
            rgT_B0 = rgLoop_B0.Brep.Trims[idx_rgT_B0_Reduced]
            #sEval='rgT_B0.TrimIndex'; print sEval+':',eval(sEval)
            
            rgT_B0_Prev_Keep = rgLoop_B0.Brep.Trims[idx_rgTs_B0_NotRemoved[(i-1) % len(idx_rgTs_B0_NotRemoved)]]
            if i != 0:
                rgT_B1_LastAdded = rgL_B1.Trims[rgL_B1.Trims.Count-1] # Get new one.
                if not rgT_B0.PointAtStart.EpsilonEquals(rgT_B1_LastAdded.PointAtEnd, fTrimEndPtTol):
                    if bDebug: print "This trim has already been added."
                    continue
            
            rgT_B0_Next_Keep = rgLoop_B0.Brep.Trims[idx_rgTs_B0_NotRemoved[(i+1) % len(idx_rgTs_B0_NotRemoved)]]
            
            """ Check whether end point of current trim is equivalent to the
            start point of the following one to keep.  This will also support
            loops that connect at the current trim end point with only one vertex."""
            if rgT_B0.PointAtEnd.EpsilonEquals(
                    rgT_B0_Next_Keep.PointAtStart, fTrimEndPtTol):
                if rgT_B0.TrimType == rg.BrepTrimType.Singular:
                    rgB1 = addSingularTrimToLoop(rgL_B1, rgT_B0.IsoStatus,
                    rgT_B0.PointAtStart, rgT_B0.PointAtEnd, bDebug)
                else:
                    rgB1 = addEdgeTrimToLoopAsDuplicateFromOtherBrep(rgL_B1, rgT_B0, bDebug)
            else: # Trim(s) to remove follow current trim.
                if xBrepTrim.isSenw(rgT_B0):
                    if (rgT_B0_Next_Keep.IsoStatus == rgT_B0.IsoStatus and
                            isParamPtAlphaBeforeParamPtBetaPerIsoStatusPerLoopDir(
                            rgL_B1.Face, rgT_B0.IsoStatus, rgT_B0.PointAtEnd,
                            rgT_B0_Next_Keep.PointAtStart)): # IsoStatus of current and next trim are the same and in that order, so merge.
                        if rgT_B0.TrimType == rg.BrepTrimType.Seam:
                            rgB1 = addSenwTrimToLoop(
                                    rgL_B1,
                                    rgT_B0.IsoStatus,
                                    rgT_B0.PointAtStart,
                                    rgT_B0.PointAtEnd,
                                    bDebug=bDebug)
                        else:
                            rgB1 = addSenwTrimToLoop(
                                    rgL_B1,
                                    rgT_B0.IsoStatus,
                                    rgT_B0.PointAtStart,
                                    rgT_B0_Next_Keep.PointAtEnd,
                                    bDebug=bDebug)
                            continue # Because latest trim added connects with next trim to keep.
                    else: # Next trim to keep is not SENW.
                        isoStats_rgT_B0_Next_Keep_StartPt = xBrepTrim.senwIsoStatusIntersectingTrimPointAtStart(
                                rgT_B0_Next_Keep)
                        
                        # TODO: Fix bug here
                        
                        if len(isoStats_rgT_B0_Next_Keep_StartPt) == 2:
                            isoStat_rgT_B0_Next_Keep_StartPt = greaterOf2ContiguousIsoStatuses(isoStats_rgT_B0_Next_Keep_StartPt)
                        else:
                            isoStat_rgT_B0_Next_Keep_StartPt = isoStats_rgT_B0_Next_Keep_StartPt[0]
                        
                        if (isoStat_rgT_B0_Next_Keep_StartPt is not None and
                                rgT_B0.IsoStatus == isoStat_rgT_B0_Next_Keep_StartPt and
                                isParamPtAlphaBeforeParamPtBetaPerIsoStatusPerLoopDir(
                                rgL_B1.Face, rgT_B0.IsoStatus, rgT_B0.PointAtEnd,
                                rgT_B0_Next_Keep.PointAtStart)):
                            # End trim to next keep trim start point.
                            rgB1 = addSenwTrimToLoop(
                                    rgL_B1,
                                    rgT_B0.IsoStatus,
                                    rgT_B0.PointAtStart,
                                    rgT_B0_Next_Keep.PointAtStart,
                                    bDebug=bDebug)
                            continue # Because latest trim added connects with next trim to keep.
                        else: # Add trim to next surface corner.
                            rgPt3d_T_corner = paramPt3dAtEndOfIsoStatusPerLoopDir(
                                    rgL_B1.Face, rgT_B0.IsoStatus)
                            rgB1 = addSenwTrimToLoop(
                                    rgL_B1,
                                    rgT_B0.IsoStatus,
                                    rgT_B0.PointAtStart,
                                    rgPt3d_T_corner,
                                    bDebug=bDebug)
                else: # Trim is not SENW.
                    if rgT_B0.TrimType == rg.BrepTrimType.Singular:
                        rgB1 = addSingularTrimToLoop(
                                rgL_B1,
                                rgT_B0.IsoStatus,
                                rgT_B0.PointAtStart,
                                rgT_B0.PointAtEnd,
                                bDebug=bDebug)
                    else:
                        rgB1 = addEdgeTrimToLoopAsDuplicateFromOtherBrep(
                                rgL_B1,
                                rgT_B0,
                                bDebug=bDebug)
                
                rgT_B1_LastAdded = rgL_B1.Trims[rgL_B1.Trims.Count-1] # Get trim last added.
                
                if rgT_B1_LastAdded.PointAtEnd.EpsilonEquals(
                        rgT_B0_Next_Keep.PointAtStart, fTrimEndPtTol):
                    continue
                
                #print "Add full length natural trims that removed trims span over."
                
                isoStats_rgT_B1_LastAdded_End = xBrepTrim.senwIsoStatusIntersectingTrimPointAtEnd(rgT_B1_LastAdded)
                if isoStats_rgT_B1_LastAdded_End is None: return
                if len(isoStats_rgT_B1_LastAdded_End) == 2:
                    isoStat_rgT_B1_LastAdded_EndPt = greaterOf2ContiguousIsoStatuses(isoStats_rgT_B1_LastAdded_End)
                else:
                    isoStat_rgT_B1_LastAdded_EndPt = isoStats_rgT_B1_LastAdded_End[0]
                
                """ First connect any missing trims if last trim endpoint is on the
                same IsoStatus as the next trim to keep start point."""
                if (rgT_B0_Next_Keep.IsoStatus == isoStat_rgT_B1_LastAdded_EndPt):
                    #print "Merge to end of next trim to keep."
                    if xSurface.isIsoStatusAtSeam(rgL_B1.Face, rgT_B0_Next_Keep.IsoStatus):
                        rgB1 = addSenwTrimToLoop(
                                rgL_B1,
                                rgT_B0_Next_Keep.IsoStatus,
                                rgT_B1_LastAdded.PointAtEnd,
                                rgT_B0_Next_Keep.PointAtStart,
                                bDebug=bDebug)
                    else:
                        rgB1 = addSenwTrimToLoop(
                                rgL_B1,
                                rgT_B0_Next_Keep.IsoStatus,
                                rgT_B1_LastAdded.PointAtEnd,
                                rgT_B0_Next_Keep.PointAtEnd,
                                bDebug=bDebug)
                        continue # Because latest trim added connects with next trim to keep.
                else:
                    isoStats_rgT_B0_Next_Keep_StartPt = xBrepTrim.senwIsoStatusIntersectingTrimPointAtStart(
                            rgT_B0_Next_Keep)
                    if len(isoStats_rgT_B0_Next_Keep_StartPt) == 2:
                        isoStat_rgT_B0_Next_Keep_StartPt = greaterOf2ContiguousIsoStatuses(isoStats_rgT_B0_Next_Keep_StartPt)
                    else:
                        isoStat_rgT_B0_Next_Keep_StartPt = isoStats_rgT_B0_Next_Keep_StartPt[0]
                    
                    if (isoStat_rgT_B0_Next_Keep_StartPt
                            == isoStat_rgT_B1_LastAdded_EndPt):
                        # Extend to start of next trim to keep.
                        rgB1 = addSenwTrimToLoop(
                                rgL_B1,
                                isoStat_rgT_B0_Next_Keep_StartPt,
                                rgT_B1_LastAdded.PointAtEnd,
                                rgT_B0_Next_Keep.PointAtStart,
                                bDebug=bDebug)
                        continue # Because latest trim added connects with next trim to keep.
                    else:
                        rgPt3d_T_corner = paramPt3dAtEndOfIsoStatusPerLoopDir(
                                rgL_B1.Face, isoStat_rgT_B1_LastAdded_EndPt)
                        rgB1 = addSenwTrimToLoop(
                                rgL_B1,
                                isoStat_rgT_B1_LastAdded_EndPt,
                                rgT_B1_LastAdded.PointAtEnd,
                                rgPt3d_T_corner,
                                bDebug=bDebug)
                
                rgT_B1_LastAdded = rgL_B1.Trims[rgL_B1.Trims.Count-1] # Get trim last added.
                
                #print "If end point of last trim added matches the start point of the next trim to keep, continue to next trim."
                if rgT_B1_LastAdded.PointAtEnd.EpsilonEquals(
                        rgT_B0_Next_Keep.PointAtStart, fTrimEndPtTol): continue
                
                #print "Using greater because if at corner, the previous point comparison should have been True."
                isoStats_rgT_B0_Next_Keep_StartPt = xBrepTrim.senwIsoStatusIntersectingTrimPointAtStart(
                        rgT_B0_Next_Keep)
                if len(isoStats_rgT_B0_Next_Keep_StartPt) == 2:
                    isoStat_rgT_B0_Next_Keep_StartPt = greaterOf2ContiguousIsoStatuses(isoStats_rgT_B0_Next_Keep_StartPt)
                else:
                    isoStat_rgT_B0_Next_Keep_StartPt = isoStats_rgT_B0_Next_Keep_StartPt[0]
                
                while (nextIsoStatus(rgT_B1_LastAdded.IsoStatus) !=
                        isoStat_rgT_B0_Next_Keep_StartPt):
                    sc.escape_test()
                    
                    rgT_B1_LastAdded = rgL_B1.Trims[rgL_B1.Trims.Count-1]
                    rgB1 = addSenwTrimToLoop(rgL_B1,
                            nextIsoStatus(rgT_B1_LastAdded.IsoStatus),
                            None,
                            None,
                            bDebug=bDebug)
                    
                    rgT_B1_LastAdded = rgL_B1.Trims[rgL_B1.Trims.Count-1] # Get new one.
                
                #print "Add trim connecting corner of natural border to next trim."
                if xBrepTrim.isSenw(rgT_B0_Next_Keep):
                    
                    rgB1 = addSenwTrimToLoop(
                            rgL_B1,
                            rgT_B0_Next_Keep.IsoStatus,
                            rgT_B1_LastAdded.PointAtEnd,
                            rgT_B0_Next_Keep.PointAtEnd,
                            bDebug=bDebug)
                else:
                    isoStats_rgT_B0_Next_Keep_StartPt = xBrepTrim.senwIsoStatusIntersectingTrimPointAtStart(
                            rgT_B0_Next_Keep)
                    if len(isoStats_rgT_B0_Next_Keep_StartPt) == 2:
                        isoStat_rgT_B0_Next_Keep_StartPt = greaterOf2ContiguousIsoStatuses(isoStats_rgT_B0_Next_Keep_StartPt)
                    else:
                        isoStat_rgT_B0_Next_Keep_StartPt = isoStats_rgT_B0_Next_Keep_StartPt[0]
                    
                    rgB1 = addSenwTrimToLoop(
                            rgL_B1,
                            isoStat_rgT_B0_Next_Keep_StartPt,
                            rgT_B1_LastAdded.PointAtEnd,
                            rgT_B0_Next_Keep.PointAtStart,
                            bDebug=bDebug)
                
        return rgB1


    def indicesOfNotSenwTrimsInSublists(rgLoop):
        """
        Returns:
            idx_rgTrims_Natural: Flat list of trim indices of natural edges
            idx_rgTrims_NotSenw_SubLoops: List of lists of trim indices of not natural edges
        """
        idx_rgTrims_Natural = []
        idx_rgTrims_NotSenw_SubLoops = []
        
        bLastTrimWasOfNaturalEdge = False
        
        for t, rgTrim in enumerate(rgLoop.Trims):
            if xBrepTrim.isSenw(rgTrim):
                bLastTrimWasOfNaturalEdge = True
            else:
                if (len(idx_rgTrims_NotSenw_SubLoops) == 0 or
                        bLastTrimWasOfNaturalEdge or
                        xBrepTrim.isStartPointOnSenw(rgTrim)):
                    idx_rgTrims_NotSenw_SubLoops.append([]) # Start new NotNaturalEdge sublist.
                idx_rgTrims_NotSenw_SubLoops[-1].append(rgTrim.TrimIndex)
                bLastTrimWasOfNaturalEdge = False
        
        if (len(idx_rgTrims_NotSenw_SubLoops) > 1 and not xBrepTrim.isStartPointOnSenw(
                rgLoop.Brep.Trims[idx_rgTrims_NotSenw_SubLoops[0][0]])):
            idx_rgTrims_NotSenw_SubLoops = (
                    idx_rgTrims_NotSenw_SubLoops[1:-1] +
                    [idx_rgTrims_NotSenw_SubLoops[-1] +
                    idx_rgTrims_NotSenw_SubLoops[0]])
        
        return idx_rgTrims_NotSenw_SubLoops


    def separateSublistsOfAWhetherAnyItemInListB(listA, listB):
        listTrue = []; listFalse = []
        for sublistInA in listA:
            if any(itemInA in sublistInA for
                    itemInA in listB):
                listTrue.append(sublistInA)
            else:
                listFalse.append(sublistInA)
        return listTrue, listFalse



    fTol = 0.0
    
    rgFace0 = rgBrep0.Faces[idx_rgFace]

    # Create brep for new face.
    rgB1 = rg.Brep()
    rgSrf = rgFace0.UnderlyingSurface()
    idx_Srf = rgB1.AddSurface(rgSrf)
    rgFace1 = rgB1.Faces.Add(idx_Srf)
    
    """ Loop through each loop of selected trims to
            Create modified outer loop, if outer loop was selected.
            Aquire loop's non-selected edges."""
    for iL, idx_rgLoop in enumerate(idx_rgLoops):
        rgLoop_B0 = rgBrep0.Loops[idx_rgLoop]
        
        """ If any of its trims are part of idx_rgTrims_ToRemove,
        custom create outer loop. """
        if rgLoop_B0.LoopType == (rg.BrepLoopType.Outer):
            rgB1 = addOuterLoop(
                    rgB1,
                    rgLoop_B0,
                    idx_rgTrims_ToRemove[iL],
                    bDebug=bDebug)
            
            idx_rgTrims_NotSenw_SubLoops = (
                    indicesOfNotSenwTrimsInSublists(rgLoop_B0))
            
            idx_rgTs_NotSenw_Remove_SubLoops, idx_rgTs_NotSenw_Keep_SubLoops = (
                    separateSublistsOfAWhetherAnyItemInListB(
                    idx_rgTrims_NotSenw_SubLoops, idx_rgTrims_ToRemove[iL]))
            
            idx_rgTs_B0_Senw_Keep = [i for sub in idx_rgTs_NotSenw_Keep_SubLoops for i in sub]
            
            rgEdges_InLoopsWithSel = [
                    rgTrim.Edge for rgTrim in rgLoop_B0.Trims if (
                    rgTrim.TrimIndex not in idx_rgTrims_ToRemove[iL]
                    and
                    rgTrim.TrimType != rg.BrepTrimType.Seam
                    and
                    not xBrepTrim.isSenw(rgTrim)
                    and
                    rgTrim.TrimIndex not in idx_rgTs_B0_Senw_Keep
                    )]
        else:
            """ For 3D curve creation, create list of edges
            not selected in loops with selected edges. """
            rgEdges_InLoopsWithSel = [
                    rgTrim.Edge for rgTrim in rgLoop_B0.Trims if (
                    rgTrim.TrimIndex not in idx_rgTrims_ToRemove[iL]
                    and
                    rgTrim.TrimType != rg.BrepTrimType.Seam
                    and
                    not xBrepTrim.isSenw(rgTrim)
                    )]
        # Uncomment to create unjoined curves:
        #map(sc.doc.Objects.AddCurve, rgEdges_InLoopsWithSel)
        
    """ Duplicate Loops with edges not part of selection from rgBrep0 to
        rgB1."""
    
    # Create list of loops that have none of the selected trims.
    idx_rgLoops_B0_NotSel = [i.LoopIndex for
            i in rgBrep0.Faces[idx_rgFace].Loops if
            i.LoopIndex not in idx_rgLoops]
    
    # Create loops and subtopology.
    for idxL_B0_NotSel in idx_rgLoops_B0_NotSel:
        
        # Add loop.
        rgL_B1 = rgB1.Loops.Add(
                rgBrep0.Loops[idxL_B0_NotSel].LoopType, rgFace1)
        
        # Create vertices, edges, trims.
        rgTrims_L_B0 = rgBrep0.Loops[idxL_B0_NotSel].Trims
        
        for rgT_B0 in rgTrims_L_B0:
            if rgT_B0.TrimType == rg.BrepTrimType.Singular:
                rgB1 = addSingularTrimToLoop(
                        rgL_B1,
                        rgT_B0.IsoStatus,
                        rgT_B0.PointAtStart,
                        rgT_B0.PointAtEnd,
                        bDebug=bDebug)
            else: # Trim is not singular.
                rgB1 = addEdgeTrimToLoopAsDuplicateFromOtherBrep(
                        rgL_B1,
                        rgT_B0,
                        bDebug=bDebug)
            # End of Trim loop
        # End of Loop loop

    if rgB1 and not rgB1.IsValid:
        bRepaired = rgB1.Repair(tolerance=sc.doc.ModelAbsoluteTolerance)
        if bDebug:
            sEval='bRepaired'; print sEval+':',eval(sEval)
            # True from Brep.Repair doesn't mean that the Brep is valid.
            sEval='rgB1.IsValid'; print sEval+':',eval(sEval)
        
        if rgB1.IsValid:
            print "Invalid brep was repaired." \
                "  **  Possible bug in xBrepFace_removeTrim"
        else:
            bBrep1IsValid, sLog = rgB1.IsValidWithLog()
            print sLog
            print '-'*40

            #bBrep1IsValid, sLog = rgB1.IsValidTopology()
            #print sLog
            
            rgB1.Dispose()
            return


    if bDebug:
        if rgB1 is None:
            print "rgB1 is None."
        else:
            sEval='rgB1.Loops.Count'; print sEval+':',eval(sEval)
            sEval='rgB1.Trims.Count'; print sEval+':',eval(sEval)
            sEval='rgB1.Edges.Count'; print sEval+':',eval(sEval)
            sEval='rgB1.Vertices.Count'; print sEval+':',eval(sEval)
            #            for t in rgB1.Loops[0].Trims:
            #                print t.Edge.StartVertex.VertexIndex, t.Edge.EndVertex.VertexIndex
            sEval='rgB1.IsValid'; print sEval+':',eval(sEval)
    
    return rgB1, rgEdges_InLoopsWithSel


def processBrepObjects(gBreps0, idx_rgTs_ToRemove_PerBs0, bEcho=True, bDebug=False):
    
    gBreps1_1F_withRemovedTrim_All = []
    gs_CrvsAdded_All = []



    def processBrepObject(gBrep0, idx_rgTs_ToRemove, bEcho=True, bDebug=False):
        #rgBreps1 = [[] for i in range(len(gBreps0))]
        rdBrep0 = sc.doc.Objects.FindId(gBrep0) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(idBrep0)
        rgBrep0 = rdBrep0.BrepGeometry
    
    
        def createOrderedLists_FacesLoopsTrims(rgBrep0, idx_rgTrims):
            """
            Per provided trim indices per brep,
            create ordered lists of their face, loop, and trim indices.
            """
        
            idx_rgFaces = []
            idx_rgLoops = []
            idx_rgTrims_ToRemove = []
        
            for idx_rgTrim in idx_rgTrims:
                rgTrim = rgBrep0.Trims[idx_rgTrim]
                idx_rgFace = rgTrim.Face.FaceIndex
                idx_rgLoop = rgTrim.Loop.LoopIndex
            
                if not idx_rgFace in idx_rgFaces:
                    idx_rgFaces.append(idx_rgFace)
                    idx_rgLoops.append([idx_rgLoop])
                    idx_rgTrims_ToRemove.append([[idx_rgTrim]])
                else:
                    iF = idx_rgFaces.index(idx_rgFace)
                    if not idx_rgLoop in idx_rgLoops[iF]:
                        idx_rgLoops[iF].append(idx_rgLoop)
                        idx_rgTrims_ToRemove[iF].append([idx_rgTrim])
                    else:
                        iL = idx_rgLoops[iF].index(idx_rgLoop)
                        idx_rgTrims_ToRemove[iF][iL].append(idx_rgTrim)
        
            return idx_rgFaces, idx_rgLoops, idx_rgTrims_ToRemove
    
    
        rc = createOrderedLists_FacesLoopsTrims(rgBrep0, idx_rgTs_ToRemove)
        idx_rgFaces_B0, idx_rgLoops_2Nest, idx_rgTrims_ToRemove_3Nest = rc
    
        attr = rdBrep0.Attributes.Duplicate()
    
        rgBs1_1F = []
        rgCrvs1_Joined = []

        # Get new face and curves of each face of selected trims.
        for iF in xrange(len(idx_rgFaces_B0)):

            rc = removeTrim(
                    rgBrep0,
                    idx_rgFaces_B0[iF],
                    idx_rgLoops_2Nest[iF],
                    idx_rgTrims_ToRemove_3Nest[iF],
                    bDebug=bDebug)
            if rc is None: continue

            rgB1_1F, rgEdges_InLoopsWithSel = rc
            if not rgB1_1F: continue

            rgBs1_1F.append(rgB1_1F)

            # For joined curves.
            if len(rgEdges_InLoopsWithSel) > 0:
                rgCrvs1_Joined.extend(rg.Curve.JoinCurves(
                        rgEdges_InLoopsWithSel))

        if len(rgBs1_1F) != len(idx_rgFaces_B0):
            print "Some of the new brep geometry is missing."
            return


        # Modify brep.

        gBreps1_1F_withRemovedTrim = []
    
        for iB in range(len(rgBs1_1F)):
            if not rgBs1_1F[iB].IsValid:

                    continue # to next Face.

            # Brep is valid.

            # Add single-face breps and create list of faces to remove from main brep.
            gBrep1 = sc.doc.Objects.AddBrep(rgBs1_1F[iB], attr)
            if gBrep1 != Guid.Empty:
                gBreps1_1F_withRemovedTrim.append(gBrep1)
            else:
                print "Could not add a face with removed trim for {}.".format(gBrep0)
                continue # to next Face.

            # No breaks (fails), so remove faces from polyface brep.
            gBreps1_B0WithFacesRemoved = xBrepObject.removeFaces(
                    gBrep0,
                    idx_rgFaces_B0)
            if gBreps1_B0WithFacesRemoved is None:
                print "Faces could not be removed from {}.".format(gBrep0)
                return


        # Add curves.
        gCrvs1_Added = []
        if len(rgCrvs1_Joined) == 0:
            print "No curves added."
        else:
            for rgCrv1_Joined in rgCrvs1_Joined:
                gCrv1_Added = sc.doc.Objects.AddCurve(rgCrv1_Joined)
                if gCrv1_Added != Guid.Empty:
                    gCrvs1_Added.append(gCrv1_Added)
            s  = "{}".format(len(gCrvs1_Added))
            s += " curve{} added.".format('s' if len(gCrvs1_Added) != 1 else '')
            if len(rgCrvs1_Joined) != len(gCrvs1_Added):
                s += "  {} curves failed to be added.".format(
                        len(rgCrvs1_Joined) - len(gCrvs1_Added))
            print s

        return gBreps1_1F_withRemovedTrim, gCrvs1_Added



    # Iterate through list of breps.
    for iB0, (gBrep0, idx_rgTs_ToRemove) in enumerate(zip(gBreps0, idx_rgTs_ToRemove_PerBs0)):
        rc = processBrepObject(
                gBrep0=gBrep0,
                idx_rgTs_ToRemove=idx_rgTs_ToRemove,
                bEcho=bEcho,
                bDebug=bDebug)
        if bDebug: print rc
        if rc is None: continue
        gBreps1_1F_withRemovedTrim, gCrvs1_Added = rc
        
        gBreps1_1F_withRemovedTrim_All.append(gBreps1_1F_withRemovedTrim)
        gs_CrvsAdded_All.append(gCrvs1_Added)
    
    return gBreps1_1F_withRemovedTrim_All, gs_CrvsAdded_All


def main():
    
    rc = getInput()
    if rc is None: return

    (
        gBreps0,
        idx_rgTs_ToRemove_PerBs0,
        bEcho,
        bDebug,
    ) = rc
    
    if bDebug:
        print "Running with Debug mode on."
        import sys
        for sModule in list(sys.modules):
            if sModule[0] == 'x':
                try:
                    reload(sys.modules[sModule])
                    print "{} reloaded.".format(sModule)
                except:
                    s  = "{} NOT reloaded.".format(sModule)
                    s += "  Does the module contain a bug"
                    s += " or was its name changed?"
                    print s
    else:
        sc.doc.Views.RedrawEnabled = False

    rc = processBrepObjects(
        gBreps0=gBreps0,
        idx_rgTs_ToRemove_PerBs0=idx_rgTs_ToRemove_PerBs0,
        bEcho=bEcho,
        bDebug=bDebug)

    gBreps1_1F_withRemovedTrim_All, gs_CrvsAdded_All = rc

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()