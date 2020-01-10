"""
180530: Due to failures in senwIsoStatusIntersectingTrimPointAtStart, changed fRh0Tol (in all functions) from Rhino.RhinoMath.ZeroTolerance to 0.1 * sc.doc.ModelAbsoluteTolerance.
190413: Created by extracting functions from xBrep.  Renamed a function.
190513: Moved a couple of functions to only script that currently uses them.
"""

import Rhino
import scriptcontext as sc


def endVertexOfEdgeTrim(rgTrim):
    if rgTrim.IsReversed():
        if rgTrim.Edge.StartVertex is not None:
            return rgTrim.Edge.StartVertex
    else:
        if rgTrim.Edge.EndVertex is not None:
            return rgTrim.Edge.EndVertex


def isSenw(rgTrim):
    return (rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.South or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.East or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.North or
            rgTrim.IsoStatus == Rhino.Geometry.IsoStatus.West)


def isStartPointOnSenw(rgTrim):
    pt = rgTrim.PointAtStart
    face = rgTrim.Face
    #fRh0Tol = Rhino.RhinoMath.ZeroTolerance
    fRh0Tol = 0.1 * sc.doc.ModelAbsoluteTolerance
    return (Rhino.RhinoMath.EpsilonEquals(pt.X, face.Domain(0).Min, fRh0Tol) or
            Rhino.RhinoMath.EpsilonEquals(pt.X, face.Domain(0).Max, fRh0Tol) or
            Rhino.RhinoMath.EpsilonEquals(pt.Y, face.Domain(1).Min, fRh0Tol) or
            Rhino.RhinoMath.EpsilonEquals(pt.Y, face.Domain(1).Max, fRh0Tol))
    return r


def previousNonSingularTrim(rgTrims_L, iTrim_Start, iTrim_Now):
    """
    Find previous non-singular trim, using recursion if necessary.
    """
    iTrim_Now = (iTrim_Now-1) % rgTrims_L.Count
    if iTrim_Now == iTrim_Start: return None
    if rgTrims_L[iTrim_Now].TrimType == Rhino.Geometry.BrepTrimType.Singular:
        return previousNonSingularTrim(rgTrims_L, iTrim_Start, iTrim_Now)
    else: return rgTrims_L[iTrim_Now]


def senwIsoStatusIntersectingTrimPointAtEnd(rgTrim):
    """
    Returns: Tuple of values found or None.
    """
    pt = rgTrim.PointAtEnd
    face = rgTrim.Face
    #fRh0Tol = Rhino.RhinoMath.ZeroTolerance
    fRh0Tol = 0.1 * sc.doc.ModelAbsoluteTolerance
    isoStatus = None
    domainV = face.Domain(1)
    if Rhino.RhinoMath.EpsilonEquals(pt.Y, domainV.Min, fRh0Tol):
        isoStatus = Rhino.Geometry.IsoStatus.South
    if Rhino.RhinoMath.EpsilonEquals(pt.Y, domainV.Max, fRh0Tol):
        isoStatus = Rhino.Geometry.IsoStatus.North
    
    domainU = face.Domain(0)
    if Rhino.RhinoMath.EpsilonEquals(pt.X, domainU.Max, fRh0Tol):
        if isoStatus is None: return (Rhino.Geometry.IsoStatus.East, )
        else: return (isoStatus, Rhino.Geometry.IsoStatus.East)
    if Rhino.RhinoMath.EpsilonEquals(pt.X, domainU.Min, fRh0Tol):
        if isoStatus is None: return (Rhino.Geometry.IsoStatus.West, )
        else: return (isoStatus, Rhino.Geometry.IsoStatus.West)
    if isoStatus is None: return
    else: return (isoStatus, )


def senwIsoStatusIntersectingTrimPointAtStart(rgTrim):
    """
    Returns: Tuple of values found or None.
    """
    pt = rgTrim.PointAtStart
    face = rgTrim.Face
    #fRh0Tol = Rhino.RhinoMath.ZeroTolerance
    fRh0Tol = 0.1 * sc.doc.ModelAbsoluteTolerance
    isoStatus = None
    domainV = face.Domain(1)
    if Rhino.RhinoMath.EpsilonEquals(pt.Y, domainV.Min, fRh0Tol):
        isoStatus = Rhino.Geometry.IsoStatus.South
    if Rhino.RhinoMath.EpsilonEquals(pt.Y, domainV.Max, fRh0Tol):
        isoStatus = Rhino.Geometry.IsoStatus.North
    
    domainU = face.Domain(0)
    if Rhino.RhinoMath.EpsilonEquals(pt.X, domainU.Max, fRh0Tol):
        if isoStatus is None: return (Rhino.Geometry.IsoStatus.East, )
        else: return (isoStatus, Rhino.Geometry.IsoStatus.East)
    if Rhino.RhinoMath.EpsilonEquals(pt.X, domainU.Min, fRh0Tol):
        if isoStatus is None: return (Rhino.Geometry.IsoStatus.West, )
        else: return (isoStatus, Rhino.Geometry.IsoStatus.West)
    if isoStatus is None: return
    else: return (isoStatus, )


def startVertexOfEdgeTrim(rgTrim):
    if rgTrim.IsReversed():
        if rgTrim.Edge.EndVertex is not None:
            return rgTrim.Edge.EndVertex
    else:
        if rgTrim.Edge.StartVertex is not None:
            return rgTrim.Edge.StartVertex
