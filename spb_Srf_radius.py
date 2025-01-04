"""
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
160419: Created.
...
210612: Rewrote the routine in constantRadiusOfSurface that checks curvature/radii at various points.
        constantRadiusOfSurface now returns a negative radius for rounds.
210617: constantRadiusOfSurface now handles surfaces with flat spots.
220613: Now states the face's shape if a cylinder, sphere, or torus.
221108: Now, fRadTol is also used for getting primitive shape.  Refactored.
230210: Now, tolerance for getting primitive shape is maximized at 0.1 * ModelAbsoluteTolerance.
230403,07: Bug fix.
240622: Added NegSign option.
241231: Bug and efficiency fixes.

TODO:
    Obtain face's shape only once.  Refer to 220613 revision.  Another script calls relevant function.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List

import time

import xBlock
import xBrepFace_tryGetPrimitiveShape
import xSurface


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fRadTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'TolForConstRad'
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bPickPt'; keys.append(key)
    values[key] = False
    names[key] = 'PostPickPtForVarRad'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bDualUnits'; keys.append(key)
    values[key] = True
    names[key] = 'DualUnitsForInchOrMm'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bAddDot'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDotHeight'; keys.append(key)
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=3)
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'AddPrefixR'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bNegSign'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'iDecPlaces'; keys.append(key)
    values[key] = sc.doc.ModelDistanceDisplayPrecision - 2
    riOpts[key] = ri.Custom.OptionInteger(values[key], setLowerLimit=True, limit=-1)
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

        if key == 'fRadTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.values[key] = cls.riOpts[key].InitialValue
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]
                return

            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.riOpts:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = cls.riOpts[key].CurrentValue
            return

        if key in cls.listValues:
            sc.sticky[cls.stickyKeys[key]] = cls.values[key] = idxList

        print("Invalid key?")


def closestFaceToPoint(rgBrep, pt3d):
    if rgBrep.Faces.Count == 1: return rgBrep.Faces[0]
    fDist_Min = sc.doc.ModelAbsoluteTolerance * 1000.0
    idxFace_Closest = None
    for f, rgFace in enumerate(rgBrep.Faces):
        bPt, u, v = rgFace.ClosestPoint(pt3d)
        if rgFace.IsPointOnFace(u, v): # Works because values are 0 for Exterior, 1 for Interior, 2 for Boundary
            ptOn = rgFace.PointAt(u, v)
            fDist = pt3d.DistanceTo(ptOn)
            if fDist < fDist_Min:
                idxFace_Closest = f
                fDist_Min = fDist
    return rgBrep.Faces[idxFace_Closest] if idxFace_Closest is not None else None


def getInput():
    """
    Get surface with optional input
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select face")

    go.GeometryFilter = (
        rd.ObjectType.Surface |
        rd.ObjectType.Brep |
        rd.ObjectType.InstanceReference)
    
    go.AcceptNumber(True, acceptZero=True)
    go.EnableHighlight(False)

    #go.SubObjectSelect = False # This was used when RC6's
    # ObjRef.SelectionPoint() was not at the absolute coordinates of the pick.
    # It was instead at a point per the block definition.
    # Fixed in Rhino 6 8/30/16 WIP.
    
    rgFace = None
    
    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while rgFace is None:
        while True:
            go.ClearCommandOptions()

            idxs_Opts.clear()

            addOption('fRadTol')
            addOption('bPickPt')
            addOption('bDualUnits')
            addOption('bAddDot')
            if Opts.values['bAddDot']:
                addOption('iDotHeight')
                addOption('AddPrefixR')
                addOption('bNegSign')
                addOption('iDecPlaces')
            addOption('bEcho')
            addOption('bDebug')

            res = go.Get()

            if res == ri.GetResult.Cancel:
                go.Dispose()
                return

            if res == ri.GetResult.Object:
                break

            if res == ri.GetResult.Number:
                key = 'fRadTol'
                Opts.riOpts[key].CurrentValue = go.Number()
                Opts.setValue(key)
                continue

            # An option was selected.
            for key in idxs_Opts:
                if go.Option().Index == idxs_Opts[key]:
                    Opts.setValue(key, go.Option().CurrentListOptionIndex)
                    break

        objref = go.Object(0)
        # Rhino 6: When sub-object selecting, go.Object(0).Geometry() is the geometry.
        rdObj = objref.Object() # Rhino 6: When selecting a block instance, this is 'ModelGeometry: (unnamed) (0)'.
        # Rhino 6: When selecting a block instance, regardless of SubObjectSelect, rdObj.ObjectType is 'InstanceReference'.
        ptPicked = objref.SelectionPoint() #; sc.doc.Objects.AddPoint(ptPicked)
        
        if go.ObjectsWerePreselected:
            if rdObj.ObjectType == rd.ObjectType.Brep:
                rgObj = rdObj.BrepGeometry
            elif rdObj.ObjectType == rd.ObjectType.Extrusion:
                rgObj = rdObj.Geometry.ToBrep()
            else:
                sc.doc.Objects.UnselectAll() # Necessary if instance reference was preselected.
                sc.doc.Views.Redraw()
                continue
            if rgObj.Faces.Count == 1:
                rgFace = rgObj.Faces[0]
            else:
                sc.doc.Objects.UnselectAll() # Necessary to avoid getting stuck in endless loop if brep has more than 1 face.
                sc.doc.Views.Redraw()
            continue
        rgObj = go.Object(0).Geometry()
        
        if (rdObj.ObjectType == rd.ObjectType.Brep or 
                rdObj.ObjectType == rd.ObjectType.InstanceReference):
            if rgObj.ObjectType == rd.ObjectType.Brep: # For Rhino 6, not 5.
                rgFace = closestFaceToPoint(rgObj, ptPicked)
            elif rgObj.ObjectType == rd.ObjectType.InstanceReference: # For both Rhino 5 & 6.
                rgFace, ptPicked = xBlock.tryPickedFaceOfBlock(rdObj, ptPicked)
        else: rgFace = objref.Face() # This also works for ExtrusionObject.
        
        if rgFace is None: sc.doc.Objects.UnselectAll() # Necessary when go.Get() is repeated.

    go.Dispose()

    return rgFace, ptPicked


def curvaturesAtNormalizedParameters(rgSrf, u_Norm, v_Norm):
    """
    Returns:
        Radius of abs(maximum curvature) (minimum radius)
        Radius of abs(maximum curvature) (minimum radius)
    """
    fDomain_u = rgSrf.Domain(0); fDomain_v = rgSrf.Domain(1)
    fRh0 = Rhino.RhinoMath.ZeroTolerance
    
    kappa_Max = kappa_Min = 0
    u = (fDomain_u[1] - fDomain_u[0]) * u_Norm + fDomain_u[0]
    v = (fDomain_v[1] - fDomain_v[0]) * v_Norm + fDomain_v[0]
    c = rgSrf.CurvatureAt(u, v)
    if c is not None: # None will happen at singularities.
        fK_Max = c.Kappa(0)
        fK_Min = c.Kappa(1)
        if abs(fK_Max) > fRh0: kappa_Max = fK_Max
        if abs(fK_Min) > fRh0: kappa_Min = fK_Min
    return kappa_Max, kappa_Min


def constantRadiusOfSurface(rgFace, fRadTol=None, bAlsoReturnShape=False, bDebug=False):
    """
    Even if surface has 2 constant radii, returns only the average of the minima.
    TODO: In the case of having 2 constant radii, return both,
    and modify code (internal and external) referencing this function.
    """

    if bDebug: print('constantRadiusOfSurface()')
    
    
    def radiusOfPrimitiveShape(rgPrimitive, bDebug=False):
        """
        """
        typeShape = rgPrimitive.GetType()
        if bDebug: print(typeShape)
        if typeShape == rg.Plane:
            return 0.0
        if typeShape == rg.Cylinder:
            rgCir = rgPrimitive.CircleAt(0.)
            return rgCir.Radius
        elif typeShape == rg.Sphere:
            return rgPrimitive.Radius
        elif typeShape == rg.Torus:
            return rgPrimitive.MinorRadius
    
    
    if fRadTol is None: fRadTol = 0.1 * sc.doc.ModelAbsoluteTolerance
    
    fTolerance_Shape = min(fRadTol, 1e-6)
    
    # Frist, check for planarity.
    rc = xBrepFace_tryGetPrimitiveShape.tryGetPrimitiveShape(
        rgFace,
        bMatchToShrunkFace=False,
        bPlane=True,
        bCylinder=False,
        bCone=True,
        bSphere=False,
        bTorus=False,
        fTolerance=1e-9,
        bDebug=bDebug)
    if rc and rc[0]:
        return
    
    # If the face is of a cylinder, sphere, or torus, obtain the shape's radius.
    rc = xBrepFace_tryGetPrimitiveShape.tryGetPrimitiveShape(
        rgFace,
        bMatchToShrunkFace=True,
        bPlane=False,
        bCylinder=True,
        bCone=False,
        bSphere=True,
        bTorus=True,
        fTolerance=fTolerance_Shape,
        bDebug=bDebug)
    if rc and rc[0]:
        rgPrimitive, fTol_Used, sShrunkOrNot = rc[0]

        radius = radiusOfPrimitiveShape(rgPrimitive)
        if radius is not None:
            kappa_Max, kappa_Min = (
                    curvaturesAtNormalizedParameters(rgFace, 0.5, 0.5))
            if rgFace.OrientationIsReversed:
                kappa_Max, kappa_Min = -kappa_Max, -kappa_Min
            if kappa_Max < 0.0:
                if bAlsoReturnShape:
                    return -radius, rgPrimitive
                else:
                    return -radius
            if bAlsoReturnShape:
                return radius, rgPrimitive
            else:
                return radius
    
    #
    # If a single singularity or short senw exists,
    # the face is probably the end of a fillet chain.
    # Only points on the border opposite the (near) singularity will be
    # checked for a constant radius.
    
    iSenwsWithSingularityOrAlmost = [side for side in (0,1,2,3) if
            rgFace.IsSingular(side)]
    rgNurbsSrf = rgFace.ToNurbsSurface()
    iSenwsWithSingularityOrAlmost.extend(
            xSurface.shortSenws(rgNurbsSrf, sc.doc.ModelAbsoluteTolerance)
    )
    rgNurbsSrf.Dispose()
    iSenwsWithSingularityOrAlmost = list(set(iSenwsWithSingularityOrAlmost))
    if len(iSenwsWithSingularityOrAlmost) == 1:
        # Create normalized parameter generator.
        iSenw = iSenwsWithSingularityOrAlmost[0]
        numShort = 2 # Number of columns of samples to take on end opposite of singularity progressing toward singularity.
        numFull = 5 # Number of rows of samples to take on end opposite of singularity traversing parallel with singularity.
        m = 1./(numFull-1) # Division multiplier
        if iSenw == 0:
            #print('South')
            uvs = ((m*u, m*v) for u in xrange(0, numFull) for
                    v in xrange(numFull-numShort, numFull))
        elif iSenw == 1:
            #print('East')
            uvs = ((m*u, m*v) for u in xrange(0, numShort) for
                    v in xrange(0, numFull))
        elif iSenw == 2:
            #print('North')
            uvs = ((m*u, m*v) for u in xrange(0, numFull) for
                    v in xrange(0, numShort))
        elif iSenw == 3:
            #print('West')
            uvs = ((m*u, m*v) for u in xrange(numFull-numShort, numFull) for
                    v in xrange(0, numFull))
        else:
            print("What happened?")
            return
        
        radii_Signed = []
        for u,v in uvs:
    #            print(u, v)
    #            sc.doc.Objects.AddPoint(rgFace.PointAt(
    #                    rgFace.Domain(0).ParameterAt(u),
    #                    rgFace.Domain(1).ParameterAt(v)))
    #            sc.doc.Views.Redraw()
            kappa_Max, kappa_Min = curvaturesAtNormalizedParameters(
                rgFace, u, v)

            if bDebug: sEval = 'kappa_Max'; print(sEval + ':', eval(sEval))
            if kappa_Max == 0: break

            if rgFace.OrientationIsReversed:
                kappa_Max, kappa_Min = -kappa_Max, -kappa_Min

            radii_Signed.append(1.0/kappa_Max)
            # Check whether radii extents are within fRadTol of a mean value.
            if not Rhino.RhinoMath.EpsilonEquals(
                    min(radii_Signed), max(radii_Signed), 2.0*fRadTol):
                break
        else:
            if bAlsoReturnShape:
                return sum(radii_Signed) / float(len(radii_Signed)), None
            else:
                return sum(radii_Signed) / float(len(radii_Signed))
    
    #
    
    #
    # Sample various points across diagonal of surface to check for constant radius.
    
    rh0 = Rhino.RhinoMath.ZeroTolerance
    
    uAndV = 0.5
    if bDebug: sEval='uAndV'; print(sEval+':',eval(sEval))
    
    
    kappa_Max, kappa_Min = curvaturesAtNormalizedParameters(rgFace, uAndV, uAndV)

    if rgFace.OrientationIsReversed:
        kappa_Max, kappa_Min = -kappa_Max, -kappa_Min


    if bDebug:
        sEval='kappa_Max'; print(sEval+':',eval(sEval),)
        sEval='1.0/kappa_Max'; print(sEval+':',eval(sEval))
        sEval='kappa_Min'; print(sEval+':',eval(sEval),)
        sEval='1.0/kappa_Min'; print(sEval+':',eval(sEval))

    if abs(kappa_Max) > rh0:
        radii_Minima_signed = [1.0/kappa_Max]
    else:
        radii_Minima_signed = []

    if abs(kappa_Min) > rh0:
        radii_Maxima_signed = [1.0/kappa_Min]
    else:
        radii_Maxima_signed = []

    #radii_Signed = list(radii_Signed_MidT)

    for uAndV in 0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9:
        # Notice that points along borders of surface are skipped to avoid creating a special case for singularities.
        
        if bDebug:
            sEval='uAndV'; print(sEval+':',eval(sEval))
        
        kappa_Max, kappa_Min = curvaturesAtNormalizedParameters(rgFace, uAndV, uAndV)

        if rgFace.OrientationIsReversed:
            kappa_Max, kappa_Min = -kappa_Max, -kappa_Min

        if bDebug:
            sEval='kappa_Max'; print(sEval+':',eval(sEval),)
            sEval='1.0/kappa_Max'; print(sEval+':',eval(sEval))
            sEval='kappa_Min'; print(sEval+':',eval(sEval),)
            sEval='1.0/kappa_Min'; print(sEval+':',eval(sEval))

        if abs(kappa_Max) > rh0:
            radii_Minima_signed.append(1.0/kappa_Max)

        if abs(kappa_Min) > rh0:
            radii_Maxima_signed.append(1.0/kappa_Min)

        if len(radii_Minima_signed) == 0 and len(radii_Maxima_signed) == 0:
            print("Flat spot detected on surface, so it has no constant radius."
                  " Should the face be split?")
            return

        if abs(max(radii_Minima_signed) - min(radii_Minima_signed)) > fRadTol:
            if radii_Maxima_signed and abs(max(radii_Maxima_signed) - min(radii_Maxima_signed)) > fRadTol:
                # No constant radius for either maxima or minima.
                return

        #        if abs(kappa_Max) > rh0 and abs(abs(1.0/kappa_Max) - abs(radii_Signed_MidT[0])) <= fRadTol:
        #            if abs(kappa_Min) > rh0 and abs(abs(1.0/kappa_Min) - abs(radii_Signed_MidT[-1])) <= fRadTol:
        #                # Both curvatures of current sample point match those at first sample point.
        #                # Record radius with largest magnitude.
        #                if abs(1.0/kappa_Max) > abs(1.0/kappa_Min):
        #                    radii_Signed.append(1.0/kappa_Max)
        #                else:
        #                    radii_Signed.append(1.0/kappa_Min)
        #            else:
        #                # Only kappa_Max matches.
        #                radii_Signed.append(1.0/kappa_Max)
        #        elif abs(kappa_Min) > rh0 and abs(abs(1.0/kappa_Min) - abs(radii_Signed_MidT[-1])) <= fRadTol:
        #            radii_Signed.append(1.0/kappa_Min)
        #        else:
        #            # No curvatures match between first and second sample points.
        #            return

    if bDebug:
        sEval='min(radii_Minima_signed)'; print(sEval+':',eval(sEval))
        sEval='max(radii_Minima_signed)'; print(sEval+':',eval(sEval))
        sEval='min(radii_Maxima_signed)'; print(sEval+':',eval(sEval))
        sEval='max(radii_Maxima_signed)'; print(sEval+':',eval(sEval))
        print

    if abs(max(radii_Minima_signed) - min(radii_Minima_signed)) <= fRadTol:
        if bAlsoReturnShape:
            return sum(radii_Minima_signed) / float(len(radii_Minima_signed)), None
        else:
            return sum(radii_Minima_signed) / float(len(radii_Minima_signed))

    return



    #    if abs(max(radii_Maxima_signed) - min(radii_Maxima_signed)) <= fRadTol:
    #        return sum(radii_Maxima_signed) / float(len(radii_Maxima_signed))

    #    # Check whether radii extents are within fRadTol of midrange.
    #    if not Rhino.RhinoMath.EpsilonEquals(
    #            min(radii_Signed), max(radii_Signed), 2.0*fRadTol
    #    ):
    #        return
    #
    #    return abs(sum(radii_Signed) / len(radii_Signed))


def main():
    
    rc = getInput()
    if rc is None: return
    rgFace0, ptPicked = rc


    fRadTol = Opts.values['fRadTol']
    bPickPt = Opts.values['bPickPt']
    bDualUnits = Opts.values['bDualUnits']
    bAddDot = Opts.values['bAddDot']
    iDotHeight = Opts.values['iDotHeight']
    AddPrefixR = Opts.values['AddPrefixR']
    bNegSign = Opts.values['bNegSign']
    iDecPlaces = Opts.values['iDecPlaces']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    def deleteTempRdObjs(gTemp):
        if not gTemp: return
        idBreps1_NETList = List[Guid](gTemp)
        sc.doc.UndoRecordingEnabled = False
        sc.doc.Objects.Delete(idBreps1_NETList, True)
        sc.doc.UndoRecordingEnabled = True


    bDocWasAlreadyModified = sc.doc.Modified
    
    rgBrep_1Face = rgFace0.DuplicateFace(True) # Used for highlighting and pick point boundary constraint.
    rgFace0.Brep.Dispose()
    rgBrep_1Face.Faces.ShrinkFaces()
    rgFace1_From1FaceB = rgBrep_1Face.Faces[0]

    if not bAddDot and ptPicked[0] != Rhino.RhinoMath.UnsetValue:
        # Alternative
        ## Add and highlight duplicate of face.
        #idBrep = sc.doc.Objects.AddBrep(rgBrep_1Face)
        #rdBrep = sc.doc.Objects.Find(idBrep)
        #rdBrep.Attributes.WireDensity = -1
        #rdBrep.Highlight(True)
        
        # Add and highlight duplicates of face border curves.
        sc.doc.UndoRecordingEnabled = False
        necs = rgBrep_1Face.DuplicateNakedEdgeCurves(nakedOuter=True, nakedInner=True)
        gCrvs_FaceBorders = [sc.doc.Objects.AddCurve(c) for c in necs]
        sc.doc.UndoRecordingEnabled = True
        
        for c in gCrvs_FaceBorders:
            co = sc.doc.Objects.FindId(c) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(c)
            co.Highlight(True),
        
        sc.doc.Views.Redraw()
        time.sleep(0.05)
    else: gCrvs_FaceBorders = None #idBrep = None
    
    bPlanar = rgFace1_From1FaceB.IsPlanar(1e-6)
    
    if bPlanar:
        s = "Surface is planar."
    else:
        fRad_Srf = None
        rv = constantRadiusOfSurface(
            rgFace=rgFace1_From1FaceB,
            fRadTol=fRadTol,
            bAlsoReturnShape=True,
            bDebug=bDebug,
            )
        if rv is not None:
            fRad_Srf, rgPrimitive = rv
            if rgPrimitive is not None:
                sPrimitive = rgPrimitive.GetType().Name
                if sPrimitive == "Cylinder":
                    s = "Cylindrical face"
                elif sPrimitive == "Sphere":
                    s = "Spherical face"
                elif sPrimitive == "Torus":
                    s = "Toric face"
            else:
                s = "Face has a constant radius."
        else:
            if bPickPt or ptPicked[0] == Rhino.RhinoMath.UnsetValue:
                ptPicked = None
                gp = ri.Custom.GetPoint()
                gp.SetCommandPrompt("Point on surface to analyze radius")
                gp.Constrain(rgBrep_1Face, -1, -1, False)
                gp.Get()
                if gp.CommandResult() == Rhino.Commands.Result.Success:
                    ptPicked = gp.Point()
                gp.Dispose()
            
            if ptPicked is not None:
                bPt, u, v = rgFace1_From1FaceB.ClosestPoint(ptPicked)
                if bPt:
                    surface_curvature = rgFace1_From1FaceB.CurvatureAt(u, v)
                    fK_Max = abs(surface_curvature.Kappa(0))
                    if fK_Max > 0.0:
                        fRad_Srf = 1.0 / fK_Max
                        s = "At picked point:"
                    else:
                        s = "Surface is planar at picked point."
            else:
                deleteTempRdObjs(gCrvs_FaceBorders)
                rgBrep_1Face.Dispose()
                return
        
        if fRad_Srf:
            if bAddDot:
                rgDot = rg.TextDot(
                    '{0}{1:.{2}f}'.format(
                        ('','R')[AddPrefixR],
                        fRad_Srf if bNegSign else abs(fRad_Srf),
                        iDecPlaces),
                    location=ptPicked)
                rgDot.FontHeight = iDotHeight
                sc.doc.Objects.AddTextDot(rgDot)
                rgDot.Dispose()

            if fRad_Srf > 0.0:
                s += "  concave (fillet)  (Hence, positive radius.)"
            else:
                s += "  convex (round)  (Hence, negative radius.)"

            if (not bDualUnits or
                (sc.doc.ModelUnitSystem != Rhino.UnitSystem.Inches
                 and
                 sc.doc.ModelUnitSystem != Rhino.UnitSystem.Millimeters)
            ):
                s += "\nDiameter = {0:.{1}f} {2}".format(
                        2.0*fRad_Srf, sc.doc.ModelDistanceDisplayPrecision, sc.doc.ModelUnitSystem)
                s += "\nRadius = {0:.{1}f} {2}".format(
                        fRad_Srf, sc.doc.ModelDistanceDisplayPrecision, sc.doc.ModelUnitSystem)
            elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Inches:
                s += "\nDiameter = {0:.{1}f} inches".format(
                        2.0*fRad_Srf, sc.doc.ModelDistanceDisplayPrecision)
                s += ' [{0:.{1}f} mm]'.format(
                        2.0*fRad_Srf*25.4, sc.doc.ModelDistanceDisplayPrecision-1)
                s += "\nRadius = {0:.{1}f} inches".format(
                        fRad_Srf, sc.doc.ModelDistanceDisplayPrecision)
                s += ' [{0:.{1}f} mm]'.format(
                        fRad_Srf*25.4, sc.doc.ModelDistanceDisplayPrecision-1)
            elif sc.doc.ModelUnitSystem == Rhino.UnitSystem.Millimeters:
                s += "\nDiameter = {0:.{1}f} mm".format(
                        2.0*fRad_Srf, sc.doc.ModelDistanceDisplayPrecision)
                s += ' [{0:.{1}f} inches]'.format(
                        2.0*fRad_Srf/25.4, sc.doc.ModelDistanceDisplayPrecision+1)
                s += "\nRadius = {0:.{1}f} mm".format(
                        fRad_Srf, sc.doc.ModelDistanceDisplayPrecision)
                s += ' [{0:.{1}f} inches]'.format(
                        fRad_Srf/25.4, sc.doc.ModelDistanceDisplayPrecision+1)
    
    print(s)
    
    # Alternative
    #if idBrep is not None: sc.doc.Objects.Delete(idBrep, True)
    
    deleteTempRdObjs(gCrvs_FaceBorders)
    
    sc.doc.Views.Redraw()
    
    if not bDocWasAlreadyModified and not bAddDot: sc.doc.Modified = False
    
    rgBrep_1Face.Dispose()


if __name__ == '__main__': main()