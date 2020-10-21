"""
160419: Created.
160427: Now supports selecting faces in block instances.
160428: Now will select shaded faces not only on isocurves or edges within block instances.
        Now text dot font height is a command option.
        Now all block objects' layers and their parents are tested for visibility.
160501: Now the correct pick point is calculated when shaded faces of breps in block instances are picked.
        This is a workaround for a bug in Rhino.DocObjects.ObjRef.SelectionPoint(), 
        which currently (Rhino 5 & 6) selects a point on the face of the first brep in the block.
        It also selects breps not visible based on the state of the objects' layer or layers' parents.
160504: Bug fix.
160601: Modularizations.
160618: More modularizations and moved some functions to a library.
160625: Fixed minor bug in face selection.
160702: Fixed bug in face selection.  Upgraded sticky code.
160714: Added sc.doc.UndoRecordingEnabled so that the undo buffer doesn't record during the highlighting occurrence, which temporarily adds objects to the model.
160721-22: Rewrote constantRadiusOfSurface(), simplifying the point sampling and now:
            First checks for a primitive shape from which the radius can be directly obtained.
            Then checks the for surfaces that end at a singularity and samples radii at oppsosite side.
            Only records minimum radii (from maximum curvature).
160819: Updated constantRadiusOfSurfacX() for revision to output of primitiveShapeTypXOfFace().
160823: For subobject filter new in Rhino 6:
            Added closestFaceToPoint()
            Modified conditional testing ObjectType of Geometry vs. RhinoObject.
        Replaced rhinoscriptsyntax functions with RC.
160824: go.SubObjectSelect = False added due to SelectionPoint() not being at
        actual picked point during sub-object selection.
160901: Now allows sub-object selection after bug fix in RC6 WIP 8/30/16.
160919: Fixed bug that locked script in endless loop when a mult-face brep was preselected.
161126: Changed default output decimal places from ModelDistanceDisplayPrecision-1 to ModelDistanceDisplayPrecision.
170123: Added sc.doc.Modified code so that document doesn't show as being
        modified due to only this script when text dots are not added.
        Commented out alias maker.
170525: Added command option for whether "R" is the prefix in the dot.
180427: Changed default dot size from 14 to 11, and 'R' prefix is added by default.
181220: closestFaceToPoint no longer crashes when ClosestPoint only finds points Exterior per IsPointOnFace.
190129: Renamed from surfaceRadius.py to Surface_Radius.py.
        Moved constantRadiusOfSurface and curvaturesAtNormalizedParameters from a module.  Added relevant imports.
        WIP: Modified constantRadiusOfSurface and curvaturesAtNormalizedParameters for improved saddle surface support.
190505: Imported a function.  Updated an import name.
190524,29: Import-related updates.
190703: Added second units when ModelUnitSystem is Inches or Millimeters.
190717: Added Opts.  Added bDualUnits.
190722-23: Import-related updates.
190810, 0903: Fixed bug in constantRadiusOfSurface.
191118: Import-related update.
191206: Bug fix.
200518: For accuracy, now uses 0.1*sc.doc.ModelAbsoluteTolerance as a tolerance to recognize primitive shapes.
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


sOpts = (
        'fRadTol',
        'bPickPt',
        'bDualUnits',
        'bAddDot',
        'iDotHeight',
        'bAddR',
        'iDecPlaces',
        'bEcho',
        'bDebug',
)

class Opts():
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    
    for key in sOpts:
        keys.append(key)
        names[key] = key[1:] # Overwrite as wanted in the following.
    
    key = 'fRadTol'
    values[key] = 10.0 * sc.doc.ModelAbsoluteTolerance
    names[key] = 'ToleranceToRecognizeConstantRadius'
    riOpts[key] = ri.Custom.OptionDouble(initialValue=values[key], setLowerLimit=True, limit=0.0)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bPickPt'
    values[key] = False
    names[key] = 'PostPickPtForVarRadius'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDualUnits'
    values[key] = True
    names[key] = 'DualUnitsForInchOrMm'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddDot'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iDotHeight'
    values[key] = 11
    riOpts[key] = ri.Custom.OptionInteger(
            values[key], setLowerLimit=True, limit=3)
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bAddR'
    values[key] = True
    names[key] = 'AddPrefixR'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'iDecPlaces'
    values[key] = sc.doc.ModelDistanceDisplayPrecision - 2
    names[key] = 'DecimalPlaces'
    riOpts[key] = ri.Custom.OptionInteger(
            values[key], setLowerLimit=True, limit=-1)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bEcho'
    values[key] = True
    names[key] = 'Echo'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    values[key] = False
    names[key] = 'Debug'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]
            else:
                # For OptionList.
                values[key] = sc.sticky[stickyKeys[key]]
    
    
    @classmethod
    def setValues(cls):
        for key in sOpts:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                pass
    
    
    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue
            else:
                # For OptionList.
                sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def dispose(listDispose): map(lambda x: x.Dispose, listDispose)


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
    
    disposeUs = []
    go = ri.Custom.GetObject()
    
    go.SetCommandPrompt("Select face")
    
    go.GeometryFilter = (rd.ObjectType.Surface |
            rd.ObjectType.Brep |
            rd.ObjectType.InstanceReference)
    
    go.AcceptNumber(True, True)
    go.EnableHighlight(False)
    #
    
    #go.SubObjectSelect = False # This was used when RC6's
    # ObjRef.SelectionPoint() was not at the absolute coordinates of the pick.
    # It was instead at a point per the block definition.
    # Fixed in Rhino 6 8/30/16 WIP.
    
    rgFace = None
    
    while rgFace is None:
        while True:
            go.AddOptionDouble(Opts.names['fRadTol'], Opts.riOpts['fRadTol'])
            go.AddOptionToggle(Opts.names['bPickPt'], Opts.riOpts['bPickPt'])
            go.AddOptionToggle(Opts.names['bDualUnits'], Opts.riOpts['bDualUnits'])
            go.AddOptionToggle(Opts.names['bAddDot'], Opts.riOpts['bAddDot'])
            if Opts.values['bAddDot']:
                go.AddOptionInteger(Opts.names['iDotHeight'], Opts.riOpts['iDotHeight'])
                go.AddOptionToggle(Opts.names['bAddR'], Opts.riOpts['bAddR'])
                go.AddOptionInteger(Opts.names['iDecPlaces'], Opts.riOpts['iDecPlaces'])
            go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
            go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
            res = go.Get()
            
            if res == ri.GetResult.Cancel:
                go.Dispose()
                return
            elif res == ri.GetResult.Object:
                break
            elif res == ri.GetResult.Number:
                Opts.riOpts['iDecPlaces'].CurrentValue = go.Number()
            
            Opts.setValues()
            Opts.saveSticky()
            go.ClearCommandOptions()
        
        objref = go.Object(0); disposeUs.append(objref)
        # Rhino 6: When sub-object selecting, go.Object(0).Geometry() is the geometry.
        rdObj = objref.Object() # Rhino 6: When selecting a block instance, this is 'ModelGeometry: (unnamed) (0)'.
        # Rhino 6: When selecting a block instance, regardless of SubObjectSelect, rdObj.ObjectType is 'InstanceReference'.
        ptPicked = objref.SelectionPoint() #; sc.doc.Objects.AddPoint(ptPicked)
        
        if go.ObjectsWerePreselected:
            if rdObj.ObjectType == rd.ObjectType.Brep:
                rgObj = rdObj.BrepGeometry; disposeUs.append(rgObj)
            elif rdObj.ObjectType == rd.ObjectType.Extrusion:
                rgObj = rdObj.Geometry.ToBrep(); disposeUs.append(rgObj)
            else:
                sc.doc.Objects.UnselectAll() # Necessary if instance reference was preselected.
                sc.doc.Views.Redraw()
                continue
            if rgObj.Faces.Count == 1: rgFace = rgObj.Faces[0]
            else:
                sc.doc.Objects.UnselectAll() # Necessary to avoid getting stuck in endless loop if brep has more than 1 face.
                sc.doc.Views.Redraw()
            continue
        rgObj = go.Object(0).Geometry(); disposeUs.append(rgObj)
        
        if (rdObj.ObjectType == rd.ObjectType.Brep or 
                rdObj.ObjectType == rd.ObjectType.InstanceReference):
            if rgObj.ObjectType == rd.ObjectType.Brep: # For Rhino 6, not 5.
                rgFace = closestFaceToPoint(rgObj, ptPicked)
            elif rgObj.ObjectType == rd.ObjectType.InstanceReference: # For both Rhino 5 & 6.
                rgFace, ptPicked = xBlock.tryPickedFaceOfBlock(rdObj, ptPicked)
        else: rgFace = objref.Face() # This also works for ExtrusionObject.
        
        if rgFace is None: sc.doc.Objects.UnselectAll() # Necessary when go.Get() is repeated.
    
    go.Dispose()
    
    return tuple([rgFace, ptPicked] + [Opts.values[key] for key in sOpts])


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


def constantRadiusOfSurface(rgFace, fRadTol=None, bEcho=False, bDebug=False):
    if bDebug: print 'constantRadiusOfSurface()'
    
    
    def radiusOfPrimitiveShape(rgPrimitive, bDebug=False):
        """
        """
        typeShape = rgPrimitive.GetType()
        if bDebug: print typeShape
        if typeShape == rg.Plane:
            return 0.
        if typeShape == rg.Cylinder:
            rgCir = rgPrimitive.CircleAt(0.)
            return rgCir.Radius
        elif typeShape == rg.Sphere:
            return rgPrimitive.Radius
        elif typeShape == rg.Torus:
            return rgPrimitive.MinorRadius
    
    
    if fRadTol is None:
        fRadTol=10.0*sc.doc.ModelAbsoluteTolerance
    
    # If the face is of a cylinder, sphere, or torus, obtain the shape's radius.
    rc = xBrepFace_tryGetPrimitiveShape.tryGetPrimitiveShape(
            rgFace,
            bPlane=False,
            bCylinder=True,
            bCone=False,
            bSphere=True,
            bTorus=True,
            fTolerance=0.1*sc.doc.ModelAbsoluteTolerance,
            bDebug=bDebug)
    if rc and rc[0]:
        rgPrimitive, fTol_Used, sShrunkOrNot = rc[0]
        radius = radiusOfPrimitiveShape(rgPrimitive)
        if radius is not None: return radius
    
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
            #print 'South'
            uvs = ((m*u, m*v) for u in xrange(0, numFull) for
                    v in xrange(numFull-numShort, numFull))
        elif iSenw == 1:
            #print 'East'
            uvs = ((m*u, m*v) for u in xrange(0, numShort) for
                    v in xrange(0, numFull))
        elif iSenw == 2:
            #print 'North'
            uvs = ((m*u, m*v) for u in xrange(0, numFull) for
                    v in xrange(0, numShort))
        elif iSenw == 3:
            #print 'West'
            uvs = ((m*u, m*v) for u in xrange(numFull-numShort, numFull) for
                    v in xrange(0, numFull))
        else:
            print "What happened?"
            return
        
        radii_Signed = []
        for u,v in uvs:
    #            print u, v
    #            sc.doc.Objects.AddPoint(rgFace.PointAt(
    #                    rgFace.Domain(0).ParameterAt(u),
    #                    rgFace.Domain(1).ParameterAt(v)))
    #            sc.doc.Views.Redraw()
            kappa_Max, kappa_Min = (
                    curvaturesAtNormalizedParameters(rgFace, u, v)
            )
            if bDebug: sEval = 'kappa_Max'; print sEval + ':', eval(sEval)
            if kappa_Max == 0: break
            radii_Signed.append(1.0/kappa_Max)
            # Check whether radii extents are within fRadTol of midrange.
            if not Rhino.RhinoMath.EpsilonEquals(
                    min(radii_Signed), max(radii_Signed), 2.0*fRadTol):
                break
        else:
            return abs(sum(radii_Signed) / len(radii_Signed))
    
    #
    
    #
    # Sample various points across diagonal of surface to check for constant radius.
    
    rh0 = Rhino.RhinoMath.ZeroTolerance
    
    numFull = 11 # Number of samples across surface if none are skipped.
    m = 1.0/(numFull-1) # Division multiplier.
    if bDebug:
        sEval='numFull'; print sEval+':',eval(sEval)
        sEval='m'; print sEval+':',eval(sEval)
    
    kappa_MaxX, kappa_MinX = curvaturesAtNormalizedParameters(rgFace, m, m)
    if bDebug:
        sEval='kappa_MaxX'; print sEval+':',eval(sEval)
        sEval='kappa_MinX'; print sEval+':',eval(sEval)
    if kappa_MaxX <= rh0 and kappa_MinX <= rh0: return # Flat spot.
    if kappa_MaxX == 0.0:
        radii_Signed_1 = 1.0/kappa_MinX,
    elif kappa_MinX == 0.0:
        radii_Signed_1 = 1.0/kappa_MaxX,
    else:
        radii_Signed_1 = 1.0/kappa_MaxX, 1.0/kappa_MinX
    
    radii_Signed = []
    for i in range(2, numFull-1): # Note that points along borders of surface are skipped to avoid creating a special case for singularities.
        kappa_MaxX, kappa_MinX = curvaturesAtNormalizedParameters(rgFace, m*i, m*i)
        if bDebug:
            sEval='kappa_MaxX'; print sEval+':',eval(sEval)
            sEval='kappa_MinX'; print sEval+':',eval(sEval)
        if kappa_MaxX <= rh0 and kappa_MinX <= rh0: return # Flat spot.
        if kappa_MaxX > rh0 and abs(1.0/kappa_MaxX - radii_Signed_1[0]) <= fRadTol:
            if kappa_MinX > rh0 and abs(1.0/kappa_MinX - radii_Signed_1[-1]) <= fRadTol:
                # Both curvatures of current sample point match those at first sample point.
                # Record radius with largest magnitude.
                if abs(1.0/kappa_MaxX) > abs(1.0/kappa_MinX):
                    if i == 2: radii_Signed.append(radii_Signed_1[0])
                    radii_Signed.append(1.0/kappa_MaxX)
                else:
                    if i == 2: radii_Signed.append(radii_Signed_1[-1])
                    radii_Signed.append(1.0/kappa_MinX)
            else:
                # Only kappa_Max matches.
                if i == 2: radii_Signed.append(radii_Signed_1[0])
                radii_Signed.append(1.0/kappa_MaxX)
        elif kappa_MinX > rh0 and abs(1.0/kappa_MinX - radii_Signed_1[-1]) <= fRadTol:
            if i == 2: radii_Signed.append(radii_Signed_1[-1])
            radii_Signed.append(1.0/kappa_MinX)
        else:
            # No curvatures match between first and second sample points.
            return
        
    # Check whether radii extents are within fRadTol of midrange.
    if not Rhino.RhinoMath.EpsilonEquals(
            min(radii_Signed), max(radii_Signed), 2.0*fRadTol):
        return
    
    return abs(sum(radii_Signed) / len(radii_Signed))


def main():
    
    rc = getInput()
    if rc is None: return
    rgFace0, ptPicked = rc[:2]
    for key, value in zip(sOpts, rc[2:]):
        exec("{} = {}".format(key, value))
    
    disposeUs = []
    
    bDocWasAlreadyModified = sc.doc.Modified
    
    rgBrep_1Face = rgFace0.DuplicateFace(True) # Used for highlighting and pick point boundary constraint.
    rgFace0.Brep.Dispose()
    rgBrep_1Face.Faces.ShrinkFaces()
    rgFace1_From1FaceB = rgBrep_1Face.Faces[0]
    
    if not bAddDot:
        # Alternative
        ## Add and highlight duplicate of face.
        #idBrep = sc.doc.Objects.AddBrep(rgBrep_1Face)
        #rdBrep = sc.doc.Objects.Find(idBrep)
        #rdBrep.Attributes.WireDensity = -1
        #rdBrep.Highlight(True)
        
        # Add and highlight duplicates of face border curves.
        sc.doc.UndoRecordingEnabled = False
        idCrvs_FaceBorders = [sc.doc.Objects.AddCurve(c) for
                c in rgBrep_1Face.DuplicateNakedEdgeCurves(True, True)]
        sc.doc.UndoRecordingEnabled = True
        
        for c in idCrvs_FaceBorders:
            co = sc.doc.Objects.FindId(c) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(c)
            co.Highlight(True),
        
        sc.doc.Views.Redraw()
        time.sleep(0.05)
    else: idCrvs_FaceBorders = None #idBrep = None
    
    bPlanar = rgFace1_From1FaceB.IsPlanar(0.1*sc.doc.ModelAbsoluteTolerance)
    
    if bPlanar:
        s = "Surface is planar."
    else:
        sRadDescr = None
        fRad_Srf = constantRadiusOfSurface(
                rgFace=rgFace1_From1FaceB,
                fRadTol=fRadTol,
                bEcho=bEcho,
                bDebug=bDebug,
        )
        if fRad_Srf is not None:
            s = "Face has constant radius."
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
                    if fK_Max > sc.doc.ModelAbsoluteTolerance:
                        fRad_Srf = 1.0 / fK_Max
                        s = "Diameter, radius at picked point"
                    else:
                        s = "Surface is planar at picked point."
            else:
                rgBrep_1Face.Dispose()
                return
        
        if fRad_Srf:
            if bAddDot:
                rgDot = rg.TextDot('{0}{1:.{2}f}'.format(
                        ('','R')[bAddR], fRad_Srf, iDecPlaces), ptPicked)
                rgDot.FontHeight = iDotHeight
                sc.doc.Objects.AddTextDot(rgDot)
                rgDot.Dispose()

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
    
    print s
    
    # Alternative
    #if idBrep is not None: sc.doc.Objects.Delete(idBrep, True)
    
    if idCrvs_FaceBorders is not None:
        idBreps1_NETList = List[Guid](idCrvs_FaceBorders)
        sc.doc.UndoRecordingEnabled = False
        sc.doc.Objects.Delete(idBreps1_NETList, True)
        sc.doc.UndoRecordingEnabled = True
    
    sc.doc.Views.Redraw()
    
    if not bDocWasAlreadyModified and not bAddDot: sc.doc.Modified = False
    
    rgBrep_1Face.Dispose()
    
    dispose(disposeUs)


if __name__ == '__main__': main()