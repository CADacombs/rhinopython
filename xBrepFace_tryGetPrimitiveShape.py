"""
190521,22: Import-related update.
190523-24: Added some functions.
190529: Removed a function.
190530-31: Improved printed output.
190605: Fixed a minor bug.
190617: Improved efficiency of a function.
190619: Added a function.
190721: Argument name change.
190722-24: Refactoring.  Changed output of some functions.
190730: More feedback is now returned from a function.
190731: Removed some debug prints.
190809-10: Corrected the default return value in a function.
190822: Modified a comment.
191118-19: Import-related update.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc

import xPlaneSurface
import xPrimitiveShape


def tryGetPlane(rgFace0, bMatchToShrunkFace=True, fTolerance=sc.doc.ModelAbsoluteTolerance, bDebug=False):
    
    if bMatchToShrunkFace:
        rgBrep_1Face_Shrunk = rgFace0.DuplicateFace(False)
        rgBrep_1Face_Shrunk.Faces.ShrinkFaces()
        rgFace_Shrunk = rgBrep_1Face_Shrunk.Faces[0]
        rgSrfs_ToTry = rgFace0.UnderlyingSurface(), rgFace_Shrunk.UnderlyingSurface()
    else:
        rgSrfs_ToTry = rgFace0.UnderlyingSurface(),
    
    for i_NotShrunk_Shrunk, rgSrf_ToTry in enumerate(rgSrfs_ToTry):
        rc = xPrimitiveShape.Surface.tryGetPlane(
                rgSrf_ToTry,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc[0]:
            if bMatchToShrunkFace: rgBrep_1Face_Shrunk.Dispose()
            rgPlane_Found, fTol_PlaneFound = rc
            if not rgPlane_Found.IsValid:
                return None, "Plane is not valid."
            bShrunkUsed = bool(i_NotShrunk_Shrunk)
            return (rgPlane_Found, fTol_PlaneFound, bShrunkUsed), None
    
    if bMatchToShrunkFace: rgBrep_1Face_Shrunk.Dispose()
    
    return None, "Plane not found within tolerance of {:.2e}.".format(fTolerance)


def tryGetRoundPrimitive(rgFace0, bMatchToShrunkFace=True, bCylinder=True, bCone=True, bSphere=True, bTorus=True, fTolerance=sc.doc.ModelAbsoluteTolerance, bDebug=False):
    """
    Round primitives are Cone, Cylinder, Sphere, and Torus.
    """
    
    if bDebug:
        reload(xPrimitiveShape)
    
    if bMatchToShrunkFace:
        rgBrep_1Face_Shrunk = rgFace0.DuplicateFace(False)
        rgBrep_1Face_Shrunk.Faces.ShrinkFaces()
        rgFace_Shrunk = rgBrep_1Face_Shrunk.Faces[0]
        rgFaces_ToTest = rgFace0, rgFace_Shrunk
        rgSrfs_ToTry = rgFace0.UnderlyingSurface(), rgFace_Shrunk.UnderlyingSurface()
    else:
        rgSrfs_ToTry = rgFace0.UnderlyingSurface(),
    
    for i_NotShrunk_Shrunk, rgSrf_ToTry in enumerate(rgSrfs_ToTry):
        rc = xPrimitiveShape.Surface.tryGetRoundPrimitive(
                rgSrf0=rgSrf_ToTry,
                bCylinder=bCylinder,
                bCone=bCone,
                bSphere=bSphere,
                bTorus=bTorus,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc[0]:
            if bMatchToShrunkFace: rgBrep_1Face_Shrunk.Dispose()
            rgShape_Found, fTol_ShapeFound = rc
            if not rgShape_Found.IsValid:
                return None, "{} is not valid.".format(rgShape_Found.GetType().Name)
            bShrunkUsed = bool(i_NotShrunk_Shrunk)
            return (rgShape_Found, fTol_ShapeFound, bShrunkUsed), None
    
    if bMatchToShrunkFace: rgBrep_1Face_Shrunk.Dispose()
    
    return None, "Round primitive not found within tolerance of {:.2e}.".format(fTolerance)


def tryGetPlaneSurface(rgFace0, bMatchToShrunkFace=True, fTolerance=0.1*sc.doc.ModelAbsoluteTolerance, bDebug=False):
    """
    Return on success:
        PlaneSurface at a size relative to the BrepFace
        fTol_Used
        bShrunkUsed
    Return on fail:
        None
        String: Log statement.
    """
    
    rgSrf_Face0 = rgFace0.UnderlyingSurface()
    
    if isinstance(rgSrf_Face0, rg.PlaneSurface):
        return None, "Surface is already a PlaneSurface."
    
    rc = tryGetPlane(
            rgFace0,
            bMatchToShrunkFace=bMatchToShrunkFace,
            fTolerance=fTolerance,
            bDebug=bDebug)
    if rc[0] is None:
        return rc
    else:
        rgPlane1, fTol_Used, bShrunkUsed = rc[0]
    
    rgSrf1 = xPlaneSurface.createFromPlaneAndObjectSize(
            rgPlane=rgPlane1,
            obj_ForSize=rgFace0,
            bDebug=bDebug)
    if rgSrf1 is None:
        return None, "createFromPlaneAndObjectSize fail."
    
    return (rgSrf1, fTol_Used, bShrunkUsed), None


def tryGetRevSurfaceOfPrimitiveShape(rgFace0, bMatchToShrunkFace=True, bCylinder=True, bCone=True, bSphere=True, bTorus=True, fTolerance=0.1*sc.doc.ModelAbsoluteTolerance, bDebug=False):
    """
    """
    rgSrf0 = rgFace0.UnderlyingSurface()
    
    rc = tryGetRoundPrimitive(
            rgFace0=rgFace0,
            bMatchToShrunkFace=bMatchToShrunkFace,
            bCylinder=bCylinder,
            bCone=bCone,
            bSphere=bSphere,
            bTorus=bTorus,
            fTolerance=fTolerance,
            bDebug=bDebug)
    if rc[0] is None: return None, "No matching primitive shape found."
   
    rgShape1, fTol_PrimitiveMatch, bFoundUsingShunkFace = rc[0]
    
    if isinstance(rgSrf0, rg.RevSurface):
        if isinstance(rgShape1, rg.Cylinder) or isinstance(rgShape1, rg.Cone):
            if rgSrf0.IsClosed(0) or rgSrf0.IsClosed(1):
                return None, "Surface is already a RevSurface of {} shape.".format(rgShape1.GetType().Name)
        
        if isinstance(rgShape1, rg.Sphere) or isinstance(rgShape1, rg.Torus):
            if rgSrf0.IsClosed(0) and rgSrf0.IsClosed(1):
                return None, "Surface is already a RevSurface of {} shape.".format(rgShape1.GetType().Name)
    
    typeShape1 = rgShape1.GetType()
    if not rgShape1.IsValid:
        return None, "{} is not valid.".format(typeShape1)
    
    if typeShape1 == rg.Plane:
        rgPlane = rgShape1
        sType = "Plane"
        rgSrf_Converted = xPlaneSurface.createFromPlaneAndObjectSize(
                rgPlane, rgFace0)
        if rgSrf_Converted is not None:
            return (rgSrf_Converted, fTol_PrimitiveMatch, bFoundUsingShunkFace), None
        else:
            sPrint = 'rgSrf_Converted'
            return None, sPrint + ':', eval(sPrint)
    
    elif typeShape1 == rg.Cylinder:
        rgCyl = rgShape1
        sType = "Cylinder"
        # Fix cylinder that is infinite (has 0 height).
        if not rgCyl.IsFinite:
            if bDebug: print "Making cylinder finite..."
            rgBbox_Face = rgFace0.GetBoundingBox(True)
            fAddLen = 1.1 * rgBbox_Face.Min.DistanceTo(rgBbox_Face.Max)
            rgCyl.Height1 = -fAddLen
            rgCyl.Height2 = fAddLen
        
        
        #    for sAttr in dir(rgCyl):
        #        if (sAttr != 'Unset' and sAttr != '__doc__' and
        #                not callable(getattr(rgCyl, sAttr))):
        #            sPrint = 'rgCyl.' + sAttr; print sPrint + ':', eval(sPrint)
    
    #        sPrint = 'dir(rgCyl)'; print sPrint + ':', eval(sPrint)
        #sPrint = 'rgCyl.Unset'; print sPrint + ':', eval(sPrint)
    #        for sAttr in dir(rgCyl):
    #            if (sAttr != 'Unset' and sAttr != '__doc__' and
    #                    not callable(getattr(rgCyl, sAttr))):
    #                sPrint = 'rgShape1.' + sAttr; print sPrint + ':', eval(sPrint)
        
        elif rgCyl.TotalHeight < 2.0*sc.doc.ModelAbsoluteTolerance:
            if bDebug:
                print "Resultant height of cylinder is {}.".format(
                    rgCyl.TotalHeight)
                print "Using 100 times model's absolute tolerance."
            rgCyl.Height1 = -100.0*sc.doc.ModelAbsoluteTolerance
            rgCyl.Height2 = 100.0*sc.doc.ModelAbsoluteTolerance
        else:
            # Increase length of cylinder so that surface is always larger than trim.
            rgBbox_Face = rgFace0.GetBoundingBox(False)
            rgLines = rgBbox_Face.GetEdges()
            fMaxLen = 0
            for rgLine in rgLines:
                if rgLine.Length > fMaxLen: fMaxLen = rgLine.Length
            fAddLen = 1.1 * fMaxLen
            #fAddLen = 1.1 * rgBbox_Face.Min.DistanceTo(rgBbox_Face.Max)
            rgCyl.Height1 = -fAddLen
            rgCyl.Height2 = fAddLen
        
        rgSrf_Converted = rgCyl.ToRevSurface()
        if rgSrf_Converted is not None:
            return (rgSrf_Converted, fTol_PrimitiveMatch, bFoundUsingShunkFace), None
        else:
            return None, "RevSurface could not be obtained from Cylinder."
        
    elif typeShape1 == rg.Cone:
        rgCone = rgShape1
        sType = "Cone"
        rgSrf_Converted = rgCone.ToRevSurface()
        if rgSrf_Converted is not None:
            return (rgSrf_Converted, fTol_PrimitiveMatch, bFoundUsingShunkFace), None
        else:
            return None, "RevSurface could not be obtained from Cone."
        
        # Increase length of cone so that surface is always larger than trim.
        xf = rg.Transform.Scale(rgCone.ApexPoint, 1.1)
        rgSrf_Converted.Transform(xf)
        
    elif typeShape1 == rg.Sphere:
        rgSphere = rgShape1
        sType = "Sphere"
        rgSrf_Converted = rgSphere.ToRevSurface()
        if rgSrf_Converted is not None:
            return (rgSrf_Converted, fTol_PrimitiveMatch, bFoundUsingShunkFace), None
        else:
            return None, "RevSurface could not be obtained from Sphere."
        
    elif typeShape1 == rg.Torus:
        rgTorus = rgShape1
        sType = "Torus"
        rgSrf_Converted = rgTorus.ToRevSurface()
        if rgSrf_Converted is not None:
            return (rgSrf_Converted, fTol_PrimitiveMatch, bFoundUsingShunkFace), None
        else:
            return None, "RevSurface could not be obtained from Torus."
    
    return None, "What happened?"
    
    s  = "  {} found".format(sType)
    s += " matching {} underlying surface".format(
            "shrunk" if bFoundUsingShunkFace else "full")
    s += " at a tolerance of {}.".format(fTol_PrimitiveMatch)
    
    return (rgSrf_Converted, fTol_PrimitiveMatch, bFoundUsingShunkFace), s


def tryGetPrimitiveShape(rgFace0, bMatchToShrunkFace=True, bPlane=True, bCylinder=True, bCone=True, bSphere=True, bTorus=True, fTolerance=1e-9, bDebug=False):
    """
    Updated in May of 2019.  Based on hasPrimitiveShape.
    
    Returns:
        On Success (regardless if primitive is found):
            tuple:
                (
                    tuple: (rg.Plane or rg.Cylinder or rg.Cone or rg.Sphere or rg.Torus,
                    float of tolerance used to find shape,
                    'shrunk' or 'not shrunk')
                ,
                    str: (Fail log)
                )
        On fail: None, None
    """
    rgSrf_Face0 = rgFace0.UnderlyingSurface()
    

    if bPlane:
        rc = tryGetPlane(
                rgFace0,
                bMatchToShrunkFace=bMatchToShrunkFace,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc[0] is not None:
            return rc

    if bCylinder or bCone or bSphere or bTorus:
        rc = tryGetRoundPrimitive(
                rgFace0=rgFace0,
                bMatchToShrunkFace=bMatchToShrunkFace,
                bCylinder=bCylinder,
                bCone=bCone,
                bSphere=bSphere,
                bTorus=bTorus,
                fTolerance=fTolerance,
                bDebug=bDebug)
        if rc[0] is not None:
            return rc
    
    return None, None

