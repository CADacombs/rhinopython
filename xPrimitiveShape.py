"""
"""

"""
190522-24: Started from another module.
...
220831: Modified tolerance range in shape searching.
241231: Bug fix in tolerances in tolerance loop after cylinder search.
        Added check to quickly eliminate some iterative searches for shapes.

Notes: Surface.IsTorus, etc., has trouble finding shapes for NurbSurfaces converted from RevSurfaces.
"""

import Rhino
import Rhino.Geometry as rg
import scriptcontext as sc


MAX_ACCEPTED_SIZE = 1e4


class AnyShape():

    @staticmethod
    def areEqual(shapes, epsilon=1e-9, bDebug=False):
        """
        """
    
        if len(shapes) < 2: return
    
        shapeA = shapes[0]
        sShapeA = shapeA.GetType().Name
    
        for i in range(1, len(shapes)):
            shapeB = shapes[i]
            sShapeB = shapeB.GetType().Name
        
            if sShapeA != sShapeB: 
                return False
        
            if sShapeA == 'Plane':
                if Plane.areEqual([shapeA, shapeB], epsilon=epsilon, bDebug=bDebug): continue
                else: return False
            elif sShapeA == 'Cylinder':
                if Cylinder.areEqual([shapeA, shapeB], epsilon=epsilon, bDebug=bDebug): continue
                else: return False
            elif sShapeA == "Cone":
                if Cone.areEqual([shapeA, shapeB], epsilon=epsilon, bDebug=bDebug): continue
                else: return False
            elif sShapeA == 'Sphere':
                if Sphere.areEqual([shapeA, shapeB], epsilon=epsilon, bDebug=bDebug): continue
                else: return False
            elif sShapeA == "Torus":
                if Torus.areEqual([shapeA, shapeB], epsilon=epsilon, bDebug=bDebug): continue
                else: return False
    
        return True


class Cylinder():

    @staticmethod
    def areEqual(cylinders, epsilon=1e-9, bDebug=False):
        if len(cylinders) < 2: return
    
        for cyl in cylinders:
            if not isinstance(cyl, rg.Cylinder): return False
    
        cylA = cylinders[0]
    
        for i in range(1, len(cylinders)):
            cylB = cylinders[i]
        
            if cylA.EpsilonEquals(other=cylB, epsilon=epsilon):
                if bDebug: print "Cylinders are EpsilonEqual."
                continue
        
            # Simple cylinder comparison did not find a match, but this doesn't mean
            # they do not.
        
            # Compare Radii.
            if abs(cylA.Radius - cylB.Radius) > epsilon:
                return False
            if bDebug: print "Cylinders' Radii are EpsilonEqual."
        
            # Compare vectors of Axes.
            if (
                    not cylA.Axis.EpsilonEquals(other=cylB.Axis, epsilon=epsilon)
                    and
                    not cylA.Axis.EpsilonEquals(other=-cylB.Axis, epsilon=epsilon)
            ):
                return False
            if bDebug: print "Cylinders' Axes are EpsilonEqual."
        
            # Compare locations of Axes.
            lineA = rg.Line(cylA.Center, cylA.Axis)
            pt_ClosestOnA = lineA.ClosestPoint(testPoint=cylB.Center, limitToFiniteSegment=False)
            if pt_ClosestOnA.DistanceTo(cylB.Center) > epsilon:
                return False
    
        return True


    @staticmethod
    def extendToObjectSize(rgCyl, rgObjForSizeRef, bDebug=False):
        if not rgCyl.IsFinite: # Fix cylinder that is infinite (has 0 height).
            if bDebug: print "Making cylinder finite..."
            rgBbox_ForSizeRef = rgObjForSizeRef.GetBoundingBox(True)
            fAddLen = 1.1 * rgBbox_ForSizeRef.Min.DistanceTo(rgBbox_ForSizeRef.Max)
            rgCyl.Height1 = -fAddLen
            rgCyl.Height2 = fAddLen
        else: # Increase length of cylinder so that it extends the reference object.
            if bDebug: print "Height of cylinder before: {}".format(rgCyl.TotalHeight)
            rgBbox_ForSizeRef = rgObjForSizeRef.GetBoundingBox(False)
            fAddLen = 1.1 * rgBbox_ForSizeRef.Min.DistanceTo(rgBbox_ForSizeRef.Max)
            if rgCyl.Height1 < rgCyl.Height2:
                rgCyl.Height1 -= fAddLen
                rgCyl.Height2 += fAddLen
            else:
                rgCyl.Height1 += fAddLen
                rgCyl.Height2 -= fAddLen
            if bDebug: print "Height of cylinder after: {}".format(rgCyl.TotalHeight)


    @staticmethod
    def isSizeAcceptable(cylinder, bDebug=False):
    
        mmPerModelUnit = 1.0 / Rhino.RhinoMath.UnitScale(
                Rhino.UnitSystem.Millimeters,
                sc.doc.ModelUnitSystem)
    
        if (cylinder.Radius * mmPerModelUnit) >= MAX_ACCEPTED_SIZE:
            if bDebug:
                print "Cylinder with {} Radius skipped.".format(
                        cylinder.Radius)
            return False
        return True


class Cone():

    @staticmethod
    def areEqual(cones, epsilon=1e-9, bDebug=False):
        """
        # Changing the Cone's Height changes
        # its BasePoint location, but not its Radius.
        # Conclusion: Do not compare Height, Radius, and BasePoint.
        # Also not comparing Planes due to needed research on how to properly
        # do this, but there are probably enough other properties being tested.
        """
    
        if len(cones) < 2: return
    
        for cone in cones:
            if not isinstance(cone, rg.Cone): return False
    
        coneA = cones[0]
    
        for i in range(1, len(cones)):
            coneB = cones[i]
        
            if coneA.EpsilonEquals(other=coneB, epsilon=epsilon):
                if bDebug: print "Cones are EpsilonEqual."
                continue

            # Simple cone comparison did not find a match, but this doesn't mean
            # they do not, so will keep testing.
        
            if not coneA.Axis.EpsilonEquals(other=coneB.Axis, epsilon=epsilon):
                return False
            if bDebug: print "Cones' Axes are EpsilonEqual."
        
            # Keep testing.
            if not coneA.ApexPoint.EpsilonEquals(
                        other=coneB.ApexPoint,
                        epsilon=10.0*epsilon
            ):
                s  = "Not EpsilonEquals because Cone.ApexPoint location difference is"
                s += " {}, ".format(
                        coneA.ApexPoint.DistanceTo(other=coneB.ApexPoint))
                s += "but the epsilon is {}.".format(epsilon)
                print s
                return False
            if bDebug: print "Cones' ApexPoints are EpsilonEqual."
        
            # Keep testing.
            fTol_Angle = sc.doc.ModelAngleToleranceDegrees
            angle_A = coneA.AngleInDegrees()
            angle_B = coneB.AngleInDegrees()
            if not abs(angle_A - angle_B) <= fTol_Angle:
                s  = "Not EpsilonEquals because Cone.AngleInDegrees() difference is"
                s += " {}, ".format(
                        abs(coneA.AngleInDegrees()-coneB.AngleInDegrees()))
                s += "but the epsilon is {}.".format(fTol_Angle)
                print s
                return False
            if bDebug: print "Cones' Radii are EpsilonEqual."
        
    
        return True


    @staticmethod
    def isSizeAcceptableAtPoint(cone, point, bDebug=False):
    
        mmPerModelUnit = 1.0 / Rhino.RhinoMath.UnitScale(
                Rhino.UnitSystem.Millimeters,
                sc.doc.ModelUnitSystem)
    
        if abs(cone.ApexPoint.MaximumCoordinate * mmPerModelUnit) >= MAX_ACCEPTED_SIZE:
            if bDebug:
                print "Cone with ApexPoint of {} skipped.".format(
                        cone.ApexPoint)
            return False
        else:
            apexToPtOnCone = cone.ApexPoint.DistanceTo(point)
            if apexToPtOnCone <= sc.doc.ModelAbsoluteTolerance:
                # Point must be at apex.
                return
        
            line_Axis = rg.Line(cone.ApexPoint, cone.BasePoint)
            pt_onAxis = line_Axis.ClosestPoint(point, limitToFiniteSegment=False)
            apexToPtOnAxis = cone.ApexPoint.DistanceTo(pt_onAxis)
            if apexToPtOnAxis <= sc.doc.ModelAbsoluteTolerance:
                # Point must be at apex.
                return
        
            radiusAtPt = point.DistanceTo(pt_onAxis)
        
            if radiusAtPt >= MAX_ACCEPTED_SIZE:
                if bDebug: print "Cone with {} Radius skipped.".format(
                        cone_radius_at_srf)
                return False
        
            if radiusAtPt < sc.doc.ModelAbsoluteTolerance:
                if bDebug: print "Cone with {} Radius skipped.".format(
                        cone_radius_at_srf)
                return False
        
            # Older method.
            #            cone_radius_at_srf = (apexToPtOnCone**2.0 - apexToPtOnAxis**2.0)**0.5
            #            if (cone_radius_at_srf * mmPerModelUnit) >= MAX_ACCEPTED_SIZE:
            #                if bDebug:
            #                    print "Cone with {} Radius skipped.".format(
            #                            cone_radius_at_srf)
            #                return False
    
        return True


class Plane():

    @staticmethod
    def areEqual(planes, epsilon=1e-9, bDebug=False):
        if len(planes) < 2: return
    
        for plane in planes:
            if not isinstance(plane, rg.Plane): return False
    
        planeA = planes[0]
    
        ptA0 = planeA.PointAt(0.0, 0.0)
        ptAU = planeA.PointAt(1.0, 0.0)
        ptAV = planeA.PointAt(0.0, 1.0)
    
        #    for pt in (ptA0, ptAU, ptAV):
        #        sc.doc.Objects.AddPoint(pt)
    
        for i in range(1, len(planes)):
            planeB = planes[i]
        
            #        if bDebug:
            #            ptB_ToA0 = planeB.ClosestPoint(ptA0)
            #            ptB_ToAU = planeB.ClosestPoint(ptAU)
            #            ptB_ToAV = planeB.ClosestPoint(ptAV)
            #            for pt in (ptB_ToA0, ptB_ToAU, ptB_ToAV):
            #                sc.doc.Objects.AddPoint(pt)
        
            ptB_ToA0 = planeB.ClosestPoint(ptA0)
            fDelta_0 = ptA0.DistanceTo(ptB_ToA0)
            if fDelta_0 > epsilon:
                #            if bDebug:
                #                s  = "No match because planes' difference at A's "
                #                s += "0,0 is {:.2e}.".format(fDelta_0)
                #                print s
                return False
        
            ptB_ToAU = planeB.ClosestPoint(ptAU)
            fDelta_U = ptAU.DistanceTo(ptB_ToAU)
            if fDelta_U > epsilon:
                #            if bDebug:
                #                s  = "No match because planes' difference at A's "
                #                s += "1,0 is {:.2e}.".format(fDelta_U)
                #                print s
                return False
        
            ptB_ToAV = planeB.ClosestPoint(ptAV)
            fDelta_V = ptAV.DistanceTo(ptB_ToAV)
            if fDelta_V > epsilon:
                #            if bDebug:
                #                s  = "No match because planes' difference at A's "
                #                s += "0,1 is {:.2e}.".format(fDelta_V)
                #                print s
                return False
            
            # False positive may exist if planeA was derived from a small BrepFace.
            # Therefore, test the planes again in the opposite way.
            for uv in (0.0, 0.0), (1.0, 0.0), (0.0, 1.0):
                ptB = planeB.PointAt(*uv)
            
                ptA_ToB = planeA.ClosestPoint(ptB)
                fDelta = ptB.DistanceTo(ptA_ToB)
                if fDelta > epsilon:
                    return False
    
        return True


class Sphere():

    @staticmethod
    def areEqual(spheres, epsilon=1e-9, bDebug=False):
        """
        """
    
        if len(spheres) < 2: return
    
        for sphere in spheres:
            if not isinstance(sphere, rg.Sphere): return False
    
        sphereA = spheres[0]
    
        for i in range(1, len(spheres)):
            sphereB = spheres[i]
        
            if not sphereA.Center.EpsilonEquals(sphereB.Center, epsilon):
                return False
            if not Rhino.RhinoMath.EpsilonEquals(sphereA.Radius, sphereB.Radius, epsilon):
                return False
    
        return True


    @staticmethod
    def isSizeAcceptable(sphere, bDebug=False):
    
        mmPerModelUnit = 1.0 / Rhino.RhinoMath.UnitScale(
                Rhino.UnitSystem.Millimeters,
                sc.doc.ModelUnitSystem)
    
        if (sphere.Radius * mmPerModelUnit) >= MAX_ACCEPTED_SIZE:
            if bDebug:
                print "Sphere with {} Radius skipped.".format(
                        tol_Try, sphere.Radius)
            return False
        return True


class Torus():

    @staticmethod
    def areEqual(tori, epsilon=1e-9, bDebug=False):
        """
        """
    
        if len(tori) < 2: return
    
        for torus in tori:
            if not isinstance(torus, rg.Torus): return False
    
        torusA = tori[0]
    
        for i in range(1, len(tori)):
            torusB = tori[i]
        
            if torusA.EpsilonEquals(other=torusB, epsilon=epsilon):
                if bDebug: print "Tori are EpsilonEqual."
                continue
            else:
                fDelta_MajRad = abs(torusA.MajorRadius - torusB.MajorRadius)
                if fDelta_MajRad <= epsilon:
                    pass
                else:
                    if bDebug:
                        s  = "No match because tori's"
                        s += " MajorRadius difference is {:.2e}.".format(fDelta_MajRad)
                        print s
                    return False
            
                fDelta_MinRad = abs(torusA.MinorRadius - torusB.MinorRadius)
                if fDelta_MinRad <= epsilon:
                    pass
                else:
                    if bDebug:
                        s  = "No match because tori's"
                        s += " MinorRadius difference is {:.2e}.".format(fDelta_MinRad)
                        print s
                    return False
            
                # Using custom function because Plane.EpsilonEquals is too strict for this application.
                if Plane.areEqual(
                        [torusA.Plane, torusB.Plane],
                        epsilon=epsilon,
                        bDebug=bDebug):
                    if bDebug: print "Tori's Planes are EpsilonEqual."
                else:
                    return False
    
        return True


    @staticmethod
    def isSizeAcceptable(torus, bDebug=False):
    
        mmPerModelUnit = 1.0 / Rhino.RhinoMath.UnitScale(
                Rhino.UnitSystem.Millimeters,
                sc.doc.ModelUnitSystem)
    
        if (torus.MajorRadius * mmPerModelUnit) >= MAX_ACCEPTED_SIZE:
            if bDebug:
                print "Torus with {} MajorRadius skipped.".format(torus.MajorRadius)
            return False
        elif (torus.MinorRadius * mmPerModelUnit) >= MAX_ACCEPTED_SIZE:
            if bDebug:
                print "Torus with {} MinorRadius skipped.".format(torus.MinorRadius)
            return False
        return True


class Surface():

    @staticmethod
    def tryGetPlane(rgSrf0, fTolerance=1e-9, bDebug=False):
        """
        Updated in May of 2019.
    
        Returns a tuple:
            On success: (rg.Plane instance), (float: tolerance actually used to obtain shape)
            On fail: None, None
        """

        exp = -52
        m_eps = 2**exp # 2.22044604925e-16

        if fTolerance <= 0.0:
            return None, None
        elif fTolerance < m_eps:
            fTolerance = m_eps

        if isinstance(rgSrf0, rg.BrepFace):
            rgSrf0 = rgSrf0.UnderlyingSurface()
    
        if isinstance(rgSrf0, rg.NurbsSurface):
            rgNurbsSrf1 = rgSrf0.Duplicate()
        else:
            rgNurbsSrf1 = rgSrf0.ToNurbsSurface()
    
        # Start tolerance to use at a low value and iterate up to input tolerance.
        fTol_Attempting = m_eps

        while fTol_Attempting <= fTolerance:
            sc.escape_test()

            b, plane = rgNurbsSrf1.TryGetPlane(fTol_Attempting)
            if b:
                rgNurbsSrf1.Dispose()
                return plane, fTol_Attempting

            if fTol_Attempting == fTolerance:
                break

            exp += 1
            fTol_Attempting = 2**exp

            if fTol_Attempting > fTolerance:
                fTol_Attempting = fTolerance
    
        return None, None


    @staticmethod
    def tryGetRoundPrimitive(rgSrf0, bCylinder=True, bCone=True, bSphere=True, bTorus=True, fTolerance=1e-9, bDebug=False):
        """
        Updated in May of 2019.
    
        Returns a tuple:
            On success: (rg.Cone, Cylinder, Sphere, or Torus instance), (float: tolerance actually used to obtain shape)
            On fail: None, None
    
        Tolerance up to fTolerance parameter is first tested for a Cylinder
        Then the other 3 shapes are tested together each tolerance trial at a time.
        """

        if not any((bCylinder, bCone, bSphere, bTorus)):
            return None, None

        if fTolerance <= 0.0:
            return None, None

        exp = -52
        m_eps = 2**exp # 2.22044604925e-16


        if isinstance(rgSrf0, rg.BrepFace):
            rgSrf0 = rgSrf0.UnderlyingSurface()
    
        if isinstance(rgSrf0, rg.NurbsSurface):
            rgNurbsSrf1 = rgSrf0.Duplicate()
        else:
            rgNurbsSrf1 = rgSrf0.ToNurbsSurface()

        # If initially True, these may later be set to False when either
        #   1. The shape cannot be found at fTolerance.
        #   2. The relative shape is found but not at an acceptable size
        bCone_WIP = bCone
        bSphere_WIP = bSphere
        bTorus_WIP = bTorus
    
        if fTolerance < m_eps:
            fTolerance = m_eps

        # Start tolerance to use at a low value and iterate up to input tolerance.
        fTol_Attempting = m_eps

        # Cylinder has preference, so all tolerances will be iterated, trying for that shape.
        if bCylinder:
            # Check max. tolerance first.
            b, cylinder = rgNurbsSrf1.TryGetCylinder(fTolerance)
            if b:
                cylinder = None
                while fTol_Attempting <= fTolerance:
                    sc.escape_test()
                    b, cylinder = rgNurbsSrf1.TryGetCylinder(fTol_Attempting)
                    if b:
                        if Cylinder.isSizeAcceptable(cylinder, bDebug=bDebug):
                            rgNurbsSrf1.Dispose()
                            return cylinder, fTol_Attempting
                        else:
                            break
            
                    if fTol_Attempting == fTolerance:
                        break
    
                    exp += 1
                    fTol_Attempting = 2**exp
    
                    #fTol_Attempting *= 10.0
    
                    if fTol_Attempting > fTolerance:
                        fTol_Attempting = fTolerance

        if bCone_WIP:
            b, cone = rgNurbsSrf1.TryGetCone(fTolerance)
            if b:
                cone = None
            else:
                bCone_WIP = False
        if bSphere_WIP:
            b, sphere = rgNurbsSrf1.TryGetSphere(fTolerance)
            if b:
                sphere = None
            else:
                bSphere_WIP = False
        if bTorus_WIP:
            b, torus = rgNurbsSrf1.TryGetTorus(fTolerance)
            if b:
                torus = None
            else:
                bTorus_WIP = False

        if not any((bCone_WIP, bSphere_WIP, bTorus_WIP)):
            return None, None

        # Reset starting tolerance for non-cylinder search.
        exp = -52
        fTol_Attempting = m_eps

        while fTol_Attempting <= fTolerance:
            sc.escape_test()
        
            if bCone_WIP:
                b, cone = rgNurbsSrf1.TryGetCone(fTol_Attempting)
                if b:
                    for pt in (
                            rgNurbsSrf1.PointAt(
                                    rgNurbsSrf1.Domain(0).T0,
                                    rgNurbsSrf1.Domain(1).T0),
                            rgNurbsSrf1.PointAt(
                                    rgNurbsSrf1.Domain(0).T1,
                                    rgNurbsSrf1.Domain(1).T1)
                    ):
                        rc = Cone.isSizeAcceptableAtPoint(
                                cone,
                                point=pt,
                                bDebug=bDebug)
                        if rc is None:
                            continue
                        elif rc is False:
                            bCone_WIP = False
                            break
                    else:
                        # Size is acceptable.
                        rgNurbsSrf1.Dispose()
                        return cone, fTol_Attempting
            if bSphere_WIP:
                b, sphere = rgNurbsSrf1.TryGetSphere(fTol_Attempting)
                if b:
                    if Sphere.isSizeAcceptable(sphere, bDebug=bDebug):
                        rgNurbsSrf1.Dispose()
                        return sphere, fTol_Attempting
                    else:
                        bSphere_WIP = False
            if bTorus_WIP:
                b, torus = rgNurbsSrf1.TryGetTorus(fTol_Attempting)
                if b:
                    if Torus.isSizeAcceptable(torus, bDebug=bDebug):
                        rgNurbsSrf1.Dispose()
                        return torus, fTol_Attempting
                    else:
                        bTorus_WIP = False
        
            if fTol_Attempting == fTolerance:
                rgNurbsSrf1.Dispose()
                return None, None

            exp += 1
            fTol_Attempting = 2**exp

            #fTol_Attempting *= 10.0

            if fTol_Attempting > fTolerance:
                fTol_Attempting = fTolerance


    @staticmethod
    def tryGetPrimitiveShape(rgSrf0, bPlane=True, bCylinder=True, bCone=True, bSphere=True, bTorus=True, fTolerance=1e-9, bDebug=False):
        """
        190722: Routines were split into 2 other functions.
    
        Returns a tuple:
            On success: (rg.Cone, Cylinder, Plane, Sphere, or Torus instance), (float: tolerance actually used to obtain shape)
            On fail: None, None
        """
    
        if isinstance(rgSrf0, rg.BrepFace):
            rgSrf0 = rgSrf0.UnderlyingSurface()
    
        if isinstance(rgSrf0, rg.NurbsSurface):
            rgNurbsSrf1 = rgSrf0
        else:
            rgNurbsSrf1 = rgSrf0.ToNurbsSurface()
    
        rc = None, None
    
        if bPlane:
            rc = Surface.tryGetPlane(
                    rgSrf0=rgSrf0,
                    fTolerance=fTolerance,
                    bDebug=bDebug)
            if rc[0] is not None:
                return rc
    
        if bCylinder or bCone or bSphere or bTorus:
            rc = Surface.tryGetRoundPrimitive(
                    rgSrf0=rgSrf0,
                    bCylinder=bCylinder,
                    bCone=bCone,
                    bSphere=bSphere,
                    bTorus=bTorus,
                    fTolerance=fTolerance,
                    bDebug=bDebug)
    
        return rc


class BrepFace():


    @staticmethod
    def tryGetPlane(rgFace0, bMatchToShrunkFace=True, fTolerance=1e-9, bDebug=False):
        """
        Return on success: tuple(
            tuple(
                rg.Plane
                float(Tolerance used)
                bool(Whether shrunk surface was used),
            None)
        Return on fail: tuple(
            None
            str(LogStatement))
        """
        if bMatchToShrunkFace:
            rgBrep_1Face_Shrunk = rgFace0.DuplicateFace(False)
            rgBrep_1Face_Shrunk.Faces.ShrinkFaces()
            rgFace_Shrunk = rgBrep_1Face_Shrunk.Faces[0]
            rgSrfs_ToTry = rgFace0.UnderlyingSurface(), rgFace_Shrunk.UnderlyingSurface()
        else:
            rgSrfs_ToTry = rgFace0.UnderlyingSurface(),
    
        for i_NotShrunk_Shrunk, rgSrf_ToTry in enumerate(rgSrfs_ToTry):
            rc = Surface.tryGetPlane(
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


    @staticmethod
    def tryGetRoundPrimitive(rgFace0, bMatchToShrunkFace=True, bCylinder=True, bCone=True, bSphere=True, bTorus=True, fTolerance=1e-9, bDebug=False):
        """
        Return on success: tuple(
            tuple(
                rg.Cone, Cylinder, Sphere, or Torus
                float(Tolerance used)
                bool(Whether shrunk surface was used),
            None)
        Return on fail: tuple(
            None
            str(LogStatement))
        """
    
        if bMatchToShrunkFace:
            rgBrep_1Face_Shrunk = rgFace0.DuplicateFace(False)
            rgBrep_1Face_Shrunk.Faces.ShrinkFaces()
            rgFace_Shrunk = rgBrep_1Face_Shrunk.Faces[0]
            rgFaces_ToTest = rgFace0, rgFace_Shrunk
            rgSrfs_ToTry = rgFace0.UnderlyingSurface(), rgFace_Shrunk.UnderlyingSurface()
        else:
            rgSrfs_ToTry = rgFace0.UnderlyingSurface(),
    
        for i_NotShrunk_Shrunk, rgSrf_ToTry in enumerate(rgSrfs_ToTry):
            rc = Surface.tryGetRoundPrimitive(
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


    @staticmethod
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
            rc = BrepFace.tryGetPlane(
                    rgFace0,
                    bMatchToShrunkFace=bMatchToShrunkFace,
                    fTolerance=fTolerance,
                    bDebug=bDebug)
            if rc[0] is not None:
                return rc

        if bCylinder or bCone or bSphere or bTorus:
            rc = BrepFace.tryGetRoundPrimitive(
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

