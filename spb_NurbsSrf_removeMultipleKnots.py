"""

"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221224: WIP: Created, starting with another script.

TODO: Add support for periodic surfaces.
"""


import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid
from System.Collections.Generic import List
from System.Diagnostics import Stopwatch

import xBrep_getDistancesBetween2
import xBrepObject


class Opts:

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'bLimitDev'; keys.append(key)
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fDevTol'; keys.append(key)
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bReplace'; keys.append(key)
    values[key] = True
    names[key] = 'DocAction'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Add', 'Replace')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bExtract'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'bSrfOnFail'; keys.append(key)
    values[key] = True
    names[key] = 'AddSrfOnFail'
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

        if key == 'fDevTol':
            if cls.riOpts[key].CurrentValue < 0.0:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue
            elif cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = 1e-6

            cls.values[key] = cls.riOpts[key].CurrentValue
            sc.sticky[cls.stickyKeys[key]] = cls.values[key]
            return

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return

        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get BrepFaces with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select faces")

    go.GeometryFilter = rd.ObjectType.Surface

    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)

    bPreselectedObjsChecked = False

    idxs_Opt = {}

    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

        go.AcceptNumber(Opts.values['bLimitDev'], acceptZero=Opts.values['bLimitDev'])

        addOption('bLimitDev')
        if Opts.values['bLimitDev']:
            addOption('fDevTol')
        addOption('bReplace')
        if Opts.values['bReplace']:
            addOption('bExtract')
        addOption('bSrfOnFail')
        addOption('bEcho')
        addOption('bDebug')


        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, ignoreUnacceptablePreselectedObjects=True)
            continue

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            if Opts.values['bLimitDev']:
                key = 'fDevTol'
                Opts.riOpts[key].CurrentValue = go.Number()
            else:
                continue

            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def multiplicityRangeOfInteriorKnots(ns, iDir):
    degree = ns.Degree(iDir)
    knots = ns.KnotsV if iDir else ns.KnotsU
    ms = []
    iK = degree
    while iK < knots.Count-degree:
        sc.escape_test()
        m = knots.KnotMultiplicity(iK)
        ms.append(m)
        iK += m
    return ms


def someInteriorKnotsHaveMultiplicityAbove1(ns):
    for iDir in (0,1):
        degree = ns.Degree(iDir)
        knots = ns.KnotsV if iDir else ns.KnotsU

        for iK in range(degree, knots.Count-degree):
            if knots.KnotMultiplicity(iK) > 1:
                return True

    return False


def processSrf(ns_In, fDevTol, **kwargs):

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script only works in Rhino V6 and above.")
        return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    if not isinstance(ns_In, rg.NurbsSurface):
        if bDebug: print("{} skipped.".format(ns_In.GetType().Name))
        return

    if ns_In.SpanCount(0) == 1 and ns_In.SpanCount(1) == 1:
        if bDebug: print("No interior knots.")
        return

    degrees = ns_In.Degree(0), ns_In.Degree(1)
    knots_In = ns_In.KnotsU, ns_In.KnotsV

    ms = (
        multiplicityRangeOfInteriorKnots(ns_In, 0),
        multiplicityRangeOfInteriorKnots(ns_In, 1),
        )

    if max(ms[0]) == 1 and max(ms[1]) == 1:
        if bDebug: print("No interior knots with multiplicity above 1.")
        return


    if bDebug: print("In removeKnots function:")


    # Duplicate NurbsSurface since RemoveKnots modifies the input surface.
    ns_Out = ns_In.Duplicate()
    
    knots_Out = ns_Out.KnotsU, ns_Out.KnotsV


    rgMeshParams = rg.MeshingParameters.QualityRenderMesh


    # NurbsCurveKnotList.RemoveMultipleKnots' minimumMultiplicity and
    # maximumMultiplicity are, respectively, the
    # minimum and maximum allowed (not to remove) knot multiplicities.

    for iDir in (0,1):
        knots_Out[iDir].RemoveMultipleKnots(
            minimumMultiplicity=0,
            maximumMultiplicity=degrees[iDir]+1,
            tolerance=Rhino.RhinoMath.UnsetValue)


    rc = xBrep_getDistancesBetween2.getDistancesBetweenBreps(
        ns_In,
        ns_Out,
        rgMeshParams,
        bCalcBrepIntersection=False)

    if rc[0]:
        rgSrf1_LastGood = ns_Out.Duplicate()
        fDev = rc[1]

    if fDevTol is None or (fDev <= fDevTol):

        iCt_KnotsRemoved = (
            knots_In[0].Count - knots_Out[0].Count +
            knots_In[1].Count - knots_Out[1].Count)

        return ns_Out, iCt_KnotsRemoved, fDev

    ns_Out.Dispose()


    print(ms[0])
    print(ms[1])

    for iDir in (0,1):
        if len(ms[(iDir+1)%2]) == 0:
            for minimumMultiplicity in range(max(ms[iDir])-1, 0, -1):
                print(minimumMultiplicity)

    return

    for minimumMultiplicity_U in range(1, degrees[0]-1):
        for minimumMultiplicity_V in range(1, degrees[1]-1):
            pass


    ns_Out = ns_In.Duplicate()






    # If surface is within deviation tolerance,
    # make a backup of the new surface for the next iteration
    # and continue.
    # Otherwise, revert rgSrf1 to rgSrf1_LastGood.

    # Duplicate rgSrf1 in case removing a knot creates an unacceptable deviation, etc.
    rgSrf1_LastGood = ns_Out.Duplicate()
    fDev_LastGood = None
    
    for iDir in 0,1:
        sc.escape_test()
        
        idxKnot1 = degrees[iDir]
        
        # Iterate through each interior knot of original NurbsSurface.
        for idxKnot0 in xrange(degrees[iDir], knots_In[iDir].Count-degrees[iDir]-1):
            sc.escape_test()
            
            if bDebug:
                s = ""
                sEval = 'iDir'; s += "{}:{}".format(sEval, eval(sEval))
                sEval = 'knots0[iDir].Count'; s += "  {}:{}".format(sEval, eval(sEval))
                sEval = 'idxKnot0'; s += "  {}:{}".format(sEval, eval(sEval))
                sEval = 'idxKnot1'; s += "  {}:{}".format(sEval, eval(sEval))
                sEval = 'knots1[iDir].Count'; s += "  {}:{}".format(sEval, eval(sEval))
                print(s)
            
            knot0Multiplicity = knots_In[iDir].KnotMultiplicity(idxKnot0)

            if (
                knot0Multiplicity >= nKnotMulti_Min and
                knot0Multiplicity <= nKnotMulti_Max
            ):
                if knots_Out[iDir].RemoveKnots(index0 = idxKnot1,
                                      index1 = idxKnot1 + 1):
                    rc = xBrep_getDistancesBetween2.getDistancesBetweenBreps(ns_In, ns_Out,
                            rgMeshParams, bCalcBrepIntersection=False)

                    # If surface is within deviation tolerance,
                    # make a backup of the new surface for the next iteration
                    # and continue.
                    # Otherwise, revert rgSrf1 to rgSrf1_LastGood.
                    if rc[0] and rc[1] <= fDevTol:
                        rgSrf1_LastGood = ns_Out.Duplicate()
                        fDev_LastGood = rc[1]
                        # Notice that since the knot was removed,
                        # the next knot index will be the same as this iteration.
        #                    if bDebug:
        #                        sEval = 'idxKnot0'; print(sEval + ':', eval(sEval)
                        continue
                    else:
                        ns_Out = rgSrf1_LastGood.Duplicate()
                        knots_Out = ns_Out.KnotsU, ns_Out.KnotsV
                        idxKnot1 += 1 # Advance knot index to try next knot.
                
                else:
                    print("knots1[iDir].RemoveKnots fail!")
                    return
    
    if bDebug:
        print("Final")
        for iDir in 0,1:
            sEval = 'iDir'; s = "{}:{}".format(sEval, eval(sEval))
            sEval = 'knots0[iDir].Count'; s += "  {}:{}".format(sEval, eval(sEval))
            sEval = 'idxKnot0'; s += "  {}:{}".format(sEval, eval(sEval))
            sEval = 'idxKnot1'; s += "  {}:{}".format(sEval, eval(sEval))
            sEval = 'knots1[iDir].Count'; s += "  {}:{}".format(sEval, eval(sEval))
            print(s)
    
    iCt_KnotsRemoved = (
        knots_In[0].Count - knots_Out[0].Count +
        knots_In[1].Count - knots_Out[1].Count)



    if iCt_KnotsRemoved > 0:
        return ns_Out, iCt_KnotsRemoved, fDev_LastGood


def processFace(rgFace_In, fDevTol, **kwargs):
    """
    Parameters:
        rgFace_In
        fDevTol
        bDebug
    Returns:
        (rg.Brep (1-face), float(deviation)), None
        (None, float(deviation needed)), sLog
        None, None
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script only works in Rhino V6 and above.")
        return



    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rgSrf_In = rgFace_In.UnderlyingSurface()


    rgNurbsSrf0 = rgSrf_In


    rc = processSrf(
        rgNurbsSrf0,
        fDevTol,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rc is None or rc[0] is None:
        return rc

    rgNurbsSrf_Converted, iCt_KnotsRemoved, srf_dev = rc

    # Success in creating NurbsSurface.  Now, create correctly trimmed brep.
    
    rgBrep0_1Face = rgFace_In.DuplicateFace(duplicateMeshes=False)
    
    rgBrep0_1Face.Faces[0].RebuildEdges(0.1*sc.doc.ModelAbsoluteTolerance, True, True)
    #rgBrep0_1Face.Faces.ShrinkFaces()

    if rgBrep0_1Face.IsSurface:
        rgBrep0_1Face.Dispose()
        rgBrep1_1Face = rgNurbsSrf_Converted.ToBrep()
        rgNurbsSrf_Converted.Dispose()
        if not rgBrep1_1Face.IsValid:
            rgBrep1_1Face.Dispose()
            return None, "Invalid brep geometry after ToBrep."
        return rgBrep1_1Face, srf_dev, None
    
    # Test areas before trimming.
    fArea_Trimmed = rgBrep0_1Face.GetArea()
    if fArea_Trimmed:
        rgBrep_Untrimmed = rgBrep0_1Face.Faces[0].UnderlyingSurface().ToBrep()
        fArea_Untrimmed = rgBrep_Untrimmed.GetArea()
        rgBrep_Untrimmed.Dispose()
        if fArea_Untrimmed:
            if abs(fArea_Trimmed - fArea_Untrimmed) <= sc.doc.ModelAbsoluteTolerance:
                rgBrep1_1Face = rgNurbsSrf_Converted.ToBrep()
                rgNurbsSrf_Converted.Dispose()
                if not rgBrep1_1Face.IsValid:
                    rgBrep1_1Face.Dispose()
                    return None, "Invalid brep geometry after ToBrep."
                return rgBrep1_1Face, srf_dev, None

    rgBrep1_1Face = xBrepFace.retrimFace(
            rgBrep0_1Face.Faces[0],
            rgSrf_Replacement=rgNurbsSrf_Converted,
            fSplitTol=1.0*sc.doc.ModelAbsoluteTolerance if fDevTol is None else fDevTol,
            bDebug=bDebug
    )
    rgNurbsSrf_Converted.Dispose()
    rgBrep0_1Face.Dispose()

    if rgBrep1_1Face is None:
        return None, "xBrepFace.createMonofaceBrep returned None."

    if not rgBrep1_1Face.IsValid:
        rgBrep1_1Face.Dispose()
        return None, "An invalid brep was skipped."

    return rgBrep1_1Face, srf_dev, None


def processBrep(rgBrep_In, idxs_rgFaces, fDevTol, **kwargs):
    """
    Returns:
        (rgBreps_1F_Mod, idxs_rgFaces_Rebuilt, srf_devs), sLogs
        None
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script only works in Rhino V6 and above.")
        return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    rgBreps_1F_Mod = []
    idxs_rgFaces_Rebuilt = []
    srf_devs = []
    srf_devs_Needed = []
    
    rgB_WIP = rgBrep_In.DuplicateBrep()
    
    sLogs = []
    
    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt
    
    idxs_AtTenths = [int(round(0.1*i*len(idxs_rgFaces),0)) for i in range(10)]
    
    for iF, idx_rgFace in enumerate(idxs_rgFaces):
        if sc.escape_test(False):
            print("Searching interrupted by user.")
            return
        
        if iF in idxs_AtTenths:
            s = sCmdPrompt0 + ", {:d}% of {} faces in current brep ...".format(
                int(100.0 * (iF+1) / len(idxs_rgFaces)), len(idxs_rgFaces))
            
            if bDebug:
                print(s)
            else:
                Rhino.RhinoApp.SetCommandPrompt(s)
        
        rgFace_In = rgB_WIP.Faces[idx_rgFace]


        rc = processFace(
            rgFace_In,
            fDevTol=fDevTol,
            bDebug=bDebug
            )

        if rc is None: return

        rgBrep_1F_Converted, srf_dev, sLog = rc

        if rgBrep_1F_Converted is None:
            if srf_dev is not None:
                srf_devs_Needed.append(srf_dev)
            if sLog is not None:
                sLogs.append(sLog)
            continue
        
        rgBreps_1F_Mod.append(rgBrep_1F_Converted)
        idxs_rgFaces_Rebuilt.append(idx_rgFace)
        srf_devs.append(srf_dev)
    
    rgB_WIP.Dispose()
    
    if srf_devs_Needed:
        s  = "Need tolerances of {:.2e} through {:.2e}".format(
                min(srf_devs_Needed), max(srf_devs_Needed))
        s += " to convert remaining convertible surfaces."
        sLogs.append(s)
    
    return (rgBreps_1F_Mod, idxs_rgFaces_Rebuilt, srf_devs), sLogs


def processBrepObjects(rhBreps, idx_Faces=None, fDevTol=None, **kwargs):
    """
    Parameters:
        rhBreps: Objrefs of brep with face components, GUIDs, rd.Breps.
    """
    
    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script only works in Rhino V6 and above.")
        return


    def getOpt(key): return kwargs[key] if key in kwargs else Opts.values[key]

    bReplace = getOpt('bReplace')
    bExtract = getOpt('bExtract')
    bSrfOnFail = getOpt('bSrfOnFail')
    bEcho = getOpt('bEcho')
    bDebug = getOpt('bDebug')


    def getRhinoObject(rhObj):
        """
        'Deleted objects cannot be found by id.'
        (https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_DocObjects_Tables_ObjectTable_FindId.htm)
        """
        rdObj = None
        if isinstance(rhObj, rd.RhinoObject):
            rdObj = rhObj
        elif isinstance(rhObj, rd.ObjRef):
            rdObj = rhObj.Object()
        elif isinstance(rhObj, Guid):
            rdObj = sc.doc.Objects.FindId(rhObj) if Rhino.RhinoApp.ExeVersion >= 6 else sc.doc.Objects.Find(rhObj)
        return rdObj


    def getBrepObject(rhObj):
        rdObj = getRhinoObject(rhObj)
        if rdObj and (rdObj.ObjectType == rd.ObjectType.Brep):
            return rdObj


    def getSortedBrepIdsAndFaces(objrefs):
        """
        Parameters:
            list(objrefs)
        Returns:
            list(Brep GUIDs)
            list(lists(integers of Face indices) per brep)
        """
        
        gBreps0 = []
        idxs_Faces_perBrep = []
    
        for o in objrefs:
            gBrep0 = o.ObjectId
            rdBrep_In = o.Object()
            rgBrep_In = o.Brep()
        
            if not rgBrep_In.IsValid:
                print("Brep {} is invalid.  Fix first.".format(gBrep0))
                rgBrep_In.Dispose()
                continue
        
            idx_CompIdx = o.GeometryComponentIndex.Index
            if idx_CompIdx == -1:
                if gBrep0 in gBreps0:
                    idxs_Faces_perBrep[gBreps0.index(gBrep0)] = range(rgBrep_In.Faces.Count)
                else:
                    gBreps0.append(gBrep0)
                    idxs_Faces_perBrep.append(range(rgBrep_In.Faces.Count))
            else:
                rgFace_Brep0 = o.Face()
                if gBrep0 in gBreps0:
                    if rgFace_Brep0 in idxs_Faces_perBrep[gBreps0.index(gBrep0)]:
                        continue
                    else:
                        idxs_Faces_perBrep[gBreps0.index(gBrep0)].append(rgFace_Brep0.FaceIndex)
                else:
                    gBreps0.append(gBrep0)
                    idxs_Faces_perBrep.append([rgFace_Brep0.FaceIndex])

        return gBreps0, idxs_Faces_perBrep


    def isKnotVectorUniform(knots):
        return (
            (knots.KnotStyle == rg.KnotStyle.Uniform) or
            (knots.KnotStyle == rg.KnotStyle.QuasiUniform) or
            (
                (knots.KnotStyle == rg.KnotStyle.PiecewiseBezier) and
                knots.Count == knots.KnotMultiplicity(0) * 2)
            )


    def getNurbsSurfaceChangeDescription(rgNurbsSrf1, rgNurbsSrf2):
        s  = "  Prop:I->O"
        s += "  {}:({}->{})x({}->{})".format(
                "Deg",
                rgNurbsSrf1.OrderU-1,
                rgNurbsSrf2.OrderU-1,
                rgNurbsSrf1.OrderV-1,
                rgNurbsSrf2.OrderV-1,
        )
        s += "  {}:({}->{})x({}->{})".format(
                "PtCt",
                rgNurbsSrf1.Points.CountU,
                rgNurbsSrf2.Points.CountU,
                rgNurbsSrf1.Points.CountV,
                rgNurbsSrf2.Points.CountV,
        )
        if Rhino.RhinoApp.ExeVersion >= 7:
            s += "  {}:({}->{})x({}->{})".format(
                    "IsUniform",
                    str(isKnotVectorUniform(rgNurbsSrf1.KnotsU))[0],
                    str(isKnotVectorUniform(rgNurbsSrf2.KnotsU))[0],
                    str(isKnotVectorUniform(rgNurbsSrf1.KnotsV))[0],
                    str(isKnotVectorUniform(rgNurbsSrf2.KnotsV))[0],
            )
        s += "  {}:{}->{}".format(
                "IsRational",
                str(rgNurbsSrf1.IsRational)[0],
                str(rgNurbsSrf2.IsRational)[0],
        )
        s += "  {}:({}->{})x({}->{})".format(
                "IsClosed",
                str(rgNurbsSrf1.IsClosed(0))[0],
                str(rgNurbsSrf2.IsClosed(0))[0],
                str(rgNurbsSrf1.IsClosed(1))[0],
                str(rgNurbsSrf2.IsClosed(1))[0],
        )
        if (    rgNurbsSrf1.IsClosed(0) or rgNurbsSrf1.IsClosed(1) or
                rgNurbsSrf2.IsClosed(0) or rgNurbsSrf2.IsClosed(1)
        ):
            s += "  {}:{}->{}".format(
                    "IsPeriodic",
                    str(rgNurbsSrf1.IsPeriodic)[0],
                    str(rgNurbsSrf2.IsPeriodic)[0],
            )
        return s


    gBreps0, idxs_rgFace_PerBrep = getSortedBrepIdsAndFaces(rhBreps)
    if not gBreps0: return

    gBs1_perB0 = []
    srf_devs_All = []
    sLogs_All = []
    
    len_gBreps0 = len(gBreps0)
    idxs_AtTenths = [int(round(0.1*i*len_gBreps0,0)) for i in range(10)]
    
    sCmdPrompt0 = Rhino.RhinoApp.CommandPrompt
    
    if len(rhBreps) == 1:
        s = sCmdPrompt0 + "Brep"
        Rhino.RhinoApp.SetCommandPrompt(s)
    
    for iB, (gBrep0, idxFaces) in enumerate(zip(gBreps0, idxs_rgFace_PerBrep)):
        rdBrep_In = getBrepObject(gBrep0)
        rgBrep_In = rdBrep_In.Geometry

        if len(rhBreps) > 1:
            if iB in idxs_AtTenths:
                s = sCmdPrompt0 + "  At {:d}% of {} breps".format(
                    int(100.0 * (iB+1) / len_gBreps0), len_gBreps0)
        
            if bDebug:
                print(s)
            else:
                Rhino.RhinoApp.SetCommandPrompt(s)


        rc = processBrep(
            rgBrep_In,
            idxFaces,
            fDevTol=fDevTol,
            bDebug=bDebug)
        if rc is None: return

        sLogs_All.extend(rc[1])
        if not rc[0]: continue

        rgBreps_1F_Mod_thisBrep, idxsFs_Mod, srf_devs = rc[0]
        
        if not rgBreps_1F_Mod_thisBrep:
            rgBrep_In.Dispose()
            continue
        
        srf_devs_All.extend(srf_devs)
        
        if not bReplace:
            gBrep1_thisBrep = []
            for rgBrep_1F_New, idxF in zip(rgBreps_1F_Mod_thisBrep, idxsFs_Mod):
                gBrep1 = sc.doc.Objects.AddBrep(rgBrep_1F_New)
                if gBrep1 != Guid.Empty:
                    gBrep1_thisBrep.append(gBrep1)
            gBs1_perB0.append(gBrep1_thisBrep)
        else:
            # bReplace==True.
            if bExtract:
                rc = xBrepObject.replaceFaces(
                        rdBrep_In,
                        idxsFs_Mod,
                        rgBreps_1F_Mod_thisBrep,
                        bExtract=True,
                        fTolerance_Join=max(srf_devs))
                if rc:
                    gBreps1_NewFaces_thisBrep, gBreps_RemainingBrep = rc
                    gBs1_perB0.append(gBreps1_NewFaces_thisBrep)
            else:
                fTols_Edges = [edge.Tolerance for edge in rgBrep_In.Edges]
                fTolerance_Join = max((
                        2.0*fDevTol if fDevTol is not None else 0.0,
                        1.1*max(fTols_Edges),
                        sc.doc.ModelAbsoluteTolerance))
                rc = xBrepObject.replaceFaces(
                        rdBrep_In,
                        idxsFs_Mod,
                        rgBreps_1F_Mod_thisBrep,
                        bExtract=False,
                        fTolerance_Join=fTolerance_Join)
                if rc:
                    gBreps_withReplacedFaces_thisBrep = rc
                    gBs1_perB0.append(gBreps_withReplacedFaces_thisBrep)

        if bDebug or bEcho and len(gBreps0) == 1 and len(idxsFs_Mod)==1:
            rgBrep_1F_Mod = rgBreps_1F_Mod_thisBrep[0]
            s  = getNurbsSurfaceChangeDescription(
                    rgBrep_In.Faces[idxsFs_Mod[0]].UnderlyingSurface(),
                    rgBrep_1F_Mod.Surfaces[0])
            s += "  Deviation: {0:.2e}".format(srf_devs[0])
            print(s)

        for brep in rgBreps_1F_Mod_thisBrep: brep.Dispose()
        rgBrep_In.Dispose()

        if (
                any(g for gs in gBs1_perB0 for g in gs) and
                len(gBreps0)==1
                and (bDebug or bEcho)
        ):
            sc.doc.Objects.UnselectAll()
            if bReplace:
                if bExtract:
                    print("{} face(s) extracted from brep and replaced.".format(len(gBreps1_NewFaces_thisBrep)))
                else:
                    print("{} face(s) replaced in brep.".format(len(idxsFs_Mod)))
            else:
                print("{} monoface brep(s) added.".format(
                    sum(len(bs) for bs in rgBreps_1F_Mod_thisBrep)))
    
    return (gBs1_perB0, srf_devs_All), sLogs_All


def formatDistance(fDistance):
    """Returns: str"""
    if fDistance is None:
        return "(None)"
    elif fDistance == 0.0:
        return "exactly 0".format(fDistance)
    elif fDistance < 10.0**(-(sc.doc.DistanceDisplayPrecision-2)):
        return "{:.1e}".format(fDistance)
    else:
        return "{:.{}f}".format(fDistance, sc.doc.ModelDistanceDisplayPrecision)


def main():

    if Rhino.RhinoApp.ExeVersion < 6:
        print("This script works only in Rhino V6 and above.")
        return

    objrefs_In = getInput()
    if objrefs_In is None: return

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    fDevTol = Opts.values['fDevTol'] if Opts.values['bLimitDev'] else None
    bReplace = Opts.values['bReplace']
    bExtract = Opts.values['bExtract']
    bSrfOnFail = Opts.values['bSrfOnFail']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']


    if not bDebug: sc.doc.Views.RedrawEnabled = False

    sc.doc.Objects.UnselectAll()


    bDocModified = False


    rc = processBrepObjects(
        rhBreps=objrefs_In,
        idxs_Faces=None,
        fDevTol=fDevTol,
        bSrfOnFail=bSrfOnFail,
        bReplace=bReplace,
        bExtract=bExtract,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rc is None: return


    (
        (gBreps1, srf_devs_All),
        sLogs,
       ) = rc


    sc.doc.Views.RedrawEnabled = True

    return

    gBs_1F_Added = []
    nAddedBreps = nAddedBrepsForUntrimmedSrfs = 0
    
    nBreps0 = len(gBreps0)


    def knotMultiplicityList(knots):
        """Returns a list."""
        i = 0
        iMulties = []
        fKnotTs_Unique = []
        while True:
            knot = knots[i]
            fKnotTs_Unique.append(knot)
            iMulti = knots.KnotMultiplicity(index=i)
            iMulties.append(iMulti)
            #print("{} at {:.4f}".format(iMulti, knot),
            i += iMulti
            if i >= knots.Count:
                break
        return iMulties


    #rgMeshParams = rg.MeshingParameters.Minimal 
    #rgMeshParams = rg.MeshingParameters.FastRenderMesh
    rgMeshParams = rg.MeshingParameters.QualityRenderMesh
    
    stopwatch = Stopwatch()
    stopwatch.Start()
    
    for ib, gBrep0 in enumerate(gBreps0):
        sc.escape_test()
        
        Rhino.RhinoApp.SetCommandPrompt (
            "Processing brep {} of {}...".format(ib+1, nBreps0))
        
        rdBrep0 = sc.doc.Objects.FindId(gBrep0)
        rgBrep0 = rdBrep0.Geometry
        
        attr = rdBrep0.Attributes
        
        idx_rgFacesWithModifiedSrfs = []
        
        timeFace = timeFacePrev = stopwatch.Elapsed.TotalSeconds
        
        sRhCmdPrompt_Brep = Rhino.RhinoApp.CommandPrompt
        
        for iF in range(rgBrep0.Faces.Count):
            sc.escape_test()
            
            timeFace = stopwatch.Elapsed.TotalSeconds
            if timeFace - timeFacePrev > 1.0:
                Rhino.RhinoApp.SetCommandPrompt(
                    sRhCmdPrompt_Brep + "    Face {} of {}...".format(
                        iF, rgBrep0.Faces.Count))
                timeFacePrev = timeFace
            
            rgFace0 = rgBrep0.Faces[iF]
            rgSrf0 = rgFace0.UnderlyingSurface()
            
            # Only process NurbsSurfaces.  Do not convert to NurbsSurfaces.
            if isinstance(rgSrf0, rg.NurbsSurface):
                pass
            elif isinstance(rgSrf0, rg.SumSurface):
                def isSumSrfReducible(ss):
                    for iDir in 0, 1:
                        c = ss.IsoCurve(iDir, ss.Domain(iDir % 2).Min)
                        if isinstance(c, rg.NurbsCurve):
                            if c.Points.Count > c.Degree + 1:
                                return True
                    return False

                if not isSumSrfReducible(rgSrf0):
                    continue

                rgSrf0 = rgSrf0.ToNurbsSurface()
            else:
                continue


            # Remove/reduce knots.
            rc = processSrf(rgSrf0,
                    fDevTol, nKnotMulti_Min, nKnotMulti_Max, rgMeshParams,
                    bEcho=bEcho, bDebug=bDebug)
            if rc is None: continue


            if bEcho and not rc and nBreps0==1 and rgBrep0.Faces.Count==1: # For when rc is False, not None.
                s  = "Input surface Deg:{}, PtCt:{},{}, IsUniform:, IsClosed:".format(
                        rgSrf0.Degree,
                        rgSrf0.Points.CountU, rgSrf0.Points.CountV)
                s += "\nA surface could not be found for entered parameters:"
                s += " SrfSrfDev_Tol={}  MinMultiplicityToRemove={}" \
                        "  MaxMultiplicityToRemove={}".format(
                        fDevTol, nKnotMulti_Min, nKnotMulti_Max)
                print(s)
                continue
            
            rgSrf1, iCt_KnotsRemoved, fSrfDev = rc
            if bEcho and not rgSrf1.IsValid and nBreps0==1:
                print("New surface for {} is not valid.".format(gBrep0))
            
            # Duplicate rgFace0 to the brep that will be modified.
            rgBrep1 = rgFace0.DuplicateFace(False)
            rgFace1 = rgBrep1.Faces[0]
            rgBrep1.AddSurface(rgSrf1)
            rgFace1.ChangeSurface(1)
            rgBrep1.Compact()
            rgFace1.RebuildEdges(1e-5, True, True)
            
            if not rgBrep1.IsValid:
                s = "New brep for {} is not valid".format(gBrep0)
                gBrep_Added = sc.doc.Objects.AddBrep(rgBrep1, attr)
                if gBrep_Added == Guid.Empty:
                    s += " and could not be added to the document."
                else:
                    gBs_1F_Added.append(gBrep_Added)
                    s += "."
                
                if bSrfOnFail and sc.doc.Objects.AddSurface(rgSrf1) != Guid.Empty:
                    s += "\nUntrimmed surface added instead."
                    print(s)
                    nAddedBrepsForUntrimmedSrfs += 1
                else:
                    s += "\nUntrimmed surface also could not be added."
                    print(s)
                    continue
            
            else:
                gBrep_Added = sc.doc.Objects.AddBrep(rgBrep1, attr)
                if gBrep_Added != Guid.Empty:
                    bDocModified = True
                    gBs_1F_Added.append(gBrep_Added)
                    idx_rgFacesWithModifiedSrfs.append(iF)
                    nAddedBreps += 1
                else:
                    print("Brep for {} could not be added to"
                          " the document.".format(gBrep0))
            
            rgBrep1.Dispose()
            
            if bEcho and nBreps0 == 1 and rgBrep0.Faces.Count==1:
                s = ""
                s += "Knots removed: {}".format(iCt_KnotsRemoved)
                s += "    Surface deviation: {}".format(
                    formatDistance(fSrfDev))
                for iDir, sDir in zip((0,1), ('U','V')):
                    s += "    {} Dir:".format(sDir)
                    s += "{}:{}".format("Degree", rgSrf0.Degree(iDir))
                    knots0 = rgSrf0.KnotsV if iDir else rgSrf0.KnotsU
                    knots1 = rgSrf1.KnotsV if iDir else rgSrf1.KnotsU
                    s += "    {}KnotMultiplicities:{}->{}".format(
                        sDir,
                        ",".join(str(i) for i in knotMultiplicityList(knots0)),
                        ",".join(str(i) for i in knotMultiplicityList(knots1)))
                    s += "    {}:{}->{}".format("PtCt",
                            (rgSrf0.Points.CountU, rgSrf0.Points.CountV)[iDir],
                            (rgSrf1.Points.CountU, rgSrf1.Points.CountV)[iDir])
            #        s += "  {}:{}-{}".format("IsRational",
            #                str(rgSrf0.IsRational)[0],
            #                str(rgSrf1.IsRational)[0])
            #        s += "  {}:{}-{}".format("IsClosed",
            #                str(rgSrf0.IsClosed)[0],
            #                str(rgSrf1.IsClosed)[0])
            #        if rgSrf0.IsClosed or rgSrf1.IsClosed:
            #            s += "  {}:{}-{}".format("IsPeriodic",
            #                    str(rgSrf0.IsPeriodic)[0],
            #                    str(rgSrf1.IsPeriodic)[0])
                print(s)
        
        nFacesAffected = len(idx_rgFacesWithModifiedSrfs)
        
        # No faces had their surfaces modified.
        if nFacesAffected == 0: continue
        
        # If bReplace, remove affected faces from brep.
        if bReplace and nFacesAffected == rgBrep0.Faces.Count:
            sc.doc.Objects.Delete(gBrep0, not bEcho)
        elif bReplace and nFacesAffected != rgBrep0.Faces.Count:
            
            idx_rgFacesWithModifiedSrfs = sorted(idx_rgFacesWithModifiedSrfs,
                    reverse=True)
            
            for iF in idx_rgFacesWithModifiedSrfs:
                rgBrep0.Faces.RemoveAt(iF)
            
            sc.doc.Objects.Replace(gBrep0, rgBrep0)
        
        # Select new faces.
        gBs_1F_Added_NETList = List[Guid](gBs_1F_Added)
        nSelected = sc.doc.Objects.Select(gBs_1F_Added_NETList)
        
        if nSelected > 0:
            print("{} monoface breps of simplified surfaces are"
                  " selected.".format(nSelected))
        
    
    stopwatch.Stop()
    if bDebug:
        print("Runtime for main loop in main: {} seconds".format(
                stopwatch.Elapsed.TotalSeconds))
    
    if bEcho:
        s  = "{} out of {} brep faces had knots removed.".format(
                nAddedBreps, nBreps0)
        if nAddedBrepsForUntrimmedSrfs:
            s += "  {} untrimmed surfaces added instead.".format(
                    nAddedBrepsForUntrimmedSrfs)
        print(s)
    
    if bDocModified:
        sc.doc.Views.Redraw()


if __name__ == '__main__': main()