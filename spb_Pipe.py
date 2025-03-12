"""
Observations of Brep.CreatePipe:
    localBlending setting doesn't seem to have a significant affect when pipes of
        a single radius is created.  Regardless _SelDup notices differences,
        possibly including control point locations.
    fitRail setting has no effect on monocurves, e.i., lines, arc, and NURBS.
    fitRail=True can simplify linear segments in polylines and polycurves.
"""

#! python 2  Must be on a line number less than 32.
from __future__ import absolute_import, division, print_function, unicode_literals

"""
171019: Created.
190830: Modified how numbers are entered at command prompt.
220821: Modified an option default value.
231030: Modified an option default value.
231106: Added printed info of result.
250107: Corrected printed output of pipe's face count.
250301: Bug fix in printed output.

TODO: Create pipes at multiple tolerances?
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc

import math

from System import Enum


def getInput():
    
    # Load sticky.
    stickyKeys = [
            'bForDraftSilhouette({})'.format(__file__),
            'fRadius({})'.format(__file__, sc.doc.Name),
            'fDraft_Deg({})'.format(__file__),
            'fDist({})({})'.format(__file__, sc.doc.Name),
            'bLocalBlending({})'.format(__file__),
            'idxPipeCapMode({})'.format(__file__),
            'bFitRail({})'.format(__file__),
            'fDistTol({})'.format(__file__, sc.doc.Name),
            'fAngTol_Deg({})'.format(__file__),
            ]
    stickyValues = [
            False,
            1.0,
            20.0,
            10.0,
            False,
            0,
            False,
            1.0*sc.doc.ModelAbsoluteTolerance,
            1.0*sc.doc.ModelAngleToleranceDegrees,
    ] # Default values.
    for i, stickyKey in enumerate(stickyKeys):
        if sc.sticky.has_key(stickyKey): stickyValues[i] = sc.sticky[stickyKey]
    (
            bForDraftSilhouette,
            fRadius,
            fDraft_Deg,
            fDist,
            bLocalBlending,
            idxPipeCapMode,
            bFitRail,
            fDistTol,
            fAngTol_Deg
    ) = stickyValues
    
    # Get curves with optional input.
    
    optionListIndices = {'PipeCapMode': None}
    
    def addOptions():
        go.AddOptionToggle('ForDraftSilhouette', optT_bForDraftSilhouette)
        if bForDraftSilhouette:
            go.AddOptionDouble('Draft', optD_fDraft_Deg)
            go.AddOptionDouble('Dist', optD_fDist)
        else:
            go.AddOptionDouble('Radius', optD_fRadius)
        go.AddOptionToggle('Blending', optT_bLocalBlending)
        optionListIndices['PipeCapMode'] = go.AddOptionList(
                "CapMode", sPipeCapModeTypes, idxPipeCapMode)
        go.AddOptionToggle('PolyCrvFitRail', optT_bFitRail)
        go.AddOptionDouble('DistTol', optD_fDistTol)
        go.AddOptionDouble('AngleTol', optD_fAngTol_Deg)
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select rail curves")
    
    go.GeometryFilter = rd.ObjectType.Curve
    
    optT_bForDraftSilhouette = ri.Custom.OptionToggle(
            bForDraftSilhouette, 'No', 'Yes')
    optD_fRadius = ri.Custom.OptionDouble(fRadius, True, 1e-6)
    optD_fDraft_Deg = ri.Custom.OptionDouble(fDraft_Deg, True,
            Rhino.RhinoMath.ZeroTolerance)
    optD_fDist = ri.Custom.OptionDouble(fDist, True, 0.0)
    optT_bLocalBlending = ri.Custom.OptionToggle(
            bLocalBlending, 'Global', 'Local')
    optT_bFitRail = ri.Custom.OptionToggle(bFitRail, 'No', 'Yes')
    optD_fDistTol = ri.Custom.OptionDouble(fDistTol, True, 0.0)
    optD_fAngTol_Deg = ri.Custom.OptionDouble(fAngTol_Deg, True,
            Rhino.RhinoMath.ZeroTolerance)
    
    # Create simplified strings of rg.PipeCapMode enumerations.
    aPipeCapModeTypes_Strings = Enum.GetNames(rg.PipeCapMode)
    sPipeCapModeTypes = []
    for idx in range(len(aPipeCapModeTypes_Strings)):
        for i, c in enumerate(aPipeCapModeTypes_Strings[idx]):
            if c.isupper(): idxLast = i
        sPipeCapModeTypes.append(aPipeCapModeTypes_Strings[idx][idxLast:])
    
    addOptions()
    
    go.AcceptNumber(True, True)
    
    go.AlreadySelectedObjectSelect = True
    go.DeselectAllBeforePostSelect = False # So objects won't be deselected on repeats of While loop.
    #go.SubObjectSelect = False
    go.EnableClearObjectsOnEntry(False) # Keep objects in go on repeats of While loop.
    go.EnableUnselectObjectsOnExit(False)
    
    bPreselectedObjsChecked = False
    
    print("Dist is the distance the pipe will be moved away from the rail.")
    
    while True:
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        # Use bPreselectedObjsChecked so that only objects before the
        # first call to go.GetMultiple is considered.
        if not bPreselectedObjsChecked and go.ObjectsWerePreselected:
            bPreselectedObjsChecked = True
            go.EnablePreSelect(False, True)
            continue
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        
        if res == ri.GetResult.Option:
            idxOptPicked = go.OptionIndex()
            
            if idxOptPicked == optionListIndices['PipeCapMode']:
                idxPipeCapMode = go.Option().CurrentListOptionIndex
            
            bForDraftSilhouette = optT_bForDraftSilhouette.CurrentValue
            fRadius = optD_fRadius.CurrentValue
            fDraft_Deg = optD_fDraft_Deg.CurrentValue
            fDist = optD_fDist.CurrentValue
            bLocalBlending = optT_bLocalBlending.CurrentValue
            bFitRail = optT_bFitRail.CurrentValue
            fDistTol = optD_fDistTol.CurrentValue
            fAngTol_Deg = optD_fAngTol_Deg.CurrentValue
        elif bForDraftSilhouette and res == ri.GetResult.Number:
            fDraft_Deg = optD_fDraft_Deg.CurrentValue = abs(go.Number())
        elif not bForDraftSilhouette and res == ri.GetResult.Number:
            fRadius = optD_fRadius.CurrentValue = abs(go.Number())
        else: break
        
        # Save sticky.
        stickyValues = (
                bForDraftSilhouette,
                fRadius,
                fDraft_Deg,
                fDist,
                bLocalBlending,
                idxPipeCapMode,
                bFitRail,
                fDistTol,
                fAngTol_Deg
        )
        for i, stickyKey in enumerate(stickyKeys):
            sc.sticky[stickyKey] = stickyValues[i]
        
        go.ClearCommandOptions()
        addOptions()
    
    # Get curve geometry.
    rgCrvs = []
    for i in range(go.ObjectCount):
        rgCrv = go.Object(i).Curve()
        if not rgCrv:
            print("Error!  Curve could not be obtained.  Will be skipped.")
            continue
        rgCrvs.append(rgCrv)
    
    go.Dispose()
    
    pipeCapMode = Enum.ToObject(rg.PipeCapMode, idxPipeCapMode)
    
    return (rgCrvs, bForDraftSilhouette, fRadius, fDraft_Deg, fDist,
            bLocalBlending, pipeCapMode, bFitRail,
            fDistTol, fAngTol_Deg)


def main():
    
    ret = getInput()
    if ret is None: return
    (rgCrvs, bForDraftSilhouette, fRadius, fDraft_Deg, fDist,
            bLocalBlending, pipeCapMode, bFitRail,
            fDistTol, fAngTol_Deg) = ret
    
    rgBreps_Pipe_All = []
    fLengths_All = []
    
    nIters = 1
    
    cw = 8
#    print "{:<{cw}}{:<{cw}}{:<{cw}}".format(
#            "Tol", "CrvCt", "TotalLn", cw=cw)
    
    fModelTol = sc.doc.ModelAbsoluteTolerance
    fAngleTolRadians = sc.doc.ModelAngleToleranceRadians
    
    fAbsoluteTols = [fModelTol * 10**(i+1) for i in xrange(nIters)]
    fAngleTolsRadians = [fAngleTolRadians * 10**(i) for i in xrange(1)]
    
    if bForDraftSilhouette:
        fRadius = math.sin(math.radians(fDraft_Deg)) * fDist
    
    for fAbsoluteTol in fAbsoluteTols:
        for fAngleTolRadians in fAngleTolsRadians:
            for rgCrv in rgCrvs:
                # Create the pipe.
                rgBreps_Pipe = rg.Brep.CreatePipe(
                    rail=rgCrv,
                    radius=fRadius,
                    localBlending=bLocalBlending,
                    cap=pipeCapMode,
                    fitRail=bFitRail,
                    absoluteTolerance=fDistTol,
                    angleToleranceRadians=Rhino.RhinoMath.ToRadians(fAngTol_Deg)
                    )
                if rgBreps_Pipe.Count > 0:
                    rgBreps_Pipe_All.extend(rgBreps_Pipe)
            #for rgCrv in rgCrvs: sc.doc.Objects.AddCurve(rgCrv)
            
    #            # Print stats of this iteration.
    #            fLength = 0.0
    #            for rgCrv in rgCrvs:
    #                fLength += rgCrv.GetLength()
    #            fLengths_All.append(fLength)
    #            print "{:<{cw}}{:<{cw}}{:<{cw}}".format(
    #                    fAbsoluteTol, rgCrvs.Count, fLength, cw=cw)

    #    i_Use = 0
    #    fLength_Use = fLengths_All[i_Use]
    #    for i in xrange(1, nIters):
    #        if fLengths_All[i] > 0.99 * fLength_Use:
    #            i_Use = i
    #            fLength_Use = fLengths_All[i]

    if not rgBreps_Pipe_All:
        print("No pipes were created.")
        return

    iCt_Faces = 0

    for rgB in rgBreps_Pipe_All:
        rgB.Faces.SplitKinkyFaces()
        #rg.Brep.Compact()
        g = sc.doc.Objects.AddBrep(rgB) # Missing splitKinkySurfaces parameter will split them.
        if g is not g.Empty:
            rdB = sc.doc.Objects.FindId(g)
            if rgB.Faces.Count != rdB.BrepGeometry.Faces.Count:
                print("Faces count of Geometry.Brep: {}".format(
                    rgB.Faces.Count))
                print("Faces count of BrepObject added to document: {}".format(
                    rdB.BrepGeometry.Faces.Count))
            iCt_Faces += rdB.BrepGeometry.Faces.Count

    print("{} breps with a total of {} faces added.".format(len(rgBreps_Pipe_All), iCt_Faces))

    sc.doc.Views.Redraw()


if __name__ == '__main__': main()