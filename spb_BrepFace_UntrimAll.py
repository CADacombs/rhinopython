"""
190724: Created to avoid having to enter some options every time command is called with an alias.
...
200617: Added support for V7's PerFaceColor.
200619: Added support for invalid brep input.
200623: Now skips untrimmed faces.
200629: Import-related update.
210312: Bug fix.
210711: Now normal direction of face is maintained.
230831: More printed feedback added.
        Moved areAllTrimsEqualToSrfSides into this script.  Added string feedback to return of areAllTrimsEqualToSrfSides.
241231: Added more printed feedback.
250121: Bug fix in printed feedback for when input breps are invalid.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Input as ri
import scriptcontext as sc

from System import Guid

import xBlock
import xBrepObject

import spb_BrepFace_UnderlyingSrf


sOpts = (
        'bShrinkUnderlying',
        'bTryForPrimitive',
        'fTryForPrimitiveTol',
        'bNurbsAllowedForPrimitives',
        'bOnlyClosedSrfsForRoundPrimitives',
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
    
    key = 'bShrinkUnderlying'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bTryForPrimitive'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'fTryForPrimitiveTol'
    values[key] = 0.1 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key])
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)
    
    key = 'bNurbsAllowedForPrimitives'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bOnlyClosedSrfsForRoundPrimitives'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    # Load sticky.
    for key in stickyKeys:
        if stickyKeys[key] in sc.sticky:
            if key in riOpts:
                riOpts[key].CurrentValue = values[key] = sc.sticky[stickyKeys[key]]


    @classmethod
    def setValues(cls):
        for key in sOpts:
            if key in cls.riOpts:
                cls.values[key] = cls.riOpts[key].CurrentValue


    @classmethod
    def saveSticky(cls):
        for key in cls.stickyKeys:
            if key in cls.riOpts:
                sc.sticky[cls.stickyKeys[key]] = cls.riOpts[key].CurrentValue


def getInput():
    """
    """
    
    # Get face with optional input.
    
    go = ri.Custom.GetObject()
    
    go.SetCommandPrompt("Select face")
    
    go.GeometryFilter = rd.ObjectType.Surface
    
    go.AcceptNumber(True, acceptZero=True)
#    go.EnableHighlight(False)
    
    idxsOpts_Main = {}

    while True:
        go.AddOptionToggle(Opts.names['bShrinkUnderlying'], Opts.riOpts['bShrinkUnderlying'])
        go.AddOptionToggle(Opts.names['bTryForPrimitive'], Opts.riOpts['bTryForPrimitive'])
        if Opts.values['bTryForPrimitive']:
            go.AddOptionDouble(Opts.names['fTryForPrimitiveTol'], Opts.riOpts['fTryForPrimitiveTol'])
        go.AddOptionToggle(Opts.names['bNurbsAllowedForPrimitives'], Opts.riOpts['bNurbsAllowedForPrimitives'])
        go.AddOptionToggle(Opts.names['bOnlyClosedSrfsForRoundPrimitives'], Opts.riOpts['bOnlyClosedSrfsForRoundPrimitives'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)
        
        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return tuple([objrefs] + [Opts.values[key] for key in sOpts])
#            sc.doc.Objects.UnselectAll() # Prepare for repeat of go.Get().
#            sc.doc.Views.Redraw()
        elif res == ri.GetResult.Cancel:
            return
        
        # An option was selected or a number was entered.
        
        res = go.Result()
        
        key = 'fTryForPrimitiveTol'
        if res == ri.GetResult.Number:
            Opts.riOpts[key].CurrentValue = go.Number()
        if Opts.riOpts[key].CurrentValue < 0.0:
            Opts.riOpts[key].CurrentValue = Opts.riOpts[key].InitialValue
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def areAllTrimsEqualToSrfSides(rgFace):

    rgF = rgFace

    if not rgF.IsSurface:
        return False, "Not IsSurface."


    if rgF.Brep.Trims.Count != 4:
        return False, "Face.Brep.Trims.Count == {}".format(rgF.Brep.Trims.Count)


    iCt_NotSingular = 0
    for rgT in rgF.Brep.Trims:
        if rgT.TrimType != Rhino.Geometry.BrepTrimType.Singular:
            iCt_NotSingular += 1
    
    idxEs = rgF.AdjacentEdges()

    iCt_FaceEdges = len(idxEs)
    
    if iCt_NotSingular != iCt_FaceEdges:
        return False, "iCt_NotSingular != iCt_FaceEdges"

    for idxE in idxEs:
        rgE = rgF.Brep.Edges[idxE]
        if rgE.Tolerance != 0.0:
            return False, "Non-0.0 edge tolerance found."

    return spb_BrepFace_UnderlyingSrf.areAllTrimsOnSrfSides(rgF), None


def main():
    """
    """
    
    rc = getInput()
    if rc is None or not rc[0]: return
    (
        objrefs,
        bShrinkUnderlying,
        bTryForPrimitive,
        fTryForPrimitiveTol,
        bNurbsAllowedForPrimitives,
        bOnlyClosedSrfsForRoundPrimitives,
        bEcho,
        bDebug,
        ) = rc
    
    for objref in objrefs:
        rdBrep_In = objref.Object()
        rgFace_In = objref.Face()

        if rdBrep_In.IsValid:
            # Test for special case that face is already untrimmed and
            # is at desired parameters.
            bTrimsEqSrfSides, sWhy = areAllTrimsEqualToSrfSides(
                rgFace_In)
            if bTrimsEqSrfSides:
                print "Face boundary is already its underlying surface sides."
                continue # to next face.
            
            if bDebug or (bEcho and len(objrefs)) == 1:
                print "Face was found to not be just its underlying surface because " + sWhy
            
            rc = spb_BrepFace_UnderlyingSrf.getOptimalSurface(
                    rgFace_In,
                    bShrinkPlanar=bShrinkUnderlying,
                    bShrinkRoundPrims=bShrinkUnderlying,
                    bShrinkOthers=bShrinkUnderlying,
                    bTryForPrimitive=bTryForPrimitive,
                    fTolForPrimitive=fTryForPrimitiveTol,
                    bOnlyClosedSrfsForRoundPrimitives=bOnlyClosedSrfsForRoundPrimitives,
                    bUseShrunkForPrimIfNeeded=False,
                    bEcho=bEcho,
                    bDebug=bDebug)
            
            if rc[0]:
                rgSrf_Res, fTol_Used, bShrunkUsed = rc[0]
        else:
            print "Input brep is invalid, so only getting UnderlyingSurface."
            rgSrf_Res = rgFace_In.UnderlyingSurface()
            fTol_Used = None
            bShrunkUsed = False

        attr = rdBrep_In.Attributes
        
        if not bDebug: sc.doc.Views.RedrawEnabled = False

        rgB_1F_Out = rgSrf_Res.ToBrep()
        rgSrf_Res.Dispose()
        if not rgB_1F_Out.IsSolid and rgFace_In.OrientationIsReversed:
            rgB_1F_Out.Flip()

        if Rhino.RhinoApp.ExeVersion >= 7:
            rgB_1F_Out.Faces[0].PerFaceColor = rgFace_In.PerFaceColor

        gBrep_Out = sc.doc.Objects.AddBrep(rgB_1F_Out, attr)
        if gBrep_Out != Guid.Empty:
            if bEcho:
                sSrfType = rgSrf_Res.GetType().Name
                s  = "{} was added to document.".format(sSrfType)
                if fTol_Used:
                    s += " Face {} to be shrunk to obtain the shape at a tolerance of {}.".format(
                        "needed" if bShrunkUsed else "did not need",
                        fTol_Used)
                print(s)
            
            if rdBrep_In.GetType() == rd.BrepObject:
                idx_GCI = objref.GeometryComponentIndex.Index
                idx_Face = 0 if idx_GCI == -1 else idx_GCI
                if xBrepObject.removeFaces(rdBrep_In, idx_Face) is not None:
                    print "Input face was removed."
            else:
                if bEcho:
                    print "New face could not be added to the document." \
                          "Face was not removed from {}.".format(rdBrep_In.GetType().Name)
        else:
            if bEcho: print "Surface could not be added to document."

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()