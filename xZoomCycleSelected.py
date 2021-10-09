"""
170929: Created to replace a RhinoScript and with more options.
171011: ...
180403: Now objects made red are removed from the selection set.
180626: Added ChangeLayer option.
180914: Fixed bug for ChangeLayer option. 
191019: Added MatchProps option.
210706: Added numerical input for jumping to an index+1.  Added ZoomIn and ZoomOut options.
211009a: Bug fix for zooming into TextDot, which ZoomSelected differently than other objects.
211009b: Added print statements reporting the zoom factor.  Commented out debug print statements.
"""

import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


def zoomToObject(gObj, fZoomFactor=1.0):
    rs.EnableRedraw(False)
    rs.UnselectAllObjects()
    rs.SelectObject(gObj)
    #print fZoomFactor, rs.ViewRadius(),
    if rs.ObjectType(gObj) == 8192:
        # TextDot
        rs.ViewCameraTarget(target=rs.TextDotPoint(gObj))
    else:
        rs.ZoomSelected()
        rs.ViewRadius(radius=rs.ViewRadius()*fZoomFactor)
    #print rs.ViewRadius()
    rs.EnableRedraw()


def main():
    
    stickyKey = 'fZoomFactor({})({})'.format(__file__, sc.doc.Name)
    fZoomFactor = sc.sticky[stickyKey] if sc.sticky.has_key(stickyKey) else 1.0
    print "Zoom factor: {}".format(fZoomFactor)
    
    gObjs = rs.GetObjects("Select objects", preselect=True)
    if gObjs is None: return
    
    idxSelected = 0
    
    go = ri.Custom.GetOption()
    
    go.AcceptNumber(True, acceptZero=True)
    
    sOptions = (
        'All',
        'Current',
        'Next',
        'Previous',
        'ZoomIn',
        'ZoomOut',
        'HighlightAll',
        'ColorRed',
        'Hide',
        'Lock',
        'Delete',
        'ChangeLayer',
        'MatchProps'
    )
    go.SetDefaultString(sOptions[2])
    for s in sOptions : go.AddOption(s)
    
    # Set first selection and view.
    sView = rs.CurrentView()
    gObj_current = gObjs[idxSelected]
    
    zoomToObject(gObj_current, fZoomFactor)
    
    while len(gObjs) > 0:
        go.SetCommandPrompt(
                "Selected object {} of {}.  Press Esc when done.".format(
                idxSelected + 1, len(gObjs)))
        
        sView = rs.CurrentView()
        
        #rs.EnableRedraw()
        res = go.Get()
        if res == ri.GetResult.Cancel:
            rs.SelectObjects(gObjs)
            return None
        elif res == ri.GetResult.Option:
            sOption = go.Option().EnglishName
        elif res == ri.GetResult.Number:
            iRes = int(go.Number())
            if 0 < iRes <= len(gObjs):
                idxSelected = iRes - 1
            sOption = None
        else:
            sOption = go.StringResult()
            print sOption
        
        #rs.EnableRedraw(False)
        
        if sOption == 'All':
            rs.SelectObjects(gObjs)
            rs.ZoomSelected(sView)
            continue
        elif sOption == 'Next':
            idxSelected += 1
        elif sOption == 'Previous':
            idxSelected -= 1
        elif sOption == 'ZoomIn':
            if rs.ObjectType(gObj_current) == 8192:
                # TextDot
                rs.ViewRadius(radius=0.5*rs.ViewRadius())
            fZoomFactor *= 0.5
            print "Zoom factor: {} -> {}".format(sc.sticky[stickyKey], fZoomFactor)
            sc.sticky[stickyKey] = fZoomFactor
        elif sOption == 'ZoomOut':
            if rs.ObjectType(gObj_current) == 8192:
                # TextDot
                rs.ViewRadius(radius=2.0*rs.ViewRadius())
            fZoomFactor *= 2.0
            print "Zoom factor: {} -> {}".format(sc.sticky[stickyKey], fZoomFactor)
            sc.sticky[stickyKey] = fZoomFactor
        elif sOption == 'HighlightAll':
            rs.SelectObjects(gObjs)
            continue
        elif sOption == 'ColorRed':
            rs.ObjectColorSource(gObj_current, 1)
            rs.ObjectColor(gObj_current, rs.coercecolor((255,0,0)))
            gObjs.remove(gObj_current)
        elif sOption == 'Hide':
            rs.HideObject(gObj_current)
            gObjs.remove(gObj_current)
        elif sOption == 'Lock':
            rs.LockObject(gObj_current)
            gObjs.remove(gObj_current)
        elif sOption == 'Delete':
            rs.DeleteObject(gObj_current)
            gObjs.remove(gObj_current)
        elif sOption == 'ChangeLayer':
            layer = rs.GetLayer()
            if layer is not None:
                rs.ObjectLayer(gObj_current, layer)
                gObjs.remove(gObj_current)
        elif sOption == 'MatchProps':
            #rs.EnableRedraw()
            gSource = rs.GetObject("Select object which to match its properties")
            if gSource:
                rs.MatchObjectAttributes(gObj_current, gSource)
                gObjs.remove(gObj_current)
            #rs.EnableRedraw(False)
        
        if len(gObjs) == 0:
            print "No more objects."
            return
        idxSelected = idxSelected % len(gObjs)
        gObj_current = gObjs[idxSelected]
        
        zoomToObject(gObj_current, fZoomFactor)


if __name__ == '__main__': main()