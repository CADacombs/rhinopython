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
211010: Now supports perspective viewports.
"""

import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


def zoomToObject(gObj, fMagFactor=1.0):
    rs.EnableRedraw(False)
    
    rs.UnselectAllObjects()
    rs.SelectObject(gObj)
    
    if rs.ObjectType(gObj) == 8192:
        # TextDot
        rs.ViewCameraTarget(target=rs.TextDotPoint(gObj))
    else:
        rs.ZoomSelected()
        viewport = sc.doc.Views.ActiveView.ActiveViewport
        viewport.Magnify(fMagFactor, mode=False)
    
    rs.EnableRedraw()


def main():
    
    stickyKey = 'fMagFactor({})'.format(__file__)
    if sc.sticky.has_key(stickyKey):
        fMagFactor = sc.sticky[stickyKey]
    else:
        fMagFactor = sc.sticky[stickyKey] = 1.0
    
    print "Magnification factor: {}".format(fMagFactor)
    
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
    gObj_current = gObjs[idxSelected]
    
    zoomToObject(gObj_current, fMagFactor)
    
    while len(gObjs) > 0:
        go.SetCommandPrompt(
                "Selected object {} of {}.  Press Esc when done.".format(
                idxSelected + 1, len(gObjs)))
        
        
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
        
        
        if sOption == 'All':
            rs.EnableRedraw()
            rs.SelectObjects(gObjs)
            rs.ZoomSelected()
            rs.EnableRedraw(False)
            continue
        elif sOption == 'Next':
            idxSelected += 1
        elif sOption == 'Previous':
            idxSelected -= 1
        elif sOption == 'ZoomIn':
            if rs.ObjectType(gObj_current) == 8192:
                # TextDot
                viewport = sc.doc.Views.ActiveView.ActiveViewport
                viewport.Magnify(2.0, mode=False)
            fMagFactor *= 2.0
            print "Magnification factor: {} -> {}".format(sc.sticky[stickyKey], fMagFactor)
            sc.sticky[stickyKey] = fMagFactor
        elif sOption == 'ZoomOut':
            if rs.ObjectType(gObj_current) == 8192:
                # TextDot
                viewport = sc.doc.Views.ActiveView.ActiveViewport
                viewport.Magnify(0.5, mode=False)
            fMagFactor *= 0.5
            print "Magnification factor: {} -> {}".format(sc.sticky[stickyKey], fMagFactor)
            sc.sticky[stickyKey] = fMagFactor
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
            gSource = rs.GetObject("Select object which to match its properties")
            if gSource:
                rs.MatchObjectAttributes(gObj_current, gSource)
                gObjs.remove(gObj_current)
        
        if len(gObjs) == 0:
            print "No more objects."
            return
        idxSelected = idxSelected % len(gObjs)
        gObj_current = gObjs[idxSelected]
        
        zoomToObject(gObj_current, fMagFactor)


if __name__ == '__main__': main()
