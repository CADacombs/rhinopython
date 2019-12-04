"""
190610: Created, starting with another script.
"""

import Rhino
import Rhino.DocObjects as rd
import Rhino.Geometry as rg
import scriptcontext as sc

import xKnotList


def main(bEcho=True, bDebug=False):
    """
    """
    
    # Custom geometry filter to only select non-NurbsCurve wires.
    def nurbsCurveGeometryFilter(rhObject, geom, compIdx):
        if not isinstance(geom, rg.Curve): return False
        return geom.GetType() == rg.NurbsCurve
    
    Rhino.RhinoApp.SetCommandPrompt("Searching ...")
    
    gObjs_Preselected = [rdObj.Id for rdObj in sc.doc.Objects.GetSelectedObjects(includeLights=False, includeGrips=False)]
    
    gCrvs_Target = []
    iter = rd.ObjectEnumeratorSettings()
    iter.NormalObjects = True
    iter.LockedObjects = False
    iter.IncludeLights = False
    iter.IncludeGrips = False
    for rdRhinoObject in sc.doc.Objects.GetObjectList(iter):
        if rdRhinoObject.Id in gObjs_Preselected:
            continue
        if rdRhinoObject.ObjectType != rd.ObjectType.Curve:
            continue
        rgCrv = rdRhinoObject.CurveGeometry
        sCrvType = rgCrv.GetType().Name
        if sCrvType == 'NurbsCurve':
            if not xKnotList.isUniform(rgCrv.Knots):
                gCrvs_Target.append(rdRhinoObject.Id)
        rgCrv.Dispose()

    ct_Selected  = 0

    sc.doc.Views.RedrawEnabled = False
    for gCrv in gCrvs_Target:
        if sc.doc.Objects.Select(gCrv):
            ct_Selected += 1
    sc.doc.Views.RedrawEnabled = True

    if bEcho:
        if gCrvs_Target:
            s  = "{} curve{} added to selection.".format(
                    ct_Selected,
                    '' if ct_Selected == 1 else 's')
            print s
        else:
            print "\nNo curves added to selection."
    
    return gCrvs_Target


if __name__ == '__main__': main(bEcho=bool(1), bDebug=bool(0))