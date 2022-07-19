"""
https://mcneel.myjetbrains.com/youtrack/issue/RH-64854

Per https://developer.rhino3d.com/api/RhinoCommon/html/M_Rhino_Geometry_Curve_CreateArcLineArcBlend.htm :

CreateArcLineArcBlend(...)
        CreateArcLineArcBlend(startPt: Point3d, startDir: Vector3d, endPt: Point3d, endDir: Vector3d, radius: float) -> Curve
        
            Creates an arc-line-arc blend curve between two curves.
                    The output 
             is generally a PolyCurve with three segments: arc, line, arc.
                    In 
             some cases, one or more of those segments will be absent because they would 
             have 0 length. 
                    If there is only a single segment, the result will 
             either be an ArcCurve or a LineCurve.
        
        
            startPt: Start of the blend curve.
            startDir: Start direction of the blend curve.
            endPt: End of the blend curve.
            endDir: End direction of the arc blend curve.
            radius: The radius of the arc segments.
            Returns: The blend curve if successful, false otherwise.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
220121: Created.
220201: Bug fix.
220501: Cleaned.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


class Opts():
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fRadius'; keys.append(key)
    values[key] = 1.0
    riOpts[key] = ri.Custom.OptionDouble(values[key])
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
        else:
            idxOpt = go.AddOptionList(
                englishOptionName=cls.names[key],
                listValues=cls.listValues[key],
                listCurrentIndex=cls.values[key])

        return idxOpt


    @classmethod
    def setValue(cls, key, idxList=None):

        if key == 'fRadius':
            if cls.riOpts[key].CurrentValue < 2.0*sc.doc.ModelAbsoluteTolerance:
                cls.riOpts[key].CurrentValue = 0.0

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curves with optional input.
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select first curve near end")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve

    go.DisablePreSelect()
    go.OneByOnePostSelect = True

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fRadius')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=2, maximumNumber=2)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fRadius'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def createCurve(rgCrv_A, rgCrv_B, bT1WorkEnd_A, bT1WorkEnd_B, fRadius, bEcho=True, bDebug=False):
    """
    Returns:
        Single PolyCurve on success.
        None on fail.
    """

    rgC_A = rgCrv_A
    rgC_B = rgCrv_B

    t_WorkEnd_A = rgCrv_A.Domain.T1 if bT1WorkEnd_A else rgCrv_A.Domain.T0
    t_WorkEnd_B = rgCrv_B.Domain.T1 if bT1WorkEnd_B else rgCrv_B.Domain.T0


    ptA = rgC_A.PointAt(t_WorkEnd_A)
    ptB = rgC_B.PointAt(t_WorkEnd_B)


    def doCurveEndsMeet():
        dist = ptA.DistanceTo(ptB)
        if dist < 100.0 * sc.doc.ModelAbsoluteTolerance:
            return True
        return False

    if doCurveEndsMeet():
        if bEcho:
            print("Ends are too close to one another for a meaningful result."
                  "Script canceled.")
        return


    dirA = rgCrv_A.TangentAt(t_WorkEnd_A)
    dirB = rgCrv_B.TangentAt(t_WorkEnd_B)

    if not bT1WorkEnd_A: dirA = -dirA 
    if bT1WorkEnd_B: dirB = -dirB 


    ####

    return rg.Curve.CreateArcLineArcBlend(
        startPt=ptA,
        startDir=dirA,
        endPt=ptB,
        endDir=dirB,
        radius=fRadius)

    ####


def createCurveObject(objrefA, objrefB, fRadius, bEcho=True, bDebug=False):

    rgC_A = objrefA.Curve()

    bSuccess, tA = rgC_A.ClosestPoint(objrefA.SelectionPoint())
    if not bSuccess:
        rgC_A.Dispose()
        return


    rgC_B = objrefB.Curve()
    bSuccess, tB = rgC_B.ClosestPoint(objrefB.SelectionPoint())
    if not bSuccess:
        rgC_A.Dispose()
        rgC_B.Dispose()
        return

    bT1WorkEnd_A = tA > rgC_A.Domain.Mid
    bT1WorkEnd_B = tB > rgC_B.Domain.Mid

    rc = createCurve(
            rgCrv_A=rgC_A,
            rgCrv_B=rgC_B,
            bT1WorkEnd_A=bT1WorkEnd_A,
            bT1WorkEnd_B=bT1WorkEnd_B,
            fRadius=fRadius,
            bEcho=bEcho,
            bDebug=bDebug,
    )
    if rc is None:
        print("Curve was not created.")
        return

    rgC_Res = rc

    gCrv_Out = sc.doc.Objects.AddCurve(rgC_Res)
    if gCrv_Out == gCrv_Out.Empty:
        if bEcho:
            print("Curve could not be added.")
        return

    if bEcho: print("{} was added.".format(rgC_Res.GetType().Name))
    sc.doc.Views.Redraw()

    return gCrv_Out


def main():
    
    rc = getInput()
    if rc is None: return
    objrefA, objrefB = rc

    fRadius = Opts.values['fRadius']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    createCurveObject(
        objrefA,
        objrefB,
        fRadius=fRadius,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()