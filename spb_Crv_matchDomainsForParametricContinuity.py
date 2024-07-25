"""
Provided a curve to modify (A) and a reference curve (R),
this script will match the domain of A so that 
D(A) = D(R) * (L(A) / L(R))
where
D is the end span domain length
L is the distance between the end control point and next control point
"""

from __future__ import absolute_import, division, print_function, unicode_literals

"""
221016: Created.
240725: Bug fix when processing multi-spanned curves.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import rhinoscriptsyntax as rs
import scriptcontext as sc


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


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

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve with optional input.
    """

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select 2 curves near their ends")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    go.GeometryAttributeFilter = ri.Custom.GeometryAttributeFilter.WireCurve

    go.OneByOnePostSelect = True
    go.DisablePreSelect()

    idxs_Opt = {}
    def addOption(key): idxs_Opt[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opt.clear()

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

        for key in idxs_Opt:
            if go.Option().Index == idxs_Opt[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def formatDistance(fDistance, iPrecision=15):
    if fDistance is None: return "(No deviation provided)"
    
    if fDistance < 0.01:
        return "{:.{}e}".format(fDistance, iPrecision)
    
    if fDistance < 0.1:
        return "{:.{}g}".format(fDistance, iPrecision+1)
    
    return "{:.{}g}".format(fDistance, iPrecision)


def wasCurvePickedCloserToT1(objref_Crv):
    rgCrv = objref_Crv.Curve()

    bSuccess, t = rgCrv.ClosestPoint(objref_Crv.SelectionPoint())
    if not bSuccess:
        return

    return t >= rgCrv.Domain.Mid


def modifyCurveDomain(rgCrv_A_In, rgCrv_R, bEvalT1End_NotT0_A, bEvalT1End_NotT0_R, bEcho=True, bDebug=False):
    """
    Returns tuple of string, string: (strShortContinuityDescription, strLongContinuityDescription)
    """

    nc_A_In = rgCrv_A_In.ToNurbsCurve()
    nc_R = rgCrv_R.ToNurbsCurve()

    deg_A = nc_A_In.Degree
    deg_R = nc_R.Degree

    if bEvalT1End_NotT0_A:
        domain_EndSpan_A = nc_A_In.SpanDomain(nc_A_In.SpanCount - 1)
        nc_EndSpan_A = nc_A_In.Trim(domain_EndSpan_A)
        t_A = nc_EndSpan_A.Domain.T1
        crvEvalSide_A = rg.CurveEvaluationSide.Below
        idxCp_Pos_A = nc_EndSpan_A.Points.Count - 1
        idxCp_Tan_A = nc_EndSpan_A.Points.Count - 2
    else:
        domain_EndSpan_A = nc_A_In.SpanDomain(0)
        nc_EndSpan_A = nc_A_In.Trim(domain_EndSpan_A)
        t_A = nc_EndSpan_A.Domain.T0
        crvEvalSide_A = rg.CurveEvaluationSide.Above
        idxCp_Pos_A = 0
        idxCp_Tan_A = 1
    pt_Pos_A = nc_EndSpan_A.Points[idxCp_Pos_A].Location
    pt_Tan_A = nc_EndSpan_A.Points[idxCp_Tan_A].Location
    fDist_CPs_A = pt_Pos_A.DistanceTo(pt_Tan_A)
    #sc.doc.Objects.AddPoint(pt_Pos_A)
    #sc.doc.Objects.AddPoint(pt_Tan_A)

    if bEvalT1End_NotT0_R:
        domain_EndSpan_R = nc_R.SpanDomain(nc_R.SpanCount - 1)
        nc_EndSpan_R = nc_R.Trim(domain_EndSpan_R)
        t_B = nc_EndSpan_R.Domain.T1
        crvEvalSide_B = rg.CurveEvaluationSide.Below
        idxCp_Pos_R = nc_EndSpan_R.Points.Count - 1
        idxCp_Tan_R = nc_EndSpan_R.Points.Count - 2
    else:
        domain_EndSpan_R = nc_R.SpanDomain(0)
        nc_EndSpan_R = nc_R.Trim(domain_EndSpan_R)
        t_B = nc_EndSpan_R.Domain.T0
        crvEvalSide_B = rg.CurveEvaluationSide.Above
        idxCp_Pos_R = 0
        idxCp_Tan_R = 1
    pt_Pos_R = nc_EndSpan_R.Points[idxCp_Pos_R].Location
    pt_Tan_R = nc_EndSpan_R.Points[idxCp_Tan_R].Location
    fDist_CPs_R = pt_Pos_R.DistanceTo(pt_Tan_R)

    nc_A_In.Dispose()
    nc_R.Dispose()

    #sc.doc.Objects.AddCurve(nc_EndSpan_A)
    #sc.doc.Objects.AddCurve(nc_EndSpan_R)
    #sc.doc.Views.Redraw()
    #1/0


    # D(A) = D(R) * (L(A) / L(R))

    m = (
        (domain_EndSpan_R.Length / domain_EndSpan_A.Length) * 
        (fDist_CPs_A / fDist_CPs_R)
        )

    if bDebug:
        sEval = "domain_EndSpan_A.Length"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))
        sEval = "domain_EndSpan_R.Length"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))
        sEval = "fDist_CPs_A"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))
        sEval = "fDist_CPs_R"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))
        sEval = "m"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))
        sEval = "m - 1.0"; print("{}: {}".format(sEval, formatDistance(eval(sEval))))

    if abs(m - 1.0) == 0.0:
        if bEcho: print("Domains already exactly match.")
        return
    if abs(m - 1.0) <= 2**-53:
        if bEcho:
            print("Domains already match within {} (machine epsilon/2.0).".format(
                formatDistance(2**-53)))
        return
    if abs(m - 1.0) <= 2**-52:
        if bEcho:
            print("Domains already match within {} (machine epsilon).".format(
                formatDistance(2**-52)))
        return
    #if abs(m - 1.0) <= Rhino.RhinoMath.ZeroTolerance:
    #    if bEcho:
    #        print("Domains already match within {}.".format(
    #            formatDistance(abs(m - 1.0))))
    #    return

    if bEcho:
        print("Domain length multiplier to apply: {}".format(
            formatDistance(m)))

    length_Domain_A_Out = rgCrv_A_In.Domain.Length * m
    domain_In = rgCrv_A_In.Domain
    rgCrv_A_In.Domain = rg.Interval(0.0, length_Domain_A_Out)
    domain_Out = rgCrv_A_In.Domain
    return domain_In != domain_Out


def processCurveObjects(rhCrv_A_In, rhCrv_R, bEvalT1End_NotT0_A, bEvalT1End_NotT0_R, bEcho=True, bDebug=False):
    """
    """

    rdCrv_A = rs.coercerhinoobject(rhCrv_A_In) # IsDocumentControlled.
    crv_A_In = rdCrv_A.Geometry # IsDocumentControlled.
    crv_R = rs.coercecurve(rhCrv_R) # IsDocumentControlled.

    rc = modifyCurveDomain(
        crv_A_In,
        crv_R,
        bEvalT1End_NotT0_A,
        bEvalT1End_NotT0_R,
        bEcho=bEcho,
        bDebug=bDebug,
        )
    if rc is None:
        return
    if not rc:
        print("Curve domain could not be modified.")
        return

    if not rdCrv_A.CommitChanges():
        print("CurveObject was NOT modified.")
        return

    s = "CurveObject was modified"

    sc.doc.Objects.UnselectAll()
    if sc.doc.Objects.Select(rdCrv_A.Id):
        s += " and is selected."
    else:
        s += " but could not be selected."

    print(s)


def main():

    objrefs = getInput()
    if objrefs is None: return

    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    if not bDebug: sc.doc.Views.RedrawEnabled = False

    bEvalT1End_NotT0_A = wasCurvePickedCloserToT1(objrefs[0])
    bEvalT1End_NotT0_R = wasCurvePickedCloserToT1(objrefs[1])

    processCurveObjects(
        rhCrv_A_In=objrefs[0],
        rhCrv_R=objrefs[1],
        bEvalT1End_NotT0_A=bEvalT1End_NotT0_A,
        bEvalT1End_NotT0_R=bEvalT1End_NotT0_R,
        bEcho=bEcho,
        bDebug=bDebug,
        )

    sc.doc.Views.RedrawEnabled = True


if __name__ == '__main__': main()