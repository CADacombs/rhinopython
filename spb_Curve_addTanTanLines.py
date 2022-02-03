"""
Creates tangent-tangent lines along curve with itself.
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
211101: Created.
220201: Added more options.
        Now, equal but opposing lines are treated as duplicates.
"""

import Rhino
import Rhino.Geometry as rg
import Rhino.Input as ri
import scriptcontext as sc


lengths_skipped = [] # Global log.


class Opts():

    keys = []
    values = {}
    names = {}
    riOpts = {}
    listValues = {}
    stickyKeys = {}


    key = 'fMinLineLength'; keys.append(key)
    values[key] = 1000.0 * sc.doc.ModelAbsoluteTolerance
    riOpts[key] = ri.Custom.OptionDouble(values[key], True, 1e-6)
    stickyKeys[key] = '{}({})({})'.format(key, __file__, sc.doc.Name)

    key = 'bKnots_NotGrevs'; keys.append(key)
    values[key] = False
    names[key] = 'SeedPts'
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'Greville', 'Knot')
    stickyKeys[key] = '{}({})'.format(key, __file__)

    key = 'fEpsilon'; keys.append(key)
    values[key] = 1e-5
    riOpts[key] = ri.Custom.OptionDouble(values[key], True, Rhino.RhinoMath.ZeroTolerance)
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

        if key == 'fMinLineLength':
            if cls.riOpts[key].CurrentValue < 1e-6:
                cls.riOpts[key].CurrentValue = cls.riOpts[key].InitialValue

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput():
    """
    Get curve and parameter with optional input.
    """

    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select curves")

    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve # Includes selection of BrepsEdges.


    def customGeometryFilter(rdObj, rgObj, compIdx):
        if not rgObj.IsPlanar(tolerance=Opts.values['fEpsilon']):
            return False

        return True

    go.SetCustomGeometryFilter(customGeometryFilter)


    go.AcceptNumber(True, acceptZero=True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('fMinLineLength')
        addOption('bKnots_NotGrevs')
        addOption('fEpsilon')
        addOption('bEcho')
        addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        if res == ri.GetResult.Number:
            key = 'fMinLineLength'
            Opts.riOpts[key].CurrentValue = go.Number()
            Opts.setValue(key)
            continue

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def uniqueKnots(nc):
    
    iK = 0
    ts_Unique = []
    
    while iK < nc.Knots.Count:
        k = nc.Knots[iK]
        m = nc.Knots.KnotMultiplicity(index=iK)

        if k < nc.Domain.T0:
            # Overlapping knot from periodic curve.
            iK += m
            continue

        if k > nc.Domain.T1:
            # Overlapping knot from periodic curve.
            break

        if m < nc.Degree:
            ts_Unique.append(k)

        iK += m

    return ts_Unique


def createLines(crv, fMinLineLength, bKnots_NotGrevs, fEpsilon, bEcho=True, bDebug=False):
    """
    """
    
    nc = crv.ToNurbsCurve()
    
    if bKnots_NotGrevs:
        ts = uniqueKnots(nc)
    else:
        ts = list(nc.GrevilleParameters())
        if not nc.IsPeriodic:
            ts = ts[1:-1]
    
    lines = [] # for duplicate testing.
    
    for a, tA in enumerate(ts):
        for tB in ts[a+1:]:
            rc = rg.Line.TryCreateBetweenCurves(
                nc, nc, tA, tB, perpendicular0=False, perpendicular1=False)
            
            if not rc[0]: continue
            
            tA_Out, tB_Out, line = rc[1:]
            
            line_reversed = rg.Line(line.To, line.From)
            
            for lineInLines in lines:
                if line.EpsilonEquals(lineInLines, epsilon=fEpsilon):
                    break
                if line_reversed.EpsilonEquals(lineInLines, epsilon=fEpsilon):
                    break

            else:
                # No match found.
                if line.Length < fMinLineLength:
                    lengths_skipped.append(line.Length)
                    continue

                lines.append(line)

    return lines


def main():

    objrefs = getInput()
    if objrefs is None: return

    fMinLineLength = Opts.values['fMinLineLength']
    bKnots_NotGrevs = Opts.values['bKnots_NotGrevs']
    fEpsilon = Opts.values['fEpsilon']
    bEcho = Opts.values['bEcho']
    bDebug = Opts.values['bDebug']

    Rhino.RhinoApp.SetCommandPrompt("Working ...")

    gLines_Out = []

    for objref in objrefs:
        crv = objref.Curve()
        
        lines = createLines(crv, fMinLineLength, bKnots_NotGrevs, fEpsilon, bEcho, bDebug)
        if not lines: continue

        gLines_Out.extend(sc.doc.Objects.AddLine(line) for line in lines)
    
    sPrint = ["{} lines added.".format(len(gLines_Out) if gLines_Out else "No")]
    
    if len(lengths_skipped) > 0:
        sPrint.append("Skipped {0:} lines of lengths {1:.{3:}f} through {2:.{3:}f}.".format(
            len(lengths_skipped),
            min(lengths_skipped),
            max(lengths_skipped),
            sc.doc.ModelDistanceDisplayPrecision))
    
    print('  '.join(sPrint))
    
    sc.doc.Views.Redraw()


if __name__ == '__main__': main()