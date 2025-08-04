"""
190708-10: Created.
190715: Refactored.
190725: Bug fixes.
200420: Increased number of decimal places in printed output.
210206: Corrected typo.
250803: Modified float display from e to G formatting.
"""

import Rhino
import Rhino.Input as ri
import scriptcontext as sc


sOpts = (
        'bAddPts',
        'bNormalize',
        'bEcho',
        'bDebug',
)


class Opts:
    
    keys = []
    values = {}
    names = {}
    riOpts = {}
    stickyKeys = {}
    
    for key in sOpts:
        keys.append(key)
        names[key] = key[1:] # Overwrite as wanted in the following.
    
    key = 'bAddPts'
    values[key] = False
    names[key] = 'CreatePoint'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bNormalize'
    values[key] = False
    names[key] = 'InterpretNumericEntryAsNormalized'
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bEcho'
    values[key] = True
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
    stickyKeys[key] = '{}({})'.format(key, __file__)
    
    key = 'bDebug'
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(initialValue=values[key], offValue='No', onValue='Yes')
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


def getInput_Curve():
    """Get curve with optional input."""
    
    go = ri.Custom.GetObject()
    go.SetCommandPrompt("Select curve")
    
    go.GeometryFilter = Rhino.DocObjects.ObjectType.Curve
    
    while True:
        go.AddOptionToggle(Opts.names['bAddPts'], Opts.riOpts['bAddPts'])
        go.AddOptionToggle(Opts.names['bNormalize'], Opts.riOpts['bNormalize'])
        go.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        go.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])
        
        res = go.Get()
        
        if res == ri.GetResult.Cancel:
            go.Dispose()
            return
        elif res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            return tuple([objref] + [Opts.values[key] for key in sOpts])
        
        Opts.setValues()
        Opts.saveSticky()
        go.ClearCommandOptions()


def getInput_Parameter(rgCrv):
    """Get parameter with optional input."""
    
    gp = ri.Custom.GetPoint()
    gp.SetCommandPrompt("Pick point or enter parameter")
    
    gp.Constrain(curve=rgCrv, allowPickingPointOffObject=False)
    
    gp.AcceptNumber(True, acceptZero=True)
    
    while True:
        gp.AddOptionToggle(Opts.names['bAddPts'], Opts.riOpts['bAddPts'])
        gp.AddOptionToggle(Opts.names['bNormalize'], Opts.riOpts['bNormalize'])
        gp.AddOptionToggle(Opts.names['bEcho'], Opts.riOpts['bEcho'])
        gp.AddOptionToggle(Opts.names['bDebug'], Opts.riOpts['bDebug'])

        res = gp.Get()
        if res == ri.GetResult.Cancel:
            gp.Dispose()
            return
        elif res == ri.GetResult.Number:
            t = gp.Number()
            if Opts.values['bNormalize']:
                t = rgCrv.Domain.ParameterAt(t)
            return tuple([t] + [Opts.values[key] for key in sOpts])
        elif res == ri.GetResult.Point:
            # Use pick point for radius.
            crv_, t = gp.PointOnCurve()
            crv_.Dispose()
            gp.Dispose()
            return tuple([t] + [Opts.values[key] for key in sOpts])
        Opts.setValues()
        Opts.saveSticky()
        gp.ClearCommandOptions()


def _formatParam(t, fPrecision=16):
    if t is None:
        return "(None)"

    if t == Rhino.RhinoMath.UnsetValue:
        return "(Infinite)"

    if t == 0.0:
        return "0.0"

    if t < 10.0**(-(fPrecision-1)):
        # For example, if t is 1e-5 and fPrecision == 5,
        # the end of this script would display only one digit.
        # Instead, this return displays 2 digits.
        return "{:.2e}".format(t)

    return "{:.{}G}".format(t, fPrecision)


def main():

    bAddPts = None
    bNormalize = None
    bEcho = None
    bDebug = None
    
    rc = getInput_Curve()
    if rc is None: return
    objref = rc[0]
    for key, value in zip(sOpts, rc[1:]):
        exec("{} = {}".format(key, value))
    
    rgCrv = objref.Curve()
    
    s  = "In domain [{},{}] with {} spans:".format(
            _formatParam(rgCrv.Domain.T0),
            _formatParam(rgCrv.Domain.T1),
            rgCrv.SpanCount)
    print s
    
    while True:
        rc = getInput_Parameter(rgCrv)
        if rc is None: return
        t = rc[0]
        for key, value in zip(sOpts, rc[1:]):
            exec("{} = {}".format(key, value))
        
        s  = "  Parameter = {}  Normalized parameter = {}".format(
                _formatParam(t),
                _formatParam(rgCrv.Domain.NormalizedParameterAt(t)))
        
        for iSpan in xrange(rgCrv.SpanCount):
            if (
                    rgCrv.SpanDomain(iSpan).T0 <=
                    t <=
                    rgCrv.SpanDomain(iSpan).T1
            ):
                s += "  SpanIndex={}".format(iSpan)
        
        print s
        
        if bAddPts:
            pt = rgCrv.PointAt(t) if bNormalize else rgCrv.PointAt(t)
            sc.doc.Objects.AddPoint(pt)
            sc.doc.Views.Redraw()


if __name__ == '__main__': main()