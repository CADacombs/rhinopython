"""
190410: Created from another script.
190510: Fixed bug in testing for invalid knot vector.

TODO: Test for invalid knot vector may need more work.
"""

import Rhino


def isUniform(knots):
    
    knots = list(knots)
    
    len_knots = len(knots)
    
    if len_knots == 2:
        return True
    elif len_knots < 2:
        return
    
    bPeriodicOrDeg1 = not Rhino.RhinoMath.EpsilonEquals(
            knots[0], knots[1], Rhino.RhinoMath.ZeroTolerance)
    
    if not bPeriodicOrDeg1 and len_knots > 2:
        if not Rhino.RhinoMath.EpsilonEquals(
                knots[-1], knots[-2], Rhino.RhinoMath.ZeroTolerance
        ):
            print "Knot vector is invalid: {}".format(knots)
            return
    
    knotVect1stSpan = knots[1] - knots[0]
    
    if not bPeriodicOrDeg1:
        for i in xrange(1, len_knots//2):
            if not Rhino.RhinoMath.EpsilonEquals(knots[i] - knots[i+1],
                    knotVect1stSpan,
                    Rhino.RhinoMath.ZeroTolerance):
                break
        else:
            print "Degree cannot be determined."
            return
        
        degree = i + 1
    else:
        degree = 1 # Even if periodic of degree > 1.
    
    # Set up for loop's range's start and stop based on whether knot vector is periodic.
    # iKnot_Start: Start of spans to check.
    # iKnot_End: End of spans to check.
    if bPeriodicOrDeg1:
        iKnot_Start = 0
        iKnot_End = len_knots - 1
        knotVect1stSpan_Not0 = knotVect1stSpan
    else:
        iKnot_Start = degree - 1
        iKnot_End = len_knots - degree
        if iKnot_Start == iKnot_End: return True # Is this also true when periodic?
        knotVect1stSpan_Not0 = knots[degree] - knots[degree-1]
    pass
    #  A uniform knot vector's interior parameter spans are all equal.
    for i in range(iKnot_Start, iKnot_End):
        span = knots[i+1] - knots[i]
        if not Rhino.RhinoMath.EpsilonEquals(
                span,
                knotVect1stSpan_Not0,
                Rhino.RhinoMath.ZeroTolerance):
            return False
    
    return True
