"""
"""

from __future__ import absolute_import, print_function, unicode_literals

"""
180919: Created.
200714: Added copy option.
220201: Streamlined the UI.  Now, multiple target matrices are allowed.
"""

import Rhino
import Rhino.DocObjects as rd
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


    key = 'bCopy'; keys.append(key)
    values[key] = False
    riOpts[key] = ri.Custom.OptionToggle(values[key], 'No', 'Yes')
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

        #if key == 'fRadius':
        #    if cls.riOpts[key].CurrentValue < 2.0*sc.doc.ModelAbsoluteTolerance:
        #        cls.riOpts[key].CurrentValue = 0.0

        if key in cls.riOpts:
            cls.values[key] = cls.riOpts[key].CurrentValue
        elif key in cls.listValues:
            cls.values[key] = idxList
        else:
            return
        sc.sticky[cls.stickyKeys[key]] = cls.values[key]


def getInput_ObjsToTransform():
    """
    """
    
    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select objects to transform")

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        go.ClearCommandOptions()

        idxs_Opts.clear()

        addOption('bCopy')
        #addOption('bEcho')
        #addOption('bDebug')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            return objrefs

        # An option was selected.
        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


def getInput_XformMatrices():
    """
    """

    sc.doc.Objects.UnselectAll()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select block instance for start transformation matrix")
    go.SetCommandPromptDefault("Identity")

    go.GeometryFilter = rd.ObjectType.InstanceReference

    go.SubObjectSelect = False

    go.AcceptNothing(True)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        idxs_Opts['Identity'] = go.AddOption('Identity')
        addOption('bCopy')

        res = go.Get()

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            geomA = None
            matrixA = rg.Transform.Identity
            break

        if res == ri.GetResult.Object:
            objref = go.Object(0)
            go.Dispose()
            geomA = objref.Geometry()
            matrixA = geomA.Xform
            break

        # An option was selected.

        if go.Option().Index == idxs_Opts['Identity']:
            go.Dispose()
            geomA = None
            matrixA = rg.Transform.Identity
            break

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


    sc.doc.Objects.UnselectAll()

    go = ri.Custom.GetObject()

    go.SetCommandPrompt("Select block instances for end transformation matrices")
    if geomA: go.SetCommandPromptDefault("Identity")

    go.GeometryFilter = rd.ObjectType.InstanceReference

    go.SubObjectSelect = False

    go.AcceptNothing(geomA)

    idxs_Opts = {}

    def addOption(key): idxs_Opts[key] = Opts.addOption(go, key)

    while True:
        if geomA: idxs_Opts['Identity'] = go.AddOption('Identity')
        addOption('bCopy')

        res = go.GetMultiple(minimumNumber=1, maximumNumber=0)

        if res == ri.GetResult.Cancel:
            go.Dispose()
            return

        if res == ri.GetResult.Nothing:
            go.Dispose()
            matricesB = [rg.Transform.Identity]
            break

        if res == ri.GetResult.Object:
            objrefs = go.Objects()
            go.Dispose()
            matricesB = [objref.Geometry().Xform for objref in objrefs]
            break

        # An option was selected.

        if go.Option().Index == idxs_Opts['Identity']:
            go.Dispose()
            matricesB = [rg.Transform.Identity]
            break

        for key in idxs_Opts:
            if go.Option().Index == idxs_Opts[key]:
                Opts.setValue(key, go.Option().CurrentListOptionIndex)
                break


    return matrixA, matricesB


def main():

    objrefs_ToXform = getInput_ObjsToTransform()
    if objrefs_ToXform is None: return


    sc.doc.Objects.UnselectAll()

    rc = getInput_XformMatrices()
    if rc is None: return

    xform_Start, xforms_Target = rc

    sc.doc.Objects.UnselectAll()

    rc, inverse = xform_Start.TryGetInverse()
    if not rc:
        print("Inverse transformation matrix could not be computed.")
        return

    bCopy = Opts.values['bCopy']


    if len(xforms_Target) == 1:
        xform_Target = xforms_Target[0]
        xform_To_apply = xform_Target * inverse
        for objref in objrefs_ToXform:
            sc.doc.Objects.Transform(
                objref,
                xform=xform_To_apply,
                deleteOriginal=not bCopy)
        sc.doc.Views.Redraw()
        return


    # Multiple target matrices.
    for xform_Target in xforms_Target:
        xform_To_apply = xform_Target * inverse

        for objref in objrefs_ToXform:
            sc.doc.Objects.Transform(
                objref,
                xform=xform_To_apply,
                deleteOriginal=False)

    if not bCopy:
        for objref in objrefs_ToXform:
            sc.doc.Objects.Delete(
                objref=objref,
                quiet=False)


    sc.doc.Views.Redraw()


if __name__ == '__main__': main()