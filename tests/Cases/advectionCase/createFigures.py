#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
advectionCasefoam = OpenFOAMReader(FileName='advectionCase.foam')
advectionCasefoam.MeshRegions = ['internalMesh']
advectionCasefoam.CellArrays = ['psi', 'psiLimitedLinear', 'psiWENO']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1063, 651]

# show data in view
advectionCasefoamDisplay = Show(advectionCasefoam, renderView1)
# trace defaults for the display properties.
advectionCasefoamDisplay.Representation = 'Surface'
advectionCasefoamDisplay.ColorArrayName = [None, '']
advectionCasefoamDisplay.OSPRayScaleArray = 'psi'
advectionCasefoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
advectionCasefoamDisplay.SelectOrientationVectors = 'None'
advectionCasefoamDisplay.ScaleFactor = 0.2
advectionCasefoamDisplay.SelectScaleArray = 'None'
advectionCasefoamDisplay.GlyphType = 'Arrow'
advectionCasefoamDisplay.GlyphTableIndexArray = 'None'
advectionCasefoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
advectionCasefoamDisplay.PolarAxes = 'PolarAxesRepresentation'
advectionCasefoamDisplay.ScalarOpacityUnitDistance = 0.08275538451783654
advectionCasefoamDisplay.GaussianRadius = 0.1
advectionCasefoamDisplay.SetScaleArray = ['POINTS', 'psi']
advectionCasefoamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
advectionCasefoamDisplay.OpacityArray = ['POINTS', 'psi']
advectionCasefoamDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(advectionCasefoamDisplay, ('CELLS', 'psiLimitedLinear'))

# rescale color and/or opacity maps used to include current data range
advectionCasefoamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
advectionCasefoamDisplay.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 'psiLimitedLinear'
psiLimitedLinearLUT = GetColorTransferFunction('psiLimitedLinear')
psiLimitedLinearLUT.RGBPoints = [-4.875717451308104e-12, 0.231373, 0.298039, 0.752941, 0.4958424866175115, 0.865003, 0.865003, 0.865003, 0.9916849732398987, 0.705882, 0.0156863, 0.14902]
psiLimitedLinearLUT.ScalarRangeInitialized = 1.0

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# hide color bar/color legend
advectionCasefoamDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 0.0, 4.568607944718666]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.05000000074505806]
renderView1.CameraParallelScale = 1.5

# save screenshot
SaveScreenshot('Figures/limitedLinear.png', renderView1, ImageResolution=[1063, 651],
    TransparentBackground=1)

# set scalar coloring
ColorBy(advectionCasefoamDisplay, ('CELLS', 'psiWENO'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(psiLimitedLinearLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
advectionCasefoamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
advectionCasefoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'psiWENO'
psiWENOLUT = GetColorTransferFunction('psiWENO')
psiWENOLUT.RGBPoints = [-0.00010480391210876405, 0.231373, 0.298039, 0.752941, 0.5046661401429446, 0.865003, 0.865003, 0.865003, 1.009437084197998, 0.705882, 0.0156863, 0.14902]
psiWENOLUT.ScalarRangeInitialized = 1.0

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 0.0, 4.568607944718666]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.05000000074505806]
renderView1.CameraParallelScale = 1.5

advectionCasefoamDisplay.SetScalarBarVisibility(renderView1, False)
# save screenshot
SaveScreenshot('Figures/WENO.png', renderView1, ImageResolution=[1063, 651],
    TransparentBackground=1)

# set scalar coloring
ColorBy(advectionCasefoamDisplay, ('CELLS', 'psi'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(psiWENOLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
advectionCasefoamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
advectionCasefoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'psi'
psiLUT = GetColorTransferFunction('psi')
psiLUT.ScalarRangeInitialized = 1.0

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 0.0, 4.568607944718666]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.05000000074505806]
renderView1.CameraParallelScale = 1.5

advectionCasefoamDisplay.SetScalarBarVisibility(renderView1, False)
# save screenshot
SaveScreenshot('Figures/analyticalResult.png', renderView1, ImageResolution=[1063, 651],
    TransparentBackground=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 0.0, 4.568607944718666]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.05000000074505806]
renderView1.CameraParallelScale = 1.4150971698348158

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
