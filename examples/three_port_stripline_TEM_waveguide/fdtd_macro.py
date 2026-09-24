from paraview.simple import *
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from paraview import servermanager
import numpy as np

paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------
# Known MaterialID scheme
# ----------------------------
NEGATIVE_ID_LABELS = {
    -5: ('IGP'                      , 0.75),
    -4: ('bounding box'             , 0.0),
    -3: ('PML layers'               , 0.15),
    -2: ('ports identifier'         , 0.75),
    -1: ('lumped ports'             , 1.0),
}

NEGATIVE_ID_COLORS = {
    -5: (0.0, 0.0, 1.0),
    -4: (0.5, 0.5, 0.5),
    -3: (1.0, 0.0, 0.0),
    -2: (1.0, 1.0, 1.0),
    -1: (0.0, 0.0, 0.0),
}

POSITIVE_ID_COLORS = [
    (1.0, 1.0, 0.0),
    (1.0, 0.0, 1.0),
    (0.0, 1.0, 1.0),
    (1.0, 0.63, 0.63),
    (0.67, 0.5, 0.33),
    (0.5, 1.0, 0.5),
]

# ----------------------------
# Get active source and view
# ----------------------------
source = GetActiveSource()
renderView1 = GetActiveViewOrCreate('RenderView')

# ----------------------------
# Show original (outline)
# ----------------------------
sourceDisplay = Show(source, renderView1, 'RectilinearGridRepresentation')
sourceDisplay.Representation = 'Outline'

renderView1.ResetCamera(False, 0.9)
renderView1.Update()

# ----------------------------
# Threshold (main display)
# ----------------------------
threshold1 = Threshold(registrationName='Threshold1', Input=source)
threshold1.LowerThreshold = -6.0
threshold1.UpperThreshold = 20.0

Hide(source, renderView1)
SetActiveSource(threshold1)

threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')
threshold1Display.Representation = 'Surface With Edges'

# Force pipeline update BEFORE doing anything else
threshold1.UpdatePipeline()
renderView1.Update()
Render()

# ----------------------------
# Detect present MaterialIDs
# ----------------------------
present_ids = set()

try:
    vtk_data = servermanager.Fetch(threshold1)
    wrapped = dsa.WrapDataObject(vtk_data)

    all_ids = []

    # Handle composite datasets (multiple blocks)
    if hasattr(wrapped, 'GetNumberOfBlocks'):
        for i in range(wrapped.GetNumberOfBlocks()):
            block = wrapped.GetBlock(i)
            if block is None:
                continue
            block_wrapped = dsa.WrapDataObject(block)
            if 'MaterialID' in block_wrapped.CellData.keys():
                ids = np.unique(block_wrapped.CellData['MaterialID'].astype(int))
                all_ids.extend(ids.tolist())
    # Handle flat datasets
    elif 'MaterialID' in wrapped.CellData.keys():
        ids = np.unique(wrapped.CellData['MaterialID'].astype(int))
        all_ids.extend(ids.tolist())
    else:
        print("WARNING: 'MaterialID' not found in CellData")

    present_ids = set(all_ids)

except Exception as e:
    print(f"WARNING: Could not detect IDs: {e}")

print("Detected IDs:", present_ids)

# ----------------------------
# Build categorical colormap
# ----------------------------
annotations = []
indexed_colors = []
indexed_opacities = []

# Negative IDs
for id_val in sorted(NEGATIVE_ID_LABELS.keys()):
    if id_val not in present_ids:
        continue
    label, opacity = NEGATIVE_ID_LABELS[id_val]
    annotations += [str(id_val), label]
    indexed_colors += list(NEGATIVE_ID_COLORS[id_val])
    indexed_opacities.append(opacity)

# Positive IDs
pos_ids = sorted(i for i in present_ids if i >= 0)

for idx, id_val in enumerate(pos_ids):
    color = POSITIVE_ID_COLORS[idx % len(POSITIVE_ID_COLORS)]
    annotations += [str(id_val), f'material {id_val}']
    indexed_colors += list(color)
    indexed_opacities.append(1.0)

print("Annotations:", annotations)

print("IndexedColors length:", len(indexed_colors))
print("IndexedOpacities length:", len(indexed_opacities))
print("Annotations:", annotations)
print("Annotations length:", len(annotations))

# ----------------------------
# Apply colormap (ParaView 6 correct way)
# ----------------------------
ColorBy(threshold1Display, ('CELLS', 'MaterialID'))
threshold1Display.SetScalarColoring('MaterialID', 1)

# Rescale FIRST before touching the LUT
threshold1Display.RescaleTransferFunctionToDataRange(True, False)
renderView1.Update()

# NOW get the LUT and configure categorical mode
materialIDLUT = GetColorTransferFunction('MaterialID')
materialIDOpacityFunc = GetOpacityTransferFunction('MaterialID')

# Disable auto-rescaling so ParaView won't override your settings
materialIDLUT.AutomaticRescaleRangeMode = "Never"

# Enable categorical mode
materialIDLUT.InterpretValuesAsCategories = 1
materialIDLUT.Annotations = annotations
materialIDLUT.IndexedColors = indexed_colors
materialIDLUT.IndexedOpacities = indexed_opacities

# Force ParaView to recognize the changes
materialIDLUT.Modified()
renderView1.Update()
Render()

# Re-apply after render (ParaView sometimes resets on first render)
materialIDLUT.InterpretValuesAsCategories = 1
materialIDLUT.Annotations = annotations
materialIDLUT.IndexedColors = indexed_colors
materialIDLUT.IndexedOpacities = indexed_opacities
materialIDLUT.Modified()

threshold1Display.MapScalars = 1
threshold1Display.LookupTable = materialIDLUT

# ----------------------------
# Axes grid
# ----------------------------
renderView1.AxesGrid.Visibility = 1

#----------------------------
# Finally, at the end and not before, change lower threshold up to -4.0, then updated the scale bar
#----------------------------
threshold1.LowerThreshold = -5.0
threshold1.UpdatePipeline()
renderView1.Update()
Render()

threshold1Display.SetScalarBarVisibility(renderView1, True)
UpdateScalarBars()
Render()

# Final update
renderView1.Update()
RenderAllViews()