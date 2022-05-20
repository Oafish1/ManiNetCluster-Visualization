# 1.5.12
- Bugfix for plot cluster coloring
- More appropriate cluster bounding given color scheme

# 1.5.11
- Fix button alignment
- Add `pseudocells.R` informational file

# 1.5.10
- Fix `Pseudocell Help` popup not popping up on second+ activation

# 1.5.9
- Add pseudocell popup instructions
- Additional notation on optional, required uploads

# 1.5.8
- Fix low dim visualization plotting bug
- Update download function

# 1.5.7
- `eventReactive` bugfix on `get_data` that caused crashing when no correspondence matrix was uploaded and `use_boma` flag was false
- Fix error on modifying clustering after changing dimension by removing dependence on `input$d`

# 1.5.6
- Small formatting fixes/tweaks

# 1.5.5
- Added more color options
- Fixed small bug for statistics with action buttons
- Styling changes
- UI changes to correspond with figure

# 1.5.4
- Added a run button
- Removed `boma` folder
- Reorganize functions
- Revise UI on uploads

# 1.5.3
- Add color theme selector
- Add higher resolution figure 1
- Add slider for global alignment KNN
- Filter pseudocell datasets
- Rearrange UI to be more intuitive

# 1.5.2
- UI changes
- More user explanations
- Dataset renaming and filtering
- BOMA ordering bugfix
- Added KNN method for local alignment
- Automatically detect genes in common between datasets

# 1.5.1
- UI changes
- Dataset uploads
- Dataset size reductions
- Common gene detection

# 1.5.0
- Add better dataset annotations
- Add figure
- Add more datasets
- Added `preprocess.R` for adding datasets more easily
- Better error handling for statistics
- More error catching for sample datasets
- Various UI changes and fixes

# 1.4.0
- Helper text revisions
- Implement pairwise distance
- Implement quick dataset loading framework
- Legend changes
- Moved distance heatmap
- Properly implement BOMA sorting method
- Properly implement DTW method
- Revised and compressed layout

# 1.3.12
- BOMA rough implementation
- Fix UI bugs

# 1.3.11
- Add BOMA UI framework and reactive elements

# 1.3.10
- Add color bar and unfinished integrated application
- Revise certain UI elements
- Add runtime warnings

# 1.3.9
- Initial direction/scale of plots same
- Switch default to 'time'
- Add 'color' annotation to series selection
- Remove medoids bar
- Add link to download base dataset
- Error catching

# 1.3.8
- Added compatibility for `ShinyApps.io`

# 1.3.7
- Uploaded default data

# 1.3.6
- Added movable plots
- Transitioned main plots to `plotly`

# 1.3.5
- Added heatmap labels
- Various UI updates
- Better accuracy barplot formatting
- Fixed cluster coloration bug
- Added BOMA checkbox, still needs full implementation

# 1.3.4
- Change UI text
- Code optimization
- Fixed legend title bug
- Fix startup errors
- Add indication for alignments in progress
- Move default data read to somewhere more appropriate
- Update helper UI icons
- Directory reorganization
- Uploaded raw svg for icons

# 1.3.3
- Removed upload/download size limitation
- Update `README.md` with minimal usage instructions

# 1.3.2
- Remove `time` default coloring
- Add dynamic selection of coloring column in metadata
- Code cleanup

# 1.3.1
- Remove `None` clustering
- Add www folder
- Add helping images
- Additional metadata error analysis and flexibility
 - Added bar plot for label transfer accuracy

# 1.3.0
- Clustering visualization fix
- Visualization QOL updates
- UI and UX changes (including terminology, coloration, and labeling)
- Added optional clustering visualization for Scatter3D
- Update downloaded information

# 1.2.0
- Add heatmap visualizing cross-modal clusters
- QOL UI improvements
- Spectral cluster framework
- Added placeholder for gene outcomes

# 1.1.1
- 3D graph label bugfix

# 1.1.0
- More QOL notes and default changes
- Enhanced explanations
- Spread out UI
- `Visualization` clustering section added with `K-Means` and `PAM` methods
- Optional `label transfer accuracy` calculation added

# 1.0.0
- Initial Release
- Main UI, dynamic parameter selection, basic upload/download
