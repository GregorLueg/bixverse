# class unions -----------------------------------------------------------------

# define class unions where methods can be safely shared across

# Tag union between SingleCells, SingleCellsSubset and MetaCells where stuff can
# be shared
ScOrMc <- S7::new_union(SingleCells, SingleCellsSubset, MetaCells)

# Tag union between SingleCells and SingleCellsSubset
ScOrScSubset <- S7::new_union(SingleCells, SingleCellsSubset)
