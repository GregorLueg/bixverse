# class unions -----------------------------------------------------------------

# define class unions where methods can be safely shared across

# Tag union between SingleCells, SingleCellsSubset and MetaCells where stuff can
# be shared
ScOrMc <- S7::new_union(SingleCells, SingleCellsSubset, MetaCells)

# Tag union between SingleCells and SingleCellsSubset
ScOrScSubset <- S7::new_union(SingleCells, SingleCellsSubset)

# Tag union of the objects that carry an ATAC store
AtacHost <- S7::new_union(SingleCellsAtac, SingleCellsMultiModal)
