# class unions -----------------------------------------------------------------

# define class unions where methods can be safely shared across

# Tag union between SingleCells, SingleCellsSubset and MetaCells where stuff can
# be shared
ScOrMc <- S7::new_union(SingleCells, SingleCellsSubset, MetaCells)

# Tag union between SingleCells and SingleCellsSubset
ScOrScSubset <- S7::new_union(SingleCells, SingleCellsSubset)

# Tag union between SpatialSpot and SpatialSpotSubset. The spatial methods
# resolve `exp_id` differently on the two but do the same thing afterwards.
SpOrSpSubset <- S7::new_union(SpatialSpot, SpatialSpotSubset)
