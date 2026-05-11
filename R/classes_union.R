# class unions -----------------------------------------------------------------

# define class unions where methods can be safely shared across

## single and meta-cells -------------------------------------------------------

# Tag union between SingleCells and MetaCells where stuff can be shared
ScOrMc <- S7::new_union(SingleCells, MetaCells)
