#=
# DataFrames extension

This module simply extends Core's `reconstruct_table` method to:
- work with DataFrames
- not copy columns unless it is necessary to do so
- allow passing through known kwargs to the constructor


In the future, if we ever end up defining `ApplyToFeatures` on a table,
then we will need to add some form of method for that...
which will likely entail adding an extra positional argument to the
reconstruct_table method, and checking whether `other_column_names`
is equal to `setdiff(Tables.columnnames(input), geometry_column_names)`.

If it is not then we will have to reconstruct the whole DataFrame from the 
GI.Feature named-tuple row table representation.
=#
module GeometryOpsDataFramesExt

import GeometryOpsCore
using DataFrames

GeometryOpsCore.used_reconstruct_table_kwargs(::DataFrames.DataFrame) = (:copycols,)

function GeometryOpsCore.reconstruct_table(
        input::DataFrames.DataFrame, geometry_column_names, geometry_columns, 
        other_column_names, args...; 
        copycols = true, kwargs...
    )
    # Create a new dataframe, let the rest be the same
    new_df = DataFrame(input; copycols)

    # Copy over the geometry columns
    for (colname, col) in zip(geometry_column_names, geometry_columns)
        new_df[!, colname] = col
    end

    # The other columns were already copied over by the constructor,
    # so we're done.
    # Metadata is set in the `_apply_table` method.
    return new_df
end

end