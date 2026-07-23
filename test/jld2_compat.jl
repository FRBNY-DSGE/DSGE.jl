# Shared JLD2-compatibility shims for the test suite.
#
# Committed reference fixtures (e.g. regime_switch_data.jld2) were written under older
# versions of DataFrames/Distributions. When the installed package's internal struct
# layout differs, JLD2 can no longer map the stored bytes onto the real type and instead
# hands back an opaque `JLD2.ReconstructedMutable{...}`, which then fails wherever a
# concrete type is expected (e.g. a `df::DataFrame` keyword).
#
# Wrapping a `load(...)` in `as_dataframe(...)` makes the read self-healing: native loads
# pass straight through, and a reconstructed frame is rebuilt from its stored fields. This
# absorbs future layout drift with no need to re-save the fixture.

using DataFrames

"""
    as_dataframe(F) -> DataFrame

Return `F` unchanged if it is already a `DataFrame`; otherwise rebuild one from a
JLD2-reconstructed DataFrame. A DataFrame serializes as fields `(:columns, :colindex)`,
where `colindex` is a `DataFrames.Index` whose `:names` field holds the column symbols.
"""
function as_dataframe(F)
    F isa DataFrame && return F
    cols = collect(getproperty(F, :columns))
    idx  = getproperty(F, :colindex)
    nms  = collect(getproperty(idx, :names))
    return DataFrame(cols, nms; copycols = false)
end
