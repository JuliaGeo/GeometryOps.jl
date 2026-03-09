module GeometryOpsTestHelpersLibGEOSExt

using GeometryOpsTestHelpers
using LibGEOS

function __init__()
    # Register LibGEOS in the test modules list
    push!(GeometryOpsTestHelpers.TEST_MODULES, LibGEOS)
end

end # module GeometryOpsTestHelpersLibGEOSExt
