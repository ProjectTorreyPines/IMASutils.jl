module IMASutils

import StaticArrays: SVector, @SVector

include("integration.jl")
export trapz, cumtrapz

include("optimization.jl")
export mirror_bound

include("contour.jl")
export contour_cache, contour_from_midplane, contour_from_midplane!

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end