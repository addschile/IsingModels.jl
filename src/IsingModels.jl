module IsingModels

using Random

include("ising.jl")
include("callbacks.jl")
include("mc.jl")
include("glauber.jl")
include("observables.jl")

export IsingModel2D,randomdata!,allup!,alldown!
export mcsweeps,mcsweep
export glaubersweeps
export calc_magnetization,calc_avg_magnetization

end
