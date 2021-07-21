using Random

abstract type IsingModel end

mutable struct IsingModel2D <: IsingModel
  h::Float64
  J::Float64
  beta::Float64
  n::Int64
  m::Int64
  data::AbstractArray
  rng::AbstractRNG
  function IsingModel2D(h::Float64,J::Float64,beta::Float64,n::Int64,m::Int64,seed::UInt64)
    data = zeros(n,m)
    rng = MersenneTwister(seed)
    new(h,J,beta,n,m,data,rng)
  end
end

# unseeded constructors
IsingModel2D(h::Float64,J::Float64,beta::Float64,n::Int64,m::Int64) = IsingModel2D(h,J,beta,n,m,rand(UInt64))
IsingModel2D(h::Float64,J::Float64,beta::Float64,n::Int64) = IsingModel2D(h,J,beta,n,n,rand(UInt64))
IsingModel2D(n::Int64,m::Int64) = IsingModel2D(1.0,1.0,1.0,n,m,rand(UInt64))
IsingModel2D(n::Int64) = IsingModel2D(1.0,1.0,1.0,n,n,rand(UInt64))

## seeded constructors -- don't need full constructor since it's the original one
#IsingModel2D(h::Float64,J::Float64,beta::Float64,n::Int64,seed::Int64) = IsingModel2D(h,J,beta,n,n,seed)
#IsingModel2D(n::Int64,m::Int64,seed::Int64) = IsingModel2D(1.0,1.0,1.0,n,m,seed)
#IsingModel2D(n::Int64,seed::Int64) = IsingModel2D(1.0,1.0,1.0,n,n,seed)

function randomdata!(model::IsingModel)
  model.data = convert(Array{Int64},bitrand((model.n,model.m)))
end

function flipbit(ind::Int64,model::IsingModel)
  return model.data[ind] == 1 ? 0 : 1
end

function flipbit!(ind::Int64,model::IsingModel)
  model.data[ind] == 1 ? model.data[ind] = 0 : model.data[ind] = 1
end

function localediff(ind::Int64,model::IsingModel)
  fbit::Int64 = flipbit(ind,model) - model.data[ind]
  # nearest-neighbor interactions
  x::Int64 = ind%model.m
  x = x != 0 ? x : model.m
  y::Int64 = ceil(Int64,ind/model.n)
  nn::Int64 = 0
  # check x-direction w/ pbc
  nn += (x+1) <= model.n ? model.data[x+1,y] : model.data[1,y]
  nn += (x-1) >= 1 ? model.data[x-1,y] : model.data[model.n,y]
  # check y-direction w/ pbc
  nn += (y+1) <= model.m ? model.data[x,y+1] : model.data[x,1]
  nn += (y-1) >= 1 ? model.data[x,y-1] : model.data[x,model.m]
  return model.h*fbit + model.J*fbit*nn
end

function defaultcb(model::IsingModel)
  return nothing
end

function mcsweeps(nsweeps::Int64,model::IsingModel,cbevery::Int64=1,cb::Function=defaultcb)
  ndata::Int64 = length(model.data)
  for i in 1:nsweeps # loop over sweeps
    if i%cbevery == 0 && !(cb == nothing)
      cb(model)
    end
    randinds = rand(model.rng,1:ndata,ndata)
    for ind in randinds # loop entries in the data
      ediff::Float64 = localediff(ind,model)
      if ediff < 0
        flipbit!(ind,model)
      else
        if rand(model.rng) < exp(-model.beta*ediff)
          flipbit!(ind,model)
        end
      end
    end
  end
end

function glaubersweeps(nsweeps::Int64,model::IsingModel,cbevery::Int64=1,cb::Function=defaultcb)
  ndata::Int64 = length(model.data)
  for i in 1:nsweeps # loop over sweeps
    if i%cbevery == 0 && !(cb == nothing)
      cb(model)
    end
    randinds = rand(model.rng,1:ndata,ndata)
    for ind in randinds # loop entries in the data
      ediff::Float64 = localediff(ind,model)
      if rand(model.rng) < 1/(1+exp(model.beta*ediff))
        flipbit!(ind,model)
      end
    end
  end
end

model1 = IsingModel2D(1.0,1.0,1.0,30)
randomdata!(model1)

mcsweeps(10,model1)
