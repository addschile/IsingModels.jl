abstract type IsingModel end

mutable struct IsingModel2D <: IsingModel
  h::Float64
  J::Float64
  T::Float64
  n::Int64
  m::Int64
  data::AbstractArray
  rng::AbstractRNG
  function IsingModel2D(h::Float64,J::Float64,T::Float64,n::Int64,m::Int64;seed::UInt64=rand(UInt64))
    data = zeros(Int64,n,m)
    rng = MersenneTwister(seed)
    new(h,J,T,n,m,data,rng)
  end
end
IsingModel2D(h::Float64,J::Float64,T::Float64,n::Int64;seed::UInt64=rand(UInt64)) = IsingModel2D(h,J,T,n,n;seed=seed)
IsingModel2D(n::Int64,m::Int64;seed::UInt64=rand(UInt64)) = IsingModel2D(1.0,1.0,1.0,n,m;seed=seed)
IsingModel2D(n::Int64;seed::UInt64=rand(UInt64)) = IsingModel2D(1.0,1.0,1.0,n,n;seed=seed)

function randomdata!(model::IsingModel)
  model.data = rand(Int64[-1,1],(model.n,model.m))
end

function allup!(model::IsingModel)
  model.data = ones(Int64,(model.n,model.m))
end

function alldown!(model::IsingModel)
  model.data = -ones(Int64,(model.n,model.m))
end

function flipbit(ind::Int64,model::IsingModel)
  return model.data[ind] == 1 ? -1 : 1
end

function flipbit!(ind::Int64,model::IsingModel)
  model.data[ind] == 1 ? model.data[ind] = -1 : model.data[ind] = 1
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
  return -(model.h*fbit + model.J*fbit*nn)
end

#function calc_energy(model::IsingModel)
#  en::Float64 = model.h*sum(mdoel.data)
#end

function calc_magnetization(model::IsingModel)
  return sum(mdoel.data)
end

function calc_avg_magnetization(model::IsingModel)
  return sum(model.data)/length(model.data)
end
