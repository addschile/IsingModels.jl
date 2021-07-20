using Random

abstract type IsingModel end

mutable struct IsingModel2D <: IsingModel
  h::Float64
  J::Float64
  n::Int64
  m::Int64
  data::AbstractArray
  function IsingModel2D(h,J,n,m)
    data = zeros(n,m)
    new(h,J,n,m,data)
  end
end

mutable struct SquareIsingModel <: IsingModel
  h::Float64
  J::Float64
  n::Int64
  m::Int64
  data::AbstractArray
  function SquareIsingModel(h,J,n)
    data = zeros(n,n)
    m = n
    new(h,J,n,m,data)
  end
end

function randomdata!(model::IsingModel)
  model.data = bitrand((model.n,model.m))
end

function mcsweeps(nsweeps::Int64)
end

function glauberdynamics(nsteps::Int64)
end

model1 = IsingModel2D(1.0,1.0,4,2)
model2 = SquareIsingModel(1.0,1.0,2)

randomdata!(model1)
randomdata!(model2)

println(model1.data)
println(model2.data)
