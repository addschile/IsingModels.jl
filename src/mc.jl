function mcsweep(model::IsingModel)
  ndata::Int64 = length(model.data)
  randinds = rand(model.rng,1:ndata,ndata)
  for ind in randinds # loop entries in the data
    ediff::Float64 = localediff(ind,model)
    if ediff < 0
      flipbit!(ind,model)
    else
      if rand(model.rng) < exp(-ediff/model.T)
        flipbit!(ind,model)
      end
    end
  end
end

function mcsweeps(nsweeps::Int64,model::IsingModel;cbevery::Int64=1,cb::Function=defaultcb,cbargs=nothing)
  ndata::Int64 = length(model.data)
  for i in 1:nsweeps # loop over sweeps
    if i%cbevery == 0 && !(cb == nothing)
      cb(model,cbargs)
    end
    mcsweep(model)
  end
end
