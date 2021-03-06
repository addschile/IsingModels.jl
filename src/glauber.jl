function glaubersweeps(nsweeps::Int64,model::IsingModel,cbevery::Int64=1,cb::Function=defaultcb)
  ndata::Int64 = length(model.data)
  for i in 1:nsweeps # loop over sweeps
    if i%cbevery == 0 && !(cb == nothing)
      cb(model)
    end
    randinds = rand(model.rng,1:ndata,ndata)
    for ind in randinds # loop entries in the data
      ediff::Float64 = localediff(ind,model)
      if rand(model.rng) < 1/(1+exp(ediff/model.T))
        flipbit!(ind,model)
      end
    end
  end
end
