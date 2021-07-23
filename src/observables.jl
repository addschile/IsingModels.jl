#function calc_energy(model::IsingModel)
#  en::Float64 = model.h*sum(mdoel.data)
#end

function calc_magnetization(model::IsingModel)
  return sum(mdoel.data)
end

function calc_avg_magnetization(model::IsingModel)
  return sum(model.data)/length(model.data)
end
