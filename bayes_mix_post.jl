using Distributions
using StatsBase
using PyPlot
using HDF5
using MAT
using JLD

nmc = 10000;
cs = [collect(10:-1:2);1.5;1;0.67;0.5];
Hs = zeros(nmc,length(cs));
ctr = 0;
for c=cs
  ctr = ctr+1;
  dat = load(string("Outputs/bayesout",replace(replace("$c","\.","_"),"_0",""),".jld"));
  H = dat["H"];
  Hs[:,ctr] = H;
end
colnames(Hs) = string(cs);
boxplot(Hs); title(L"Adjacency matrix $L_0$ error"); xticks(1:13,cs); xlabel("1/c");
savefig("Figures/bayes_mix_hamming_box.png")
