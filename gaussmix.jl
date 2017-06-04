
using Distributions
using StatsBase
using PyPlot
using HDF5
using MAT
using JLD
include("fgaussmix.jl")

tau = 9;
N = 100; c=0.67; K=100;
sigma = 1/(c.*N);
ztr = sample(1:N,WeightVec((1/N).*ones(N)),N);
nu = (1/N).*ones(N);

mutr = zeros(N);
mutr = (1/N).*[1:N];

y = zeros(N);
[y[j] = rand(Normal(mutr[ztr[j]],sigma)) for j in 1:N];


nmc = 10000;
z = deepcopy(ztr);
sy = zeros(K); sy2 = zeros(K); n = zeros(K);
[n[j] = sum(1*(z.==j)) for j=1:K];
[sy[j] = sum(y[z.==j]) for j=1:K];
[sy2[j] = sum(y[z.==j].^2) for j=1:K];
mlik = zeros(K);
H = zeros(nmc);

Atr = zeros(N,N);
A = zeros(N,N);
[Atr[j1,j2] = 1*(z[j1]==z[j2]) for j1=1:N, j2=1:N];
Delt = zeros(nmc);

for t=1:nmc
  if mod(t,100)==0
    #print(t);
    #plot(1:(t-1),H[1:(t-1)]);
    print(t);
    print("\n");
    print("Hamming: ")
    print(mean(H[1:(t-1)]));
    print("\n")
    print("Class movement:")
    print(mean(Delt[1:(t-1)]));
    print("\n")
  end
  zlast = deepcopy(z);
  for i=1:N
      sy,sy2,n = remove_y(z[i],y[i],sy,sy2,n)
      for k=1:K
          sy, sy2, n = add_y(k,y[i],sy,sy2,n)
          mlk = lmlik(n,K,sy,sy2,tau,sigma,nu)
          mlik[k] = mlk;
          sy, sy2, n = remove_y(k,y[i],sy,sy2,n)
      end
      mlik = exp(mlik-maximum(mlik))
      wts = mlik./sum(mlik);
      zt = sample(1:K,WeightVec(wts));
      #print(zt)
      z[i] = zt;
      sy,sy2,n = add_y(z[i],y[i],sy,sy2,n)
  end
  [A[j1,j2] = 1*(z[j1]==z[j2]) for j1=1:N, j2=1:N];
  Ad = 1.*abs(A-Atr);
  H[t] = sum(Ad[:]);
  Delt[t] = sum(1.*(zlast.!=z));

end

save(string("/Users/jamesjohndrow/Dropbox (Personal)/Projects/Microclustering/Code/Outputs/bayesout",replace("$c","\.","_"),".jld"),"H",H,"Delt",Delt)
file = matopen(string("/Users/jamesjohndrow/Dropbox (Personal)/Projects/Microclustering/Code/Outputs/bayesout",replace("$c","\.","_"),".mat"), "w");
write(file, "H", H); write(file,"Delt",Delt);
close(file)
