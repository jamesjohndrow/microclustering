function lmlik(n,K,sy,sy2,tau,sigma,nu)
  lpr = sum(n.*log(nu)); tau2 = tau^2; sigma2 = sigma^2;
  ll = sum(log(sigma)-n.*log(sqrt(2*pi)*sigma)-.5.*log(n.*tau2+sigma2)-sy2./(2*sigma2)+(tau2*sy.^2)./(2*sigma2*(n.*tau2+sigma2)));
  lml = ll+lpr;
  return(lml)
end

function remove_y(k,yi,sy,sy2,n)
  sy[k] = sy[k]-yi; sy2[k] = sy2[k]-yi.^2; n[k] = n[k] - 1;
  return(sy,sy2,n)
end

function add_y(k,yi,sy,sy2,n)
  sy[k] = sy[k]+yi; sy2[k] = sy2[k]+yi.^2; n[k] = n[k] + 1;
  return(sy,sy2,n)
end
