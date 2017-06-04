function m = nmmoments( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

issmall = N<=10;
m = zeros(length(N),1);

m(N==1) = 1;
for j=2:10
   m(N==j) = exactsum(j);
end

m(~issmall) = gammainc(1,N(~issmall),'upper');



end

function s = exactsum(k)
    s = sum( (1:k)./(factorial(1:k).*factorial(k-(1:k))).*derangements(k-(1:k)));
end

function d = derangements(q)
    d = ones(1,length(q));
    d(q>0) = floor(factorial(q(q>0))./exp(1)+1/2); 
end