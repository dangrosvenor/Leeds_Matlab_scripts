%calculates the confluent hypergeometric function of the first order

function [sum,n]=conhypgeo1(a,b,z,tol);

sum=1;
X=1;
n=1;

while ( abs(X)>tol | n>1000 )

    n=n+1;
    X=X*(a+n-2)/(b+n-2)*z*factorial(n-2)/factorial(n-1); %facorial(0) = 1
    sum=sum + X;
    
end
    
    