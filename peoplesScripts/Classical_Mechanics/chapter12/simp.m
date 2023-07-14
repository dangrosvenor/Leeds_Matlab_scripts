%simp.m
function simpu=simp(f,inc)   
%Simpson's rule for numerical integration
%f is an odd array of evaluated functions in steps inc
ip=length(f);       %must be an odd number
s1=sum(f(2:2:ip-1));%sums all even terms
s2=sum(f(3:2:ip-2));%sums all odd term does not include f(1) and f(ip)
simpu=(4.*s1+2.*s2+f(1)+f(ip))*inc/3.0;%finally add f(1) and f(ip)