%trap.m - uses the native trapz
function trapez=trap(f,inc)   
%trapezoid method for numerical integration
%f is an array of evaluated functions in steps inc
trapez=inc*trapz(f);