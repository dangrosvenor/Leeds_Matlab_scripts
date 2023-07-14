%inert_el2.m
function inerel2=inert_el2(m,x,y,z)    
%Inertia function integrands in cartesian coords
%if m=7 it does the volume integrand
if     m==1 inerel2=y^2+z^2;
elseif m==2 inerel2=-x*y;
elseif m==3 inerel2=-x*z;
elseif m==4 inerel2=x^2+z^2;
elseif m==5 inerel2=-y*z;
elseif m==6 inerel2=x^2+y^2;
elseif m==7 inerel2=1.0;
else
    disp ' only 7 integrands are needed '
    return
end