%inert_el.m
function inerel=inert_el(m,x,y,z)    
%Inertia function integrands in cartesian coords
if m==1 inerel=y^2+z^2;
elseif m==2 inerel=-x*y;
elseif m==3 inerel=-x*z;
elseif m==4 inerel=x^2+z^2;
elseif m==5 inerel=-y*z;
elseif m==6 inerel=x^2+y^2;
else
    disp ' only 6 elements are needed '
    return
end