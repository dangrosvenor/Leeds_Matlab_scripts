%binary2_der.m: returns the derivatives for binary mass system
function derivs = binary1_der( t, w, flag,c,m1,m2)
% c=4*pi^2, m1=mass of 1st body, m2=mass of 2nd body
% Entries in the vector of dependent variables are:
% w(1,2,...8)=x1,v1x,y1,v1y,x2,v2x,y2,v2y
wr=sqrt((w(5)-w(1)).^2+(w(7)-w(3)).^2);
derivs = [w(2); c*m2*(w(5)-w(1))./wr.^3; w(4); c*m2*(w(7)-w(3))./wr.^3;...
          w(6); -c*m1*(w(5)-w(1))./wr.^3;w(8); -c*m1*(w(7)-w(3))./wr.^3];
