%ecoll_2dfun0.m finds the zeros of the momentum, energy equations
%in the the 2d elastic collision in order to find the unknowns.
%Based on the use of the fsolve.m function in the optimization toolbox
%
%function froot = ecoll_2dfun0(X,v1i,v2i,m1,m2,tp10,tp20,tp1)%use this line without "global" stement
function froot =ecoll_2dfun0(X)%use this line with a "global" statement
global v1i v2i m1 m2 tp10 tp20 tp1
%X(1)=v1f, X(2)=v2f, X(3)=tp2
%energy equation: m1*v1i^2+m2*v2i^2-m1*v1f^2-m2*v2f^2=0
%x-momentum: m1*v1i*cos(tp10)+m2*v2i*cos(tp20)-m1*v1f*cos(tp1)-m2*v2f*cos(tp2)=0
%y-momentum: -m1*v1i*sin(tp10)+m2*v2i*sin(tp20)-m1*v1f*sin(tp1)+m2*v2f*sin(tp2)=0
froot = [m1.*v1i.^2+m2.*v2i.^2-m1.*X(1).^2-m2.*X(2).^2;...
         m1.*v1i.*cos(tp10)+m2.*v2i.*cos(tp20)-m1.*X(1).*cos(tp1)-m2.*X(2).*cos(X(3));...
        -m1.*v1i.*sin(tp10)+m2.*v2i.*sin(tp20)-m1.*X(1).*sin(tp1)+m2.*X(2).*sin(X(3));];
