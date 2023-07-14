%det_soln2_2d.m - for a 2 dimensional solid
%numerically finds roots of a matrix determinant as well as
%eigenvectors, eigenvalues of the rectangle's inertia tensor
clear; format compact;
a=1;b=1;                          %solid plate sides
M=1; rho=M/(a*b);                 %solid mass, density
ax=0; ay=0;                       %origin at corner as lower limits
bx=a+ax; Nx=11; dx=(bx-ax)/(Nx-1);%x upper limit, x spacing
by=b+ay; Ny=11; dy=(by-ay)/(Ny-1);%y upper limit, y spacing
ie=6;                             %matrix elements needed
x=[ax:dx:bx];y=[ay:dy:by];        %x,y grid
%below, trap is MATLAB's trapezoid rule
%simp is the more accurate simpson rule
    for j=1:Ny
        for i=1:Nx
            for m=1:ie
                fx(m,i)=rho*inert_el(m,x(i),y(j),0);%z=0 in 2d
            end
        end
        for m=1:ie
            %fy(m,j)=trap(fx(m,:),dx);%trapezoidal rule
            fy(m,j)=simp(fx(m,:),dx); %Simpson rule
        end
    end
%finally integrate over the y coord to get the moments
for m=1:ie
    %ff=trap(fy(m,:),dy);%trapezoidal
    ff=simp(fy(m,:),dy); %Simpson
       if     m==1 Ixx=ff;
       elseif m==2 Ixy=ff;
       elseif m==3 Ixz=ff;
       elseif m==4 Iyy=ff;
       elseif m==5 Iyz=ff;
       elseif m==6 Izz=ff;
       end
%fprintf('m= %2i, The integral is %4.3f\n',m,ff)  
end
Iyx=Ixy; Izx=Ixz; Izy=Iyz;   %use symmetry for the rest
A=[[Ixx,Ixy,Ixz];[Iyx,Iyy,Iyz];[Izx,Izy,Izz]];%inertia tensor
disp 'Inertia Tensor'
disp(rats(A))                     %display A in string fraction form
I= i*triu(tril(ones(3)));         %create a unit matrix 3x3
M=A-I;                            %whose determinant we seek
[P,Q]=eig(A);                     %P=column eigenvectors,Q=eigenvalues
ac=round(100*acos(P)*180/pi)/100; %angles related to eigenvectors to 2 dec.
disp 'Direction cosines matrix'
fprintf('cos(%4.2d) cos(%4.2d) cos(%4.2d)\n',ac)%direction cosines
disp 'Principal moments'
S(1:3)=diag(Q);disp(rats(S)) %principal moments in string fraction form
%can check for orthogonaly of P1 and P2
%dot(P(1:3,1),P(1:3,2)),dot(P(1:3,1),P(1:3,3)),dot(P(1:3,2),P(1:3,3))
%can check that the eigenvectors add to 1
%sum(P(1:3,1).^2), sum(P(1:3,2).^2),sum(P(1:3,3).^2)
