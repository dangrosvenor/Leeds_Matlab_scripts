%det_soln.m - symbolically finds roots of a matrix determinant as well as
%eigenvectors, eigenvalues of the cube's inertia tensor
%type "help symbolic" within MATLAB'S command line for help on symbolic functions
clear; format compact;
syms c i a x y z M r;   %symbolic variables
r=M/a^3;                %cube density
%matrix elements in units of c use 3d symbolic integration
%diagonal elements - think of c=M*a^2 - our inertia unit
Ixx=r*int(int(int((y^2+z^2),x,0,a),y,0,a),z,0,a)*c/(M*a^2);
Iyy=r*int(int(int((x^2+z^2),x,0,a),y,0,a),z,0,a)*c/(M*a^2);
Izz=r*int(int(int((x^2+y^2),x,0,a),y,0,a),z,0,a)*c/(M*a^2);
%products of inertia
Ixy=-r*int(int(int((x*y),x,0,a),y,0,a),z,0,a)*c/(M*a^2);
Ixz=-r*int(int(int((x*z),x,0,a),y,0,a),z,0,a)*c/(M*a^2);
Iyz=-r*int(int(int((y*z),x,0,a),y,0,a),z,0,a)*c/(M*a^2);
Iyx=Ixy; Izx=Ixz; Izy=Iyz;   %use symmetry for the rest
A=[[Ixx,Ixy,Ixz];[Iyx,Iyy,Iyz];[Izx,Izy,Izz]]%inertia tensor
I= i*triu(tril(ones(3)));    %create a unit matrix 3x3
M=A-I;                       %whose determinant we seek
R=simplify(factor(det(M)))   %to see the roots
[P,Q]=eig(A)                 %P=column eigenvectors,Q=eigenvalues
