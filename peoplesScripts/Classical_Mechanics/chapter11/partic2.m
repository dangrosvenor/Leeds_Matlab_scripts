%partic2.m - given the mass and initial positions, velocities, and
%accelerations of 3 particles, this script calculates their linear
%and angular momenta, energies, forces, and torques 
%help symbolic             %type this if need help on symbolic functions
clear;
format compact;            %Suppress extra line-feeds on outputs
u=[1,1,1];                 %vector used for component sum
m=[1,2,3];                 %mass vector
M=sum(m);                  %total mass
r =[5,4,0;-1,2,0;1,-3,0];  %position vectors
v =[4,0,0;-1,2,0;0,-6,0];  %velocity vectors
a =[4,0,0;2,0,0;0,-6,0];   %acceleration vectors
rcm =[1.0,-0.1667,0];      %center of mass position vector
vcm =[0.3333,-2.3333,0];   %center of mass velocity vector
acm =[1.3333,-3.0000,0];   %center of mass acceleration vector          
for i=1:3
    p(i,:)=m(i)*v(i,:);                %momentum calculation
    l(i,:)=cross(r(i,:),p(i,:));       %angular momentum calculation
    e(i)=dot(p(i,:),p(i,:))/2/m(i);    %kinetic energies
    f(i,:)=m(i)*a(i,:);                %forces
    tau(i,:)=cross(r(i,:),f(i,:));     %torques
end
p                                      %obtained momenta
P=u*p                                  %net momentum 1st way
P=M*vcm                                %P center of mass way
l                                      %obtained angular momenta
L=u*l                                  %net angular momentum
e                                      %energies
E=sum(e(:))                            %total kinetic energy
f                                      %obtained forces
F=u*f                                  %net force 1st way
F=M*acm                                %F center of mass way
tau                                    %obtained torques
Tau=u*tau                              %Net torque