%partic1.m - given the mass and initial time-dependent positions of 3 particles,
%this script calculates the velocities, accelerations and their center of mass 
%values. Type "help symbolic" within the command line 
%if need help on symbolic functions
clear;
format compact;                        %Suppress extra line-feeds on outputs
syms t;                                %declare t as a symbol
t = sym('t','real');                   %let t, and the masses be real
m=[1,2,3]                              %mass vector
M=sum(m)                               %total mass
r=[3+2*t^2,4,0;-2+1/t,2*t,0;1,-3*t^2,0]%given position vectors 
v=diff(r,t,1)                          %velocity vector
a=diff(r,t,2)                          %acceleration vector
rcm=m*r/M                              %center of mass coordinate
vcm=m*v/M                              %center of mass velocity
acm=m*a/M                              %center of mass acceleration
%================== evaluate =====================
t=1;                                   %set the time at which to evaluate expressions
r=eval(r)
v=eval(v)
a=eval(a)
rcm=eval(rcm)
vcm=eval(vcm)
acm=eval(acm)
