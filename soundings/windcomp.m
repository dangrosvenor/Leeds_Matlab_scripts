function [u,v]=windcomp(speed,dir,dircomp)
% function [u,v]=windcomp(speed,dircomp)
% calculates components of vector of mag speed and direction dir if dircomp is the domain orientation
% v becomes the direction parallel to domain direction, so use v for 2-d LEM input

u=-speed.*sin(dir*pi/180-dircomp*pi/180);
v=-speed.*cos(dir*pi/180-dircomp*pi/180); %take component of speed at dircomp degrees
