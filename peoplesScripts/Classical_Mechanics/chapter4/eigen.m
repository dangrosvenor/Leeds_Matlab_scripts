%eigen.m
clear;
M=[[1,-1];[-1,1]]  %input matrix M
[P,L]=eig(M)       %eigen vectors, P, and eigenvalues, L
PI=inv(P)          %inverse of P
P*PI               %the unit matrix should result
PI*M*P             %check eigenvectors. Should get eigenvalues back