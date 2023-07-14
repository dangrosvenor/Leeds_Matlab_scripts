%theta_max.m - plots the maximum scattering angle theta_1 versus the
%m2/m1 ratio
clear;
x=[0.0:0.0125:1];
y=acos(sqrt(1-x.^2))/pi;
plot(x,y,'b.')
xlabel('m_2/m_1','FontSize',13)
ylabel('\theta_{1max} (\pi)','FontSize',13)
str=cat(2,'case of m_2 < m_1: Plot of \theta_{1max} versus m_2/m_1 ratio');
title(str,'FontSize',13)
