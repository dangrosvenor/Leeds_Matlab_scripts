function [hemm,iwczt]=readEMMf(fname);
% [temm,Temm,arr]=readEMMf(fdir,fname,col);
% col is the column in which the time is stored
iwczt=dlmread([fname],' ');
%hemm=iwczt(