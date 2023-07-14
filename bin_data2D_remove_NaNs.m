function [X2,Y2]=bin_data2D_remove_NaNs(X,Y)
%function [X2,Y2]=bin_data2D_remove_NaNs(X,Y)
%removes NaN values from two corresponding vectors X and Y

a=isnan(Y);
b=find(a==0);

Y3=Y(b);
X3=X(b);  %ignore all the NaN data in Y

a=isnan(X3);
b=find(a==0);

Y2=Y3(b);
X2=X3(b);  %ignore all the NaN data in X

