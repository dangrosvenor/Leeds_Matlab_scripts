function [X_new,Y_new,Z_new]=bin_data3D_remove_NaNs(X,Y,Z)
%function [X2,Y2,Z2]=bin_data3D_remove_NaNs(X,Y,Z)
%removes NaN values from corresponding vectors X, Y and Z

a=isnan(Y);
b=find(a==0);

Y2=Y(b);
X2=X(b);  %ignore all the NaN data in Y
Z2=Z(b);

a=isnan(X2);
b=find(a==0);

Y3=Y2(b);
X3=X2(b);  %ignore all the NaN data in X
Z3=Z2(b);

a=isnan(Z3);
b=find(a==0);

Y_new=Y3(b);
X_new=X3(b);  %ignore all the NaN data in X
Z_new=Z3(b);

