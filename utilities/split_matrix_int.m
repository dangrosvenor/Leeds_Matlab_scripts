function mat_new=split_matrix_int(mat,N)
%takes a matrix mat and "splits" each element into N elements
%interpolates in between

siz=size(mat);

x=1:1/(N):siz(1);
y=1:1/(N):siz(2);

[xi,yi]=meshgrid(x,y); %

X=[1:siz(1)];
Y=[1:siz(2)];

mat_new=interp2(Y,X,mat,yi,xi,'linear')';