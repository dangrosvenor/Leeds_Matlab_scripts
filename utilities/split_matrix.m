function mat_new=split_matrix(mat,N)
%takes a matrix mat and "splits" each element into N elements of the same
%value. I.e. makes a higher resolution version but doesn't interpolate the
%values - just takes the value for the cell within which the new point
%falls

siz=size(mat);

x=1:siz(1);
x2=repmat(x,[N 1]);
x3=x2(:);

y=1:siz(2);
y2=repmat(y,[N 1]);
y3=y2(:);

[xi,yi]=meshgrid(x3,y3); %creates a matrix that references the original matrix
%N times in the correct positions

inds=sub2ind(size(mat),xi(:),yi(:));
mat_new = mat(inds);
mat_new = reshape(mat_new,size(xi));
mat_new = mat_new';

