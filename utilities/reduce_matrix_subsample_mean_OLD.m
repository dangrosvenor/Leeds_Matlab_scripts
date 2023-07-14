function mat_new=reduce_matrix_subsample_mean(mat,N,M)
%mat_new=reduce_matrix_subsample_mean(mat,N,M)
%takes a matrix mat and reduces it in size by averaging over a certain
%window of [N M]
%Think don't have all the N and M right... so at the moment it only works for N=M... sort this out one day!

siz=size(mat);

psiz=prod(siz);

nblocks = psiz / (N*M);

L1 = siz(1);
L2 = siz(2);

%taking the approach where we create a list of linear indices that
%reference the array in such a way that consequetive M*N elements will 
%contain the data for each block
%Will do this by creating i and j indices for this referencing and then
%doing sub2ind to convert to a big list of linear indices to reorder the
%array

%so for i will need e.g. 1111 2222 3333 .... L1L1L1L1 (an [M*L1 1] matrix) repeated L2/N times 
%Start by creating 1111 2222 3333 .... L1L1L1L1 first. Imagine this
%array re-arranged to be [1111; 2222; ...] i.e. a [N L1] matrix than each
%row is just 1:L1 and repeat this M times
i = [1:L1];
i2= repmat(i,[M 1]);
i3 = i2(:);
%at this stage have e.g. 1111 2222 3333....
%replicate for all the rows
i4 = repmat(i3,[1 L2/M]);

%for j need e.g. [1234 1234 .....;
%                 5678 5678 ....
%make each block of [1234; 5678] first by just re-arraning a linear array
sj = size(i4);
j=[1:sj(2)*M];
j2=reshape(j,[N sj(2)]);
%now just repeat this
j3 = repmat(j2,[L1 1]);

%make the linear indices for the array
inds = sub2ind(size(mat),i4(:),j3(:));

%make the new array with the changed order
mat_new=mat(inds);
%reshape so that we have all M*N elements for each block in each row
mat_new = reshape(mat_new,[M*N psiz/(M*N)]);
%do non-NaN mean over each row
mat_new = meanNoNan(mat_new,1);
%reshape into the original size divided by N and M
mat_new = reshape(mat_new,[L1/N L2/M]);
% 


