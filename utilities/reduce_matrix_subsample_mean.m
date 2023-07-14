function [mat_new,N_vals,std1_new,mat2_new,std2_new,cov_new]=reduce_matrix_subsample_mean(mat,N,M,options,mat2,thresh_N)
%mat_new=reduce_matrix_subsample_mean(mat,N,M)
%
%takes a matrix mat and reduces it in size by averaging over a certain
%window of [N M]
%NaNs are ignored for the means for each block.
%N_vals is the number of non-NaN values that went into the average for each
%block.

icalc_cov=0;
ithresh_N=0;

if exist('options')
    if length(findstr(options,'covariance')>0)
        icalc_cov=1;
    end
    if length(findstr(options,'N_threshold')>0)
       ithresh_N=1;
    end

end

siz=size(mat);
psiz=prod(siz);

L1 = siz(1);
L2 = siz(2);

N1=L1/N;
N2=L2/M;

%We need the array to be divisible exactly by N and M, so chop some off if it
%isn't
L1_new = floor(N1)*N;
L2_new = floor(N2)*M;

%Amount to chop off - divide bewteen each side (so /2)
chop_N = (L1 - L1_new)/2;
chop_M = (L2 - L2_new)/2;
%Account for odd numbers
chop_N1 = floor(chop_N);
chop_N2 = ceil(chop_N);
chop_M1 = floor(chop_M);
chop_M2 = ceil(chop_M);

%Chop it off the array
mat = mat(1+chop_N1:end-chop_N2,1+chop_M1:end-chop_M2);
%Re-calculate array properties
L1=L1_new;
L2=L2_new;
N1=L1/N;
N2=L2/M;
siz=size(mat);
psiz=prod(siz);


%


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
%Each block of 1234... needs to be L1 long
j3=[];
for x=1:L2/M
    j=[(x-1)*M+1:(x*M)]';
%    j2=repmat(j,[L1/N 1]);
    j2=repmat(j,[L1 1]);    
    j3=cat(1,j3,j2);
end



%make the linear indices for the array
inds = sub2ind(size(mat),i4(:),j3(:));

%make the new array with the changed order
mat_new=mat(inds);
%reshape so that we have all M*N elements for each block in each row
mat_new = reshape(mat_new,[M*N psiz/(M*N)]);




if icalc_cov==1
    mat2_new=mat2(inds);
    mat2_new = reshape(mat2_new,[M*N psiz/(M*N)]);


    for icov=1:size(mat_new,2)
        A = mat_new(:,icov);
        B = mat2_new(:,icov);
        inanA=isnan(A); inanB=isnan(B); 
        inan=find(inanA==1 | inanB==1);
        A(inan)=[];
        B(inan)=[];
        if length(A)>0
            cov_new2 = cov(A,B);
            cov_new(icov,1) = cov_new2(1,2);
        else
            cov_new(icov,1)=NaN;
        end
    end
        %Do for second array if providing two arrays (e.g. co-variance)
        [mat2_new, N_vals2, std2_new] = meanNoNan(mat2_new,1);
end

%do non-NaN mean and std. dev over each row
[mat_new, N_vals, std1_new] = meanNoNan(mat_new,1);

%reshape into the original size divided by N and M
mat_new = reshape(mat_new,[L1/N L2/M]);
std1_new = reshape(std1_new,[L1/N L2/M]);
N_vals = reshape(N_vals,[L1/N L2/M]);
% 
if icalc_cov==1
    mat2_new = reshape(mat2_new,[L1/N L2/M]);
    std2_new = reshape(std2_new,[L1/N L2/M]);
    cov_new = reshape(cov_new,[L1/N L2/M]);
    N_vals2 = reshape(N_vals2,[L1/N L2/M]);
end

if ithresh_N==1
   inan = find(N_vals<thresh_N);   
   mat_new(inan)=NaN;
   std1_new(inan)=NaN;
   
   if icalc_cov==1
      inan2 = find(N_vals2<thresh_N);
      mat2_new(inan2)=NaN;
      std2_new(inan2)=NaN;
      
      inan12 = find(N_vals<thresh_N | N_vals2<thresh_N);
      cov_new(inan12)=NaN;
   end    
   
end


