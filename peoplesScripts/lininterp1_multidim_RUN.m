function YI = lininterp1_multidim_RUN(X,Y,XI,Xdefault,dim)

sX = size(X);
sY = size(Y);

if length(sX)==2 & min(sX)==1 % X is a vector, assume Y is the array
    Y = shiftdim(Y,-1); %make sure the array has at least 3 dimensions
    sY2 = size(Y);

    ndims=[1:length(sY2)];
    %create a vector with all of the dimensions except the required one
    ndims2=ndims;
    ndims2(dim+1)=[]; % +1 because we added an extra dimension

    Y = permute(Y,[dim ndims2]);

end


YI = lininterp1f_multidim(X,Y,XI,9e99);

