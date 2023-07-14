function out=load_am3(var,ind_str,nc,gcm_idays,cosp_flag,iLON)
%function out=load_am3(var,ind_str,nc,gcm_idays)
%Swaps the lon and lat dimensions in AM3 arrays to be ordered lat,lon
%as in CAM5. Also flips the lat dimension to make it the same as CAM5.
%(lat increases with increasing array index)

eval_str_dat = ['nc{''' var '''}' ind_str ';'];
%note the use of gcm_idays is determined in ind_str
eval(eval_str_dat);
dat=ans;

siz=size(dat);
N=length(siz);

if length(gcm_idays)==1 & siz(2)~=1
    dat=shiftdim(dat,-1); %shifts the matrix dimensions to the right giving a singleton matrix
    %at the left - this ensures that the read function works when only have
    %one time dimension. But only when have more than 2 non-singleton dims (otherwise
    %have e.g [40 1] arrays)
end

%size may have changed now
siz=size(dat);
N=length(siz);

%want to switch the lat and lon and flip one of the dimensions
switch N
    case 2
        perm = [2 1];
        iflipdir = 1;
    case 3
        perm = [1 3 2];
        iflipdir = 2;
    case 4
        perm = [1 2 4 3];
        iflipdir = 3;        
end

dat = permute(dat,perm);

out = flipdim(dat,iflipdir);

if nargin>4
    if cosp_flag==1
        out(out<-1e29)=NaN;
    end
end

