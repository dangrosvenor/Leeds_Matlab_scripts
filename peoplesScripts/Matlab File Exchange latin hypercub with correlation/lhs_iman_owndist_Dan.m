function z=lhs_iman(corr_in,nsample,x_pdf,f_pdf,ntry)
% z=lhs_iman(xmean,xsd,corr_in,nsample,nloop)
% LHS with correlation, normal distribution 
% method of Iman & Conover
% Iman, R. L., and W. J. Conover. 1982. A Distribution-free Approach to Inducing Rank Correlation 
%      Among Input Variables. Communications in Statistics B 11:311-334
%
% Input:
%   xmean   : mean of data (1,nvar)
%   xsd     : std.dev of data (1,nvar)
%   corr_in    : correlation matrix of the variables (nvar,nvar)
%   nsample : no. of samples
%   ntry    : optional, no of trial to get a close correlation matrix
% Output:
%   z       : random sample (nsample,nvar)
%   Budiman (2004)

nvar=length(x_pdf);

if(nargin==4), ntry=1; end;

% induce data with correlation
P = chol(corr_in);
P=P';

%Note - origininal lhs_iman was passing means of zeros and std of ones here 
% - producing a standard normalised distrubtion
%xm=zeros(1,nvar);
%xs=ones(1,nvar);
%R=latin_hs(xm,xs,nsample,nvar);
R2=latin_hs_owndist_Dan(nsample,nvar,x_pdf,f_pdf);

%Have to convert to a zero based distribution with a std dev of
%one. Seems like the correlation introduction works even if the base distribution R is not gaussian
xmean = mean(R2,1);
xsd = std(R2,1);
    for j=1:nvar 
        R(:,j) = (R2(:,j) - xmean(j))/xsd(:,j);
    end

T = corrcoef(R);
Q=chol(T);
Q=Q';
    
S = P * inv(Q);
RB= R*S';

amin=realmax;
for il=1:ntry
    for j=1:nvar    
        % rank RB
        [r,id]=ranking(RB(:,j));
        % sort R
        [RS,id]=sort(R(:,j));
        % permute RS so has the same rank as RB
        z(:,j) = RS(r).*xsd(j)+xmean(j); 
%        z(:,j) = RS(r) + xmean(j);         
%        z(:,j) = RS(r);
    end
    ae=sum(sum(abs(corrcoef(z)-corr_in)));
    if(ae<amin),
        zb=z;
        amin=ae;
    end;
end

z=zb;