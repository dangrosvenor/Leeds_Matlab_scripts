function Z = PDFdistZ(x,f,p)
% function Z = PDFdistZ(x,f,p)=
% returns the value x=Z value for which integral( P(x) )[0 Z] = p
% IMPORTANT - x are the bin EDGES! Requried for a correct cumulative PDF
% p runs from 0 to 1
% f can be normalised or not



if size(x,2)~=1
    x=x';
end

if size(f,2)~=1
    y=y';
end

%cumulative sum of frequencies
cum = cumsum(f)./sum(f);
%In f we have the freqencies for the bins given by bin edges in x. So have
%N-1 f values for N x values
cum = cat(1,0,cum); %add a zero to the start since the first value in the cumsum would be the total for
%bin 1. But at xc=0 we want P=0 for the cumulative PDF P(xc)

%N=100;
%Xbins=make_PDF_bins(cum,N);
%X = 0.5 * (Xbins(1:end-1) + Xbins(2:end));

%qh = ndHistc_run([cum], Xbins);

%Z = interp1(cum,x,p); %can't do this since the cum values are not unique

%Use this c script that I wrote
Z = lininterp1f_multidim(cum,x,p,9e99);


