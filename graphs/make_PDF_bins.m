function Ybins=make_PDF_bins(Y,nYpdf,min_ovr,max_ovr)

%nYpdf=100;

if nargin==4
    minY=min_ovr;
    maxY=max_ovr;
else
    minY=minALL(Y);
    maxY=maxALL(Y);
end

dY=(maxY-minY)/nYpdf;

Ybins=[minY:dY:maxY];