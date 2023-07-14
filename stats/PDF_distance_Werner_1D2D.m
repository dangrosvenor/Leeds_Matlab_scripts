function [Dtot,D] = PDF_distance_Werner_1D2D(X_bins,F1,Y_bins,F2)

if nargin==4
    nD_in=2;
elseif nargin==2
    nD_in=1;
end

%Some checks on the input data
% Overall frequencies should be the same for each PDF
S1 = sum(F1(:));
S2 = sum(F2(:));
if S1~=S2
    error('The total frequencies of the two PDFs must be the same! I.e. sum(F1(:)) == sum(F2(:)).');
end

nD = length(size(F1));
if nD>2
    error('Only set up for 1D or 2D PDFs at the moment. Should be extendable to nD - just need to write the roll out routine for nD.');
end

if nD~=nD_in
    error(['Number of input arguments does not match the dimensionality of the PDFs. PDFs are '...
        num2str(nD) 'D, whereas input vars are as for ' num2str(nD_in) 'D']);
end

%Normalise the X and Y bins to lie between zero and one - otherwise the X
%and Y distances will be weighted differently since they are different
%physical quantities.
%To do this will just use the range provided for now, but perhaps shoudl
%calculate e.g. the limit between the 1% and 99% percentiles or something?
dX = X_bins(end) - X_bins(1);
X_bins_norm = (X_bins - X_bins(1)) / dX;

if nD==2
    dY = Y_bins(end) - Y_bins(1);
    Y_bins_norm = (Y_bins - Y_bins(1)) / dY;
end

if nD==1
    [X1] = expand_PDF(X_bins_norm,F1);
    [X2] = expand_PDF(X_bins_norm,F2);

elseif nD==2
    %Roll out the PDF to give (paired up) vectors for X and Y for both PDFs
    [X1,Y1] = expand_PDF_2D(X_bins_norm,Y_bins_norm,F1);
    [X2,Y2] = expand_PDF_2D(X_bins_norm,Y_bins_norm,F2);
end

D = NaN*ones([length(X1) 1]);
for i=1:length(X1)
    D(i) = abs( X1(i) - X2(i) ) + abs( Y1(i) - Y2(i) );
end

Dtot = sum(D,1);

