function [Dtot,D] = PDF_distance_Werner_2D(X_bins,Y_bins,F1,F2,N_norm)

if nargin<5
    N_norm = 1e6;
end

%Some checks on the input data
% Overall frequencies should be the same for each PDF
% But need the values in actual counts rather than normalised to less than
% one since are rolling out each value F(i,j) times
S1 = sum(F1(:));
S2 = sum(F2(:));
if S1~=S2 | (S1+S2)<2.01   %the latter is to avoid keeping PDFs that are normalised to one - we will want to normalise these to N_norm
    %but it may be the case that we are comparing low count PDFs, so will
    %just do this if <2.01
    
    %    error('The total frequencies of the two PDFs must be the same! I.e. sum(F1(:)) == sum(F2(:)).');
    
    %If the overall counts are not the same then normalise to a set number
    % = N_norm. For accuracy this should be quite high since need to round
    % the counts to an integer.
    % 
    %Normalise frequencies so that they sum to N_norm and are integers
    F1 = round( N_norm * F1 / S1 );
    F2 = round( N_norm * F2 / S2 );
end



%Normalise the X and Y bins to lie between zero and one - otherwise the X
%and Y distances will be weighted differently since they are different
%physical quantities.
%To do this will just use the range provided for now, but perhaps shoudl
%calculate e.g. the limit between the 1% and 99% percentiles or something?
dX = X_bins(end) - X_bins(1);
X_bins_norm = (X_bins - X_bins(1)) / dX;

dY = Y_bins(end) - Y_bins(1);
Y_bins_norm = (Y_bins - Y_bins(1)) / dY;

%Roll out the PDF to give (paired up) vectors for X and Y for both PDFs
[X1,Y1] = expand_PDF_2D(X_bins_norm,Y_bins_norm,F1);
[X2,Y2] = expand_PDF_2D(X_bins_norm,Y_bins_norm,F2);

%Do the calculation
D = abs( X1 - X2 ) + abs( Y1 - Y2 );
Dtot = mean(D,1);

