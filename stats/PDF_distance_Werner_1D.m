function [Dtot,D,N_norm] = PDF_distance_Werner_2D(X_bins,F1_in,F2_in,N_norm_in,method)
%F2_in can be consist of multiple vectors of data to be compared
%to F1_in. Needs to be ordered as [nvectors,nbins]

if nargin<4
    N_norm_in = 1e6;   
end

if ~exist('method')
    method = 'multi-dim';
end


N_norm = 1; %Value in case no normalisation is needed.
        


%Some checks on the input data
% Overall frequencies should be the same for each PDF 
% But need the values in actual counts rather than normalised to less than
% one since are rolling out each value F(i,j) times
S1 = sum(F1_in(:));
S2 = sum(F2_in(:));
if S1~=S2 | (S1+S2)<2.01   %the latter is to avoid keeping PDFs that are normalised to one - we will want to normalise these to N_norm
    %but it may be the case that we are comparing low count PDFs, so will
    %just do this if <2.01

%     error('The total frequencies of the two PDFs must be the same! I.e. sum(F1(:)) == sum(F2(:)). Could normalise them.');
% end

%If the overall counts are not the same then normalise to a set number
    % = N_norm. For accuracy this should be quite high since need to round
    % the counts to an integer.
    % 
    %Normalise frequencies so that they sum to N_norm and are integers
    %However, if N_norm is too low then the total number of elements won't
    %be the same, so might need to increase.
    sum_F1=0; sum_F2=-1; %dummy values that just need to be different
    N_norm = N_norm_in;
    while sum_F1~=sum_F2

        F1 = round( N_norm * F1_in / S1 );
        F2 = round( N_norm * F2_in / S2 );
        sum_F1 = sum(F1);
        sum_F2 = sum(F2);
        N_norm=N_norm+N_norm_in; %Increase N_norm until we have the same number of points
        
    end
    
    N_norm = sum_F1; %Since N_norm was incremented at the end of the loop and there may be rounding differences

end





switch method

    case 'cumsum'
        %This method is quicker since it avoids the roll-out - gives a
        %slightly different answer to the multi-dim, though...
        %Perhaps to do with the Xbin normalization?

        X1 = cumsum(F1/N_norm,2); %Will ingnore the initial zero since this always matches between PDFs.
        X2 = cumsum(F2/N_norm,2);
        X1 = X1'; X2 = X2';

    case 'multi-dim'
        %Normalise the X and Y bins to lie between zero and one - otherwise the X
        %and Y distances will be weighted differently since they are different
        %physical quantities.
        %To do this will just use the range provided for now, but perhaps shoudl
        %calculate e.g. the limit between the 1% and 99% percentiles or something?
        dX = X_bins(end) - X_bins(1);
        X_bins_norm = (X_bins - X_bins(1)) / dX;

        %Roll out the PDF to give (paired up) vectors for X and Y for both PDFs
        [X1] = expand_PDF(X_bins_norm,F1);
        [X2] = expand_PDF(X_bins_norm,F2);

end

%Do the calculation
D = abs( X1 - X2 );
Dtot = mean(D,1);



