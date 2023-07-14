
for isza=1:length(Re_all_pdfs_lowCF)
    Re_vals(isza) = eval(['mean( (Re_all_pdfs_' low_or_highCF '{isza}(idat).dat).^-2.5 );']);
    Tau_vals(isza) = mean( (Tau_all_pdfs{isza}(itau).dat).^0.5 );
    %as would have been used in the actual calc of mean
    %Nd
%    prod_vals(isza) = eval(['mean( Re_all_pdfs_' low_or_highCF '{isza}(idat).dat.^-2.5 .* Tau_all_pdfs{isza}(itau).dat.^0.5 );']);
end

Nall = eval(['Nd_vs_sza_' low_or_highCF '(idat).x'';']);
N1 = Nall(1);
N2 = Nall(end);

%const = N1 / prod_vals(1); %const from N1 = k * mean(Tau^0.5 .* Re.^-2.5)

%scale_fac = Re_vals(1)*Tau_vals(1) / prod_vals(1);

%change scale_fac to be a vector (i.e. different for each sza)
%scale_fac = Re_vals.*Tau_vals ./ prod_vals;

%this is now a multiplication factor to be applied (combines the scale
%factor and the k constant)
scale_fac = Nall ./ (Re_vals.*Tau_vals);
