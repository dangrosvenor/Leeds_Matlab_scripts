function [prcs_out,imed] = percentiles_from_PDF(x_PDF,PDF,prcs,dim)
%[prcs_out,imed] = percentiles_from_PDF(x_PDF,PDF,prcs,dim)
%calculate the prcs percentile values from the PDF of x_PDF
%PDF can be a matrix of pdf vectors - dim tells it which dimension is the
%PDF vector dimension


%x_PDF2 = 0.5 * (x_PDF(1:end-1) + x_PDF(2:end));


prcs=prcs/100;
imed=find(prcs==0.5);





%Nd_prcs=prctile(Nd_PDF_multi,prcs,6);

ndims=[1:length(size(PDF))];
%create a vector with all of the dimensions except the required one
ndims2=ndims;
ndims2(dim)=[]; 

%put the required dimension first
PDF = permute(PDF,[dim ndims2]);

%Ntots = sum(PDF,dim); %total number of datapoints for normalizing each PDF
Ntots = sum(PDF,1); %total number of datapoints for normalizing each PDF

size_PDF = size(PDF);

if size_PDF(1)>1 %added this bit - needs testing. Trying to prevent repmat when only have a 1D PDF
%    Ntots = repmat(Ntots,[ones([1 length(ndims)-1]) size(PDF,dim)]);
    Ntots = repmat(Ntots,[size(PDF,1) ones([1 length(ndims)-1])]);
end


%csum = cumsum(PDF./Ntots,dim); %cumulative sum of the PDFs, so that 0.2
%corresponds to 20th percentile, etc.

%change the order of the array so that the PDF dimension is first
% - this is required for the interp1 function - perhaps change so that the
% dimension can be specified for the c version
%PDF = permute(cumsum(PDF./Ntots,dim),[dim ndims2]);
PDF = cumsum(PDF./Ntots,1);
%N.B. NaN values will be created where Ntots=0 (no data in the PDF)
%PDF(isnan(PDF))=0;

sPDF2=size(PDF);
z=zeros([1 sPDF2(2:end)]);
%add zeros for the start of the cumsum (interpolate using bin edges)
PDF=cat(1,z,PDF);

%x_PDF3 = repmat(x_PDF2',[1 size_PDF(2:end)]);

disp('Running the interp routine');
%Now need to interpolate along the cumulative PDF for the desired percentiles
%x-vector is is the cumulative PDF
%y-vector the bin values
%prcs_out = interp1(PDF,x_PDF2,prcs); %won't work as x needs to be
%distinct (i.e. the same for all of the PDFs) and also needs to be regularly
%spaced, which it isn't since it's the cumulative distribution

%so, run a C based mex version of this 1D interpolation routine (whihc is also much faster than interp1)
% - also, works on mutli-dimensional PDF arrays
%prcs_out = lininterp1f_multidim(PDF,x_PDF2,prcs,9e99);
%prcs_out = lininterp1f(PDF,x_PDF2,prcs,[]); 

prcs_out = lininterp1f_multidim_run(PDF,x_PDF,prcs,1);

fprintf(1,'\nDone\n');
prcs_out(prcs_out>8.9e99)=NaN;



