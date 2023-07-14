function p=eta_get_p(nc)
%calculate pressure given the eta values and the surface pressure
%function p=eta_get_p(nc)


psfc=nc{'PSFC'}(:);
ptop=nc{'P_TOP'}(:);
eta=nc{'ZNU'}(:);
%eta=nc{'ZNW'}(:);

psfc_rep=permute(repmat(psfc,[1 1 length(eta)]) , [3 1 2]); %make 3D matrix consisting of all the horizontal slices replicated in height
eta_rep=permute(repmat(eta,[size(psfc,1) 1 size(psfc,2)]) , [2 1 3]); %make a 3D matrix of columns of eta

p = ptop + eta_rep.*(psfc_rep-ptop); %calculate pressure given the eta values and the surface pressure