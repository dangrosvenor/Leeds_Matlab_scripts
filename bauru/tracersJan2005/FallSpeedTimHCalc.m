clear twoD ice snow graupel icenc snownc graupelnc

npes=npess2(idir);

ice=icediagsALL(idir).i(:,dumprange,42)/npes; %mean ice MR time/height
snow=icediagsALL(idir).i(:,dumprange,40)/npes; %mean snow MR time/height
graupel=icediagsALL(idir).i(:,dumprange,41)/npes; %mean graupel MR time/height

icenc=icediagsALL(idir).i(:,dumprange,43)/npes; %mean ice NC time/height

% if idir==2
%     %scale number for idir=2 so that the average mass is the same as in idir=1
%     icenc=icediagsALL(2).i(:,dumprange,42)/npes .* icediagsALL(1).i(:,dumprange,43)/npes ./ icediagsALL(1).i(:,dumprange,42)/npes;
% end

snownc=icediagsALL(idir).i(:,dumprange,45)/npes; %mean snow NC time/height
graupelnc=icediagsALL(idir).i(:,dumprange,44)/npes; %mean graupel NC time/height


% 	switch it
% 	case 1 %ice
%         lab='Ice';
%         im=6;
%         in=7;
% 	case 2 %snow
%         lab='Snow';
%         im=4;
%         in=9;
% 	case 3 %graupel
%         lab='Graupel';
%         im=5;
%         in=8;
% 	end
    
twoD.Q(:,:,6)=ice;
twoD.Q(:,:,7)=icenc;
twoD.Q(:,:,4)=snow;
twoD.Q(:,:,9)=snownc;
twoD.Q(:,:,5)=graupel;
twoD.Q(:,:,8)=graupelnc;


clear ice_flux snow_flux graupel_flux tot_fallflux

const=[1 1 1;1 1 2;1 1 3;1 1 4;1 1 5]; %graupel schemes
%const=[1 1 1;1 2 1;1 3 1;1 4 1;1 5 1;1 6 1]; %snow schemes
%const=[1 1 1;2 1 1;3 1 1]; %ice

%for i=1:size(const,1)
    
    for i=1:1
	Vi=FallSpeed(GridDan(idir).RHON,twoD,'ice',const(i,:));
	Vs=FallSpeed(GridDan(idir).RHON,twoD,'snow',const(i,:));
	Vg=FallSpeed(GridDan(idir).RHON,twoD,'graupel',const(i,:));
	
	rho=repmat(GridDan(idir).RHON,[1 size(Vi,2)]);
	
	ice_flux(i).i=rho.*ice.*Vi;
	snow_flux(i).i=rho.*snow.*Vs;
	graupel_flux(i).i=rho.*graupel.*Vg;
	
	tot_fallflux(i).i=ice_flux(i).i+snow_flux(i).i+graupel_flux(i).i;

end



'done fall speed calc'