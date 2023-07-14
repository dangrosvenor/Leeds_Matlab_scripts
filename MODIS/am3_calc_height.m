function [h_full,h_half]=am3_calc_height(Pfull,Tfull,Phalf)
%function [h_full,h_half]=am3_calc_height(Pfull,Tfull,Phalf) 
%N.B. using the hydrostatic approach here tends to fail at 
%P lower than ~50hPa. So get NaNs for above here.

warning('off','MATLAB:interp1:NaNinY');

h_full = NaN*ones(size(Pfull));
h_half = NaN*ones(size(Phalf));

for it=1:size(Pfull,1)
    it
    for ilat=1:size(Pfull,3)
        for ilon=1:size(Pfull,4)


            [P,h]=hyd_solve(Tfull(it,:,ilat,ilon),Pfull(it,:,ilat,ilon),[Pfull(it,end,ilat,ilon) Pfull(it,1,ilat,ilon)]);

            h_full(it,:,ilat,ilon) = interp1(P,h,Pfull(it,:,ilat,ilon));
            h_half(it,:,ilat,ilon) = interp1(P,h,Phalf(it,:,ilat,ilon));

        end
    end
end

h_half(:,end,:,:) = 0;

% --- can't do this because h will be a different size each time....
%use the new multidim interpolation (uses C) to speed up this part -
%should also sort out an ODE solver from C too
%   h_full = lininterp1f_multidim_RUN(P,h,Pfull,2);
%   h_half = lininterp1f_multidim_RUN(P,h,Phalf,2);
  
  
            
            
            
