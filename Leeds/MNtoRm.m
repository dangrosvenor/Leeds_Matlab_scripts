function [Rm]=MNtoRm(M, N, density, sigma)
%    real(wp), intent(in) :: M, N, density, sigma
%    real(wp) :: MNtoRm

global ccn_tidy

%    if (N==0 | M==0) % shouldn't really be here
    if (N<1e-10 | M<1e-20) % shouldn't really be here        
%    if (N<ccn_tidy | M<ccn_tidy*1e-18) % shouldn't really be here        
      Rm=0.0;
    else
      Rm=(3.0*M*exp(-4.5*log(sigma).^2)/(4.0*N*pi*density)).^(1.0/3.0);
%      Rm=(3.0*M*exp(-4.5*log(sigma).^2)/(4.0*N*3.00*density)).^(1.0/3.0);      
    end