function [dm1, dm2, dn1, dn2, dm_ratio_accum, dn_ratio_accum,r_new_accum,r_new_coarse,rm,r_ratio,r1,r2] = which_mode(dm, dn, method)
%    real(wp), intent(in) :: dm  % total mass increment
%     real(wp), intent(in) :: dn  % total number increment
%
%     real(wp), intent(in) :: r1_in  % mean radius of mode 1
%     real(wp), intent(in) :: r2_in  % mean radius of mode 2
%
%     real(wp), intent(in) :: density % density of aerosol, assumed the same across all modes
%
%     real(wp), intent(out) :: dm1  % mass increment to mode 1
%     real(wp), intent(out) :: dn1  % number increment to mode 1
%     real(wp), intent(out) :: dm2  % mass increment to mode 2
%     real(wp), intent(out) :: dn2  % number increment to mode 2

% Orig form - %function which_mode(dm, dn, r1_in, r2_in, density, dm1, dm2,dn1, dn2)

% method = 'Dan Robin';
% method = 'Robin';
%method = 'lognormal dN/lnD ratio';
%method = 'UM';

switch method
    case 'UM pre bug-fix'
        i_use_MNtoRm=0;
    otherwise
        i_use_MNtoRm=1;
end

global ccn_tidy

density = 1777; %kg/m3 - fixed aerosol density (used for all modes)
sigma = 1.5; %default for aerosol

ccn_tidy = 0.1e+6; %in thresholds.F90

%Values in the namelist (initial values)
mass_accum_nml = 4.5e-9; 
N_accum_nml = 3.8e8;

f_init = 0.5;
%f_init = 1e-5;
mass_accum = mass_accum_nml*f_init; %Assume that currently have half the mass and number of the accumluation mode
N_accum = N_accum_nml*f_init; %This will give the same size.


%Will calc these based on the namelist values (as done initially in code)
%instead of using input values.
r1_in = MNtoRm(mass_accum, N_accum , density, sigma);
r2_in = MNtoRm(0, 0 , density, sigma);

% r1_in = 1.71593852907412857E-7;
% r2_in = 0.0;
% 
% r1_in = 1.7113646216939575E-7;
% r2_in = 1.7113646216939575E-7;



max_accumulation_mean_radius = 0.25e-6;  %orig
%max_accumulation_mean_radius = 0.26e-6;
min_coarse_mean_radius = 1.0e-6;  %orig
%min_coarse_mean_radius = 0.25e-6;
%min_coarse_mean_radius = 2.0e-6;

default_coarse_mean_radius = 1.0e-6;


dm1=0.0;
dm2=0.0;
dn1=0.0;
dn2=0.0;

l_aeroproc_no_coarse=false;
%l_aeroproc_no_coarse=true;


%    if (l_ukca)
%        imethod=iukca_method;
%    end

if (dm*dn > 0.0)
    % dm and dn should be positive and of the same sign
    r1=min(max_accumulation_mean_radius, r1_in)
    
    if r2_in<1e-30
        r2 = default_coarse_mean_radius
    else
        r2 = r2_in
    end    
%    r2=max(min_coarse_mean_radius, r2_in)

    ftpi=3.14159*4./3.0;
    rftpi=1./ftpi;
    
    if i_use_MNtoRm==0
        rm=(rftpi*dm/dn/density).^(1.0/3.0)  %Should use MNtoRm(mode_M,mode_N,density,aerophys(k)%sigma(imode)) here I think for consistency.
    else
        rm = MNtoRm(dm, dn, density, sigma)
    end
    
    %      select case(imethod)
    %      case (iukca_method)
    % Call ukca subroutine, or add UKCA code here?
    %      case default
    
%    if (rm > r2 && ~l_aeroproc_no_coarse)

                r1_3=r1*r1*r1;
                r2_3=r2*r2*r2;


    if (rm > r2)
        dm1=0.0;
        dn1=0.0;
        dm2=dm;
        dn2=dn;
    elseif (rm < r1 | l_aeroproc_no_coarse)
        %Force into here if l_aeroproc_no_coarse - i.e. no processing into the
        %coarse mode
        dm1=dm;
        dn1=dn;
        dm2=0.0;
        dn2=0.0;
    else
        switch method
            case {'UM','UM pre bug-fix'}

                gamma=(rm*rm*rm-r1_3)/(r2_3-r1_3);
                %        gamma=(rm-r1)/(r2-r1); %This makes it worse...
                dn2=dn*(gamma);
                dn1=dn - dn2;
                if i_use_MNtoRm==0
                    dm2=ftpi*density*r2_3*dn2;
                else
                    dm2=ftpi*density*r2_3*dn2 * exp(4.5*log(sigma).^2); %4/3 * pi * rho * integral(N * r.^3)
                end
                dm1=dm-dm2;

            case 'lognormal dN/lnD ratio'
                
                %Work out the ratio of dn/dlnD of the accum and coarse
                %modes to the total dn/dlnD at the mean radius of the
                %aerosol to be added and apportion number according to
                %that.
                %Then do the same for dV/dlnD and do the same for the mass
                
                dNdlnD_accum = lognormal(rm*2,dn,sigma,r1*2);
                dNdlnD_coarse = lognormal(rm*2,dn,sigma,r2*2);  
                gamma = dNdlnD_coarse ./ (dNdlnD_accum + dNdlnD_coarse)
                
                dVdlnD_accum = lognormal_volume(rm*2,dn,sigma,r1*2);
                dVdlnD_coarse = lognormal_volume(rm*2,dn,sigma,r2*2); 
                gamma_vol = dVdlnD_coarse ./ (dVdlnD_accum + dVdlnD_coarse)
                
                %Note, that these ratios are actually the same since if we
                %use the same N (=dn), sigma and diameter (of evaporating
                %paricles) the stuff in front of the expoentials cancels
                %leaving only the exponentials, whihc are the same for
                %dV/dlnD and for dN/dlnD. 
                % But prob should use the numbers in each dist for N rather
                % than dn - but problematic if we have no current coarse
                % aerosol!
                 

%                 dn2=dn*(gamma);
%                 dn1=dn - dn2;  %This is wrong - we do not want to conserve
%                 if i_use_MNtoRm==0
%                     dm2=ftpi*density*r2.^3*dn2;
%                 else
%                     dm2=ftpi*density*r2_3*dn2 * exp(4.5*log(sigma).^2); %4/3 * pi * rho * integral(N * r.^3)
%                 end
%                 dm1=max(0.0,dm-dm2); % N.B. - we shouldn't really expect to conserver the number of the evaporating aerosol since if we transfer to the coarse mode
%                                      % then we are increasing it's size!
%                                      %Could base the ratio on dV/dlnD and
%                                      %then calculate the numbers?

                   dn2 = dn * gamma;
                   dn1 = dn - dn2;
                   dm2 = dm * gamma_vol;
                   dm1 = dm - dm2;
                   
            case 'Robin'
                
                gamma_num = (rm - r1)/(r2 - r1)                
                gamma_vol = (rm*rm*rm-r1_3)/(r2_3-r1_3)
                
                dn2 = dn * gamma_num;
                dn1 = dn - dn2;
                dm2 = dm * gamma_vol;
                dm1 = dm - dm2;
                
            case 'Robin2'
                
                gamma_vol = (rm - r1)/(r2 - r1);                
                gamma_num = (rm*rm*rm-r1_3)/(r2_3-r1_3);
                
                dn2 = dn * gamma_num;
                dn1 = dn - dn2;
                dm2 = dm * gamma_vol;
                dm1 = dm - dm2;
                
            case 'Dan Robin'
                
%                gamma_num = (rm - r1)/(r2 - r1)                
                gamma_vol = (rm*rm*rm-r1_3)/(r2_3-r1_3)
                gamma_num = gamma_vol;
                
                dn2 = dn * gamma_num;
                dn1 = dn - dn2;
                dm2 = dm * gamma_vol;
                dm1 = dm - dm2;   
                
                %If we want the new aerosol to be added whilst preserving
                %the size of the evaporating aerosol then would just add
                %the mass and number according to the same ratio?
                
            case 'Dan Robin2'

                gamma_num = (rm - r1)/(r2 - r1);
                %gamma_vol = (rm*rm*rm-r1_3)/(r2_3-r1_3)
                gamma_vol = gamma_num;

                dn2 = dn * gamma_num;
                dn1 = dn - dn2;
                dm2 = dm * gamma_vol;
                dm1 = dm - dm2;

                %If we want the new aerosol to be added whilst preserving
                %the size of the evaporating aerosol then would just add
                %the mass and number according to the same ratio?
                
            case 'midway all or nothing'
                if rm>(r1+r2)/2
                    dm2=dm;
                    dn2=dn;
                else
                    dm1=dm;
                    dn1=dn;
                end

        end
    end
    %      end select
    
    dm_ratio_accum = dm1./(dm1+dm2);
    dn_ratio_accum = dn1./(dn1+dn2);
    
    r_ratio = (rm-r1)/(r2-r1)
    gamma_old = (rm*rm*rm-r1_3)/(r2_3-r1_3)
    
    M_new_accum = mass_accum + dm1;  %Work out new mass based on half the original namelist values for mass and number
    N_new_accum = N_accum + dn1;    %Should prob do all this outside the function
    
    M_new_coarse = 0 + dm2;
    N_new_coarse = 0 + dn2;  
    
    r_new_accum = MNtoRm(M_new_accum, N_new_accum, density, sigma);
    r_new_coarse = MNtoRm(M_new_coarse, N_new_coarse, density, sigma);
    
    r_new_ratio_accum = r_new_accum ./r1
    r_new_ratio_coarse = r_new_coarse ./r2
    

end


function [Rm]=MNtoRm(M, N, density, sigma)
%    real(wp), intent(in) :: M, N, density, sigma
%    real(wp) :: MNtoRm

global ccn_tidy

%    if (N==0 | M==0) ! shouldn't really be here
    if (N<1e-10 | M<1e-20) ! shouldn't really be here        
%    if (N<ccn_tidy | M<ccn_tidy*1e-18) ! shouldn't really be here        
      Rm=0.0;
    else
      Rm=(3.0*M*exp(-4.5*log(sigma).^2)/(4.0*N*pi*density)).^(1.0/3.0);
%      Rm=(3.0*M*exp(-4.5*log(sigma).^2)/(4.0*N*3.00*density)).^(1.0/3.0);      
    end