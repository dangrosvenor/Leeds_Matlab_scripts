function run_ukca_ss_flux()

u=[1:50];
for iu=u
[aer_mas_primss,aer_num_primss,aer_numflux_primss] = ukca_prim_ss(u(iu));
flux(iu)=sum(aer_numflux_primss);
end

figure
plot(u,flux); grid on
xlabel('U10 (m s^{-1})');
ylabel('Total SS number flux from surface (# m^{-2} s{-1})');


function [aer_mas_primss,aer_num_primss,aer_numflux_primss] = ukca_prim_ss(u_scalar_10m)
%function [aer_mas_primss,aer_num_primss,aer_numflux_primss] = ukca_prim_ss(u_scalar_10m,P,T)
%function [aer_mas_primss,aer_num_primss] = ukca_prim_ss(mass,u_scalar_10m,aird)

%Replaced the aird input as used in the UKCA code with air Pressure and Temp. (P and T)
%aird was the number of air molecules per cc, i.e. aird = Nair/V.
% So,  P = aird *kT and so aird = P/kT
% Also, don't need mass (the grid box air mass)

%REAL, INTENT(IN)    :: mass(row_length,rows, model_levels)
% mass of each cell
% For Matlab can just pass this as one value - Fortran only uses the lowest vertical level for this
% Guess will need to work out from the dz used in my UM setup?
%REAL, INTENT(IN)    :: u_scalar_10m(row_length,rows)
% Scalar 10m wind
%REAL, INTENT(IN)    :: aird(row_length,rows, model_levels)
% No density (/cm^3)
%REAL, INTENT(INOUT) :: aer_mas_primss(row_length,rows,            ...
%model_levels,nmodes)
% mass ems flux
%REAL, INTENT(INOUT) :: aer_num_primss(row_length,rows,            ...
%model_levels,nmodes)
% number ems flux




% Local variables
%INTEGER :: jl
%INTEGER :: jv
%INTEGER :: modemt
%INTEGER, PARAMETER :: nbins=20
%REAL    :: mbmid(nbins)
%REAL    :: mblo(nbins)
%REAL    :: mbhi(nbins)
%REAL    :: drmid(nbins)
%REAL    :: deltar(nbins)
%REAL    :: drhi(nbins)
%REAL    :: drlo(nbins)
%REAL    :: flux(row_length,rows)
%REAL    :: boxflux(row_length,rows)
%REAL    :: deln(row_length,rows)
%REAL    :: deln_mflux(row_length,rows)
%REAL    :: bet
%REAL    :: agong
%REAL    :: s98flux(row_length,rows)
%REAL    :: m86flux(row_length,rows)
%REAL    :: combflux(row_length,rows)
%REAL    :: dfdr(row_length,rows)
%REAL    :: dum

nbins=20;
mbsmall = 156.0;
mblarge = 4.6e13;

% changed settings above to match bin -- %REALised that difference in
% sea-salt emissions between bin and mode was due to bin upper-limit
% for sea-salt emissions being ~ 7 microns, whereas mode is ~14 microns
% Setting MBSMALL,MBLARGE as above results in equivalent emissions
% grid being used for GLOMAP-mode as in the standard GLOMAP_bin run.
%
% .. use 35 nm cut-off for r at rh=80 (approx 17.5 nm cutoff for dryr)

cutoff=0.0175e-6;

% Molar mass of dry air (kg/mol)
%REAL :: mm_da  % =avc*zboltz/ra

%INTEGER, PARAMETER :: i_method_smith = 1
%INTEGER, PARAMETER :: i_method_monahan = 2
%INTEGER, PARAMETER :: i_method_combined = 3

i_method_smith = 1;
i_method_monahan = 2;
i_method_combined = 3;


avc = 6.022e23; %Avogadro's constant
% Boltzmanns constant (J K-1)
zboltz = 1.3804e-23;
ra = 287.05; %Gas constant for dry air

%% From ukca_mode_setup.F90
cp_cl=4;  % Index to store NaCl   cpt
nmodes=4; %actually set to 7 in ukca_mode_setup, but only use indices 3 and 4 here.

% Mass density of components (kg m^-3)
rhocomp=[1769.0,1500.0,1500.0,1600.0,2650.0,1500.0];

% Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168];

% Specify size limits of geometric mean radius for each mode
% Set dlim34 here to be 500nm to agree with bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6];
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5];


%% End of setup constants, etc.

aer_mas_primss = zeros([size(u_scalar_10m,1) size(u_scalar_10m,2) size(u_scalar_10m,3) nmodes]);
aer_num_primss = zeros(size(aer_mas_primss));
aer_numflux_primss = zeros(size(aer_mas_primss));

%aird = P./(zboltz.*T) / 1e6; %no. air moleculues per cm3

% Molar mass of dry air (kg/mol)
mm_da = avc*zboltz/ra;

seaice_frac = zeros(size(u_scalar_10m));
land_fraction = zeros(size(u_scalar_10m));

% .. Set emission method

i_ems_method = i_method_monahan;

% DEPENDS ON: ukca_ingridg
%CALL ukca_ingridg(nbins,mbsmall,mblarge,mbmid,mblo,mbhi)
[mbmid,mblo,mbhi] = ukca_ingridg(nbins,mbsmall,mblarge);

dum=3.0*mm(cp_cl)/(avc*rhocomp(cp_cl)*4.0*pi);
for jv=1:nbins
    deltar(jv)=((dum*mbhi(jv))^(1.0/3.0)-                          ...
        (dum*mblo(jv))^(1.0/3.0))*1.0e6; % in microns
    drmid(jv)=(dum*mbmid(jv))^(1.0/3.0);         % in m
    drlo(jv) =(dum*mblo(jv ))^(1.0/3.0);         % in m
    drhi(jv) =(dum*mbhi(jv ))^(1.0/3.0);         % in m
end

% .. Loop over sea-salt emission bins
for jv = 1:nbins
    if (drhi(jv) > cutoff)       % if bin size is > cutoff for flux
        %                   (cutoff is lowest dry radius for the Gong-Monahan scheme)
        
        modemt=0;
        if (drmid(jv) <  0.5*ddplim0(4)); modemt=3; end    % soluble accum. mode
        if (drmid(jv) >= 0.5*ddplim0(4)); modemt=4; end    % soluble coarse mode
        
        % .. now emit according to widths as set in DDPLIM0 rather than
        % .. as hard-coded at 1 micron dry diameter (equivalent for usual setup)
        
        if (modemt > 0)
            
            if (i_ems_method == i_method_smith |                     ...
                    i_ems_method == i_method_combined)
                s98flux(:,:) = 0.0;
                %        WHERE (land_fraction < 0.999 .AND. seaice_frac < 0.999)
                % .. Smith et al. (1998)
                s98flux(:,:) = 1.4367*0.2*(u_scalar_10m(:,:)^3.5)*     ...
                    exp(-1.5*(log(2.0*drmid(jv)/2.088e-6 ))^2) +         ...
                    1.4367 * 6.8e-3*(u_scalar_10m(:,:)^3  )*                    ...
                    exp(-1.0*(log(2.0*drmid(jv)/20.88e-6))^2);
                %        END WHERE
            end
            
            if (i_ems_method == i_method_monahan |                     ...
                    i_ems_method == i_method_combined)
                m86flux(:,:) = 0.0;
                bet = (0.433 - log10(drmid(jv)*2.0e6))/0.433;
                agong = 4.7*(1.0 + 30.0*(drmid(jv)*2.0e6))^(-0.017*      ...
                    (drmid(jv)*2.0e6)^(-1.44));
                %        WHERE (land_fraction < 0.999 .AND. seaice_frac < 0.999)
                % ..  Gong (2003) ---- extension from Monahan et al. (1986)
                %          (factor of 2.0 to convert dF/dr(80) to dF/dr(dry) )
                % DPG - Uses Eqn. 2 of Gong (2003). Gives dF/dr, i.e. number flux per radius
                % interval of size dist.
                m86flux(:,:) = 2.0*1.373*(u_scalar_10m(:,:)^3.41)*     ...
                    (drmid(jv)*2.0e6)^(-agong)*             ...
                    (1.0+0.057*((2.0e6*drmid(jv))^(3.45)))* ...
                    10.0^(1.607*exp(-bet*bet));
                %        END WHERE
            end
            
            
            
            % .. Choose Smith, Monahan, or combined flux
            if (i_ems_method == i_method_smith)
                dfdr(:,:) = s98flux(:,:);
                %        if (verbose >= 2) CALL umPrint(                        ...
                %            'Smith sea-salt flux scheme selected',             ...
                %            src='ukca_prim_ss')
            elseif (i_ems_method == i_method_monahan)
                dfdr(:,:) = m86flux(:,:);
                %        if (verbose >= 2) CALL umPrint(                        ...
                %            'Monahan sea-salt flux scheme selected',           ...
                %            src='ukca_prim_ss')
            elseif (i_ems_method == i_method_combined)
                % .. If using a combination of Smith and Monahan, then the recommendation
                %     is to use the larger of the two fluxes at any radius
                %        if (verbose >= 2) CALL umPrint(                        ...
                %            'Combined sea-salt flux scheme selected',          ...
                %            src='ukca_prim_ss')
                %        WHERE (s98flux > m86flux)
                %          combflux = s98flux
                %        ELSEWHERE
                %          combflux = m86flux
                %        END WHERE
                iwhere = s90flux>m86flux;
                combflux(iwhere) = s98flux;
                combflux(~iwhere) = m86flux;
                dfdr(:,:) = combflux(:,:);
            else
                %        cmessage = ' No method for sea-salt flux specified'
                error(' No method for sea-salt flux specified');
                %        errcode = 1
                %        CALL ereport('ukca_prim_ss',errcode,cmessage)
            end
            
            % .. Calculate sea-salt emission flux in particles m-2 s-1
            if (drlo(jv) > cutoff)
                flux(:,:) = dfdr(:,:)*deltar(jv);
            else
                flux(:,:) = dfdr(:,:)*(drhi(jv) - cutoff)*1.0e6;
            end
            
            % .. Calculate  sea-salt aerosol source (particles m-2 s-1)
            boxflux(:,:) = flux(:,:)*(1.0-seaice_frac(:,:))*            ...
                (1.0-land_fraction(:,:));
            
            % .. Calculate change in grid box aerosol number density (particles/cc m-2 s-1)
            % DPG - boxflux should be in #/m2/s at this stage (as in Gong paper)
            % So, here we just divide by the volume of the air (in cm3) = Vair.
            % Vair = Nair / aird (from aird = Nair/V, since aird is the number density of air molecues per cc)
            % Nair = no. moles * avc (avc=Avogadro)
            % Nair = mass * avc / mm_da  (since no. moles = mass/molecular mass)
%DPG - deln is not acutally used, though, so don't need the air mass...
            %deln(:,:) = boxflux(:,:)*aird(:,:,1)*                       ...
            %    mm_da/(mass(:,:,1)*avc);
            
            % .. Calculate change in grid box aerosol mass density (#/m2/s)
            deln_mflux(:,:) = boxflux(:,:)*mm_da/avc;
            
            % .. sum up each emitted mass into mode in kg/m2/s
            aer_mas_primss(:,:,1,modemt) =                              ...
                aer_mas_primss(:,:,1,modemt) + deln_mflux*mbmid(jv);
            % .. sum up each emitted number into mode in equiv-kg/m2/s
            aer_num_primss(:,:,1,modemt) =                              ...
                aer_num_primss(:,:,1,modemt) + deln_mflux(:,:);

             % .. sum up each emitted boxflux number into mode in #/m2/s
            aer_numflux_primss(:,:,1,modemt) =                              ...
                aer_numflux_primss(:,:,1,modemt) + boxflux(:,:);

                      
        end    % if modemt > 0
    end     % if drhi(jv) > cutoff
end      % jv = 1,nbins







function [mbmid,mblo,mbhi] = ukca_ingridg(nbins,mbsmall,mblarge)

lgmrange=log(mblarge)-log(mbsmall);
for jv=1:nbins+1
    lgmgrid = log(mbsmall) + lgmrange * (jv-1)/(nbins);
    if (jv < (nbins+1)); mblo(jv)=exp(lgmgrid); end
    if (jv > 1); mbhi(jv-1)=exp(lgmgrid); end
end

for jv=1:nbins
    mbmid(jv)=exp(0.5*log(mblo(jv)*mbhi(jv))); % geometric mean
end









