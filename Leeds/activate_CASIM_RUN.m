%calls activate_CASIM.m

nrun=1; %simulate nrun cloud activations (run activate nrun times)


dt=1; %doesn't really matter about the timestep - set to one since answers get divided by this

fact=0.75;


fevap = 0.98; %fraction of droplets that evporate each cycle for testing.
fevap=0;

cloud_mass = 0.2e-3; %kg/kg, approx 0.2 g/m3
cloud_number_init = 0; % per kg
w = 0.5; %m/s
rho = 0.8; %kg/m3
T = 282; %K
p = 800e2; %hPa


%Values in the namelist (initial values)
mass_accum_nml = 4.5e-9;  %kg per kg
N_accum_nml = 3.8e8; %per kg, approx 380 per cc

%initial defaults
dnccn_all=0;
dmac_all=0;
dnumber=0;
% interstitial aerosol number conc per kg
aerophys.N(1) = 0;  %Aitken mode
aerophys.N(3) = 0; % Coarse mode (approx 380 per cc)

% interstitial aerosol mass, kg per kg
aerophys.Minit(1) = 0;
aerophys.Minit(2) = mass_accum_nml;
aerophys.Minit(3) = 0;

aerophys.Ninit(1) = 0;  %Aitken mode
aerophys.Ninit(2) = N_accum_nml; % Accum mode (approx 380 per cc)
aerophys.Ninit(3) = 0; % Coarse mode (approx 380 per cc)

aerophys.N = aerophys.Ninit;
aerophys.M = aerophys.Minit;
cloud_number = cloud_number_init;

%simulate nrun cloud activations (run activate nrun times)
for irun=1:nrun

    aerophys.N = aerophys.N - dnccn_all;  %Aitken mode, accum and coarse modes
    aerophys.M = aerophys.M - dmac_all;  %Aitken mode, accum and coarse modes
    cloud_number = (1 - fevap) * (cloud_number + dnumber)


    aeroact = []; %don't think this does anything.. used as optional for ARG scheme
    dustchem=[]; %not used either
    dustliq=[]; %ditto

    dustphys.N(1)=0;
    dustphys.M(1)=0;

    for imode=1:length(aerophys.N)
        aerochem.density(imode) = 1777; %kg/m3 - fixed aerosol density (used for all modes)
        aerophys.sigma(imode) = 1.5; %default for aerosol
        aerophys.rd(imode) = MNtoRm(aerophys.M(imode), aerophys.N(imode), aerochem.density(imode), aerophys.sigma(imode));
    end


    [dnumber, dmac, dnccn_all, dmac_all, dnumber_d, dmass_d, active] = activate_CASIM(fact,dt, cloud_mass, cloud_number, w, rho,  T, p, ...
        aerophys, aerochem, aeroact, dustphys, dustchem, dustliq);

    act_frac_mass(irun,:) = dmac_all./aerophys.M;    
    
    dnccn_all_irun(irun,:) = dnccn_all;
    
    %Amount of aerosol as proportion of the initial value that has
    %activated cumulatively over all activation cycles
    act_frac_mass_cum(irun,:) = (aerophys.Minit - (aerophys.M - dmac_all) )./aerophys.Minit;
    
    act_frac_num_cum(irun,:) = (aerophys.Ninit - (aerophys.N - dnccn_all) )./aerophys.Ninit;    
end

act_frac_mass(1:nrun,:)
act_frac_mass_cum(1:nrun,:)
act_frac_num_cum(1:nrun,:)
dnccn_all_irun(1:nrun,:)
