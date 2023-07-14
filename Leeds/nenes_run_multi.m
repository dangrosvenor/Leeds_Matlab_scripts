clear nccn_all_Nenes 

M_aero0 = 1.e-7; % Total aerosol mass mixing ratio kg/kg

dw = 0.1;
w = [dw:dw:1]; % m/s

df=0.1;
factor = [df:df:1];

%df=1;
%factor = [df:df:10];

for k=1:length(factor)  %in range(0,len(factor)):
    %    nccn=zeros([length(w) nmodes]);

    switch N_type
        case 'prescribed number'
            %Vary N_aero, but keep M_aero0 the same
            M_aero = M_aero0;
            N_aero = factor(k)*N_aero0
            % N_aero set outside
            tit_str = 'for varying number (constant mass)';

        case 'prescribed single aerosol mass'
            %Vary M_aero and keep r_single (translates to mass for monomodal dist as used in CASIM code to get number) the same
            %This single mass is then used to get he number, so number will
            %scale with the mass (i.e. by factor(k)).
            M_aero = factor(k)*M_aero0;
            N_aero = 1; %dummy value - will be calculated in the function based on M_aero and r_single (size)
            tit_str = 'for varying mass and number';

        case 'prescribed median radius (rd)'

    end

    %N.B. - think that the Nenes function needs w as cm/s, so multiplying
    %by 100 here
    [nccn_all_Nenes(k,:),Dg,N_aero_out,r_single] = nenes_run_func(N_aero,M_aero,w*100,r_single,N_type);

end

%% Make plots

%    legend(leg);    
    figure
    
    w2=[w(1)-dw w w(end)+dw];
    w2=0.5*[w2(2:end) + w2(1:end-1)];
        

%    num = factor.*Maero0;
    f2=[factor(1)-df factor factor(end)+df];    
    f2=0.5*[f2(2:end) + f2(1:end-1)];  
    
    dpcolor(w2,f2,nccn_all_Nenes/1e6); %Converted to per mg (same as per cc for air density=1)
    shading flat %Note shading interp messes up the axes and gives a blank row and column...
    colorbar
    xlabel('W (m s^{-1})','fontsize',15);
%    ylabel('Nc (cm^{-3})','fontsize',15);
    ylabel(['factor ' tit_str],'fontsize',15);
    title(['rmin=' num2str(r_single*1e9,'%.0f') ' nm, rdmin=' num2str(Dg/2*1e9,'%.0f') ' nm, N_{aero MAX}=' num2str(N_aero_out(1)/1e6,'%.0f') ' mg^{-1}'],'fontsize',15);
