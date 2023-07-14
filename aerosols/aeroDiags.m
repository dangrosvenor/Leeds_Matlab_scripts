%%run multipop,  Combined_Sdist4 & massbin

clear NaeroDetQ MaeroDetQ flux Mfluxes

iload=0;
%exdir='c:/cygwin/home/login/runs/aeroruns/aerodiags2/results/AeroDiags.mat';


dx=1000;
jjp=500;

if iload==1
    load(exdir); %loads diag().dg and times arrays
end



%MaeroDet=detrain(diag,'MaeroF ',times,dgstrDan); %mass of material detrained at this height for all partitions
%NaeroDet=detrain(diag,'NaeroF ',times,dgstrDan); %number of aerosol detrained at this height


%Mloss=dgrate(diag,'Mloss ',times,dx,jjp,dgstrDan,GridDan); %looks for Mloss and multiplies each rate by the time between diags to get the total for each time
%Nloss=dgrate(diag,'Nloss ',times,dx,jjp,dgstrDan,GridDan);


for iq=16:54 %54
    dgs=strcat('ALL_WQ',num2str(iq));
    [NaeroDetQ(:,:,iq-15),flux(:,:,iq-15)]=detrain(diag,dgs,times,dgstrDan);
    MaeroDetQ(:,:,iq-15)=NaeroDetQ(:,:,iq-15)*mmav2(iq-15);
    Mfluxes(:,:,iq-15)=flux(:,:,iq-15)*mmav2(iq-15);
end

sumNdet=sum(NaeroDetQ,3);
Nflux=sum(flux,3);
sumMdet=sum(MaeroDetQ,3);
Mflux=sum(Mfluxes,3);


