dirpath='C:\Documents and Settings\Login\My Documents\samis_model\icemodel1\';

 out_folder{1}='output_1st_run\';
 out_folder{2}='output\';
 out_folder{3}='output_HM_20\';

 runs=[1 2 3];
 
iread=0;
iread_ice=1;

fig=3;

%
plot_RH_T=1;    % piirret��nk� RH ja l�mp�tila ajan funktiona
plot_cgas=1;    % piirret��nk� kaasujen pitoisuudet ajan funktiona
plot_distrib=1; % piirret��nk� hiukkasten kokojakaumat
plot_kohler=1;  % piirret��nk� K�hler k�yr�
plot_final_c=1; % piirret��nk� kokojakauma lopussa
%

if iread==1
    fprintf(1,'\nReading in cmodel data from %s%s.....',dirpath,out_folder);
	[n_time,nbin,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry]=read_data(dirpath,out_folder);
    fprintf(1,'\nFinshed reading.');
	n_bins=sum(nbin);    
end

if iread_ice==1    
    for iruns=1:length(runs)
        irun=runs(iruns);
        fprintf(1,'\nReading in cmodel data from %s%s.....',dirpath,out_folder{irun});
        [n_time(irun).dat,n_ice(irun).dat,n_drops(irun).dat,time(irun).dat,r_ice(irun).dat,c_ice(irun).dat,r_drops(irun).dat,c_drops(irun).dat,...
                RH(irun).dat,c_ice_tot(irun).dat,sum_c(irun).dat,c_end(irun).dat,sum_cy1(irun).dat,sum_cy2(irun).dat,aa(irun).dat]=read_data_ice(dirpath,out_folder{irun});
        fprintf(1,'\nFinshed reading.');
    end
end

% n_time                    Aika-arvojen lukum��r�
% n_bins                    Binien lukum��r�
% n_gases                   Kaasujen lukum��r�
% time(n_time)              Aika (s)
% r_bin(n_time,n_bins)      Hiukkasten s�de ajan funktiona
% RH(n_time)                RH ajan funktiona
% T(n_time)                 L�mp�tila (K) ajan funktiona
% c_gases(n_time,n_gases)   Kaasujen pitoisuus ajan funktiona (mol/cm^3)
% c_bin(n_bins)             Hiukkasten lukum��r�pitoisuus alussa (1/cm^3)
% r_dry(n_bins)             Hiukkasten s�de ennen alkutasapainoa (m)
%
if iread_ice~=1
	fprintf(1,'\nAbout to plot....');
	f=cmodel_plot_dan(fig,n_time,nbin,n_gases,time,r_bin,RH,T,c_gases,c_bin,r_dry,n_bins);
end