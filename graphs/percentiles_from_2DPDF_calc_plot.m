%deals with overlaying percentiles from 2D PDFs onto other 2D PDFs
% two modes of action :-
% 1) Calculate the percentiles from a 2D PDF that has just been drawn and save
% them
% 2) Load in the percentiles and plot then over the current 2D PDF plot
% Set up the desired variable in pdf2d first

ino_polder=0;   %flag to say whether to avoid plotting POLDER percentiles - is a problem in that the PDF cuts off
%at 20um. So, means and median will be underestimated.

plot_pcolor='yes';
plot_pcolor='no';

action = 'calc and save';
action = 'load and plot';

mod35='no';
mod35='yes';


savedir_prctiles = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/';



%% load and plot
switch action
    case 'load and plot'
        leg_loc='NorthWest';
        
        switch plot_pcolor
            case 'yes'
                plotTimeHeightVap3
            case 'no'  %just do the x-y plot with no colors
                plotTimeHeightVap3
                posit2=posit*0.7;
                figure('position',[posit2(1) posit2(2) posit2(3)*1.1 posit2(4)],'name',[figlab ' ' xlabelstr ' ' ylabelstr]);
                xlabel(xlabelstr);
                ylabel(ylabelstr);
                title(titlenam);
                fontsize_figure(gcf,gca,16);
                grid minor
        end

        switch y_axis_vals
            case {'CF COSP-MODIS GCM','Max low cloud (no screening)'}
                
               
                
                %GCM
                %load the latest filename from the appropriate FILENAME .mat file
                savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_' gcm_str '_' xlabelstr '_' ylabelstr '.mat'];
                savefilenames_prctiles = remove_problem_chars(savefilenames_prctiles);
                %Below will load ['savefile_prctiles_' gcm_str]
                load(savefilenames_prctiles);
                
                %Obs
                %load the latest AMSRE filename from the appropriate FILENAME .mat file
                switch mod35
                    case 'no'                        
                        savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_MODIS_Longitude_Cloud_Fraction.mat'];
                        leg_str1='CFliq';
                    case 'yes'
                        savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_MODIS_Longitude_MOD35_Daytime_Cloud_Fraction.mat'];
                        leg_str1='MOD35';
                end

                load(savefilenames_prctiles);
                savefile_prctiles_OBS=savefile_prctiles_MODIS;  
                
                
            case {'CDR Polder2','Re COSP GCM','R_{eff 2.1 \mum} (\mum) minus 20%'}
                leg_str1='POLDER';
                
                %modis minus 20%, CF>0.8, NP>50
                savefile_prctiles_STORE{1}=[savedir_prctiles 'saved_prctiles_2Dpdf_Longitude_R_{eff_2.1_mum}_(mum)_minus_20pct_20130208T141551.mat'];
                %MODIS minus 20%, no screening
                savefile_prctiles_STORE{1}=[savedir_prctiles 'saved_prctiles_2Dpdf_Longitude_R_{eff_2.1_mum}_(mum)_minus_20pct_20130208T182040.mat'];

                %POLDER
                %savefile_prctiles_STORE{2}=[savedir_prctiles 'saved_prctiles_2Dpdf_Longitude_CDR_(mum)_20130208T152720.mat'];
                savefile_prctiles_STORE{2}=[savedir_prctiles 'saved_prctiles_2Dpdf_Longitude_CDR_(mum)_20130208T175744.mat'];
                % with PDFs normalised along y
                savefile_prctiles_OBS=[savedir_prctiles 'saved_prctiles_2Dpdf_POLDER_Longitude_CDR_(mum)_20130209T112259.mat'];

                %GCMs
                %2deg spacing
                %CAM5_CLUBBv2
                savefile_prctiles_CAM5_CLUBBv2_COSP=[savedir_prctiles 'saved_prctiles_2Dpdf_CAM5_CLUBBv2_COSP_Longitude_COSP_Reff_20130210T200056.mat'];
                %CAM5
                savefile_prctiles_CAM5_COSP=[savedir_prctiles 'saved_prctiles_2Dpdf_CAM5_COSP_Longitude_COSP_Reff_20130209T113715.mat'];
                %CAM-CLUBBv1
                savefile_prctiles_CAM5_CLUBB_COSP=[savedir_prctiles 'saved_prctiles_2Dpdf_CAM5_CLUBB_COSP_Longitude_COSP_Reff_20130210T195533.mat'];

                
                %load the latest GCM filename from the appropriate FILENAME .mat file
                savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_' gcm_str '_' xlabelstr '_' ylabelstr '.mat'];
                savefilenames_prctiles = remove_problem_chars(savefilenames_prctiles);
                load(savefilenames_prctiles);
                
                leg_loc = 'NorthEast';
                
                
            case {'Nd from grid vals timeseries3','Nd GCM'}
                leg_str1='MODIS';
                %modis, CF>0.8, NP>50
                savefile_prctiles_OBS=[savedir_prctiles 'saved_prctiles_2Dpdf_MODIS_Longitude_N_d_(cm^{-3})_20130210T210959.mat'];

                %GCMs
                %2deg spacing
                %CAM5_CLUBBv2
                savefile_prctiles_CAM5_CLUBBv2_COSP=[savedir_prctiles 'saved_prctiles_2Dpdf_CAM5_CLUBBv2_COSP_Longitude_Nd_(cm^{-3})_20130210T210635.mat'];
                %CAM5
                savefile_prctiles_CAM5_COSP=[savedir_prctiles 'saved_prctiles_2Dpdf_CAM5_COSP_Longitude_Nd_(cm^{-3})_20130210T211850.mat'];
                %CAM-CLUBBv1
                savefile_prctiles_CAM5_CLUBB_COSP=[savedir_prctiles 'saved_prctiles_2Dpdf_CAM5_CLUBB_COSP_Longitude_Nd_(cm^{-3})_20130210T211622.mat'];
                %AM3-CLUBB
                savefile_prctiles_AM3_CLUBB=[savedir_prctiles 'saved_prctiles_2Dpdf_AM3_CLUBB_Longitude_Nd_(cm^{-3})_20130210T214059.mat'];


                %load the latest GCM filename from the appropriate FILENAME .mat file
                savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_' gcm_str '_' xlabelstr '_' ylabelstr '.mat'];
                savefilenames_prctiles = remove_problem_chars(savefilenames_prctiles);
                load(savefilenames_prctiles);
                


                
                
            case {'TLWP DAYTIME GCM, grid-box average','TLWP NIGHTTIME GCM, grid-box average','LWP DAYTIME GCM, grid-box average','LWP NIGHTTIME GCM, grid-box average'}
                leg_str1='AMSRE';
                
                %load the latest filename from the appropriate FILENAME .mat file
                savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_' gcm_str '_' xlabelstr '_' ylabelstr '.mat'];
                savefilenames_prctiles = remove_problem_chars(savefilenames_prctiles);
                %Below will load ['savefile_prctiles_' gcm_str]
                load(savefilenames_prctiles);
                
                switch y_axis_vals
                    case {'TLWP DAYTIME GCM, grid-box average','LWP DAYTIME GCM, grid-box average'}
                        %load the latest AMSRE filename from the appropriate FILENAME .mat file
                        savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_AMSRE_Longitude_Daytime_AMSRE_TLWP_(g_m^{-2}).mat'];

                    case {'TLWP NIGHTTIME GCM, grid-box average','LWP NIGHTTIME GCM, grid-box average'}
                        %load the latest AMSRE filename from the appropriate FILENAME .mat file
                        savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_AMSRE_Longitude_Nighttime_AMSRE_TLWP_(g_m^{-2}).mat'];
                end
                
                load(savefilenames_prctiles);
                savefile_prctiles_OBS=savefile_prctiles_AMSRE;                    

        end




        switch gcm_str
            case 'POLDER'
                %load and plot MODIS
                savefile_prctiles=savefile_prctiles_STORE{1};
                load(savefile_prctiles);
                hold on
                %errorbarYY plots an error bar from a given point with the
                %positive and negative errors the same - so plot from the
                %midpoint and give it half the amplitude
                mid_errorbar_pos = 0.5*(prcY_2Dpdf(1,:)+prcY_2Dpdf(3,:));
                L_errorbar = (prcY_2Dpdf(3,:)-prcY_2Dpdf(1,:))/2;
                spanX = mid_Xbins(end) - mid_Xbins(1);
                errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'k','none',5,0.015,spanX);
                errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'w','none',2,0.015,spanX);
                xsave=mid_Xbins; ysave=prcY_2Dpdf(2,:);

                if ino_polder==0
                    %load and plot POLDER
                    savefile_prctiles=savefile_prctiles_STORE{2};
                    load(savefile_prctiles);
                    hold on

                    mid_errorbar_pos = 0.5*(prcY_2Dpdf(1,:)+prcY_2Dpdf(3,:));
                    L_errorbar = (prcY_2Dpdf(3,:)-prcY_2Dpdf(1,:))/2;
                    spanX = mid_Xbins(end) - mid_Xbins(1);
                    errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'k','none',5,0.010,spanX);
                    errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'g','none',2,0.010,spanX);

                end

                plot(xsave,ysave,'ks','markerfacecolor','w','markersize',10);

                if ino_polder==0
                    plot(mid_Xbins,prcY_2Dpdf(2,:),'kd','markerfacecolor','g','markersize',10);
                end

            otherwise
                %load and plot POLDER (reff) or MODIS (Nd)
%                savefile_prctiles=savefile_prctiles_STORE{2};
                savefile_prctiles=savefile_prctiles_OBS;                
                load(savefile_prctiles);
                hold on

                per10 = prcY_2Dpdf(1,:);               
                mid_errorbar_pos = 0.5*(per10+prcY_2Dpdf(3,:));
                L_errorbar = (prcY_2Dpdf(3,:)-per10)/2;
                spanX = mid_Xbins(end) - mid_Xbins(1);
%                errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'k','none',5,0.015,spanX);
                errorbarYY('vert filled bar',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'g','none',2,0.015,spanX);                
                xsave=mid_Xbins;
                ysave=prcY_2Dpdf(2,:);
                Y_mean_save = Y_mean;
              

                %load and plot GCM percentiles
                savefile_prctiles=eval(['savefile_prctiles_' gcm_str]);
                load(savefile_prctiles);
                hold on

                per10 = prcY_2Dpdf(1,:);
                per10(isnan(per10))=0;
                mid_errorbar_pos = 0.5*(per10+prcY_2Dpdf(3,:));
                L_errorbar = (prcY_2Dpdf(3,:)-per10)/2;
                spanX = mid_Xbins(end) - mid_Xbins(1);
                errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'k','none',5,0.015,spanX);
                errorbarYY('vert',mid_Xbins,mid_errorbar_pos,L_errorbar,gca,'w','none',2,0.015,spanX);

                
                h1=plot(xsave,ysave,'kd','markerfacecolor','g','markersize',16);
                h2=plot(mid_Xbins,prcY_2Dpdf(2,:),'kd','markerfacecolor','w','markersize',14);
                h3=plot(xsave,Y_mean_save,'k^','markerfacecolor','g','markersize',16);
                h4=plot(mid_Xbins,Y_mean,'k^','markerfacecolor','w','markersize',14);
                
                leg_str2=remove_character(gcm_str,'_','-');

        end

        
        switch y_axis_vals
            case {'Nd from grid vals timeseries3','Nd GCM'}
                set(gca,'ylim',[0 300]);
            case {'TLWP DAYTIME GCM, grid-box average','TLWP NIGHTTIME GCM, grid-box average','LWP DAYTIME GCM, grid-box average','LWP NIGHTTIME GCM, grid-box average'}
                set(gca,'ylim',[-20 250]);
                leg_loc = 'NorthWest';
            case {'CF COSP-MODIS GCM','CF GCM'}
                set(gca,'ylim',[0 1.6]);
                leg_loc = 'NorthWest';
        end

%        legend([h1 h2 h2 h4 0],{leg_str1,leg_str2,'Medians','Means','Bars: 10th & 90th percentiles'},'location',leg_loc);
        legend([h1 h2 h2 h4],{leg_str1,leg_str2,'Medians','Means'},'location',leg_loc);    
        grid

%% calc and save
    case 'calc and save'
        %note that PlotTimeHeightVap3 gives X_mean and Y_mean values (and
        %mid_Xbins and mid_Ybins, X_mode, Y_mode, std_dev_Y and rms_Y

        savefile_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_' gcm_str '_' xlabelstr '_' ylabelstr '_' datestr(now,30) '.mat'];

        %         savefile_prctiles=remove_character(savefile_prctiles,'*','');
        %         savefile_prctiles=remove_character(savefile_prctiles,'\','');
        %         savefile_prctiles=remove_character(savefile_prctiles,' ','_'); %replace all spaces with underscores - latex can handle single
        %         savefile_prctiles=remove_character(savefile_prctiles,'<','.LT.'); %replace all spaces with underscores - latex can handle single
        %         savefile_prctiles=remove_character(savefile_prctiles,'>','.GT.'); %replace all spaces with underscores - latex can handle single
        %         savefile_prctiles=remove_character(savefile_prctiles,'%','pct');
        savefile_prctiles = remove_problem_chars(savefile_prctiles); %this removes the above

        %below is the savefile that stores the latest filename for
        %savefile_prctiles to save copying and pasting the name
        savefilenames_prctiles = [savedir_prctiles 'saved_prctiles_2Dpdf_FILENAME_' gcm_str '_' xlabelstr '_' ylabelstr '.mat'];
        savefilenames_prctiles = remove_problem_chars(savefilenames_prctiles);
        savefile_gcmstr = ['savefile_prctiles_' gcm_str];
        eval([savefile_gcmstr '=savefile_prctiles';]);
        save(savefilenames_prctiles,savefile_gcmstr);

        prcs = [10 50 90];
        [prcY_2Dpdf,imed]=percentiles_from_PDF(Ybins,qh(1:end-1,1:end-1),[10 50 90],1);
        [prcX_2Dpdf,imed]=percentiles_from_PDF(Xbins,qh(1:end-1,1:end-1),[10 50 90],2);


        save(savefile_prctiles,'Xbins');
        save(savefile_prctiles,'Ybins','mid_Xbins','mid_Ybins','X_mean','Y_mean','X_mode','Y_mode','std_dev_X','std_dev_Y','rms_X','rms_Y','prcX_2Dpdf','prcY_2Dpdf','-APPEND');


end


