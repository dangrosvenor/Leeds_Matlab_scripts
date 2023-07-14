imdp=1;
modis_data_plot_strs{imdp}='Number of droplets cell values time mean - specific days'; mdp_str{imdp}='Nd'; imdp=imdp+1;     
modis_data_plot_strs{imdp}='LWP cell values time mean (grid-box mean) - specific days'; mdp_str{imdp}='LWP'; imdp=imdp+1;     
modis_data_plot_strs{imdp}='Cloud depth cell values time mean - specific days'; mdp_str{imdp}='H'; imdp=imdp+1;     

clear modis_var
     isave=1;
     modis_var{isave}='Ndays_multi'; isave=isave+1;
     modis_var{isave}='Ndata_multi'; isave=isave+1;     
                        
     for imdp2=1:imdp-1
         eval_str = ['modis_var{isave}=''' mdp_str{imdp2} '_multi'';']; eval(eval_str); isave=isave+1;
         eval_str = ['modis_var{isave}=''' mdp_str{imdp2} '_sq_multi'';']; eval(eval_str); isave=isave+1;
         eval_str = ['modis_var{isave}=''' mdp_str{imdp2} '_PDF_multi'';']; eval(eval_str); isave=isave+1;
     end
     
     
     %ignore the above if clearing here
%     clear modis_var; isave=1; 
%     eval_str = ['modis_var{isave}=''' mdp_str{imdp2} '_PDF_multi'';']; eval(eval_str); isave=isave+1;
     
     
     '';