%this saves the timeseries of the variables in modis_var
%(i.e. the name with _tim added)
fprintf(1,'\n saving L2 data in modis_var array...');

make_mockL3_variables_script_name = 'make_mockL3_variables_Puijo_CDP_save_from_scatter'; %all variables
make_mockL3_variables_script_name = 'make_variables_Irshad_scatter'; %just the ones that Irshad wants (cat I)
%make_mockL3_variables_script_name = 'make_variables_Irshad_scatter_cat2'; %just the ones that Irshad wants

eval(make_mockL3_variables_script_name);


tag = '_no_std_filtering';
tag = '_no_screening';
%tag = '_std37_filtering_screening';
%tag = '_std37_filtering_screening_FINAL';
%tag = '_std37_filtering_screening_FINAL_Irshad_data_only';
%tag = '_std37_filtering_screening_FINAL_Irshad_data_only_cat2';
tag = '0.5x1deg_cat3_CF80';
tag = '5x5km_cat2_CF99_CTT265';
tag = '5x5km_cat1_CF99_CTT273_ceilometer_screening';

savevarname = ['/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/5km_matches_Dan' tag '.mat'];
savevarname_text = ['/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/ground_station_matches' tag '.txt'];
fid=fopen(savevarname_text,'wt');

app_str=')';

save_type = 'all';
save_type = 'scatter_plot values';





switch save_type
    case 'all'
        inotnan = find(isnan(Cloud_Fraction_Liquid.timeseries3)==0);
         fprintf(fid,'NO FILTERING (ALL data)\n');
    case 'scatter_plot values'
        inotnan = find(isnan(Cloud_Fraction_Liquid.timeseries3(ihtot(istd)))==0);
        inotnan = ihtot(istd(inotnan));
        
        fprintf(fid,'SCREENING APPLIED : %s\n',thresh_str);
        
end


fprintf(fid,'%s ','Year');
fprintf(fid,'%s ','Day');
fprintf(fid,'%s ','Month');
fprintf(fid,'%s ','Hour');
fprintf(fid,'%s ','Minute');

eval_str = ['sdat = size(squeeze(' modis_var{1} '.timeseries3(inotnan) ) );'];
eval(eval_str);


for iread_modis=1:length(modis_var)
    
    MODIS_varname = modis_var{iread_modis};
    MODIS_varname2 = remove_character(MODIS_varname,' ','_');
    eval_str = ['fprintf(fid,''%s '',''' MODIS_varname2 ''');']; eval(eval_str);

end

fprintf(fid,'\n');

for i=1:sdat(1)
    
    [Y,M,D,H,MM,S] = datevec(Date_Time_Swath.timeseries3(inotnan(i)));
    
    fprintf(fid,'%f ',Y);
    fprintf(fid,'%f ',D);
    fprintf(fid,'%f ',M);
    fprintf(fid,'%f ',H);
    fprintf(fid,'%f ',MM);


    for iread_modis=1:length(modis_var)

%        fprintf(1,'\n saving %d out of %d',iread_modis,length(modis_var));

        MODIS_varname = modis_var{iread_modis};

        %double quotes tells matlab to insert a ' in the string (without
        %closing the quote for the sring prcessing)


       


        if i==1
            eval_str = ['save_L2_Irshad_func(' MODIS_varname '.timeseries3(inotnan),savevarname,MODIS_varname,app_str);'];
            eval(eval_str);

            app_str=',''-APPEND'')';

        end

        if strcmp(MODIS_varname,'MODIS_swath_filename')
            eval_str = ['fprintf(fid,''%s '',char(' MODIS_varname '.timeseries3{inotnan(i)}) );'];
            eval(eval_str);
        else
            eval_str = ['fprintf(fid,''%f '',' MODIS_varname '.timeseries3(inotnan(i)) );'];
            eval(eval_str);

        end

    end

    fprintf(fid,'\n');



end

fclose(fid);


% save(savevarname,'timLAT','-APPEND');
% save(savevarname,'timLON','-APPEND');
%MLAT = timLAT;
%MLON = timLON;
save(savevarname,'MLAT','-APPEND');
save(savevarname,'MLON','-APPEND');
%save(savevarname,'scantime_matlab_tim','-APPEND');

fprintf(1,'\nDone\n');


