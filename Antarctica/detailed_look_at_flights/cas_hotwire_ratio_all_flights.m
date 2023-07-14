%load all flights and calculate and save the CAS to hotwire LWC ratios
filedirA='Y:\BAS_flights\';



flight_choose=[99 100 101 102 104 105 108 113 117 120 122 123];
%119 and 129 CAS data doesn't seem to load in properly??

data_case = 'LWC';
data_case = 'Number of droplets';

switch data_case
    case 'LWC'


        for ilwc_ratio=1:length(flight_choose)
            flight=flight_choose(ilwc_ratio);
            filedir=[filedirA 'flight' num2str(flight) '\'];

            call_from_cas_hotwire_ratio_all_flights=1; %overwrites filedir in multi_read_plot_CAS
            %Need to set this flag after every run of multi_plot_CAS
            multi_read_plot_CAS

            scatter_plot;  %select 'CAS hotwire matches'
            close(gcf);

            %store data in a structure
            lwc_ratio_dat(ilwc_ratio).flight=flight;
            lwc_ratio_dat(ilwc_ratio).median_ratio=med_vals;
            lwc_ratio_dat(ilwc_ratio).median_ratio2=med_vals2;
            lwc_ratio_dat(ilwc_ratio).median_ratio3=med_vals3;

            lwc_ratio_dat(ilwc_ratio).nvals=nvals;
            lwc_ratio_dat(ilwc_ratio).nvals2=nvals2;
            lwc_ratio_dat(ilwc_ratio).nvals3=nvals3;

        end

        isave=0;
        if isave==1
            save_file='LWC_ratio_matlab_data.mat';
            save([filedirA save_file],'lwc_ratio_dat');
        end


    case 'Number of droplets'
         for ilwc_ratio=2:length(flight_choose)
            flight=flight_choose(ilwc_ratio);
            filedir=[filedirA 'flight' num2str(flight) '\'];

            call_from_cas_hotwire_ratio_all_flights=1; %overwrites filedir in multi_read_plot_CAS
            %Need to set this flag after every run of multi_plot_CAS
            multi_read_plot_CAS

            scatter_plot;  %select 'Number of droplets vs. mode diameter'
            close(gcf);

            %store data in a structure
            lwc_ratio_dat(ilwc_ratio).flight=flight;
            lwc_ratio_dat(ilwc_ratio).median_ratio=med_vals;
            lwc_ratio_dat(ilwc_ratio).median_ratio2=med_vals2;
            lwc_ratio_dat(ilwc_ratio).median_ratio3=med_vals3;
            lwc_ratio_dat(ilwc_ratio).median_ratio4=med_vals4;
            lwc_ratio_dat(ilwc_ratio).median_ratio5=med_vals5;
            lwc_ratio_dat(ilwc_ratio).median_ratio6=med_vals6;
            
            
            lwc_ratio_dat(ilwc_ratio).nvals=nvals;
            lwc_ratio_dat(ilwc_ratio).nvals2=nvals2;
            lwc_ratio_dat(ilwc_ratio).nvals3=nvals3;
            lwc_ratio_dat(ilwc_ratio).nvals4=nvals4;
            lwc_ratio_dat(ilwc_ratio).nvals5=nvals5;
            lwc_ratio_dat(ilwc_ratio).nvals6=nvals6;            

        end

        isave=0;
        if isave==1
            save_file='Number_of_droplets_vs_mode_matlab_data.mat';
            save([filedirA save_file],'lwc_ratio_dat');
        end

        




end

%load([filedirA save_file]);


