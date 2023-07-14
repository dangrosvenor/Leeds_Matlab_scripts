LAT_start=-90.5;  %old
LAT_start=-92.5; %SH
LAT_end=0.5;

%NH
LAT_start=0.5;
LAT_end=92.5;

LAT_jump=5;
LAT_jump=1;

LAT_width=0; %gives 1 degree actual widths
LAT_width=1; %gives 2 degree actual widths
LAT_width=2; %gives 3 degree widths

LAT_vals_sza = [LAT_start:LAT_jump:LAT_end];

idays_sza=[1:365];
%idays_sza=[1:182];
%idays_sza=[183:365];

%nstart = 355; %21st Dec 2005

clear SZA_sza grad_sza mid_sza LAT_sza med_Nd med_grad

%create arrays first to save memory - and will tell us if there is enough
%memory early on
Nd_sza = NaN * ones([201 93 length(idays_sza)]);
Np_sza = NaN * ones([201 93 length(idays_sza)]);

for iday_sza=1:length(idays_sza)
    
    fprintf(1,'\niday_sza=%d of %d',iday_sza,length(idays_sza));

    nstart=idays_sza(iday_sza);


    for ilat_sza=1:length(LAT_vals_sza)

        ioverride_pdf=1;

        minXbins=0;
        maxXbins=3000;

        minYbins=0;
        maxYbins=90;

        %select the number of PDF bins
        nXpdf=600;
        nYpdf=25;
        nYpdf=35;
        nYpdf=200;
        %    nYpdf=15;



        LAT_val = [LAT_vals_sza(ilat_sza) LAT_vals_sza(ilat_sza)+LAT_width];
        plotTimeHeightVap3
        LAT_sza(:,ilat_sza)=[LAT(ilat(1)) LAT(ilat(end))]; %the actual latitudes used
        %these refer to the latitudes as stored in the D3 grids, which go
        %from -89.5 to 89.5 and so represent the mid points of -90 to -89,
        %etc.
        man_choose_water_graph=1;
        graph=96;
        time_highlight_path=[];

        iytick_relabel=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
        y_axis_type=''; %default
        x_axis_type='';
        i_set_dateticks=0;
        iadd_nums_above=0;

        xlims=0;
        fsize=14;

        idatetick=0; %flag to put times as HH:MM instead of decimal time
        noplot=1;
        waterVapourMay2005  %                        

        Nd_sza(:,ilat_sza,iday_sza)=xdat(1).x;
%        SZA_sza(:,ilat_sza,iday_sza)=ydat(1).y;
%        grad_sza(:,ilat_sza,iday_sza) = diff(xdat(1).x)./diff(ydat(1).y');
%       
        Np_sza(:,ilat_sza,iday_sza)=sum(qh,2); %number of points that go into the average for each SZA

        close all

    end

%     for isza=1:size(mid_sza,1)
% 
%         med_dat = grad_sza(isza,:,iday_sza);
%         med_dat(isnan(med_dat))='';
%         med_grad(isza) = median(med_dat);
%         if length(med_dat)==0
%             med_grad(isza)=NaN;
%         end
% 
%     end
% 
% 
%     for isza=1:size(Nd_sza,1)
%         med_dat = Nd_sza(isza,:,iday_sza);
%         med_dat(isnan(med_dat))='';
%         med_Nd(isza) = median(med_dat);
%         if length(med_dat)==0
%             med_Nd(isza)=NaN;
%         end
%     end
% 



end

mid_sza = 0.5*(ydat(1).y(1:end-1)+ydat(1).y(2:end));
%LAT=0.5*(LAT_sza(1,:)+LAT_sza(2,:));

clear save_vars
isv=1;
save_vars{isv}='Nd_sza'; isv=isv+1;
%save_vars{isv}='SZA_sza'; isv=isv+1;
%save_vars{isv}='grad_sza'; isv=isv+1;
save_vars{isv}='mid_sza'; isv=isv+1;
save_vars{isv}='Np_sza'; isv=isv+1;
%save_vars{isv}='med_dat'; isv=isv+1;
save_vars{isv}='LAT_sza'; isv=isv+1;


savedir_var='~/modis_work/saved_data_L3/'
savevarname = [savedir_var 'sza_multiple_lats_data_2006_days' num2str(idays_sza(1)) '-' num2str(idays_sza(end)) '_' datestr(now,30)];

fprintf(1,'\n saving Nd_timeseries etc.');

for isv=1:length(save_vars)
    if isv==1
        eval(['save(''' savevarname ''',''' save_vars{isv} ''')']);
    else
        eval(['save(''' savevarname ''',''' save_vars{isv} ''',''-APPEND'')']);
    end
end


