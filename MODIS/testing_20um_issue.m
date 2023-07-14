
%% No confidence screening, but with phase, sunglint and phase outcome.
%  Also, no screening for re>re_min or max abs difference

SD_id = hdfsd('start',filename_h5,'read'); %open the file
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO);
[re_diff,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius_Difference',SD_id,INFO);
status = hdfsd('end',SD_id);  %close the file


%No confidence screending:-
 ind=find( mask_5km(1,:,:) < 1 | mask_5km(2,:,:) > 0 | qapq_5km(7,:,:) < 1);  
 

 re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
    
% Additional screening for liquid phase, successful outcome and no sunglint   
ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];
    

re_21 = re;
re_21(ihtot)=NaN;
%re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

bins=[0:0.2:30];
nd_16=ndhistc_run(re_16(:),bins);
nd_21=ndhistc_run(re_21(:),bins);
nd_37=ndhistc_run(re_37(:),bins);

figure
plot(bins(1:end-1),nd_16,'k-'); hold on;
plot(bins(1:end-1),nd_21,'b-');
plot(bins(1:end-1),nd_37,'r-');

me_16 = meanNoNan(re_16(:),1);
me_21 = meanNoNan(re_21(:),1);
me_37 = meanNoNan(re_37(:),1);

n_16 = sum(nd_16);
n_21 = sum(nd_21);
n_37 = sum(nd_37);

legend({['re16 N=' num2str(n_16) '  me=' num2str(me_16)],['re21 N=' num2str(n_21) '  me=' num2str(me_21)],['re37 N=' num2str(n_37) '  me=' num2str(me_37)]});

%% As above, but with confidence screening 
SD_id = hdfsd('start',filename_h5,'read'); %open the file
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO);
[re_diff,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius_Difference',SD_id,INFO);
status = hdfsd('end',SD_id);  %close the file


%Confidence screending:-
 ind=find( mask_5km(1,:,:) < 1 | mask_5km(2,:,:) > 0 | qapq_5km(7,:,:) < 3);  
 

 re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
    
% Additional screening for liquid phase, successful outcome and no sunglint   
ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];   
    
    

re_21 = re;
re_21(ihtot)=NaN;
%re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

bins=[0:0.2:30];
nd_16=ndhistc_run(re_16(:),bins);
nd_21=ndhistc_run(re_21(:),bins);
nd_37=ndhistc_run(re_37(:),bins);

figure
plot(bins(1:end-1),nd_16,'k-'); hold on;
plot(bins(1:end-1),nd_21,'b-');
plot(bins(1:end-1),nd_37,'r-');


me_16 = meanNoNan(re_16(:),1);
me_21 = meanNoNan(re_21(:),1);
me_37 = meanNoNan(re_37(:),1);

n_16 = sum(nd_16);
n_21 = sum(nd_21);
n_37 = sum(nd_37);

legend({['re16 N=' num2str(n_16) '  me=' num2str(me_16)],['re21 N=' num2str(n_21) '  me=' num2str(me_21)],['re37 N=' num2str(n_37) '  me=' num2str(me_37)]});



%% As above, but with the re>re_min  screen
SD_id = hdfsd('start',filename_h5,'read'); %open the file
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO);
[re_diff,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius_Difference',SD_id,INFO);
status = hdfsd('end',SD_id);  %close the file


%Confidence screending:-
 ind=find( mask_5km(1,:,:) < 1 | mask_5km(2,:,:) > 0 | qapq_5km(7,:,:) < 3);  
 

 re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
    
% Additional screening for liquid phase, successful outcome and no sunglint   
ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];   
    
    

re_21 = re;
re_21(ihtot)=NaN;
%re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

%Screenings for re<20, re>re_min and re_diff_max
re_min=0.5;
re_diff_max=15;

%ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min | re_16>20);
ihtot2 = find (re_16<re_min | re_16>20);
%ihtot2 = find (re_16>20);
re_16(ihtot2)=NaN;

%ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min | re_37>20);
ihtot2 = find (re_37<re_min | re_37>20);
%ihtot2 = find (re_37>20);
re_37(ihtot2)=NaN;

bins=[0:0.2:30];
nd_16=ndhistc_run(re_16(:),bins);
nd_21=ndhistc_run(re_21(:),bins);
nd_37=ndhistc_run(re_37(:),bins);

figure
plot(bins(1:end-1),nd_16,'k-'); hold on;
plot(bins(1:end-1),nd_21,'b-');
plot(bins(1:end-1),nd_37,'r-');


me_16 = meanNoNan(re_16(:),1);
me_21 = meanNoNan(re_21(:),1);
me_37 = meanNoNan(re_37(:),1);

n_16 = sum(nd_16);
n_21 = sum(nd_21);
n_37 = sum(nd_37);

legend({['re16 N=' num2str(n_16) '  me=' num2str(me_16)],['re21 N=' num2str(n_21) '  me=' num2str(me_21)],['re37 N=' num2str(n_37) '  me=' num2str(me_37)]});
title ('With re>re min');



    
%% As above, but with the re>re_min and max abs of difference
SD_id = hdfsd('start',filename_h5,'read'); %open the file
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO);
[re_diff,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius_Difference',SD_id,INFO);
status = hdfsd('end',SD_id);  %close the file


%Confidence screending:-
 ind=find( mask_5km(1,:,:) < 1 | mask_5km(2,:,:) > 0 | qapq_5km(7,:,:) < 3);  
 

 re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
    
% Additional screening for liquid phase, successful outcome and no sunglint   
ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];   
    
    

re_21 = re;
re_21(ihtot)=NaN;
%re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

%Screenings for re<20, re>re_min and re_diff_max
re_min=0.5;
re_diff_max=15;

ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min | re_16>20);
%ihtot2 = find (re_16<re_min | re_16>20);
%ihtot2 = find (re_16>20);
re_16(ihtot2)=NaN;

ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min | re_37>20);
%ihtot2 = find (re_37<re_min | re_37>20);
%ihtot2 = find (re_37>20);
re_37(ihtot2)=NaN;

bins=[0:0.2:30];
nd_16=ndhistc_run(re_16(:),bins);
nd_21=ndhistc_run(re_21(:),bins);
nd_37=ndhistc_run(re_37(:),bins);

figure
plot(bins(1:end-1),nd_16,'k-'); hold on;
plot(bins(1:end-1),nd_21,'b-');
plot(bins(1:end-1),nd_37,'r-');


me_16 = meanNoNan(re_16(:),1);
me_21 = meanNoNan(re_21(:),1);
me_37 = meanNoNan(re_37(:),1);

n_16 = sum(nd_16);
n_21 = sum(nd_21);
n_37 = sum(nd_37);

legend({['re16 N=' num2str(n_16) '  me=' num2str(me_16)],['re21 N=' num2str(n_21) '  me=' num2str(me_21)],['re37 N=' num2str(n_37) '  me=' num2str(me_37)]});
title ('With re>re min and max abs diff');


%% As above, but using only same points (making others NaN when one is NaN)
SD_id = hdfsd('start',filename_h5,'read'); %open the file
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO);
[re_diff,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius_Difference',SD_id,INFO);
status = hdfsd('end',SD_id);  %close the file


%Confidence screending:-
 ind=find( mask_5km(1,:,:) < 1 | mask_5km(2,:,:) > 0 | qapq_5km(7,:,:) < 3);  
 

 re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
    
% Additional screening for liquid phase, successful outcome and no sunglint   
ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];   
    
    

re_21 = re;
re_21(ihtot)=NaN;
%re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

%Screenings for re<20, re>re_min and re_diff_max
re_min=0.5;
re_diff_max=15;

ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min | re_16>20);
%ihtot2 = find (re_16<re_min | re_16>20);
%ihtot2 = find (re_16>20);
re_16(ihtot2)=NaN;

ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min | re_37>20);
%ihtot2 = find (re_37<re_min | re_37>20);
%ihtot2 = find (re_37>20);
re_37(ihtot2)=NaN;

re_16(isnan(re_21) | isnan(re_37))=NaN;
re_21(isnan(re_16) | isnan(re_37))=NaN;
re_37(isnan(re_16) | isnan(re_21))=NaN;

bins=[0:0.2:30];
nd_16=ndhistc_run(re_16(:),bins);
nd_21=ndhistc_run(re_21(:),bins);
nd_37=ndhistc_run(re_37(:),bins);

figure
plot(bins(1:end-1),nd_16,'k-'); hold on;
plot(bins(1:end-1),nd_21,'b-');
plot(bins(1:end-1),nd_37,'r-');


me_16 = meanNoNan(re_16(:),1);
me_21 = meanNoNan(re_21(:),1);
me_37 = meanNoNan(re_37(:),1);

n_16 = sum(nd_16);
n_21 = sum(nd_21);
n_37 = sum(nd_37);

legend({['re16 N=' num2str(n_16) '  me=' num2str(me_16)],['re21 N=' num2str(n_21) '  me=' num2str(me_21)],['re37 N=' num2str(n_37) '  me=' num2str(me_37)]});
title ('With re>re min and max abs diff, consistent points');



%% Using same points, but with no confidence screening (and no re<20
%% screening)
SD_id = hdfsd('start',filename_h5,'read'); %open the file
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO);
[re_diff,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius_Difference',SD_id,INFO);
status = hdfsd('end',SD_id);  %close the file


% Confidence screending:-
 ind=find( mask_5km(1,:,:) < 1 | mask_5km(2,:,:) > 0 | qapq_5km(7,:,:) < 1);  %No screening for <1
 

 re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
    
% Additional screening for liquid phase, successful outcome and no sunglint   
ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];   
    
    

re_21 = re;
re_21(ihtot)=NaN;
%re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

%Screenings for re<20, re>re_min and re_diff_max
re_min=0.5;
re_diff_max=15;

%ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min | re_16>20);
ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min);
%ihtot2 = find (re_16<re_min | re_16>20);
%ihtot2 = find (re_16>20);
re_16(ihtot2)=NaN;

%ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min | re_37>20);
ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min);
%ihtot2 = find (re_37<re_min | re_37>20);
%ihtot2 = find (re_37>20);
re_37(ihtot2)=NaN;

re_16(isnan(re_21) | isnan(re_37))=NaN;
re_21(isnan(re_16) | isnan(re_37))=NaN;
re_37(isnan(re_16) | isnan(re_21))=NaN;

bins=[0:0.2:30];
nd_16=ndhistc_run(re_16(:),bins);
nd_21=ndhistc_run(re_21(:),bins);
nd_37=ndhistc_run(re_37(:),bins);

figure
plot(bins(1:end-1),nd_16,'k-'); hold on;
plot(bins(1:end-1),nd_21,'b-');
plot(bins(1:end-1),nd_37,'r-');


me_16 = meanNoNan(re_16(:),1);
me_21 = meanNoNan(re_21(:),1);
me_37 = meanNoNan(re_37(:),1);

n_16 = sum(nd_16);
n_21 = sum(nd_21);
n_37 = sum(nd_37);

legend({['re16 N=' num2str(n_16) '  me=' num2str(me_16)],['re21 N=' num2str(n_21) '  me=' num2str(me_21)],['re37 N=' num2str(n_37) '  me=' num2str(me_37)]});
title ('Consistent points, no confidence screening and no re<20um screening');



    
    
    
        
    
    
    
    
    
        
    
    
    
    
    
        
    
    
    
    
    
        
    
    