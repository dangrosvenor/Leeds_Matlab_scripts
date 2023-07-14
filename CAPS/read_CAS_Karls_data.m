%read Karl's CAS spreadsheets - after having copied the CAS tab to another spreadsheet and saving
%as .txt

file_choose='Karl Flight 553';

switch file_choose
    case 'Karl Flight 553'
        filedir='Y:\BAS_flights\Karls_CDP_CAS_comparisons\';
        fileName=[filedir 'Flight_553_CAS.txt'];
        
        
end

clear CAS_bins_Karl
CAS_bins_Karl_mid=[0.575	0.645	0.715	0.785	0.855	0.925	0.9950	1.0650	1.135	1.21	1.375	1.75	2.25	2.75	3.25	3.75	4.5	5.75	6.85	7.55	9.05	11.35	13.75	17.5	22.5	27.5	32.5	37.5	42.5	47.5];

Dupper=50;
CAS_bins_Karl(30) = Dupper;
for i=length(CAS_bins_Karl_mid):-1:2
    CAS_bins_Karl(i-1) = 2*CAS_bins_Karl_mid(i) - Dupper;
    Dupper=CAS_bins_Karl(i-1);
end
%just checking that we had the same size bins as before - we do! Karl
%was giving the mid-points in his spreadsheet.




clear CAS_per_cc  CAS_time_Karl CAS_tot_0_06to50 CAS_tot_5to50 CAS_LWC_3to50 ...
CAS_meanVolDiam CAS_effective_diam
 

fid=fopen(fileName,'rt');

%skip some lines to get to the data
for i=1:8
    textline=fgetl(fid);
end

i=0;
go=1;
while go==1
i=i+1;
 [vals_CAS]=textscan(fid,...   %
        '%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1); 
 
    if size(vals_CAS{1},1)==0
        break %have reached the end
    end

 CAS_time_Karl(i)=vals_CAS{1}; 
     
     
 if isnan(vals_CAS{3})==1
     CAS_tot_0_06to50(i)=0;
     CAS_tot_5to50(i)=0;
     CAS_LWC_3to50(i)=0;
     CAS_meanVolDiam(i)=0;
     CAS_effective_diam(i)=0;

     for idat=8:37
         CAS_per_cc(i,idat-7)=0;
     end
 else

     CAS_tot_0_06to50(i)=double(vals_CAS{3});
     CAS_tot_5to50(i)=double(vals_CAS{4});
     CAS_LWC_3to50(i)=double(vals_CAS{5});
     CAS_meanVolDiam(i)=double(vals_CAS{6});
     CAS_effective_diam(i)=double(vals_CAS{7});

     for idat=8:37
         CAS_per_cc(i,idat-7)=double(vals_CAS{idat});
     end
 
 
 end
 
    
end


disp('Finished reading CAS data');