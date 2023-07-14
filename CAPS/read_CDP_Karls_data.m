%read Karl's CDP spreadsheets - after having copied the CDP tab to another spreadsheet and saving
%as .txt

file_choose='Karl Flight 553';
%file_choose=1;

switch file_choose
    case 'Karl Flight 553'
        filedir='Y:\BAS_flights\Karls_CDP_CAS_comparisons\';
        fileName=[filedir 'Flight_553_CDP.txt'];
        
        
end



clear CDP_per_cc  CDP_time_Karl CDP_tot_2to50 CDP_tot_3to50 CDP_tot_5to50 ...
CDP_LWC_3to50 CDP_meanVol CDP_effective_diam


CDP_bins_Karl_mid=[2.5	3.5	4.5	5.5	6.5	7.5	8.5	9.5	10.5	11.5	12.5	13.5	15	17	19	21	23	25	27	29	31	33	35	37	39	41	43	45	47	49];
Dupper=50;
CDP_bins_Karl(31)=Dupper;
for i=length(CDP_bins_Karl_mid):-1:1
    CDP_bins_Karl(i) = 2*CDP_bins_Karl_mid(i) - Dupper;
    Dupper=CDP_bins_Karl(i);
end



fid=fopen(fileName,'rt');

%skip some lines to get to the data
for i=1:8
    textline=fgetl(fid);
end

i=0;
go=1;
while go==1
i=i+1;
 [vals_CDP]=textscan(fid,...   %
        '%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1); 
 
    if size(vals_CDP{1},1)==0
        break %have reached the end
    end

 CDP_time_Karl(i)=vals_CDP{1}; 
 
 if isnan(double(vals_CDP{3}))==1
     CDP_tot_2to50(i)=0;
     CDP_tot_3to50(i)=0;
     CDP_tot_5to50(i)=0;
     CDP_LWC_3to50(i)=0;
     CDP_meanVol(i)=0;
     CDP_effective_diam(i)=0;

     for idat=9:38
         CDP_per_cc(i,idat-8)=0;
     end

 else
 
     CDP_tot_2to50(i)=double(vals_CDP{3});
     CDP_tot_3to50(i)=double(vals_CDP{4});
     CDP_tot_5to50(i)=double(vals_CDP{5});
     CDP_LWC_3to50(i)=double(vals_CDP{6});
     CDP_meanVol(i)=double(vals_CDP{7});
     CDP_effective_diam(i)=double(vals_CDP{8});

     for idat=9:38
         CDP_per_cc(i,idat-8)=double(vals_CDP{idat});
     end
 
 end
 
    
end


disp('Finished reading CDP data');