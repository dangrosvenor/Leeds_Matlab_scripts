%% Start


%low low = lowCF, low sensor

%% Tau

%relative increases in mean tau
rel_incs_low_low([2])   = 100*(mean(Tau_all_pdfs{12}(1).dat) / mean(Tau_all_pdfs{1}(1).dat) - 1); 
rel_incs_low_high([2])   = 100*(mean(Tau_all_pdfs{12}(2).dat) / mean(Tau_all_pdfs{1}(2).dat) - 1); 
rel_incs_high_low([2])   = 100*(mean(Tau_all_pdfs{12}(3).dat) / mean(Tau_all_pdfs{1}(3).dat) - 1); 
rel_incs_high_high([2])   = 100*(mean(Tau_all_pdfs{12}(4).dat) / mean(Tau_all_pdfs{1}(4).dat) - 1);

%% Re
%; rel_incs([3 1 5]);
rel_incs_low_low([4 7 10])   = [100*(mean(Re_all_pdfs_lowCF{12}(3).dat) / mean(Re_all_pdfs_lowCF{1}(3).dat) - 1),... 
    100*(mean(Re_all_pdfs_lowCF{12}(1).dat) / mean(Re_all_pdfs_lowCF{1}(1).dat) - 1),...
    100*(mean(Re_all_pdfs_lowCF{12}(5).dat) / mean(Re_all_pdfs_lowCF{1}(5).dat) - 1),...
                                ];

%rel_incs_low_high([4 7 10])   = rel_incs([4 2 6]);
rel_incs_low_high([4 7 10])   = [100*(mean(Re_all_pdfs_lowCF{12}(4).dat) / mean(Re_all_pdfs_lowCF{1}(4).dat) - 1)...
                                100*(mean(Re_all_pdfs_lowCF{12}(2).dat) / mean(Re_all_pdfs_lowCF{1}(2).dat) - 1)...
                                100*(mean(Re_all_pdfs_lowCF{12}(6).dat) / mean(Re_all_pdfs_lowCF{1}(6).dat) - 1)...
                                ];


% rel_incs_high_low([4 7 10])   = rel_incs([3 1 5]);
rel_incs_high_low([4 7 10])   = [100*(mean(Re_all_pdfs_highCF{12}(3).dat) / mean(Re_all_pdfs_highCF{1}(3).dat) - 1)...
                                100*(mean(Re_all_pdfs_highCF{12}(1).dat) / mean(Re_all_pdfs_highCF{1}(1).dat) - 1)...
                                100*(mean(Re_all_pdfs_highCF{12}(5).dat) / mean(Re_all_pdfs_highCF{1}(5).dat) - 1)...
                                ];
 
 
%rel_incs_high_high([4 7 10])   = rel_incs([4 2 6]);
rel_incs_high_high([4 7 10])   = [100*(mean(Re_all_pdfs_highCF{12}(4).dat) / mean(Re_all_pdfs_highCF{1}(4).dat) - 1)...
                                100*(mean(Re_all_pdfs_highCF{12}(2).dat) / mean(Re_all_pdfs_highCF{1}(2).dat) - 1)...
                                100*(mean(Re_all_pdfs_highCF{12}(6).dat) / mean(Re_all_pdfs_highCF{1}(6).dat) - 1)...
                                ];
                            
%% Nd                            
%rel_incs_low_low([1 6 9])   = rel_incs([3 1 5]);  %Nd
rel_incs_low_low([1 6 9])   = [100*(mean(Nd_all_pdfs_lowCF{12}(3).dat) / mean(Nd_all_pdfs_lowCF{1}(3).dat) - 1)...
                                100*(mean(Nd_all_pdfs_lowCF{12}(1).dat) / mean(Nd_all_pdfs_lowCF{1}(1).dat) - 1)...
                                100*(mean(Nd_all_pdfs_lowCF{12}(5).dat) / mean(Nd_all_pdfs_lowCF{1}(5).dat) - 1)...
                                ];


%rel_incs_low_high([1 6 9])   = rel_incs([4 2 6]);
rel_incs_low_high([1 6 9])   = [100*(mean(Nd_all_pdfs_lowCF{12}(4).dat) / mean(Nd_all_pdfs_lowCF{1}(4).dat) - 1)...
                                100*(mean(Nd_all_pdfs_lowCF{12}(2).dat) / mean(Nd_all_pdfs_lowCF{1}(2).dat) - 1)...
                                100*(mean(Nd_all_pdfs_lowCF{12}(6).dat) / mean(Nd_all_pdfs_lowCF{1}(6).dat) - 1)...
                                ];
                            
                            
%rel_incs_high_low([1 6 9])   = rel_incs([3 1 5]);
rel_incs_high_low([1 6 9])   = [100*(mean(Nd_all_pdfs_highCF{12}(3).dat) / mean(Nd_all_pdfs_highCF{1}(3).dat) - 1)...
                                100*(mean(Nd_all_pdfs_highCF{12}(1).dat) / mean(Nd_all_pdfs_highCF{1}(1).dat) - 1)...
                                100*(mean(Nd_all_pdfs_highCF{12}(5).dat) / mean(Nd_all_pdfs_highCF{1}(5).dat) - 1)...
                                ];
                            
%rel_incs_high_high([1 6 9])   = rel_incs([4 2 6]); 
rel_incs_high_high([1 6 9])   = [100*(mean(Nd_all_pdfs_highCF{12}(4).dat) / mean(Nd_all_pdfs_highCF{1}(4).dat) - 1) ...
                                100*(mean(Nd_all_pdfs_highCF{12}(2).dat) / mean(Nd_all_pdfs_highCF{1}(2).dat) - 1) ...
                                100*(mean(Nd_all_pdfs_highCF{12}(6).dat) / mean(Nd_all_pdfs_highCF{1}(6).dat) - 1) ...
                                ]; 
                            
                            
%% Increases due to tau and Re increases


%% dN due to dtau (keeping Re constant)
       
 itau=1; %low CF, low sensor
 idat=3; %2.1 um
 low_or_highCF='lowCF';
 make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
 N12 = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
 rel_incs_low_low([3]) = (N12 - N1) ./ (N2 - N1);
 
 itau=2; %low CF, high sensor
 idat=4; %2.1 um
 low_or_highCF='lowCF';
 make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
 N12 = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
 rel_incs_low_high([3]) = (N12 - N1) ./ (N2 - N1);
 
 itau=3; %high CF, low sensor
 idat=3; %2.1 um
 low_or_highCF='highCF';
 make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
 N12 = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
 rel_incs_high_low([3]) = (N12 - N1) ./ (N2 - N1);
 
 itau=4; %high CF, high sensor
 idat=4; %2.1 um
 low_or_highCF='highCF';
 make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
 N12 = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
 rel_incs_high_high([3]) = (N12 - N1) ./ (N2 - N1);
                    
                    
%% dN due to dRe_2.1 (keeping Tau constant)
itau=1; %low CF, low sensor
idat=3; %2.1 um
low_or_highCF='lowCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_low_low([5]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_low_low([12]) = (N12_tau - N1) ./ (N2 - N1);

itau=2; %low CF, high sensor
idat=4; %2.1 um
low_or_highCF='lowCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_low_high([5]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_low_high([12]) = (N12_tau - N1) ./ (N2 - N1);


itau=3; %high CF, low sensor
idat=3; %2.1 um
low_or_highCF='highCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_high_low([5]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_high_low([12]) = (N12_tau - N1) ./ (N2 - N1);


itau=4; %high CF, high sensor
idat=4; %2.1 um
low_or_highCF='highCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_high_high([5]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_high_high([12]) = (N12_tau - N1) ./ (N2 - N1);
 


%% dN due to dRe_1.6 (keeping Tau constant)
itau=1; %low CF, low sensor
idat=1; %1.6 um
low_or_highCF='lowCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_low_low([8]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_low_low([13]) = (N12_tau - N1) ./ (N2 - N1);

itau=2; %low CF, high sensor
idat=2; %1.6 um
low_or_highCF='lowCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_low_high([8]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_low_high([13]) = (N12_tau - N1) ./ (N2 - N1);

itau=3; %high CF, low sensor
idat=1; %1.6 um
low_or_highCF='highCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_high_low([8]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_high_low([13]) = (N12_tau - N1) ./ (N2 - N1);

itau=4; %high CF, high sensor
idat=2; %1.6 um
low_or_highCF='highCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_high_high([8]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_high_high([13]) = (N12_tau - N1) ./ (N2 - N1);



       
       
%% dN due to dRe_3.7 (keeping Tau constant)
itau=1; %low CF, low sensor
idat=5; %3.7 um
low_or_highCF='lowCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_low_low([11]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_low_low([14]) = (N12_tau - N1) ./ (N2 - N1);

itau=2; %low CF, high sensor
idat=6; %3.7 um
low_or_highCF='lowCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_low_high([11]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_low_high([14]) = (N12_tau - N1) ./ (N2 - N1);

itau=3; %high CF, low sensor
idat=5; %3.7 um
low_or_highCF='highCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_high_low([11]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_high_low([14]) = (N12_tau - N1) ./ (N2 - N1);

itau=4; %high CF, high sensor
idat=6; %3.7 um
low_or_highCF='highCF';
make_sza_rel_increases_table_vals_subfunc %sub function/script - calcs N1, N2, Re_vals, Tau_vals and scale_fac
N12 = Re_vals(end) .* Tau_vals(1) .* scale_fac(end) ;
rel_incs_high_high([11]) = (N12 - N1) ./ (N2 - N1);
N12_tau = Re_vals(1) .* Tau_vals(end) .* scale_fac(end) ;
rel_incs_high_high([14]) = (N12_tau - N1) ./ (N2 - N1);


%% print the table to file in Latex format
% fid = fopen('/home/disk/eos1/d.grosvenor/matlab/work/MODIS/Nd_attr_table.txt','wt');
% 
% fprintf(fid,'low CF, low sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_low_low(1),rel_incs_low_low(2),rel_incs_low_low(3),rel_incs_low_low(4),rel_incs_low_low(5),rel_incs_low_low(6),rel_incs_low_low(7),rel_incs_low_low(8),rel_incs_low_low(9),rel_incs_low_low(10),rel_incs_low_low(11));
% fprintf(fid,'low CF, high sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_low_high(1),rel_incs_low_high(2),rel_incs_low_high(3),rel_incs_low_high(4),rel_incs_low_high(5),rel_incs_low_high(6),rel_incs_low_high(7),rel_incs_low_high(8),rel_incs_low_high(9),rel_incs_low_high(10),rel_incs_low_high(11));
% fprintf(fid,'high CF, low sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_high_low(1),rel_incs_high_low(2),rel_incs_high_low(3),rel_incs_high_low(4),rel_incs_high_low(5),rel_incs_high_low(6),rel_incs_high_low(7),rel_incs_high_low(8),rel_incs_high_low(9),rel_incs_high_low(10),rel_incs_high_low(11));
% fprintf(fid,'high CF, high sensor & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f & %3.1f\\%% & %3.1f\\%% & %1.2f \\\\ \n',rel_incs_high_high(1),rel_incs_high_high(2),rel_incs_high_high(3),rel_incs_high_high(4),rel_incs_high_high(5),rel_incs_high_high(6),rel_incs_high_high(7),rel_incs_high_high(8),rel_incs_high_high(9),rel_incs_high_high(10),rel_incs_high_high(11));
% 
% fclose(fid);


fid = fopen('/home/disk/eos1/d.grosvenor/matlab/work/MODIS/Nd_attr_table.txt','wt');

fprintf(fid,'low CF, low sensor & %3.1f\\%% & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f \\\\ \n',rel_incs_low_low(2),rel_incs_low_low(1),rel_incs_low_low(4),rel_incs_low_low(5),rel_incs_low_low(12),rel_incs_low_low(6),rel_incs_low_low(7),rel_incs_low_low(8),rel_incs_low_low(13),rel_incs_low_low(9),rel_incs_low_low(10),rel_incs_low_low(11),rel_incs_low_low(14));
fprintf(fid,'low CF, high sensor & %3.1f\\%% & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f \\\\ \n',rel_incs_low_high(2),rel_incs_low_high(1),rel_incs_low_high(4),rel_incs_low_high(5),rel_incs_low_high(12),rel_incs_low_high(6),rel_incs_low_high(7),rel_incs_low_high(8),rel_incs_low_high(13),rel_incs_low_high(9),rel_incs_low_high(10),rel_incs_low_high(11),rel_incs_low_high(14));
fprintf(fid,'high CF, low sensor & %3.1f\\%% & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f \\\\ \n',rel_incs_high_low(2),rel_incs_high_low(1),rel_incs_high_low(4),rel_incs_high_low(5),rel_incs_high_low(12),rel_incs_high_low(6),rel_incs_high_low(7),rel_incs_high_low(8),rel_incs_high_low(13),rel_incs_high_low(9),rel_incs_high_low(10),rel_incs_high_low(11),rel_incs_high_low(14));
fprintf(fid,'high CF, high sensor & %3.1f\\%% & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f & & %3.1f\\%% & %3.1f\\%% & %1.2f & %1.2f \\\\ \n',rel_incs_high_high(2),rel_incs_high_high(1),rel_incs_high_high(4),rel_incs_high_high(5),rel_incs_high_high(12),rel_incs_high_high(6),rel_incs_high_high(7),rel_incs_high_high(8),rel_incs_high_high(13),rel_incs_high_high(9),rel_incs_high_high(10),rel_incs_high_high(11),rel_incs_high_high(14));
                                                       
fclose(fid);
                            
                            
                            
                            
                            