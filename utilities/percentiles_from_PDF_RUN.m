Nd_bins_mon = [0:25:200 300:100:1400 1500:500:5000];  %29 bins

Nd_bins_mon2 = 0.5 * (Nd_bins_mon(1:end-1) + Nd_bins_mon(2:end) );

prcs=[5 10 20 30 40 50 60 70 80 90 95]; 

[prcs_out,imed] = percentiles_from_PDF(Nd_bins_mon,Nd_PDF_multi,prcs,6);

%calculate the percentile values from the PDF of Nd
