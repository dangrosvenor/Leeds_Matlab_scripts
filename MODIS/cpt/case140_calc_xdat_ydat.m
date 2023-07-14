
%code snippet for case 140 of watervap


xdat(idat).x=eval(['midXbins_' pdf_str '_' data_str ';']); %
%ydat(idat).y=eval(['PDF_' pdf_str '_' data_str './sum(PDF_' pdf_str '_' data_str ');']); %

%this data should already be bin normalised and divided by the total N - so
%sum(PDF * bin_widths) = 1
ydat(idat).y=eval(['PDF_' pdf_str '_' data_str ';']); %
mean_pdf(idat) = sum(ydat(idat).y.*xdat(idat).x);

        ylab='Bin Normalised Frequency';

switch pdf_plot_type
    case 'Cumulative'
        bins = eval(['Xbins_' pdf_str '_' data_str ';']);
        bin_widths = diff(bins);
        ytemp = cumsum( ydat(idat).y .* bin_widths );
        ydat(idat).y = cat(1,[0 ytemp]);
        %bin edges
        xdat(idat).x = bins;
        
        ylab='Cumulative Frequency';
                
                
end