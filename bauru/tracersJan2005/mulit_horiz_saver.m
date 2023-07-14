[sx sy st]=size(pressure_hslice(1).dat);

idir=1;

%for supersaturation run timeseries.m first to calculate si array
for it=7:44  %st           
    wrap_slice;
    multisaveplot;
    close(gcf);
end

