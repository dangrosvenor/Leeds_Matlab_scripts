loc_labs='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

for i=1:length(ilat)
    fprintf(1,'\n%c %f',loc_labs(i),pdat(1).p(ilat(i),ilon(i)))
end
