ih_start=findheight((GridDan(1).Z+620)/1000 , 15);
ih_end=findheight((GridDan(1).Z+620)/1000 , 19);
for ih=ih_start:4:ih_end
    
    wrap_slice;
    multisaveplot;
end
    