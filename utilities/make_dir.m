function make_dir(direc)

istr = strfind(direc,'/');
i0=1;
for i=2:length(istr)
    inds = i0:istr(i);
    if exist([direc(inds)])~=7
        eval(['!mkdir ' direc(inds)]); %make it
    end
end