ttt=[38:83];
exdir='g:/runs/sdlavap/results/MPCcomp/iceNC_15-22_gwave/'

for ianim=1:length(ttt)
    tt=ttt(ianim);
    plotTimeHeightVap2;
    set(gcf,'Paperpositionmode','auto');
    exname=[exdir 'iceMR_' time2d '.emf'];
    print(gcf,'-dmeta',exname);
    close(gcf);
end