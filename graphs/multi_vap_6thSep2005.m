dire='c:/documents and settings/tom/my documents/dan/';

direc(1).dir='dmi1715_5ppmv_25km_3';
direc(2).dir='dmidamp_2';
direc(3).dir='damp_ccn960';
direc(4).dir='dmi1715_5ppmv';
direc(5).dir='damp_inx10';

plotcase=2;

switch plotcase
case 1
	for idir=2:2
        exname=[dire 'Temp_start_' direc(idir).dir '.emf'];
        watervapourmay2005;
        set(gcf,'paperpositionmode','auto');
        print(gcf,exname,'-dmeta');
        close(gcf);
	end

case 2
	for idir=1:3
        %exname=[dire 'MicroRateTimeHeight_3_' direc(idir).dir '.emf'];
        %exname=[dire 'MicroRate_5_' direc(idir).dir '.emf'];
        exname=[dire 'dq_dehyd_3.67_2_' direc(idir).dir '.emf'];
        plottimeheightvap2;
        set(gcf,'paperpositionmode','auto');
        print(gcf,exname,'-dmeta');
        close(gcf);
	end

end