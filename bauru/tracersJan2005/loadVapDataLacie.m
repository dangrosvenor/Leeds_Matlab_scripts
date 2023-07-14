clear direc

%direc(1).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv_25km_3/';
%direc(2).dir='c:/documents and settings/tom/my documents/dan/dmidamp_2/';
direc(1).dir='f:/les/dmidamp_2/';
direc(2).dir='f:/les/damp_ccn960/';
%direc(4).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv/';
direc(3).dir='f:/les/damp_inx10/';
direc(4).dir='f:/les/500mres/';

endt={'65' '88' '65' '88'};
endtALL={'65' '65' '65' '88'};

%saveflag=[0 0 0 0 0 0 0];
loadflag=0;
loadflag=zeros([1 10]);
loadflag([8:10])=1;
%loadflag([3])=1;


clear GridDan
for j=1:length(direc)
	exdir=[direc(j).dir 'results/diags/'];
	exname=[exdir 'gridDan'];
	load(exname);
	GridDan(j)=gridDan;
	
	if loadflag(1)==1
	exname=[exdir 'icediag4_5thSep2005'];
	load(exname);
	icediag4(j)=icediags4;
    clear icediags4;
	end
    
    if loadflag(2)==1 %microphysical process rates
	exname=[exdir 'icediag_5thSep2005'];
	load(exname);
	icediag(j)=icediags;
	end
	
	if loadflag(3)==1
        exname=[exdir 'maxW_1-65'];
        load(exname);
        MaxW(j)=maxW;
        
	end
    
    if loadflag(4)==1
		exname=[exdir 'vap_prctiles'];
		load(exname);
		vap_prc(j)=vap_prctiles;
	end
    
     if loadflag(5)==1
		%exname=[exdir 'DQdehydration_10thSep2005_inc44-86'];
        exname=[exdir 'DQdehydration_10thSep2005_inc44-86'];
        load(exname,'dq_vap','dq_tots');
        dq_tot(j)=dq_tots;
        if exist('dq_vap'); dq_vaps(j)=dq_vap; end
	end
    
    if loadflag(6)==1
		exname=[exdir 'vap+tot_prcs_1-65'];
		load(exname,'tot_prc','vap_prc');
        vap_prctiles(j)=vap_prc;
        tot_prctiles(j)=tot_prc;
	end
    
    if loadflag(7)==1
		exname=[exdir 'nn_10thSep2005_1-86'];
		load(exname,'n');
        nn(j)=n;
	end
    
    if loadflag(8)==1
		exname=[exdir 'icediagsALL_1-' endtALL{j}];
		load(exname,'icediagALL');
        icediagsALL(j)=icediagALL;
	end
    
    if loadflag(9)==1
		exname=[exdir 'icediag_1-' endt{j}];
		load(exname,'icediags');
        icediag(j)=icediags;
	end
    
    if loadflag(10)==1
		exname=[exdir 'fluxes_1-65'];
        load(exname,'flux');
		fluxes(j)=flux;
	end



end

'loadVapData finished'

