clear direc

comp='lacieLap';
%comp='uni';

ndir=20;

loadflag=0;
loadflag=zeros([1 12]); 
%loadflag([11])=1;
%loadflag([3])=1;

%loadflag([5])=1; %diags to load

% endt{1}='58'; %direc(1).dir='h:/5ppmv_qnums2/';
% endtALL{1}='58';



switch comp
	case 'uni'
        topdir(1).dir='c:/cygwin/home/user/runs/';
        topdir(2).dir='c:/cygwin/home/user/runs/';
        topdir(3).dir='c:/cygwin/home/user/runs/';
        topdir(4).dir='g:/runs/';
        topdir(5).dir='h:/';
        topdir(6).dir='g:/runs/';
        topdir(7).dir='h:/';
        topdir(8).dir='g:/runs/';
        topdir(9).dir='';
        topdir(10).dir='';
        topdir(11).dir='';
        topdir(12).dir='';
        topdir(13).dir='h:/';
        topdir(14).dir='h:/';
        topdir(15).dir='g:/runs/';
        topdir(16).dir='c:/cygwin/home/user/runs/';
        topdir(17).dir='c:/cygwin/home/user/runs/';
        topdir(18).dir='c:/cygwin/home/user/runs/';
        topdir(19).dir='c:/cygwin/home/user/runs/';
        topdir(20).dir='c:/cygwin/home/user/runs/';
        
    case 'lacieLap'
        for i=1:ndir
            topdir(i).dir='e:/les/';
        end
        
end


loadselect=[1:4]; %files to load

        %direc(1).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv_25km_3/';
%direc(2).dir='c:/documents and settings/tom/my documents/dan/dmidamp_2/';
direc(1).dir=[topdir(1).dir 'dmidamp_2/']; %zdmp=22.858km
direc(2).dir=[topdir(2).dir 'damp_ccn960/']; %zdmp=22.858km  
%direc(4).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv/';
direc(3).dir=[topdir(3).dir 'damp_inx10/']; %zdmp=22.858km
direc(4).dir=[topdir(4).dir '500mres/']; %zdmp=22.858km

direc(5).dir=[topdir(5).dir '5ppmv_qnums2/']; %xc for the hotpool =0 so think didn't run properly **zdmp=27.946km

%direc(1).dir='c:/cygwin/home/user/runs/5ppmv/';
direc(6).dir=[topdir(6).dir '5ppmv/']; %no icediag_nums   **zdmp=27.946km

direc(7).dir=[topdir(7).dir '5ppmv_1000km_3/']; %has icediag_nums  **zdmp=27.946km

direc(8).dir=[topdir(8).dir '5ppmv_ccn960_2/']; %Bezier runs  **zdmp=27.946km
direc(9).dir=[topdir(9).dir '5ppmv_indiv10/']; %  **zdmp=27.946km
direc(10).dir=[topdir(10).dir '5ppmv_th26.8_2/']; % **zdmp=27.946km
direc(11).dir=[topdir(11).dir '4ppmv_th26.8_2/']; % **zdmp=27.946km
 
direc(12).dir=[topdir(12).dir '4ppmv/']; % zdmp=27.946km 

direc(13).dir=[topdir(13).dir 'lowE/']; % **zdmp=27.946km

direc(14).dir=[topdir(14).dir '250mres_1000km/']; % **zdmp=22.858km
direc(15).dir=[topdir(15).dir '3dNewton/']; 
direc(16).dir=[topdir(16).dir '500mNewt/']; %1000km 
direc(17).dir=[topdir(17).dir '1kmNewt/']; %1000km
direc(18).dir=[topdir(18).dir '250m13.4Newt/']; %1000km

direc(19).dir=[topdir(19).dir '13.02.12utc/']; %13th Feb case study - started at 15UTC (not 12utc but 12 LT). 250km domain but warmpool put at 
                                                %250km so will only get half of it in domain as was using uncorrected code
direc(20).dir=[topdir(19).dir '13.02_2/']; %as above but with warmpool at x=10km so all of it will be within domain - both done on PC

machineSTORE=[1 1 1 1 1 1 1 2 2 2 2 1 1 3 3 3 3 3 1 1]; %array containing the machine code for the file list above - 1=PC, 2=Bezier, 3+Newton



endt={'_1-65' '88' '65' '88' '58' '88' '55' '60' '' '' '' '88' '' '' };
endtALL={'_1-65' '_1-65' '_1-65' '_1-88' '_1-58' '_1-88' '_1-55' '_1-60' '' '' '' '_1-88' '' '' '' '' '' '' '' ''};
endtWmax={'88' '88' '88' '88' '88' 'xx' '55' '60' '' '' '' '88' '' ''};


%loadselect=[8];
loadflag([6 5 8])=1;

gridflag=1; %flag to say whether to load in Grid
clear GridDan
for jdir=1:length(loadselect)   %length(direc)
    
    j=loadselect(jdir);
    direcDan(jdir).dir=direc(j).dir;
    machine(jdir)=machineSTORE(j);
    
    exdir=[direc(j).dir 'results/diags/'];
    
    if gridflag==1
	exname=[exdir 'gridDan.mat'];
	load(exname);
	GridDan(jdir)=gridDan;
    end
	
	if loadflag(1)==1
	exname=[exdir 'icediag4_5thSep2005'];
	load(exname);
	icediag4(jdir)=icediags4;
    clear icediags4;
	end
    
    if loadflag(2)==1 %microphysical process rates
	exname=[exdir 'icediag_5thSep2005'];
	load(exname);
	icediag(jdir)=icediags;
	end
	
	if loadflag(3)==1
        exname=[exdir 'maxW_1-' endtWmax{j}];
        load(exname);
        MaxW(jdir)=maxW;
        
	end
    
    if loadflag(4)==1
		exname=[exdir 'vap_prctiles'];
		load(exname);
		vap_prc(jdir)=vap_prctiles;
	end
    
     if loadflag(5)==1
		%exname=[exdir 'DQdehydration_10thSep2005_inc44-86'];
        %exname=[exdir 'dq_dehyd_1-' endtALL{j}];
        exname=[exdir 'dq_dehyd' endtALL{j}];
        load(exname,'dq_vap','dq_tots');
        dq_tot(jdir)=dq_tots;
        if exist('dq_vap'); dq_vaps(jdir)=dq_vap; end
	end
    
    if loadflag(6)==1
		exname=[exdir 'vap+tot_prcs' endtALL{j}];
		load(exname,'tot_prc','vap_prc');
        vap_prctiles(jdir)=vap_prc;
        tot_prctiles(jdir)=tot_prc;
	end
    
    if loadflag(7)==1
		exname=[exdir 'nn_10thSep2005_1-86'];
		load(exname,'n');
        nn(jdir)=n;
	end
    
    if loadflag(8)==1
		exname=[exdir 'icediagsALL' endtALL{j} '.mat'];
		load(exname,'icediagALL');
        icediagsALL(jdir)=icediagALL;
	end
    
    if loadflag(9)==1
		exname=[exdir 'icediag_1-' endt{j}];
		load(exname,'icediags');
        icediag(jdir)=icediags;
	end
    
    if loadflag(10)==1
		exname=[exdir 'fluxes_1-65'];
        load(exname,'flux');
		fluxes(jdir)=flux;
	end
    
    if loadflag(11)==1
		exname=[exdir 'icediag_nums_1-' endt{j}];
		load(exname,'icediags_nums');
        icediag_nums(jdir)=icediags_nums;
	end
    
%     if loadflag(12)==1
%         exname=[exdir 'W_prctiles'];
%         load(exname,'w');
%         w_prctiles(jdir)=w;	
% 	end
    
    if loadflag(12)==1
        %exname=[exdir 'maxW_1-' endt{j}];
        clear wd wb w
        exname=[exdir 'W_prc+dist'];
        load(exname,'w','wb','wd');
        
        w_prctiles(jdir)=w;
        wdist(jdir)=wd;
        wbb(jdir)=wb;
        
        
        
    end


end

'loadVapData finished'

