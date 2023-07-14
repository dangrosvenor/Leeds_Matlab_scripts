saveflag=1;

clims=[];

clear cat

cat{1}='Cloudwater';
cat{2}='Condensation Freezing';
cat{3}='Contact Nucleation';
cat{4}='Homogeneous Nucleation';
cat{5}='Splinters 1';
cat{6}='Frozen Rain 1';
cat{7}='Splinters 2';
cat{8}='RF Splinters 1';
cat{9}='Frozen Rain 2';
cat{10}='Splinters 3';
cat{11}='RF Splinters 2';
%cat{12}='Rainwater';

cols=[1 4:13]; %rainwater is not output in iwc_up file but in rwc_up


zemm=hemm(:,1);

files={'naero' 'ssat' 'lwc' 'ncw' 'rwc' 'nr'};

comp='pc';
comp='laptop';

switch comp
case 'laptop'
	emmdir=['c:/documents and settings/g/my documents/'];
case 'pc'
    emmdir=['c:/documents and settings/login/my documents/'];
end

typ='mr';
typ='nc';
%typ='other';
%typ='w';

switch typ
	case 'mr'
		savedir=[emmdir 'emm_0012/' outdir '/pics/iwc/'];
		
		for i=1:length(cat)
            plotemmtimHf(temm,zemm,iwczt,cols(i),savedir,[cat{i} '_Z'],saveflag);
		end
	
	case 'nc'
		savedir=[emmdir '/emm_0012/' outdir '/pics/inc/'];
		
		for i=1:length(cat)
            plotemmtimHf(temm,zemm,inczt,cols(i),savedir,[cat{i} '_NC_Z'],saveflag);
		end
    
    case 'other'
        savedir=[emmdir '/emm_0012/' outdir '/pics/'];
        for i=1:length(files)
            if i==2
                clims=[0 5];
            else
                clims=[];
            end
            
            comm=['plotemmtimHf(temm,zemm,' files{i} ',1,savedir,[files{i} ''_Z''],saveflag,clims);' ];
            eval(comm);
        end
    case 'w'
        tgrid=W(:,3);
		tgrid=reshape(tgrid,[s1 s2]); %time grid
        ltemm=length(temm);
        plotemmtimHf(tgrid(:,1:ltemm),hemm(:,1:ltemm),wemm(:,1:ltemm),1,savedir,'_W_Z',saveflag)
    
end
