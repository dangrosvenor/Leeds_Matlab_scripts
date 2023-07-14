%June 2007 - this is the one to use as of this date
idir=1;

typ='mr';
%typ='nc';
typ='other';
%typ='w';



idir2=runs(idir);

saveflag=1;

if ~exist('savedir'); savedir='dir';end
clims=[];


if ~exist('idir'); idir=1; end


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


files={'naero' 'ssat' 'lwc' 'ncw' 'rwc' 'nr','qup'};

files={files{[1:7]}};
%files={files{[2]}};

comp='pc';
comp='laptop';
%comp='lacie';
comp='readEMM2';
comp='ydrive';

switch comp
case 'laptop'
	emmdir=['c:/documents and settings/g/my documents/emm_0012/'];
case 'pc'
    emmdir=['c:/documents and settings/login/my documents/emm_0012/'];
case 'lacie'
    emmdir=['e:/emm_0012/'];
case 'readEMM2'
    emmdir=rootdir;
case 'ydrive'
    emmdir='y:/';    
    
end



tlims=[0 20.6];
zlims=[0 14];

tlims=[];
zlims=[];

tag='2030_UTC-2';
tag='1';

zflag=1;

if zflag==1
    zdat=vec(idir).z;
else
    zdat=Temm;
end

temm=vec(idir).time;
Temm=vec(idir).temp;



switch typ
	case 'mr'
%		savedir=[emmdir 'emm_0012/' outdir{idir2} '/pics/iwc/'];
		savedir=['y:/' outdir{idir2} '/pics/iwc/'];
        
		for i=1:length(cat)
            plotemmtimHf(temm,zdat,emmdat(idir).iwczt,cols(i),savedir,[tag '-' cat{i}],saveflag,[],tlims,zlims,zflag);
		end
	
	case 'nc'
		savedir=[emmdir outdir{idir2} '/pics/inc/'];
		
		for i=1:length(cat)

            plotemmtimHf(temm,zdat,emmdat(idir).inczt,cols(i),savedir,[tag '-' cat{i} 'NC'],saveflag,clims,tlims,zlims,zflag);
		end
    
    case 'other'
%        savedir=[emmdir '/emm_0012/' outdir{idir2} '/pics/'];
        savedir=[emmdir outdir{idir2} '/pics/'];
        
        for i=1:length(files)
            if strmatch(files{i},'ssat')==1
                clims=[0 6];
                clims=[];
            elseif strmatch(files{i},'lwc')==1
                clims=[0 3.2];
            elseif strmatch(files{i},'rwc')==1
                clims=[0 3.2];    
            else
                clims=[];
            end
            
            comm=['plotemmtimHf(temm,zdat,emmdat(idir).' files{i} ',1,savedir,[tag ''-'' files{i}],saveflag,clims,tlims,zlims,zflag);' ];
            eval(comm);
        end
    case 'w'
      tlims=[vec(idir).time(1) vec(idir).time(end)];
      zlims=[vec(idir).z(1) vec(idir).z(end)];
        
        savedir=[emmdir '/emm_0012/' outdir{idir2} '/pics/'];
        plotemmtimHf(temm2(1,:),Zemm(:,1),wemm(:,:),1,savedir,'_W',saveflag,clims,tlims,zlims,zflag)
    
end
