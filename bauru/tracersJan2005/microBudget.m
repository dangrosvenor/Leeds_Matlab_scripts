     dgs={'PGDEP','PGMLT', ...  %1-2
          'PRAUT',   'PGSHD',   'PRACW',   'PSMLT', ... %3-6
          'PIMLT',   'PSAUT',   'PSDEP',   'PIACR_G', ... %7-10
          'PSACI',   'PGACI',   'PGACW',   'PGACS', ... %11-14
          'PGACR',   'PSACR',   'PRACS',   'PSACW', ... %15-18
          'PGAUT',   'PGFR',    'PRACI_G', 'PGWET', ... %19-22
          'PGDRY',   'PGSUB',   'PSSUB',   'PREVP', ... %23-26
          'PISUB',   'DQI  ',   'PIHAL',   'PIPRM', ... %27-30
          'PIDEP',   'PIACW',   'PICNT',   'PIFRW', ... %31-34
          'PIACR_S', 'PRACI_S', 'PRACI',   'PIACR', ... %35-38
          'RIACR_S', 'RSACR',   'RIACR_G', 'RGFR', ...  %39-42
          'RGACS',   'RGACR',   'RSAUT',   'RIACI', ... %43-46
          'RSACS',   'RSBRK'  ...                       %47-48
      , 'RGmlt'		  ... %49
      , 'RSmlt'		   ...%50
      , 'RImlt'		  ...
      , 'RSaci'		  ...
      , 'RGaci'		  ... 
      , 'RGacw'		  ...  %??? - not used?
      , 'RRacs'		  ... %55
      , 'RGaut'		  ...
      , 'RGsub'		  ... 
      , 'RSsub'		  ...
      , 'RIcnt'		  ...
      , 'RIhal'		  ... %60
      , 'RIprm'		  ... 
      , 'RIfrw'		  ... 
      , 'RRaci_S'	  ... 
      , 'RRaci_G'	  ...
      , 'RSaut_I'	  ... %65
      , 'RSaut_S'	  ...
      , 'RGaut_S'	  ... 
      , 'RGaut_G'	  ...
      , 'RSacr_S'	  ...
      , 'RSacr_G' }   ;   %70


all=1; %flag to say if want to sum in icediag arrays over all times or whether are using values from a single dump )=0)


%ice number sources ***********************************************************

insour=[59 60 61 62];
insink=[53 52 63 64 46 65];

imi0=[29 30 33 34]; %all sources
imass=[12 11 36 21]; %all sinks
inorm=[45 46]; %all sinks
        
%    dqNI(J)=RIhal(J)+RIprm(J)+RIcnt(J)+RIfrw(J)
%     &         -(RGaci(J)+RSaci(J)+
%     &         RRaci_S(J)+RRaci_G(J))
%     &         -RIaci(J)-RSaut(J)

%           dqNI(J)=MAX( (1.-AQNIjk)*RDT,
%     &         MIN( (rmax_conc-AQNIjk)*RDT,
%     &         (PIhal(J)+PIprm(J)+PIcnt(J)+RHOMOG*PIfrw(J))/RMI0
%     &         -(PGaci(J)+PSaci(J)+
%     &         PRaci_S(J)+PRaci_G(J))*Rmass_I
%     &         -RIaci(J)-RSaut(J) ))



jc=1;
lz=length(GridDan(jc).Z);

source=zeros([lz 1]);
sink=zeros([lz 1]);

mi0=zeros([lz 1]);
mass=zeros([lz 1]);
norm=zeros([lz 1]);

% siz=size(icediag_num(jc).i);
% 
% source=zeros([siz(1) siz(2)]);
% sink(1:length(GridDan(jc).Z))=0;
% 
% mi0(1:length(GridDan(jc).Z))=0;
% mass(1:length(GridDan(jc).Z))=0;
% norm(1:length(GridDan(jc).Z))=0;


for i=1:length(insour)
    if all==0
        source=source+getDGAVsNoArea(['ALL_' dgs{insour(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        source=source+sum(icediag_nums(jc).i(:,:,insour(i)),2);
    end
end
for i=1:length(insink)
    if all==0
        sink=sink+getDGAVsNoArea(['ALL_' dgs{insink(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        sink=sink+sum(icediag_nums(jc).i(:,:,insink(i)),2);
    end
end
for i=1:length(imi0)  
    if all==0
        mi0=mi0+getDGAVsNoArea(['ALL_' dgs{imi0(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        mi0=mi0+sum(icediag_nums(jc).i(:,:,imi0(i)),2);
    end
end
for i=1:length(imass)  
    if all==0
        mass=mass+getDGAVsNoArea(['ALL_' dgs{imass(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        mass=mass+sum(icediag_nums(jc).i(:,:,imass(i)),2);
    end
end
for i=1:length(inorm)    
    if all==0
        norm=norm+getDGAVsNoArea(['ALL_' dgs{inorm(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        norm=norm+sum(icediag_nums(jc).i(:,:,inorm(i)),2);
    end
end


total=source-sink;

qi=getDGAVsNoArea('ALL_Q06',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1)';
ni=getDGAVsNoArea('ALL_Q07',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1)';

RmassIa=ni./qi;
RmassI=RmassIa;
RmassI(RmassIa==Inf)=0;
RmassI(ni==0)=0;
RmassI=RmassI';

total2 = mi0/1e-15 - mass .* RmassI - norm;



%snow number sources *****************************************************************

!     &        dqNS(J)=MAX( (1.-AQNSjk)*RDT ,
!     &               RIacr_S(J)+RSaut(J)+RSbrk(J)
!     &             -(PSsub(J)+PSmlt(J)+PGaut(J))*Rmass_S
!     &             -RGacs(J)-RSacr(J)-RSacs(J) )

% dqNS(J) = RIacr_S(J)+RSaut(J)+RSbrk(J)
%      &             -(RSsub(J)+RSmlt(J)+RGaut(J))
%      &             -RGacs(J)-RSacr(J)-RSacs(J)
    

issour=[39 66 48];
issink=[43 69 47 58 50 67];

isnormsour=[39 45 48];
isnormsink=[43 40 47];
ismass=[25 6 49]; %a sink


ssource=zeros([lz 1]);
ssink=zeros([lz 1]);

smass=zeros([lz 1]);
snormsour=zeros([lz 1]);
snormsink=zeros([lz 1]);

for i=1:length(issour)
    if all==0
        ssource=ssource+getDGAVsNoArea(['ALL_' dgs{issour(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
	else
        ssource=ssource+sum(icediag_nums(jc).i(:,:,issour(i)),2);
    end
end
for i=1:length(issink)
    if all==0
        ssink=ssink+getDGAVsNoArea(['ALL_' dgs{issink(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        ssink=ssink+sum(icediag_nums(jc).i(:,:,issink(i)),2);
    end
end
for i=1:length(ismass)
    if all==0
		smass=smass+getDGAVsNoArea(['ALL_' dgs{ismass(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        smass=smass+sum(icediag_nums(jc).i(:,:,ismass(i)),2);
    end
end
for i=1:length(isnormsour)  
    if all==0
        snormsour=snormsour+getDGAVsNoArea(['ALL_' dgs{isnormsour(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        snormsour=snormsour+sum(icediag_nums(jc).i(:,:,isnormsour(i)),2);
   end
end
for i=1:length(isnormsink)  
    if all==0
        snormsink=snormsink+getDGAVsNoArea(['ALL_' dgs{isnormsink(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        snormsink=snormsink+sum(icediag_nums(jc).i(:,:,isnormsink(i)),2);
    end
end


stotal=ssource-ssink;

qs=getDGAVsNoArea('ALL_Q04',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1)';
ns=getDGAVsNoArea('ALL_Q09',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1)';

RmassSa=ns./qs;
RmassS=RmassSa;
RmassS(RmassSa==Inf)=0;
RmassS(ns==0)=0;
RmassS=RmassS';

stotal2 = snormsour - snormsink - smass .* RmassS;



%graupel number sources *****************************************************************8


        
%               dqNG(J)=RGaut(J)
%      &               +RSacr(J)+RIacr_G(J)+RGfr(J)
%      &             -(RGsub(J)+RGmlt(J))
% 
%     
%     &        dqNG(J)=MAX( (1.-AQNGjk)*RDT ,
%     &                PGaut(J)*Rmass_S 
%     &               +RSacr(J)+RIacr_G(J)+RGfr(J)
%     &             -(PGsub(J)+PGmlt(J))*Rmass_G )

igsour=[68 70 41 42];
igsink=[57 49];

igmass_g=[22 2]; %all sinks
igmass_s=[19];
ignorm=[38 41 42]; %all sinks


jc=1;

gsource=zeros([lz 1]);
gsink=zeros([lz 1]);

gmass_s=zeros([lz 1]);
gmass_g=zeros([lz 1]);
gnorm=zeros([lz 1]);


for i=1:length(igsour)
    if all==0
        gsource=gsource+getDGAVsNoArea(['ALL_' dgs{igsour(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        gsource=gsource+sum(icediag_nums(jc).i(:,:,igsour(i)),2);
    end
end
for i=1:length(igsink)
    if all==0
        gsink=gsink+getDGAVsNoArea(['ALL_' dgs{igsink(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        gsink=gsink+sum(icediag_nums(jc).i(:,:,igsink(i)),2);
    end
end
for i=1:length(igmass_s)    
    if all==0
        gmass_s=gmass_s+getDGAVsNoArea(['ALL_' dgs{igmass_s(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        gmass_s=gmass_s+sum(icediag_nums(jc).i(:,:,igmass_s(i)),2);
    end
end
for i=1:length(igmass_g)    
    if all==0
        gmass_g=gmass_g+getDGAVsNoArea(['ALL_' dgs{igmass_g(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        gmass_g=gmass_g+sum(icediag_nums(jc).i(:,:,igmass_g(i)),2);
    end
end
for i=1:length(ignorm)    
    if all==0
        gnorm=gnorm+getDGAVsNoArea(['ALL_' dgs{ignorm(i)}],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    else
        gnorm=gnorm+sum(icediag_nums(jc).i(:,:,ignorm(i)),2);
    end
end


gtotal=gsource-gsink;

qg=getDGAVsNoArea('ALL_Q05',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1)';
ng=getDGAVsNoArea('ALL_Q08',dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1)';

RmassGa=ni./qi;
RmassG=RmassGa;
RmassG(RmassGa==Inf)=0;
RmassG(ng==0)=0;
RmassG=RmassG';

gtotal2 = gmass_s .* RmassS - gmass_g .* RmassG + gnorm;


if all==1
	dqni=sum(icediagALL(1).i(:,:,34),2);
	dqns=sum(icediagALL(1).i(:,:,36),2);
	dqng=sum(icediagALL(1).i(:,:,35),2);
else
    dqni=getDGAVsNoArea(['ALL_DQ07'],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    dqns=getDGAVsNoArea(['ALL_DQ09'],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
    dqng=getDGAVsNoArea(['ALL_DQ08'],dgstrDan(jc).dg,TimeAvDan(jc).DGAV,1,1);
end


%snow mr sources

    ismrsour=[8 9 18 11 36 35];
    ismrsink=[25 14 17 19 6];




