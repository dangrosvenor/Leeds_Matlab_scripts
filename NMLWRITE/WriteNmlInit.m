function [RH,U,V,TH,IndicesRH,IndicesU,IndicesV,IndicesTH]=WriteNmlInit(FileName,data);

% function to write the nml file from trmm data
[U,V,TH,RH,IndicesRH,IndicesU,IndicesV,IndicesTH]=NmlGrid(data);

maxind=max(max(IndicesRH,IndicesU));
maxind=max(max(maxind,IndicesV));
Z=data(:,1)-data(1,1);
% open the file
fid=fopen(FileName,'w');

fprintf(fid,'&CNTRL NN=4,NNDIAG=4,NNDUMP=3,ISTART=1 &END\n');
fprintf(fid,'&TIMENML NTMHALT=1,TIMHALT(1)=1000.,IPRTDG=1,NSTEPMAX=10000 &END\n');
fprintf(fid,' &JOBINFO FILEA=''filea'',FILEB=''fileb'',FILEZ=''filez'',\n');
fprintf(fid,' MSGFILE=''mesg'',USERID=''MAP4.HS'' &END\n');
NRUN=input('Enter the run number');
fprintf(fid,'&INPUT NRUN=%d,Z0=0.0001,Z0TH=0.00001,\n',NRUN);
fprintf(fid,' PSF=%.1f,SHFLX_SEN=600.0,SHFLX_LAT(1)=600.0 &END\n',100.*data(1,2));
fprintf(fid,'&GRID NSMTH=5,HGD(1)=%.1f,KGD(1)=100,,\n',Z(maxind));
fprintf(fid,' ZZTOP=%.1f &END\n',Z(maxind));
fprintf(fid,'&THPROF\n');

size(IndicesTH);
imax=ans(2);
for i=1:imax
    fprintf(fid,' ZNREF_READ(%d)=%.1f, THREF_READ(%d)=%.1f,\n',i,Z(IndicesTH(i)),i,TH(IndicesTH(i)));
end
fprintf(fid,'THREF0=%.1f\n',TH(IndicesTH(1)));
for i=1:imax
    fprintf(fid,' ZNINIT_READ(%d)=%.1f, THINIT_READ(%d)=%.1f,\n',i,Z(IndicesTH(i)),i,TH(IndicesTH(i)));
end
for i=1:imax
    fprintf(fid,' UINIT(%d)=%.1f, ZUINIT(%d)=%.1f,\n',i,U(IndicesU(i)),i,Z(IndicesU(i)));
end
for i=1:imax
    fprintf(fid,' VINIT(%d)=%.1f, ZVINIT(%d)=%.1f,\n',i,V(IndicesV(i)),i,Z(IndicesV(i)));
end
for i=1:imax
    fprintf(fid,' ZHINIT(%d)=%.1f, HINIT(%d)=%.3f,\n',i,Z(IndicesRH(i)),i,RH(IndicesRH(i)));
end
fprintf(fid,'T_INIT=1.0,ZT_INIT=1300.0,W_INIT=1500.0,NT_INIT=1 &END\n');
fprintf(fid,' &SUBMODEL SUBB=40.0,SUBC=16.0,SUBG=1.2,SUBH=0.0,\n');
fprintf(fid,' SUBP=1.,SUBQ=1.,SUBR=4.,ATH2_N=0.3,A2_N=0.23,PR_N=0.7,\n');
fprintf(fid,'  RIC=0.25,RMLMAX=23. &END\n');
fprintf(fid,' &DIAGNOST\n');
fprintf(fid,' IDGU=3,IDGV=3,IDGW=3,IDGTH=3,IDGQ(1)=3,\n');
fprintf(fid,' IDGCL=2,IDGAV=1,IDGTL=2,IDGP=3,IDGPV=0,IDGPD(1)=1,\n');
fprintf(fid,' &END\n');
fprintf(fid,'&DYNAMICS UG0=0.,VG0=0.,DUGDZ=0.,DVGDZ=0.,FCORIOL=0.0 &END\n');
fprintf(fid,'&NUMERICS DXX=100.,DYY=100.,DTM=1.0,DTMMAX=1.,DTMMIN=0.01 &END\n');
fprintf(fid,'&PHYSICS &END\n');
fprintf(fid,'&DAMPNML DMPTIM=0.001,ZDMP=%.1f,HDMP=2000. &END\n',Z(maxind)-2000.);
fprintf(fid,'&OVRIDE1 &END');


fclose(fid);