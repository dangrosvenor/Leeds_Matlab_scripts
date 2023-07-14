%imports n files directly from .DG00xx files from n directories in direc(n).text
clear SER TwoD TwoDDan
fnall=input('Enter the Diagnostic file number, with directories in direcDan(i).dir: ');
ndir=input('Enter the no. of directories to import :');

size(fnall);
fnsiz=ans(2);

for i=1:ndir;
    
    if fnsiz>1
        fn=fnall(i);
    else
        fn=fnall;
    end
    
if fn<10;
    fn=strcat('0',int2str(fn));
else;
    fn=int2str(fn);
end;


    FileName=strcat('c:\cygwin\home\user\',direcDan(i).dir,'\RUN0001.DG00',fn);
    DIAG2_3_v2;
    GridDan(i)=Grid;
    TwoDDan(i)=TwoD;
    TwoDDan(i).Q(:,1,1)=0;
    SerDan(i).SER=SER;
    textdataDan(i).text=textdat;
    TimeAvDan(i)=TimeAv;
end;
[DG,SERSTR]=SORTHEADER(HEADER2);