clear Surf Surf2
direcDan(1).dir='forcecons';
direcDan(2).dir='forceconsccn720';


load 'C:\matlabR12\work\bauru\casestudy\forcecons\diag\Surf_2-169';
load 'C:\matlabR12\work\bauru\casestudy\forcecons\diag\echotops_2-169';
load 'C:\matlabR12\work\bauru\casestudy\forcecons\diag\cappi_2-169';
Surf.echotops=echotops./1000;
Surf.cappi=cappi;
Surf.cappi2=[];
Surf2(2)=Surf;

load 'C:\matlabR12\work\bauru\casestudy\forceconsccn720\diag\Surf-2-169';
load 'C:\matlabR12\work\bauru\casestudy\forceconsccn720\diag\echotops-2-169';
load 'C:\matlabR12\work\bauru\casestudy\forceconsccn720\diag\cappi-2-169';
Surf(1).echotops=echotops./1000;
Surf(1).cappi=cappi;
Surf.cappi2=[];
Surf2(1)=Surf;  %Surf2(1)=ccn720 2=normal


%Surf(2)=Surf2;

direcDan(1).dir='c1000-1200s';   %114
direcDan(2).dir='c125km';    %121
direcDan(3).dir='c250km';    %121
direcDan(4).dir='c500km';    %121
direcDan(5).dir='c2000-500m';    %120

for i=1:length(direcDan)
    dire=strcat('c:/matlabr12/work/bauru/casestudy/',direcDan(i).dir,'/diag/');
    d=dir(dire);
    for j=3:length(d);
        if strcmp(d(j).name(1:4),'Surf')==1
            surfname=d(j).name;
        end
    end
    load(strcat(dire,surfname));
    try
        Surf2(2+i).MaxW=Surf.MaxW;
        Surf2(2+i).RainMR=Surf.RainMR;
        Surf2(2+i).prec=Surf.prec;
        Surf2(2+i).pmax=Surf.pmax;
        Surf2(2+i).instant=Surf.instant;
       % Surf2(2+i).cappi2=Surf.cappi2;
    catch
    end
        
        Surf2(2+i).echotops=Surf.echotops./1000;
        Surf2(2+i).cappi=Surf.cappi;
    end


for i=1:3     %length(Surf2)
    Surf2(i).cappi2=Surf2(i).cappi;
    Surf2(i).cappi=[];
    temp=Surf2(i).echotops;
    Surf2(i).echotops=[];
    
    for j=1:size(temp,2)
        Surf2(i).echotops(:,1,j)=temp(:,j);
    end
    
    for j=1:size(Surf2(i).cappi2,3)
        Surf2(i).cappi(1:size(Surf2(i).cappi2,2),1,j)=Surf2(i).cappi2(2,:,j);
    end
end

