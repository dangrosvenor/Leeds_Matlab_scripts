clear Surf Surf2

load 'C:\matlabR12\work\bauru\casestudy\forceconsccn720\diag\Surf-2-169';
load 'C:\matlabR12\work\bauru\casestudy\forceconsccn720\diag\echotops-2-169';
load 'C:\matlabR12\work\bauru\casestudy\forceconsccn720\diag\cappi-2-169';
Surf.echotops=echotops./1000;
Surf.cappi=cappi;
Surf2(1)=Surf;

clear Surf;
load 'C:\matlabR12\work\bauru\casestudy\forcecons\diag\Surf_2-169';
load 'C:\matlabR12\work\bauru\casestudy\forcecons\diag\echotops_2-169';
load 'C:\matlabR12\work\bauru\casestudy\forcecons\diag\cappi_2-169';
Surf.echotops=echotops./1000;
Surf.cappi=cappi;
Surf2(2)=Surf;





for i=1:length(Surf2)
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
        