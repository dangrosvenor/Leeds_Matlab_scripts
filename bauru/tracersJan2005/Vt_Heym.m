clear vt

n=2.25; %column
A=0.26;
k=0.49;

n=2.1; %6 bullet rosette
A=0.05;
k=0.97*0.81;

n=2.52; %8 bullet rosette
A=0.32;
k=0.52;

n=2.38; %planar crystals P1a-P1b-P1d
A=0.79;
k=0.084;

g=9.81*100;
rhoA=1.2/1000;
T=213;

vis=airprop2(T,'my')*10; %*10 to convert to cgs

D=[0:1000e-9:1000e-6]*100; %think D should be in cm

X=4/3 * g*k*A^(n-1)*D.^3 / rhoA /vis^2;

lims(1).l=[0.01 10];
lims(2).l=[10 585];
lims(3).l=[585 1.56e5];
lims(4).l=[1.56e5 1e99];

afs=[0.04394 0.06049 0.2072 1.0865];
bfs=[0.970 0.831 0.638 0.499];

for j=1:length(lims)
	
    i=find(X>lims(j).l(1) & X<=lims(j).l(2));
	af=afs(j);
	bf=bfs(j);
	vt(i)=Vt(af,bf,vis,D(i),A,k,rhoA,n,g)/100; %in m/s

end