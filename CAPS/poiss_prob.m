%calculate the Poisson probability for >=n events
%i.e. prob of getting n or more events when sampling L litres of air with mean conc. of C

n=[1:25];

C=0.0; %concentration per L
N=0;
%L=100.4; %volume of air sampled
%L=50;
L=N/C;
L=10.4*10;

a=L*C; %expected number of ice crystals in L litres


for j=1:length(n)

    sum_p=0;
    for i=0:n(j)-1
        sum_p = sum_p + poisspdf(i,a);
    end

    fprintf(1,'\nFor n>=%d (%f per L) prob = %f %%',n(j),n(j)/L,100*(1-sum_p));

end

fprintf(1,'\n');

Cs=[0.1 0.2 0.4 0.6 0.8 1.6 3.2 10]*10.4/L;



p=90; %percentage confidence interval required

dp=(1-p/100)/2;
for i=1:length(Cs)

    %for the range of concentrations over which 90 % of the particles fall we can use
    C=Cs(i);
    Ncrystals=poissinv([dp 1-dp],C*L);
    conc_range = Ncrystals/L;
    perror= 100*( poissinv([1-dp],C*L)/L - C )/C;   %=56.25% for this case
    %the error would then be +/- from the mean (C)
    
    fprintf('\nNmean= %f N=%f to %f conc=%f to %f error=%f %%error=%f',C*L,Ncrystals,conc_range,perror*C/100,perror);

end
fprintf(1,'\n');


