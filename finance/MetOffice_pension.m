ages=[66:70]; 
inflation_percent=0;
current_age=42;
clear leg_str

salary=45; 
clear t2; 
iage=0; 
for age=ages
    iage=iage+1;
    t=0; 
    for i=1:age-current_age;
        t=t+salary*2.32/100;
        t=t*(1+inflation_percent/100);
    end;
    t2(iage)=t;
end

figure('color','w');
plot(ages,t2,'bo-');
leg_str{1} = ['Salary = £' num2str(salary) 'k'];
xlabel('Retirment age');
ylabel('Annual income on retirement (k£)');
%title(['Based on a salary of £' num2str(salary) 'k']);

grid on
hold on


salary=35; 
clear t2; 
iage=0; 
for age=ages
    iage=iage+1;
    t=0; 
    for i=1:age-current_age;
        t=t+salary*2.32/100;
        t=t*(1+inflation_percent/100);
    end;
    t2(iage)=t;
end

plot(ages,t2,'ro-');
leg_str{2} = ['Salary = £' num2str(salary) 'k'];

legend(leg_str);




