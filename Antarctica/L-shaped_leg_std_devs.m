clear std_leg mean_leg
for i=1:length(xdat)
    iend=find(xdat(1).x<-66.8));
    std_dev(i) = std(ydat(i).y(1:iend));
    mean_dev(i) = mean(ydat(i).y(1:iend));
end

mean_leg = mean_dev;
std_leg = std_dev;    