function [savings] = saver_func(nyears,monthly,rate)

f=12*100;

savings=0;
for i=1:nyears*f
    savings = savings + savings*rate/100/f;
    savings = savings + monthly*12/f;
end
