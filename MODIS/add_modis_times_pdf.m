range_X = max(X)-min(X);
range_Y = max(Y)-min(Y);

for i=1:length(ihtot)
    
    text(X_orig(i)+range_X/100,Y_orig(i)+range_Y/100,datestr(scantime_matlab_tim_pdf(ihtot(i)),31));
    fprintf(1,'\n%s',datestr(scantime_matlab_tim_pdf(ihtot(i)),31));
    
end