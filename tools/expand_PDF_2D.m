function [xdat,ydat] = expand_PDF_2D(x,y,F)
%function dat = expand_PDF(x,F)
%expands a 2D PDF into a two linear arrays of elements representing all of the
%individual data points (and with the x and y matching up)
%i.e. takes all of the x (mid-point) values and makes F_i occurences of x_i
%for all i ( =1:length(x) )


Ndat=sum(sum(F));
if ~isnan(Ndat)
    xdat=NaN*ones([Ndat 1]);
    ydat=NaN*ones([Ndat 1]);
else
    xdat=NaN;
    ydat=NaN;
    return
end


n=1; %counter for the linear array to be created
for i=1:length(x)
    for j=1:length(y)

        N=F(i,j);
        xdat(n:n+N-1)=x(i);
        ydat(n:n+N-1)=y(j);        
        n=n+N;

    end
end
    