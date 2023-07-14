function [diffdan]=diffdan(y,x,n)
%central differencing routine
%[diffdan]=diffdan(y,x,n)
%finds dy/dx along dimension n
%first answer is [y(3)-y(1)]/[(x(3) -x(1)] and result will be of original length along n minus 2

if ~exist('n'); n=1; end

[diffy,err]=shiftmat(y,n);
if err==1
    break
end

[diffx,err]=shiftmat(x,n);
if err==1
    break
end


diffdan=diffy./diffx;