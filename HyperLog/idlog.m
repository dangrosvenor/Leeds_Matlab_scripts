function d=idlog(x,xc)
%inverse of dlog

i=find(x>=log10(xc));
d(i)=10.^(x(i));

i=find(x<log10(xc)-1);
d(i)=-10.^(-(x(i)+1-2*log10(xc)));

i=find(x<log10(xc) & x>=log10(xc)-1);
a=log10(xc);
b=a-1;
d(i)=-2*xc*(x(i)-a)/(b-a) + xc;

d=reshape(d,size(x));

% for i=1:length(x)
%     
% 	if x(i)>=log10(xc)
%         d(i)=10^(x(i));
%     elseif x(i)<log10(xc)-1
% 		d(i)=-10^(-(x(i)+1-2*log10(xc)));
% 	else
%         a=log10(xc);
%         b=a-1;
%         d(i)=-2*xc*(x(i)-a)/(b-a) + xc;
% 	end
% 
% end