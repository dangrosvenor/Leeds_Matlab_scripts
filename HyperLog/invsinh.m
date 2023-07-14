function y=invsinh(x)

for i=1:length(x)
	y(i)=fzero(@invsinhROOT,0,[],x(i));
end

