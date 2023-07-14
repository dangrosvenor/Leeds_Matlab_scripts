clear wd wb



colors={'b-','r-.','k--','g'};
%figure;

for i=1:4
    ntim=size(wbb(i).w,2);
    
	step=wbb(i).w(2,end)-wbb(i).w(1,end);
	minw=minALL(wbb(i).w);
	maxw=maxALL(wbb(i).w);
	
	wb(i).w=minw:step:maxw;
	
	n=0;
	for j=2:length(wb(i).w)
        [a b]=find(wbb(i).w==wb(i).w(j));
        a=a-1; %values stored in one less than the wbb index
        b(a==0)=[];
        a(a==0)=[];
        [x]=sub2ind(size(wdist(i).w),a,b);
        n=n+length(x);
        wd(i).w(j)=sum(wdist(i).w(x))/250/ntim; %sum up all the values in this bin and normalise to number of points
	end
	
    [a]=find(wd(i).w>0);
    wd(i).w(a(end)+1:end)=[];
    wb(i).w(a(end)+1:end)=[];
    
    
	
	%plot(wb(i).w,log10(wd(i).w),colors{i});
	%hold on;
    
end

%	xlabel('Vertical Velocity (m/s)');
%	ylabel('log10 of normalised dn/dw where n is no. grid points');