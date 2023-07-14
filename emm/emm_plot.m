
isave=0;
plotcase=4;
x=100;



for plotcase=1:4

figure;    
switch plotcase

case 1
    x=103;
	plot(D1(4,end-x:end),D1(3,end-x:end),'x-');
	xlabel('droplet_N (/kg)');
    ylabel('z (km)');
    picname=['c:/cygwin/home/user/emm/debug_graphs/droplet_N'];
    
case 2
	plot(N1(3,end-x:end),N1(8,end-x:end),'x-');
	xlabel('N-aero (/kg)');
    ylabel('z (km)');
    picname=['c:/cygwin/home/user/emm/debug_graphs/N_aero'];
    
case 3
	plot(N1(4,end-x:end),N1(8,end-x:end),'x-');
	xlabel('N-w-source (/kg)');
    ylabel('z (km)');
    picname=['c:/cygwin/home/user/emm/debug_graphs/N_w_source2'];    
    
case 4
	plot(N1(6,end-x:end),N1(8,end-x:end),'x-');
	xlabel('s (%)');
    ylabel('z (km)');
    picname=['c:/cygwin/home/user/emm/debug_graphs/s_l'];    
    
end  
    
    
	
	if isave==1
        set(gcf,'paperpositionmode','auto');
		print(gcf,[picname '.emf'],'-dmeta');
	end
    
   
 
end
