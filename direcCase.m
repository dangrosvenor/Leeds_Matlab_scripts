flag=8;


if flag==8
    direcDan(1).dir='g:casestudy/c1000-1200s';   %114
    direcDan(2).dir='g:casestudy/c1000NormWind'; %121
    direcDan(3).dir='g:casestudy/c125-1200s';    %121
    direcDan(4).dir='g:casestudy/c250-1200s';    %121
    direcDan(5).dir='g:casestudy/c500-1200s';    %121
    direcDan(6).dir='g:casestudy/c2000-500m';    %120
    direcDan(7).dir='g:casestudy/c1000gal2';     %123
    direcDan(8).dir='g:casestudy/cold1000kmspin';%147
end

if flag==7
    direcDan(1).dir='g:casestudy/hot243+20';   %97
    direcDan(2).dir='g:casestudy/hot243-20';   %97
    direcDan(3).dir='g:casestudy/hot243-ccn960';   %97
    direcDan(4).dir='g:casestudy/hot243sslowjetspin';   %97
    direcDan(5).dir='g:casestudy/cold1000kmspin';   %91
    direcDan(6).dir='g:casestudy/hot234s-lowjet';   %49
    direcDan(7).dir='g:casestudy/hot243-ref';   %23
    direcDan(8).dir='g:casestudy/hot243-20mins';   %14
end

if flag==6
    direcDan(1).dir='g:casestudy/hotpool5000-200pts'; %30
    direcDan(2).dir='g:casestudy/coldpool';   %13
    direcDan(3).dir='g:casestudy/coldpool2';   %13
    direcDan(4).dir='g:casestudy/coldpool2test';   %16
    direcDan(5).dir='g:casestudy/coldpoolspin5000';   %30
    direcDan(6).dir='g:casestudy/coldpoolspin2';   %5
    direcDan(7).dir='g:casestudy/hotpoolspin2';   %5
end


if flag==4
	direcDan(1).dir='baurufield/13.02.18alt';   %71
    direcDan(2).dir='baurufield/13.02.18';      %48
    direcDan(3).dir='g:casestudy/03.03.12alt-50'; %51
    direcDan(4).dir='g:casestudy/24.02.2015altSpin'; %45
    direcDan(5).dir='g:casestudy/24.02.2015alt18km'; %44
    direcDan(6).dir='g:casestudy/24.02.2015alt50vert'; %30
    direcDan(7).dir='g:casestudy/24.02.2015alt5deg'; %30
	direcDan(8).dir='g:casestudy/24.02.2015alt';  %30
    direcDan(9).dir='g:casestudy/24.02.2015alt2deg';   %29
	direcDan(10).dir='g:casestudy/24.02.2015alt4deg'; %23
    direcDan(11).dir='g:casestudy/24.02.2015'; %21
    direcDan(12).dir='baurufield/03.03.12';      %73
    direcDan(13).dir='g:casestudy/24.02alt4degccn480'; %10
    
    
	
	%direcDan(5).dir='g:casestudy/24.02.2015alt';
end

if flag==5
    direcDan(1).dir='g:casestudy/24.02.2015alt5deg'; %30
    direcDan(2).dir='g:casestudy/24.02.2015alt';  %30
    direcDan(3).dir='g:casestudy/24.02.2015alt2deg';   %29
	direcDan(4).dir='g:casestudy/24.02.2015alt4deg'; %23
    direcDan(5).dir='g:casestudy/24.02.2015'; %21
    
    	%direcDan(5).dir='g:casestudy/24.02.2015alt';
end


if flag==1
	direcDan(1).dir='g:fluxes/flux2';      %63
	direcDan(2).dir='g:fluxes/fluxlowres';   %60
	direcDan(3).dir='g:fluxes/fluxlowres2000'; %172
	%direcDan(4).dir='baurufield/24.02.2015alt'; %73
	%direcDan(5).dir='baurufield/03.03.12';      %73
end


if flag==2
    direcDan(1).dir='baurufield/03.03.12';      %73
    direcDan(2).dir='g:casestudy/03.03.12alt-50'; %51
	direcDan(3).dir='baurufield/13.02.18';      %48
	direcDan(4).dir='baurufield/13.02.18alt';   %71
	direcDan(5).dir='baurufield/24.02.2015';    %45 - problem with 45 - think will have to start again if want to continue
	direcDan(6).dir='baurufield/24.02.2015alt'; %73
	
end

if flag==3
	direcDan(1).dir='g:fluxes/onesineinv';      %48
	direcDan(2).dir='baurufield/13.02.18alt';   %71
	direcDan(3).dir='baurufield/24.02.2015';    %45 - problem with 45 - think will have to start again if want to continue
	direcDan(4).dir='baurufield/24.02.2015alt'; %73
	direcDan(5).dir='baurufield/03.03.12';      %73
end