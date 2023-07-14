comp='pc';
switch comp
case 'laptop'
	savedir='c:/documents and settings/g/my documents/emm_0012/output/pics/';
case 'pc'
    savedir='c:/documents and settings/login/my documents/emm_0012/output/pics/';
end
plotemmtimHf(temm,Temm,lwc,1,savedir,'lwc');