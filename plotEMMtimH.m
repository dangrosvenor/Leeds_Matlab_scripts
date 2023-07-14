cat{1}='Cloudwater';
cat{2}='Condensation Freezing';
cat{3}='Contact Nucleation';
cat{4}='Homogeneous Nucleation';
cat{5}='Splinters 1';
cat{6}='Frozen Rain 1';
cat{7}='Splinters 2';
cat{8}='RF Splinters 1';
cat{9}='Frozen Rain 2';
cat{10}='Splinters 3';
cat{11}='RF Splinters 2';
cat{12}='Rainwater';

cols=[1 4:14];

typ='nc';
switch typ
	case 'mr'
	savedir='c:/documents and settings/g/my documents/emm_0012/output/pics/iwc/';
	
	for i=1:length(cat)
        plotemmtimHf(temm,Temm,iwczt,cols(i),savedir,cat{i});
	end
	
	case 'nc'
	savedir='c:/documents and settings/g/my documents/emm_0012/output/pics/inc/';
	
	for i=1:length(cat)
        plotemmtimHf(temm,Temm,inczt,cols(i),savedir,[cat{i} 'NC']);
	end

end
