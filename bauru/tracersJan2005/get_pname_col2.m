function imphys=get_pname_col2(dgcell)

pnames; %puts process rate strings into pname(i).p
            
for idg=1:length(dgcell)            
	for imp=1:length(pname)
		if strcmp(upper(pname(imp).p),upper(dgcell{idg}))==1
            imphys(idg)=imp;
            break
		end
	end
end