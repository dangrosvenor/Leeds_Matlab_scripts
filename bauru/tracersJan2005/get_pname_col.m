function imphys=get_pname_col(dgcell)

pnames; %puts process rate strings into pname(i).p
            
for idg=1:length(dgcell)            
	for imp=1:length(pname)
		if strcmp(pname(imp).p,dgcell{idg})==1
            imphys(1).p(idg)=imp;
            break
		end
	end
end