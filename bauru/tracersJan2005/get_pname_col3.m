function imphys=get_pname_col3(dgcell,idir,pname_all)
%imphys=get_pname_col2(dgcell,idir,pname_all)
%dgecell can be more than one string - e.g. = {'PRAUT','PGDEP'}

            
for idg=1:length(dgcell)            
	for imp=1:length(pname_all(idir).pname)
		if strcmp(upper(pname_all(idir).pname(imp).p),upper(dgcell{idg}))==1
            imphys(idg)=imp;
            break
		end
	end
end