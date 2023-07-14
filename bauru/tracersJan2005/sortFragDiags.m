idir=1;
nfiles=72;
names={'wind','vap','pressure','potemp','icemr','icenc','snowmr','snownc','graupelmr','graupelnc'};
members={'V';'v';'p';'p';'i';'i';'i';'i';'i';'i'};
members{1,2}='W';

[a b]=size(members);
altdir='c:/temp/250m/';

for inames=9:9 %:length(names)
    fprintf(1,'%d',inames);
    %clear( upper(names{inames}) );  %prob best to manually clear variables in case want to append to an array
    for ib=1:b 
       if length(members{inames,ib})>0
            for jj=1:nfiles
%                savename=[direcDan(idir).dir 'results/diags/' names{inames} '_' num2str(jj)];
                savename=[altdir names{inames} '_' num2str(jj)];

                load(savename);
                comm=[upper(names{inames}) '(1).' members{inames,ib} '(:,:,jj)=' names{inames} '(1).' members{inames,ib} ';']; 
                eval(comm);
            end
        
            savename=[direcDan(idir).dir 'results/diags/' upper(names{inames}) '_' members{inames,ib}];
            save(savename,upper(names{inames}));
            clear(upper(names{inames}));
        end
    end
end

'done'

