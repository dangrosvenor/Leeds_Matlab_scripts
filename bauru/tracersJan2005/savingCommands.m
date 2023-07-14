%saving commands
savecom='save(exname';

for i=1:length(varlist)
        savestr{i}=[varlist{i} 'SAVE'];
        comm=[ varlist{i} 'SAVE=' varlist{i} '(jdir);'];
        eval(comm);     
		savecom=[savecom ',''' savestr{i} '''']  ;      
    end

savecom=[savecom ');'];

%fprintf(1,'\n%s file=%s\n',savecom,exname);
%qans=input('You sure you want to write to file? (type ''y'' with quotes)');

% if strcmp('y',qans)
 	eval(savecom);
%     fprintf(1,'*** files saved ***');
% else
%     fprintf(1,'*** files NOT saved ***');
% end
