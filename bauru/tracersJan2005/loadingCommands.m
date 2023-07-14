%loading commands

savecom='load(exname';
for i=1:length(varlist)
    savestr{i}=[varlist{i} 'SAVE'];
	savecom=[savecom ',''' savestr{i} ''''];
end
savecom=[savecom ');'];        
eval(savecom);

for i=1:length(varlist)
    comm=[ varlist{i} '(jdir)=' varlist{i} 'SAVE;'];
    eval(comm);     
end
        
