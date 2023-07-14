%smiles parser

str1='CC(C(CC(=O)C(O)=O)CO)(C)O';

[str_out,str_in,nleft,nleft2]=smiles_bracket(str1)




%could have nested structure like this - but would be quite complicated to actually disseminate!
br{1}='CC(C(CC(=O)C(O)=O)CO)(C)O';

[a,b]=smiles_bracket(br{1});
br{1}.br{1}=a;
br{1}.br{2}=b;

[a21,b21]=smiles_bracket(a);
br{1}.br{1}.br{1}=a21;
br{1}.br{1}.br{2}=b21;

[a22,b22]=smiles_bracket(b);
br{1}.br{2}.br{1}=a22;
br{1}.br{2}.br{2}=b22;

%perhaps better to branch from the carbon chain

carbons(1).br{1}='';
carbons(2).br{1}='C(CC(=O)C(O)=O)CO'; %=str_in - would then break this down more
carbons(2).br{2}='C';
carbons(2).br{3}='O';

carbons(2).br{1}.chain='CCO';  %i.e. have two C atoms that need to label
carbons(2).br{1}.chain.carbons(1).br{1}='CC(=O)C(O)=O'; %would become CCC(O)=O and =O
carbons(2).br{1}.chain.carbons(2).br{1}='O'; 

carbons(2).br{1}.br{2}='CC(=O)C(O)=O';


str_out2='CCCC(OCC)(CCC)C(C)CC';  %say that had this carbon atom
%in order to count the number of C atoms in the chain can do something like

nleft=9999;
while nleft~=0
    [str_out2,str_in2,nleft,nleft2]=smiles_bracket(str_out2); %basically keeps reducing down until are left with a chain of standalone carbons
end

ncarbons=length(findstr('C',str_out2));

%or function written to do this

% ncarbons=smiles_ncarbons(str_out2);
% so would then know how many carbon atoms to label for the chain


