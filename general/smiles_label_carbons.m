function [text]=smiles_label_carbons(str1)
[ncarbons,posC,chain]=smiles_ncarbons(str1);

[str_out,str_in,nleft,nleft2,npos]=smiles_bracket(str1);        
    nleft2_save=nleft2;
    str_in_save=str_in;
    npos_save=npos;

%chain is the main chain without the bracketed items
%sort out non-bracketed atoms from main chain

		iC=findstr('C',chain);
        iall=[1:length(chain)];
        iall(iC)=[]; %i.e. indices that aren't C
        inonC=iall;
        
%        str1(inonC)='';
        
        for i=1:length(iall)
            carbon(iall(i)-1).br{1}=chain(inonC(i)); %set as the non C atom
            iall=iall-1; %reduce the indices by one as have removed one non C atom
        end
        

%npos is the position of the bracket in terms of C atoms     

for icarb=1:ncarbons
    
    if npos==icarb    
        carbon(icarb).br{2}=str_in;
        %carbon(icarb).br{1}=smiles_label_carbons(str_in);
        nleft_save=nleft;
        for ileft=1:nleft_save;
            [str_out,str_in,nleft,nleft2,npos]=smiles_bracket(str_out);
            carbon(icarb).br{ileft+2}=str_in;
        end    
        

            
                        
    else
        carbon(icarb).br{1}='';
    end
    
end

        
        
        if nleft2_save>0
            carbon(npos_save).br{1}.br=smiles_label_carbons(str_in_save);
        end        
        

text=carbon;
    
    