function [carbon]=smiles_label_carbons2(str1)
% creates the stucture to describe the molecules
% usage: e.g. carbons=smiles_label_carbons2('C(CC(=O)C(O)=O)CO')
% would return e.g. carbons(1..n).br
% then carbons(i).br{ib} is the ib'th branch on carbon one
% if this can be broken down further then carbons(i).br{ib+1}.carbons(1...n)
% would contain the breakdown for that etc. etc.


[ncarbons]=smiles_ncarbons(str1); %counts the number of carbons in the main chain
        
%trying to strip down to the main chain
%so first of all remove all the bracketed items

nbr(1:ncarbons)=0; %an array to index the number of branches off each carbon

nleft=9999;
str_out=str1;
while nleft~=0    %keep going until have broken down all the brackets

    [str_out,str_in,nleft,nleft2,npos,npos2]=smiles_bracket(str_out); %break off the first complete bracket
    
    if npos~=0   %if had a bracket there then assign it to the correct carbon atom
        nbr(npos)=nbr(npos)+1;  %increase the branch index for this carbon
        carbon(npos).br{nbr(npos)}=str_in; %assign to the carbon
        
        if nleft2~=0  %if this has further brackets then break these down further with another call to this function
            nbr(npos)=nbr(npos)+1;   %increase branch index so that the chain is not overwritten (could remove this if don't want to display it)
            carbon(npos).br{nbr(npos)}.carbons=smiles_label_carbons2(str_in); %break down further
        end
    end
    
end

%now are left with the main chain (str_out = e.g. 'COCCOC')

%from this might want to check whether Os are in the chain or are at the end (acid)
%could be done within the for loop below

%now remove and assign the non-bracketed atoms (e.g. O atoms) from main chain

		iC=findstr('C',str_out);  %find position of all C atoms
        iall=[1:length(str_out)]; %index of all positions
        iall(iC)=[]; %i.e. remove indices that are C so are left with indices of non C atoms
        inonC=iall;  %duplicate     
        
        for i=1:length(inonC)  %go though all non C atoms and assign
            nc=inonC(i)-1;     %carbon atom is one before position of the non C atom
            nbr(nc)=nbr(nc)+1; %increase branch index
            carbon(nc).br{nbr(nc)}=str_out(iall(i)); %assign the non C atom
            inonC=inonC-1; %reduce the indices by one as have removed one non C atom
        end

        
    
    