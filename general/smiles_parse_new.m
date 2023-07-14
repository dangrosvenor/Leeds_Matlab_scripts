function [base]=smiles_parse_new(str1)
% creates the stucture to describe the molecules
% usage: e.g. bases=smiles_label_bases2('C(CC(=O)C(O)=O)CO')
% would return e.g. bases(1..n).br
% then bases(i).br{ib} is the ib'th branch on base one
% if this can be broken down further then bases(i).br{ib+1}.bases(1...n)
% would contain the breakdown for that etc. etc.


%[nbase]=smiles_nbase(str1); %counts the number of bases in the main chain
 
str_out=str1;
nleft=9999;
while nleft~=0
    [str_out,str_in,nleft,nleft2,npos,npos_start]=smiles_bracket(str_out); 
    %basically keeps removing complete brackets until are left with a chain of standalone carbons and other non-bracketed atoms (e.g. 'COCO')
end

nbase=length(str_out);
nbr(1:nbase)=0; %an array to index the number of branches off each base

if strcmp(str_out(end-1),'=')==1
    nbase=nbase-2;
end

for i=1:nbase
    base(i).atom=str_out(i); %label all the base atoms
end


if strcmp(str_out(end-1),'=')==1
    nbr(nbase)=nbr(nbase)+1;
    base(nbase).br{nbr(nbase)}=str_out(end-1:end);
    str_out(end-1:end)=[];
end

%ieq=findstr('=',str_out);
%str_out(ieq)=[];
 

%function [str_0th_brack,str_1st_brack,nleft,nleft2,ipos,ipos_end,ipos_start]=smiles_bracket(str1)
 
%trying to strip down to the main chain
%so first of all remove all the bracketed items



nleft=9999;
str_out=str1;
while nleft~=0    %keep going until have broken down all the brackets

    [str_out,str_in,nleft,nleft2,npos_C,npos2,npos]=smiles_bracket(str_out); %break off the first complete bracket   
    
    if npos~=0   %if had a bracket there then assign it to the correct base chain atom
        nbr(npos)=nbr(npos)+1;  %increase the branch index for this base atom
        base(npos).br{nbr(npos)}=str_in; %assign to the base atom 
        
        if nleft2~=0  %if this has further brackets then break these down further with another call to this function
            %nbr(npos)=nbr(npos)+1;   %increase branch index so that the chain is not overwritten (could remove this if don't want to display it)
            base(npos).next{nbr(npos)}.base=smiles_parse_new(str_in); %break down further
        end
    end
    
end

%now are left with the main chain (str_out = e.g. 'COCCOC')

%from this might want to check whether Os are in the chain or are at the end (acid)
%could be done within the for loop below

%now remove and assign the non-bracketed atoms (e.g. O atoms) from main chain
% 
% 		iC=findstr('C',str_out);  %find position of all C atoms
%         iall=[1:length(str_out)]; %index of all positions
%         iall(iC)=[]; %i.e. remove indices that are C so are left with indices of non C atoms
%         inonC=iall;  %duplicate     
%         
%         for i=1:length(inonC)  %go though all non C atoms and assign
%             nc=inonC(i)-1;     %base atom is one before position of the non C atom
%             nbr(nc)=nbr(nc)+1; %increase branch index
%             base(nc).br{nbr(nc)}=str_out(iall(i)); %assign the non C atom
%             inonC=inonC-1; %reduce the indices by one as have removed one non C atom
%         end

        
    
    