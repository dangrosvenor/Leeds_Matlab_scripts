function [str_0th_brack,str_1st_brack,nleft,nleft2,ipos,ipos_end,ipos_start]=smiles_bracket(str1)
%breaks off one complete bracket and returns what is outside (str_0th_brack) and inside (str_1st_brack) the bracket
%plus the number of brackets left outside (nleft) and inside the bracket (nleft2)
%and the position of the bracket in terms of carbons atoms in the main chain (ipos)
%and the position of the end of the brackets in terms of characters of full string (ipos_end)
%and the position of the end of the brackets in terms of characters of full string minus one (ipos_start)

nopen=findstr('(',str1);  %finds indices of opening brackets
nclose=findstr(')',str1); %and of closing ones

if length(nopen)==0  %if aren't any then return null values
    nleft=0;
    nleft2=0;
    str_0th_brack=str1;
    ipos=0;
    break
end

ipos_start=nopen(1)-1;
ipos_end=nclose(end);

before=str1(1:ipos_start);
iC=findstr('C',before); %indices of all C atoms before bracket       
ipos=length(iC); %position of the first branch in terms of C atoms (e.g. 1st carbon atom)

str_0th_brack=str1(1:ipos_start); %first bit of string outside the first opening bracket

 go=1;
 i=0;
 while go==1
    i=i+1;  %iterates through all the closing brackets
    a=length(find(nopen<nclose(i)));  %looking to see if the number opening brackets matches number of closing
    if a==i
        go=0;  %if it does then have a complete closed set of brackets and can output this
    end
    
end

str_0th_brack = [str_0th_brack str1(nclose(i)+1:end)];  %add the text to the right of the bracket

str_1st_brack = str1(nopen(1)+1:nclose(i)-1);  %the text inside the bracket

nleft=length( findstr('(',str_0th_brack) );   %determine how many brackets we have left to process
nleft2=length( findstr('(',str_1st_brack) );