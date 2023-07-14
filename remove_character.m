function str_out=remove_character(str_in,char_to_remove,str_replace)
% function str_out=remove_character(str_in,char_to_remove,str_replace)
% Replaces all instances of char_to_remove from str_in and replaces with str_replace.
% str_replace can be longer than the number of the characters removed.

LS = length(char_to_remove);
str_orig=str_in;
itext=findstr(char_to_remove,str_in);
itext=[1-LS itext];
str_out='';

if nargin<3
    str_replace='';
end

for i=2:length(itext)
    str_out = [str_out str_orig(itext(i-1)+LS:itext(i)-1) str_replace];
end
str_out = [str_out str_orig(itext(end)+LS:end)];


