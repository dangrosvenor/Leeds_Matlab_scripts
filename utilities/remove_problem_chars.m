function out = remove_problem_chars(in)

out=in;

out=remove_character(out,'*','');
out=remove_character(out,'\','');
out=remove_character(out,' ','_'); %replace all spaces with underscores - latex can handle single spaces
out=remove_character(out,'<','.LT.'); 
out=remove_character(out,'>','.GT.'); 
out=remove_character(out,'%','pct');