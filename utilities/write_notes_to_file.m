function [filename_out]=write_notes_to_file(filename_in,str)
%Opens a text file and appends a string to it if it already exists (on a new line). If not
%then it creates a new file with a unique filename (using the date and
%time) and returns the name for future use

if exist(filename_in)==2 %if the file exists
    fid=fopen(filename_in,'at');
else
    filename_out = [filename_in '_' datestr(now,30) '.txt'];
    fid = fopen(filename_out,'wt');
end

fprintf(fid,'%s\n',str);
fclose(fid);
