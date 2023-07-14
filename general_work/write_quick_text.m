function write_quick_text(text)

fid=fopen('/home/disk/eos1/d.grosvenor/text.txt','wt');
fprintf(fid,text);
fclose(fid);