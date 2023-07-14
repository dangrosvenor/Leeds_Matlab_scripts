filename_in = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/rose-app.conf';
filename_out = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/rose-app.conf_minus_upukca';

fid=fopen(filename_in,'rt');
fid_out = fopen(filename_out,'wt');

go=1;
i=0;
while go==1
    i=i+1;
    line = fgetl(fid);    
    if line==-1 %eof
        break
    end    
    
    if length(strfind(line,'[namelist:umstash_streq'))>0
        dat_temp{1}=line;
        for j=2:7
           line = fgetl(fid);         
           dat_temp{j}=line;
        end
        if length(strfind(line,'use_name=''UPUKCA'''))==0
            for j=1:7
                fprintf(fid_out,'%s\n',dat_temp{j});                
            end
        else
            '';
        end
                    
        
    else
       fprintf(fid_out,'%s\n',line); 
    end
                
    
end

fclose(fid);
fclose(fid_out);

fprintf(1,'\nFinished rm UPUKCA\n');