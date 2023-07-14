infile='c:/cygwin/home/user/emm/EMM1-dan2.cpp';
outfile='c:/cygwin/home/user/emm/EMM_nofprintf.cpp';

fid=fopen(infile);
fid2=fopen(outfile,'wt');

noread=0;
t='go'; %any char
ifflag=0;
printif=1;
count=0;
while ischar(t)
    if noread==0
        t=fgetl(fid);
        count=count+1;
    else
        t=t2;
    end
    noread=0;
    
    if count==2558
        '';
    end
    
    f=findstr(lower(t),'fprintf');
    f2=findstr(lower(t),'fopen');
    f3=findstr(lower(t),'fclose');
    f4=findstr(lower(t),'if(');
    f44=findstr(lower(t),'if (');
    f5=findstr(lower(t),'{');
    f6=findstr(lower(t),'else');
    f7=findstr(lower(t),'for');
    
    if (length(f44)>0 & length(f5)==0) | (length(f4)>0 & length(f5)==0) | (length(f6)>0 & length(f5)==0) | (length(f7)>0 & length(f5)==0 ) %if there is an if with no curly brackets on previous line
        if ifflag==1
            fprintf(fid2,'%s\n',told); %in case have two ifs in a row
        end
        ifflag=1; %then will need to delete the if too (or could perhaps replace with a printf or something
        told=t; %store line for possible later printing
        printif=0;
    end
    
    if length(f)>0 | length(f2)>0 | length(f3)>0
        ifflag=0; %reset ifflag so next line can be printed
        
        ff=findstr(t,'{'); 
        ff2=findstr(t,'}'); 
        if length(ff)>0 & lenght(ff2)==0
            fprintf(fid2,'%s\n','{');
        elseif length(ff)>0 & lenght(ff2>0)
            fprintf(fid2,'%s\n','{}');
        end
        
        
            
        
        ff=findstr(t,')'); 
        if length(ff)==0  %if the matching bracket is not found then look for it    
            flag=1;
            while flag==1
                t2=fgetl(fid);
                ft2=findstr(t2,')');
                if length(ft2)>0
                    ff=findstr(t,'}'); 
                    if length(ff)>0
                        fprintf(fid2,'%s\n','}');
                    end
                    flag=0; %if we find it then stop looking 
                end
            end
        end
        
    else
        if ifflag==0
            fprintf(fid2,'%s\n',t);
        elseif printif==1
            fprintf(fid2,'%s\n',told);
            fprintf(fid2,'%s\n',t);
            ifflag=0;
        end
        printif=1;
    end
    
    
end
    

fclose(fid);
fclose(fid2);
    
'done'  








