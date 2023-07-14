lesfile='c:/cygwin/home/login/parallel_Bezier.f';
outfile='c:/cygwin/home/login/parallel_Bezier_done.f';

fid=fopen(lesfile);
fid2=fopen(outfile,'wt');

noread=0;
t='go';
while ischar(t)
    if noread==0
        t=fgetl(fid);
    else
        t=t2;
    end
    noread=0;
    
    f=findstr(upper(t),'CALL GC');
    if length(f)>0
        ff=findstr(t,'!');
        if length(ff)>0
            if ff(1)==1
                t(1)=' ';
            end
            fprintf(fid2,'%s\n',t);
            
            flag=1;
            while flag==1
                t2=fgetl(fid);
                ft2=findstr(t2,'&');
                if length(ft2)>0
                    if ft2(1)==7
                        fft2=findstr(t2,'!');
                        if length(fft2)>0
                            if fft2(1)==1
                                t2(1)='';
                            end
                            fprintf(fid2,'%s\n',t2);
                        end   
                    else
                        flag=0;
                        noread=1;
                    end
                else
                    flag=0;
                    noread=1;    
                end
            end
         
        else
            fprintf(fid2,'%s\n',t);
        end
    else
        fprintf(fid2,'%s\n',t);
    end
    
    
end
    

fclose(fid);
fclose(fid2);
    
'done'    