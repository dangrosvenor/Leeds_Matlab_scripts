%function to read in matrix of values from given file - input the first few data points so can find beginning and also width of table


function diracd=readtable(pat,firstvals,nocols)

fid=fopen(pat,'rt');

 cont=1;
 while (cont==1)
     for i=1:length(firstvals)
         inum(i)=0;
         if isstr(firstvals(i)) | iscell(firstvals(i))
             temp=fscanf(fid,'%s',[1]);
             if length(findstr(char(firstvals(i)),temp))>0; 
                 cont=0; 
             else
                 cont=1;
                 break
                 %fscanf(fid,'%g',[1]);
             end
         else
             inum(i)=1;
            temp=fscanf(fid,'%g',[1]);
            if abs(temp-firstvals(i))<1e-10;
                diracd(i,1)=temp;
                cont=0;
            else
                cont=1;
                %break
                if length(temp)==0;
                    fscanf(fid,'%s',[1]);
                end
            end
         end
    end
 end
 
i=sum(inum);
 if nocols>1
     diracd(i+1:nocols,1)=fscanf(fid,'%g',[nocols-i]);
 end

 dirac2=fscanf(fid,'%g',[nocols Inf]);
 diracd(:,2:size(dirac2,2)+1)=dirac2;    
 fclose(fid);
 


