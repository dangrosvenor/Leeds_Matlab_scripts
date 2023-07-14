%sorts header out so that lists the order of the diagnotic profiles in dgstr

clear dgstr dgstr2 rem
a=findstr('RNDGS&',char(HEADER2));
b=findstr('&',char(HEADER2(a:end)));
bend=a+b(2);

i=1;
j=1;
flag=0;

[st rem]=strtok(char(HEADER2(a+9:end)),','); %finds next 'word'
dgstr{1}=st;
dgstr2{1}=st;
%while ( strcmp(st(1),'&') ~=1 )
while (size(st,1)>0)
    x=findstr(st,'ALL');
    
    if size(x,2)>0
        flag=flag+1;
        j=j+1;
        dgstr2{j}=st;
    end
    
    if (flag==0)
        j=j+1;
        dgstr2{j}=st;
    end
        
    i=i+1;
    [st rem]=strtok(rem,',');
    
    x=findstr(st,'&');
    if size(x,1)==0
        dgstr{i}=st;
    else
        %dgstr{i}=strtok(st,'&');
        %st2=st;
        %st='';
        
    end
    
end

clear rem

    
%see findhead to find specific 'word' within dgstr