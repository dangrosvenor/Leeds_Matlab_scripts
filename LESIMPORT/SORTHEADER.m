function [DG,SERSTR]=SORTHEADER(HEADER)
% Function reads header and makes useful


% Allocate memory
K(1:20)=zeros(1,20);
RNDGS_START=findstr(char(HEADER),'RNDGS');
[tok,REM]=strtok(char(HEADER(RNDGS_START:end)),'&');
REM=REM(2:end);
[TOK,REM]=strtok((REM(2:end)),'&');
%TOK(1:20)
% Now TOK contains all the RNDGS definitions

i=1;
REM=TOK(1:end-1);
while(max(size(REM))>1)
    [TOK,REM]=strtok((REM),',');
    switch(TOK(1:3))
    case {'ALL'}
        j=1;
        K(j)=K(j)+1;
    case {'ALu'}
        j=2;
        K(j)=K(j)+1;
    case {'ALd'}
        j=3;
        K(j)=K(j)+1;
    case {'CLu'}
        j=4;
        K(j)=K(j)+1;
    case {'CLd'}
        j=5;
        K(j)=K(j)+1;
    case {'ACu'}
        j=6;
        K(j)=K(j)+1;
    case {'ACd'}
        j=7;
        K(j)=K(j)+1;
    case {'BYu'}
        j=8;
        K(j)=K(j)+1;
    case {'BCu'}
        j=9;
        K(j)=K(j)+1;
    case {'PPd'}
        j=10;
        K(j)=K(j)+1;
    case {'MOu'}
        j=11;
        K(j)=K(j)+1;
    case {'BMu'}
        j=12;
        K(j)=K(j)+1;
    case {'M_1'}
        j=13;
        K(j)=K(j)+1;
    case {'M1u'}
        j=14;
        K(j)=K(j)+1;
    case {'M1d'}
        j=15;
        K(j)=K(j)+1;
    case {'M_2'}
        j=16;
        K(j)=K(j)+1;
    case {'M2u'}
        j=17;
        K(j)=K(j)+1;
    case {'M2d'}
        j=18;
        K(j)=K(j)+1;
    case {'M_3'}
        j=19;
        K(j)=K(j)+1;
    otherwise
        j=20;
        K(j)=K(j)+1;
    end     
    %DGAVSTR(j).DG{K(j)}=TOK;
    %fprintf(1,'DG(%d).DG.a%s=%d',j,TOK,i);
    %eval(sprintf('DG(%d).DG.a%s=%d',j,TOK,i));
    INDS=find(TOK~=' ');
    DG{i}=TOK(INDS);
    i=i+1;
end
SER_START=findstr(char(HEADER),'TIMES');
[tok,REM]=strtok(char(HEADER(SER_START:end)),'&');

[TOK,REM]=strtok((REM(1:end)),',');
TOK=TOK(1:end-4);
%TOK(1:20)
% Now TOK contains all the SER definitions
SERSTR{1}='TIMES';
i=2;
REM=TOK;
while(REM)
    [TOK,REM]=strtok((REM),'&');
    SERSTR{i}=TOK;
    %fprintf(1,'DG.a%s=%d',DGAVSTR{i},i);
    %eval(sprintf('DG.a%s=%d',DGAVSTR{i},i));
    i=i+1;
end
