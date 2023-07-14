function [DG,SERSTR,PARTITIONS]=SORTHEADER(HEADER)
% Function reads header and makes useful


% Allocate memory
K(1:20)=zeros(1,20);
RNDGS_START=findstr(char(HEADER),'RNDGS');
[tok,REM]=strtok(char(HEADER(RNDGS_START:end)),'&');

[TOK,REM]=strtok((REM(2:end)),'&');
%TOK(1:20)
% Now TOK contains all the RNDGS definitions
i=1;
REM=TOK;
PARTITIONS = {'ALL','ALu','ALd','CLu','CLd','ACu','ACd','BYu','BCu','PPd','MOu','BMu','M_1','M1u','M1d','M_2','M2u','M2d','M_3'};
doPart=zeros(1,19);
while(max(size(REM))>=2)
    [TOK,REM]=strtok((REM),',');
    char(TOK);
    switch(TOK(1:3))
    case {'ALL'}
        j=1;
        K(j)=K(j)+1;
        doPart(1) = 1;
    case {'ALu'}
        j=2;
        K(j)=K(j)+1;
        doPart(2) = 2;
    case {'ALd'}
        j=3;
        K(j)=K(j)+1;
        doPart(3) = 3;
    case {'CLu'}
        j=4;
        K(j)=K(j)+1;
        doPart(4) = 4;
    case {'CLd'}
        j=5;
        K(j)=K(j)+1;
        doPart(5) = 5;
    case {'ACu'}
        j=6;
        K(j)=K(j)+1;
        doPart(6) = 6;
    case {'ACd'}
        j=7;
        K(j)=K(j)+1;
        doPart(7) = 7;
    case {'BYu'}
        j=8;
        K(j)=K(j)+1;
        doPart(8) = 8;
    case {'BCu'}
        j=9;
        K(j)=K(j)+1;
        doPart(9) = 9;
    case {'PPd'}
        j=10;
        K(j)=K(j)+1;
        doPart(10) = 10;
    case {'MOu'}
        j=11;
        K(j)=K(j)+1;
        doPart(11) = 11;
    case {'BMu'}
        j=12;
        K(j)=K(j)+1;
        doPart(12) = 12;
    case {'M_1'}
        j=13;
        K(j)=K(j)+1;
        doPart(13) = 13;
    case {'M1u'}
        j=14;
        K(j)=K(j)+1;
        doPart(14) = 14;
    case {'M1d'}
        j=15;
        K(j)=K(j)+1;
        doPart(15) = 15;
    case {'M_2'}
        j=16;
        K(j)=K(j)+1;
        doPart(16) = 16;
    case {'M2u'}
        j=17;
        K(j)=K(j)+1;
        doPart(17) = 17;
    case {'M2d'}
        j=18;
        K(j)=K(j)+1;
        doPart(18) = 18;
    case {'M_3'}
        j=19;
        K(j)=K(j)+1;
        doPart(19) = 19;
    otherwise
        j=20;
        K(j)=K(j)+1;
    end     
    eval(sprintf('DG(%d).DG.a%s=%d;',j,TOK,i));
    i=i+1;
end
DG1 = DG(doPart(find(doPart)));  % sort the data
DG1(end+1) = DG(end);clear DG;
DG = DG1;

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
    char(TOK);
    if(~length(TOK))
        break;
    end
    SERSTR{i}=TOK;
    i=i+1;
end

PARTITIONS = PARTITIONS(doPart(find(doPart)));% sort the data
