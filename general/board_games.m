%test run of Doodle for games night with weighting

%Could have an array that is 5 (no. days in week) * N_people * N_weeks
%0 could represent negative choice, 1 positive
%Also need a weighting array of length N_people*N_weeks (allows to keep
%track of weighting)

choices = NaN*ones([5 25 104]); %Have the default as =1 so that people who don't vote don't get their weight increased.
choice_override=0;
        
weights = zeros([25 104]);
weights(:,1)=1; %start the intial values at one for carrying through

peeps={'Dan     ','Joey    ','Tim     ','Chetan  ','Leighton','Robin   ','Luis    ','Ross    ','Douglas ','Ilkka   ','Peggy   '};
days={'Monday   ','Tuesday  ','Wednesday','Thursday ','Friday   '};

%Data entry :-

%% Historical data
maj_system=1; %set to one if were using the majority choice system (historical data)

%%



%%
iweek=2; %for week begnning 25th Jan - week 1 is the intial values

ipers=1; %Dan
choices(:,ipers,iweek) = [1 1 0 1 1];
ipers=2; %Joey
choices(:,ipers,iweek) = [0 0 1 0 1];
% ipers=3; %Tim
% choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=4; %Chetan
choices(:,ipers,iweek) = [0 0 0 0 1];
% ipers=5; %Leighton
% choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=6; %Robin
choices(:,ipers,iweek) = [0 0 0 0 1];
ipers=7; %Luis
choices(:,ipers,iweek) = [0 0 1 1 0];
% ipers=8; %Ross
% choices(:,ipers,iweek) = [1 1 1 1 1];
% ipers=9; %Douglas
% choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=10; %Ilkka
choices(:,ipers,iweek) = [0 0 0 0 1];

% Do the calculations
board_games_compute

%%
iweek=3; %for week begnning 1st Feb - week 1 is the intial values

ipers=1; %Dan
choices(:,ipers,iweek) = [0 1 1 1 1];
% ipers=2; %Joey
% choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=3; %Tim
choices(:,ipers,iweek) = [0 0 1 1 1];
ipers=4; %Chetan
choices(:,ipers,iweek) = [0 0 0 1 0];
ipers=5; %Leighton
choices(:,ipers,iweek) = [0 0 0 1 1];
ipers=6; %Robin
choices(:,ipers,iweek) = [1 1 0 1 1];
ipers=7; %Luis
choices(:,ipers,iweek) = [0 1 1 1 0];
ipers=8; %Ross
choices(:,ipers,iweek) = [0 1 0 1 1];
ipers=9; %Douglas
choices(:,ipers,iweek) = [1 1 1 1 0];

% Do the calculations
board_games_compute

%%
iweek=4; %for week begnning 8th Feb

ipers=1; %Dan
choices(:,ipers,iweek) = [1 0 0 1 1];
ipers=2; %Joey
choices(:,ipers,iweek) = [1 0 0 0 0];
ipers=3; %Tim
choices(:,ipers,iweek) = [1 1 1 0 0];
ipers=4; %Chetan
choices(:,ipers,iweek) = [0 0 1 1 0];
ipers=5; %Leighton
choices(:,ipers,iweek) = [0 0 0 1 1];
ipers=6; %Robin
choices(:,ipers,iweek) = [1 0 1 1 1];
ipers=7; %Luis
choices(:,ipers,iweek) = [0 1 0 1 0];


% Do the calculations
board_games_compute


%%
iweek=5; %for week begnning 15th Feb - week 1 is the intial values

% ----- weights started
maj_system=0; %Using weights


ipers=1; %Dan
choices(:,ipers,iweek) = [1 1 0 1 1];
%ipers=2; %Joey
%choices(:,ipers,iweek) = [1 1 1 0 1];
ipers=3; %Tim
choices(:,ipers,iweek) = [1 1 0 1 0];
% ipers=4; %Chetan
% choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=5; %Leighton
choices(:,ipers,iweek) = [0 0 0 0 1];
ipers=6; %Robin
choices(:,ipers,iweek) = [1 1 1 1 1];
% ipers=7; %Luis
% choices(:,ipers,iweek) = [1 1 1 1 1];

choice_override = 4; %Thursday was chosen

% Do the calculations
board_games_compute


%%
iweek=iweek+1; %for week begnning 22nd Feb

% ----- weights started
maj_system=0; %Using weights


ipers=1; %Dan
choices(:,ipers,iweek) = [1 1 0 1 1];
ipers=2; %Joey
choices(:,ipers,iweek) = [1 1 0 0 1];
ipers=3; %Tim
choices(:,ipers,iweek) = [1 0 1 0 0];
ipers=4; %Chetan
choices(:,ipers,iweek) = [0 0 0 0 1];
ipers=5; %Leighton
choices(:,ipers,iweek) = [0 0 0 1 1];
ipers=6; %Robin
choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=7; %Luis
choices(:,ipers,iweek) = [0 0 0 0 0];
% ipers=8; %Ross
% choices(:,ipers,iweek) = [1 1 1 1 1];
ipers=9; %Douglas
choices(:,ipers,iweek) = [0 1 1 0 0];
%ipers=10; %Ilkka
%choices(:,ipers,iweek) = [0 0 0 0 1];
ipers=11; %Peggy
choices(:,ipers,iweek) = [0 1 0 0 0];

%choice_override = 4; %Thursday was chosen

% Do the calculations
board_games_compute


%%
iweek=iweek+1; %for week begnning 29th Feb

% ----- weights started
maj_system=0; %Using weights


ipers=1; %Dan
choices(:,ipers,iweek) = [1 1 0 1 1];
ipers=2; %Joey
choices(:,ipers,iweek) = [1 1 1 0 0];
%ipers=3; %Tim
%choices(:,ipers,iweek) = [1 0 1 0 0];
%ipers=4; %Chetan
%choices(:,ipers,iweek) = [0 0 0 0 1];
ipers=5; %Leighton
choices(:,ipers,iweek) = [0 0 0 1 0];
ipers=6; %Robin
choices(:,ipers,iweek) = [1 1 1 1 1];
%ipers=7; %Luis
%choices(:,ipers,iweek) = [0 0 0 0 0];
% ipers=8; %Ross
% choices(:,ipers,iweek) = [1 1 1 1 1];
%ipers=9; %Douglas
%choices(:,ipers,iweek) = [0 1 1 0 0];
%ipers=10; %Ilkka
%choices(:,ipers,iweek) = [0 0 0 0 1];
%ipers=11; %Peggy
%choices(:,ipers,iweek) = [0 1 0 0 0];

choice_override = 1; %Chose Monday

% Do the calculations
board_games_compute


%%
iweek=iweek+1; %for week begnning 7th March

% ----- weights started
maj_system=0; %Using weights


ipers=1; choices(:,ipers,iweek) = [1 1 0 1 0]; %Dan
ipers=2; choices(:,ipers,iweek) = [1 1 1 0 0]; %Joey
%ipers=3; choices(:,ipers,iweek) = [1 0 1 0 0]; %Tim
%ipers=4; choices(:,ipers,iweek) = [0 0 0 0 1]; %Chetan
ipers=5; choices(:,ipers,iweek) = [0 0 1 0 0]; %Leighton
ipers=6; choices(:,ipers,iweek) = [1 0 1 1 1]; %Robin
ipers=7; choices(:,ipers,iweek) = [0 0 1 1 0]; %Luis
%ipers=8; choices(:,ipers,iweek) = [1 1 1 1 1]; %Ross
%ipers=9; choices(:,ipers,iweek) = [0 1 1 0 0]; %Douglas
%ipers=10; choices(:,ipers,iweek) = [0 0 0 0 1]; %Ilkka
ipers=11; choices(:,ipers,iweek) = [0 0 1 0 0]; %Peggy

%choice_override = 1; 

% Do the calculations
board_games_compute


%%
iweek=iweek+1; %for week begnning 14th March

% ----- weights started
maj_system=0; %Using weights


%ipers=1; choices(:,ipers,iweek) = [1 1 0 1 0]; %Dan
%ipers=2; choices(:,ipers,iweek) = [1 1 1 0 0]; %Joey
%ipers=3; choices(:,ipers,iweek) = [1 0 1 0 0]; %Tim
ipers=4; choices(:,ipers,iweek) = [0 0 0 0 0]; %Chetan
ipers=5; choices(:,ipers,iweek) = [0 0 0 0 0]; %Leighton
ipers=6; choices(:,ipers,iweek) = [0 0 0 0 0]; %Robin
ipers=7; choices(:,ipers,iweek) = [0 0 0 0 0]; %Luis
%ipers=8; choices(:,ipers,iweek) = [1 1 1 1 1]; %Ross
%ipers=9; choices(:,ipers,iweek) = [0 1 1 0 0]; %Douglas
%ipers=10; choices(:,ipers,iweek) = [0 0 0 0 1]; %Ilkka
%ipers=11; choices(:,ipers,iweek) = [0 0 1 0 0]; %Peggy

choice_override = 1;  %manually choose the day when there is a tie

% Do the calculations
board_games_compute