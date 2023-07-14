function [xar,yar,zar,zh,hclass,slatr,slongr,adate,ahhmm,latitude,longitude,radClass]=radarClassARMPaul(fnam)


radClass = {'unclassified','drizzle',...
        'rain','dry low density snow', ...
        'dry high density snow','melting snow', ...
        'dry graupel','wet graupel','small hail < 2cms',...
        'large hail > 2 cms','rain hail mix'};

fid = fopen(fnam,'r');

%
cline1  = fscanf(fid,'%14c',1);     
adate   = (cline1(1,1:8));   % date
ahhmm   = (cline1(1,10:13)); % time

%
cline2 = fscanf(fid,'%79c',1);
slatr = str2num(cline2(1,1:8));     % Radar latitude
slongr = str2num(cline2(1,10:18));  % Radar longitude
xmina = str2num(cline2(1,20:26));   % Minimum x value of radar grid
xmaxa = str2num(cline2(1,28:34));   % Maximum x value of radar grid
nx = str2num(cline2(36:38));        % Number of x points
delx=(xmaxa-xmina)/(nx-1);          % x resolution
ymina = str2num(cline2(1,40:46));   % Minimum y value of radar grid
ymaxa = str2num(cline2(1,48:54));   % Maximum y value of radar grid
ny = str2num(cline2(56:58));        % Number of y points
dely=(ymaxa-ymina)/(ny-1);          % y resolution
zmin = str2num(cline2(1,60:66));    % Minimum z value of radar grid
zmax = str2num(cline2(1,68:74));    % Maximum z value of radar grid
nz = str2num(cline2(76:78));        % Number of z points
delz=(zmax-zmin)/(nz-1);            % z resolution

% Read in all data.
vals = fscanf(fid,'%f ',[242, inf]);
fclose(fid);

% The data is arranged as follows x(zh,zclass),y,z
xar = linspace(xmina,xmaxa,nx);
yar = linspace(ymina,ymaxa,ny);
zar = linspace(zmin,zmax,nz);

% Now arrange the data
zh      = permute(reshape(vals(1:2:242,:),[nx ny nz]),[2 1 3]);
hclass  = permute(reshape(vals(2:2:242,:),[nx ny nz]),[2 1 3]);
zh(find(zh(:)==-99))=NaN;
hclass(find(hclass(:)==-9))=NaN;


latitude = slatr + yar./1.852./60;
Re=6378.1e3;
longitude = slongr + 1./(2.*pi.*Re.*cos(latitude.*pi./180)./360) .*xar.*1000;
% subplot(2,1,1)
% 
% pcolor(xar,yar,zh(:,:,iht)')
% 
% shading flat
% 
% colormap(jet)

