function ncvar = ncl4matlab(nclfun,varargin)

% ========================================================================
% 
% By:      Gabe Kooperm, SIO, UCSD, 2011.
% Purpose: Use NCL functions in Matlab.
% Depends: NCAR Command Language (NCL).
% Notes:   1. The first argument 'nclfun' must be a string name of the
%             NCL function followed by the input arguments NCL expects.
%          2. True and False must be entered as strings.
% Example:
%      ZPL = ncl4matlab('vinth2p',Z,hyam,hybm,plevs,PS,2,1000,1,'False');
% 
% ========================================================================
path2ncl = '/Users/mikep/Library/ncl_6.0.0';
setenv('NCARG_ROOT',path2ncl); % Path to NCL

wkdir = '/Users/mikep/Desktop'; % Workspace

% ========================================================================

if exist(sprintf('%s/ncl4matlab.nc',wkdir),'file');
    unix(sprintf('rm %s/ncl4matlab.nc',wkdir));
end

if exist(sprintf('%s/ncl4matlab.ncl',wkdir),'file');
    unix(sprintf('rm %s/ncl4matlab.ncl',wkdir));
end

if exist(sprintf('%s/ncl4matlab2.nc',wkdir),'file');
    unix(sprintf('rm %s/ncl4matlab2.nc',wkdir));
end

ncid = netcdf.create(sprintf('%s/ncl4matlab.nc',wkdir),'CLOBBER');
for ii = 1:length(varargin);
    if iscellstr(varargin(ii));
        eval(sprintf('arg%02d = char(varargin(ii));',ii));
    elseif numel(cell2mat(varargin(ii))) == 1;
        eval(sprintf('arg%02d = num2str(cell2mat(varargin(ii)));',ii));
    else
        varin = cell2mat(varargin(ii));
        varin = permute(varin,ndims(varin):-1:1);
        dims = size(varin);
        kk = 0;
        for jj = 1:numel(dims);
            if dims(jj) > 1;
                kk = kk + 1;
                dimname = sprintf('arg%02ddim%02d',ii,kk);
                dimid(kk) = netcdf.defDim(ncid,dimname,dims(jj));
            end
        end
        varid = netcdf.defVar(ncid,sprintf('arg%02d',ii),'double',dimid);
        netcdf.endDef(ncid);
        netcdf.putVar(ncid,varid,varin);
        netcdf.reDef(ncid);
    end
    clear dimid varid;
end
netcdf.endDef(ncid);
netcdf.close(ncid);
clear ncid;

fid = fopen(sprintf('%s/ncl4matlab.ncl',wkdir),'w');
fprintf(fid,'begin\n');
fprintf(fid,'load "%s/lib/ncarg/nclscripts/csm/contributed.ncl"\n',path2ncl);
fprintf(fid,'load "%s/lib/ncarg/nclscripts/csm/gsn_code.ncl"\n',path2ncl);
fprintf(fid,'load "%s/lib/ncarg/nclscripts/csm/gsn_csm.ncl"\n',path2ncl);
fprintf(fid,'fin = addfile("%s/ncl4matlab.nc","r")\n',wkdir);
fprintf(fid,'fout = addfile("%s/ncl4matlab2.nc","c")\n',wkdir);
funstr = sprintf('ncvar=%s(',nclfun);
for ii = 1:length(varargin);
    if exist(sprintf('arg%02d',ii),'var') == 1;
        eval(sprintf('funstr = strcat(funstr,arg%02d);',ii));
        funstr = strcat(funstr,',');
    else
        fprintf(fid,'arg%02d=fin->arg%02d\n',ii,ii);        
        funstr = sprintf('%sarg%02d,',funstr,ii);
    end
end
funstr(end) = ')';
fprintf(fid,'%s\n',funstr);
fprintf(fid,'fout->ncvar=ncvar\n');
fprintf(fid,'end\n');
fclose(fid);

unix(sprintf('%s/bin/ncl %s/ncl4matlab.ncl',path2ncl,wkdir));
if ~exist(sprintf('%s/ncl4matlab2.nc',wkdir),'file');
    error('NCL script failed.');
end

ncid = netcdf.open(sprintf('%s/ncl4matlab2.nc',wkdir),'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'ncvar');
ncvar = netcdf.getVar(ncid,varid);
ncvar(ncvar == -9999) = nan;
ncvar = permute(ncvar,ndims(ncvar):-1:1);
netcdf.close(ncid);

unix(sprintf('rm %s/ncl4matlab.nc',wkdir));
unix(sprintf('rm %s/ncl4matlab.ncl',wkdir));
unix(sprintf('rm %s/ncl4matlab2.nc',wkdir));
