UM_var = 'LWP';
UM_var = 'LowCF';
UM_var = 'LWP_Adrian';

idat=0;

filenames={...
    'xkjksa_pg200811120030_LWP_Adrian_TIMSER.txt'...    
    'xkjkda_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkjkwa_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkjkga_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkjkha_pg200811120030_LWP_Adrian_TIMSER.txt'...      
    'xkjkia_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkjkka_pg200811120030_LWP_Adrian_TIMSER.txt'...    
    };

filenames={...
    'xkjksa_pg200811120030_LWP_Adrian_TIMSER.txt'...    
    'xkjkga_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkcvoa_pd200811120000_LWP_TIMSER.txt'...    
    'xkjkia_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkcvva_pd200811120000_LWP_TIMSER.txt'...        
    };

%Effect of varying aerosol
filenames={...
    'xkjkia_pg200811120030_LWP_Adrian_TIMSER.txt'...
    'xkjkxa_pg200811120030_LWP_Adrian_TIMSER.txt'...    
    'xkjkya_pg200811120030_LWP_Adrian_TIMSER.txt'...     
    };

filedir = '~/UM/xkjk/';


filenames={...
    'xkqkfa_pd200810251200_LWP_TIMSER.txt'...
    'xkqkha_pd200810251259_LWP_TIMSER.txt'...   
    'xkqkja_pd200810251259_LWP_TIMSER.txt'...
    };

filedir = '~/UM/xkqk/';


for idat=1:length(filenames)
    
    %filename=['~/UM/xjqcha_pg0000.pp_' UM_var '_TIMSER.txt'];
    filename=[filedir filenames{idat}];
    fid=fopen(filename,'rt');
    tmp=fgetl(fid);
    dat=textscan(fid,'%f %f');
    fclose(fid);


    time_UM{idat} = dat{1}(:);
    eval([UM_var '_UM{idat} = dat{2}(:);']);
    filename_UM{idat} = filename;

end

fprintf(1,'\n Done UM read\n');

% % new run starts here
% filename=['~/UM/xjoiha_pg0000.pp_' UM_var '_TIMSER.txt'];
% fid=fopen(filename,'rt');
% tmp=fgetl(fid);
% dat=textscan(fid,'%f %f');
% fclose(fid);
% 
% idat=idat+1;
% time_UM{idat} = dat{1}(:);
% eval([UM_var '_UM{idat} = dat{2}(:);']);
% filename_UM{idat} = filename;
% 
% % new run starts here
% filename=['~/UM/xjoiia_pg0000.pp_' UM_var '_TIMSER.txt'];
% fid=fopen(filename,'rt');
% tmp=fgetl(fid);
% dat=textscan(fid,'%f %f');
% fclose(fid);
% 
% idat=idat+1;
% time_UM{idat} = dat{1}(:);
% eval([UM_var '_UM{idat} = dat{2}(:);']);
% filename_UM{idat} = filename;
% 
% 
% % new run starts here
% filename=['~/UM/xjqlha_pg0000.pp_' UM_var '_TIMSER.txt'];
% fid=fopen(filename,'rt');
% tmp=fgetl(fid);
% dat=textscan(fid,'%f %f');
% fclose(fid);
% 
% idat=idat+1;
% time_UM{idat} = dat{1}(:);
% eval([UM_var '_UM{idat} = dat{2}(:);']);
% filename_UM{idat} = filename;
% 
% 
% % new run starts here
% filename=['~/UM/xjqcka_pg0000.pp_' UM_var '_TIMSER.txt'];
% fid=fopen(filename,'rt');
% tmp=fgetl(fid);
% dat=textscan(fid,'%f %f');
% fclose(fid);
% 
% idat=idat+1;
% time_UM{idat} = dat{1}(:);
% eval([UM_var '_UM{idat} = dat{2}(:);']);
% filename_UM{idat} = filename;