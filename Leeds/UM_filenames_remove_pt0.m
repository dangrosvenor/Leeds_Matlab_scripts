

   for idat_UM=1:length(fileUM)
%for idat_UM=1:length(fileUM)
    idat_driver = idat_driver+1;
    
    if iscell(dirUM)==1
        dirUM_i = dirUM{idat_UM};
    else
        dirUM_i = dirUM;
    end 
    
     filename = [dirUM_i fileUM{idat_UM}];
     
     ifind = strfind(filename,'/');
     dir_name = filename(1:ifind(end));
    
    
    search_str=remove_character(filename,'VAR_NAME','*');
    files=dir([search_str '*.txt']);
    
    for ifile=1:length(files)
       
        old_name = [dir_name files(ifile).name];
        ifind = strfind(files(ifile).name,'z1500.0');
        if length(ifind)>0
            new_name = [dir_name remove_character(files(ifile).name,'z1500.0','z1500')];
            eval_str = ['!mv ' old_name ' ' new_name];
            eval(eval_str);
        end
        
        
    end
    
   end