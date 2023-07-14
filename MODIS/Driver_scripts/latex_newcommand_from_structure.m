function []=latex_newcommand_from_structure(dat,var_str,filename_save,iappend)

if iappend==0
    fid = fopen([filename_save '_STATS.tex'],'wt');
else
    fid = fopen([filename_save '_STATS.tex'],'at');
end

%Convert all of the variable names in the input structure to actual names
    %for ease of use
    name_struc='dat'; %The name of the structure
    names = eval(['fieldnames(' name_struc ');']);
    for i=1:length(names)
        eval_str = [names{i} ' = ' name_struc '.' names{i} ';']; 
        eval(eval_str);
    end
    
   for i=1:length(names)     
        if length(strmatch(names{i},'aux'))==0  %Don't print the values in the aux variable
            if length(strfind(names{i},'_precision'))>0
               continue %go to next iteration if it's the precision variable
            end
            str2 = remove_character(names{i},'_','');
            var_str = remove_character(var_str,'_','');
            var_str = remove_character(var_str,'-','');
            val = eval(names{i});
            %names{i};
            if exist([names{i} '_precision'])
                prec = eval([names{i} '_precision']);
                val_str = num2str(val,prec);
            else
                if abs(val)<1
                    val_str = num2str(val,'%.2f');
                else
                    val_str = num2str(val,'%.1f');
                end
            end
            fprintf(fid,['\\newcommand{\\' var_str str2 '}{' val_str '}\n']);
            %'\newcommand{\var_stat}{number}';
        end
    end
    
    fclose(fid);

