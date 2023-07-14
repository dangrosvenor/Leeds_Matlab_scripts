

% ivar=1; clear vars_to_save
%
% %mean cloud fraction with the lat/lon cell
% vars_to_save{ivar}='CF_mockL3'; ivar=ivar+1;
%
% vars_to_save{ivar}='Nice_mockL3'; ivar=ivar+1;
% vars_to_save{ivar}='Nreject_mockL3'; ivar=ivar+1;


        for ivar=1:length(vars_to_save)
            fprintf(1,'Saving %d of %d\n',ivar,length(vars_to_save));

            eval_str = [vars_to_save{ivar} tag '=' vars_to_save{ivar} ';'];
            eval(eval_str);

            if ~(exist(filename_savevars)==2)
                eval_str = ['save(filename_savevars,''' vars_to_save{ivar} tag ''',''-V7.3'');'];
            else
                eval_str = ['save(filename_savevars,''' vars_to_save{ivar} tag ''',''-V7.3'',''-APPEND'');'];
            end

            eval(eval_str);

        end

  














