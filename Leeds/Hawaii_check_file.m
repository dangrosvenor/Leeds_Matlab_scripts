
    for ifile=1:length(files)
        i=0;
        for iex=1:length(exceptions)                    
           if length(strfind(files(ifile).name,exceptions{iex})) == 0
               i = i + 1;
           else
               break
           end
        end
        if i == length(exceptions)
            break
        end
    end
    