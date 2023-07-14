function ctick_text=ctick_strings_make(cticks)

for j=1:length(cticks)
    te=num2str(cticks(j));    
    ctick_text(j,1:length(te))=te;    
end
