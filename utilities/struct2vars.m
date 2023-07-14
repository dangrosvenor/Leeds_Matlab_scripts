function struct2vars(struct)

names = fieldnames(struct);
for i=1:length(names)
    eval_str = [names{i} ' = struct.' names{i} ';'];
    eval(eval_str);
end