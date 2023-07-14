%Run the function below and put the variables in the current workspace

if ~exist('iadd_umid_label')
    iadd_umid_label=1;
end

[UM_case_out]=UM_case_select_runs(UM_cases,iadd_umid_label);
%Convert all of the variable names in the input structure to actual names
%for ease of use
names = fieldnames(UM_case_out);
for i=1:length(names)
    eval_str = [names{i} ' = UM_case_out.' names{i} ';'];
    eval(eval_str);
    eval_str2 = ['UM_map.' names{i} ' = ' names{i} ';'];
    eval(eval_str2);    
end