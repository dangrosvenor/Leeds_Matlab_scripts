function [savestr,savecopy]=savevar(savestr,newvar)

savestr{end+1}=[newvar 'SAVE'];
comm=['savecopy=' newvar];
eval(comm);