           
function [stoic,chems] = str2stoic(str)

name='[a-zA-Z_$][a-zA-Z_$0-9:\-]*';
empty_or_num='[ \f\t\v\.0-9]*';

split = regexp(str,['(' empty_or_num ')(' name ')'],'tokens');
split=[split{:}];
stoic = cellfun(@(x) empty2one(str2double(x{1})), split);
chems = cellfun(@(x) x{2}, split,'UniformOutput',0);

end