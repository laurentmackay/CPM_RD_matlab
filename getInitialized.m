function assigned = getVars(f)
str=fileread(f);


str=regexprep(str,"%[ \w\=\(\)\+\;\:\_\.\*\,\]\[\-]+[\n\r ]",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

assigned=regexp(str,'([a-zA-Z_$][a-zA-Z_$0-9]*) *\=',"tokens"); %get variables
assigned=cellfun(@(x) x(1), assigned);
assigned=unique(assigned);

mask = cellfun(@(x) exist(x,'builtin')~=5,assigned);
mask = mask & cellfun(@(x) exist(x,'file')~=2,assigned);

assigned=assigned(mask);
mask=true(1,length(assigned));
%remove assignments that reference the variable
for i = 1:length(assigned)
    var=assigned{i};
    assignments=regexp(str,['[ \s]+' var ' *\=([^\n\r\;\=]+)[\n\r\;]{1}'],"tokens");
    assignments=cellfun(@(x) x(1), assignments);
    for ass = assignments
        numerical = regexp(ass,' ?[0-9]\.?[0-9]?');
        if ~isempty(numerical{1}) & numerical{1}==1
            break; %do not worry about variables that are assigned a numerical value at the beginning of the file
        end
        update=regexp(ass,[var '[^\(]'],'match');
        if ~isempty(update{1})
            mask(i)=false;
        end
    end
end
assigned=assigned(mask);
end

