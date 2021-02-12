function [vars, starts ] = getVars(f)
str=fileread(f);

str=regexprep(str,"%[^\n]*(\n)?","$1"); %remove block comments
str=regexprep(str,'\.\.\.(\n|\r\n)',""); %remove elipses
str=regexprep(str,"\'[A-Za-z0-9]+\'",""); %remove simple hardcoded strings with single quotes
str=regexprep(str,"(?<=(?:\,| |\[|\(|\{))\'[^\'\n\r]*\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[A-Za-z0-9]+\"',""); %remove simple hardcoded strings with double quotes
str=regexprep(str,'\"[^\"\n\r]*\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition



ws=' \t\f';
code='\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';

name='([a-zA-Z][a-zA-Z_0-9]*)';
namelist=['((?:[ ]*' name '[ ]*\,[ ]*)?(?:[ ]*' name '[ ]*))'];
nameref=['(?:^|[' code '\n]+[' ws ']*|[' ws ']+|\n)' name ];

%remove local variables from anonymous functions as they are not defined in
%the scope of the script/function file
match_anon=regexp(str,['@\(?' namelist '\)?[^\n]+'],'match');
anon_local=regexp(str,['@\(?' namelist '\)?[^\n]+'],'tokens');
if ~isempty(anon_local)
    anon_local=regexp([anon_local{:}],'[ ]*\,[ ]*','split');
    anon_rep=cellfun(@(l,v) regexprep(l,strcat('(?<![A-Za-z0-9_])',v,'(?![A-Za-z0-9])'),repmat("",size(v))),match_anon,anon_local,'UniformOutput',false);
    for i=1:length(match_anon)
        str=strrep(str,match_anon{i},anon_rep{i});
    end
end

vars=regexp(str,nameref,"tokens"); %get variables
vars=unique([vars{:}],'stable');


if ~isempty(vars)
    
    bi=cellfun(@(x) exist(x,'builtin'),vars);
    fi= cellfun(@(x) exist(x,'file'),vars);
    mask = bi~=5 & (fi~=2 & fi~=6);
    
    %remove those lazy func calls with strings
    funcs=vars(~mask);
    funcs=funcs(cellfun(@(v) ~any(strcmp(v,{'if','while','for'})), funcs));%dont mistake if/while/for calls as a function
    funcs_match=cellfun(@(f) ['\n[' ws ']*' f '[' ws ']+[^\n]+\n'],funcs,'UniformOutput',0);
    
    str = regexprep(str,funcs_match,'\n');
    
    [vars,starts]=regexp(str,nameref,"tokens","start"); %get variables
    vars=[vars{:}];
    [~,iu]=unique(vars,'stable');
    vars=vars(iu);
    starts=starts(iu);
    if ~isempty(vars)
        
        bi=cellfun(@(x) exist(x,'builtin'),vars);
        fi= cellfun(@(x) exist(x,'file'),vars);
        mask = bi~=5 & (fi~=2 & fi~=6);
        
        vars=vars(mask);
        starts=starts(mask);
        
    end
    
    
    
    
    
    %     vars=unique([vars ],'stable');
    
    
    isgca=cellfun(@(x) strcmp(x,'gca'),vars);
    isgcf=cellfun(@(x) strcmp(x,'gcf'),vars);
    
    vars=vars(~( isgca | isgcf));
    
else
    vars={};
    
end
% mask = mask & cellfun(@(x) exist(x,'file')~=2 && exist(x,'file')~=6,vars);

% vars=unique([  vars_bare ]);% vars_bare
end

