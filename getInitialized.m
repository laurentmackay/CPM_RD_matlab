function assigned = getInitialized(f, ignore_code_ref)
str=fileread(f);

if nargin==1
    ignore_code_ref=false;
end


str=regexprep(str,"%[^\n]*(\n)?","$1"); %remove comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,'~([^=])','$1'); %remove tildes
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition
str=regexprep(str,"\'",""); %remove transposes

code='0-9 \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\:\|\{\}\^\.\@';
name='[a-zA-Z_\^][a-zA-Z_$0-9]*';
assignment='(?:[ \t\f]*|[a-zA-Z_$0-9])\=[ \t\f]*';
% assignment='(?<!~)\=[ \t\f]*';
% assignment='(?:([ \t\f]*)\=[ \t\f]*)';
seps='(:?\,| )';



list_lhs=['((:?' seps '*' name '[ \t\f]*?' seps ')+(:?' name '[ ]*))'];
% list_rhs=['((:?[' code ']*' ref '[' code ']*(?:' code ')*)?(:?[' code ']*' ref '[' code ']*))'];
list_rhs=['((:?[' code ']*' name '[' code ']+)*(:?[' code ']*' name '[' code ']*?))'];

% assigned=regexp(str,['('  name  ') *\=[^\(\)\n]+\n'],"tokens"); %get simply assigned variables
% assigned=cellfun(@(x) x(1), assigned);
% assigned=unique(assigned);

[vars,i_vars]=regexp(str,['(' name  ')'],"tokens","start");
[vars,ind]=unique([vars{:}],'stable');
i_vars=i_vars(ind);


[unpack,i_unpack]=regexp(str,['\[' list_lhs '\]' assignment name '\(' list_rhs '\)' ],"tokens","start");

rhs= cellfun(@(x) regexp([x{2}],['(?<!\d)' name],"match"), unpack,'UniformOutput',0);
lhs=cellfun(@(x) regexp(strtrim(regexprep(x{1},[seps '+'],' ')),' ',"split"), unpack,'UniformOutput',0);

[assign,i_assign]=regexp(str,[ '(' name ')' assignment '(:?(:?[' code ']*)?(' name ')(:?[' code ']*)?)*' ],"tokens","start");
% io3=regexp(str,[ '(' name ')' assignment '(:?(:?[' code ']*)?(' name ')(:?[' code ']*)?)*' ],"match");
lhs2=cellfun(@(x) x{1}, assign,'UniformOutput',0);
rhs2=cellfun(@(x) regexp(x{2}, ['(?<!\d)' name] ,"match"), assign,'UniformOutput',0);

i_tot=[i_unpack, i_assign];
rhs=[rhs rhs2];
lhs=[lhs lhs2];

[i_tot,sorted]=sort(i_tot);
rhs = rhs(sorted);
lhs = lhs(sorted);

assigned=cellfun(@(r,l) setdiff(l,r),rhs,lhs,'UniformOutput',0);

if ~isempty(i_tot)
    i_tot=repelem(i_tot,cellfun(@(a) numel(a),assigned));
    [assigned,inds]=unique([assigned{:}],'stable');
    i_tot=i_tot(inds);
    if ~ignore_code_ref
        inds=arrayfun(@(a,i) i_vars(strcmp(a,vars))>=i,assigned,i_tot);
        assigned=assigned(inds);
    end
    
end
errs=regexp(str,[ 'catch[ \t\f]*(' name ')'],'tokens');
assigned = [assigned errs{:}];
mask=cellfun(@(a) numel(a)~=0,assigned);
assigned=assigned(mask);
if ~isempty(assigned)
    assigned=unique(assigned,'stable');
end



end

