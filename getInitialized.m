function assigned = getInitialized(f)
str=fileread(f);


str=regexprep(str,"%[^\n]*(\n)?","$1"); %remove comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition
str=regexprep(str,"\'",""); %remove transposes

code='0-9 \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\:\|\{\}\^\.';
name='[a-zA-Z_$][a-zA-Z_$0-9]*';
assignment='(?:[ \t\f]*\=[ \t\f]*)';
seps='(:?\,| )';
array_access='(:?\([^\n]*\))?';
ref=[name array_access];
list=['((:?[ ]*' ref '[ ]*' seps ')?(:?' ref '))[ ]*'];
% list2=['(:?(:?[ ]*' ref '[ ]*' seps ')?(:?' ref '))[ ]*'];

assigned=regexp(str,['('  name  ') *\=[^\(\)\n]+\n'],"tokens"); %get simply assigned variables
assigned=cellfun(@(x) x(1), assigned);
assigned=unique(assigned);


inout=regexp(str,['\[' list '\]' assignment name '\(' list '\)' ],"tokens"); %get variables assigned by a function with nargout>1

  
% assigned2=cellfun(@(x) split([x{:}],',')', assigned2,'UniformOutput',0);
rhs=cellfun(@(x) regexp([x{2}],['[ \t\f]*' seps '[ \t\f]*'],"split"), inout,'UniformOutput',0);
lhs=cellfun(@(x) regexp([x{1}],['[ \t\f]*' seps  '[ \t\f]*'],"split"), inout,'UniformOutput',0);

io2=regexp(str,[ '(' name ')' assignment '(:?(:?[' code ']*)?(' name ')(:?[' code ']*)?)*' ],"tokens"); 
% io3=regexp(str,[ '(' name ')' assignment '(:?(:?[' code ']*)?(' name ')(:?[' code ']*)?)*' ],"match"); 
lhs2=cellfun(@(x) x{1}, io2,'UniformOutput',0);
rhs2=cellfun(@(x) regexp(x{2}, name ,"match"), io2,'UniformOutput',0);


assigned2=cellfun(@(r,l) setdiff(l,r),[rhs rhs2],[lhs lhs2],'UniformOutput',0);
assigned2=unique([assigned2{:}]);
% assigned2=setdiff([out{:}],[in{:}]);
% assigned=cellfun(@(x) regexp([x{2}],['[ \t\f]*' seps '[ \t\f]*'],"split"), inout,'UniformOutput',0)

if ~isempty(assigned2)
    assigned2=unique(assigned2(:))';
end

% mask = cellfun(@(x) exist(x,'builtin')~=5,assigned);
% mask = mask & cellfun(@(x) exist(x,'file')~=2,assigned);

% % assigned=assigned(mask);
% mask=true(1,length(assigned));
% number=' ?[0-9]+((?:\.|[edfEDF][+\-]?[0-9][0-9\.]*)?';
% %remove assignments that reference the variable
% for i = 1:length(assigned)
%     var=assigned{i};
%     assignments=regexp(str,['[ ^\n\r\]*' var ' *\=([^\n\r\;\=]+)[\n\r\;]'],"tokens");
%     assignments=cellfun(@(x) x(1), assignments);
%     for ass = assignments
%         numerical = regexp(ass,number);
%         if ~isempty(numerical{1}) & numerical{1}==1
%             break; %do not worry about variables that are assigned a numerical value at the beginning of the file
%         end
%         update=regexp(ass,[var '[^\(]'],'match');
%         if ~isempty(update{1})
%             mask(i)=false;
%         end
%     end
% end
% assigned=assigned(mask);
% assigned=[assigned, assigned2];
assigned=assigned2;
end

