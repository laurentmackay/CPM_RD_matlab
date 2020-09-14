function vars = getVars(f)
str=fileread(f);

str=regexprep(str,"%[^\n]*(\n)?","$1"); %remove block comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
numer=['0-9' code(2:end)];
name=['([a-zA-Z][a-zA-Z_0-9]*)'];

assignment='(?:[ \t\f]*\=[ \t\f]*)';
% array_access='(:?\([^\n]*\))?';
% ref=[name array_access];


lhs=[ '(?:[^\n\r' numer '_]+)??(' name ')[' numer ']*\=' ];
lhs=[ '(?:[' numer ']+)??(' name ')[' numer ']*\=' ];
array_inds=[ '\([' numer ']*(' name ')[' numer ']*\)' ];
% rhs=[ '\=[' code ']*(' name ')' ];

rhs=[ '=([^\n\;\r]+)'];
% vars00=regexp(str,rhs,"match")
vars0=regexp(str,rhs,"tokens");
name2=['(?:^|[' code '])' name];
vars0=regexp([vars0{:}], name2 ,"tokens");
vars0=[vars0{:}];
% vars=regexp(str,lhs,"tokens"); %get variables
% vars2=regexp(str,rhs,"tokens"); %get variables
% vars3=regexp(str,array_inds,"tokens"); %get variables
vars=regexp(str,{lhs,array_inds},"tokens");%cellfun(@(x) ,,'UniformOutput',0); %check for the pre-defined variable referencing scenarios
% assigned=vars{1};
vars=feval(@(x,y,z) [x{:} y{:} z{:}],vars{:},vars0); %flatten the cell array of cell arrays
vars0=unique([vars0{:}],'stable');
% vars=[vars{:} vars2{:} vars3{:}];
% vars0=regexp(str,name,"match");
vars=unique(vars,'stable');

vars_bare=regexp(str,['(?:^|[' code '])(' name ')[\*\+\-\/\,]'],"tokens"); %get variables

vars_bare=unique([vars_bare{:}],'stable');




bi=cellfun(@(x) exist(x,'builtin'),vars);
fi= cellfun(@(x) exist(x,'file'),vars);
mask = bi~=5 & (fi~=2 & fi~=6);

vars=vars(mask);


vars=unique([vars vars_bare ],'stable');
% mask = mask & cellfun(@(x) exist(x,'file')~=2 && exist(x,'file')~=6,vars);

% vars=unique([  vars_bare ]);% vars_bare 
end

