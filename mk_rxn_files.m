function  mk_rxn_files(f,save_dir)
if nargin<2
    save_dir='.'
end
simp=false;
predef_spatial = {'bndry_mask'}; %pre-defined spatial variables that will be provided at runtime
elementwise_operations = @(str) regexprep(str,'(?<pre>[^\.])(?<op>\/|\*|\^)','$<pre>.$<op>');  %make multiplication, division, and exponentiation element-wise operations


[chems,S_,rate_constants,fast_chems,fast_pair,fast_affinity, S_cat, species_fast,stoic_fast, S_fast, S_cat_fast,r,p,r_fast,p_fast] = getChemRxns(f);

init=getInitialConditions(f,chems);



vars=getInitialized(f,true);
chems=regexprep(chems,':','');
chems=regexprep(chems,'-','_');

str=fileread(f);

[~,f,~]=fileparts(f);

str=regexprep(str,"%[^\n]*(\n)?",'$1'); %remove block comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,"\'[^\'\n\r]+\'",''); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',''); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',''); %remove function definition




code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
numer=['0-9' code(2:end)];
name=['([a-zA-Z][a-zA-Z_0-9]*)'];





%figure out which quantities are spatially variable
ass=cellfun(@(x) ['(?<![a-zA-Z_0-9])' x '(?:[ \t\f]*)?=([^\n\r\;]+)'],vars,'UniformOutput',0);


par_defs = regexp(str,'(?:\n|^)(?:p|par|param)[ \t\f]([^\n]+)','tokens');
str = regexprep(str,'(?:\n|^)(?:p|par|param)[ \t\f]([^\n]+)','');

model_defs = regexp(str,ass,'match');
model_defs = [model_defs{:}];
if ~isempty(model_defs)
    model_split = regexp(model_defs,[name '[ \t\f]*=([^=]+)'],'tokens');
    model_split=[model_split{:}];
    model_pars = cellfun(@(x) x{2}, model_split, 'UniformOutput', false);
    model_pars = regexp(model_pars,name,'tokens');
    model_pars = [model_pars{:}];
    model_pars = [model_pars{:}];
else
    model_defs = {};
    model_split = {};
    model_pars = {};
end



rate_con_vars = regexp(rate_constants,name,'tokens');
rate_con_vars = [rate_con_vars{:}];
rate_con_vars = [rate_con_vars{:}];

fast_params = regexp(fast_affinity,name,'tokens');
fast_params = [fast_params{:}];
if ~isempty(fast_params)
    fast_params = [fast_params{:}];
end


model_vars = cellfun(@(x) x{1}, model_split, 'UniformOutput', false);
model_var_defs = cellfun(@(x) x{2}, model_split, 'UniformOutput', false);

model_pars = unique([model_pars fast_params rate_constants],'stable');
model_pars = setdiff(model_pars, [model_vars chems],'stable');
model_par_vals = cellstr(repmat('0.0',1,length(model_pars)));









if ~isempty(par_defs)
    par_defs=[par_defs{:}];
    
    par_default = regexp(par_defs,'\*(?:[ \t\f]*)?=([^\,\;]+)', 'tokens' );
    par_default = par_default{1};
    
    if ~isempty(par_default)
        model_par_vals = cellstr(repmat(string(par_default{1}{1}),1,length(model_pars)));
    end
    
    pars= regexp(par_defs,[name '(?:[ \t\f]*)?='], 'tokens' );
    pars=[pars{:}];
    pars=[pars{:}];
    
    par_vals =      regexp(par_defs,[name(2:end-1) '(?:[ \t\f]*)?=([^\,\;\n\r]+)'], 'tokens' );
    par_vals = [par_vals{:}];
    par_vals = [par_vals{:}];
    
    if ~isempty(pars)
        i_par = cellfun(@(p) find(strcmp(p,model_pars),1),pars,'UniformOutput',0);
        matched = cellfun(@(x) ~isempty(x) ,i_par);
        model_par_vals([i_par{:}]) = par_vals(matched);
    end
    vars=setdiff(vars,pars,'stable');
    
    
end

%% re-ordering model_pars to the same order that parameter values are specified in
[model_pars_tot, ia, ib] = union(pars,model_pars,'stable');
model_par_vals_tot = [par_vals(ia) model_par_vals(ib)];
model_pars = intersect(model_pars_tot,model_pars,'stable');
if size(model_pars,1)~=1
    model_pars=model_pars';
end

if ~isempty(pars)
    i_par = cellfun(@(p) find(strcmp(p,model_pars),1),pars,'UniformOutput',0);
    matched = cellfun(@(x) ~isempty(x) ,i_par);
    unknown = setdiff(1:length(model_pars),[i_par{:}]);
    if ~isempty(unknown)
        error(['Please provide values for the following parameters: ' strjoin(model_pars(unknown),', ')])
    end
    model_par_vals([i_par{:}]) = par_vals(matched);
end


if ~isempty(model_pars)
    par_refs = strcat(['(?<=(?:^|[' code ']))('],model_pars_tot,[')(?=(?:$|[' code ']))'])';
    par_reps = string(model_par_vals_tot);
    par_reps = regexprep(par_reps,par_refs,par_reps);
    model_par_vals = regexprep(model_par_vals,par_refs,par_reps);
else
    par_reps={};
    par_refs={};
end
spatial=false(size(vars));
tmp=~spatial;

rhs=regexp(str,ass,'tokens');
rhs=[rhs{:}];
if ~isempty(rhs)
    rhs=cellfun(@(x) x{1},rhs,'UniformOutput',false);
end
spatial_vars=[chems predef_spatial];
ref=@(x) ['(:?[' code '\n]|^)(?<var>' x ')(:?[' code '\n]|$)'];

spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);
if ~isempty(rhs)
    tmp=cell2mat(cellfun(@(x) any(cellfun(@(y) ~isempty(y),regexp(x,spatial_ref,'tokens'))),rhs,'UniformOutput',0));
    
    while any(tmp~=spatial)
        
        spatial_vars=unique([spatial_vars vars(tmp)],'stable');
        spatial(tmp&~spatial)=1;
        
        spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);
        tmp=cell2mat(cellfun(@(x) any(cellfun(@(y) ~isempty(y),regexp(x,spatial_ref,'tokens'))),rhs,'UniformOutput',0));
        
    end
    
end

spatial_vars=spatial_vars(length(chems)+1:end);
ref=@(x) ['(?<pre>[' code '\n]|^)(?<var>' x ')(?<post>[' code '\n]|$)'];
spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);

model_pars=model_pars(~cellfun(@isempty,model_pars));

sym_names = unique([model_pars model_vars chems],'stable');
sym_str= strjoin(strcat(sym_names,"=sym('",sym_names,"','real');"),newline);


pair_inds = arrayfun(@(i) strcmp(chems,fast_pair{i})',1:length(fast_pair),'UniformOutput',false);
fast_inds = arrayfun(@(i) strcmp(chems,fast_chems{i})',1:length(fast_chems),'UniformOutput',false);


S_tot = [S_ S_fast];

lcon = null(S_tot','r');
lcon=rref(lcon')';

N_con=size(lcon,2);

rate_str = @(r) regexprep(strjoin(strcat(chems(r>0),cellstr(repmat('.^',[nnz(r(r>0)),1]))',cellstr(string(r(r>0)))),'.*'),'\.\^1(?=($|\.))','');
rates_fast=cell(length(fast_affinity),1);
for i=1:length(fast_affinity)
    r__=S_fast(:,2*i-1:2*i);
    rates_fast{i} = strcat(rate_str(r__(:,1)'),'*',fast_affinity{i},'-',rate_str(r__(:,2)'));
end


%% user-specified model definitions

lines=[model_defs arrayfun(@(i) ['alpha_chem(vox+' num2str(i-1) '*sz)=(' rate_constants{i} ')' strjoin(arrayfun(@(j) chemRateString(chems{j},-S_(j,i)),1:size(S_,1),'UniformOutput',0),'')],1:size(S_,2),'UniformOutput',0)];

lines = regexprep(lines,spatial_ref,'$<pre>$<var>\(vox\)$<post>');
% fast_affinity_0=fast_affinity;
% fast_affinity = regexprep(fast_affinity,spatial_ref,'$<pre>$<var>\(vox\)$<post>');

chem_ref=cellfun(@(x) ref(x),chems,'UniformOutput',0);
chem_rep=arrayfun(@(i) ['$<pre>x\(vox+' num2str(i-1) '*sz\)$<post>'],1:length(chems),'UniformOutput',0);

lines=regexprep(lines,chem_ref,chem_rep);


lines=cellfun(@(x) [x ';'], lines,'UniformOutput',0);

%% initializing some symbolics stuff



% chem_str2 = strjoin(chems(1:N_slow),',');
if ~isempty(model_defs)
    def_str=strjoin(strcat(model_defs,";"),newline);
else
    def_str='';
end

eval(strcat(sym_str,';',newline, def_str))

% eval([ 'assume([' strjoin([[model_vars']' model_pars chems]) ']>0)'])
% eval([ 'assume([' strjoin([[model_vars']' model_pars chems]) ']>0)'])
% S_fast  = [-[pair_inds{:}]+[fast_inds{:}] [pair_inds{:}]-[fast_inds{:}]]; %only simple reversible complexation is supported right now :(
%this gives you all the slow variables
is_fast=~any(S_',1)&any(S_fast',1);

if any(ismissing(init(~is_fast)))
    error('Please specify some initial condtions for the system')
end

[rates_naive, rxn_naive] = rate_strings(chems,r,p,rate_constants);


%% ENUMERATING CONSTRAINTS
QSSA_defs = regexp(str,['(?:\n|^)QSS(?:A?)\(([ \t\f]*' name(2:end-1) '[ \t\f]*[\,])*' name '(?=\))'],'tokens');
QSSA_defs = [QSSA_defs{:}];
if ~isempty(QSSA_defs)
    QSSA_defs = QSSA_defs(cellfun(@(x) ~isempty(x),QSSA_defs));
    i_QSSA = cellfun(@(x) find(strcmp(x,chems),1),QSSA_defs);
else
    i_QSSA=[];
end
fast_chems = [QSSA_defs chems(is_fast)'];

is_fast(i_QSSA)=true;


% elim =fast_chems;
% elim_0=elim;

fast_eqns = (0 == str2sym([rates_naive{i_QSSA} rates_fast]));
fast_eqns_0 = fast_eqns;
con_str=[];
N_chem = length(chems);
N_fast = nnz(is_fast);
N_fast_eqn = length(fast_eqns);
N_slow = length(chems)-N_fast;
% N_elim_0 = length(elim);

p_elim_con={};
for i=1:N_con
    
    l=lcon(:,i);
    inds=l~=0;
    con_str = [con_str strjoin(strcat(string(l(inds)),'*',chems(inds)'),'+')];
    
    p_elim_con{end+1}=chems(inds);
    
end



consrv_nm = strcat('cnsrv_',string(1:N_con)) ;
eval(['syms ' char(strjoin(consrv_nm,' '))])

eval(['assume([ ' char(strjoin(consrv_nm,', ')) ']>0)'])
consrv_eqns = eval(['[' char(strjoin(strcat(consrv_nm,'==',con_str),';')) ']']);


eqns=sym2cell([consrv_eqns; fast_eqns ]);

%% ELIMINATING VARIABLES
p_elim_fast = regexp(cellstr(string(fast_eqns)),name,'tokens');
p_elim_fast=cellfun(@(x) [x{:}],p_elim_fast,'UniformOutput',0);
p_elim_fast=cellfun(@(x) intersect(x,chems),p_elim_fast,'UniformOutput',0);

elim_con={};
elim_con_eqn={};
% 
% if N_elim_0<N_con %under specified
%     
%     
%     eqns=fast_eqns;
%     del=true;
%     for e=fast_chems'
%         p_elim=[p_elim_fast; p_elim_con'];
%         ind=cellfun(@(pe) strcmp(e,pe),p_elim, 'UniformOutput', 0);
%         contains = cellfun(@any,ind);
%         take = find(contains,1);
%         
%         del=true;
%         if isempty(setdiff(p_elim{take},e, 'stable'))
%             del=false;
%         end
%         p_elim_fast(contains(1:N_fast)) = cellfun(@(pe) setdiff(pe,e, 'stable'), p_elim_fast(contains(1:N_fast)),'UniformOutput',0);
%         p_elim_con(contains(N_fast+1:end)) = cellfun(@(pe) setdiff(pe,e, 'stable'), p_elim_con(contains(N_fast+1:end)),'UniformOutput',0);
%         if take <= N_fast
%             p_elim_fast(take)=[];
%             eqns(take)=[];
%         end
%     end
%     
% %     p_elim = [p_elim_con'; p_elim_fast]; %p_elim_fast should be empty....
% % if numel(p_elim)>1
% 
% % else
% %     p_elim=1;
% % end
% 
% else %over specified
%     elim = elim(1:N_fast);
% end

    elim_fast={};
    p_elim_fast = cellfun(@(x) intersect(x,fast_chems),p_elim_fast,'UniformOutput',0);
    while length(elim_fast)~=length(fast_eqns)
        e=p_elim_fast{1}{1};
        p_elim=[p_elim_fast; p_elim_con'];
        ind=cellfun(@(pe) strcmp(e,pe),p_elim, 'UniformOutput', 0);
        contains = cellfun(@any,ind);
        take = find(contains,1);
        
        del=true;
        if isempty(setdiff(p_elim{take},e, 'stable'))
            del=false;
        end
        tmp = length(p_elim_fast);
        p_elim_fast(contains(1:tmp)) = cellfun(@(pe) setdiff(pe,e, 'stable'), p_elim_fast(contains(1:tmp)),'UniformOutput',0);
        p_elim_con(contains(tmp+1:end)) = cellfun(@(pe) setdiff(pe,e, 'stable'), p_elim_con(contains(tmp+1:end)),'UniformOutput',0);
        
        if take <= N_fast_eqn
            p_elim_fast(take)=[];
            eqns(take)=[];
        end
        
        elim_fast{end+1}=e;
        
    end



    p_elim = cellfun(@(pe) intersect(pe,chems(~is_fast),'stable'),p_elim_con,'UniformOutput',0);
    p_elim = p_elim_con;
    eqns = sym2cell(consrv_eqns);
    
    while length(elim_con)~=N_con
        e=p_elim{1}{1};
        p_elim = cellfun(@(pe) setdiff(pe,e, 'stable'), p_elim,'UniformOutput',0);
        p_elim(1)=[];
        
%         elim{end+1} = e ;
        elim_con_eqn{end+1} = eqns{1};
        elim_con{end+1} = e;
        
        eqns(1)=[];
        
    end
    consrv_eqns_elim = eval(['[' char(strjoin(string(elim_con_eqn),' ')) ']']);






if size(elim_fast,1)==1
    elim_fast=elim_fast';
end

N_elim = length(elim_fast);
N_wp =  N_chem - N_elim;


warning ('off','symbolic:solve:SolutionsDependOnConditions');

if ~isempty(fast_chems)
    
    sol_fast_0 = eval(['solve( fast_eqns,[' strjoin(elim_fast,' ')  '])']);
    if length(fast_chems)~=1
        sol_fast_0 = struct2cell(sol_fast_0);
    else
        sol_fast_0 = {sol_fast_0};
    end
else
    elim_defs={};
    
end

%SOLVE THE CONSERVATION EQUATIONS
sol_consrv=eval(['solve( consrv_eqns_elim,[' strjoin(elim_con,' ')  '])']);
if length(elim_con)==1
%     elim_defs =  strcat(elim_con{1},'=',string(sol_consrv));
    sol_rhs =  {sol_consrv'};
    consrv_subs =  sol_consrv;
else
%     elim_defs = cellfun(@(e) strcat(e,'=',string(sol_consrv.(e))),elim_con)';
    sol_rhs = struct2cell(sol_consrv);
    consrv_subs = cell2sym(struct2cell(sol_consrv));
    sol_consrv = struct2cell(sol_consrv);
end

% fast_eqns

% fast_refs = strcat('(?<=(?:^|[\+\- \*\(]))(',chems(is_fast),')(?=(?:$|[\+\- \*\.\)]))')';

% fast_reps =regexprep(string(sol_rhs{1:N_fast}),slow_refs, slow_init_reps);
% consrv_arr = eval(strcat('[',strjoin(string(sol_consrv)),']'));
% fast_arr = eval(strcat('[',strjoin(string(fast_chems)),']'));
if ~isempty(fast_chems)
%     fast_eqns = subs(fast_eqns, str2sym(elim_con'), consrv_subs);
    sol_fast=eval(['solve( fast_eqns,[' strjoin(elim_fast,' ')  '])']);
    if length(fast_chems)==1
        elim_defs = cellstr(strcat(fast_chems{1},'=',string(sol_fast)));
        sol_rhs = [{sol_fast}; sol_rhs];
        sol_fast={sol_fast};
    else
        elim_defs = [cellfun(@(e) strcat(e,'=',string(sol_fast.(e))),elim_fast)];
        sol_rhs =  [struct2cell(sol_fast); sol_rhs];
        sol_fast = struct2cell(sol_fast);
    end
else
    sol_fast={};
    %     elim_defs={};
    %      sol_rhs = {};
end


f_chems = cellfun(@(c) ['f_' c], chems, 'UniformOutput', false);


eval([ ' syms ' strjoin(f_chems,' ') ';']);
eval([ 'assume([' strjoin([[model_vars']' model_pars]) ']>0)']);
eval(strcat(' assume([', strjoin(f_chems), "],'real')"));
eval(strcat(' assume([', strjoin(chems), "]>0)"));



ec0=subs(consrv_eqns_elim,str2sym(elim_fast),cell2sym(sol_fast));
sol_consrv=eval(['solve( ec0,[' strjoin(elim_con,' ')  '])']);
if length(elim_con)==1
%     elim_defs =  strcat(elim_con{1},'=',string(sol_consrv));
   sol_rhs(N_fast_eqn+1:end) =  {sol_consrv'};
    consrv_subs =  sol_consrv;
else
%     elim_defs = cellfun(@(e) strcat(e,'=',string(sol_consrv.(e))),elim_con)';
    sol_rhs(N_fast_eqn+1:end) = struct2cell(sol_consrv);
    consrv_subs = cell2sym(struct2cell(sol_consrv));
    sol_consrv = struct2cell(sol_consrv);
end




slow_refs = strcat('(?<=(?:^|[',code,']))(',chems(~is_fast),')(?=(?:$|[',code,']))')';
fast_refs = strcat('(?<=(?:^|[',code,']))(',chems(is_fast),')(?=(?:$|[',code,']))')';
consrv_rhs = regexprep(string(consrv_eqns),'^.*==','');


slow_init_reps = cellstr(string(init(~is_fast)));
con_fast = cellfun(@(x) any(strcmp(x,fast_chems)),elim_con);
if N_fast>0
    
    sol_fast_1 = sym2cell(subs(cell2sym(sol_fast_0), str2sym(elim_con'), cell2sym(sol_consrv)));
    fast_init_reps = regexprep(cellstr(string([sol_fast_1; sol_consrv(con_fast)])),slow_refs, slow_init_reps);
    fast_init_reps = cellstr(string(str2sym(fast_init_reps)));
    fast_init_reps = regexprep(fast_init_reps, par_refs, par_reps);
    init(is_fast) = fast_init_reps;
else
    fast_init_reps = {};
end


consrv_defs_init = regexprep(consrv_rhs,[slow_refs; fast_refs],[slow_init_reps'; fast_init_reps]');
consrv_defs_init = cellstr(string(str2sym(consrv_defs_init)));
consrv_def_vals = consrv_defs_init;

consrv_defs = regexprep(string(consrv_eqns),'==','=');

par_eqn = str2sym(model_pars');
par_val_eqn = str2sym(model_par_vals);
consrv_nm_eqn = str2sym(consrv_nm');
consrv_val_eqn = str2sym(consrv_def_vals);
chem_eqn = str2sym(chems(~is_fast));
init_eqn = str2sym(init(~is_fast));

for i=1:length(sol_consrv)
if length(sol_consrv{i})>1
    sol = sol_consrv{i};
    sols_sub = subs(sol,[par_eqn; consrv_nm_eqn; chem_eqn'],[par_val_eqn; consrv_val_eqn; init_eqn']);
    positive = isAlways(sols_sub>0);
    
    if ~any(positive)
        error(['No positive solution for the conservation of matter was found:' ec0{i}]);
    else
        sol_consrv{i}=sol(find(positive,1));
    end
    
end
end



warning ('on','symbolic:solve:SolutionsDependOnConditions');

consrv_arr = eval(strcat('[',strjoin(string(sol_consrv)),']'));
fast_chem_arr = eval(strcat('[',strjoin(string(elim_fast)),']'));
fast_arr = eval(strcat('[',strjoin(string(sol_fast)),']'));
if iscell(sol_consrv)
    sol_consrv = cell2sym(sol_consrv);
end
tt=subs(sol_consrv,fast_chem_arr,fast_arr);
sol_rhs(N_fast_eqn+1:end) =  sym2cell(tt);

cellfun(@(e) disp(['Eliminating: ' e]),elim_defs)
elim_str = strjoin(elim_defs,newline);
if ~isempty(S_fast)
    is_slow = ~any(S_fast');
else
    is_slow = true(size(chems));
    is_slow(i_QSSA)=false;
end
is_elim = cellfun(@(c) any(strcmp(c,elim_fast)),chems);



% f_slow = cellfun(@(c) ['f_' c], {chems{~is_elim}}, 'UniformOutput', false);

f_chems = cellfun(@(c) ['f_' c], chems, 'UniformOutput', false);


% eval([ ' syms ' strjoin(f_chems,' ') ';']);
% eval([ 'assume([' strjoin([[model_vars']' model_pars]) ']>0)']);
% eval(strcat(' assume([', strjoin(f_chems), "],'real')"));
% eval(strcat(' assume([', strjoin(chems), "],'real')"));

%% Computing GAMMA

% Gamma = eval([ ' [' strjoin( cellfun(@(aff, pair ) ['simplify(' aff '*' pair ',"Steps",10)'],fast_affinity, fast_pair,'UniformOutput',0),';') ']']   );
if length(fast_eqns)>=1
    Gamma=eval(strcat("[",strjoin(string([sol_fast_0; sol_consrv(con_fast)])),"]"))';
else
    Gamma=[];
end
Gamma_0=Gamma;
Gamma_simp = Gamma;
i=0;
% nm = ['sigma0__' int2str(i)];
% [Gamma_simp,sigma__rep] = subexpr(Gamma_simp,nm);

%% COMPUTING J_GAMMA
J_gamma=arrayfun(@(gi) cellfun(@(j) simplify(diff(gi,j),'Steps',10),{chems{~is_fast}},'UniformOutput',0),Gamma_0,'UniformOutput',false)';
if ~isempty(J_gamma)
    J_gamma = [J_gamma{:}];
    J_gamma = [J_gamma{:}];
    J_gamma=reshape(J_gamma,[N_slow,size(Gamma,1)])';
end
J_gamma_0 = J_gamma;
%
%


big = [Gamma_0 J_gamma];
i=0;
nm = ['subs__' int2str(i)];

sigma__reps={};
if ~isempty(big)
    [big_simp,sigma__rep] = subexpr(big,nm);
    if ~isempty(sigma__rep)
        
        while ~isempty(sigma__rep)
            
            
            sigma__reps{end+1}=[nm  ' = ' char(simplify(sigma__rep,'Steps',10))];
            
            i=i+1;
            nm = ['subs__' int2str(i)];
            [big_simp, sigma__rep] = subexpr(big_simp,nm);
            
            
        end
        
    end
    
    big_simp = simplify(big_simp);
    if simp
        J_gamma = big_simp(:,2:end);
        Gamma = big_simp(:,1);
    end
end



% fast_pair_sym = eval(['[' strjoin(fast_pair) ']'])';
% aff_simp =  Gamma ./ fast_pair_sym;



%% setting up the liner system to determine reactions rates

lcon_fast_0 = lcon(is_fast,:);

i_con_fast = any(lcon_fast_0 > 0,1);

N_con_fast = nnz(i_con_fast);

% m0 = sym('m0', [N_con_fast N_slow]);

j_gam = sym('J_gamma_', [N_fast N_slow]);
if ~isempty(J_gamma)
    %     m0( J_gamma == 0 )=0;
    
    j_gam( J_gamma == 0) =0;
end


m11=j_gam;

lcon_fast=lcon_fast_0(:,i_con_fast);
lmix = lcon(~is_fast,i_con_fast)';
lslow_0 = lcon(~is_fast,~i_con_fast)';
% for i=1:N_con_fast
%     con_fast=lcon_fast(:,i)~=0;
%     if any(con_fast)
%         m0(i,:)=lmix(i,:)+sum(J_gamma(con_fast,:),1);
%     end
% end

l_sf=sym(lmix(1:N_con_fast,:));
for i=1:N_con_fast
    l_sf(i,:)=l_sf(i,:)+sum((lcon_fast(:,i)'.*j_gam')',1);
end
% l_sf=(lcon_fast'.*j_gam)+lmix(1:N_con_fast,:);
% b__0 = eval(['[' strjoin( f_slow,'; ') ']'] );
% for i=1:size(lslow,1)
% % irep =  find(any(m10(mask,:).*lslow(i,:),2),1);
% % mask
% % m10(irep,:) = lslow(i,:)
% end

% is_slow = ~is_elim;
is_solo = any(lmix( sum(abs(lmix),2)==1,:),1);
i_non_fast = find(~is_fast);
ind = is_slow&~is_fast;

eye_slow_0= diag(ind,0);
eye_slow_0 = eye_slow_0(ind,~is_fast);
% if ~isemtpty(i_non_fast(is_solo))
eye_slow_0 = eye_slow_0(~any(i_non_fast(is_solo) == find(ind)',2),:);
%finding a way to express our system that isnt undetermined
A_mat=[lslow_0; eye_slow_0]';
[A_rref,p_A]=rref(A_mat);
A_li_rows = A_mat(:,p_A)';
if size(A_li_rows,1)+N_con_fast>N_slow
    A_li_rows=A_li_rows(1:N_slow-N_con_fast,:);
elseif size(A_li_rows,1)+N_con_fast<N_slow
    error('Underdetermined slow sub-system, try removing some fast reactions.')
end

mask_slow_0 = sum(A_li_rows~=0,2)==1;

% i_non_fast = find(~is_fast);
mask_non_elim = any(i_non_fast==find(~is_elim)',1);
% i_slow_0 = i_non_fast(mask_slow_0);

i_0 = mod(find(A_li_rows(mask_slow_0,:)')-1,N_slow)'+1;


ij_slow = mod(find(A_li_rows(mask_slow_0,:)')-1,N_slow)+1;
f__0 = eval([ '[' strjoin(strcat('f_', chems )) ']']);
b__0 =  eval([ '[' strjoin(strcat('f_', chems(~is_fast) )) ']']);

b__0(:) =  0;
inds= i_non_fast(ij_slow);
b__0(mask_slow_0) = f__0(inds);

project=true;
if project

f_mix = [A_li_rows;l_sf]\[b__0]';
% f_mix2 = eval(['[m10;m0]\[' strjoin( f_slow,'; ') ';' strjoin(cellstr(string(zeros(1,size(J_gamma,1)))),'; ')  ']']);

else
    f_mix = f__0(~is_fast)';

end

%
f_fast =  simplify(j_gam*f_mix,'Steps',20);
m100 = A_li_rows(mask_slow_0,:);
f_mix_simp = simplify([f_mix; f_fast],'Steps',20);
i=0;
nm = ['subs2__' int2str(i)];
[f_mix_simp,sigma2__rep] = subexpr(f_mix_simp,nm);
sigma2__reps={};
if ~isempty(sigma2__rep)
    
    while ~isempty(sigma2__rep)
        sigma2__reps{end+1}=[nm ' = ' char(simplify(sigma2__rep)) ];
        i=i+1;
        nm = ['subs2__' int2str(i)];
        [f_mix_simp, sigma2__rep] = subexpr(f_mix_simp,nm);
        
    end
    
end
% sigma2__reps=sigma2__reps';
% if simp
    f_tot0=[f_mix; f_fast];
    f_tot=f_mix_simp;
% end
% f_tot = sym('f_tot',[N_chem,1]);
% f_tot(~is_fast)=f_mix_simp;
% f_tot(is_fast)=f_fast;

% f_fast_simp = simplify(f_fast,'Steps',20);
% i=0;
% nm = ['sub3__' int2str(i)];
% [f_fast_simp,sigma3__rep] = subexpr(f_fast_simp,nm);
% sigma3__reps={};
% if ~isempty(sigma3__rep)
%     
%     while ~isempty(sigma3__rep)
%         sigma3__reps{end+1}=[nm ' = ' char(simplify(sigma3__rep)) ];
%         i=i+1;
%         nm = ['subs3__' int2str(i)];
%         [f_fast_simp, sigma3__rep] = subexpr(f_fast_simp,nm);
%         
%     end
%     
% end
% % sigma2__reps=sigma2__reps';
% f_fast=f_fast_simp;

warning ('off','symbolic:sym:isAlways:TruthUnknown');
non_zero = ~arrayfun(@(m) isAlways(j_gam(m)==0) || isAlways(j_gam(m)==1),1:numel(j_gam));
warning ('on','symbolic:sym:isAlways:TruthUnknown');

m1_gamma=j_gam(non_zero);
m0=J_gamma(non_zero);

predefs = arrayfun(@(i) [char(m1_gamma(i)) '=' char(m0(i))],1:numel(m1_gamma),'UniformOutput',false);

% predefs2 = arrayfun(@(i) ['affinity__' int2str(i) '=' char(aff_simp(i))],1:length(aff_simp),'UniformOutput',false);

rates=cell(size(S_,2),1);
for i=1:size(S_,2)
%     r=-S_(:,i);
%     r(r<0)=0;
    inds = r(:,i)>0;
    terms = strcat(chems(inds),cellstr(repmat('^',[nnz(inds),1]))',cellstr(string(r(r(:,i)>0,i)))');
    terms = regexprep(terms,'\^1$',''); %remove exponent 1's from rate computation
    rate_str = [ rate_constants{i} '*' strjoin(terms,'*')];
%     if any(S_cat(:,i))
%         rate_str = [ rate_str '.*' strjoin(chems(S_cat(:,i)),'.*') ];
%     end
    rates{i} = rate_str;
end
elim_refs = nameref(elim_fast);
elim_reps = strcat('(',string(sol_fast),')');

%  fast_refs = strcat('(?<=(?:^|[\+\- \*]))(',fast_chems,')(?=(?:$|[\+\- \*.]))')';
%  elim_reps = regexprep(elim_reps,fast_refs,string(sol_fast))
rates_naive = rates;

if ~isempty(elim_fast)
    rates_elim = regexprep(rates,elim_refs,elim_reps);
    rates = rates_elim;
end












ind_slow = find(is_slow);
N_pair = floor(size(S_,2)/2);

paired = nan(size(S_,2),1);
for i=1:size(S_,2)
    if isnan(paired(i))
        s_vec=S_(:,i);
        pairing = find(all(s_vec==-S_),1);
        if ~isempty(pairing) && isnan(paired(pairing))
            paired(i)=pairing;
            paired(pairing)=i;
        end
    end
end
c=0;


consrv_refs = strcat('(?<=(?:^|[\+\- \*\(]))(',elim_con,')(?=(?:$|[\+\- \*\.\)\^]))')';
consrv_reps = cellfun(@(x) strcat("(",x,")"),  string(sol_consrv));
rates_consrv = regexprep(rates_naive,consrv_refs,consrv_reps);
 
fid=fopen(strcat(save_dir,filesep, f,'.rxns'),'w');

for i=1:length(paired)
    if paired(i)~=-1
        if ~isnan(paired(i))
            sl=['R' int2str(c) '=' rates_naive{i} ' - ' rates_naive{paired(i)} ';' newline];
            paired(paired(i))=-1;
        else
            sl=['R' int2str(c) '=' rates_naive{i} ';' newline];
        end
        sl=regexprep(sl,'\.\*','*');
        sl=regexprep(sl,'\_','');
        fwrite(fid,sl,'char');
        c=c+1;
    else
        skip=false;
    end
    
end

fwrite(fid,[newline newline],'char');
fwrite(fid,strcat(elim_str,newline,strjoin(consrv_defs,newline)),'char');

fwrite(fid,[newline newline],'char');
c=0;
rates_paired=cell(size(S_,2),1);
for i=1:length(paired)
    if paired(i)~=-1
        if ~isnan(paired(i))
            sl=['R' int2str(c) '=' rates_consrv{i} ' - ' rates_consrv{paired(i)} ';' newline];
            rates_paired{paired(i)}='';
            paired(paired(i))=-1;
        else
            sl=['R' int2str(c) '=' rates_consrv{i} ';' newline];
        end
        rates_paired{i}=['R' int2str(c)];
        sl=regexprep(sl,'\.\*','*');
        sl=regexprep(sl,'\_','');
        fwrite(fid,sl,'char');
        c=c+1;
    else
        skip=false;
    end
    
end

fwrite(fid,[newline newline],'char');

% f_slow={};
for chem_ind = find(~is_fast)
    %     chem_ind = find(m10(i,:),1);
    S_chem = S_(chem_ind,:);
    inds = find(S_chem);
    inds=inds(cellfun(@(x) ~isempty(x),rates_paired(inds)));
%     pairs = paired(inds);
%     pairs = pairs(pairs~=-1);
    
    terms = strcat(cellstr(num2str(S_chem(inds)')),'*',rates_paired(inds),'');
    terms = regexprep(terms,'(?<![0-9\.\-\+A-Za-z_])1\*',''); %remove multiplication by 1's from rate computation
    terms = regexprep(terms,'(?<=^)-1\*','-'); %remove multiplication by -1's from rate computation
    sum_terms = strjoin(terms,'+'); % sum terms
    sum_terms = regexprep(sum_terms,'+[ \+]*+','+');
    sum_terms = regexprep(sum_terms,'+-','-');% remove +-, replace with -
    fwrite(fid, char(strcat("d",chems(chem_ind),"dt=",sum_terms, string(newline))) ,'char');
%     f_slow{end+1} = sum_terms;
    
end

fwrite(fid,[newline newline],'char');
fwrite(fid,strcat('$Assumptions = {',strjoin(strcat([model_pars, chems],'>0'),', '),'};'),'char');
fwrite(fid,[newline newline],'char');

fclose(fid);


rates_slow={};
f_slow={};
for chem_ind = i_non_fast(ij_slow)
    %     chem_ind = find(m10(i,:),1);
    S_chem = S_(chem_ind,:);
    inds = find(S_chem);
    terms = strcat(cellstr(num2str(S_chem(inds)')),repmat('*(',[length(inds),1]),rates(inds),repmat(')',[length(inds),1]));
     terms = strcat(cellstr(num2str(S_chem(inds)')),'*(',rates(inds),')');

    terms = regexprep(terms,'(?<![0-9\.\-\+A-Za-z_])1\*',''); %remove multiplication by 1's from rate computation
    terms = regexprep(terms,'(?<=^)-1\*','-'); %remove multiplication by -1's from rate computation
    sum_terms = strjoin(terms,'+'); % sum terms
    sum_terms = regexprep(sum_terms,'+-','-');% remove +-, replace with -
    if length(sum_terms)==0
        sum_terms='0';
    end
    rates_slow{end+1} = ['f_' chems{chem_ind} ' = ' sum_terms];
    f_slow{end+1} = sum_terms;
    
end



rates_chem={};
% f_slow={};
for chem_ind = find(~is_fast)
    %     chem_ind = find(m10(i,:),1);
    S_chem = S_(chem_ind,:);
    inds = find(S_chem);
    terms = strcat(cellstr(num2str(S_chem(inds)')),repmat('*(',[length(inds),1]),rates(inds),repmat(')',[length(inds),1]));
    terms = regexprep(terms,'(?<![0-9\.\-\+A-Za-z_])1\*',''); %remove multiplication by 1's from rate computation
    terms = regexprep(terms,'(?<=^)-1\*','-'); %remove multiplication by -1's from rate computation
    sum_terms = strjoin(terms,'+'); % sum terms
    sum_terms = regexprep(sum_terms,'+-','-');% remove +-, replace with -
    if length(sum_terms)==0
        sum_terms='0';
    end
    rates_chem{end+1} = ['f_' chems{chem_ind} ' = ' sum_terms];
%     f_slow{end+1} = sum_terms;
    
end

f_slow = simplify(str2sym(f_slow));

f_non_elim={};
for chem_ind = find(~is_elim)
    %     chem_ind = find(m10(i,:),1);
    S_chem = S_(chem_ind,:);
    inds = find(S_chem);
    terms = strcat(cellstr(num2str(S_chem(inds)')),repmat('*(',[length(inds),1]),rates(inds),repmat(')',[length(inds),1]));
    terms = regexprep(terms,'(?<![0-9\.\-\+A-Za-z_])1\*',''); %remove multiplication by 1's from rate computation
    terms = regexprep(terms,'(?<=^)-1\*','-'); %remove multiplication by -1's from rate computation
    sum_terms = strjoin(terms,'+'); % sum terms
    sum_terms = regexprep(sum_terms,'+-','-');% remove +-, replace with -
    if length(sum_terms)==0
        sum_terms='0';
    end
%     rates_slow{end+1} = ['f_' chems{chem_ind} ' = ' sum_terms];
    f_non_elim{end+1} = sum_terms;
    
end


f_non_elim = simplify(str2sym(f_non_elim));

disp('Reaction Rates (naive):')
disp(strjoin(strcat('d',chems(~is_elim),'/dt=',string(f_non_elim)),newline))

chem_rep_FVM=arrayfun(@(i) ['$<pre>u\(:,' num2str(i) '\)$<post>'],1:length(chems),'UniformOutput',0);
spatial_rep_FVM=repmat({'$<pre>$<var>\(:\)$<post>'},size(spatial_vars));

model_defs = regexprep(model_defs,elim_refs,elim_reps);

rxn_body = [strjoin([model_defs'; rates_slow'; sigma__reps'; predefs'; sigma2__reps';],[';' newline])...
    ';' newline newline ['Rx = [' char(strjoin(string(f_tot),[',...' newline])) '];']];
rxn_body = regexprep(rxn_body,chem_ref ,chem_rep_FVM);% reshape the chemical names

rxn_body = elementwise_operations(rxn_body); %make multiplication, division, and exponentiation element-wise operations


% fast_affinity_FVM = regexprep(fast_affinity_0,chem_ref ,chem_rep_FVM);
% fast_affinity_FVM=regexprep(fast_affinity_FVM,spatial_ref,'$<pre>$<var>\(:\)$<post>')
fid=fopen(strcat(save_dir,filesep,'eval_Rx.m'),'w');
fwrite(fid,rxn_body,'char');
fclose(fid);




if ~isempty(model_defs)
    model_body = [strjoin(model_defs',[';' newline]) ';'];
else
    model_body='';
end
chem_rep_full=arrayfun(@(i) ['$<pre>x\(:,:,' num2str(i) '\)$<post>'],1:length(chems),'UniformOutput',0);%this only works for 2d
model_body = regexprep(model_body,chem_ref ,chem_rep_full);% reshape the chemical names
model_body = elementwise_operations(model_body);
fid=fopen(strcat(save_dir,filesep,'eval_model.m'),'w');
fwrite(fid,model_body,'char');
fclose(fid);


% init = sym('init', )








consrv_assgn = cellfun(@(nm,v) [nm '=' v],cellstr(consrv_nm), cellstr(consrv_defs_init'),'UniformOutput',false);

consrv_defs_init = strcat("par ",consrv_nm, " = " , consrv_defs_init', string(newline));

%  if simp
% ode_body = [strjoin([model_defs'; rates_slow'; sigma__reps'; predefs'; sigma2__reps';], newline)...
%     ';' newline newline char(strjoin(strcat('d',chems(~is_fast),'/dt=',string(f_mix)'),newline))];
%  else
%  end

% var_refs = nameref(model_vars)

% ode_body_elim = {model_defs', cellstr(elim_defs), rates_slow', sigma__reps', predefs', sigma2__reps'};


f_mix_elim = subs(f_mix,m1_gamma,m0);
 f_tot_explicit = subs(f_tot0,m1_gamma,m0);
% rates_new = eval(['[' strjoin(rates_elim) ']']);

f_mix_elim = simplify(subs(f_mix_elim,f__0(i_non_fast(ij_slow)),f_slow),'steps',10);
f_tot_explicit = subs(f_tot_explicit,f__0(i_non_fast(ij_slow)),f_slow);
f_mix_explicit = f_mix_elim ;
% f_tot_explicit = subs(f_tot_explicit,str2sym(elim),cell2sym(sol_fast));
f_tot_explicit = subs(f_tot_explicit,str2sym(model_vars),subs(str2sym(model_var_defs)));
f_tot_explicit = subs(f_tot_explicit,str2sym(elim_fast),subs(cell2sym(sol_fast),str2sym(model_vars),subs(str2sym(model_var_defs))));
anon_str = strjoin(regexprep(string(f_tot_explicit),...
                   nameref(chems),...
                   arrayfun(@(i) ['u(' int2str(i) ')'],1:length(chems),'UniformOutput',false)),', ');


fid=fopen(strcat(save_dir,filesep,'model_anon.m'),'w');
fwrite(fid,['rhs_anon = @(u) [' char(anon_str) '];'],'char');
fclose(fid);

% f_mix_elim = simplify(subs(f_mix_elim),'steps',10);
% elim_con_elim = subs(sol_consrv,str2sym(elim),cell2sym(sol_fast))
% f_mix_elim = simplify(subs(f_mix_elim,str2sym([elim; elim_con']),cell2sym(sol_rhs)),'steps',10);
% ode_body_elim = cellfun(@(x) strjoin(x,newline),ode_body_elim(~cellfun(@isempty,ode_body_elim)), 'UniformOutput',0);
arrayfun(@(e) disp(strcat("Eliminating: ",e)),strcat(elim_con'," = ",string(sol_consrv)));

% cellfun(@(e) disp(['Eliminating: ' e]),elim_defs)

f_mix_ode = subs(f_mix_explicit, str2sym(elim_con'), sol_consrv);
% f_mix_elim_ode = subs(f_mix_elim, str2sym(elim_con'), sol_consrv);
is_elim_con0 = cellfun(@(c) any(strcmp(c,elim_con)),chems);
is_elim_con = is_elim_con0(~is_fast);

ode_body_elim = char(strjoin(strcat('d',chems(~is_elim_con0 & ~is_fast),'/dt=',string(f_mix(~is_elim_con))'),newline));

ode_body_elim=regexprep(ode_body_elim,'\.\*','\*');

% ode_body_elim = regexprep(ode_body_elim,elim_refs,elim_reps);

consrv_str = strjoin(strcat(strcat('#', consrv_defs), string(newline), consrv_defs_init'), newline);
consrv_str = char(consrv_str);


init_str = ['init ' char(strjoin(strcat(chems(~is_elim_con0 & ~is_fast)', "=", string(init(~is_elim_con0 & ~is_fast)')), ', '))];

rates_0 = char(strjoin( strcat(string(f__0(i_non_fast(ij_slow))),'=',string(f_slow)), newline));

if ~isempty(par_defs)
    model_par_vals = compose('%.9f',eval(['[' strjoin(model_par_vals) ']']));
    par_str0 = [char(strjoin(strcat(model_pars,"=",model_par_vals),[';' newline])) newline];
    par_str = ['par ' char(strjoin(strcat(model_pars,"=",model_par_vals),', ')) newline];
else
    par_str=newline;
end

par_str =  breaklines(par_str,255,['\']);
init_str =  breaklines(init_str,255,['\']);

model_str = strjoin(model_defs,newline);


  elim_con_defs = regexprep(string(sol_consrv),nameref(consrv_nm),consrv_def_vals);
%   elim_con_defs = string(eval(['[' char(strjoin(elim_con_defs)) ']'])');
  
 elim_con_defs = char(strjoin(strcat(elim_con','=',elim_con_defs), newline));

 
%  elim_con_eqn = [cellfun(@(e) strcat(e,'=',string(sol_fast.(e))),fast_chems)];
ode_str=[ode_body_elim newline newline ...
    init_str newline newline ...
    elim_con_defs newline newline ...
    model_str newline newline ...
    rates_0 newline newline ...
    par_str newline newline ];



nms0 = regexp(ode_str,['(?<=(?:d))' name '(?=/dt)'], 'tokens')';
nms0=[nms0{:}];


nms = regexp(ode_str,['(?<=(?:^|[' code ']|\n|d{1}))' name '(?=(?:$|[' code ']|\n|\r))'], 'tokens')';
nms=unique([nms{:}],'stable');

is_ddt = any(cell2mat(cellfun(@(nm) strcmp(['d' nm],nms),nms0,'UniformOutput',false)'));

nms(is_ddt)=[];

too_long=cellfun(@length,nms)>9;

if nnz(too_long)>0
    
    
    too_long_refs = strcat(['(?<=(?:^|[' code ']|\n|d))('],nms(too_long),[')(?=(?:$|[' code ']|\n|\r))'])';
    too_long_reps = arrayfun(@(i) string(['r_' int2str(i) ]),1:nnz(too_long));
    
    ode_str = regexprep(ode_str, too_long_refs, too_long_reps);
    
    ode_str = [ode_str '#THE FOLLOWING REPLACEMENTS HAVE BEEN MADE IN THIS FILE' newline ...
        char(strjoin(strcat("# ",nms(too_long)," -> ",too_long_reps), newline)) ];
    
end
fid=fopen(strcat(save_dir,filesep,f, '.ode'),'w');
fwrite(fid,ode_str,'char');
fclose(fid);













chem_ref_elim=cellfun(@(x) ref(x),chems(~is_elim_con0 & ~is_fast),'UniformOutput',0);
% chem_rep_elim=arrayfun(@(i) ['$<pre>x\(vox+' num2str(i-1) '*sz\)$<post>'],1:length(chems(~is_elim_con)),'UniformOutput',0);


%  fast_eqns = subs(fast_eqns, str2sym(elim_con'), consrv_subs);

chem_rep_FVM_elim=arrayfun(@(i) ['$<pre>u\(:,' num2str(i) '\)$<post>'],1:length(chems(~is_elim_con0 & ~is_fast)),'UniformOutput',0);

rxn_body_elim = [strjoin([model_defs'; rates_chem'; sigma__reps'; predefs'; sigma2__reps';],[';' newline])...
    ';' newline newline ['Rx = [' char(strjoin(string(f_mix_ode(~is_elim_con)),[',...' newline])) '];']];
rxn_body_elim = regexprep(rxn_body_elim,elim_refs,elim_reps);
rxn_body_elim = regexprep(rxn_body_elim,consrv_refs,consrv_reps);
rxn_body_elim = regexprep(rxn_body_elim,chem_ref_elim,chem_rep_FVM_elim);% reshape the chemical names

rxn_body_elim = regexprep(rxn_body_elim,'(?<pre>[^\.])(?<op>\/|\*|\^)','$<pre>.$<op>'); %make multiplication, division, and exponentiation element-wise operations



% fast_affinity_FVM = regexprep(fast_affinity_0,chem_ref ,chem_rep_FVM);
% fast_affinity_FVM=regexprep(fast_affinity_FVM,spatial_ref,'$<pre>$<var>\(:\)$<post>')
fid=fopen(strcat(save_dir,filesep,'eval_Rx_elim.m'),'w');
fwrite(fid,rxn_body_elim,'char');
fclose(fid);


% ode_body=regexprep(ode_body,'\.\*','\*');










pars_tot=[model_pars pars];
[tmp, i ] = unique([model_pars pars],'stable');
tmp2 = [model_par_vals par_vals];
par_vals_tot = tmp2(i);
pars_tot = pars_tot(i);
if ~isempty(pars_tot)
    fid=fopen(strcat(save_dir,filesep,'model_params.m'),'w');
    fwrite(fid, strjoin(strcat([pars_tot cellstr(consrv_nm)],'=',[par_vals_tot consrv_def_vals'],';'),newline),'char');
    fclose(fid);
end


% nms_elim = regexp(string(f_mix_elim(mask_non_elim)),name, 'tokens')';
% nms_elim=[nms_elim{:}];

params = inline_script('model_params');
eval_rhs_str = inline_script('eval_Rx_elim');

rhs_str =['function Rx = rhs_fun(t,u)' newline 'u=transpose(u);' newline params newline eval_rhs_str newline 'Rx=transpose(Rx);' newline 'end'];

fid=fopen(strcat(save_dir,filesep,'rhs_fun.m'),'w');
fwrite(fid,rhs_str,'char');
fclose(fid);

init_vals = eval(strcat("[",strjoin(regexprep(init(~ismissing(init)),par_refs,par_reps)),"]"));


fid=fopen(strcat(save_dir,filesep,'model_ic.m'),'w');
fwrite(fid,strcat("ic = [", strjoin(string(init_vals)), "];"),'char');
fclose(fid);

is_ode = ~is_elim_con0 & ~is_fast;

% ic_ode=init_vals(~is_fast);
ic_ode=init_vals(is_ode);
tol=1e-14;
fprintf(['Integrating to find fixed point (abstol = ' num2str(tol) ')...']);
t0=tic();
T_vec=0;
while ~all(abs(rhs_fun(0,ic_ode'))<tol)
    [T_vec,Y_vec] = ode15s(@ rhs_fun,T_vec(end)+[0 1e3],ic_ode,odeset('NonNegative',1:length(ic_ode)));
    ic_ode = Y_vec(end,:);
end


fp=init;
fp(~is_elim_con0 & ~is_fast)=num2str(ic_ode',12);
fp(is_elim_con0) = subs(str2sym(elim_con)',str2sym(elim_con)' ,subs(subs(sol_consrv),consrv_nm_eqn,consrv_val_eqn));
fp(is_elim_con0) = regexprep(fp(is_elim_con0), par_refs, par_reps);
fp(is_elim_con0) = regexprep(fp(is_elim_con0),  nameref(chems(is_ode)), string(ic_ode));
fp(is_fast) = regexprep(string(subs(cell2sym(sol_fast))),nameref(chems(~is_fast)),fp(~is_fast));
fp(is_fast) = regexprep(fp(is_fast), par_refs, par_reps);
fp = eval(strcat('[',strjoin(fp),']'));
t1=toc(t0);
fprintf(['Done (' num2str(t1) ' seconds elapsed, final time = ' num2str(T_vec(end)) ')' newline]);

fid=fopen(strcat(save_dir,filesep,'model_fp.m'),'w');
fwrite(fid,['fp = [' num2str(fp,12) '];'],'char');
fclose(fid);

init_relaxed = ic_ode;
disp(strjoin(strcat(chems(is_ode),'=',string(init_relaxed)),', '));
disp(strjoin(strcat(chems(~is_ode),'=',string(fp(~is_ode))),', '));

u_read = cellfun(@(c,i) [c '=U(' int2str(i) ')'],chems(is_ode),num2cell(1:nnz(is_ode)),'UniformOutput',false);

if length(model_pars) > 10
    valid_ind = true(1,min(length(model_pars)+1,36));
    valid_ind(11) = false;
else
    valid_ind = true(1,length(model_pars));
end
% valid_ind(24) = false;
par_read = cellfun(@(p,i) [p '=PAR(' int2str(i) ')'],model_pars(1:nnz(valid_ind)),num2cell(find(valid_ind)),'UniformOutput',false);
par_read = [par_read ...
    cellfun(@(p,v) [p '=' v],model_pars(nnz(valid_ind)+1:end),model_par_vals(nnz(valid_ind)+1:end),'UniformOutput',false)];

par_assgn = cellfun(@(v,i,j) ['PAR(' int2str(i) ')=' v ' !' model_pars{j} ],model_par_vals(1:nnz(valid_ind)),num2cell(find(valid_ind)), num2cell(1:nnz(valid_ind)),'UniformOutput',false);
u_assgn = cellfun(@(v,i) ['U(' int2str(i) ')=' v],cellstr(num2str(init_relaxed',12))',num2cell(1:nnz(~is_elim_con)),'UniformOutput',false);

F_assgn = arrayfun(@(f,i)[  'F(' int2str(i) ')=' char(string(f))],f_mix_ode(~is_elim_con),(1:nnz(~is_elim_con))','UniformOutput',false);


% i_00 = setdiff(i_0,find(is_elim_con0));

func_str=[strjoin(u_read,newline) newline newline strjoin(par_read,newline) newline newline...
    strjoin(consrv_assgn,newline) newline elim_con_defs newline newline newline strjoin(model_defs, newline)... 
    newline newline strjoin(predefs, newline) newline newline char(strjoin(strcat(string(f__0(i_0)),'=',string(f_slow)),newline))...
    newline newline strjoin(F_assgn, newline)];
func_str = fortran_subroutine('FUNC','NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP',...
    ['IMPLICIT NONE' newline 'INTEGER NDIM, IJAC, ICP(*)' newline ...
    'DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)'],...
    unique([model_pars, chems(~is_elim) cellstr(consrv_nm), model_vars, string(m1_gamma), cellstr(string(f__0(i_0)))],'stable'),...
func_str);


stpnt_str = fortran_subroutine('STPNT','NDIM,U,PAR,T',...
    ['IMPLICIT NONE' newline 'INTEGER NDIM' newline ...
    'DOUBLE PRECISION U(NDIM), PAR(*), T'],{},...
    [ strjoin(par_assgn,newline) newline newline strjoin(u_assgn,newline)]);

auto_str = [ func_str stpnt_str fortran_subroutine('BCND') fortran_subroutine('ICND') fortran_subroutine('FOPT') fortran_subroutine('PVLS')  ];
auto_str = breaklines(auto_str,80,'&');
auto_str = strrep(auto_str,"^","**");

fid=fopen(strcat(save_dir,filesep,f,'.f90'),'w');
fwrite(fid,auto_str,'char');
fclose(fid);

fid=fopen(strcat(save_dir,filesep,'c.',f),'w');
fwrite(fid,...
['NDIM=   ' int2str(length(ic_ode)) ', NPAR=   ' int2str(length(valid_ind)) ', IPS =   1, IRS =   0, ILP =   1' newline ...
'parnames = {'  char(strjoin(strcat(int2str(find(valid_ind)'),": '",model_pars(1:nnz(valid_ind))',"'"),','))  '}' newline ...
'unames = {' char(strjoin(strcat(int2str((1:nnz(~is_elim_con0 & ~is_fast))'),": '",chems(~is_elim_con0 & ~is_fast)',"'"),','))  '}' newline...
'ICP =  [1]' newline ...
'NTST=  50, NCOL=   4, IAD =   3, ISP =   2, ISW = 1, IPLT= 0, NBC= 0, NINT= 0' newline...
'NMX= 10000000, NPR=  10000, MXBF=  10, IID =   2, ITMX= 8, ITNW= 5, NWTN= 3, JAC= 0' newline...
'EPSL= 1e-7, EPSU = 1e-7, EPSS = 1e-05' newline ...
'DS  =   0.001, DSMIN= 0.0000001, DSMAX=   0.05, IADS=   1' newline ...
'THL =  {11: 0.0}, THU =  {}' newline ...
'UZSTOP={}' newline ...
'UZR={}']...
,'char');
fclose(fid);




preamble={'if length(vox)>1'...
    '[tmp,tmp2]=meshgrid(ir0,vox);'...
    '    I_rx=tmp+tmp2;'...
    'else'...
    '    I_rx=vox+ir0;'...
    'end'...
    'a_c_0=alpha_chem(I_rx);' '' '' };


out=strjoin( [preamble lines {''  'alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);'}],'\n');


fid=fopen(strcat(save_dir,filesep,'update_alpha_chem.m'),'w');
fwrite(fid,out,'char');
fclose(fid);


% tmp=arrayfun(@(i) ['x(vox+[' num2str(find(S_(:,i))'-1) ']*sz)' ],1:length(rate_constants),'UniformOutput',0);
% tmp=arrayfun(@(i) [tmp{i} '=' tmp{i} '+[' num2str(S_(S_(:,i)~=0,i)') '];'] ,1:length(tmp),'UniformOutput',0);
% 
% slow_rx=['if rx==1' newline '    ' tmp{1} newline ....
%     strjoin(arrayfun(@(i) ['elseif rx==' int2str(i)  newline '    ' tmp{i}],2:length(tmp),'UniformOutput',0),newline)...
%     newline 'end' newline];



D=getDiffusionRates(f,chems);



fid=fopen(strcat(save_dir,filesep,'initialize_chem_params.m'),'w');
fwrite(fid,['N_species = ' int2str(length(chems)) ';' newline...
    'N_rx = ' int2str(length(rate_constants)) ';' newline...
    'D = [' num2str(D) '];'],'char');
fwrite(fid,[newline 'N_slow = ' int2str(N_slow) ';' newline...
    ['chems={' char(strjoin(strcat("'",chems,"'"),',')) '};']],'char');
fwrite(fid,[repelem(newline,4)  newline], 'char');

fclose(fid);

extra_files=dir('rxn_files/*.m');
extra_files = {extra_files.name};
if ~isempty(extra_files)
addpath(genpath('rxn_files'));

M=struct();
M.chems = chems;
M.is_fast = is_fast;

for file=extra_files
    file = regexprep(file{1},'\..*','');
    eval([file '(str,M)']);
end

rmpath(genpath('rxn_files'))
end

end

