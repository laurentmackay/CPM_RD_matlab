function  mk_rxn_files(f)
predef_spatial = {'bndry_mask'}; %pre-defined spatial variables that will be provided at runtime



[chems,S,rate_constants,fast_chems,fast_pair,fast_affinity] = getChemRxns(f);
vars=getInitialized(f,true);
chems=regexprep(chems,':','');
chems=regexprep(chems,'-','_');

str=fileread(f);

str=regexprep(str,"%[^\n]*(\n)?",'$1'); %remove block comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,"\'[^\'\n\r]+\'",''); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',''); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',''); %remove function definition




code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
numer=['0-9' code(2:end)];
name=['([a-zA-Z][a-zA-Z_0-9]*)'];


par_def = regexp(str,'(?:\n|^)(?:p|par|param)[ \t\f]([^\n]+)','tokens');

if ~isempty(par_def)
    
    par_def=[par_def{:}];
    
    pars= regexp(par_def,[name '(?:[ \t\f]*)?='], 'tokens' );
    pars=[pars{:}];
    pars=[pars{:}];
    
    par_vals =      regexp(par_def,[name(2:end-1) '(?:[ \t\f]*)?=([^\,\;])+'], 'tokens' );
    par_vals = [par_vals{:}];
    par_vals = [par_vals{:}];
    
    [~,i_vars]=setdiff(vars,pars);
    
    vars = vars(sort(i_vars));
end
spatial=false(size(vars));
tmp=~spatial;



%figure out which quantities are spatially variable
ass=cellfun(@(x) ['(?<![a-zA-Z_0-9])' x '(?:[ \t\f]*)?=([^\n\r\;]+)'],vars,'UniformOutput',0);




model_defs = regexp(str,ass,'match');
if ~isempty(model_defs)
    model_split = regexp([model_defs{:}],[name '[ \t\f]*=([^=]+)'],'tokens');
    model_params = cellfun(@(x) x{1}{2}, model_split, 'UniformOutput', false);
    model_params = regexp(model_params,name,'tokens');
    model_params = [model_params{:}];
    model_params = [model_params{:}];
else
    
    model_split = {};
    model_params = {};
end

fast_params = regexp(fast_affinity,name,'tokens');
fast_params = [fast_params{:}];
if ~isempty(fast_params)
    fast_params = [fast_params{:}];
end

model_params = unique([model_params fast_params rate_constants],'stable');

model_vars = cellfun(@(x) x{1}{1}, model_split, 'UniformOutput', false);









rhs=regexp(str,ass,'tokens');
rhs=cellfun(@(x) x{1},rhs);
spatial_vars=[chems predef_spatial];
ref=@(x) ['(:?[' code ']+|^)(?<var>' x ')(:?[' code ']+|$)'];

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
ref=@(x) ['(?<pre>[' code ']+|^)(?<var>' x ')(?<post>[' code ']+|$)'];
spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);



sym_str= ['syms ' strjoin(unique([model_params model_vars chems],'stable'),' ') ];


pair_inds = arrayfun(@(i) strcmp(chems,fast_pair{i})',1:length(fast_pair),'UniformOutput',false);
fast_inds = arrayfun(@(i) strcmp(chems,fast_chems{i})',1:length(fast_chems),'UniformOutput',false);

S_fast  = [-[pair_inds{:}]+[fast_inds{:}] [pair_inds{:}]-[fast_inds{:}]]; %only simple reversible complexation is supported right now :(
%this gives you all the slow variables
if ~isempty(S_fast)
    is_slow = ~any(S_fast');
else
    is_slow = true(size(chems));
end

N_fast = length(fast_chems);
N_slow = length(chems)-N_fast;

m10= diag(is_slow(1:N_slow),0);
m10 = m10(is_slow,:);
f_slow = cellfun(@(c) ['f_' c], {chems{is_slow}}, 'UniformOutput', false);


chem_str = strjoin(chems(1:N_slow));

chem_str2 = strjoin(chems(1:N_slow),',');
if ~isempty(model_defs)
    def_str=strjoin([model_defs{:}],[';' newline]);
else
    def_str='';
end
    
eval([sym_str ' ' chem_str ' ' strjoin(f_slow,' ') newline def_str ';'])
eval([ 'assume([' strjoin([[model_vars']' model_params]) ']>0)'])

Gamma = eval([ ' [' strjoin( cellfun(@(aff, pair ) ['simplify(' aff '*' pair ',"Steps",10)'],fast_affinity, fast_pair,'UniformOutput',0),';') ']']   );
Gamma_0=Gamma;
Gamma_simp = Gamma;
i=0;
% nm = ['sigma0__' int2str(i)];
% [Gamma_simp,sigma__rep] = subexpr(Gamma_simp,nm);

J_gamma=arrayfun(@(gi) cellfun(@(j) simplify(diff(gi,j),'Steps',10),{chems{:,1:N_slow}},'UniformOutput',0),Gamma_0,'UniformOutput',false)';
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
    J_gamma = big_simp(:,2:end);
    Gamma = big_simp(:,1);
end



fast_pair_sym = eval(['[' strjoin(fast_pair) ']'])';
aff_simp =  Gamma ./ fast_pair_sym;


lines=[model_defs{:} ...
    arrayfun(@(i) ['alpha_chem(vox+' num2str(i-1) '*sz)=(' rate_constants{i} ')' strjoin(arrayfun(@(j) chemRateString(chems{j},-S(j,i)),1:size(S,1),'UniformOutput',0),'')],1:length(rate_constants),'UniformOutput',0)];

lines = regexprep(lines,spatial_ref,'$<pre>$<var>\(vox\)$<post>');
fast_affinity_0=fast_affinity;
fast_affinity = regexprep(fast_affinity,spatial_ref,'$<pre>$<var>\(vox\)$<post>');

chem_ref=cellfun(@(x) ref(x),chems,'UniformOutput',0);
chem_rep=arrayfun(@(i) ['$<pre>x\(vox+' num2str(i-1) '*sz\)$<post>'],1:length(chems),'UniformOutput',0);

lines=regexprep(lines,chem_ref,chem_rep);


lines=cellfun(@(x) [x ';'], lines,'UniformOutput',0);



lcon = logical(null([S S_fast]','r'));
%not super sure how to proceed, I was somewhat unsure at the theoreitcal
%level as well. couldnt express things in a nice way.




lcon_fast = lcon(N_slow+1:end,:);
i_con_fast = any(lcon_fast);
N_con_fast = nnz(i_con_fast);

% N_con = size(lcon,2);
m0 = sym('m0', [N_con_fast N_slow]);
m00 = sym('m00', [N_con_fast N_slow]);
m1 = sym('J_gamma_', [N_con_fast N_slow]);
if ~isempty(J_gamma)
    m0( J_gamma == 0 )=0;
    m00(lcon(1:N_slow,i_con_fast)' == 0 )=0;
    m1( J_gamma == 0) =0;
end


m11=m1;

lcon_fast=lcon_fast(:,i_con_fast);
lmix = lcon(1:N_slow,i_con_fast)';
lslow = lcon(1:N_slow,~i_con_fast)';
for i=1:N_con_fast
    con_fast=lcon_fast(:,i);
    if any(con_fast)
        m0(i,:)=double(lmix(i,:))+J_gamma(con_fast,:);
    end
end

m1_gamma = m1;
m1=m1+double(lmix(1:N_con_fast,:));
b__0 = eval(['[' strjoin( f_slow,'; ') ']'] );
% for i=1:size(lslow,1)
% % irep =  find(any(m10(mask,:).*lslow(i,:),2),1);
% % mask
% % m10(irep,:) = lslow(i,:)
% end

%finding a way to express our system that isnt undetermined
A_mat=[double(lslow(:,is_slow)); m10(:,is_slow)]';
[~,p_A]=rref(A_mat);
A_li_rows = A_mat(:,p_A)';
m10 = zeros(size(m10));
m10(:,is_slow) = A_li_rows;
mask_slow_raw = sum(m10,2)==1;
N_slow2 = length(mask_slow_raw);
i_slow_raw = mod(find(m10(mask_slow_raw,is_slow)')-1,N_slow2)+1;
chems_slow = chems(is_slow);


b__0(mask_slow_raw) =  eval([ '[' strjoin(strcat('f_',chems_slow(i_slow_raw)')) ']']);
b__0(~mask_slow_raw) =  0;

f_mix = eval(['[m10;m1]\[b__0;' strjoin(cellstr(string(zeros(1,size(J_gamma,1)))),'; ')  ']']);
% f_mix2 = eval(['[m10;m0]\[' strjoin( f_slow,'; ') ';' strjoin(cellstr(string(zeros(1,size(J_gamma,1)))),'; ')  ']']);


%
f_fast =  m11*f_mix;
m100 = m10(mask_slow_raw,is_slow);
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
f_mix=f_mix_simp;



f_fast_simp = simplify(f_fast,'Steps',20);
i=0;
nm = ['sub3__' int2str(i)];
[f_fast_simp,sigma3__rep] = subexpr(f_fast_simp,nm);
sigma3__reps={};
if ~isempty(sigma3__rep)
    
    while ~isempty(sigma3__rep)
        sigma3__reps{end+1}=[nm ' = ' char(simplify(sigma3__rep)) ];
        i=i+1;
        nm = ['subs3__' int2str(i)];
        [f_fast_simp, sigma3__rep] = subexpr(f_fast_simp,nm);
        
    end
    
end
% sigma2__reps=sigma2__reps';
f_fast=f_fast_simp;

warning ('off','symbolic:sym:isAlways:TruthUnknown');
non_zero = ~arrayfun(@(m) isAlways(m1_gamma(m)==0) || isAlways(m1_gamma(m)==1),1:numel(m1));
warning ('on','symbolic:sym:isAlways:TruthUnknown');

m1_gamma=m1_gamma(non_zero);
m0=m0(non_zero);

predefs = arrayfun(@(i) [char(m1_gamma(i)) '=' char(m0(i))],1:numel(m1_gamma),'UniformOutput',false);
predefs2 = arrayfun(@(i) ['affinity__' int2str(i) '=' char(aff_simp(i))],1:length(aff_simp),'UniformOutput',false);


rates=cell(size(S,2),1);
for i=1:size(S,2)
    r=-S(:,i);
    r(r<0)=0;
    inds = r>0;
    terms = strcat(chems(inds),cellstr(repmat('.^',[nnz(inds),1]))',cellstr(string(r(inds)))');
    terms = regexprep(terms,'\.\^1$',''); %remove exponent 1's from rate computation
    rate_str = [ rate_constants{i} '.*' strjoin(terms,'.*')];
    rates{i} = rate_str;
end


rates_slow={};
ind_slow = find(is_slow);
for chem_ind = ind_slow(i_slow_raw)
    %     chem_ind = find(m10(i,:),1);
    S_chem = S(chem_ind,:);
    inds = find(S_chem);
    terms = strcat(cellstr(num2str(S_chem(inds)')),repmat('*(',[length(inds),1]),rates(inds),repmat(')',[length(inds),1]));
    terms = regexprep(terms,'[^0-9\.\-\+]*1\*',''); %remove multiplication by 1's from rate computation
    sum_terms = strjoin(terms,'+'); % sum terms
    sum_terms = regexprep(sum_terms,'+-','-');% remove +-, replace with -
    if length(sum_terms)==0
        sum_terms='0';
    end
    rates_slow{end+1} = ['f_' chems{chem_ind} ' = ' sum_terms];
    
end

chem_rep_FVM=arrayfun(@(i) ['$<pre>u\(:,' num2str(i) '\)$<post>'],1:length(chems),'UniformOutput',0);
spatial_rep_FVM=repmat({'$<pre>$<var>\(:\)$<post>'},size(spatial_vars));

rxn_body = [strjoin([[model_defs{:}]'; rates_slow'; sigma__reps'; predefs2'; predefs'; sigma2__reps';],[';' newline])...
    ';' newline newline ['Rx = [' char(strjoin(string(f_mix),[',...' newline])) '];']];
rxn_body = regexprep(rxn_body,chem_ref ,chem_rep_FVM);% reshape the chemical names

rxn_body = regexprep(rxn_body,'(?<pre>[^\.])(?<op>\/|\*|\^)','$<pre>.$<op>'); %make multiplication, division, and exponentiation element-wise operations


% fast_affinity_FVM = regexprep(fast_affinity_0,chem_ref ,chem_rep_FVM);
% fast_affinity_FVM=regexprep(fast_affinity_FVM,spatial_ref,'$<pre>$<var>\(:\)$<post>')
fid=fopen('eval_Rx.m','w');
fwrite(fid,rxn_body,'char');
fclose(fid);

if ~isempty(model_defs)
    model_body = [strjoin([model_defs{:}]',[';' newline]) ';'];
else
    model_body='';
end
chem_rep_full=arrayfun(@(i) ['$<pre>x\(:,:,' num2str(i) '\)$<post>'],1:length(chems),'UniformOutput',0);%this only works for 2d
model_body = regexprep(model_body,chem_ref ,chem_rep_full);% reshape the chemical names

fid=fopen('eval_model.m','w');
fwrite(fid,model_body,'char');
fclose(fid);


ode_body = [strjoin([[model_defs{:}]'; rates_slow'; sigma__reps'; predefs2'; predefs'; sigma2__reps';], newline)...
    ';' newline newline char(strjoin(strcat('d',chems,'/dt=',string(f_mix)'),newline))];



if ~isempty(par_def)
    str_par = ['par' strjoin(strcat(pars,'=',par_vals,';'),', ') newline]
else
    str_par=newline;
end
ode_body=regexprep(ode_body,'\.\*','\*');

fid=fopen([f '.ode'],'w');
fwrite(fid,[ode_body newline newline str_par newline newline],'char');
fclose(fid);




preamble={'if length(vox)>1'...
    '[tmp,tmp2]=meshgrid(ir0,vox);'...
    '    I_rx=tmp+tmp2;'...
    'else'...
    '    I_rx=vox+ir0;'...
    'end'...
    'a_c_0=alpha_chem(I_rx);' '' '' };


out=strjoin( [preamble lines {''  'alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);'}],'\n');


fid=fopen('update_alpha_chem.m','w');
fwrite(fid,out,'char');
fclose(fid);


tmp=arrayfun(@(i) ['x(vox+[' num2str(find(S(:,i))'-1) ']*sz)' ],1:length(rate_constants),'UniformOutput',0);
tmp=arrayfun(@(i) [tmp{i} '=' tmp{i} '+[' num2str(S(S(:,i)~=0,i)') '];'] ,1:length(tmp),'UniformOutput',0);

slow_rx=['if rx==1' newline '    ' tmp{1} newline ....
    strjoin(arrayfun(@(i) ['elseif rx==' int2str(i)  newline '    ' tmp{i}],2:length(tmp),'UniformOutput',0),newline)...
    newline 'end' newline];


i_fast=cellfun(@(x) find(strcmp(chems,x)),fast_chems)-1;
i_pair=cellfun(@(x) find(strcmp(chems,x)),fast_pair)-1;
rx_complex=arrayfun(@(i) find( S(i,:)),i_pair,'UniformOutput',false);

fast_rx=arrayfun(@(i)['if ' strjoin(arrayfun(@(j) ['rx==' int2str(j)],rx_complex{i},'UniformOutput',false),"||")...
    newline '    xtot=sum(x(vox+[' int2str([i_pair(i) i_fast(i)]) ']*sz));' newline ...
    '    x(vox+' int2str(i_pair(i)) '*sz)=round(xtot/(1+' fast_affinity{i} '));' newline '    x(vox+' int2str(i_fast(i)) '*sz)=xtot-x(vox+' int2str(i_pair(i)) '*sz);' newline 'end' ],1:length(rx_complex),'UniformOutput',false);

out=strjoin({slow_rx,strjoin(fast_rx,newline)},newline);



fid=fopen('perform_rx.m','w');
fwrite(fid,out,'char');
fclose(fid);


project = arrayfun(@(i) ['utot = u(:,' int2str(i_pair(i)+1) ')+u(:,' int2str(i_fast(i)+1) ');' newline...
    'u(:,' int2str(i_pair(i)+1) ') = utot/(1+affinity__' int2str(i) ');' newline...
    'u(:,' int2str(i_fast(i)+1) ') = utot - u(:,' int2str(i_pair(i)+1)  ');'], 1:length(fast_pair), 'UniformOutput',0);


project = regexprep(project,'(?<pre>[^\.])(?<op>\/|\*|\^)','$<pre>.$<op>'); %convert to element-wise operations


fid=fopen('project_fast.m','w');
fwrite(fid,strjoin(project,newline),'char');
fclose(fid);




D=getDiffusionRates(f,chems);

fid=fopen('initialize_chem_params.m','w');
fwrite(fid,['N_species = ' int2str(length(chems)) ';' newline...
    'N_rx = ' int2str(length(rate_constants)) ';' newline...
    'D = [' num2str(D) '];'],'char');
fwrite(fid,[newline 'N_slow = ' int2str(N_slow) ';' newline...
    ['chems={' char(strjoin(strcat("'",chems,"'"),',')) '};']],'char');
fwrite(fid,[repelem(newline,4)  newline], 'char');

fclose(fid);


if ~isempty(par_def)
    fid=fopen('model_params.m','w');
    fwrite(fid, strjoin(strcat(pars,'=',par_vals,';'),newline),'char');
    fclose(fid);
end

end

