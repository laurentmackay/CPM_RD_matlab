function [N_fast,N_slow,S_,chem_rates_slow,chem_ref,chems,code,consrv_eqns,...
consrv_nm,elim,elim_refs,elim_reps,elim_str,f,f__0,f_mix,f_slow,i_0,ij_non_fast,...
is_elim,is_fast,m0,m1_gamma,mask_non_elim,model_defs,model_par_vals,model_pars,...
model_var_defs,model_vars,name,par_defs,par_refs,par_reps,predefs,...
rate_constants,ref,rxn_rates,sigma2__reps,sigma__reps,simplify,slow_refs,...
sol_fast_0,sol_rhs,spatial_vars] = consume_fmix_fun(N_fast,N_slow,S_,...
chem_rates_slow,chem_ref,chems,code,consrv_eqns,consrv_nm,elim,elim_refs,elim_reps,...
elim_str,f,f__0,f_mix,f_slow,i_0,ij_non_fast,is_elim,is_fast,m0,m1_gamma,...
mask_non_elim,model_defs,model_par_vals,model_pars,model_var_defs,model_vars,name,...
par_defs,par_refs,par_reps,predefs,rate_constants,ref,rxn_rates,sigma2__reps,...
sigma__reps,simplify,slow_refs,sol_fast_0,sol_rhs,spatial_vars)

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


fid=fopen([f '.rxns'],'w');
fwrite(fid,elim_str,'char');
for i=1:length(paired)
    if paired(i)~=-1
        if ~isnan(paired(i))
            sl=['R' int2str(c) '=' rxn_rates{i} ' - ' rxn_rates{paired(i)} ';' newline];
            paired(paired(i))=-1;
        else
            sl=['R' int2str(c) '=' rxn_rates{i} ';' newline];
        end
        sl=regexprep(sl,'\.\*','*');
        sl=regexprep(sl,'\_','');
        fwrite(fid,sl,'char');
        c=c+1;

    end
    
end
fclose(fid);



chem_rep_FVM=arrayfun(@(i) ['$<pre>u\(:,' num2str(i) '\)$<post>'],1:length(chems),'UniformOutput',0);
spatial_rep_FVM=repmat({'$<pre>$<var>\(:\)$<post>'},size(spatial_vars));

model_defs = regexprep(model_defs,elim_refs,elim_reps);

rxn_body = [strjoin([model_defs'; chem_rates_slow'; sigma__reps'; predefs'; sigma2__reps';],[';' newline])...
    ';' newline newline ['Rx = [' char(strjoin(string(f_mix),[',...' newline])) '];']];
rxn_body = regexprep(rxn_body,chem_ref ,chem_rep_FVM);

rxn_body = regexprep(rxn_body,'(?<pre>[^\.])(?<op>\/|\*|\^)','$<pre>.$<op>'); 


fid=fopen('eval_Rx.m','w');
fwrite(fid,rxn_body,'char');
fclose(fid);



if ~isempty(model_defs)
    model_body = [strjoin(model_defs',[';' newline]) ';'];
else
    model_body='';
end
chem_rep_full=arrayfun(@(i) ['$<pre>x\(:,:,' num2str(i) '\)$<post>'],1:length(chems),'UniformOutput',0);
model_body = regexprep(model_body,chem_ref ,chem_rep_full);

fid=fopen('eval_model.m','w');
fwrite(fid,model_body,'char');
fclose(fid);


init=getInitialConditions(f,chems);

slow_init_reps = cellstr(string(init(~is_fast)));
if N_fast>0
    fast_init_reps = regexprep(cellstr(string(sol_fast_0(1:N_fast))),slow_refs, slow_init_reps);
    fast_init_reps = regexprep(fast_init_reps, par_refs, par_reps);
else
    fast_init_reps = {};
end

consrv_rhs = regexprep(string(consrv_eqns),'^.*==','');
slow_refs = strcat('(?<=(?:^|[',code,']))(',chems(~is_fast),')(?=(?:$|[',code,']))')';
fast_refs = strcat('(?<=(?:^|[',code,']))(',chems(is_fast),')(?=(?:$|[',code,']))')';



consrv_defs_init = regexprep(consrv_rhs,[slow_refs; fast_refs],[slow_init_reps'; fast_init_reps]');
consrv_defs_init = string(eval(['[' char(strjoin(consrv_defs_init)) ']'])');
consrv_def_vals = consrv_defs_init;






consrv_assgn = cellfun(@(nm,v) [nm '=' v],cellstr(consrv_nm), cellstr(consrv_defs_init'),'UniformOutput',false);

consrv_defs_init = strcat("par ",consrv_nm, " = " , consrv_defs_init', string(newline));





f_mix_elim = subs(f_mix,m1_gamma,m0);


f_mix_elim = simplify(subs(f_mix_elim,f__0(ij_non_fast),f_slow),'steps',10);

f_mix_elim = simplify(subs(f_mix_elim,str2sym(model_vars),str2sym(model_var_defs)),'steps',10);
f_mix_elim = simplify(subs(f_mix_elim),'steps',10);
f_mix_elim = simplify(subs(f_mix_elim,str2sym(elim),cell2sym(sol_rhs)),'steps',10);

ode_body_elim = char(strjoin(strcat('d',chems(~is_elim),'/dt=',string(f_mix_elim(mask_non_elim))'),newline));


chem_ref_elim=cellfun(@(x) ref(x),chems(~is_elim),'UniformOutput',0);
chem_rep_elim=arrayfun(@(i) ['$<pre>x\(vox+' num2str(i-1) '*sz\)$<post>'],1:length(chems(~is_elim)),'UniformOutput',0);
chem_rep_FVM_elim=arrayfun(@(i) ['$<pre>u\(:,' num2str(i) '\)$<post>'],1:length(chems(~is_elim)),'UniformOutput',0);

rxn_body_elim = [strjoin([model_defs'; chem_rates_slow'; sigma__reps'; predefs'; sigma2__reps';],[';' newline])...
    ';' newline newline ['Rx = [' char(strjoin(string(f_mix(mask_non_elim)),[',...' newline])) '];']];
rxn_body_elim = regexprep(rxn_body_elim,elim_refs,elim_reps);
rxn_body_elim = regexprep(rxn_body_elim,chem_ref_elim,chem_rep_FVM_elim);

rxn_body_elim = regexprep(rxn_body_elim,'(?<pre>[^\.])(?<op>\/|\*|\^)','$<pre>.$<op>'); 



fid=fopen('eval_Rx_elim.m','w');
fwrite(fid,rxn_body_elim,'char');
fclose(fid);

if ~isempty(par_defs)
    model_par_vals = string(eval(['[' strjoin(model_par_vals) ']']));
    par_str0 = [char(strjoin(strcat(model_pars,"=",model_par_vals),[';' newline])) newline];
    par_str = ['par ' char(strjoin(strcat(model_pars,"=",model_par_vals),', ')) newline];
else
    par_str=newline;
end

init_str = ['init ' char(strjoin(strcat(chems(~is_elim)', "=", string(init(~is_elim)')), ', '))];


ode_body_elim=regexprep(ode_body_elim,'\.\*','\*');


consrv_defs = regexprep(string(consrv_eqns),'==','=');
consrv_str = strjoin(strcat(strcat('#', consrv_defs), string(newline), consrv_defs_init'), newline);
consrv_str = char(consrv_str);


ode_str=[ode_body_elim newline newline ...
    init_str newline newline ...
    consrv_str newline newline ...
    par_str newline newline ];

if ~isempty(par_defs)
    fid=fopen('model_params.m','w');
    fwrite(fid, strjoin(strcat([model_pars cellstr(consrv_nm)],'=',[model_par_vals consrv_def_vals'],';'),newline),'char');
    fclose(fid);
end


nms_elim = regexp(string(f_mix_elim(mask_non_elim)),name, 'tokens')';
nms_elim=[nms_elim{:}];

params = inline_script('model_params');
eval_rhs_str = inline_script('eval_Rx_elim');

str =['function Rx = rhs_fun(t,u)' newline 'u=transpose(u);' newline params newline eval_rhs_str newline 'Rx=transpose(Rx);' newline 'end'];

fid=fopen('rhs_fun.m','w');
fwrite(fid,str,'char');
fclose(fid);


ic_ode=init(~is_fast);
ic_ode=init(~is_elim);
tol=1e-14;
fprintf(['Integrating to find fixed point (abstol = ' num2str(tol) ')...']);
t0=tic();
T_vec=0;
while ~all(abs(rhs_fun(0,ic_ode'))<tol)
    [T_vec,Y_vec] = ode15s(@ rhs_fun,T_vec(end)+[0 1e3],ic_ode,odeset('NonNegative',1:length(ic_ode)));
    ic_ode = Y_vec(end,:);
end

t1=toc(t0);
fprintf(['Done (' num2str(t1) ' seconds elapsed, final time = ' num2str(T_vec(end)) ')' newline])
init_relaxed = ic_ode;
disp(strjoin(strcat(chems(~is_elim),'=',string(init_relaxed)),', '));

u_read = cellfun(@(c,i) [c '=U(' int2str(i) ')'],chems(~is_elim),num2cell(1:nnz(~is_elim)),'UniformOutput',false);

if length(model_pars) > 10
    valid_ind = true(1,min(length(model_pars)+1,36));
    valid_ind(11) = false;
else
    valid_ind = true(1,length(model_pars));
end
par_read = cellfun(@(p,i) [p '=PAR(' int2str(i) ')'],model_pars(1:nnz(valid_ind)),num2cell(find(valid_ind)),'UniformOutput',false);
par_read = [par_read ...
    cellfun(@(p,v) [p '=' v],model_pars(nnz(valid_ind)+1:end),model_par_vals(nnz(valid_ind)+1:end),'UniformOutput',false)];

par_assgn = cellfun(@(v,i) ['PAR(' int2str(i) ')=' v],model_par_vals(1:nnz(valid_ind)),num2cell(find(valid_ind)),'UniformOutput',false);
u_assgn = cellfun(@(v,i) ['U(' int2str(i) ')=' v],cellstr(string(init_relaxed)),num2cell(1:nnz(~is_elim)),'UniformOutput',false);

F_assgn = arrayfun(@(f,i)[  'F(' int2str(i) ')=' char(string(f))],f_mix(mask_non_elim),(1:nnz(mask_non_elim))','UniformOutput',false);


func_str=[strjoin(u_read,newline) newline newline strjoin(par_read,newline) newline newline...
    strjoin(consrv_assgn,newline) newline newline strjoin(model_defs, newline)... 
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

fid=fopen([f '.f90'],'w');
fwrite(fid,auto_str,'char');
fclose(fid);


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
fid=fopen([f '.ode'],'w');
fwrite(fid,ode_str,'char');
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


tmp=arrayfun(@(i) ['x(vox+[' num2str(find(S_(:,i))'-1) ']*sz)' ],1:length(rate_constants),'UniformOutput',0);
tmp=arrayfun(@(i) [tmp{i} '=' tmp{i} '+[' num2str(S_(S_(:,i)~=0,i)') '];'] ,1:length(tmp),'UniformOutput',0);

slow_rx=['if rx==1' newline '    ' tmp{1} newline ....
    strjoin(arrayfun(@(i) ['elseif rx==' int2str(i)  newline '    ' tmp{i}],2:length(tmp),'UniformOutput',0),newline)...
    newline 'end' newline];



D=getDiffusionRates(f,chems);



fid=fopen('initialize_chem_params.m','w');
fwrite(fid,['N_species = ' int2str(length(chems)) ';' newline...
    'N_rx = ' int2str(length(rate_constants)) ';' newline...
    'D = [' num2str(D) '];'],'char');
fwrite(fid,[newline 'N_slow = ' int2str(N_slow) ';' newline...
    ['chems={' char(strjoin(strcat("'",chems,"'"),',')) '};']],'char');
fwrite(fid,[repelem(newline,4)  newline], 'char');

fclose(fid);
end
