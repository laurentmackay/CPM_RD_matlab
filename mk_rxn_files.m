function  mk_rxn_files(f)
predef_spatial = {'bndry_mask'}; %pre-defined spatial variables that will be provided at runtime



[chems,S,rates,fast_chems,fast_pair,fast_affinity] = getChemRxns(f);
vars=getInitialized(f);


str=fileread(f);

str=regexprep(str,"%[^\n]*(\n)?",'$1'); %remove block comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,"\'[^\'\n\r]+\'",''); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',''); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',''); %remove function definition

code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
numer=['0-9' code(2:end)];
name=['([a-zA-Z][a-zA-Z_0-9]*)'];

spatial=false(size(vars));
tmp=~spatial;



%figure out which quantities are spatially variable
ass=cellfun(@(x) [x '(?:[ \t\f]*)?=([^\n\r\;]+)'],vars,'UniformOutput',0);
rhs=regexp(str,ass,'tokens');
rhs=cellfun(@(x) x{1},rhs);
spatial_vars=[chems predef_spatial];
ref=@(x) ['(:?[' code ']+|^)(?<var>' x ')(:?[' code ']+|$)'];
while any(tmp~=spatial)
    spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);
    tmp=cell2mat(cellfun(@(x) any(cellfun(@(y) ~isempty(y),regexp(x,spatial_ref,'tokens'))),rhs,'UniformOutput',0));
   
    spatial_vars=unique([spatial_vars vars(tmp)],'stable');
    spatial(tmp&~spatial)=1;

    spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);
    tmp=cell2mat(cellfun(@(x) any(cellfun(@(y) ~isempty(y),regexp(x,spatial_ref,'tokens'))),rhs,'UniformOutput',0));

end

spatial_vars=spatial_vars(length(chems)+1:end);
ref=@(x) ['(?<pre>[' code ']+|^)(?<var>' x ')(?<post>[' code ']+|$)'];
spatial_ref=cellfun(@(x) ref(x),spatial_vars,'UniformOutput',0);

model_defs = regexp(str,ass,'match');
model_split = regexp([model_defs{:}],[name '[ \t\f]*=([^=]+)'],'tokens');
model_params = cellfun(@(x) x{1}{2}, model_split, 'UniformOutput', false);
model_params = regexp(model_params,name,'tokens');
model_params = [model_params{:}];
model_params = [model_params{:}];

model_vars = cellfun(@(x) x{1}{1}, model_split, 'UniformOutput', false);

sym_str= ['syms ' strjoin(unique([model_params model_vars],'stable'),' ') ];

eval([sym_str newline strjoin([model_defs{:}],[';' newline]) ';'])

% model_params = 

lines=[model_defs{:} ...
 arrayfun(@(i) ['alpha_chem(vox+' num2str(i-1) '*sz)=(' rates{i} ')' strjoin(arrayfun(@(j) chemRateString(chems{j},-S(j,i)),1:size(S,1),'UniformOutput',0),'')],1:length(rates),'UniformOutput',0)];

lines=regexprep(lines,spatial_ref,'$<pre>$<var>\(vox\)$<post>');

chem_ref=cellfun(@(x) ref(x),chems,'UniformOutput',0);
chem_rep=arrayfun(@(i) ['$<pre>x\(vox+' num2str(i-1) '*sz\)$<post>'],1:length(chems),'UniformOutput',0);
lines=regexprep(lines,chem_ref,chem_rep);


lines=cellfun(@(x) [x ';'], lines,'UniformOutput',0);


pair_inds = arrayfun(@(i) strcmp(chems,fast_pair{i})',1:length(fast_pair),'UniformOutput',false);
fast_inds = arrayfun(@(i) strcmp(chems,fast_chems{i})',1:length(fast_chems),'UniformOutput',false);

S_fast  = [-[pair_inds{:}]+[fast_inds{:}] [pair_inds{:}]-[fast_inds{:}]]; %only simple reversible complexation is supported right now :(


lcon = logical(null([S S_fast]','r'));
%not super sure how to proceed, I was somewhat unsure at the theoreitcal
%level as well. couldnt express things in a nice way.

%this gives you all the slow variables
is_slow = ~any(S_fast');

mat = diag(is_slow(1:6),0);
mat = eye(6);

J_gamma=zeros(2,6);
J_gamma(1,2) = 1;
J_gamma(1,6) = 0.25;
J_gamma(2,2) = 0.5;
J_gamma(2,6) = 2;


N_fast = length(fast_chems);
N_slow = length(chems)-N_fast;
N_con = size(lcon,2);
m0 = double(lcon(1:N_slow,:)');

for i=1:N_con
    con_fast=lcon(N_slow+1:end,i);
    if any(con_fast)
         m0(i,:)=m0(i,:)+J_gamma(con_fast,:);
    end
end

m1=diag(lcon(1:N_slow,~any(lcon(N_slow+1:end,:))));
m1 = double(m1(any(m1,2),:));


%now we need to repeat this but with a symbolic matrix
null([m1;m0],'r')

rx_lines=cell(size(S,1),1);
for i=1:size(S,1)
    line=['Rx(:,' int2str(i) ')='];
    for j=1:size(S,2)
        if S(i,j)~=0
            rate = [rates{j} strjoin(arrayfun(@(k) chemRateString(['u(' int2str(k) ',:)'],-S(k,j)),1:size(S,1),'UniformOutput',0)) ];
            if S(i,j)>0
                line=[ line '+' int2str(S(i,j)) '*' rate];
            else
                line=[ line int2str(S(i,j)) '*' rate ];
            end
        end
        
    end
    rx_lines{i} = line;
end


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


 tmp=arrayfun(@(i) ['x(vox+[' num2str(find(S(:,i))'-1) ']*sz)' ],1:length(rates),'UniformOutput',0);
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

D=getDiffusionRates(f,chems);

fid=fopen('initialize_chem_params.m','w');
fwrite(fid,['N_species = ' int2str(length(chems)) ';' newline 'N_rx = ' int2str(length(rates)) ';' newline 'D = [' num2str(D) '];'],'char');
fclose(fid);

end

