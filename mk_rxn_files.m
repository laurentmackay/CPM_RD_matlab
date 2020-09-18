function  mk_update_alpha_chem(f)
[chems,S,rates] = getChemRxns(f);
vars=getInitialized(f);


str=fileread(f);

str=regexprep(str,"%[^\n]*(\n)?","$1"); %remove block comments
str=regexprep(str,"\.\.\.\n",""); %remove elipses
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
numer=['0-9' code(2:end)];
name=['([a-zA-Z][a-zA-Z_0-9]*)'];

spatial=false(size(vars));
tmp=~spatial;



%figure out which 
ass=cellfun(@(x) [x '(?:[ \t\f]*)?=([^\n\r\;]+)'],vars,'UniformOutput',0);
rhs=regexp(str,ass,'tokens');
rhs=cellfun(@(x) x{1},rhs);
spatial_vars=chems;
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

lines=regexp(str,ass,'match');
i=1;
strjoin(arrayfun(@(j) chemRateString(chems{j},-S(j,i)),1:size(S,1),'UniformOutput',0),'')

lines=[lines...
 arrayfun(@(i) ['alpha_chem(vox+' num2str(i-1) '*sz)=(' rates{i} ')' strjoin(arrayfun(@(j) chemRateString(chems{j},-S(j,i)),1:size(S,1),'UniformOutput',0),'')],1:length(rates),'UniformOutput',0)];

lines=regexprep([lines{:}],spatial_ref,'$<pre>$<var>\(vox\)$<post>');

chem_ref=cellfun(@(x) ref(x),chems,'UniformOutput',0);
chem_rep=arrayfun(@(i) ['$<pre>x\(vox+' num2str(i-1) '*sz\)$<post>'],1:length(chems),'UniformOutput',0);
lines=regexprep(lines,chem_ref,chem_rep);

lines=cellfun(@(x) [x ';'], lines,'UniformOutput',0);

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

out=['if rx==1' newline '    ' tmp{1} newline ....
strjoin(arrayfun(@(i) ['elseif rx==' int2str(i)  newline '    ' tmp{i}],2:length(tmp),'UniformOutput',0),newline)...
 newline 'end' newline];

fid=fopen('perform_rx.m','w');
fwrite(fid,out,'char');
fclose(fid);

fid=fopen('initialize_chem_params.m','w');
fwrite(fid,['N_species = ' int2str(length(chems)) ';' newline 'N_rx = ' int2str(length(rates)) ';'],'char');
fclose(fid);

end

