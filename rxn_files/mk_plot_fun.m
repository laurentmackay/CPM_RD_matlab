function mk_plot_fun(str,M,save_dir)

code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
numer=['0-9' code(2:end)];
name=['[a-zA-Z][a-zA-Z_0-9]*'];


plot_defs = regexp(str,['(?<=plot\()((?:' name '[ \t\f]*,[ \t\f]*)' name '[ \t\f]*)+\)'],'tokens');

slow_chems = M.chems(~M.is_fast);
if isempty(plot_defs)
    plot_defs = slow_chems;
else
    plot_defs = regexp(plot_defs{1},'[ \t\f]*,[ \t\f]*','split');
    plot_defs = [plot_defs{:}];
end

N_plot = length(plot_defs);

m = ceil(sqrt(N_plot ));
n = ceil(N_plot/m);
nl = string(newline);

init_str = ['pic_fig=figure(1);clf();' newline strjoin(cellstr(strcat('panel',string((1:N_plot)'),'=subplot(',int2str(m),',',int2str(n),',',int2str((1:N_plot)'),');')), newline)];
fid = fopen(strcat(save_dir,'/initialize_pic.m'),'w');


fwrite(fid,init_str,'char');

fclose(fid);

inds = (1:N_plot)';

plot_str=compose(strcat('plotCellIm(panel%i,reshape(%s,shape),cell_mask,i0,j0);', nl,...
              "caxis(panel%i,'auto');",nl,...
              "colorbar(panel%i);",nl,...
              "title(panel%i,'%s', 'Fontsize', 24);",newline),inds,string(plot_defs'),inds,inds,inds,string(plot_defs'));

plot_str = regexprep(plot_str,nameref(slow_chems),cellstr(strcat('u(:,',int2str((1:length(slow_chems))'),')'))); %replace literal chem names by their state in the `u` array
          
fid = fopen(strcat(save_dir,'/pic.m'),'w');
fwrite(fid,['if plotting' newline newline],'char');
fwrite(fid,['tp__0=tic;' newline newline],'char');
fwrite(fid,strjoin(plot_str,nl+nl),'char');
if exist('sgtitle')~=0
fwrite(fid,...    
       nl+nl+"sgtitle(pic_fig,['t=' num2str(time) ', t_{plot}=' num2str(double(tic-tp__0)*1e-6), ', t_{sim}=' num2str(toc)], 'Fontsize', 10,'FontWeight','bold')",...
       'char');
end

fwrite(fid,[newline newline 'end'],'char');
fclose(fid);
    
          
end

