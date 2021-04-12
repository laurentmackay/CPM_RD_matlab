% file='../_chem_Rx_Pax_Kathy/results/random_walk_B_sweep/final_B_2.25_copy1.mat';
% file='final_B_6.2586_copy2.mat';
deploy_model('chem_Rx_Pax_Asheesh')

set_experiment('random_walk_B_sweep_200k2')



%% find the results and their B values
B_vals = {};
results = ls_results()';
for result = results
    B_str = regexp(result.name,'B_([^\_]+)','tokens');
    B_vals{end+1}=str2num(B_str{1}{1});
end

B_vals_res = cell2mat(B_vals);
[B_vals, ~ ,inds] = unique(B_vals_res);
disp([ int2str(length(inds)) ' results total, ' int2str(length(B_vals)) ' unique B values'])




%
thresh=1.5;
% percent_moving=cell(length(B_vals),1);
% lambda=cell(length(B_vals),1);
% persistence_time=cell(length(B_vals),1);
figure(1);

for i=1:length(inds)

[acf_halflife, t_acf] = get_acf_halflife(results(i).name);
[t_polarize, ~] = get_polarization_time(results(i).name);

i_polarize=find(t_acf>=t_polarize,1);
t_acf(i_polarize)
clf();
plot(t_acf,acf_halflife)
xlabel('Time');
ylabel(['\lambda (B=' num2str(B_vals(inds(i))) ')'])
drawnow
    if ~isempty(i_polarize)
        N_up=nnz(acf_halflife(i_polarize:end)>thresh);
%         percent_moving{inds(i)}{end+1}=N_up/(length(acf_halflife)-i_polarize);
%         lambda{inds(i)}{end+1}=acf_halflife(i_polarize:end);
    end
end


%%
% 
% block_length(b) 
% 

lambda{cellfun(@isempty, lambda)}={};
figure(2);clf();
    hold on
    N_bin=100;

    counts=zeros(N_bin,length(B_vals));
    edges=zeros(N_bin+1,length(B_vals));
    up_times=cellfun(@(x) cellfun(@(y)  get_down_times(y,thresh), x, 'UniformOutput',0),lambda, 'UniformOutput', 0);
    up_times=[up_times{:}];
    up_times=[up_times{:}];
    edges_0=linspace(3.5, 500.5, N_bin+1);
%     [~,edges_0]=histcounts(up_times, );
    
    %     cellfun(@(x) max(cellfun(@(y) max( get_up_times(y,thresh)), x)),lambda, 'UniformOutput', 0)
for i=1:length(lambda)
  

    
    up_times=cellfun(@(x) get_down_times(x,1.5), lambda{i},'UniformOutput',0);
    [counts(:,i),edges(:,i)]=histcounts([up_times{:}], edges_0,'normalization','pdf');
%     plot(, counts)
    
    xlim([0, 200])
%     for j=lambda{i}
%         
%         
%     
%     plot(j{1});
%     
%     
%     end

%     drawnow
end
hold off
[x,y]=meshgrid(B_vals,(edges_0(1:end-1)+edges_0(1:end-1))/2);
pcolor(x,y,counts)
%%
percent_moving{cellfun(@isempty, percent_moving)}={};
p_move = cellfun(@(x) mean([x{:}]), percent_moving)
figure(1);clf();
plot(p_move);



