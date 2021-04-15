
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


%% compute the timecourse of polarization for each result
max_diff=cell(length(B_vals),1);
time_to_polarize=cell(length(B_vals),1);

for i=1:length(inds)
    
    [t,y]=timecourse(results(i).name,"[min(min(x(inds))), max(max(x(inds)))]",...
        "i_rac = find(strcmp(chems,'Rac')); inds=cell_inds(1:A)+sz*(i_rac-1); ");
    
    delta=abs( y(:,1)-y(:,2));
    delta_max=max(delta);
    max_diff{inds(i)}{end+1} = delta_max;
    
    if delta_max>0.01
        time_to_polarize{inds(i)}{end+1}=t(find(delta>=delta_max/2,1));
    else
        time_to_polarize{inds(i)}{end+1}=Inf;
    end
    disp(['processed: ' int2str(i) '/' int2str(length(inds)) ', B=' num2str(B_vals(inds(i)))  '->' num2str(delta_max>0.01)])
end
disp('done running timecourse')

%% plotting
e2=cellfun(@(x) std([x{:}]),time_to_polarize);
t2=cellfun(@(x) mean([x{:}]),time_to_polarize);
figure(1);clf();
h=errorbar(B_vals,t2,e2);

% hold (gca,'on')
% figure(1);
% h=errorbar(B_vals,t2,e2);
% hold (gca,'off')



d2=cellfun(@(x) mean([x{:}]),max_diff);
figure(2);clf();
plot(B_vals,d2)
% set(gca,'YScale','Log')