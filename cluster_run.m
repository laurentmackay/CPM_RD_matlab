clearvars all

B=linspace(0,40,10);
Nrep=2e1;


mkdir 'results'
delete results/*.mat

mk_fun('main')
% try 
parpool()
% catch e
% c = parcluster('local');
% d= c.JobStorageLocation; 
% delete(d); 
% parpool()
% end
parfor k=1:numel(B)*Nrep
    i=floor((k-1)/Nrep)+1;
    j=mod((k-1),Nrep)+1;
%     
    main_fun(B(i),j)
% disp([i,j])
end