function [chem_rates,rxn_rates] = rate_strings(chems,r,p,rate_constants)

rxn_rates=cell(size(r,2),1);

for i=1:size(r,2)
%     r=r(:,i);
%     r(r<0)=0;
    inds = r(:,i)>0;
    terms = strcat(chems(inds),cellstr(repmat('.^',[nnz(inds),1]))',cellstr(string(r(inds,i)))');
    terms = regexprep(terms,'\.\^1$',''); %remove exponent 1's from rate computation
    rate_str = [ rate_constants{i} '.*' strjoin(terms,'.*')];
    rxn_rates{i} = rate_str;
end
S=p-r;
chem_rates=cell(size(S,1),1);
for chem_ind = 1:length(chem_rates)
    %     chem_ind = find(m10(i,:),1);
    S_chem = S(chem_ind,:);
    inds = find(S_chem);
    terms = strcat(cellstr(num2str(S_chem(inds)')),repmat('*(',[length(inds),1]),rxn_rates(inds),repmat(')',[length(inds),1]));
    terms = regexprep(terms,'[^0-9\.\-\+]*1\*',''); %remove multiplication by 1's from rate computation
    sum_terms = strjoin(terms,'+'); % sum terms
    sum_terms = regexprep(sum_terms,'+-','-');% remove +-, replace with -
    if length(sum_terms)==0
        sum_terms='0';
    end
    chem_rates{chem_ind} = sum_terms;
    
end



end

