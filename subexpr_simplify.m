function [simp,reps] = subexpr_simplify(eq,base_nm)
i=0;
nm = [char(base_nm) int2str(i)];
[simp,rep] = subexpr(eq,nm);
reps={};
while ~isempty(rep)
    reps{end+1}=[nm ' = ' char(simplify(rep)) ];
    i=i+1;
    nm = [char(base_nnm) int2str(i)];
    [simp, rep] = subexpr(simp,nm);
end

end

