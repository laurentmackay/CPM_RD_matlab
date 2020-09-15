%wow diffusing
function x =  Alg3(x,dt,jump,pT0,pi,cell_inds,A,inds)

sz=size(x,1)*size(x,2);

m=zeros(A,1);
vox=cell_inds(1:A);

for chem=1:length(inds)
    if inds(chem)
        ic0=(chem-1)*sz;
        
        pT_curr=pT0(vox,chem);
        x0_curr=x(vox+ic0);
        inds_diff=pT_curr>0 & x0_curr>0;
        
        %determine number of molecules to be transferred
        m(inds_diff)=Alg5(x0_curr(inds_diff),pT_curr(inds_diff)*dt(chem));
        mi=sample_mi(m(inds_diff),pi(:,vox(inds_diff)));
        vox_diff=vox(inds_diff);
        
        for j=1:length(vox_diff)
            neighbors=jump(:,vox_diff(j))+ic0;
            x(neighbors)=x(neighbors)+mi(:,j);
        end
        % the following if check is necessary when there exists null
        % molecule in a square
        if ~isempty(vox_diff)
            x(vox_diff+ic0)=x(vox_diff+ic0)-m(inds_diff);
        end

    end
end
end