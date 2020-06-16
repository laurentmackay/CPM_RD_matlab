function x =  Alg3(x0,dt,D,h,jump,diffuse_mask,pT,pi,cell_inds,A)

sz=size(x0,1)*size(x0,2);
x=x0;


for chem=1:length(D)
    ic0=(chem-1)*sz;
    for i=1:A
        vox=cell_inds(i);
        pT_curr=pT(vox,chem);
        x0_curr=x0(vox+ic0);
        neighbors=jump(:,vox)+ic0;
        if pT_curr>0 && x0_curr>0
            m=Alg5(x0_curr,pT_curr,rand());
            mi=sample_mi(m,pi(:,vox)')';
            x(neighbors)=x(neighbors)+mi;
            x(vox+ic0)=x(vox+ic0)-m;
        end
    end    
    end
end