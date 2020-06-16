function x =  Alg3(x0,dt,D,h,jump,diffuse_mask,pT,pi,cell_inds,A,nrep)


sz=size(x0,1)*size(x0,2);

% for i=1:A
%     vox=cell_inds(i);
%     pT(vox,:) = num_vox_diff(vox)*D*dt/(h^2);
% end

% dx=zeros(size(x0));
x=x0;
for i=1:nrep
for chem=1:size(x0,3)
    ic0=(chem-1)*sz;
    for i=1:A
        vox=cell_inds(i);
        pT_curr=pT(vox,chem);
        x0_curr=x0(vox+ic0);
         jump_curr=jump(:,vox)+ic0;
        if pT_curr>0 && x0_curr>0
            m=Alg5(x0_curr,pT_curr,rand());
            mi=sample_mi(m,pi(:,vox)')';
            x(jump_curr)=x(jump_curr)+mi;
            x(vox+ic0)=x(vox+ic0)-m;
        end
    end    
end
x0=x;
end

end