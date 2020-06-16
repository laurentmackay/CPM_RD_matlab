function x =  Alg3(x0,dt,D,h,jump,diffuse_mask,pT,pi,cell_inds,A)

% pT = zeros(size(num_vox_diff,2),size(x0,3));
sz=size(x0,1)*size(x0,2);



% dx=zeros(size(x0));
x=x0;

Nchem=size(x0,3);
ic0=((1:Nchem)-1)*sz;

    
    for i=1:A
        vox=cell_inds(i);
        jump_vox=jump(:,vox);
        for chem=1:Nchem
            ic=ic0(chem);
            pT_curr=pT(vox,chem);
            x0_curr=x0(vox+ic);
            if pT_curr>0 && x0_curr>0
                m=Alg5(x0_curr,pT_curr,rand());
                mi=sample_mi(m,pi(:,vox)')';
                x(jump_vox+ic)=x(jump_vox+ic)+mi;
                x(vox+ic)=x(vox+ic)-m;
            end
        end


    end
end