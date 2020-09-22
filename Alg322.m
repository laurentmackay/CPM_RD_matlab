function x =  Alg322(x,dt,D,h,jump,diffuse_mask,pT0,pi,cell_inds,A,inds)

sz=size(x,1)*size(x,2);
% x=x0;
m=zeros(A,1);
vox=cell_inds(1:A);

for chem=1:length(inds)
    if inds(chem)
        ic0=(chem-1)*sz;
        
        pT_curr=pT0(vox,chem);
        x0_curr=x(vox+ic0);
        
        
        inds_diff=pT_curr>0 & x0_curr>0;
%         vox_diff=vox(inds_diff);

        

        
        
        m(inds_diff)=Alg52(x0_curr(inds_diff),pT_curr(inds_diff)*dt(chem));
        mi=sample_mi2(m(inds_diff),pi(:,vox(inds_diff)));
        vox_diff=vox(inds_diff);
        
        for j=1:length(vox_diff)
            neighbors=jump(:,vox_diff(j))+ic0;
            x(neighbors)=x(neighbors)+mi(:,j);
%             x(vox_diff(j)+ic0)=x(vox_diff(j)+ic0)-m(j);
        end
%         try
        x(vox_diff+ic0)=x(vox_diff+ic0)-m(inds_diff);
%         catch e
%             throw(e);
%         end
        %update x
        
        
    end
end
% disp(m_avg/counter);
end