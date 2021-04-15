function x =  Alg3(x,dt,D,h,jump,diffuse_mask,pT0,pi,cell_inds,A,inds)

sz=size(x,1)*size(x,2);
% x=x0;

m_avg=0;
counter=0;

for chem=1:length(inds)
    if inds(chem)
        ic0=(chem-1)*sz;
        for i=1:A
            
            vox=cell_inds(i);
            pT_curr=pT0(vox,chem);
            x0_curr=x(vox+ic0);
            neighbors=jump(:,vox)+ic0;
            
            if pT_curr>0 && x0_curr>0
                m=Alg5(x0_curr,pT_curr*dt(chem),rand());
                %           disp([pT_curr*dt*x0_curr, m])
                m_avg=m_avg+m;
                counter=counter+1;
                mi=sample_mi(m,pi(:,vox)')';
                x(neighbors)=x(neighbors)+mi;
                x(vox+ic0)=x(vox+ic0)-m;
                
            end
        end
    end
end
% disp(m_avg/counter);
end