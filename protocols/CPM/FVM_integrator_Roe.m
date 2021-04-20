%NUMERICS TAKEN FROM:
%J. Math. Biol. (1995) 34:148-176 Journel of
%Mathematical
%Mology
%© Springer-Verlag 1995
%Implicit-explicit methods for reaction-diffusion
%problems in pattern formation 
%%


u = x(cell_inds(1:A) + i_chem_0);

interior=~bndrys&cell_mask(:);

vox0=repmat((1:sz)',1,N_dim);
vox=vox0(interior);


%%


row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
J_diff_0 = sparse(i,j,Delta,numel(interior),sz); %the flux into/out of cell i due to transport across face j
% if exist('reacted','var')
%     if reacted
    
%         vel_1=ones(size(i0))*1;
%         vel_2=zeros(size(i0));
%         vel=zeros([shape 2]);
%         vel(:,:,1)=ones(size(i0))*vmax;
%         vel(:,:,2)=zeros(size(i0));
        vel2=repelem(vcpm,1,1,2);
%         vel_1(bndry_mask)=0;

        opp_dir=reshape(flipud(reshape(1:4,[2,N_dim/2])),[1,N_dim]);
        dim=repelem(1:N_dim/2,2);
        jump_opp=jump(:,opp_dir);
        
        mask_upw = -(1-mod(1:N_dim,2)*2).*vel2(jump+((0:N_dim-1))*sz)<0 & interior;
        row_upw = find(mask_upw);
%         [i_upw,~, ~] =find(mask_upw);
        i_upw = mod(row_upw-1,sz)+1;

        
        j_upw=jump(mask_upw);
        v_upw = vel2(mask_upw);
% 
% %     end  
    J_adv = sparse([row_upw; row_upw],[i_upw; j_upw],[-v_upw; v_upw]*h,numel(interior),sz); %add in the advective flux
%     J_adv = sparse([row_upw; row_upw],[i_upw; j_upw; ],[-v_upw; v_upw]*h,numel(interior),sz); %add in the advective flux

%     
%     vcpm_1 = vcpm(:,:,1);
%     vcpm_2 = vcpm(:,:,2);
%     vcpm(find((vcpm_1(:)>0 & bndrys(:,1)) | (vcpm_1(:)<0 & bndrys(:,2))))=0;
%     vcpm(find((vcpm_2(:)>0 & bndrys(:,3)) | (vcpm_2(:)<0 & bndrys(:,4)))+sz)=0;
    
%     lambda_plus=(vcpm(:)+abs(vcpm(:)))/2;
%     lambda_minus=(vcpm(:)-abs(vcpm(:)))/2;
    
%     lambda_plus_1=(vcpm_1(:)+abs(vcpm_1(:)))/2;
%     lambda_plus_2=(vcpm_2(:)+abs(vcpm_2(:)))/2;
%     lambda_minus_1=(vcpm_1(:)-abs(vcpm_1(:)))/2;
%     lambda_minus_2=(vcpm_2(:)-abs(vcpm_2(:)))/2;
% 
%     
%     u_plus_1 = sparse(1:sz , jump(1:sz,1)',1);
%     u_plus_2 = sparse(1:sz , jump(1:sz,3)',1);
%     u_minus_1 = sparse(1:sz , jump(1:sz,2)',1);
%     u_minus_2 = sparse(1:sz , jump(1:sz,4)',1);
%     
%     D_plus_1 =  (u_plus_1-speye(sz));
%     D_plus_2 =  (u_plus_2-speye(sz));
%     D_minus_1 =  (u_minus_1-speye(sz));
%     D_minus_2 =  (u_minus_2-speye(sz));
% 
%     
%     f_plus_1 = u_plus_1.*vcpm_1(jump(:,1));
%     f_plus_2 = u_plus_2.*vcpm_2(jump(:,3));
%     
%     f_minus_1 = u_minus_1.*vcpm_1(jump(:,2));
%     f_minus_2 = u_minus_2.*vcpm_2(jump(:,4));
% 
%     f_1=spdiags(vcpm_1(:),0,sz,sz);
%     f_2=spdiags(vcpm_2(:),0,sz,sz);
%     
%     f_star_plus_1 = (f_plus_1+f_1)/2-(D_plus_1.*sign(lambda_plus_1))/2;
%     f_star_plus_2 = (f_plus_2+f_2)/2-(D_plus_2.*sign(lambda_plus_2))/2;
%     
%     f_star_minus_1 = (f_minus_1+f_1)/2-(D_minus_1.*sign(lambda_minus_1))/2;
%     f_star_minus_2 = (f_minus_2+f_2)/2-(D_minus_2.*sign(lambda_minus_2))/2;
% %     f_star_minus = (f_minus+f)/2-(spdiags(lambda_minus,0,sz,sz)*D_plus)/2
%     
%     J_adv_2 = -((f_star_plus_1 - f_star_minus_1)+(f_star_plus_2 - f_star_minus_2))/(2*h);
%     %     [ii2,~,~]=find();
%     [ii,~,~]=find(~interior);
% %     J_adv(interior,ii)
%     intersect(ii2,ii)

%     J_adv(:)=0;
    
    
% end



[i2,j2,J_vals ] = find(J_diff_0);
i2=mod(i2-1,sz)+1;
u_xx = sparse(i2,j2,J_vals/h,sz,sz);
u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));
eye = speye(A);

% if exist('reacted','var')
    [i20,j2,J_vals ] = find(J_adv);
    i2=mod(i20-1,sz)+1;
    vu_x = sparse(i2,j2,J_vals/h,sz,sz);
%     vu_x(jump(mask_upw),:) = vu_x(jump(mask_upw),:) - vu_x(jump(jump_opp(mask_upw)),:);
%     vu_x(jump(unique(i20)),:)=vu_x(jump(unique(i20)),:)-vu_x(jump_opp(unique(i20)),:);
%     vu_x=vu_x(vox',vox');
% vu_x=J_adv_2;
% sum(abs(vu_x),'all')
    vu_x=vu_x(cell_inds(1:A),cell_inds(1:A));
%     sum(abs(vu_x),'all')
%     vu_x = -(vu_x-vu_x');

    MAT_list = arrayfun(@(Di)eye+(vu_x+Di*u_xx)*dt,D(1:N_slow),'UniformOutput',0); 
%      MAT_list = arrayfun(@(Di)eye+(vu_x)*dt,D(1:N_slow),'UniformOutput',0); 
% else
%     MAT_list = arrayfun(@(Di)eye+Di*u_xx*dt,D(1:N_slow),'UniformOutput',0); 
% end


eval_Rx_slow
u_prev = u(:,1:N_slow);

t0=time;


% MAT_list = arrayfun(@(Di)( eye/dt-Di*u_xx*9/16),D,'UniformOutput',0); %MCNAB
if any(u(:)<0)
    disp('negatory pig pen')
end

while time-t0<T_integration
    

    
%     Rx_prev=Rx;
    eval_Rx_slow
%     u_curr = u(:,1:N_slow);
%     b_=(2*u_curr-u_prev/2)/dt + (2*Rx-Rx_prev); % 2-SBDF
%     b=u/dt + 3/2*Rx-Rx_prev/2 + (D'.*(u_xx*(3*u/8 + u_prev/16))')'; %MCNAB
%     u_prev=u_curr;
    sumo=sum(u,1);
%      sum(u_xx*u)
    sum(vu_x*u)
    for i = 1:N_slow
        u(:,i) = MAT_list{i}*u(:,i)+Rx(:,i)*dt;
        
%         if any(u(:,i)<0)
%             disp(['wild ass over here -> ' int2str(i)])
%         end
%         if any(any(isnan(u)))
%             disp('welp, look whats for breakfast, its nan bread')
%         end
    end

    
%    [ (sumo-sum(u,1))/dt sum( (sumo-sum(u,1))/dt)]
    
    project_fast
    

    time=time+dt;


end

    if any(u(:)<0)
        error('Negative solutions detected, please use a smaller timestep dt')
    end

%     if all(all(u==0))
%expanding back in to the full system
x(cell_inds(1:A) + i_chem_0) = u(:);
u = reshape(x,[sz ,size(x,3)]);
eval_Rx %re-evaulate on full matrix
