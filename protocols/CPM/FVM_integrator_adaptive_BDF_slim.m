%NUMERICS TAKEN FROM:
%J. Math. Biol. (1995) 34:148-176

%Implicit-explicit methods for reaction-diffusion
%problems in pattern formation 
%%


u = x(cell_inds(1:A) + i_chem_0);

interior=~bndrys&cell_mask(:);

vox=repmat((1:sz)',1,N_dim);
vox=vox(interior)';


%%


row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox(:)';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
u_x = sparse(i,j,Delta,numel(interior),sz);

[i2,j2,v ] = find(u_x);
i2=mod(i2-1,sz)+1;
u_xx = sparse(i2,j2,v/h,sz,sz);

%system reduction
u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));




eval_Rx_slow
u_prev = u(:,1:N_slow);

t0=time;

eye = speye(A);

is_1BDF = (max(D)*dt/(h^2))>1/2;

if is_1BDF
    MAT_list = arrayfun(@(Di)( eye/dt-Di*u_xx),D(1:N_slow),'UniformOutput',0); % 1-SBDF
else
    MAT_list = arrayfun(@(Di)( eye*3/(2*dt)-Di*u_xx),D(1:N_slow),'UniformOutput',0); % 2-SBDF
end

% MAT_list = arrayfun(@(Di)( eye/dt-Di*u_xx*9/16),D,'UniformOutput',0); %MCNAB
if any(u(:)<0)
    disp('negatory pig pen')
end

while time-t0<T_integration
    

    if ~is_1BDF
        Rx_prev=Rx;
    end
    eval_Rx_slow
    u_curr = u(:,1:N_slow);
    if is_1BDF
        b_=u_curr/dt + Rx; % 1-SBDF   
    else
        b_=(2*u_curr-u_prev/2)/dt + (2*Rx-Rx_prev); % 2-SBDF  
    end
%     b=u/dt + 3/2*Rx-Rx_prev/2 + (D'.*(u_xx*(3*u/8 + u_prev/16))')'; %MCNAB
    if ~is_1BDF
        u_prev=u_curr;
    end

    for i = 1:N_slow
        u(:,i) = MAT_list{i}\b_(:,i);
        if any(u(:,i)<0)
            disp(['wild ass over here -> ' int2str(i)])
        end
    end
    
    project_fast
    

    time=time+dt;


end

    if any(u(:)<0)
        disp('wild ass over here')
    end
%expanding back in to the full system
x(cell_inds(1:A) + i_chem_0) = u(:);
u = reshape(x,[sz ,size(x,3)]);
eval_Rx %re-evaulate on full matrix
