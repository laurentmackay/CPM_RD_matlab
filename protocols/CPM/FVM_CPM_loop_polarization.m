T_integration = cpmstep;
keep_running=true;

i_rac = find(strcmp(chems,'Rac')); 
inds=cell_inds(1:A)+sz*(i_rac-1);
d_Rac_0=max(max(x(inds)))-min(min(x(inds)));

while time<Ttot && keep_running
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask); %all cell sites padded with 0s (initially)
    
    while (time-last_time)<Ttot

       
        FVM_integrator_2SBF_slim
 
        
        
        if time>=lastcpm+cpmstep
            
            for kk=1:(2*Per)/cpmsteps %itterates CPM step Per times
%                 disp('doing CPM')
                try
                    CPM_step
                catch err
                    rethrow(err)
                    break;
                end
                
                detect_bndrys
                

            end
            
            enumerate_diffusion %recalcluates diffusable sites
            lastcpm=time;
            cpmcounter=cpmcounter+1;
        end
                
        if time>=lastplot+picstep || time==lastcpm % takes video frames
 
            if cpmcounter==cpmsteps*cpm_wait
              u = reshape(x,[sz ,size(x,3)]);
              pic
              lastplot=time; 
            
                i_rac = find(strcmp(chems,'Rac')); 
                inds=cell_inds(1:A)+sz*(i_rac-1);
                d_Rac=max(max(x(inds)))-min(min(x(inds)));
                if plotting
                    gif
                end
%                 if ~isempty(getCurrentTask())
                    disp([num2str(copyNum) ': B=' num2str(B) ', t=' num2str(time) ', delta_Rac=' num2str(d_Rac)])
%                 end
                
                save_results
                cpmcounter=0;
                
                i_rac = find(strcmp(chems,'Rac')); 
                inds=cell_inds(1:A)+sz*(i_rac-1);
%                 if d_Rac/d_Rac_0<1e-4
%                    keep_running=false;
%                    break;
%                 end
                
            end
            
            
        end
        
    end
    
    last_time=time;
    
    
    
    enumerate_diffusion %recaucluates diffusable sites
end
