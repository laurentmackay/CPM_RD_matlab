diffusing_species=1:N_species; %only the first 6 species diffuse


for drx=1:size(jump,2) %itterating over all possible directions
    diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); %finding cites where diffusion is possible
    
end

for vox=1:size(diffuse_mask,2)
    num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); %total number of pssoible diffusion reactions in a voxel
end

%pre-compute single molecule diffusion probabilities for the timestep of
%length dt, as well as the directional probabilities
%this must be re-run everytime cell-geometry changes
for i=1:A
    vox=cell_inds(i);
    pT0(vox,:) = num_vox_diff(vox)*D/(h^2);
    pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
end