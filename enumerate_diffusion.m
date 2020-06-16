diffusing_species=1:6; %only the first 6 species diffuse


for drx=1:size(jump,2) %itterating over all possible directions
    diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); %finding cites where diffusion is possible
    num_diffuse(drx)=length(ij0(diffuse_mask(drx,:))); %total number of pssoible diffusion reactions in a direction
    ij_diffuse(drx,:)=[ij0(diffuse_mask(drx,:))' zeros(1,(N)*(N)-num_diffuse(drx))];
    %storing diffusible species in an array padded with zeros so its of a constant length
    
    
end
for vox=1:size(diffuse_mask,2)
    num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); %total number of pssoible diffusion reactions in a voxel
end

id0=(diffusing_species-1)*sz;

for drx=1:size(jump,2)
    diffusing_species_sum(drx,:)=sum(x(id0 + ij_diffuse(drx,1:num_diffuse(drx))')); 
    %sums the number of proteins that can diffuse in each direction per latice square 
end