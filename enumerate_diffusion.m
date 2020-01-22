diffusing_species=1:6;

ij_diffuse=cell([1,size(jump,2)]);


for drx=1:size(jump,2)
diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:);
ij_diffuse{drx}=ij0(diffuse_mask(drx,:));
end
