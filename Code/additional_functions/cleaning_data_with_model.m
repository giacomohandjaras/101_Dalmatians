%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_voxels_cleaned=cleaning_data_with_model(data_voxels,model_matrix)

subjects=size(data_voxels,1);
voxels=size(data_voxels,2);
timepoints=size(data_voxels,3);

timepoints_model=size(model_matrix,1);
features_model=size(model_matrix,2);

data_voxels_cleaned=nan(subjects,voxels,timepoints);

for sub=1:subjects
    temp=squeeze(data_voxels(sub,:,:));
    temp=temp';
    beta=model_matrix\temp;
    temp_cleaned = temp-(model_matrix*beta);
    data_voxels_cleaned(sub,:,:)=temp_cleaned';
end


end
