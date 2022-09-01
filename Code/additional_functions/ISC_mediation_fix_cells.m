%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results_fixed=ISC_mediation_fix_cells(results,voxels_split)

[permutations,chunks]=size(results);
temp_data=results{1,1};
[subjects_pairings]=size(temp_data,1);
voxel_count=max(voxels_split)-1;
voxels_split_num=numel(voxels_split);

results_fixed=cell(permutations,1);

for p=1:permutations
    temp_data=nan(subjects_pairings,voxel_count);
    for chunk=1:voxels_split_num-1
        voxel_begin=voxels_split(chunk);
        voxel_end=voxels_split(chunk+1)-1;
        temp_data(:,voxel_begin:voxel_end)=results{p,chunk};
    end
    results_fixed{p}=temp_data;
    clear temp_data
end


end
