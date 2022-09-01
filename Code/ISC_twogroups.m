%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

addpath('NIfTI_tools/');
addpath('additional_functions/');

filename_mask='../Results/mask_grey_matter.nii.gz'; %%%voxels to be considered

%%%%%Files to be saved
save_file='../Results/Results_ISC_twogroups_SD.nii'; %the nifti
save_file_mat='../Results/Results_ISC_twogroups_SD.mat'; %the workspace

groups='SD'; %%%groups in our experiment: 'TD' for Tipically developed, 'SD' for sensory deprived

switch groups
    case 'TD'
        %%%%Video-only
        filenames_groupA={'../fMRI_101_Dalmatians/ctrlV/sub-020.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-021.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-023.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-024.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-025.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-026.nii.gz','../fMRI_101_Dalmatians/ctrlV/sub-028.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-029.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-030.nii.gz', '../fMRI_101_Dalmatians/ctrlV/sub-031.nii.gz'};
        %%%%Audio-only
        filenames_groupB={'../fMRI_101_Dalmatians/ctrlA/sub-003.nii.gz'    , '../fMRI_101_Dalmatians/ctrlA/sub-004.nii.gz'    , '../fMRI_101_Dalmatians/ctrlA/sub-005.nii.gz'  , '../fMRI_101_Dalmatians/ctrlA/sub-006.nii.gz'  , '../fMRI_101_Dalmatians/ctrlA/sub-007.nii.gz', '../fMRI_101_Dalmatians/ctrlA/sub-008.nii.gz', '../fMRI_101_Dalmatians/ctrlA/sub-009.nii.gz', '../fMRI_101_Dalmatians/ctrlA/sub-010.nii.gz', '../fMRI_101_Dalmatians/ctrlA/sub-011.nii.gz', '../fMRI_101_Dalmatians/ctrlA/sub-027.nii.gz'};
        
    case 'SD'
        %%%%Deaf
        filenames_groupA={'../fMRI_101_Dalmatians/deaf/sub-044.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-045.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-046.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-047.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-048.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-049.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-050.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-051.nii.gz', '../fMRI_101_Dalmatians/deaf/sub-052.nii.gz' };
        %%%%Blind
        filenames_groupB={'../fMRI_101_Dalmatians/blind/sub-033.nii.gz', '../fMRI_101_Dalmatians/blind/sub-035.nii.gz', '../fMRI_101_Dalmatians/blind/sub-036.nii.gz', '../fMRI_101_Dalmatians/blind/sub-038.nii.gz', '../fMRI_101_Dalmatians/blind/sub-039.nii.gz', '../fMRI_101_Dalmatians/blind/sub-041.nii.gz', '../fMRI_101_Dalmatians/blind/sub-042.nii.gz', '../fMRI_101_Dalmatians/blind/sub-043.nii.gz', '../fMRI_101_Dalmatians/blind/sub-053.nii.gz' };
end


timepoints_stimulus=1614;
randomA=load('permutation_schemas/permutation_schema_final_1000perm.mat');
randomB=load('permutation_schemas/permutation_schema_final_1000perm_alt.mat');


%%%%% number of CPU cores to be used %%%%%
CPUs=4;

%%%%% The maximum number of NaN across subject's pairings allowed in a voxels is no more than 1/tollerance. NaN depends on voxels with a timeseries of zeros.%%%%%
tollerance=3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's open the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s open the mask...'));
disp(sprintf('######################################################'));

disp(sprintf('Open nifti mask %s', filename_mask));
mask=load_nii(filename_mask);

voxel_size=mask.hdr.dime.pixdim([2:4]);
disp(sprintf('Voxel size %d x %d x %d mm', voxel_size(1), voxel_size(2), voxel_size(3)));

x_size=size(mask.img,1);
y_size=size(mask.img,2);
z_size=size(mask.img,3);
disp(sprintf('Matrix size %d x %d x %d', x_size, y_size, z_size));


%%%Retrieve voxels of interests
voxel_count=0;
coordinates=[];
for x=1:x_size
    for y=1:y_size
        for z=1:z_size
            if mask.img(x,y,z)>0
                voxel_count=voxel_count+1;
                coordinates(voxel_count,:)=[x,y,z];
            end
        end
    end
end

disp(sprintf('Voxels in the mask %d', voxel_count));

permutations=size(randomA.permutation_tps,1);
disp(sprintf('Number of permutations %d', permutations));

subjects_groupA=numel(filenames_groupA);
disp(sprintf('Number of subjects group A %d', subjects_groupA));

subjects_groupB=numel(filenames_groupB);
disp(sprintf('Number of subjects group B %d', subjects_groupB));


%%%%%Let's open all the subjects
disp(sprintf('\n######################################################'));
disp(sprintf('Let''s perform the analysis for group %s', groups));
disp(sprintf('######################################################'));

temp_data_groupA=nan(voxel_count,subjects_groupA,timepoints_stimulus);
temp_data_groupB=nan(voxel_count,subjects_groupB,timepoints_stimulus);

for sub=1:subjects_groupA
    current_sub=filenames_groupA{sub};
    disp(sprintf('\nOpen subject of group A %d: %s', sub, current_sub));
    
    current_sub_data=load_nii(current_sub);
    t_size=size(current_sub_data.img,4);
    disp(sprintf('Timepoints subject %d', t_size));
    
    for voxel=1:voxel_count
        temp_data_groupA(voxel,sub,:)=squeeze(current_sub_data.img(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),:));
    end
    
    clear current_sub_data
end


for sub=1:subjects_groupB
    current_sub=filenames_groupB{sub};
    disp(sprintf('\nOpen subject of group B %d: %s', sub, current_sub));
    
    current_sub_data=load_nii(current_sub);
    t_size=size(current_sub_data.img,4);
    disp(sprintf('Timepoints subject %d', t_size));
    
    for voxel=1:voxel_count
        temp_data_groupB(voxel,sub,:)=squeeze(current_sub_data.img(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),:));
    end
    
    clear current_sub_data
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's start with ISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s start computations...'));
disp(sprintf('######################################################'));

data_groupA=matrix2cell(temp_data_groupA);
data_groupB=matrix2cell(temp_data_groupB);


results_isc_groupAB=cell(voxel_count,1);
results_isc_groupAB_perm=cell(voxel_count,1);
results_isc_se_groupAB_perm=cell(voxel_count,1);

coordinates_mask_groupA=cell(voxel_count,1);
coordinates_mask_groupB=cell(voxel_count,1);

permutation_maskA=repmat([0:subjects_groupA-1],t_size,1);
permutation_maskA=permutation_maskA.*t_size;
permutation_maskB=repmat([0:subjects_groupB-1],t_size,1);
permutation_maskB=permutation_maskB.*t_size;


disp(sprintf('\nMeasure ISC between groups....'));
warning off
cpu_in_parallel=parpool('local',CPUs);


tic
parfor voxel=1:voxel_count
    warning off
    
    temp_dataA=data_groupA{voxel};
    distances=1-pdist(temp_dataA,'correlation');
    if sum(isnan(distances))<=numel(distances)/tollerance
        coordinates_mask_groupA{voxel}=1;
    else
        coordinates_mask_groupA{voxel}=0;
    end
    
    temp_dataB=data_groupB{voxel};
    distances=1-pdist(temp_dataB,'correlation');
    if sum(isnan(distances))<=numel(distances)/tollerance
        coordinates_mask_groupB{voxel}=1;
    else
        coordinates_mask_groupB{voxel}=0;
    end
    
    if(coordinates_mask_groupB{voxel}==1 & coordinates_mask_groupA{voxel}==1)
        distances=1-pdist2(temp_dataA,temp_dataB,'correlation');
        results_isc_groupAB{voxel}=distances;
        [distances_perm_mean,distances_perm_se]=ISC_twogroups_function(temp_dataA,randomA,permutation_maskA,temp_dataB,randomB,permutation_maskB);
        results_isc_groupAB_perm{voxel}=distances_perm_mean;
        results_isc_se_groupAB_perm{voxel}=distances_perm_se;
    end
    
end
toc

delete(cpu_in_parallel);
warning on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's do the statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s do statistics...'));
disp(sprintf('######################################################'));

voxel_count_final_groupA=sum(cell2mat(coordinates_mask_groupA));
disp(sprintf('\nNumber of voxels processed in group A %d', voxel_count_final_groupA));

voxel_count_final_groupB=sum(cell2mat(coordinates_mask_groupB));
disp(sprintf('Number of voxels processed in group B %d', voxel_count_final_groupB));

coordinates_mask=cell2mat(coordinates_mask_groupA).*cell2mat(coordinates_mask_groupB);
voxel_count_final=sum(coordinates_mask);
disp(sprintf('Final number of voxel between groups %d', voxel_count_final));


results_isc_matrix_groupAB=cell2matrix(results_isc_groupAB);
results_isc_matrix_groupAB_perm=cell2matrix(results_isc_groupAB_perm);
results_isc_se_matrix_groupAB_perm=cell2matrix(results_isc_se_groupAB_perm);


%%%%%%Prepare the space to save results
RESULTS_3D=zeros(x_size,y_size,z_size,4);

results_isc_matrix_groupAB_perm_good=nan(voxel_count_final,permutations);
pdistro_between=nan(voxel_count_final,1);
effectdistro_between=nan(voxel_count_final,1);

voxel_good=1;
for voxel=1:voxel_count
    
    if(coordinates_mask(voxel)==1)
        raw_between=nanmean(results_isc_matrix_groupAB(voxel,:));
        se_between=nanstd(results_isc_matrix_groupAB(voxel,:))./sqrt(numel(results_isc_matrix_groupAB(voxel,:)));
        effect_between=raw_between./se_between;
        
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),1)=raw_between;
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),2)=effect_between;
        
        %%%%%%Get raw p-values
        null_distro_raw=squeeze(results_isc_matrix_groupAB_perm(voxel,:));
        null_distro_raw_se=squeeze(results_isc_se_matrix_groupAB_perm(voxel,:));
        null_distro=null_distro_raw./null_distro_raw_se;
        
        null_distro_cat=cat(2,null_distro,effect_between);
        null_distro_cat_sort=sort(null_distro_cat,'descend');
        effect_pos=find(null_distro_cat_sort==effect_between);
        p_between=effect_pos(end)/(permutations+1);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),3)=1-p_between;
        
        effectdistro_between(voxel_good)=effect_between;
        pdistro_between(voxel_good)=p_between;
        results_isc_matrix_groupAB_perm_good(voxel_good,:)=null_distro;
        voxel_good=voxel_good+1;
    end
end


%%%%%%%FWE correction
max_between=nanmax(results_isc_matrix_groupAB_perm_good,[],1);

for voxel=1:voxel_count
    if(coordinates_mask(voxel)==1)
        raw_between=nanmean(results_isc_matrix_groupAB(voxel,:));
        se_between=nanstd(results_isc_matrix_groupAB(voxel,:))./sqrt(numel(results_isc_matrix_groupAB(voxel,:)));
        effect_between=raw_between./se_between;
        [pvalue,critical_value_at_p]=pareto_right_tail(max_between(:), effect_between, 0.05);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),4)=abs(log10(pvalue));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's save the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s save nifti & workspace...'));
disp(sprintf('######################################################'));

results_final=make_nii(RESULTS_3D, voxel_size, [0 0 0]);
save_nii(results_final, save_file);

disp(sprintf('If you are using AFNI/FSL, please fix the nifti header using AFNI'));
disp(sprintf('3drefit -newid -view tlrc -space MNI -duporigin %s %s',filenames_groupA{1}, save_file));
disp(sprintf('3drefit -relabel_all_str ''r_between t_between 1-p fwe_log_p'' %s', save_file));


%%%%Save the null distro for further analyses
save(save_file_mat,'coordinates','results_isc_matrix_groupAB','results_isc_matrix_groupAB_perm','results_isc_se_matrix_groupAB_perm','coordinates_mask','RESULTS_3D','max_between');

