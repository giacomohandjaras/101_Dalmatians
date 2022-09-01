%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

addpath('NIfTI_tools/');
addpath('additional_functions/');

filename_mask='../Results/mask_conjunction_TD.nii.gz';  %%%voxels to be considered

%%%%%Files to be saved
save_file='../Results/Results_mediation_SD_low.nii'; %the nifti
save_file_mat='../Results/Results_mediation_SD_low.mat'; %the workspace

groups='SD'; %%%groups in our experiment: 'TD' for Tipically developed, 'SD' for sensory deprived
models='low'; %%%models to be removed from fMRI data prior to ISC: 'low' for low-level auditory and visual, 'high-w2v' for high-level word2vec semantic, 'high-gpt3' for high-level GPT-3 semantic, 'editing' for the editing model (scene cuts etc...)

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


switch models
    case 'low'
        %%%This is not a mistake! Here we want to remove fine-grained correlations of the visual model with fMRI elicited during auditory stimulation, and viceversa
        model_matrix_groupB=csvread('../Computational_Models/Reduced_Models/lowlevel_visual.1D');
        model_matrix_groupA=csvread('../Computational_Models/Reduced_Models/lowlevel_acoustic.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_groupA,1),1),model_matrix_groupA);
        model_matrix_groupB=cat(2,ones(size(model_matrix_groupB,1),1),model_matrix_groupB);
        
    case 'high-w2v'
        model_matrix_group=csvread('../Computational_Models/Reduced_Models/highlevel_word2vec.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);
        model_matrix_groupB=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);

    case 'high-gpt3'
        model_matrix_group=csvread('../Computational_Models/Reduced_Models/highlevel_gpt3.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);
        model_matrix_groupB=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);
                
    case 'editing'
        model_matrix_group=csvread('../Computational_Models/Reduced_Models/editing.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);
        model_matrix_groupB=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);
end




timepoints_stimulus=1614;
randomA=load('permutation_schemas/permutation_schema_final_1000perm.mat');
randomB=load('permutation_schemas/permutation_schema_final_1000perm_alt.mat');

%%%%% number of CPU cores to be used %%%%%
CPUs=4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's open the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Loading models...'));
disp(sprintf('######################################################'));
disp(sprintf('Models loaded for %s condition', models));
disp(sprintf('Models related to group A contains %d columns', size(model_matrix_groupA,2)));
disp(sprintf('Models related to group B contains %d columns', size(model_matrix_groupB,2)));


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
%%%%%%%Let's clean the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\nLet''s clean the data form the models....'));
tic

data_groupA=permute(temp_data_groupA,[2,1,3]);
data_groupB=permute(temp_data_groupB,[2,1,3]);

data_groupA_cleaned=cleaning_data_with_model(data_groupA,model_matrix_groupA);
data_groupB_cleaned=cleaning_data_with_model(data_groupB,model_matrix_groupB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's start with ISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s start computations...'));
disp(sprintf('######################################################'));

%%%%%Let's define all the subjects' pairings
subjects_rdm=ones(subjects_groupA,subjects_groupB);
[subjects_pos(:,1),subjects_pos(:,2)]=ind2sub(size(subjects_rdm),find(subjects_rdm>0));
subjects_pairings=size(subjects_pos,1);


disp(sprintf('\nMeasure ISC between groups....'));

results_isc_raw=zeros(subjects_pairings,voxel_count);
results_isc_model=zeros(subjects_pairings,voxel_count);

tic
for pairing=1:subjects_pairings
    results_isc_raw(pairing,:)=fast_corr(squeeze(data_groupA(subjects_pos(pairing,1),:,:))',squeeze(data_groupB(subjects_pos(pairing,2),:,:))');
    results_isc_model(pairing,:)=fast_corr(squeeze(data_groupA_cleaned(subjects_pos(pairing,1),:,:))',squeeze(data_groupB_cleaned(subjects_pos(pairing,2),:,:))');
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's do the permutation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voxels_split=[1:250:voxel_count];  %%%%let's split voxels in small chuncks to reduce CPU cache/memory bottleneck
if voxels_split(end)<voxel_count
    voxels_split(end+1)=voxel_count+1;
end
voxels_split_num=numel(voxels_split);

disp(sprintf('Let''s do the permutation....'));

permutation_maskA=repmat([0:subjects_groupA-1],t_size,1);
permutation_maskA=permutation_maskA.*t_size;
permutation_maskB=repmat([0:subjects_groupB-1],t_size,1);
permutation_maskB=permutation_maskB.*t_size;

results_isc_raw_null=cell(permutations,voxels_split_num-1);
results_isc_model_null=cell(permutations,voxels_split_num-1);

warning off
cpu_in_parallel=parpool('local',CPUs);

p_display=1;
tic

for p=1:permutations
    
    if(p_display==10)
        disp(sprintf('Iteration: %d/%d....', p,permutations));
        p_display=0;
    end
    
    data_groupA_perm=nan(subjects_groupA,voxel_count,t_size);
    data_groupB_perm=nan(subjects_groupB,voxel_count,t_size);
    permutation_matrixA=squeeze(randomA.permutation_tps(p,1:subjects_groupA,:))'+permutation_maskA;
    permutation_matrixB=squeeze(randomB.permutation_tps(p,1:subjects_groupB,:))'+permutation_maskB;
    
    for voxel=1:voxel_count
        temp_groupA=squeeze(data_groupA(:,voxel,:));
        temp_groupA_flip=temp_groupA';
        temp_groupA_flip_perm=temp_groupA_flip(permutation_matrixA);
        data_groupA_perm(:,voxel,:)=temp_groupA_flip_perm';
        
        temp_groupB=squeeze(data_groupB(:,voxel,:));
        temp_groupB_flip=temp_groupB';
        temp_groupB_flip_perm=temp_groupB_flip(permutation_matrixB);
        data_groupB_perm(:,voxel,:)=temp_groupB_flip_perm';
    end
    
    clear temp_groupA temp_groupA_flip_perm temp_groupB temp_groupB_flip_perm
    
    
    parfor chunk=1:voxels_split_num-1
        warning off
        voxel_begin=voxels_split(chunk);
        voxel_end=voxels_split(chunk+1)-1;
        
        %%%cleaning shuffled data
        temp_data_groupA_cleaned=cleaning_data_with_model(data_groupA_perm(:,voxel_begin:voxel_end,:),model_matrix_groupA);
        temp_data_groupB_cleaned=cleaning_data_with_model(data_groupB_perm(:,voxel_begin:voxel_end,:),model_matrix_groupB);
        
        temp_results_isc_raw_null=nan(subjects_pairings,voxel_end-voxel_begin+1);
        temp_results_isc_model_null=nan(subjects_pairings,voxel_end-voxel_begin+1);
        
        %%%perform ISC with raw and cleaned data
        for pairing=1:subjects_pairings
            temp_results_isc_raw_null(pairing,:)=fast_corr(squeeze(data_groupA_perm(subjects_pos(pairing,1),voxel_begin:voxel_end,:))',squeeze(data_groupB_perm(subjects_pos(pairing,2),voxel_begin:voxel_end,:))');
            temp_results_isc_model_null(pairing,:)=fast_corr(squeeze(temp_data_groupA_cleaned(subjects_pos(pairing,1),:,:))',squeeze(temp_data_groupB_cleaned(subjects_pos(pairing,2),:,:))');
        end
        
        results_isc_raw_null{p,chunk}=temp_results_isc_raw_null;
        results_isc_model_null{p,chunk}=temp_results_isc_model_null;
        
    end %chunk end parfor
    
    p_display=p_display+1;
end % permutation

toc

delete(cpu_in_parallel);
warning on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%MERGE CELLS ACROSS VOXEL CHUNKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear data_*
clear temp_*
clear permutation_*
clear model_*

results_isc_raw_null_fixed=ISC_mediation_fix_cells(results_isc_raw_null,voxels_split);
results_isc_raw_null=results_isc_raw_null_fixed;
clear results_isc_raw_null_fixed;

results_isc_model_null_fixed=ISC_mediation_fix_cells(results_isc_model_null,voxels_split);
results_isc_model_null=results_isc_model_null_fixed;
clear results_isc_model_null_fixed;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's do the statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s do statistics...'));
disp(sprintf('######################################################'));

%%%%%%Get raw p-values
results_model=get_pvalue_for_mediation(results_isc_raw,results_isc_model,results_isc_raw_null,results_isc_model_null, permutations);


%%%%%%%FWE correction
max_distro=nanmax(results_model.null_distro_model_tstat,[],1);
results_model.effect_model_pvalue_fwe=nan(voxel_count,1);

for voxel=1:voxel_count
    [pvalue,critical_value_at_p]=pareto_right_tail(max_distro(:), results_model.effect_model_tstat(voxel), 0.05);
    results_model.effect_model_pvalue_fwe(voxel)=pvalue;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's save the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s save nifti & workspace...'));
disp(sprintf('######################################################'));

RESULTS_3D=zeros(x_size,y_size,z_size,5);

for voxel=1:voxel_count
    if (isnan(nanmean(results_isc_raw(:,voxel)))~=1)
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),1)=nanmean(results_isc_raw(:,voxel));
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),2)=nanmean(results_model.effect_model(voxel,:));
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),3)=results_model.effect_model_tstat(voxel);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),4)=abs(log10(results_model.effect_model_pvalue(voxel)));
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),5)=abs(log10(results_model.effect_model_pvalue_fwe(voxel)));
    end
end

results_final=make_nii(RESULTS_3D, voxel_size, [0 0 0]);
save_nii(results_final, save_file);
disp(sprintf('If you are using AFNI/FSL, please fix the nifti header using AFNI'));
disp(sprintf('3drefit -newid -view tlrc -space MNI -duporigin %s %s',filenames_groupA{1}, save_file));
disp(sprintf('3drefit -relabel_all_str ''ISC_raw ISC_gain_model ISC_gain_model_tstat ISC_gain_model_logp ISC_gain_model_fwe'' %s', save_file));

%%%%Save the null distro for further analyses
save(save_file_mat,'coordinates','filenames_groupA','filenames_groupB','permutations','results_isc_raw','results_model','models','groups','save_file','voxel_count','voxel_size','voxels_split','x_size','y_size','z_size');

