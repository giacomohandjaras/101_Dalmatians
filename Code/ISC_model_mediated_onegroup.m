%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

addpath('NIfTI_tools/');
addpath('additional_functions/');

filename_mask='../Results/mask_conjunction_TD.nii.gz';  %%%voxels to be considered

%%%%%Files to be saved
save_file='../Results/Results_mediation_AV_low-auditory.nii'; %the nifti
save_file_mat='../Results/Results_mediation_AV_low-auditory.mat'; %the workspace

%%%%%fMRI files to open
filenames_groupA={'../fMRI_101_Dalmatians/ctrlAV/sub-012.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-013.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-014.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-015.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-016.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-017.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-018.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-019.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-022.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-032.nii.gz'};

models='low-auditory'; %%%models to be removed from fMRI data prior to ISC: 'low-visual' for low-level visual, 'low-auditory' for low-level auditory, 'high-w2v' for high-level semantic largely related to word2vec, 'high-gpt3' for high-level semantic largely related to gpt3, 'editing' for the editing model (scene cuts etc...)

switch models
    case 'low-visual'
        model_matrix_groupA=csvread('../Computational_Models/Reduced_Models/lowlevel_visual.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_groupA,1),1),model_matrix_groupA);
        
    case 'low-auditory'
        model_matrix_groupA=csvread('../Computational_Models/Reduced_Models/lowlevel_acoustic.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_groupA,1),1),model_matrix_groupA);
        
    case 'high-w2v'
        model_matrix_group=csvread('../Computational_Models/Reduced_Models/highlevel_word2vec.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);

    case 'high-gpt3'
        model_matrix_group=csvread('../Computational_Models/Reduced_Models/highlevel_gpt3.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_group,1),1),model_matrix_group);
        
    case 'editing'
        model_matrix_groupA=csvread('../Computational_Models/Reduced_Models/editing.1D');
        model_matrix_groupA=cat(2,ones(size(model_matrix_groupA,1),1),model_matrix_groupA);
end



timepoints_stimulus=1614;
randomA=load('permutation_schemas/permutation_schema_final_1000perm.mat');

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


%%%%%Let's open all the subjects
disp(sprintf('\n######################################################'));
disp(sprintf('Let''s perform the analysis...'));
disp(sprintf('######################################################'));

temp_data_groupA=nan(voxel_count,subjects_groupA,timepoints_stimulus);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's clean the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\nLet''s clean the data form the models....'));
tic

data_groupA=permute(temp_data_groupA,[2,1,3]);
data_groupA_cleaned=cleaning_data_with_model(data_groupA,model_matrix_groupA);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Let's start with ISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s start computations...'));
disp(sprintf('######################################################'));

%%%%%Let's define all the subjects' pairings
subjects_rdm=ones(subjects_groupA);
subjects_rdm(triu(subjects_rdm,1)==0)=0;
[subjects_pos(:,1),subjects_pos(:,2)]=ind2sub(size(subjects_rdm),find(subjects_rdm>0));
subjects_pairings=size(subjects_pos,1);


disp(sprintf('\nMeasure ISC....'));

results_isc_raw=zeros(subjects_pairings,voxel_count);
results_isc_model=zeros(subjects_pairings,voxel_count);

tic
for pairing=1:subjects_pairings
    results_isc_raw(pairing,:)=fast_corr(squeeze(data_groupA(subjects_pos(pairing,1),:,:))',squeeze(data_groupA(subjects_pos(pairing,2),:,:))');
    results_isc_model(pairing,:)=fast_corr(squeeze(data_groupA_cleaned(subjects_pos(pairing,1),:,:))',squeeze(data_groupA_cleaned(subjects_pos(pairing,2),:,:))');
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
    permutation_matrixA=squeeze(randomA.permutation_tps(p,1:subjects_groupA,:))'+permutation_maskA;
    
    for voxel=1:voxel_count
        temp_groupA=squeeze(data_groupA(:,voxel,:));
        temp_groupA_flip=temp_groupA';
        temp_groupA_flip_perm=temp_groupA_flip(permutation_matrixA);
        data_groupA_perm(:,voxel,:)=temp_groupA_flip_perm';
    end
    
    clear temp_groupA temp_groupA_flip_perm
    
    
    parfor chunk=1:voxels_split_num-1
        warning off
        voxel_begin=voxels_split(chunk);
        voxel_end=voxels_split(chunk+1)-1;
        
        %%%cleaning shuffled data
        temp_data_groupA_cleaned=cleaning_data_with_model(data_groupA_perm(:,voxel_begin:voxel_end,:),model_matrix_groupA);
        
        temp_results_isc_raw_null=nan(subjects_pairings,voxel_end-voxel_begin+1);
        temp_results_isc_model_null=nan(subjects_pairings,voxel_end-voxel_begin+1);
        
        %%%perform ISC with raw and cleaned data
        for pairing=1:subjects_pairings
            temp_results_isc_raw_null(pairing,:)=fast_corr(squeeze(data_groupA_perm(subjects_pos(pairing,1),voxel_begin:voxel_end,:))',squeeze(data_groupA_perm(subjects_pos(pairing,2),voxel_begin:voxel_end,:))');
            temp_results_isc_model_null(pairing,:)=fast_corr(squeeze(temp_data_groupA_cleaned(subjects_pos(pairing,1),:,:))',squeeze(temp_data_groupA_cleaned(subjects_pos(pairing,2),:,:))');
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
save(save_file_mat,'coordinates','filenames_groupA','permutations','results_isc_raw','results_model','models','save_file','voxel_count','voxel_size','voxels_split','x_size','y_size','z_size');

