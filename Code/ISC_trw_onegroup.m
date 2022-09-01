%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

addpath('NIfTI_tools/');
addpath('additional_functions/');

filename_mask='../Results/mask_conjunction_TD.nii.gz';

%%%%%Files to be saved
save_file='../Results/crap_Results_trw_onegroups_AV.nii'; %the nifti
save_file_mat='../Results/crap_Results_trw_onegroups_AV.mat'; %the workspace

%%%%%fMRI files to open
filenames_groupA={'../fMRI_101_Dalmatians/ctrlAV/sub-012.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-013.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-014.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-015.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-016.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-017.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-018.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-019.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-022.nii.gz', '../fMRI_101_Dalmatians/ctrlAV/sub-032.nii.gz'};

timepoints_stimulus=1614;
randomA=load('permutation_schemas/permutation_schema_final_1000perm.mat');

%%%%% number of CPU cores to be used %%%%%
CPUs=4;

%%%%% The maximum number of NaN across subject's pairings allowed in a voxels is no more than 1/tollerance. NaN depends on voxels with a timeseries of zeros.%%%%%
tollerance=3;

%%%%% The definition of the temporal windows, in timepoints%%%%%
windows_sizes=[1:1:120];


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
%%%%%%%Let's start with ISC/TRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n######################################################'));
disp(sprintf('Let''s start computations...'));
disp(sprintf('######################################################'));

data_groupA=matrix2cell(temp_data_groupA);

results_isc_groupA=cell(voxel_count,1); %%%it will retain the best temporal window
results_isc_groupA_corr=cell(voxel_count,1); %%%it will contain the ISC for the best temporal window
results_isc_groupA_raw=cell(voxel_count,1); %%%it will contain the ISC for all the temporal windows
results_isc_groupA_tscore=cell(voxel_count,1); %%%it will contain the t-score for the best temporal window
results_isc_groupA_tscore_perm=cell(voxel_count,1); %%%it will contain the t-score  null distribution for the best temporal window
results_isc_groupA_perm=cell(voxel_count,1); %%%it will contain the null ISC for the best temporal window

coordinates_mask_groupA=cell(voxel_count,1);

permutation_maskA=repmat([0:subjects_groupA-1],t_size,1);
permutation_maskA=permutation_maskA.*t_size;


disp(sprintf('\nMeasure ISC TRWs within group....'));
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
    
    if(coordinates_mask_groupA{voxel}==1)
        
        data_isc_temp=nan(numel(windows_sizes),subjects_groupA*(subjects_groupA-1)/2);
        for w=1:numel(windows_sizes)
            temp_dataA_res=downsampling_signal(temp_dataA,windows_sizes(w));
            distances=1-pdist(temp_dataA_res,'correlation');
            data_isc_temp(w,:)=distances(:);
        end
        
        temporal_tscore=nanmean(data_isc_temp,2)./(nanstd(data_isc_temp,[],2)./sqrt(subjects_groupA*(subjects_groupA-1)/2));
        temporal_corr=nanmean(data_isc_temp,2);
        window_best=find(temporal_corr==nanmax(temporal_corr));  %%%%let's find the best window
        
        results_isc_groupA{voxel}=windows_sizes(window_best);
        results_isc_groupA_corr{voxel}=temporal_corr(window_best);
        results_isc_groupA_tscore{voxel}=temporal_tscore(window_best);
        results_isc_groupA_raw{voxel}=data_isc_temp;
        
        %%%%%Let's do the permutation on the best window
        distances_perm_tscore=nan(permutations,1);
        distances_perm=nan(permutations,1);
        
        for perm=1:permutations
            permutation_matrixA=squeeze(randomA.permutation_tps(perm,:,:));
            
            temp_dataA_perm=nan(size(temp_dataA));
            for sub=1:subjects_groupA
                temp_dataA_perm(sub,:)=temp_dataA(sub,squeeze(permutation_matrixA(sub,:)));
            end
            
            temp_dataA_perm_res=downsampling_signal(temp_dataA_perm,windows_sizes(window_best));
            distances=1-pdist(temp_dataA_perm_res,'correlation');
            distances_perm=nanmean(distances(:));
            distances_perm_tscore(perm)=nanmean(distances(:))./(nanstd(distances(:))./sqrt(subjects_groupA*(subjects_groupA-1)/2));
        end
        
        results_isc_groupA_tscore_perm{voxel}=distances_perm_tscore;
        results_isc_groupA_perm{voxel}=distances_perm;
        
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

coordinates_mask=cell2mat(coordinates_mask_groupA);
voxel_count_final=sum(coordinates_mask);
disp(sprintf('Final number of voxel %d', voxel_count_final));

results_isc_matrix_groupA=cell2matrix(results_isc_groupA);
results_isc_matrix_groupA_corr=cell2matrix(results_isc_groupA_corr);
results_isc_matrix_groupA_tscore=cell2matrix(results_isc_groupA_tscore);
results_isc_matrix_groupA_tscore_perm=cell2matrix(results_isc_groupA_tscore_perm);
results_isc_matrix_groupA_perm=cell2matrix(results_isc_groupA_perm);

%%%%%%Prepare the space to save results
results_isc_matrix_groupA_perm_good=nan(voxel_count_final,permutations);
RESULTS_3D=zeros(x_size,y_size,z_size,5);

voxel_good=1;

for voxel=1:voxel_count
    if(coordinates_mask(voxel)==1)
        effect=results_isc_matrix_groupA(voxel,:);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),1)=effect;
        effect_corr=results_isc_matrix_groupA_corr(voxel,:);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),2)=effect_corr;
        effect_tscore=results_isc_matrix_groupA_tscore(voxel,:);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),3)=effect_tscore;
        
        %%%%%%Get raw p-values
        distro_nulla=results_isc_matrix_groupA_tscore_perm(voxel,:);
        distro_nulla_cat=cat(2,distro_nulla,effect_tscore);
        distro_nulla_cat_sort=sort(distro_nulla_cat,'descend');
        effect_pos=find(distro_nulla_cat_sort==effect_tscore);
        p_tscore=effect_pos(end)/(permutations+1);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),4)=1-p_tscore;
        
        results_isc_matrix_groupA_perm_good(voxel_good,:)=distro_nulla;
        voxel_good=voxel_good+1;
    end
end


%%%%%%%FWE correction
max_tscore=nanmax(results_isc_matrix_groupA_perm_good,[],1);

for voxel=1:voxel_count
    if(coordinates_mask(voxel)==1)
        effect_tscore=results_isc_matrix_groupA_tscore(voxel,:);
        [pvalue,critical_value_at_p]=pareto_right_tail(max_tscore(:), effect_tscore, 0.05);
        RESULTS_3D(coordinates(voxel,1),coordinates(voxel,2),coordinates(voxel,3),5)=abs(log10(pvalue));
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
disp(sprintf('3drefit -relabel_all_str ''delay_TRs delay_r delay_tscore 1-p fwe_log_p'' %s', save_file));


%%%%Save the null distro for further analyses
save(save_file_mat,'coordinates','results_isc_matrix_groupA','results_isc_matrix_groupA_corr','results_isc_groupA_raw','results_isc_matrix_groupA_tscore','results_isc_matrix_groupA_tscore_perm','results_isc_matrix_groupA_perm','max_tscore','coordinates_mask','RESULTS_3D');

