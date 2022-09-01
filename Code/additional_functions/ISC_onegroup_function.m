%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distances_perm_mean,distances_perm_se]=ISC_onegroup_function(temp_dataA,randomA,permutation_maskA)

permutations=size(randomA.permutation_schema,1);
distances_perm_mean=nan(permutations,1);
distances_perm_se=nan(permutations,1);
subjects_groupA=size(temp_dataA,1);

temp_dataA_flip=temp_dataA';
se_num=sqrt(subjects_groupA*(subjects_groupA-1)/2);

for perm=1:permutations
    
    permutation_matrixA=squeeze(randomA.permutation_tps(perm,1:subjects_groupA,:))'+permutation_maskA;
    
    temp_dataA_perm=temp_dataA_flip(permutation_matrixA);
    
    distances_perm=1-pdist(temp_dataA_perm','correlation');
    
    distances_perm_mean(perm)=nanmean(distances_perm(:));
    distances_perm_se(perm)=nanstd(distances_perm(:))./se_num;
end


end
