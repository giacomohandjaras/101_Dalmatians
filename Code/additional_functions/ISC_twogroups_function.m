%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distances_perm_mean,distances_perm_se]=ISC_twogroups_function(temp_dataA,randomA,permutation_maskA,temp_dataB,randomB,permutation_maskB)

permutations=size(randomA.permutation_schema,1);
distances_perm_mean=nan(permutations,1);
distances_perm_se=nan(permutations,1);
subjects_groupA=size(temp_dataA,1);
subjects_groupB=size(temp_dataB,1);

temp_dataA_flip=temp_dataA';
temp_dataB_flip=temp_dataB';
se_num=sqrt(subjects_groupA*subjects_groupB);

for perm=1:permutations
    permutation_matrixA=squeeze(randomA.permutation_tps(perm,1:subjects_groupA,:))'+permutation_maskA;
    permutation_matrixB=squeeze(randomB.permutation_tps(perm,1:subjects_groupB,:))'+permutation_maskB;
    
    temp_dataA_perm=temp_dataA_flip(permutation_matrixA);
    temp_dataB_perm=temp_dataB_flip(permutation_matrixB);
    
    distances_perm=1-pdist2(temp_dataA_perm',temp_dataB_perm','correlation');
    distances_perm_mean(perm)=nanmean(distances_perm(:));
    distances_perm_se(perm)=nanstd(distances_perm(:))./se_num;
end


end
