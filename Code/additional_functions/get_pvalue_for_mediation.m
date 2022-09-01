%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results=get_pvalue_for_mediation(isc_raw,isc_model,isc_raw_null,isc_model_null,permutations)

dof=size(isc_raw,1);
voxel_count=size(isc_raw,2);

results.null_distro_model=nan(voxel_count,permutations);
results.null_distro_model_tstat=nan(voxel_count,permutations);
results.effect_model=nan(voxel_count,dof);
results.effect_model_tstat=nan(voxel_count,1);
results.effect_model_pvalue=nan(voxel_count,1);

for voxel=1:voxel_count
    results.effect_model(voxel,:)=isc_raw(:,voxel)-isc_model(:,voxel); %%%This is the model gain
    effect_model_mean=nanmean(results.effect_model(voxel,:));
    effect_model_se=nanstd(results.effect_model(voxel,:))./sqrt(dof-1);
    results.effect_model_tstat(voxel,:)=effect_model_mean/effect_model_se;
    
    null_distro_mean=nan(permutations,1);
    null_distro_std=nan(permutations,1);
    null_distro_tstat=nan(permutations,1);
    for p=1:permutations
        temp_null=isc_model_null{p};
        isc_raw_temp=isc_raw_null{p};
        null_distro_mean(p)=nanmean(isc_raw_temp(:,voxel)-temp_null(:,voxel),1);
        null_distro_std(p)=nanstd(isc_raw_temp(:,voxel)-temp_null(:,voxel))./sqrt(dof-1);
        null_distro_tstat(p)=null_distro_mean(p)./null_distro_std(p);
    end
    
    results.null_distro_model(voxel,:)=null_distro_mean;
    results.null_distro_model_tstat(voxel,:)=null_distro_tstat;
    
    null_distro_cat=cat(1,null_distro_tstat,results.effect_model_tstat(voxel,:));
    null_distro_cat_sort=sort(null_distro_cat,'descend');
    effetto_pos=find(null_distro_cat_sort==results.effect_model_tstat(voxel,:));
    results.effect_model_pvalue(voxel)=effetto_pos(end)/(permutations+1);
end

end
