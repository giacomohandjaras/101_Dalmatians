%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

seed=1969; %random seed for reproducibility. You need to change it if you want to perform a between groups analysis
rand ('state', seed);

%%%%Definition of runs in TR
timeserie=cat(2,ones(1,15),ones(1,15)*2,ones(1,15)*3,ones(1,15)*4,ones(1,15)*5,ones(1,15)*6,ones(1,15)*7,ones(1,15)*8,ones(1,15)*9,ones(1,15)*10,ones(1,15)*11,ones(1,15)*12,ones(1,15)*13,ones(1,15)*14,ones(1,15)*15,ones(1,15)*16,ones(1,15)*17,ones(1,10)*18);

subjects=10;
permutations=1000;

timeserie_tps=[1:numel(timeserie)];
chunks=max(timeserie);
chunks_sign=cat(2,ones(1,chunks/2)*1,ones(1,chunks/2)*-1);

permutation_schema=nan(permutations,subjects,chunks);
permutation_tps=nan(permutations,subjects,numel(timeserie));


for perm=1:permutations 
    for sub=1:subjects
        perm_seed=randperm(chunks);
        perm_sign=chunks_sign(randperm(chunks));
        permutation_schema(perm,sub,:)=perm_seed.*perm_sign;
        
        timeserie_perm=[];
        for c=1:chunks
            tps=find(timeserie==perm_seed(c));
            if(perm_sign(c)<0)
                tps=flip(tps);
            end
            timeserie_perm=cat(2,timeserie_perm,timeserie_tps(tps));
        end
        
        permutation_tps(perm,sub,:)=timeserie_perm(:);
    end
end


save('permutation_schema_ctrl_final_1000perm.mat','permutation_tps','permutation_schema','-v7.3');


