%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

seed=2001; %random seed for reproducibility
rand ('state', seed);

%%%%Definition of runs in TR
timeserie=cat(2,ones(1,89),ones(1,89)*2,ones(1,90)*3,ones(1,75)*4,ones(1,75)*5,ones(1,75)*6,ones(1,106)*7,ones(1,107)*8,ones(1,107)*9,ones(1,108)*10,ones(1,108)*11,ones(1,109)*12,ones(1,78)*13,ones(1,79)*14,ones(1,79)*15,ones(1,80)*16,ones(1,80)*17,ones(1,80)*18);

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


save('permutation_schema_final_1000perm.mat','permutation_tps','permutation_schema','-v7.3');


