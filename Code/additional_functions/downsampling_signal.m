%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal_down=downsampling_signal(signal,window)


if window>1
    subs=size(signal,1);
    tps=size(signal,2);
    %step=ceil(window/3);
    step=1;
    
    bins=[1:step:tps];
    bins_mask=bins<(tps-(window-step));
    bins(bins_mask==0)=[];
    
    signal_down=nan(subs,numel(bins)-1);
    
    for bin=1:(numel(bins)-1)
        temp_signal=signal(:,bins(bin):(bins(bin)+window-1));
        temp_signal_avg=nanmean(temp_signal,2);
        signal_down(:,bin)=temp_signal_avg;
    end
    
else
    signal_down=signal;
end

end
