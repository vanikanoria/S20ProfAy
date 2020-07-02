function [isTrue] = sustOscillations(mh1) %checks for sustained oscillations in eval mRNA
    peaks=0;
    troughs=0;
    index_peaks=zeros(2,1);
    index_troughs=zeros(2,1);
    value_peaks=zeros(2,1);
    value_troughs=zeros(2,1);
    
    mh1_cell_1=(mh1(1,:));
    timesteps_num = size(mh1(1,:));
    for z=2:timesteps_num-1
        if mh1_cell_1(z-1)<mh1_cell_1(z)&&mh1_cell_1(z)>mh1_cell_1(z+1)
            peaks=peaks+1;
            index_peaks(peaks)=z;
            value_peaks(peaks)=mh1_cell_1(z);
        end
        if mh1_cell_1(z-1)>mh1_cell_1(z)&&mh1_cell_1(z)<mh1_cell_1(z+1)
            troughs=troughs+1;
            index_troughs(troughs)=z;
            value_troughs(troughs)=mh1_cell_1(z);
        end
    end
    if peaks>=3 %Should fail if fewer than 3 peaks in 5 days.
        half_peaks=index_peaks(index_peaks<84);%MUST CHANGE IF CHANGE NUMBER OF HOURS. should be hours/2.
        mid=size(half_peaks,1);
        if mid == 0
            mid=1;
        end
        for k=1:length(value_peaks)
            if value_peaks(k)<0.1 %Should never have peak w/ amplitude less than 0.1 (indicates noise in the system).
                isTrue=false;
                return
            end
        end
        a=find(index_peaks>=126, 1);
        b=find(index_troughs>=126,1);
        if isempty(a) || isempty(b) %Should still have oscillations 3/4 of the way through the 168-hour run.
            isTrue=false;
            return
        end
        if peaks<troughs %Case 1: Fewer peaks than troughs
            peakToTrough=zeros(peaks,1);
            for z=1:peaks
                peakToTrough(z)=value_peaks(z)/((value_troughs(z)+value_troughs(z+1))/2);
            end
            ratioEnd = peakToTrough(peaks);%peak-to-trough ratio at end of run.
            if mid<=length(peakToTrough)
                ratioMid = peakToTrough(mid); %peak-to-trough ratio at middle of run.
            else
                ratioMid = peakToTrough(peaks);
            end
        elseif peaks>troughs %Case 2: More peaks than troughs.
            peakToTrough = zeros(troughs,1);
            for z=1:troughs
                peakToTrough(z)=((value_peaks(z)+value_peaks(z+1))/2)/value_troughs(z);   
            end
            ratioEnd = peakToTrough(troughs);
            if mid<=length(peakToTrough)
                ratioMid = peakToTrough(mid);
            else
                ratioMid = peakToTrough(troughs);
            end
        elseif peaks==troughs %Case 3: Same number of peaks as troughs.
            peakToTrough = zeros(peaks,1);
            for z=1:troughs
                peakToTrough(z)=value_peaks(z)/(value_troughs(z));   
            end
            ratioEnd = peakToTrough(troughs);
            if mid<=length(peakToTrough)
                ratioMid = peakToTrough(mid);
            else
                ratioMid = peakToTrough(troughs);
            end
        end
        amplitudeMid=(value_peaks(mid)-value_troughs(mid))/2;%Calculate amplitude at middle point.
        if ratioEnd >= 1.5 && ratioMid >= 1.5 && ratioMid/ratioEnd <= 1.02 && ratioMid/ratioEnd>=0.98 && amplitudeMid>=0.15 %Apply all conditions.
            isTrue=true;
            return;
        else
            isTrue=false;
            return;
        end
    else
        isTrue=false;
        return;
    end
end