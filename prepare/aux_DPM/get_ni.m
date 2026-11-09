function ni = get_ni(time_visit,numIDs)
    % This function generates a vector with the number of images for each
    % subject.
    % time_visit = time since baseline (bl =0)
    % numIDs = number of subjects
    % ni = number of images for each subject
    
    %% Calculating number of images
    baseline=find(time_visit==0);
    numScans=length(time_visit);
    ni=baseline([2:end,end])-baseline;
    ni(end)=numScans-baseline(end)+1;
    
    %% Checking errors
    if(length(ni)~=numIDs)
        warning('Error calculating ni: There are more scans in baseline than subjects');
    end
    if(sum(ni)~=numScans)
        warning('Error calculating ni');
    end    
end