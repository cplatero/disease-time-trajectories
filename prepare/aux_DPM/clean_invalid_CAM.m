function [new_CAM] = clean_invalid_CAM(CAM,ni,time_visit)
new_CAM = CAM;
numPatients=length(ni);
numImg=1;
%fprintf('\n');
for i=1:numPatients
    
    % if there is any element which is nan:
    tmp_CAM = CAM(numImg:numImg+ni(i)-1,:);
    invalid_visits=sum(isnan(tmp_CAM),2);
    if sum(invalid_visits)>0
        % Inside this if statement the affected visit must be modified
        % or deleted in order to avoid Nan values. If deleted, function
        % input must be changed in order to include all the features.
        
        %fprintf('subject %d had some nan values.\n',i);
        
        for f=1:size(tmp_CAM,2)
            % if the feature has a nan element, interpolation/extrapolation
            % is executed and saved into the output matrix.
            if sum(isnan(tmp_CAM(:,f)))>0
                
                % First the nan values are separated:
                times = time_visit(numImg:numImg+ni(i)-1);
                valid_idx = find(~isnan(tmp_CAM(:,f)));
                valid_visits = tmp_CAM(valid_idx,f);
                valid_times = times(valid_idx);                
%                try
%                 new_CAM(numImg:numImg+ni(i)-1,f) = interp1(valid_times,valid_visits,times,'linear','extrap');
%                catch
%                    new_CAM(:,f) = valid_visits;
%                end
                no_valid_idx = find(isnan(tmp_CAM(:,f)));
                no_valid_times = times(no_valid_idx);
                if(length(valid_idx)>1)
                    new_CAM(numImg+no_valid_idx-1,f)=interp1(valid_times,...
                        valid_visits,no_valid_times,'linear','extrap');
                else
                    if(isempty(valid_visits))
                        new_CAM(numImg+no_valid_idx-1,f)=mean(CAM(:,f),'omitnan');
                        %new_CAM(numImg+no_valid_idx-1,f)=nan;
                    else
                        new_CAM(numImg+no_valid_idx-1,f)= valid_visits;
                    end
                end
            end
        end
        
        
    end
    numImg = numImg+ni(i);
end

end
