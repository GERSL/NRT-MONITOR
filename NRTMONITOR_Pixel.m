function rec_cg = NRTMONITOR_Pixel(rec_pre, clrx, clry, pid, adj_rmse, T_cg, t_span, preLength,iFFactor)
%% 
%     This function use a recursive Forgetting Factor approach to detect
% changes continuiously. 




    % Defining variables (constants)
    num_fc = 1;        % initialize NUM of Functional Curves
    rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'pos',[],...
        'rmse',[], 'rmse_cnt',[],'num_obs',[],'change_prob',[],...
        'category',[],'magnitude',[],'initial',[],'obs_last',[]);
    
    rpre_initial = rec_pre(end).initial;
    if rpre_initial == 1   % Last segment is fitted by a model (with parameters).       
        %% Three scenarios with model parameters:
        %  S1-2 --- No saved obs_last (Paper Fig. 6a-b)
        %  S4 --- Saved 12 latest obs: usually change probability < 100% (Paper Fig. 6d)
        
        % Extract parameters and RMSE LUT
        rpre_coefs = rec_pre(end).coefs(:,1:6);
        rpre_rmse = rec_pre(end).rmse(:,1:6);
        rpre_rmse_cnt = rec_pre(end).rmse_cnt;
        
        % Merge historical and new observations
        rpre_cp = rec_pre(end).change_prob;
        if rpre_cp > 0 && rpre_cp < 1 % S4
            rpre_obs_last = rec_pre(end).obs_last;
            clrx = [rpre_obs_last(:,7);clrx];
            clry = [rpre_obs_last(:,1:6);clry];
        end           
    else % Last segment isn't fitted by a model (without parameters).        
        %% One scenarios without model parameters:
        %  S3 --- Saved obs of last segment: change confirmed but not initial (Paper Fig. 6c)
        
        % Merge historical and new observations
        rpre_obs_last = rec_pre(end).obs_last;
        if ~isempty(rpre_obs_last)
            clrx = [rpre_obs_last(:,7);clrx];
            clry = [rpre_obs_last(:,1:6);clry];
        end
        
        % Set ZERO parameters and RMSE LUT
        rpre_coefs = zeros(8,6);
        rpre_rmse = zeros(46,6);
        rpre_rmse_cnt = zeros(46,1);    
    end
    
    % Remove duplicated observations (L8 and S2 both have observations for that date)
    [clrx,uniq_id,~] = unique(clrx);
    tmp_y = zeros(length(clrx),6);
    for i = 1:6
        tmp_y(:,i) = clry(uniq_id,i);
    end
    clry = tmp_y;
    
    % -------------------------- Start monitoring --------------------------
    iTotal = length(clrx);    % record the number of observations
    if iTotal < 12 && rpre_initial == 0
        %% Case0: save the record with very few obs.
        % Only one time segment in the entire time series, and it has
        % very few observations. Directly save the record, and a time 
        % model will be fitted if there are at least 6 observations.
        iSaveObs = 1;
        rec_temp = SaveRecord_TotalFewObs(pid, clrx, clry,iSaveObs);
        rec_cg(num_fc) = rec_temp(1);
    else
        i_start = 1;              % the first observation for TSFit 
        pos_PreBreak = 1;
        iLastSeg = 0;
        % The while loope starts if the remaining time segment has at 
        % least 12 observations, or the time series has initial model 
        % parameters.                        
        while i_start <= iTotal-11 | rpre_initial > 0
            % Function: Change detection.
            % If a break is detected, the remianing time seires will be
            % used to repeat the change detection. Note that even 
            % though i_start is smaller than iTotal, the last obs will 
            % still be used in the change detection.
            if rpre_initial == 1
                [f_break, pos_break, rec_temp] = MOLDmonitoring_initial(pid,...
                    clrx(i_start:end), clry(i_start:end,:),T_cg, t_span, iFFactor,preLength,...
                    adj_rmse,rpre_coefs,rpre_rmse,rpre_rmse_cnt);                
            else
                [f_break, pos_break, rec_temp] = MOLDmonitoring(pid,...
                    clrx(i_start:end), clry(i_start:end,:),T_cg, t_span, iFFactor,...
                    pos_PreBreak,preLength,adj_rmse);
            end
                       
            if f_break == 0 
                %% Case1: save the record with no change.
                % The first time segment is the last.
                rec_cg(num_fc) = rec_temp(1);
                iLastSeg = 1;
                break;      
            else					
                % Confirm a change
                rec_cg(num_fc) = rec_temp(1);
                num_fc = num_fc + 1; 
                i_start = pos_break + i_start - 1;
                pos_PreBreak = pos_break;	
                rpre_initial = 0;
            end % if		                   
        end  % while
        
        if i_start > iTotal-11 && rpre_initial == 0 && iLastSeg == 0
            %% Case3: save obs in last time segment.
            % There are at least 6  but less than 12 observations in last time 
            % segment.
            iSaveObs = 1;
            rec_temp = SaveRecord_TotalFewObs(pid, clrx(i_start:end), ...
                clry(i_start:end,:), iSaveObs);
            rec_cg(num_fc) = rec_temp(1);
        end % Case3
        
    end % else
    
% function end
end







