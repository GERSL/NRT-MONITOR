function rec_cg = NRTMOLD_Pixel(rec_pre, clrx, clry, pid, adj_rmse, T_cg, t_span, preLength,iFFactor)
%% 
%     This function use a recursive Forgetting Factor approach to detect
% changes continuiously. Developed by Rong Shang, 2020/08/25.
%
% Function Inputs:
%    rec_pre --- Previous detected record;
%    clrx --- Date of obs;
%    clry --- Reflectance data;
%    pid --- Index for a pixel;
%    adj_rmse --- minimum RMSE;
%    T_cg --- Change detection threshold based on chi-squared distribution;
%    t_span --- Change detection minimum time interval;
%    preLength --- Length of created points to stabilize recursive steps;
%    iFFactor --- Forgetting parameter.
%
% Fuction output:
%    rec_cg --- The ouput structure.



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


function [f_break, pos_break, rec_cg] = MOLDmonitoring_initial(pid, clrx, clry,...
    T_cg, t_span, iFFactor,preLength,adj_rmse,rpre_coefs,rpre_rmse,rpre_rmse_cn)
%% 
%     This function conducts MOLD monitoring with initial parameters to 
% detect a change. When a change is detected, this function will return the 
% position of change in the time series. Noticed that this function only
% detect one change. Developed by Rong Shang, 2020/08/25.
%       
% Function Inputs:
%    pid --- Index for a pixel;
%    clrx --- Date of input reflectance;
%    clry --- Reflectance values for the detected bands;
%    T_cg --- Change detection threshold based on chi-squared distribution;
%    t_span --- Change detection minimum time interval;
%    iFFactor --- Forgetting parameter; 
%    preLength --- Length of created points to stabilize recursive steps;
%    adj_rmse --- minimum RMSE;
%    rpre_coefs --- Initial parameters;
%    rpre_rmse --- Sum of difference square for each LUT bin;
%    rpre_rmse_cn --- Count for each LUT bin. 
%
% Fuction output (rec_cg):
%    f_break --- Whether a change has been detected (1 is for change);
%    pos_break --- Posistion of a change; 
%    rec_cg --- Output struct of a change record.

    % ----------------- Fuction Start ----------------%
    p_start = 1;    % the first observation for TSFit.  
	num_fc =1;
    f_break =0;    % record whether a change has been detected.
    pos_break = p_start;    % record the posistion of change.
    min_num_c = 4;    % Number of coefficient required.
    mid_num_c = 6;
    max_num_c = 8;
    n_times = 3;  
    n_Model = 8;
    B_detect = 2:6;
    nbands = 7;
    conse = 4;
	def_pT_cg = T_cg;
	def_conse = conse;
    fit_cft = rpre_coefs;
    LUT_sum = rpre_rmse;
    LUT_cnt = rpre_rmse_cn;    
    rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'pos',[],...
        'rmse',[], 'rmse_cnt',[],'num_obs',[],'change_prob',[],...
        'category',[],'magnitude',[],'initial',[],'obs_last',[]);


    %% Recursive estimation of Forgetting Factor.
    clrx_ff = clrx;
    clry_ff = clry;
    cur_mt = datevec(clrx_ff);
    clrx_doy = clrx_ff - datenum(cur_mt(:,1),1,0);  
    clrx_bin = floor(clrx_doy/8)+1;
    nLength = length(clrx_ff);
    theta0 = zeros(n_Model,nbands-1);    % variable for initial step
    rmse_FF = zeros(1, nbands-1);    % predicted rmse  
    yhat = zeros(nLength, nbands-1);    % estimates of FF                       
    v_dif_FF = zeros(nLength,nbands-1);    % diference of estimates  
    for i_B = B_detect
        % Get all estimates from forgetting factor
        theta0(:,i_B) = fit_cft(1:n_Model,i_B);
        [yhat(:,i_B), v_dif_FF(:,i_B), rmse_FF(1,i_B)] = ...
            autoFFsrFit(clrx_ff, clry_ff(:,i_B),theta0(:,i_B),...
            n_Model, iFFactor, preLength);
    end
        

    %% Detect a change using the while loope
    i = p_start;    
    mini_rmse = adj_rmse;
    adj_conse = def_conse;
    while i <= nLength - def_conse                
        % Adjusted conse according to change length           
        currdate = clrx_ff(i)+t_span;
        if currdate > clrx_ff(nLength)
            break;    % Date exceedes the time seires!!
        end       
        index_delta = find(clrx_ff(i:nLength) >= currdate);
        if length(index_delta)<1
            break;    % No points!!
        end
        adj_conse = index_delta(1);   % First point > minimum time interval.    
        if adj_conse < def_conse
            adj_conse = def_conse;
        end          
        if i+ adj_conse - 1> nLength
            break;    % Last point exceedes the time seires!!
        end
            
        % Adjust T_cg based on conse and adj_conse
        T_cg = chi2inv(def_pT_cg,length(B_detect));             
        if adj_conse > def_conse
            pT_cg = 1-(1-def_pT_cg)^(def_conse/adj_conse);
            T_cg = chi2inv(pT_cg,length(B_detect));
        end

        % Calculate RMSE from the LUT
        mini_rmse = ExtractRMSE_LUTbin(clrx_bin(i),LUT_sum,LUT_cnt);

        % Calculate change metrix for adj_conse observations
        vec_mag = zeros(adj_conse,1);
        v_dif = zeros(adj_conse,nbands-1);  
        for i_conse = 1:adj_conse           
            for i_B = B_detect
                if sum(i_B == B_detect)   
                    v_dif(i_conse,i_B) = abs(v_dif_FF(i+i_conse-1,i_B))/max(adj_rmse(i_B),mini_rmse(i_B));
                end  % if i_B
            end  % for i_B
            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
        end  % for i_conse

        % Detect a change
        if min(vec_mag) > T_cg  
            % Updating output information
            rec_cg(num_fc).t_break = clrx_ff(i);    
            rec_cg(num_fc).t_start = clrx_ff(p_start);                                                                       
            rec_cg(num_fc).pos = pid;                
            rec_cg(num_fc).change_prob = 1;      
            rec_cg(num_fc).initial = 1; 
            rec_cg(num_fc).rmse = LUT_sum;   
            rec_cg(num_fc).rmse_cnt = LUT_cnt;
            rec_cg(num_fc).obs_last = 0;
            
            if i-1<p_start
                rec_cg(num_fc).t_end = 0;  
                rec_cg(num_fc).num_obs = 0;    
                rec_cg(num_fc).coefs = rpre_coefs;   
                rec_cg(num_fc).category = 8;
            else
                rec_cg(num_fc).t_end = clrx_ff(i-1);   
                rec_cg(num_fc).num_obs = length(p_start:i-1);  
                            
                % Determine the time series model
                IDs = p_start:i-1; 
                if length(IDs)>=6
                    i_span = i-1-p_start+1;              
                    update_num_c = update_cft(i_span,n_times,min_num_c,...
                        mid_num_c,max_num_c,max_num_c);
                    outfit_cft = zeros(8,nbands-1); 
                    outrmse = zeros(1, nbands-1); 
                    for i_B = 1:6
                        [outfit_cft(:,i_B),outrmse(1,i_B),~] = autoTSFit(...
                            clrx_ff(IDs), clry_ff(IDs, i_B), update_num_c);
                    end         
                    rec_cg(num_fc).coefs = outfit_cft;   
                    rec_cg(num_fc).category =  update_num_c;
                else
                    rec_cg(num_fc).coefs = rpre_coefs;   
                    rec_cg(num_fc).category = 8;
                end
            end

            % chagne vector magnitude
            v_dif_mag = zeros(conse,nbands-1);
            for i_conse = 1:conse
                for i_B = 1:nbands-1
                    % absolute difference
                    v_dif_mag(i_conse,i_B) = clry_ff(i+i_conse,i_B)-...
                        autoTSPred(clrx_ff(i+i_conse), rec_cg(num_fc).coefs(:,i_B));
                end
            end
            rec_cg(num_fc).magnitude = median(v_dif_mag,1);  

            f_break = 1;
            pos_break = i+p_start-1;
            return;  % Return the change record!!! 
        end  % if T_cg 
            
        % Update the LUT
        curM = clrx_bin(i);
        LUT_cnt(curM) = LUT_cnt(curM) + 1;     
        for i_B = B_detect
            LUT_sum(curM,i_B) = LUT_sum(curM,i_B) + (v_dif_FF(i, i_B))*(v_dif_FF(i, i_B));                    
        end
        i = i+1;
    end  % end WHILE


    %% Process the end
    id_last = nLength;
    for i_conse = nLength:-1:i    
        for i_B = B_detect
            v_dif(i_conse,i_B) = abs(v_dif_FF(i_conse,i_B))/max(adj_rmse(i_B),mini_rmse(i_B));
        end
        vec_mag = norm(v_dif(i_conse,B_detect))^2; 
        if vec_mag < T_cg
            % the last stable id
            id_last = i_conse;
            break;
        end
    end
    
    % No change is detected for the last time segment.
    if f_break == 0  
        % Set output information
        IDs = p_start:id_last; 
        rec_cg(num_fc).t_break = 0;    
        rec_cg(num_fc).t_start = clrx_ff(p_start);             
        rec_cg(num_fc).t_end = clrx_ff(id_last);                                                  
        rec_cg(num_fc).pos = pid;   
        rec_cg(num_fc).num_obs = length(IDs);      
        rec_cg(num_fc).initial = 1; 
        rec_cg(num_fc).rmse = LUT_sum;   
        rec_cg(num_fc).rmse_cnt = LUT_cnt;
        
        % Determine the time series model       
        if length(IDs)>=6
            i_span = id_last-p_start+1;              
            update_num_c = update_cft(i_span,n_times,min_num_c,...
                mid_num_c,max_num_c,max_num_c);
            outfit_cft = zeros(8,nbands-1); 
            outrmse = zeros(1, nbands-1); 
            for i_B = 1:6
                [outfit_cft(:,i_B),outrmse(1,i_B),~] = autoTSFit(...
                    clrx_ff(IDs), clry_ff(IDs, i_B), update_num_c);
            end         
            rec_cg(num_fc).coefs = outfit_cft;   
            rec_cg(num_fc).category =   update_num_c;
        else
            rec_cg(num_fc).coefs = rpre_coefs; 
            rec_cg(num_fc).category = 8;
        end
        rec_cg(num_fc).magnitude = 0;  
        rec_cg(num_fc).change_prob= 0;
        rec_cg(num_fc).obs_last = 0;

        if id_last < nLength    % potential change at the end
            rec_cg(num_fc).t_break = clrx_ff(id_last+1); 
            rec_cg(num_fc).initial = 0; 
            
            cp_conse = (nLength-id_last+1)/adj_conse;
            cp_tspan = (clrx(end)-clrx(id_last))/t_span;
            rec_cg(num_fc).change_prob=min(cp_conse,cp_tspan);

            if nLength>=12
                obs_last = zeros(12,nbands-1);
                obs_last(:,1:6) = clry_ff(end-11:end,1:6);
                obs_last(:,7) = clrx_ff(end-11:end);
            else
                obs_last = zeros(nLength,nbands-1);
                obs_last(:,1:6) = clry_ff(:,1:6);
                obs_last(:,7) = clrx_ff;
            end
            rec_cg(num_fc).obs_last = obs_last;
        end % if id_last
        return;
    end  % if f_break

 % function end     
end       

function [f_break, pos_break, rec_cg] = MOLDmonitoring(pid,clrx, clry,...
     T_cg, t_span, iFFactor,pos_PreBreak,preLength,adj_rmse)
%% 
%     This function conducts MOLD monitoring without initial parameters to 
% detect a change. When a change is detected, this function will return the 
% position of change in the time series. Noticed that this function only
% detect one change. Developed by Rong Shang, 2020/08/25.
%       
% Function Inputs:
%    pid --- Index for a pixel;
%    clrx --- Date of input reflectance;
%    clry --- Reflectance values for the detected bands;
%    T_cg --- Change detection threshold based on chi-squared distribution;
%    t_span --- Change detection minimum time interval;
%    iFFactor --- Forgetting parameter; 
%    pos_PreBreak --- Index for previous change;
%    preLength --- Length of created points to stabilize recursive steps;
%    adj_rmse --- minimum RMSE.
%
% Fuction output (rec_cg):
%    f_break --- Whether a change has been detected (1 is for change);
%    pos_break --- Posistion of a change; 
%    rec_cg --- Output struct of a change record.

    % ----------------- Fuction Start ----------------%
    p_start = 1;    % the first observation for TSFit.  
	num_fc =1;
    f_break =0;    % record whether a change has been detected.
    pos_break = p_start;    % record the posistion of change.
    min_num_c = 4;    % Number of coefficient required.
    mid_num_c = 6;
    max_num_c = 8;
    n_times = 3;  
    n_Model = 8;
    B_detect = 2:6;
    nbands = 7;
    conse = 4;
	def_pT_cg = T_cg;
	def_conse = conse;
    
    rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'pos',[],...
        'rmse',[], 'rmse_cnt',[],'num_obs',[],'change_prob',[],...
        'category',[],'magnitude',[],'initial',[],'obs_last',[]);
                            
    % Time Model initialization
    nInitial = 24;
    [fit_cft,rec_Flag,rec_cg_minit,dStart,dStable,BL_train] = TimeModelInitialize(...
        clrx, clry, pid, adj_rmse, conse, B_detect, T_cg,pos_PreBreak, nInitial);


    if BL_train == 0  &&  rec_Flag == 0     
        %% Case1: very noisy obs.
        % The time seires only has one time segment, and it is so noisy that
        % can not pass the stable test.
        iSaveObs = 1;
        rec_temp = SaveRecord_TotalFewObs(pid, clrx, clry,iSaveObs);
        rec_cg(num_fc) = rec_temp(1);
        return;
    end

    if BL_train == 0  && rec_Flag == 1    
        %% Case2: very noisy obs.
        % break detected in model initialization, but the remaining time seires 
        % is so noisy that can not pass the stable test.
        f_break = 1;      
        pos_break = find(clrx == rec_cg_minit(1).t_break);
        rec_cg(num_fc) = rec_cg_minit(1);
        num_fc = num_fc+1;
        iSaveObs = 1;
        i_start = find(clrx == dStart);
        rec_temp = SaveRecord_TotalFewObs(pid, clrx(i_start:end),... 
            clry(i_start:end,:),iSaveObs);
        rec_cg(num_fc) = rec_temp(1);    
        return;
    end

    
    if BL_train == 1  && rec_Flag == 1    
        %% Case3:  break detected in model initialization.
        f_break = 1;      
        pos_break = find(clrx == rec_cg_minit(1).t_break);
        rec_cg(num_fc) = rec_cg_minit(1);
        return;
    end

    if BL_train == 1  && rec_Flag == 0    % break not detected in model initialization
        
        %% Update Monitoring with model initialization
        % Find the start of stable period
        i_start = find(clrx == dStart);
        i_stable = find(clrx == dStable);

        % Recursive steps started from the stable start, and the index of
        % new stable period will be 1: i_stable-i_start+1.
        clrx_ff = clrx(i_start:end);
        clry_ff = clry(i_start:end,:);
        cur_mt = datevec(clrx_ff);
        clrx_doy = clrx_ff - datenum(cur_mt(:,1),1,0);  
        clrx_bin = floor(clrx_doy/8)+1;
        nLength = length(clrx_ff);
        theta0 = zeros(n_Model,nbands-1);    % variable for initial step
        rmse_FF = zeros(1, nbands-1);    % predicted rmse  
        yhat = zeros(nLength, nbands-1);    % estimates of FF                       
        v_dif_FF = zeros(nLength,nbands-1);    % diference of estimates  
        for i_B = B_detect
            % Get all estimates from forgetting factor
            theta0(:,i_B) = fit_cft(1:n_Model,i_B);
            [yhat(:,i_B), v_dif_FF(:,i_B), rmse_FF(1,i_B)] = ...
                autoFFsrFit(clrx_ff, clry_ff(:,i_B),theta0(:,i_B),...
                n_Model, iFFactor, preLength);
        end

        % Initial look-up table of RMSE
        LUT_sum = zeros(46,nbands-1);
        LUT_cnt = zeros(46,1);   
        for iDs = 1:i_stable-i_start+1
            curM = clrx_bin(iDs);
            LUT_cnt(curM) = LUT_cnt(curM) + 1;     
            for i_B = B_detect
                LUT_sum(curM,i_B) = LUT_sum(curM,i_B) + ...
                    (v_dif_FF(iDs, i_B))*(v_dif_FF(iDs, i_B));                            
            end
        end
              

        %% Detect a changeLog using the while loope
        % Start from the point following the stable end point.
        i = i_stable-i_start+2; 
        mini_rmse = adj_rmse;
        adj_conse = def_conse;
        %% while loop to detect change
        while i <= nLength - def_conse                
            % Adjusted conse according to change length           
            currdate = clrx_ff(i)+t_span;
            if currdate > clrx_ff(nLength)
                break;    % Date exceedes the time seires!!
            end       
            index_delta = find(clrx_ff(i:nLength) >= currdate);
            if length(index_delta)<1
                break;    % No points!!
            end
            adj_conse = index_delta(1);   % First point > minimum time interval.    
            if adj_conse < def_conse
                adj_conse = def_conse;
            end          
            if i+ adj_conse - 1> nLength
                break;    % Last point exceedes the time seires!!
            end
            
            % Adjust T_cg based on conse and adj_conse
            T_cg = chi2inv(def_pT_cg,length(B_detect));             
            if adj_conse > def_conse
                pT_cg = 1-(1-def_pT_cg)^(def_conse/adj_conse);
                T_cg = chi2inv(pT_cg,length(B_detect));
            end


            % Calculate RMSE from the LUT
            mini_rmse = ExtractRMSE_LUTbin(clrx_bin(i),LUT_sum,LUT_cnt);

            % Calculate change metrix for adj_conse observations
            vec_mag = zeros(adj_conse,1);
            v_dif = zeros(adj_conse,nbands-1);  
            for i_conse = 1:adj_conse           
                for i_B = B_detect
                    if sum(i_B == B_detect)   
                        v_dif(i_conse,i_B) = abs(v_dif_FF(i+i_conse-1,i_B))/max(adj_rmse(i_B),mini_rmse(i_B));
                    end  % if i_B
                end  % for i_B
                vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
            end  % for i_conse
               
            
            % Detect a change
            if min(vec_mag) > T_cg  
                % Updating output information
                rec_cg(num_fc).t_break = clrx_ff(i);    
                rec_cg(num_fc).t_start = clrx_ff(p_start);             
                rec_cg(num_fc).t_end = clrx_ff(i-1);                                                  
                rec_cg(num_fc).pos = pid;   
                rec_cg(num_fc).num_obs = length(p_start:i-1);    
                rec_cg(num_fc).change_prob = 1;      
                rec_cg(num_fc).initial = 1; 
                rec_cg(num_fc).rmse = LUT_sum;   
                rec_cg(num_fc).rmse_cnt = LUT_cnt;
                rec_cg(num_fc).obs_last = 0;

                % Determine the time series model
                IDs = p_start:i-1; 
                i_span = i-1-p_start+1;              
                update_num_c = update_cft(i_span,n_times,min_num_c,...
                    mid_num_c,max_num_c,max_num_c);
                outfit_cft = zeros(8,nbands-1); 
                outrmse = zeros(1, nbands-1); 
                for i_B = 1:6
                    [outfit_cft(:,i_B),outrmse(1,i_B),~] = autoTSFit(...
                        clrx_ff(IDs), clry_ff(IDs, i_B), update_num_c);
                end         
                rec_cg(num_fc).coefs = outfit_cft;   
                rec_cg(num_fc).category =   update_num_c;

                % change vector magnitude
                v_dif_mag = zeros(conse,nbands-1);
                for i_conse = 1:conse
                    for i_B = 1:nbands-1
                        % absolute difference
                        v_dif_mag(i_conse,i_B) = clry_ff(i+i_conse,i_B)-...
                            autoTSPred(clrx_ff(i+i_conse), outfit_cft(:,i_B));
                    end
                end
                rec_cg(num_fc).magnitude = median(v_dif_mag,1);  

                f_break = 1;
                pos_break = i+i_start-1;
                return;  % Return the change record!!! 
            end  % if T_cg 

            % Update the LUT
            curM = clrx_bin(i);
            LUT_cnt(curM) = LUT_cnt(curM) + 1;     
            for i_B = B_detect
                LUT_sum(curM,i_B) = LUT_sum(curM,i_B) + (v_dif_FF(i, i_B))*(v_dif_FF(i, i_B));                    
            end
            i = i+1;
        end  % end WHILE


        %% Process the end
        id_last = nLength;
        for i_conse = nLength:-1:i    
            for i_B = B_detect
                v_dif(i_conse,i_B) = abs(v_dif_FF(i_conse,i_B))/max(adj_rmse(i_B),mini_rmse(i_B));
            end
            vec_mag = norm(v_dif(i_conse,B_detect))^2; 
            if vec_mag < T_cg
                % the last stable id
                id_last = i_conse;
                break;
            end
        end

        % No change is detected for the last time segment.
        if f_break == 0  
            % Set output information
            IDs = p_start:id_last; 
            rec_cg(num_fc).t_break = 0;    
            rec_cg(num_fc).t_start = clrx_ff(p_start);             
            rec_cg(num_fc).t_end = clrx_ff(id_last);                                                  
            rec_cg(num_fc).pos = pid;   
            rec_cg(num_fc).num_obs = length(IDs);      
            rec_cg(num_fc).initial = 1; 
            rec_cg(num_fc).rmse = LUT_sum;   
            rec_cg(num_fc).rmse_cnt = LUT_cnt;

            % Determine the time series model       
            i_span = id_last-p_start+1;              
            update_num_c = update_cft(i_span,n_times,min_num_c,...
                mid_num_c,max_num_c,max_num_c);
            outfit_cft = zeros(8,nbands-1); 
            outrmse = zeros(1, nbands-1); 
            for i_B = 1:6
                [outfit_cft(:,i_B),outrmse(1,i_B),~] = autoTSFit(...
                    clrx_ff(IDs), clry_ff(IDs, i_B), update_num_c);
            end         
            rec_cg(num_fc).coefs = outfit_cft;   
            rec_cg(num_fc).category =   update_num_c;
            rec_cg(num_fc).magnitude = 0;  
            rec_cg(num_fc).change_prob= 0;
            rec_cg(num_fc).obs_last = 0;

            if id_last < nLength    % potential change at the end
                rec_cg(num_fc).t_break = clrx_ff(id_last+1); 
                rec_cg(num_fc).initial = 0; 

                cp_conse = (nLength-id_last+1)/adj_conse;
                cp_tspan = (clrx(end)-clrx(id_last))/t_span;
                rec_cg(num_fc).change_prob=min(cp_conse,cp_tspan);

                if nLength>=12
                    obs_last = zeros(12,nbands-1);
                    obs_last(:,1:6) = clry_ff(end-11:end,1:6);
                    obs_last(:,7) = clrx_ff(end-11:end);
                else
                    obs_last = zeros(nLength,nbands-1);
                    obs_last(:,1:6) = clry_ff(:,1:6);
                    obs_last(:,7) = clrx_ff;
                end
                rec_cg(num_fc).obs_last = obs_last;
            end % if id_last
            return;
        end  % if f_break 
    % BL_train end
    end 
 % function end     
end       



function rmse = ExtractRMSE_LUTbin(doy_bin,LUT_sum,LUT_cnt)
%% ExtractRMSE_LUTbin
% This function extracts the RMSE from the LUT at a given DOY bin.
%                      --- Developed by Rong Shang, 2020/08/25.

    % Set the parameters
    nLeast = 24;
    nBin = 46;
    fBin = floor(nBin/2);
    nBands = 7;
    doy_bin = doy_bin+fBin;
    
    % Expand the LUT
    sumLUT = zeros(nBin*2,nBands-1);
    sumLUT(1:fBin,:)=LUT_sum(fBin+1:nBin,:);
    sumLUT(fBin+1:nBin+fBin,:)=LUT_sum;
    sumLUT(nBin+fBin+1:nBin*2,:)=LUT_sum(1:fBin,:);    
    cntLUT = zeros(nBin*2,1);
    cntLUT(1:fBin)=LUT_cnt(fBin+1:nBin,1);
    cntLUT(fBin+1:nBin+fBin)=LUT_cnt(:,1);
    cntLUT(nBin+fBin+1:nBin*2)=LUT_cnt(1:fBin,1);
    
    % Search the LUT index
    final_index = fBin+1:fBin+nBin;    
    if sum(cntLUT) > nLeast  % total > numLeast
        if cntLUT(doy_bin) >= nLeast
            final_index = doy_bin;
        else
            step_s = 1;
            cur_index = doy_bin-step_s:doy_bin+step_s;
            cur_sum = sum(cntLUT(cur_index));
            while cur_sum < nLeast && step_s<fBin-1
                step_s = step_s+1;
                cur_index = doy_bin-step_s:doy_bin+step_s;
                cur_sum = sum(cntLUT(cur_index));
            end
                final_index = cur_index;
        end
    end
 
    % Calculate RMSE
    rmse = zeros(6,1);
    for i_B = 1:6
        rmse(i_B) =  sqrt(sum(sumLUT(final_index,i_B))/sum(cntLUT(final_index)));
    end     
    
% function end
end


function [fit_cft, rec_Flag, rec_cg, dateStart, dateStable,BL_train] = TimeModelInitialize(...
    clrx, clry, pid, adj_rmse, conse, B_detect, T_cg, pos_PreBreak,lengthInitial)
%% TimeModelInitialize:
%     This function aims to find the stable start date of the time series 
% for COLDFF. Except for returning the fitting parameters and stable 
% start date, it also returns the stuct of rec_cg if a breakpoint was found.
%
% Function Inputs:
%    clrx --- Date of input reflectance;
%    clry --- Reflectance values for the detected bands;
%    pid --- Index for a pixel;
%    conse --- Number of continuous observations used for change detection;
%    B_detect --- Bands used for change dettection(e.g. 2:5);
%    adj_rmse --- Minimum RMSE used for change detection;
%    lengthInitial --- Length for Time Model Initialization.
%
% Fuction output (rec_cg):
%    fit_cft --- Fitting parameters from the stable start
%    rec_Flag --- Whether a change has been detected (1 is for change);
%    rec_cg --- Output struct of a change record;
%    dateStart --- Date of the start point;
%    dateStable --- Date of the stable point.
  

%% ChangeLog 11/26/2019:
% Add "lengthInitial" to this function.
% Use adjusted number of paramters for the time model.



    nid = pid;
    dateStable = clrx(1);
    dateStart = clrx(1);
	i_Stable_Flag = 0;
	rec_Flag = 0;
    nbands = 8;
	min_num_c = 4;
	mid_num_c = 6;
	max_num_c = 8;
	n_times = 3;	
	num_yrs = 365.25;    % number of days per year	
    fit_cft = zeros(max_num_c,nbands-1);
	
    rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'pos',[],...
        'rmse',[], 'rmse_cnt',[],'num_obs',[],'change_prob',[],...
        'category',[],'magnitude',[],'initial',[],'obs_last',[]);
    
	
	% adjust threshold based on chi-squared distribution
	def_pT_cg = T_cg;
	T_cg = chi2inv(def_pT_cg,length(B_detect));

	% initializing variables
	% the first observation for TSFit
	i_start = 1;
	% record the start of the model initialization (0=>initial;1=>done)
	BL_train = 0;
	% identified and move on for the next curve
	num_fc = 1; % NUM of Fitted Curves (num_fc)
	% initialize i_dense
	i_dense = 1;
	
	% While loope: process till the last clear observation
    % i starts with the miminum numbers of clear obs
	i = lengthInitial;
	while i<= length(clrx)
		% span of "i"
		i_span = i-i_start+1;	
		% span of time (num of years)
		time_span = clrx(i)-clrx(i_start);
		% max time difference
		time_dif = max(clrx(i_start+1:i) - clrx(i_start:i-1));
		
		% basic requrirements: 1) enough observations; 2) enough time
        if i_span < lengthInitial || time_span < num_yrs
            i=i+1;
            continue;
        end
        
        % Initializing model
        if BL_train == 0
            % check max time difference
            if time_dif > num_yrs
                i = i+1;
                i_start = i_start+1;
                % i that is dense enough
                i_dense = i_start;
                continue
            end % time_dif

            % Step 2: model fitting
            % defining computed variables
            fit_cft = zeros(max_num_c,nbands-1);          
            rmse = zeros(nbands-1,1);
            v_dif = zeros(nbands-1,1);
            rec_v_dif = zeros(i-i_start+1,nbands-1);

            i_span = i_start:i;
            update_num_c = update_cft(i_span,n_times,min_num_c,...
                mid_num_c,max_num_c,max_num_c);

            % normalized to z-score
            for i_B = B_detect
                % initial model fit
                [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)]=...
                    autoTSFit(clrx(i_start:i),clry(i_start:i,i_B),update_num_c);

                % minimum rmse
                mini_rmse = max(adj_rmse(i_B),rmse(i_B));

                % compare the first clear obs
                v_start = rec_v_dif(1,i_B)/mini_rmse;
                % compare the last clear observation
                v_end = rec_v_dif(end,i_B)/mini_rmse;
                % anormalized slope values
                v_slope = fit_cft(2,i_B)*...
                    (clrx(i)-clrx(i_start))/mini_rmse;

                % differece in model intialization
                v_dif(i_B) = abs(v_slope) + max(abs(v_start),abs(v_end));
            end
            v_dif = norm(v_dif(B_detect))^2;
						
            % Stable Test
            if v_dif > T_cg    % Not Stable
                % start from next clear obs
                i_start = i_start + 1;
                % move forward to the i+1th clear observation
                i = i + 1;
                % keep all data and move to the next obs
                continue
            else     % Stable
                % model ready!
                BL_train = 1;

                 % Record the position of stable end
                i_stable = i;                
                dateStable = clrx(i_stable);

                % Find the previous break point
                if pos_PreBreak == 1
                    i_break = 1;
                else
                    i_break = pos_PreBreak;
                end

                % Backword direction: change detection
                if i_start > i_break                               
                    % model fit at the beginning 
                    for i_ini = i_start-1:-1:i_break
                        if i_start - i_break < conse
                            ini_conse = i_start - i_break;
                        else
                            ini_conse = conse;
                        end
                        % value of difference for conse obs
                        v_dif = zeros(ini_conse,nbands-1);
                        % record the magnitude of change
                        v_dif_mag = v_dif;
                        % chagne vector magnitude
                        vec_mag = zeros(ini_conse,1);

                        for i_conse = 1:ini_conse
                            for i_B = 1:6
                                % absolute difference
                                v_dif_mag(i_conse,i_B) = ...
                                    clry(i_ini-i_conse+1,i_B)-...
                                    autoTSPred(clrx(i_ini-i_conse+1),fit_cft(:,i_B));
                                % normalized to z-scores
                                if sum(i_B == B_detect)
                                    % minimum rmse
                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));

                                    % z-scores
                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/mini_rmse;
                                end
                            end
                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                        end % i_conse


                        % change detection
                        if min(vec_mag) > T_cg 
                            % the points before current i_start can not be merged.
                            i_Stable_Flag = 1;                                        
                            break;    % break from for i_ini = i_start-1:-1:i_break
                        end

                        % update i_start (merge previous point and move back i_start)
                        i_start = i_ini;   
                    end
                end % i_start
				dateStart= clrx(i_start);
				

                % If there are more than 4 obs, record a change event in first 
				% model initialization of a pixel
                if pos_PreBreak == 1 && i_Stable_Flag == 1 && i_start > 4
                    rec_Flag = 1;		
                    % record time of curve start
                    rec_cg(num_fc).t_start = clrx(1);
                    % record time of curve end
                    rec_cg(num_fc).t_end = clrx(i_start-1);								
                    % record break time
                    rec_cg(num_fc).t_break = clrx(i_start);
                    % record postion of the pixel
                    rec_cg(num_fc).pos = nid;
                    % record change probability
                    rec_cg(num_fc).change_prob = 1;
                    % record number of observations
                    rec_cg(num_fc).num_obs = i_start - i_dense;
                    % record change magnitude
                    rec_cg(num_fc).magnitude = - median(v_dif_mag,1);

                    if i_start - i_dense >= 4								
                        % defining computed variables
                        fit_cft_pre = zeros(max_num_c,nbands-1);
                        % rmse for each band
                        rmse = zeros(nbands-1,1);

                        for i_B=1:6
                            [fit_cft_pre(:,i_B),rmse(i_B)] = ...
                                autoTSFit(clrx(i_dense:i_start-1),...
                                clry(i_dense:i_start-1,i_B),min_num_c);
                        end
                        cntLength =length(i_dense:i_start-1);
                        % record fitted coefficients
                        rec_cg(num_fc).coefs = fit_cft_pre;
                        % record rmse of the pixel
                        rec_cg(num_fc).rmse = rmse;
                        rec_cg(num_fc).rmse_cnt = cntLength;
                        % record fit category
                        rec_cg(num_fc).category =   min_num_c;
                        rec_cg(num_fc).initial = 1;
                    end						

                    % identified and move on for the next functional curve
                    num_fc = num_fc + 1;
                end % pos_PreBreak

            end
        end


		%  Model initialized, and obtain a stable period	
        if BL_train == 1 && i_Stable_Flag == 1
            dateStart= clrx(i_start);
            i_span = i_stable-i_start+1;           
            if i_span >= 6 % Estimtae a time model from at least 6 observations
                update_num_c = update_cft(i_span,n_times,min_num_c,...
                    mid_num_c,max_num_c,max_num_c);
                % Ountput: fit_cft will be used as initial model paramters of
                % Forgetting Factor.
                fit_cft = zeros(max_num_c,nbands-1);
                rmse = zeros(nbands-1,1);
                for i_B=1:6
                    [fit_cft(:,i_B),rmse(i_B)] = autoTSFit(clrx(i_start:i_stable),...
                        clry(i_start:i_stable,i_B),update_num_c);
                end
            end % i_span
            break;    % Break from the while loope
        end % BL_train

		% move forward to the i+1th clear observation
		i=i+1;
	end % end of while iterative       	
end


function rec_cg = SaveRecord_TotalFewObs(pid, clrx, clry,iSaveObs)
%% SaveRecord_TotalFewObs
% Save the record with very few observations in the time series. A time
% model will be fitted if there are at least 6 observations.

    rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'pos',[],...
        'rmse',[], 'rmse_cnt',[],'num_obs',[],'change_prob',[],...
        'category',[],'magnitude',[],'initial',[],'obs_last',[]);
    
    p_start = 1;   
	num_fc =1;
    iTotal = length(clrx);    
    numLowest = 6;    
    min_num_c = 4;    
    mid_num_c = 6;
    max_num_c = 8;
    n_times = 3;  
    nbands = 7;
    
    rec_cg(num_fc).t_break = 0;   
    rec_cg(num_fc).t_start=clrx(p_start);                
    rec_cg(num_fc).t_end=clrx(iTotal);                                          
    rec_cg(num_fc).pos = pid;     
    rec_cg(num_fc).rmse_cnt =0;
    rec_cg(num_fc).category = 0;

    % Determine the time series model
    IDs = p_start:iTotal;  
    i_span = iTotal-p_start+1;       
    if i_span < numLowest
        rec_cg(num_fc).coefs = zeros(8,nbands-1); 
        rec_cg(num_fc).rmse = 0;      
    end
    
    if i_span >= numLowest
        outfit_cft = zeros(8,nbands-1); 
        outrmse = zeros(1, nbands-1);  
        update_num_c = update_cft(i_span,n_times,min_num_c,...
            mid_num_c,max_num_c,max_num_c);
        for i_B = 1:6
            [outfit_cft(:,i_B), outrmse(1,i_B),~] = autoTSFit(...
                clrx(IDs), clry(IDs, i_B), update_num_c);
        end         
        rec_cg(num_fc).coefs = outfit_cft;                
        rec_cg(num_fc).rmse = outrmse(1, :) ;   
        rec_cg(num_fc).category =   update_num_c;
    end

    
    rec_cg(num_fc).num_obs = iTotal;      
    rec_cg(num_fc).change_prob = 0; 
    rec_cg(num_fc).initial = 0;   
    rec_cg(num_fc).magnitude = zeros(1,nbands-1);  
    if iSaveObs == 1
        olast = zeros(iTotal,6);
        olast(:,7) = clrx;
        olast(:,1:6) = clry(:,:);
        rec_cg(num_fc).obs_last = olast;
    end
end