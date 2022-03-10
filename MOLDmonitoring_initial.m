function [f_break, pos_break, rec_cg] = MOLDmonitoring_initial(pid, clrx, clry,...
    T_cg, t_span, iFFactor,preLength,adj_rmse,rpre_coefs,rpre_rmse,rpre_rmse_cn)
%% 
%     This function conducts MOLD monitoring with initial parameters to 
% detect a change. When a change is detected, this function will return the 
% position of change in the time series. Noticed that this function only
% detect one change. 



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