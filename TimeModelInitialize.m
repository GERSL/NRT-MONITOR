function [fit_cft, rec_Flag, rec_cg, dateStart, dateStable,BL_train] = TimeModelInitialize(...
    clrx, clry, pid, adj_rmse, conse, B_detect, T_cg, pos_PreBreak,lengthInitial)
%% TimeModelInitialize:
%     This function aims to find the stable start date of the time series 
% for COLDFF. Except for returning the fitting parameters and stable 
% start date, it also returns the stuct of rec_cg if a breakpoint was found.
%
 


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