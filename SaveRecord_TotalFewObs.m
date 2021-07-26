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
