function main_NRT_MONITOR()
%% main_NRT_MONITOR
% The main function to run NRT-MONITOR based on previous monitoring results
% ended in the year 2014.

    % ******* Define the input and output *******
    maindir = 'D:\NRT-MONITOR-main\';
    input_Path = fullfile(maindir,'Sample_Data.mat'); 
    pre_Path = fullfile(maindir,'Sample_previousMonitoring_end2014.mat'); 
    output_Path = fullfile(maindir,'Sample_results_updatefrom2015.mat'); 
    pid = 1; 
   
    
    % ******* Default thresholds *******
    T_cg = 0.99999;
    t_span = 30;
    preLength = 120;
    iFFactor = 0.999;


    % ******* Read previous monitoring results ended in 2014 ******* 
    rec_pre = load(pre_Path);
    rec_pre = rec_pre.rec_pre;
    rec_pre = rec_pre(end);
    adj_rmse = rec_pre.rmse;
    
    
    % ******* Read the TRA-adjusted and outlier-removed observations ******* 
    obs_cut = load(input_Path);
    obs_TRAclear = obs_cut.obs_TRAclear;    
    clrx = obs_TRAclear(:,1); 
    % Select observations start from 2015
    obs_TRAclear = obs_TRAclear(clrx>datenum(2015,1,0),:);   
    clrx = obs_TRAclear(:,1); 
    clry= obs_TRAclear(:,2:7); 


    % ******* Remove duplicated observations (same observing date) *******
    [clrx,uniq_id,~] = unique(clrx);
    tmp_y = zeros(length(clrx),6);
    for i = 1:6
        tmp_y(:,i) = clry(uniq_id,i);
    end
    clry = tmp_y;
    

    % ******* Run NRT-MONITOR and save results *******
    rec_cg = NRTMONITOR_Pixel(rec_pre, clrx, clry, pid, adj_rmse, T_cg, t_span, preLength,iFFactor);
    save(output_Path,'rec_cg');
end
