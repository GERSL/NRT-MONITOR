function rmse = ExtractRMSE_LUTbin(doy_bin,LUT_sum,LUT_cnt)
%% ExtractRMSE_LUTbin
% This function extracts the RMSE from the LUT at a given DOY bin.


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

