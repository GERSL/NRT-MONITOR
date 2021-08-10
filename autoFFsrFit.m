function [yhat, v_dif, rmse] = autoFFsrFit(x, y, theta0, df, iff, preLength)
%% 
% This function gets the etimations using the recursive approach (forgetting factor). 
% It is an important step for detecting a change. 
%

% --------------------------------- Fuction Start ------------------------------------%
    obj = recursiveLS(df, theta0,'InitialParameterCovariance',2);              
    obj.ForgettingFactor = iff;                

    n=length(x);                          
    w=2*pi/365.25;                    
    yhat = zeros(n,1);                 
    theta = zeros(n,8); 
    
    iFlag_PreCreate = 1;
    if iFlag_PreCreate==1
        i_start_date = x(1);
        xx0 = i_start_date-preLength*8:8:i_start_date-1;
        yy0 = autoTSPred(xx0',theta0);

        clrx_new = [xx0'; x(1:end)];
        clry_new = [yy0; y(1:end)];

        yhat_new = zeros(length(clrx_new),1);


        for i = 1:length(clrx_new)%n
            cos_x = cos(w*clrx_new(i));
            sin_x = sin(w*clrx_new(i));

            if(df==4)  % 4 parameters
                H=[1 clrx_new(i) cos_x sin_x];
            end

            if(df==6)   % 6 parameters
                cos_2x = cos(2*w*clrx_new(i));
                sin_2x = sin(2*w*clrx_new(i));
                H=[1 clrx_new(i) cos_x sin_x cos_2x sin_2x];
            end

            if(df==8)   % 8 parameters
                cos_2x = cos(2*w*clrx_new(i));
                sin_2x = sin(2*w*clrx_new(i));
                cos_3x = cos(3*w*clrx_new(i));
                sin_3x = sin(3*w*clrx_new(i));
                H=[1 clrx_new(i) cos_x sin_x cos_2x sin_2x cos_3x sin_3x];
%                 H=[1 0 cos_x sin_x cos_2x sin_2x cos_3x sin_3x];
            end

            % Recursive step
            [theta(i,:), EstimatedOutput] = obj(clry_new(i),H);
            yhat_new(i)=EstimatedOutput;
        end
        yhat = yhat_new(preLength+1:end);
    end
    
    if iFlag_PreCreate ==0    
        for i = 1:n
            cos_x = cos(w*x(i));
            sin_x = sin(w*x(i));

            if(df==4)  % 4 parameters
                H=[1 x(i) cos_x sin_x];
            end

            if(df==6)   % 6 parameters
                cos_2x = cos(2*w*x(i));
                sin_2x = sin(2*w*x(i));
                H=[1 x(i) cos_x sin_x cos_2x sin_2x];
            end

            if(df==8)   % 8 parameters
                cos_2x = cos(2*w*x(i));
                sin_2x = sin(2*w*x(i));
                cos_3x = cos(3*w*x(i));
                sin_3x = sin(3*w*x(i));
                H=[1 x(i) cos_x sin_x cos_2x sin_2x cos_3x sin_3x];
            end

            % Recursive step
            [theta(i,:), EstimatedOutput] = obj(y(i),H);
            yhat(i)=EstimatedOutput;
        end
    end
    v_dif = y-yhat;                                   
    rmse=norm(v_dif)/sqrt(n-df);        
    release(obj);               
end