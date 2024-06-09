function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    Ct = [eye(3) zeros(3) zeros(3) zeros(3) zeros(3)
        zeros(3) eye(3) zeros(3) zeros(3) zeros(3)]; %C matrix from Extended Kalman Filter
    R = eye(6)*0.0001;  % Noise in the measurement (assumed)
    Kt = covarEst*Ct'/(Ct*covarEst*Ct'+R); % Kalman Gain
    uCurr = uEst + Kt*(z_t - uEst(1:6)); %Updating Current mean by using measurement from the vision pose from the camera on the quadrotor.
    covar_curr = covarEst - Kt*Ct*covarEst; %Updating the current covariance.
    
end

