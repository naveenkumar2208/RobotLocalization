function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively
Ct = [eye(3) zeros(3) zeros(3) zeros(3) zeros(3)
      zeros(3) eye(3) zeros(3) zeros(3) zeros(3)];
R = eye(6)*0.0001;
% R = normrnd(0,0.001,[6,6]);
Kt = covarEst*Ct'/(Ct*covarEst*Ct'+R);
uCurr = uEst + Kt*(z_t - uEst(1:6));
covar_curr = covarEst - Kt*Ct*covarEst;

end