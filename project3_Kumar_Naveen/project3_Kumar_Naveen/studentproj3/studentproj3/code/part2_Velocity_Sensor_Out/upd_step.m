function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
%% Parameter Definition
%z_t - is the sensor data at the time step
%covarEst - estimated covar of the  state
%uEst - estimated mean of the state
Rot = [cosd(45) -sind(45) 0
    -sind(45) -cosd(45) 0  % Rotation matrix reprenting the body in Camera frame.
    0 0 -1];

Trans = [-0.04
    0        % Translation vector from Camera to Body Frame
    -0.03];

rbcb = -Rot'*Trans; %Translation vector from body to camera frame
Srbcb = [0 -rbcb(3) rbcb(2)
    rbcb(3) 0 -rbcb(1)
    -rbcb(2) rbcb(1) 0]; % Skew-symmatric matrix of translation vector from body to camera frame

n = 15;
alpha = 0.001;
beta = 2;
kappa = 1;
lambda = ((alpha^2)*(n+kappa)) - n;

cholpaug = chol(covarEst,"lower");
% Compute Sigma Points
X = zeros(15,31);
for i = 1:31
    if i == 1
        X(:,i) = uEst;
    elseif i >= 2 && i <= 16
        X(:,i) = uEst + (sqrt(n + lambda)*cholpaug(:,i-1));
    elseif i > 16 && i <= 31
        X(:,i) = uEst - (sqrt(n + lambda)*cholpaug(:,rem(i,16)));
    end
end

gfunc = zeros(3,31);
% Propogating the Sigma Points
for i = 1:31
    gfunc(:,i) = Rot*eul2rotm([X(6,i) X(5,i) X(4,i)])'*[X(7,i);X(8,i);X(9,i)] - Rot*Srbcb*Rot'*z_t(4:6);
end

Res = zeros(3,1);
% Computing the predicted mean, predicted covariance of the measurement and
% predicteing the cross covariance.
for i = 1:31
    if i == 1
        W = lambda/(lambda+n);
        Res = Res+(W*gfunc(:,i));
    elseif i >=2 && i <= 55
        W = 1/(2*(lambda+n));
        Res = Res+(W*gfunc(:,i));
    end
end
C = [];
for i = 1:31
    if i == 1
        W = lambda/(lambda+n);
        C = W*(X(:,i) - uEst)*(gfunc(:,i) - Res)';
    elseif i >=2 && i <= 31
        W = 1/(2*(lambda+n));
        C = C + W*(X(:,i) - uEst)*(gfunc(:,i) - Res)';
    end
end
R = 0.0001*eye(3);
S = [];
for i = 1:31
    if i == 1
        W = (lambda / (n + lambda)) + (1 - alpha^2 + beta);
        S = W*(gfunc(:,i) - Res)*(gfunc(:,i) - Res)' + R;
    elseif i >=2 && i <= 55
        W = 1/(2*(lambda+n));
        S = S + W*(gfunc(:,i) - Res)*(gfunc(:,i) - Res)' + R;
    end
end
Kt = C/S; % Kalman Gain
uCurr = uEst + Kt*(z_t(1:3)-Res);%Updating Current mean by using measurement from the Velocity of the camera frame with respect to the world frame expressed in the camera frame.
covar_curr = covarEst - Kt*S*Kt';

end

