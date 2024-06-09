function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
%% Parameter Definition
% uPrev - is the mean of the prev state
%covarPrev - covar of the prev state
%angVel - angular velocity input at the time step
%acc - acceleration at the timestep
%dt - difference in time

ng = [0.2
    0.001
   0.4];
na = [0.0002    % Values of different noise set after different hit and trials and comparision of vicon data.
    0.001
    0.001];
nbg = [0.01
    0.01
    0.1];
nba = [0.0015   % Values of different bias set after different hit and trials and comparision of vicon data.
   0.0001
   0.001];

n = 15; %number of variables in the uPrev
nq = 12; %number of elements in q noise
ndash = n+nq;

alpha = 0.001;
beta = 2;
kappa = 1;

lambda = ((alpha^2)*(ndash+kappa)) - ndash;

Q = [diag(ng) zeros(3,3) zeros(3,3) zeros(3,3)
    zeros(3,3) diag(na) zeros(3,3) zeros(3,3)
    zeros(3,3) zeros(3,3) diag(nbg) zeros(3,3)
    zeros(3,3) zeros(3,3) zeros(3,3) diag(nba)];

Paug = [covarPrev zeros(15,12)
    zeros(12,15) Q];
uAugPrev = [uPrev
    zeros(12,1)];

X = zeros(27,55);

cholpaug = chol(Paug,"lower");

% Computing Sigma Points
for i = 1:55
    if i == 1
        X(:,i) = uAugPrev;
    elseif i >= 2 && i <= 28
        X(:,i) = uAugPrev + (sqrt(ndash + lambda)*cholpaug(:,i-1));
    elseif i > 28 && i <= 55
        X(:,i) = uAugPrev - (sqrt(ndash + lambda)*cholpaug(:,rem(i,28)));
    end
end
%Propogation of each of the Sigma points
F = zeros(15,55);

for i = 1:55

    x3 = [X(7,i);X(8,i);X(9,i)];
    x4 = [X(10,i);X(11,i);X(12,i)];
    x5 = [X(13,i);X(14,i);X(15,i)];

    roll = X(4,i); %Phai
    pitch = X(5,i); %theta
    yaw = X(6,i); %yaw

     R = eul2rotm([yaw pitch roll]); % Rotation matrix defines the body in the world frame

    G = [cos(pitch)*cos(yaw), -sin(yaw), 0;
    cos(pitch)*sin(yaw), cos(yaw), 0;       %inverse of G Matrix
    -sin(pitch), 0, 1]\R;

    f = [x3
        G*(angVel-x4-X(16 : 18, i)) 
        R*(acc-x5-X(19 : 21, i))-[0;0;9.81] % Process model
        X(22 : 24, i) 
        X(25 : 27, i)]*dt; 

    F(:,i) = X(1:15,i) + f; %Predicted data from differen sigma points.
end
Res = zeros(15,1);

for i = 1:55
    if i == 1
        W = lambda/(lambda+ndash);
        Res = Res+(W*F(:,i)); %Computing the final predicted mean by multiplying respective weight of each sigma points and summing them up.
    elseif i >=2 && i <= 55
        W = 1/(2*(lambda+ndash));
        Res = Res+(W*F(:,i));
    end
end
uEst = Res;
ResC = zeros(15,15);

for i = 1 : 55
    if (i == 1)
        Wc = (lambda / (ndash + lambda)) + (1 - alpha^2 + beta);
        ResC = Wc * (F(:, i) - uEst) * transpose((F(:, i) - uEst));
    else  %%Computing the final predicted Covariance by multiplying respective weight of each sigma points and summing them up.
        Wc = 1/(2 * (ndash + lambda));
        ResC = ResC + (Wc * (F(:, i) - uEst) * transpose((F(:, i) - uEst)));
    end
end
covarEst = ResC;
end

