clear; % Clear variables
addpath('../data')
datasetNum = 4; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime,proj2Data] = init(datasetNum);

Z = sampledVicon(1:6,:);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.1*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0; %last time step in real time
pos = proj2Data.position;
pose = proj2Data.angle;
for i = 1:length(sampledTime)
    %% Fill in the FOR LOOP
    angVel = sampledData(i).omg; % Angular Velocity from IMU
    acc = sampledData(i).acc; %Linear accelaration from IMU
    if i == 1
        dt = sampledTime(1)-0;
    else
        dt = sampledTime(i)-sampledTime(i-1);
    end
    [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt); %Prediction Step
    z_t = [pos(i,:)';pose(i,:)'];
    [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst); %Update Step

    uPrev = uCurr;          %Updating the variables for the new iteration.
    covarPrev = covar_curr;

    savedStates(:,i) = uCurr;
end

plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);