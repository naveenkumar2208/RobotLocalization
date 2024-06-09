%% PROJECT 2 VELOCITY ESTIMATION
tic
close all;
clear all;
clc;
addpath('../data')

%Change this for both dataset 1 and dataset 4. Do not use dataset 9.
datasetNum = 4;

[sampledData, sampledVicon, sampledTime] = init(datasetNum);

%% INITIALIZE CAMERA MATRIX AND OTHER NEEDED INFORMATION
CamCallib = [311.0520        0        201.8724;

0         311.3885    113.6210;

0            0           1   ];

Rot = [cosd(45) -sind(45) 0 %Rotation from camera to body
    -sind(45) -cosd(45) 0
    0 0 -1];

Trans = [-0.04 %Translation from Camera to Body
    0
    -0.03];

Tcb = zeros(4,4);
Tcb(1:3,1:3) = Rot;
Tcb(1:3,4) = Trans; %Transformation from Camera to Body frame
Tcb(4,4) = 1;

for n = 2:length(sampledData)
    %% Initalize Loop load images
    Corntminus1 = detectFASTFeatures(sampledData(n-1).img); % Corner detection in both of the images
    CornerCurr = detectFASTFeatures(sampledData(n).img);
    %% Detect good points
    a = min(length(Corntminus1),length(CornerCurr));
    a = floor(a*(60/140));
    Corntminus1 = Corntminus1.selectStrongest(a); % Chooseing Good points for the computation
    CornerCurr = CornerCurr.selectStrongest(a);
    %% Initalize the tracker to the last frame.
    pointTracker = vision.PointTracker('MaxBidirectionalError',3); %Tracking using Computer Vision Toolbox
    initialize(pointTracker,Corntminus1.Location,sampledData(n-1).img)
    %% Find the location of the next points;
    [points,point_validity] = pointTracker(sampledData(n).img); % Finding the location of next points in new image
    %% Calculate velocity
    % Use a for loop
    Corntminus1(point_validity == 0) = [];  
    points(point_validity == 0,:) = [];
    point_validity(point_validity == 0) = [];  %Deleting those points which has point validity = 0

    velocity = zeros(length(point_validity),3);

    dt = (sampledData(n).t-sampledData(n-1).t);

    for i = 1:length(point_validity)
        a = CamCallib\[points(i,:)';1];   
        b = CamCallib\[Corntminus1.Location(i,:)';1];
        velocity(i,:) = (a' - b')/dt; % Pixel Velocity in Camera Frame
    end


    %% Calculate Height
    t = n-1;
    [position, orientation, R_c2w] = estimatePose(sampledData, t);

    Twb = zeros(4,4);
    Twb(1:3,1:3) = eul2rotm(orientation);
    Twb(1:3,4) = position;
    Twb(4,4) = 1; %Transformation Matrix from Body to World Frame.

    Tcw = Tcb/Twb; %Transformation Matrix from World to Camera Frame.

    Twc = Twb/Tcb; %Transformation Matrix from Camera to World Frame
    XYZ = CamCallib*Tcw(1:3,4);

    Z = XYZ(3); % Z depth of the middle pixel





    %% Thereshold outputs into a range.
    % Not necessary

    %% Fix the linear velocity
    % Change the frame of the computed velocity to world frame
    vel = zeros(2*length(point_validity),1);
    functi = zeros(2*length(point_validity),6);
    for i = 1:length(point_validity)
        a = 2*i-1;
        b = 2*i;
        xyz = CamCallib\[Corntminus1.Location(i,:) 1]';
        x = xyz(1);
        y = xyz(2);
        vel(a:b,1) = velocity(i,1:2)';
        functi(a,:) = [(-1/Z) 0 (x/Z) x*y -(1+(x*x)) y];
        functi(b,:) = [0 (-1/Z) (y/Z) (1+(y*y)) -x*y -x];  
    end
    velo = pinv(functi)*vel; %Calculating Camera velocity in Camera Frame.

    Adj=[ Twc(1:3,1:3) zeros(3,3); zeros(3,3)  Twc(1:3,1:3)];
    Vel = Adj*velo; %Calculating Camera Velocity in the World Frame.    




    %% ADD SOME LOW PASS FILTER CODE
    % Not neceessary but recommended

    % estimatedV(:,n) = Vel;

    %% STORE THE COMPUTED VELOCITY IN THE VARIABLE estimatedV AS BELOW
    estimatedV(:,n) = Vel; % Feel free to change the variable Vel to anything that you used.
    % Structure of the Vector Vel should be as follows:
    % Vel(1) = Linear Velocity in X
    % Vel(2) = Linear Velocity in Y
    % Vel(3) = Linear Velocity in Z
    % Vel(4) = Angular Velocity in X
    % Vel(5) = Angular Velocity in Y
    % Vel(6) = Angular Velocity in Z
end


plotData(estimatedV, sampledData, sampledVicon, sampledTime, datasetNum)
toc