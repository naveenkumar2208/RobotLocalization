function [Vel] = velocityRANSAC(optV,optPos,Z,R_c2w,e)
%% CHANGE THE NAME OF THE FUNCTION TO velocityRANSAC
%% Input Parameter Description
% optV = The optical Flow
% optPos = Position of the features in the camera frame
% Z = Height of the drone
% R_c2w = Rotation defining camera to world frame
% e = RANSAC hyper parameter
CamCallib = [311.0520        0        201.8724;

0         311.3885    113.6210;

0            0           1   ];

psucess = 0.99; % probability of hitting at least one inlier set, i.e. the probability of success
m = 3;

k = ceil(log(1-psucess)/log(1-e^3));
s = struct([]);
%  indexk = zeros(k,7);
indexk = zeros(k,1);

randomindex = zeros(k,3);

for in=1:k
    rn = zeros(1,3);
    for i=1:m
        rn(i) = randi(length(optPos)); %Three Random Points
    end
    randomindex(in,:) = rn;
    vel = zeros(6,1);
    functi = zeros(6,6);
    for i = 1:3
        a = 2*i-1;
        b = 2*i;
        xyz = CamCallib\[optPos.Location(rn(i),:) 1]';
        x = xyz(1);
        y = xyz(2);
        vel(a:b,1) = optV(rn(i),1:2)';
        functi(a,:) = [(-1/Z) 0 (x/Z) x*y -(1+(x*x)) y];
        functi(b,:) = [0 (-1/Z) (y/Z) (1+(y*y)) -x*y -x];
    end
    velo = pinv(functi)*vel; %Velosity from three random points of camera in camera frame.

    vel = zeros(2*length(optV),1);
    functi = zeros(2*length(optV),6);
    for i = 1:length(optV)
        a = 2*i-1;
        b = 2*i;
        xyz = CamCallib\[optPos.Location(i,:) 1]';
        x = xyz(1);
        y = xyz(2);
        vel(a:b,1) = optV(i,1:2)';
        functi(a,:) = [(-1/Z) 0 (x/Z) x*y -(1+(x*x)) y];
        functi(b,:) = [0 (-1/Z) (y/Z) (1+(y*y)) -x*y -x];
    end
    Hi = functi;
    Hivelo = Hi*velo; %Velocity of camera from all the points in Camera frame
    c = 0;


    indexnumber = [];
    for i = 1:length(optV)
        a = 2*i-1;
        b = 2*i;
        del = norm([Hivelo(a) - vel(a); ...
            Hivelo(b) - vel(b)].^2);
        if del<10^(-4)                         %RANSAC
            c = c+1;
            indexnumber(end+1) = i;
        end
    end

    s(in).c = c;
    s(in).indexnumber = indexnumber;

    
    indexk(in,1) = c;
   
end


%% Output Parameter Description
% Vel = Linear velocity and angualr velocity vector
[~,I] = max(indexk(:,1));
indexnumber = s(I).indexnumber;
vel = zeros(2*length(indexnumber),1);
functi = zeros(2*length(indexnumber),6);

for i = 1:length(indexnumber)
    a = 2*i-1;
    b = 2*i;
    xyz = CamCallib\[optPos.Location(indexnumber(i),:) 1]';
    x = xyz(1);
    y = xyz(2);
    vel(a:b,1) = optV(indexnumber(i),1:2)';
    functi(a,:) = [(-1/Z) 0 (x/Z) x*y -(1+(x*x)) y];
    functi(b,:) = [0 (-1/Z) (y/Z) (1+(y*y)) -x*y -x]; % RANSAC with the selected points which are less than particular threshold.
end

veloci = pinv(functi)*vel; % Final Velocities in Camera Frame.

Adj=[ R_c2w' zeros(3,3); zeros(3,3)  R_c2w];
Vel = Adj*veloci; % Final Velocities in World Frame.
end