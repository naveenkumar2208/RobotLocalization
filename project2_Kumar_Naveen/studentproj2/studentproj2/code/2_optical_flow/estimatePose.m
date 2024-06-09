function [position, orientation, R_c2w] = estimatePose(data, t)
%% CHANGE THE NAME OF THE FUNCTION TO estimatePose
% Please note that the coordinates for each corner of each AprilTag are
% defined in the world frame, as per the information provided in the
% handout. Ideally a call to the function getCorner with ids of all the
% detected AprilTags should be made. This function should return the X and
% Y coordinate of each corner, or each corner and the centre, of all the
% detected AprilTags in the image. You can implement that anyway you want
% as long as the correct output is received. A call to that function
% should made from this function.
id = data(t).id;
res = getCorner(id);
    %% Input Parameter Defination
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset
    %% Output Parameter Defination
    
    % position = translation vector representing the position of the
    % drone(body) in the world frame in the current time, in the order XYZ
    
    % orientation = euler angles representing the orientation of the
    % drone(body) in the world frame in the current time, in the order XYZ
    
    %R_c2w = Rotation which defines camera to world frame
    %% Input Parameter Defination
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset
    
    %% Output Parameter Defination
    A = zeros(length(id)*10,9);
    for i=1:length(id)
        a = 10*(i-1)+1;
        b = 10*i;
        A(a:b,:) = [res(1,i) res(2,i)  1 0        0          0 -data(t).p0(1,i)*res(1,i) -data(t).p0(1,i)*res(2,i)   -data(t).p0(1,i)
                    0        0         0 res(1,i) res(2,i)   1 -data(t).p0(2,i)*res(1,i) -data(t).p0(2,i)*res(2,i)   -data(t).p0(2,i)
                    res(3,i) res(4,i)  1 0        0          0 -data(t).p1(1,i)*res(3,i) -data(t).p1(1,i)*res(4,i)   -data(t).p1(1,i)
                    0        0         0 res(3,i) res(4,i)   1 -data(t).p1(2,i)*res(3,i) -data(t).p1(2,i)*res(4,i)   -data(t).p1(2,i)
                    res(5,i) res(6,i)  1 0        0          0 -data(t).p2(1,i)*res(5,i) -data(t).p2(1,i)*res(6,i)   -data(t).p2(1,i)
                    0        0         0 res(5,i) res(6,i)   1 -data(t).p2(2,i)*res(5,i) -data(t).p2(2,i)*res(6,i)   -data(t).p2(2,i)
                    res(7,i) res(8,i)  1 0        0          0 -data(t).p3(1,i)*res(7,i) -data(t).p3(1,i)*res(8,i)   -data(t).p3(1,i)
                    0        0         0 res(7,i) res(8,i)   1 -data(t).p3(2,i)*res(7,i) -data(t).p3(2,i)*res(8,i)   -data(t).p3(2,i)
                    res(9,i) res(10,i) 1 0        0          0 -data(t).p4(1,i)*res(9,i) -data(t).p4(1,i)*res(10,i)  -data(t).p4(1,i)
                    0        0         0 res(9,i) res(10,i)  1 -data(t).p4(2,i)*res(9,i) -data(t).p4(2,i)*res(10,i)  -data(t).p4(2,i)];

    end
   [~,~,V] = svd(A);
   h = zeros(3,3);
   h(1,:) = V(1:3,9);
   h(2,:) = V(4:6,9);
   h(3,:) = V(7:9,9);

   CamCalib = [311.0520        0        201.8724;
 
                  0         311.3885    113.6210;
 
                  0            0           1   ];
   h = h/h(3,3);

   RRT = pinv(CamCalib)*h;

   n = zeros(3,3);
   n(1:3,1) = RRT(1:3,1);
   n(1:3,2) = RRT(1:3,2);
   n(1:3,3) = cross(RRT(1:3,1),RRT(1:3,2));

   [U,~,V] = svd(n);
R = U*[1 0 0
       0 1 0
       0 0 det(U*V')]*V';
R_c2w = R;
T = RRT(1:3,3)/norm(RRT(1:3,1));


Rot = [cosd(45) -sind(45) 0
       -sind(45) -cosd(45) 0
       0 0 -1];

Trans = [-0.04
         0
         -0.03];

Tcb = zeros(4,4);
Tcb(1:3,1:3) = Rot;
Tcb(1:3,4) = Trans;
Tcb(4,4) = 1;
% Tbc = pinv(Tbc);

Tcw = zeros(4,4);
Tcw(1:3,1:3) = R;
Tcw(1:3,4) = T;
Tcw(4,4) = 1;


Twb = pinv(Tcw)*Tcb;

position = Twb(1:3,4);
orientation = rotm2eul(Twb(1:3,1:3));
end