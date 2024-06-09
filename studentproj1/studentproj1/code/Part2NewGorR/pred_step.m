function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

ng = normrnd(0,0.001,[3,1]);
na = normrnd(0,0.001,[3,1]);
nbg = normrnd(0,0.001,[3,1]);
nba = normrnd(0,0.001,[3,1]);
% G = [0 -sin(uPrev(6)) cos(uPrev(5))*cos(uPrev(6))
%      0  cos(uPrev(6)) cos(uPrev(5))*sin(uPrev(6))
%      1  0             -sin(uPrev(5))];
invG = [1 sin(uPrev(4))*sin(uPrev(5))/cos(uPrev(5)) -cos(uPrev(4))*sin(uPrev(5))/cos(uPrev(5)) 
        0 cos(uPrev(4)) sin(uPrev(4))
        0 -sin(uPrev(4))/cos(uPrev(5)) cos(uPrev(4))/cos(uPrev(5)) ];
f = [uPrev(7)
     uPrev(8)
     uPrev(9)
     invG*(angVel-[uPrev(10);uPrev(11);uPrev(12)])
     eul2rotm([uPrev(4) uPrev(5) uPrev(6)],'XYZ')*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na)-[0;0;9.81]
     nbg
     nba];
uEst = uPrev + dt*f;


At = zeros(15,15);
At(1:3,7:9) = eye(3);
At(4:6,5) = -invG*invG*[0, 0,         cos(uPrev(5))
0, 0,  sin(uPrev(4))*sin(uPrev(5))
0, 0, -cos(uPrev(4))*sin(uPrev(5))]*(angVel-[uPrev(10);uPrev(11);uPrev(12)]-ng);
At(4:6,4) = -invG*invG*[0,       0,              0
0, -sin(uPrev(4)), -cos(uPrev(4))*cos(uPrev(5))
0,  cos(uPrev(4)), -cos(uPrev(5))*sin(uPrev(4))]*(angVel-[uPrev(10);uPrev(11);uPrev(12)]-ng);
At(4:6,10:12) = -invG;

R = [                       cos(uPrev(5))*cos(uPrev(6)),                       -cos(uPrev(5))*sin(uPrev(6)),         sin(uPrev(5))
cos(uPrev(4))*sin(uPrev(6)) + cos(uPrev(6))*sin(uPrev(4))*sin(uPrev(5)), cos(uPrev(4))*cos(uPrev(6)) - sin(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)), -cos(uPrev(5))*sin(uPrev(4))
sin(uPrev(4))*sin(uPrev(6)) - cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5)), cos(uPrev(6))*sin(uPrev(4)) + cos(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)),  cos(uPrev(4))*cos(uPrev(5))];

dRdx = [                                   0,                                      0,              0
cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5)) - sin(uPrev(4))*sin(uPrev(6)), - cos(uPrev(6))*sin(uPrev(4)) - cos(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)), -cos(uPrev(4))*cos(uPrev(5))
cos(uPrev(4))*sin(uPrev(6)) + cos(uPrev(6))*sin(uPrev(4))*sin(uPrev(5)),   cos(uPrev(4))*cos(uPrev(6)) - sin(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)), -cos(uPrev(5))*sin(uPrev(4))];

dRdy = [       -cos(uPrev(6))*sin(uPrev(5)),         sin(uPrev(5))*sin(uPrev(6)),         cos(uPrev(5))
 cos(uPrev(5))*cos(uPrev(6))*sin(uPrev(4)), -cos(uPrev(5))*sin(uPrev(4))*sin(uPrev(6)),  sin(uPrev(4))*sin(uPrev(5))
-cos(uPrev(4))*cos(uPrev(5))*cos(uPrev(6)),  cos(uPrev(4))*cos(uPrev(5))*sin(uPrev(6)), -cos(uPrev(4))*sin(uPrev(5))];

dRdz = [                      -cos(uPrev(5))*sin(uPrev(6)),                         -cos(uPrev(5))*cos(uPrev(6)), 0
cos(uPrev(4))*cos(uPrev(6)) - sin(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)), - cos(uPrev(4))*sin(uPrev(6)) - cos(uPrev(6))*sin(uPrev(4))*sin(uPrev(5)), 0
cos(uPrev(6))*sin(uPrev(4)) + cos(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)),   cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5)) - sin(uPrev(4))*sin(uPrev(6)), 0]; 

At(7:9,4) = dRdx*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na);
At(7:9,5) = dRdy*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na);
At(7:9,6) = dRdz*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na);
At(7:9,13:15) = -R;

Vt = zeros(15,12);
Vt(4:6,1:3) = -invG;
Vt(7:9,4:6) = -R;
Vt(10:12,7:9) = eye(3);
Vt(13:15,10:12) = eye(3);

% Q = [ng;na;nbg;nba]*[ng',na',nbg',nba'];
Q = eye(12);
Qd = Q*dt;

Ft = eye(15) + dt*At;
covarEst = Ft*covarPrev*Ft' + Vt*Qd*Vt';

end

