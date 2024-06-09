function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

ng = [0.014
   -0.0002724
    0.86];
na = [-0.0009
   -0.0017
    0.0004];
nbg = [-0.0021
   -0.0005
   -0.0010];
nba = [0.0001
    0.0012
    0.0002];
% G = [0 -sin(uPrev(6)) cos(uPrev(5))*cos(uPrev(6))
%      0  cos(uPrev(6)) cos(uPrev(5))*sin(uPrev(6))
%      1  0             -sin(uPrev(5))];
invG = [cos(uPrev(6))*sin(uPrev(5))/cos(uPrev(5)) sin(uPrev(5))*sin(uPrev(6))/cos(uPrev(5)) 1
        -sin(uPrev(6)) cos(uPrev(6)) 0
        cos(uPrev(6))/cos(uPrev(5)) sin(uPrev(6))/cos(uPrev(5)) 0];
f = [uPrev(7)
     uPrev(8)
     uPrev(9)
     invG*(angVel-[uPrev(10);uPrev(11);uPrev(12)]-ng)
     eul2rotm([uPrev(6) uPrev(5) uPrev(4)])*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na)-[0;0;9.81]
     nbg
     nba];
uEst = uPrev + dt*f;


At = zeros(15,15);
At(1:3,7:9) = eye(3);
At(4:6,5) = -invG*invG*[0 0 -cos(uPrev(6))*sin(uPrev(5))
                            0 0 -sin(uPrev(6))*sin(uPrev(5))
                            0 0 -cos(uPrev(5))]*(angVel-[uPrev(10);uPrev(11);uPrev(12)]-ng);
At(4:6,6) = -invG*invG*[0 -cos(uPrev(6)) -sin(uPrev(6))*cos(uPrev(5))
                            0 -sin(uPrev(6)) cos(uPrev(6))*cos(uPrev(5))
                            0 0 0]*(angVel-[uPrev(10);uPrev(11);uPrev(12)]-ng);
At(4:6,10:12) = -invG;

R = [cos(uPrev(5))*cos(uPrev(6)), cos(uPrev(6))*sin(uPrev(4))*sin(uPrev(5)) - cos(uPrev(4))*sin(uPrev(6)), sin(uPrev(4))*sin(uPrev(6)) + cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5))
     cos(uPrev(5))*sin(uPrev(6)), cos(uPrev(4))*cos(uPrev(6)) + sin(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)), cos(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)) - cos(uPrev(6))*sin(uPrev(4))
      -sin(uPrev(5)),                        cos(uPrev(5))*sin(uPrev(4)),                        cos(uPrev(4))*cos(uPrev(5))];

dRdx = [0, sin(uPrev(4))*sin(uPrev(6)) + cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5)),   cos(uPrev(4))*sin(uPrev(6)) - cos(uPrev(6))*sin(uPrev(4))*sin(uPrev(5))
        0, cos(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)) - cos(uPrev(6))*sin(uPrev(4)), - cos(uPrev(4))*cos(uPrev(6)) - sin(uPrev(4))*sin(uPrev(5))*sin(uPrev(6))
        0,                        cos(uPrev(4))*cos(uPrev(5)),                         -cos(uPrev(5))*sin(uPrev(4))];

dRdy = [-cos(uPrev(6))*sin(uPrev(5)), cos(uPrev(5))*cos(uPrev(6))*sin(uPrev(4)), cos(uPrev(4))*cos(uPrev(5))*cos(uPrev(6))
        -sin(uPrev(5))*sin(uPrev(6)), cos(uPrev(5))*sin(uPrev(4))*sin(uPrev(6)), cos(uPrev(4))*cos(uPrev(5))*sin(uPrev(6))
        -cos(uPrev(5)),       -sin(uPrev(4))*sin(uPrev(5)),       -cos(uPrev(4))*sin(uPrev(5))];

dRdz = [-cos(uPrev(5))*sin(uPrev(6)), - cos(uPrev(4))*cos(uPrev(6)) - sin(uPrev(4))*sin(uPrev(5))*sin(uPrev(6)), cos(uPrev(6))*sin(uPrev(4)) - cos(uPrev(4))*sin(uPrev(5))*sin(uPrev(6))
         cos(uPrev(5))*cos(uPrev(6)),   cos(uPrev(6))*sin(uPrev(4))*sin(uPrev(5)) - cos(uPrev(4))*sin(uPrev(6)), sin(uPrev(4))*sin(uPrev(6)) + cos(uPrev(4))*cos(uPrev(6))*sin(uPrev(5))
                      0,                                      0,                                    0]; 

At(7:9,4) = dRdx*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na);
At(7:9,5) = dRdy*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na);
At(7:9,6) = dRdz*(acc-[uPrev(13);uPrev(14);uPrev(15)]-na);
At(7:9,13:15) = -R;

Vt = zeros(15,12);
Vt(4:6,1:3) = -invG;
Vt(7:9,4:6) = -R;
Vt(10:12,7:9) = eye(3);
Vt(13:15,10:12) = eye(3);


Q = eye(12);
Qd = Q*dt;

Ft = eye(15) + dt*At;
covarEst = Ft*covarPrev*Ft' + Vt*Qd*Vt';

end

