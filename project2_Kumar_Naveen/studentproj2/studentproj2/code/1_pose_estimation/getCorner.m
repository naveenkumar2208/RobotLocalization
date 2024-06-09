function res = getCorner(id)
%% CHANGE THE NAME OF THE FUNCTION TO getCorner
    %% Input Parameter Description
    % id = List of all the AprilTag ids detected in the current image(data)
    res = zeros(10,length(id));
    %% Output Parameter Description
    % res = List of the coordinates of the 4 corners (or 4 corners and the
    % centre) of each detected AprilTag in the image in a systematic method
    for i = 1:length(id)
        if id(i) >= 0 && id(i) <=35
            c(1) =  2*rem(id(i),12)*0.152 + (0.152/2); % Center % 1 means x 
            c(2) =  0.152/2 + 2*0.152*floor(id(i)/12);                         % 2 means y
            bl(1) = 0.152 + 2*rem(id(i),12)*0.152; %Bottom Left Corner
            bl(2) = 2*0.152*floor(id(i)/12);
            br(1) = 0.152 + 2*rem(id(i),12)*0.152; %Bottom Right Corner
            br(2) = 0.152 + 2*0.152*floor(id(i)/12);
            tr(1) = 2*0.152*rem(id(i),12); % Top Right Corner
            tr(2) = 0.152 + 2*0.152*floor(id(i)/12);
            tl(1) = 2*0.152*rem(id(i),12); % Top Left Corner
            tl(2) = 0 + 2*0.152*floor(id(i)/12);
              
        elseif id(i) >= 36 && id(i) <= 71
            c(1) =  2*rem(id(i),12)*0.152 + (0.152/2); % Center % 1 means x 
            c(2) =  0.152/2 + 2*0.152*floor(id(i)/12) + (0.178-0.152);                         % 2 means y
            bl(1) = 0.152 + 2*rem(id(i),12)*0.152; %Bottom Left Corner
            bl(2) = 2*0.152*floor(id(i)/12) + (0.178-0.152);
            br(1) = 0.152 + 2*rem(id(i),12)*0.152; %Bottom Right Corner
            br(2) = 0.152 + 2*0.152*floor(id(i)/12) + (0.178-0.152);
            tr(1) = 2*0.152*rem(id(i),12); % Top Right Corner
            tr(2) = 0.152 + 2*0.152*floor(id(i)/12) + (0.178-0.152);
            tl(1) = 2*0.152*rem(id(i),12); % Top Left Corner
            tl(2) = 0 + 2*0.152*floor(id(i)/12) + (0.178-0.152);

        elseif id(i) >= 72 && id(i) <= 107
            c(1) =  2*rem(id(i),12)*0.152 + (0.152/2); % Center % 1 means x 
            c(2) =  0.152/2 + 2*0.152*floor(id(i)/12) + 2*(0.178-0.152);                         % 2 means y
            bl(1) = 0.152 + 2*rem(id(i),12)*0.152; %Bottom Left Corner
            bl(2) = 2*0.152*floor(id(i)/12) + 2*(0.178-0.152);
            br(1) = 0.152 + 2*rem(id(i),12)*0.152; %Bottom Right Corner
            br(2) = 0.152 + 2*0.152*floor(id(i)/12) + 2*(0.178-0.152);
            tr(1) = 2*0.152*rem(id(i),12); % Top Right Corner
            tr(2) = 0.152 + 2*0.152*floor(id(i)/12) + 2*(0.178-0.152);
            tl(1) = 2*0.152*rem(id(i),12) + (0.178-0.152); % Top Left Corner
            tl(2) = 0 + 2*0.152*floor(id(i)/12) + 2*(0.178-0.152);
        end
        res(:,i) = [c';bl';br';tr';tl'];
    end
end