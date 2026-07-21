% ------------------------------------------------------------------------
% Author: Sophie A. Liu
% Purpose: preliminary compartment modeling
% a bit irritating to hard code parameters
% ------------------------------------------------------------------------

% --- assumptions -------
% conservation of volume in each compartment
% homeostasis will be reached
% once labeled, it cannot be unlabeled 
% (however, label only lasts 5 min and we have 10 min increments...
% ------------------------

function dydt = splitComps(t,y,p)
    b = y(1);

% --- continues down similarly, many vars ----
    l1 = y(2);        % bone marrow
    l2a = y(3);       % both red and white pulp spleen
    l2b = y(4);

    l3a = y(5);      % lymph axillary
    l3b = y(6);      % lymph iliac
    l3c = y(7);      % lymph inguinal
    l3d = y(8);      % lymph mesenteric

    l4 = y(9);       % peyer's patches

% -----------------
    u1 = y(10);    
    u2a = y(11);
    u2b = y(12);

    u3a = y(13);
    u3b = y(14);
    u3c = y(15);
    u3d = y(16);

    u4 = y(17);


% --- reverse rates ---
    kneg1 = p(1);             % bm   
    kneg2a = p(2);            % spleen   
    kneg2b = p(3);

    kneg3a = p(4);            % LNs     
    kneg3b = p(5);                
    kneg3c = p(6);                
    kneg3d = p(7);

    kneg4 = p(8);             % PP

% --- forward rates ---
    k1 = p(9);

    k2a = p(10); 
    k2b = p(11);

    k3a = p(12);
    k3b = p(13);
    k3c = p(14);
    k3d = p(15);

    k4 = p(16); 


% --- system -------------
dydt = [ 
    % --- all exiting the blood ---
    -(k1+ k2a+k2b+ k3a+k3b+k3c+k3d+ k4)*b + ...
        kneg1*(l1+u1) + ...
        kneg2a*(l2a+u2a) + kneg2b*(l2b+u2b) + ...
        kneg3a*(l3a+u3a) + kneg3b*(l3b+u3b) + kneg3c*(l3c+u3c) + kneg3d*(l3d+u3d) + ...
        kneg4*(l4 + u4); 
   
    % --- total labeled flux for each tissue compartment ---
        k1*b - kneg1*l1;
        k2a*b - kneg2a*l2a;
        k2b*b - kneg2b*l2b;

        k3a*b - kneg3a*l3a;
        k3b*b - kneg3b*l3b;
        k3c*b - kneg3c*l3c;
        k3d*b - kneg3d*l3d;

        k4*b - kneg4*l4;

    % --- unlabeled exiting each tissue compartment --------
        -kneg1*u1;
        -kneg2a*u2a;
        -kneg2b*u2b;
        -kneg3a*u3a;
        -kneg3b*u3b;
        -kneg3c*u3c;
        -kneg3d*u3d;
        -kneg4*u4;];

end

% --- initial guess for params based on R. magnitude not direction 
p = [0.002370606 , ...                   % bm reverse
    0.02488078, 0.001255718, ...         % spleen red, white
    0.0009416563,0.001221415, 0.0008115359, 0.000860843, ...   % LNs: ax, il, ing, mes               
    0.001053737 ,                        % peyer's

    0.003289281, ...                     % bm forward
    0.01185686 , 0.02052329  , ...       % spleen red, white
    0.02520642, 0.005542072 , 0.01259731 , 0.03447049 , ...   % LNs: ax, il, ing, mes
    0.008135678 ];                       % peyer's

% fraction ICs are less skewed ---
y0 = [1, ...              % blood
    0, ...                % bm
    0, 0, ...             % spleen
    0, 0, 0, 0, ...       % LNs
    0, ...                % PPs
    0, ...
    1, 1, ... 
    1, 1,1,1, ...
    1];      

% --- they can only take the experiment out to 3 days (72 hrs, 4320 sec.)
tspan = [0, 600];    % rates etc. in 1/min unit

% --- solving
[t, y] = ode45(@(t, y)splitComps(t, y, p), tspan, y0);
% plot(t,y)
% legend('blood', ...
%     'labeled bone marrow','labeled spleen', ...
%     'labeled lymph axillary', 'L LN iliac', 'L LN inguinal', 'L LN mesenteric',...
%     'labeled peyer', ...
%     'UL BM','UL spleen', ...
%     'UL LN axillary', 'UL LN iliac', 'UL LN inguinal', 'UL LN mesenteric',...
%     'UL PPs')
% 
% UL = sum(y(:,10:17),2);
% 
% figure
% plot(t,UL,'LineWidth',2)
% ylabel('Total unlabeled')
% xlabel('Time')


% --- labeled will approach completion as t -> inf ---
N_BM  = y(:,2);
N_SRP = y(:,3);
N_SWP = y(:,4);
N_AX  = y(:,5);
N_IL  = y(:,6);
N_IG  = y(:,7);
N_ME  = y(:,8);
N_PP  = y(:,9);

plot(t,[N_BM N_SRP, N_SWP, N_AX N_IL N_IG N_ME N_PP])
legend('BM','Red pulp', 'White pulp','Ax','Iliac','Inguinal','Mes','PP')
