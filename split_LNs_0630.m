function dydt = splitComps(t,y,p)
b = y(1);

% -- continues down similarly, many vars ----
l1 = y(2);       % bone marrow
l2 = y(3);       % both red and white pulp spleen

l3a = y(4);      % lymph axillary
l3b = y(5);      % lymph iliac
l3c = y(6);      % lymph inguinal
l3d = y(7);      % lymph mesenteric

l4 = y(8);       % peyer's patches

% --
u1 = y(9);    
u2 = y(10);

u3a = y(11);
u3b = y(12);
u3c = y(13);
u3d = y(14);

u4 = y(15);


% ----------------------------

kneg1 = p(1);                % reverse direction ksub-1
kneg2 = p(2);                % forward direction k1

kneg3a = p(3);                
kneg3b = p(4);                
kneg3c = p(5);                
kneg3d = p(6);

kneg4 = p(7); 

k1 = p(8);                
k2 = p(9); 

k3a = p(10);
k3b = p(11);
k3c = p(12);
k3d = p(13);

k4 = p(14); 


dydt = [-(k1+k2+ k3a+k3b+k3c+k3d+ k4)*b + kneg1*(l1+u1) + kneg2*(l2+u2) + ...
    kneg3a*(l3a+u3a) + kneg3b*(l3b+u3b) + kneg3c*(l3c+u3c) + + kneg3d*(l3d+u3d) + ...
    kneg4*(l4 + u4); 

        kneg1*u1;
        k2*b - kneg2*l2;

        k3a*b - kneg3a*l3a;
        k3b*b - kneg3b*l3b;
        k3c*b - kneg3c*l3c;
        k3d*b - kneg3d*l3d + k4*l4;

        k4*b - kneg4*l4;

        -kneg1*u1;
        -kneg2*u2;

        -kneg3a*u3a;
        -kneg3b*u3b;
        -kneg3c*u3c;
        -kneg3d*u3d;

        -kneg4*u4];

end

% initial guess for params. magnitude not direction
% could fit a least squares estimate but not sure what that entails
p = [0.004, 0.022, ...
    0.0006, 0.0014, 0.00049, 0.00063, ...
    0.0012, ...
    % now k forward
    0.0022, 0.021, ...
    0.0055, 0.0025, 0.0026, 0.014, ...
    0.0042];      

% ICs. check matrix carefully
% y0 = [34969, ... 
%     0, 0, 0, 0, 0, 0, 0, 0, ... 
%     19086, 277000, ...
%     301675, 63998, 186545, ...
%     128766];

% fracs easier and less skewed
y0 = [1, ... 
    0, 0, 0, 0, 0, 0, 0, ... 
    1, 1, ... 
    1, 1,1,1, ...
    1];      
tspan = [0, 5000];

% solving
[t, y] = ode45(@(t, y)splitComps(t, y, p), tspan, y0);
plot(t,y)
legend('blood', ...
    'labeled bone marrow','labeled spleen', ...
    'labeled lymph axillary', 'L LN iliac', 'L LN inguinal', 'L LN mesenteric',...
    'labeled peyer', ...
    'UL BM','UL spleen', ...
    'UL LN axillary', 'UL LN iliac', 'UL LN inguinal', 'UL LN mesenteric',...
    'UL PPs')

UL = sum(y(:,9:15),2);

figure
plot(t,UL,'LineWidth',2)
ylabel('Total unlabeled')
xlabel('Time')

N_BM  = y(:,2);
N_SPL = y(:,3);
N_AX  = y(:,4);
N_IL  = y(:,5);
N_IG  = y(:,6);
N_ME  = y(:,7);
N_PP  = y(:,8);

plot(t,[N_BM N_SPL N_AX N_IL N_IG N_ME N_PP])
legend('BM','Spleen','Ax','Iliac','Inguinal','Mes','PP')