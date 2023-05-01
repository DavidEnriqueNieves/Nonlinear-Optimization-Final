% Constants
G_over_c = 10000 / 3;
C2G = 1024;
C5G = 2048;


% Define initial input values
P2G = 4;
P5G = 3;
D2G = 100;
D5G = 30;
I2G = 0.5;
I5G = 0.8;

% Use MATLAB's optimization functions to find optimal values
x0 = [P2G, P5G, D2G, D5G]

% Define packet_loss function
function pl = packet_loss(distance)
    pl = 0.10 * distance;
end

% Define speed function (replace with actual implementation)
function s = speed(c, i)
    s = 0; % Replace with actual implementation
end

% Define objective function
function y = f(P2G, P5G, D2G, D5G, I2G, I5G)
    total_distance = D2G*P2G + D5G*P5G;
    number_of_phones = P2G*P5G;
    free_space_path_loss = log(4*pi*2.4*G_over_c*D2G) + log(4*pi*5*G_over_c*D5G);
    bytes_lost = C2G * packet_loss(D2G) + C5G * packet_loss(D5G);
    transfer_speed = speed(C2G,I2G) + speed(C5G,I5G);
    power_savings = I2G*P2G + I5G*P5G;
    
    y = total_distance - number_of_phones - free_space_path_loss - bytes_lost + transfer_speed + power_savings;
end

% Define clip function
function x = clip(x)
    x(1) = min(max(x(1), 0), 8);
    x(2) = min(max(x(2), 0), 8);
    x(3) = min(max(x(3), 40), 200);
    x(4) = min(max(x(4), 20), 60);
    x(5) = min(max(x(5), 0), 1);
    x(6) = min(max(x(6), 0), 1);
end

