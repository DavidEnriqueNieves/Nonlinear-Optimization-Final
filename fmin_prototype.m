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
x0 = [P2G, P5G, D2G, D5G, I2G, I5G];
lb = [0, 0, 40, 20, 0, 0];
ub = [8, 8, 200, 60, 1, 1];
% [x, fval, exitflag, output] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub, @(x) clip(x), options);

% [x, fval, exitflag, output] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub, @(x) identity(x) , options);

% THIS WORKS, AT LEAST
% [x, fval, exitflag, output] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub);

%options = optimoptions('fmincon', 'MaxIterations', 1000);
% options = optimoptions('fmincon', 'MaxIterations', 1000, 'Algorithm' , 'sqp');
% options = optimoptions('fmincon', 'MaxIterations', 1000, 'Algorithm' , 'interior-point');
options = optimoptions('fmincon', 'MaxIterations', 5000, 'Algorithm' , 'active-set', "Display", "iter");
[x, fval, exitflag, output, lambda, grad, hessian] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub, @(x) nothing(x), options);

% Display results
fprintf('Optimal Values:\n');
fprintf('P2G: %.2f\n', x(1));
fprintf('P5G: %.2f\n', x(2));
fprintf('D2G: %.2f\n', x(3));
fprintf('D5G: %.2f\n', x(4));
fprintf('I2G: %.2f\n', x(5));
fprintf('I5G: %.2f\n', x(6));


function [c,ceq] = nothing(x)
    c(1) = 0;
    c(2) = 0;
    ceq = [];
end


function id = identity(x)
        id = x;
end


% Define packet_loss function
function pl = packet_loss(distance)
    pl = 0.10 * distance;
end

% Define speed function (replace with actual implementation)
function s = average_speed(c, i,distance, freq)

  % Calculate the transmission speed of a phone given the chunk size, idle time, distance and frequency
  
  % The 2.4GHz connection starts at 150MB/s at 40m, it is assumed to decrease logarithmically to 50Mb/s at 200m.
  % For the 5GHz connection, the speed starts at 1500MB/s at 20m, it is assumed to decrease logarithmically to 250Mb/s at 60m.
  
  % The chunk size is between 512 and 4098 bytes, the idle time is between 0 and 1 second (0: completely idle, 1: 100% busy)
  % The speed increases linearly with respect to the chunk size and the idle time.
  
  if freq == 2.4
      start_speed = 250
      end_speed = 50
      start_distance = 40
      end_distance = 200
  elseif  freq == 5
      start_speed = 1500
      end_speed = 250
      start_distance = 20
      end_distance = 60
  else
	  error("Frequency must be 2.4 or 5")
  end
  
  % Scale the distance to be between 0 and 1
  fin_distance = (distance - start_distance) / (end_distance - start_distance)
  
  % Make it negative exponential
  s = end_speed + (1 -exp(-5) * (exp(5*fin_distance))) * (start_speed - end_speed)
  %    speed = speed * (chunk_size / 4096) * (idle_time)
end

% Define objective function
function y = f(P2G, P5G, D2G, D5G, I2G, I5G)
    max_power = 8;
    min_2G_path_loss = 3.60435984234;
    min_5G_path_loss = 3.6220886093;
    max_5G_path_loss = 4.09920986402;
    max_2G_path_loss = 4.30332984668;

    G_over_c = 10^9 / 299792458;
    C2G = 1024;
    C5G = 2048;
    total_distance = D2G*P2G / (200 * 8) + D5G*P5G / (60 * 8);
    number_of_phones = P2G*P5G / (8 * 8);
    free_space_path_loss =  (  log(4*pi*2.4*G_over_c*D2G) - min_2G_path_loss ) / ( max_2G_path_loss - min_2G_path_loss )  
    		         +  (  log(4*pi*5*G_over_c*D5G) - min_5G_path_loss ) / ( max_5G_path_loss - min_5G_path_loss );

    % bytes_lost = C2G * packet_loss(D2G) + C5G * packet_loss(D5G);
    bytes_lost = 0

    max_2G_speed = 250;
    max_5G_speed = 1500;

    min_2G_speed = 50;
    min_5G_speed = 250;
    

    transfer_speed = ( average_speed(C2G,I2G,D2G, 2.4) - min_2G_speed ) / (max_2G_speed - min_2G_speed) + ( average_speed(C5G,I5G,D5G,  5) - min_5G_speed ) / (max_5G_speed - min_5G_speed);

    
    power_savings = I2G*P2G / max_power + I5G*P5G / max_power;

    power_savings_coeff = 1;
    speed_coeff = 1;
    path_loss_coeff = 1;
    bytes_lost_coeff = 1;
    number_of_phones_coeff = 1;
    total_distance_coeff = 1;
y = power_savings_coeff * power_savings + ...
    speed_coeff * transfer_speed + ...
    path_loss_coeff * free_space_path_loss + ...
    bytes_lost_coeff * bytes_lost + ...
    number_of_phones_coeff * number_of_phones + ...
    total_distance_coeff * total_distance;
end

% Define clip function
function x = clip(x)
    x(1) = min(max(x(1), 0), 8);
    x(2) = min(max(x(2), 0), 8);
    x(3) = min(max(x(3), 40), 200);
    x(4) = min(max(x(4), 20), 60);
    x(5) = min(max(x(5), 0), 1);
    x(6) = min(max(x(6), 0), 1);
    x = x;
end

