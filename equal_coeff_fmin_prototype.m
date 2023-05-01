% Constants
% Define initial input values
G_over_c = 10000 / 3;
% between 512 and 4096
% These are NOT optimizable
C2G = 1024;
C5G = 2048;

D2G = 100;
D5G = 30;
I2G = 0.5;
I5G = 0.8;
P2G = 4;
P5G = 3;

% Use MATLAB's optimization functions to find optimal values

x0 = [ D2G, D5G, I2G, I5G, P2G, P5G];
lb = [ 40, 20,  0.5, 0.5,  0, 0];
ub = [ 200, 60, 1, 1, 8, 8];
f_coeffs = [ 1, 1, 1, 1, 1, 1];
% [x, fval, exitflag, output] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub, @(x) clip(x), options);

% [x, fval, exitflag, output] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub, @(x) identity(x) , options);

% THIS WORKS, AT LEAST
% [x, fval, exitflag, output] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6)), x0, [], [], [], [], lb, ub);

%options = optimoptions('fmincon', 'MaxIterations', 1000);
% options = optimoptions('fmincon', 'MaxIterations', 1000, 'Algorithm' , 'sqp');
% options = optimoptions('fmincon', 'MaxIterations', 1000, 'Algorithm' , 'interior-point');
options = optimoptions('fmincon', 'MaxIterations', 5000, 'Algorithm' , 'active-set', "Display", "iter");
[x, fval, exitflag, output, lambda, grad, hessian] = fmincon(@(x) -f(x(1), x(2), x(3), x(4), x(5), x(6), f_coeffs), x0, [], [], [], [], lb, ub, @(x) nothing(x), options);

% Display results
fprintf('Optimal Values:\n');
fprintf('D2G: %.2f\n', x(1));
fprintf('D5G: %.2f\n', x(2));
fprintf('I2G: %.2f\n', x(3));
fprintf('I5G: %.2f\n', x(4));
fprintf('P2G: %.2f\n', x(5));
fprintf('P5G: %.2f\n', x(6));

suffix = "";
file_name = strcat("~/Documents/Semester10/MTH5335/Project/results(algorithm="+ options.Algorithm+ "+iterations="+ options.MaxIterations+  "f=" + fval+ ")" + suffix + ".csv" );


write_results_to_csv(file_name, x, fval, x0, ub, lb, exitflag, output, lambda, grad, hessian, C2G, C5G, f_coeffs);
	
	% Define function to write results to a CSV file
	function write_results_to_csv(file_name, x_opt, fval, x0, ub, lb, exitflag, output, lambda, grad, hessian, C2G, C5G, f_coeffs)
	
	% Create cell array to store results
	results = cell(21, 8);
	
	results{1, 1} = "TABLE";
	% Fill in results cell array
	results{2, 1} = "D2G";
	results{3, 1} = "D5G";
	results{4, 1} = "I2G";
	results{5, 1} = "I5G";
	results{6, 1} = "P2G";
	results{7, 1} = "P5G";
	
	results{1, 2} = 'x'
	results{2, 2} = x_opt(1);
	results{3, 2} = x_opt(2);
	results{4, 2} = x_opt(3);
	results{5, 2} = x_opt(4);
	results{6, 2} = x_opt(5);
	results{7, 2} = x_opt(6);
	
	results{1, 3} = 'x_0'
	results{2, 3} = x0(1);
	results{3, 3} = x0(2);
	results{4, 3} = x0(3);
	results{5, 3} = x0(4);
	results{6, 3} = x0(5);
	results{7, 3} = x0(6);
	
	results{1, 4} = 'UB'
	results{2, 4} = ub(1);
	results{3, 4} = ub(2);
	results{4, 4} = ub(3);
	results{5, 4} = ub(4);
	results{6, 4} = ub(5);
	results{7, 4} = ub(6);
	
	results{1, 5} = 'LB'
	results{2, 5} = lb(1);
	results{3, 5} = lb(2);
	results{4, 5} = lb(3);
	results{5, 5} = lb(4);
	results{6, 5} = lb(5);
	results{7, 5} = lb(6);
	
	
	results{9, 1} = 'F-value';
	results{10, 1} = fval;
	
	results{9, 2} = 'Iterations';
	results{10, 2} = output.iterations;
	
	results{9, 3} = 'Algorithm';
	results{10, 3} = output.algorithm;
	
	results{9, 4} = 'FuncCount';
	results{10, 4} = output.funcCount;
	
	results{9, 5} = 'StepSize';
	results{10, 5} = output.stepsize;
	
	
	% if eigenvalues of hessian are greater than or equal to 0
	
	results{9, 6} = "H >= 0";
	results{10, 6} = string(all(eig(hessian) >= -0.001));
	results{9, 7} = "H > 0 ";
	results{10, 7} = string(all(eig(hessian) >= 0));
	
	results{12, 2} = "Chunk Size 2G";
	results{13, 2} = C2G;
	
	results{12, 3} = "Chunk Size 5G";
	results{13, 3} = C5G;
	
	results{14, 1} = "Grad";
	results{15, 1} = grad(1);
	results{16, 1} = grad(2);
	results{17, 1} = grad(3);
	results{18, 1} = grad(4);
	results{19, 1} = grad(5);
	results{20, 1} = grad(6);
	
	results{1, 5} = "$f(x)$ Coeffs";
	results{2, 5} = f_coeffs(1);
	results{3, 5} = f_coeffs(2);
	results{4, 5} = f_coeffs(3);
	results{5, 5} = f_coeffs(4);
	results{6, 5} = f_coeffs(5);
	results{7, 5} = f_coeffs(6);
	% Write results to CSV file
	% writematrix( results,file_name);
	cell2csv(file_name,results)
	% csvwrite(file_name, results);
	
end

function cell2csv(filename,cellArray,delimiter)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(filename,cellArray,delimiter)
%
% filename      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% delimiter = seperating sign, normally:',' (default)
%
% by Sylvain Fiedler, KA, 2004
% modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
if nargin<3
    delimiter = ',';
end

datei = fopen(filename,'w');
for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)

        var = eval(['cellArray{z,s}']);

        if size(var,1) == 0
            var = '';
        end

        if isnumeric(var) == 1
            var = num2str(var);
        end

        fprintf(datei,var);

        if s ~= size(cellArray,2)
            fprintf(datei,[delimiter]);
        end
    end
    fprintf(datei,'\n');
end
fclose(datei);

end

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
  s = ( end_speed + (1 -exp(-5) * (exp(5*fin_distance))) * (start_speed - end_speed) ) * i 
  %    speed = speed * (chunk_size / 4096) * (idle_time)
end

% Define objective function
% function y = f(P2G, P5G, D2G, D5G, I2G, I5G, f_coeffs)
function y = f( D2G, D5G, I2G, I5G,P2G, P5G, f_coeffs)
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
    
    disp("D2G");
    disp(D2G);
    disp("D5G");
    disp(D5G);

    disp(" max_5G_path_loss - min_5G_path_loss ");
    disp(max_5G_path_loss - min_5G_path_loss);

    disp(" max_2G_path_loss - min_2G_path_loss ");
    disp(max_2G_path_loss - min_2G_path_loss);

    disp(" log(4*pi*2.4*G_over_c*D2G) - min_2G_path_loss ");
    disp(log(4*pi*2.4*G_over_c*D2G) - min_2G_path_loss);

    disp(" log(4*pi*5*G_over_c*D5G) - min_5G_path_loss ");
    disp(log(4*pi*5*G_over_c*D5G) - min_5G_path_loss);

    free_space_path_loss =  (  log10(4*pi*2.4*G_over_c*D2G) - min_2G_path_loss ) ./ ( max_2G_path_loss - min_2G_path_loss )  
    		         +  (  log10(4*pi*5*G_over_c*D5G) - min_5G_path_loss ) ./ ( max_5G_path_loss - min_5G_path_loss );


    % bytes_lost = C2G * packet_loss(D2G) + C5G * packet_loss(D5G);
    bytes_lost = 0

    max_2G_speed = 250;
    max_5G_speed = 1500;

    min_2G_speed = 50;
    min_5G_speed = 250;
    

    transfer_speed = ( average_speed(C2G,I2G,D2G, 2.4) - min_2G_speed ) / (max_2G_speed - min_2G_speed) + ( average_speed(C5G,I5G,D5G,  5) - min_5G_speed ) / (max_5G_speed - min_5G_speed);

    
    % disp("I2G");
    % disp(I2G);
    % disp("I5G");
    % disp(I5G);
    % disp("P2G");
    % disp(P2G);
    % disp("P5G");
    % disp(P5G);
    % disp("max_power");
    % disp(max_power);

    power_savings = ( (I2G*P2G) / max_power ) + ( (I5G*P5G) / max_power );

    power_savings_coeff = f_coeffs(1);
    speed_coeff = f_coeffs(2);
    path_loss_coeff = f_coeffs(3);
    bytes_lost_coeff = f_coeffs(4);
    number_of_phones_coeff = f_coeffs(5);
    total_distance_coeff = f_coeffs(6);

    disp("total_distance")
    disp(total_distance)
    disp("number_of_phones")
    disp(number_of_phones)
    disp("free_space_path_loss")
    disp(free_space_path_loss)
    disp("bytes_lost")
    disp(bytes_lost)
    disp("transfer_speed")
    disp(transfer_speed)
    disp("power_savings")
    disp(power_savings)
y = -1 * power_savings_coeff * power_savings + ...
  speed_coeff * transfer_speed - ...
  path_loss_coeff * free_space_path_loss + ...
  bytes_lost_coeff * bytes_lost + ...
  number_of_phones_coeff * number_of_phones + ...
  total_distance_coeff * total_distance;
  disp("y");
  disp(y);
end

