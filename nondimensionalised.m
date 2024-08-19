tic

% THE EQUATIONS:
% C1' = C1*(r1(M) - a11*C1 - a12*C2)
% C2' = C2*(r2(M) - a21*C1 - a22*C2)
% M' = i - M - p*M*(C1 + C2)



% THINGS TO DO WHEN RUNNING A NEW BIFURCATION DIAGRAM
% 1 - Change your equations for r1 and r2
% 2 - Change your equations for r1' and r2'
% 3 - Set your parameters and hide the ones you want to plot
% 4 - Choose the length and width of your plot
% 5 - Choose the number of points to plot
% 6 - Place the 2 parameters you want to plot into the for loops
% 
% These instructions are noted with numbers throughout the code



% 1a - ALTER YOUR EQUATION FOR r1(M)
function r1result = r1(M)
    r1result = sech(M-4) - 0.2;
end


% 1b - ALTER YOUR EQUATION FOR r2(M)
function r2result = r2(M)
    r2result = sech(M-6) - 0.2;
end


% 2a - ALTER YOUR EQUATION FOR r1'(M)
function r1dashresult = r1dash(M)
    r1dashresult = -tanh(M-4)*sech(M-4);
end


% 2b - ALTER YOUR EQUATION FOR r2'(M)
function r2dashresult = r2dash(M)
    r2dashresult = -tanh(M-6)*sech(M-6);
end



% Jacobian
function matrix = Jacobian(C1,C2,M,A11,A12,A21,A22,P)
    matrix = [r1(M)-2*A11*C1-A12*C2, -A12*C1, C1*r1dash(M);
        -A21*C2, r2(M)-A21*C1-2*A22*C2, C2*r2dash(M);
        -P*M, -P*M, -1-P*C1-P*C2];
end

% To solve M1 Survival
function m1 = m1solve(M,A11,I,P)
    m1 = M*(1+(r1(M)*P/A11)) - I;
end

% To solve M2 Survival
function m2 = m2solve(M,A22,I,P)
    m2 = M*(1+(r2(M)*P/A22)) - I;
end

% To solve M3 (Coexistence)
function m3 = m3solve(M,A11,A12,A21,A22,I,P)
    m3 = M*(1 + (r1(M)*(A22-A21) + r2(M)*(A11-A12))*P/(A11*A22 - A12*A21)) - I;
end

tolerance = 1e-5;



% 3a - CHOOSE YOUR PARAMETERS
% 3b - HIDE THE TWO YOU WANT TO PLOT
A11 = 1;
A12 = 0;
A21 = 0;
A22 = 5.2;
%I = 21.6;
%P = 10;



% 4 - CHOOSE THE LENGTH AND WIDTH OF THE PLOTS
ymin = 0.1;         % This is where you want the y-axis to begin
ymax = 15;          % and end
xmin = 0.01;        % The same for the x-axis
xmax = 4;



% 5 - CHOOSE THE NUMBER OF POINTS TO PLOT
num_points = 100;   % This is the number of points along each axis (i.e. an nxn grid of points)

colours = zeros(num_points);
% This is a matrix containing the 'colours' we want to plot for each point
% This will be filled with various numbers by the end

options = optimoptions('fsolve', 'TolX', 1e-10);
% These are the tolerances we want for fsolve



% The following loops find a 'number' for each point which
% indicates the type of stability at that point. It then places it
% into the matrix of zeros we created above. A breakdown of the number
% that we assign can be found below

ynum=0;
% 6a - PLACE THE PARAMETER YOU WANT ON THE Y-AXIS HERE
for I = linspace(ymin,ymax,num_points)
    ynum = ynum + 1;
    xnum = 0;
    % 6b - PLACE THE PARAMETER YOU WANT ON THE X-AXIS HERE 
    % (Final instruction)
    for P = linspace(xmin,xmax,num_points)
        xnum = xnum + 1;
        
        colour=[0,0,0,0]; 
        % To find the relevant colour, we first start by assigning each
        % point a 4D vector. Each number represents the number of stable
        % equilibria for a given equilibria type
        % colour = [Coexistence, M2 survival, M1 survival, Extinction]
        
        
        %First we search for Extinction equilibria%
    
        J_num = Jacobian(0,0,I,A11,A12,A21,A22,P); % Calculate the Jacobian
        eigenvalues = real(eig(J_num)); % Find the real parts of the eigenvalues
        %disp(eigenvalues);
        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 % Check they are all negative
            colour(4) = colour(4) + 1; %If so, increase the number of stable Coexistence equilibria by 1
        end


        %Next we do M1 survival equilibria

        G = @(M) m1solve(M,A11,I,P);
        for i = linspace(4,6.2,3)
            % We know that the maximum number of possible equilibria is 3
            % So we test three different starting values of fsolve
            solved = fsolve(G, i, options);
            %disp(solved)
            check = G(solved); % We check it is actually a root
            %disp(check)
            if -tolerance < check && check < tolerance % If it is as root
                equil = [r1(solved)/A11,0,solved]; % Find the equilibrium's co-ordinate
                equil = round(equil, 3); % Round to deal with numerical errors
                %disp(equil)
                if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 % If it is feasible
                    J_num = Jacobian(equil(1),equil(2),equil(3),A11,A12,A21,A22,P); % Calculate the Jacobian
                    eigenvalues = real(eig(J_num)); % Find the real eigenvalues
                    %disp(eigenvalues)
                    if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 % Check real parts are negative
                        colour(3) = colour(3) + 1; %If so, increase the number of stable M1 survival equilibria by 1
                    end
                end
            end
            G = @(M) G(M)/((M-solved)^2); 
            % Once a root has been found, we alter the equation so that
            % the root can no longer be found
        end


        % We then do the same for M2 survival

        G = @(M) m2solve(M,A22,I,P);
        for i = linspace(6,8.2,3)
            % We know that the maximum number of possible equilibria is 3
            % So we test three different starting values of fsolve
            solved = fsolve(G, i, options);
            %disp(solved)
            check = G(solved);  % We check it is actually a root
            %disp(check)
            if -tolerance < check && check < tolerance % If it is as root
                equil = [0,r2(solved)/A22,solved];  % Find the equilibrium's co-ordinate
                equil = round(equil, 3); % Round to deal with numerical errors
                %disp(equil)
                if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 % If it is feasible
                    J_num = Jacobian(equil(1),equil(2),equil(3),A11,A12,A21,A22,P);  % Calculate the Jacobian
                    eigenvalues = real(eig(J_num)); % Find the real eigenvalues
                    %disp(eigenvalues)
                    if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 % Check real parts are negative
                        colour(2) = colour(2) + 1; %If so, increase the number of stable M2 survival equilibria by 1
                    end
                end
            end
            G = @(M) G(M)/((M-solved)^2);
            % Once a root has been found, we alter the equation so that
            % the root can no longer be found
        end


        % We then do the same for Coexistence

        G = @(M) m3solve(M,A11,A12,A21,A22,I,P);
        for i = linspace(4,8.2,4)
            % This time, we know that the maximum number of possible equilibria is 4
            % So we test 4 different starting values of fsolve
            solved = fsolve(G, i);
            %disp(solved)
            check = G(solved); % We check it is actually a root
            %disp(check)
            if -tolerance < check && check < tolerance % If it is as root
                equil = [((A22*r1(solved))-(A12*r2(solved)))/(A11*A22-A12*A21),...
                ((A11*r2(solved))-(A21*r1(solved)))/(A11*A22-A12*A21),...
                solved]; % Find the equilibrium's co-ordinate
                equil = round(equil, 3); % Round to deal with numerical errors
                %disp(equil)
                if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 % If it is feasible
                    J_num = Jacobian(equil(1),equil(2),equil(3),A11,A12,A21,A22,P); % Calculate the Jacobian
                    eigenvalues = real(eig(J_num)); % Find the real eigenvalues
                    %disp(eigenvalues)
                    if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 % Check real parts are negative
                        colour(1) = colour(1) + 1; %If so, increase the number of stable coexistence equilibria by 1
                    end
                end
            end
            G = @(M) G(M)/((M-solved)^2);
            % Once a root has been found, we alter the equation so that
            % the root can no longer be found
        end


        %disp(colour)


        % The following color scheme is based on the idea that
        % the maximum number of stable equilibria for each type
        % of equilibria are:
        % Extinction = 1
        % M1 = 2
        % M2 = 2
        % Coexistence = 1 (this one might not be correct)

        if colour == [0,0,0,0] % Black
            finalcolour = 0;
        end
        if colour == [0,0,0,1] % Red
            finalcolour = 1;
        end
        if colour == [0,0,1,0] | colour == [0,0,2,0] % Light Green
            finalcolour = 2;
        end
        if colour == [0,1,0,0] | colour == [0,2,0,0] % Dark Green
            finalcolour = 3;
        end
        if colour == [1,0,0,0] % Blue
            finalcolour = 4;
        end
        if colour == [0,0,1,1] | colour == [0,0,2,1] % Light Yellow
            finalcolour = 5;
        end
        if colour == [0,1,0,1] | colour == [0,2,0,1] % Dark Yellow
            finalcolour = 6;
        end
        if colour == [1,0,0,1] % Magenta
            finalcolour = 7;
        end
        if colour == [0,1,1,0] | colour == [0,2,2,0] % Green
            finalcolour = 8;
        end
        if colour == [1,0,1,0] | colour == [1,0,2,0] % Light cyan
            finalcolour = 9;
        end
        if colour == [1,1,0,0] | colour == [1,2,0,0] % Dark cyan
            finalcolour = 10;
        end
        if colour == [0,1,1,1] | colour == [0,1,2,1] | colour == [0,2,1,1] | colour == [0,2,2,1] %Yellow
            finalcolour = 11;
        end
        if colour == [1,0,1,1] | colour == [1,0,2,1] % Light grey
            finalcolour = 12;
        end
        if colour == [1,1,0,1] | colour == [1,2,0,1] %Dark grey
            finalcolour = 13;
        end
        if colour == [1,1,1,0] | colour == [1,1,2,1] | colour == [1,2,1,0] | colour == [1,2,2,0] % Cyan
            finalcolour = 14;
        end
        if colour == [1,1,1,1] | colour == [1,2,1,1] | colour == [1,1,2,1] | colour == [1,2,2,1] % White
            finalcolour = 15;
        end
    
        colours(ynum,xnum) = finalcolour;
        
    end
end

%disp(colours)


% Below is the custom Color Map created with RGB Triplets%


customColorMap = [
    0, 0, 0;        % 0 = black
    1, 0, 0;        % 1 = red
    0.67, 1, 0;     % 2 = light green
    0, 0.5, 0;      % 3 = dark green
    0, 0, 1;        % 4 = blue
    1, 1, 0.67;     % 5 = light yellow
    0.75, 0.75, 0;  % 6 = dark yellow
    1, 0, 1;        % 7 = magenta
    0, 1, 0;        % 8 = green
    0.63, 1, 1;     % 9 = light cyan
    0, 0.74, 0.74;  % 10 = dark cyan
    1, 1, 0;        % 11 = yellow
    0.8, 0.8, 0.8;  % 12 = light grey
    0.4, 0.4, 0.4;  % 13 = dark grey
    0, 1, 1;        % 14 = cyan
    1, 1, 1;        % 15 = white
    ];


% This then plots the points using imagesc

yspace = linspace(ymin,ymax,num_points);
xspace = linspace(xmin,xmax,num_points);

colormap(customColorMap)
imagesc(xspace,yspace,colours)
set(gca, 'CLim', [0 length(customColorMap) - 1]);
axis xy
colorbar;
toc