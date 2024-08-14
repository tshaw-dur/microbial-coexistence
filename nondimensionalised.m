tic

syms C1 C2 M k1 k2 h1 h2 a11 a12 a21 a22 i p
f1 = C1*(sech(M-k1)-h1-a11*C1-a12*C2);
f2 = C2*(sech(M-k2)-h2-a21*C1-a22*C2);
f3 = i-M*(1+p*(C1+C2));
F=[f1,f2,f3];
vars=[C1;C2;M];

%Jacobian matrix
J = jacobian(F,vars);

%M1
function F = m1solve(x,K1,H1,A11,I,P)
    M = x(1);
    F(1) = M*(1+((sech(M-K1)-H1)*P/A11)) - I;
end

%M2
function G = m2solve(x,K2,H2,A22,I,P)
    M = x(1);
    G(1) = M*(1+((sech(M-K2)-H2)*P/A22)) - I;
end

%Coexistence
function H = m3solve(x,K1,K2,H1,H2,A11,A12,A21,A22,I,P)
    M = x(1);
    det = A11*A22 - A12*A21;
    H(1) = M*(1 + ((sech(M-K1)-H1)*(A22-A21) + (sech(M-K2)-H2)*(A11-A12))*P/det) - I;
end

% CHOOSE YOUR PARAMETERS %
% HIDE THE TWO YOU WANT TO PLOT %
K1 = 5;
K2 = 6;
H1 = 0.2;
H2 = 0.2;

A11 = 1;
A12 = 0;
A21 = 0;
A22 = 2;
%I = 21.6;
%P = 10;

% CHOOSE THE LENGTH AND WIDTH OF THE PLOTS %
% CHOOSE THE PARAMETERS YOU WANT TO PLOT%
ymin = 0.01;
ymax = 20;
xmin = 0.01;
xmax = 2;
num_points = 100;

colours = zeros(1,num_points);

for I = linspace(ymin,ymax,num_points) %y-axis
    Icolours = [];
    for P = linspace(xmin,xmax,num_points) %x-axis
        
        colour=[0,0,0,0]; %[co,m2,m1,ext]
        
        
        %EXTINCTION%
        J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {0, 0, I, K1,K2,H1,H2,A11,A12,A21,A22,I,P});
        J_num = double(J_num); %Jacobian$
        eigenvalues = real(eig(J_num)); %Real parts of the eigenvalues
        %disp(eigenvalues);
        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check they are all negative
            colour(4) = colour(4) + 1;
        end
        %disp(colour)



        %M1 SURVIVAL
        solved1 = fsolve(@(x) m1solve(x,K1,H1,A11,I,P), 4);
        equil = [(sech(solved1-K1)-H1)/A11,0,solved1]; %Equilibrium point
        check = m1solve(solved1,K1,H1,A11,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(3) = colour(3) + 1;
                end
            end
        end

        solved2 = fsolve(@(x) m1solve(x,K1,H1,A11,I,P)/(x-solved1)^2, 4);
        equil = [(sech(solved2-K1)-H1)/A11,0,solved2]; %Equilibrium point
        check = m1solve(solved2,K1,H1,A11,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(3) = colour(3) + 1;
                end
            end
        end

        solved3 = fsolve(@(x) m1solve(x,K1,H1,A11,I,P)/(((x-solved1)^2)*((x-solved2)^2)) , 4.01);
        equil = [(sech(solved3-K1)-H1)/A11,0,solved3]; %Equilibrium point
        check = m1solve(solved3,K1,H1,A11,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(3) = colour(3) + 1;
                end
            end
        end




        %M2 SURVIVAL
        solved1 = fsolve(@(x) m2solve(x,K2,H2,A22,I,P), 6);
        equil = [0,(sech(solved1-K2)-H2)/A22,solved1]; %Equilibrium point
        check = m2solve(solved1,K2,H2,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(2) = colour(2) + 1;
                end
            end
        end

        solved2 = fsolve(@(x) m2solve(x,K2,H2,A22,I,P)/(x-solved1)^2, 7);
        equil = [0,(sech(solved2-K2)-H2)/A22,solved2]; %Equilibrium point
        check = m2solve(solved2,K2,H2,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(2) = colour(2) + 1;
                end
            end
        end

        solved3 = fsolve(@(x) m2solve(x,K2,H2,A22,I,P)/(((x-solved1)^2)*((x-solved2)^2)), 9);
        equil = [0,(sech(solved3-K2)-H2)/A22,solved3]; %Equilibrium point
        check = m2solve(solved3,K2,H2,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(2) = colour(2) + 1;
                end
            end
        end




        %COEXISTENCE
        solved1 = fsolve(@(x) m3solve(x,K1,K2,H1,H2,A11,A12,A21,A22,I,P), 2);
        equil = [((A22*(sech(solved1-K1)-H1))-(A12*(sech(solved1-K2)-H2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved1-K2)-H2))-(A21*(sech(solved1-K1)-H1)))/(A11*A22-A12*A21),...
        solved1]; %Equilibrium point
        check = m3solve(solved1,K1,K2,H1,H2,A11,A12,A21,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end

        solved2 = fsolve(@(x) m3solve(x,K1,K2,H1,H2,A11,A12,A21,A22,I,P)/(x-solved1)^2, 4.5);
        equil = [((A22*(sech(solved2-K1)-H1))-(A12*(sech(solved2-K2)-H2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved2-K2)-H2))-(A21*(sech(solved2-K1)-H1)))/(A11*A22-A12*A21),...
        solved2]; %Equilibrium point
        check = m3solve(solved2,K1,K2,H1,H2,A11,A12,A21,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end

        solved3 = fsolve(@(x) m3solve(x,K1,K2,H1,H2,A11,A12,A21,A22,I,P)/(((x-solved1)^2)*((x-solved2)^2)), 7);
        equil = [((A22*(sech(solved3-K1)-H1))-(A12*(sech(solved3-K2)-H2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved3-K2)-H2))-(A21*(sech(solved3-K1)-H1)))/(A11*A22-A12*A21),...
        solved3]; %Equilibrium point
        check = m3solve(solved3,K1,K2,H1,H2,A11,A12,A21,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end
        
        solved4 = fsolve(@(x) m3solve(x,K1,K2,H1,H2,A11,A12,A21,A22,I,P)/(((x-solved1)^2)*((x-solved2)^2)*((x-solved3)^2)), 9.4);
        equil = [((A22*(sech(solved4-K1)-H1))-(A12*(sech(solved4-K2)-H2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved4-K2)-H2))-(A21*(sech(solved4-K1)-H1)))/(A11*A22-A12*A21),...
        solved4]; %Equilibrium point
        check = m3solve(solved4,K1,K2,H1,H2,A11,A12,A21,A22,I,P); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,k1,k2,h1,h2,a11,a12,a21,a22,i,p}, {equil(1), equil(2), equil(3), K1,K2,H1,H2,A11,A12,A21,A22,I,P});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end

        disp(colour)

        % The following color scheme is based on the idea that
        % the maximum number of stable equilibria for each type
        % of equilibria are:
        % Extinction = 1
        % M1 = 2
        % M2 = 2
        % Coexistence = 1 (this one might not be correct)

        if colour == [0,0,0,0] %black
            finalcolour = 0;
        end
        if colour == [0,0,0,1] %red
            finalcolour = 1;
        end
        if colour == [0,0,1,0] | colour == [0,0,2,0] %light green
            finalcolour = 2;
        end
        if colour == [0,1,0,0] | colour == [0,2,0,0] %dark green
            finalcolour = 3;
        end
        if colour == [1,0,0,0] %blue
            finalcolour = 4;
        end
        if colour == [0,0,1,1] | colour == [0,0,2,1] %light yellow
            finalcolour = 5;
        end
        if colour == [0,1,0,1] | colour == [0,2,0,1] %dark yellow
            finalcolour = 6;
        end
        if colour == [1,0,0,1] %magenta
            finalcolour = 7;
        end
        if colour == [0,1,1,0] | colour == [0,2,2,0] %green
            finalcolour = 8;
        end
        if colour == [1,0,1,0] | colour == [1,0,2,0] %light cyan
            finalcolour = 9;
        end
        if colour == [1,1,0,0] | colour == [1,2,0,0] %dark cyan
            finalcolour = 10;
        end
        if colour == [0,1,1,1] | colour == [0,1,2,1] | colour == [0,2,1,1] | colour == [0,2,2,1]%yellow
            finalcolour = 11;
        end
        if colour == [1,0,1,1] | colour == [1,0,2,1] %light grey
            finalcolour = 12;
        end
        if colour == [1,1,0,1] | colour == [1,2,0,1]%dark grey
            finalcolour = 13;
        end
        if colour == [1,1,1,0] | colour == [1,1,2,1] | colour == [1,2,1,0] | colour == [1,2,2,0] %cyan
            finalcolour = 14;
        end
        if colour == [1,1,1,1] | colour == [1,2,1,1] | colour == [1,1,2,1] | colour == [1,2,2,1] %white
            finalcolour = 15;
        end
    
        Icolours = [Icolours, finalcolour];
        
    end
    %disp(Icolours)
    colours=[colours;Icolours];
end

colours(1,:) = [];
disp(colours)


customColorMap = [
    0, 0, 0;        %0 = black
    1, 0, 0;        %1 = red
    0.67, 1, 0;     %2 = light green
    0, 0.5, 0;      %3 = dark green
    0, 0, 1;        %4 = blue
    1, 1, 0.67;     %5 = light yellow
    0.75, 0.75, 0;  %6 = dark yellow
    1, 0, 1;        %7 = magenta
    0, 1, 0;        %8 = green
    0.63, 1, 1;     %9 = light cyan
    0, 0.74, 0.74;  %10 = dark cyan
    1, 1, 0;        %11 = yellow
    0.8, 0.8, 0.8;  %12 = light grey
    0.4, 0.4, 0.4;  %13 = dark grey
    0, 1, 1;        %14 = cyan
    1, 1, 1;        %15 = white
    ];

%disp(length(customColorMap))

yspace = linspace(ymin,ymax,num_points);
xspace = linspace(xmin,xmax,num_points);

%disp(Dspace)
%disp(Ispace)

colormap(customColorMap)
imagesc(xspace,yspace,colours)
set(gca, 'CLim', [0 15]);
axis xy
colorbar;
toc