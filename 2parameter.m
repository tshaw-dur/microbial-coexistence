syms C1 C2 M a11 a12 a21 a22 i d p1 p2
f1 = C1*(sech(M-5)-0.2-a11*C1-a12*C2);
f2 = C2*(sech(M-6.5)-0.2-a21*C1-a22*C2);
f3 = i-d*M-p1*M*C1-p2*M*C2;
F=[f1,f2,f3];
vars=[C1;C2;M];

%Jacobian matrix
J = jacobian(F,vars);

%M1
function F = m1solve(x,A11,I,D,P1)
    M = x(1);
    F(1) = M*(D+((sech(M-5)-0.2)*P1/A11)) - I;
end

%M2
function G = m2solve(x,A22,I,D,P2)
    M = x(1);
    G(1) = M*(D+((sech(M-6.5)-0.2)*P2/A22)) - I;
end

%Coexistence
function H = m3solve(x,A11,A12,A21,A22,I,D,P1,P2)
    M = x(1);
    det = A11*A22 - A12*A21;
    H(1) = M*(D*det + (sech(M-5)-0.2)*(A22*P1-A21*P2) + (sech(M-6.5)-0.2)*(A11*P2-A12*P1)) - I*det;
end

% CHOOSE YOUR PARAMETERS %
% HIDE THE TWO YOU WANT TO PLOT %
A11 = 2;
A12 = 0;
A21 = 0;
A22 = 1.4;
%I = 21.6;
%D = 2.5;
P1 = 10;
P2 = 2;

% CHOOSE THE LENGTH AND WIDTH OF THE PLOTS %
% CHOOSE THE PARAMETERS YOU WANT TO PLOT%
ymin = 1;
ymax = 5;
xmin = 1;
xmax = 40;
num_points = 100;

colours = zeros(1,num_points);

for D = linspace(ymin,ymax,num_points) %y-axis
    Icolours = [];
    for I = linspace(xmin,xmax,num_points) %x-axis
        
        colour=[0,0,0,0]; %[co,m2,m1,ext]
        
        
        %EXTINCTION%
        J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {0, 0, I/D, A11,A12,A21,A22,I,D,P1,P2});
        J_num = double(J_num); %Jacobian$
        eigenvalues = real(eig(J_num)); %Real parts of the eigenvalues
        %disp(eigenvalues);
        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check they are all negative
            colour(4) = colour(4) + 1;
        end




        %M1 SURVIVAL
        solved1 = fsolve(@(x) m1solve(x,A11,I,D,P1), 4);
        equil = [(sech(solved1-5)-0.2)/A11,0,solved1]; %Equilibrium point
        check = m1solve(solved1,A11,I,D,P1); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(3) = colour(3) + 1;
                end
            end
        end

        solved2 = fsolve(@(x) m1solve(x,A11,I,D,P1)/(x-solved1)^2, 4);
        equil = [(sech(solved2-5)-0.2)/A11,0,solved2]; %Equilibrium point
        check = m1solve(solved2,A11,I,D,P1); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(3) = colour(3) + 1;
                end
            end
        end

        solved3 = fsolve(@(x) m1solve(x,A11,I,D,P1)/(((x-solved1)^2)*((x-solved2)^2)) , 4);
        equil = [(sech(solved3-5)-0.2)/A11,0,solved3]; %Equilibrium point
        check = m1solve(solved3,A11,I,D,P1); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(3) = colour(3) + 1;
                end
            end
        end




        %M2 SURVIVAL
        solved1 = fsolve(@(x) m2solve(x,A22,I,D,P2), 6);
        equil = [0,(sech(solved1-6.5)-0.2)/A22,solved1]; %Equilibrium point
        check = m2solve(solved1,A22,I,D,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(2) = colour(2) + 1;
                end
            end
        end

        solved2 = fsolve(@(x) m2solve(x,A22,I,D,P2)/(x-solved1)^2, 7);
        equil = [0,(sech(solved2-6.5)-0.2)/A22,solved2]; %Equilibrium point
        check = m2solve(solved2,A22,I,D,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(2) = colour(2) + 1;
                end
            end
        end

        solved3 = fsolve(@(x) m2solve(x,A22,I,D,P2)/(((x-solved1)^2)*((x-solved2)^2)), 9);
        equil = [0,(sech(solved3-6.5)-0.2)/A22,solved3]; %Equilibrium point
        check = m2solve(solved3,A22,I,D,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(2) = colour(2) + 1;
                end
            end
        end




        %COEXISTENCE
        solved1 = fsolve(@(x) m3solve(x,A11,A12,A21,A22,I,D,P1,P2), 2);
        equil = [((A22*(sech(solved1-5)-0.2))-(A12*(sech(solved1-6.5)-0.2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved1-6.5)-0.2))-(A21*(sech(solved1-5)-0.2)))/(A11*A22-A12*A21),...
        solved1]; %Equilibrium point
        check = m3solve(solved1,A11,A12,A21,A22,I,D,P1,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end

        solved2 = fsolve(@(x) m3solve(x,A11,A12,A21,A22,I,D,P1,P2)/(x-solved1)^2, 4.5);
        equil = [((A22*(sech(solved2-5)-0.2))-(A12*(sech(solved2-6.5)-0.2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved2-6.5)-0.2))-(A21*(sech(solved2-5)-0.2)))/(A11*A22-A12*A21),...
        solved2]; %Equilibrium point
        check = m3solve(solved2,A11,A12,A21,A22,I,D,P1,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end

        solved3 = fsolve(@(x) m3solve(x,A11,A12,A21,A22,I,D,P1,P2)/(((x-solved1)^2)*((x-solved2)^2)), 7);
        equil = [((A22*(sech(solved3-5)-0.2))-(A12*(sech(solved3-6.5)-0.2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved3-6.5)-0.2))-(A21*(sech(solved3-5)-0.2)))/(A11*A22-A12*A21),...
        solved3]; %Equilibrium point
        check = m3solve(solved3,A11,A12,A21,A22,I,D,P1,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end
        
        solved4 = fsolve(@(x) m3solve(x,A11,A12,A21,A22,I,D,P1,P2)/(((x-solved1)^2)*((x-solved2)^2)*((x-solved3)^2)), 9.5);
        equil = [((A22*(sech(solved4-5)-0.2))-(A12*(sech(solved4-6.5)-0.2)))/(A11*A22-A12*A21),...
        ((A11*(sech(solved4-6.5)-0.2))-(A21*(sech(solved4-5)-0.2)))/(A11*A22-A12*A21),...
        solved4]; %Equilibrium point
        check = m3solve(solved4,A11,A12,A21,A22,I,D,P1,P2); %Make sure it actually is
        %disp(equil)
        %disp(check)
        if -1e-4 < check && check < 1e-4
            equil = round(equil, 3);
            if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 %Check feasibility
                J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                J_num = double(J_num); %Jacobian
                eigenvalues = real(eig(J_num)); %Eigenvalues
                %disp(eigenvalues)
                if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check real parts are negative
                    colour(1) = colour(1) + 1;
                end
            end
        end

        disp(colour)



        if colour == [0,0,0,0] %black
            finalcolour = 0;
        end
        if colour == [0,0,0,1] %red
            finalcolour = 1;
        end
        if colour == [0,0,1,0] %light green
            finalcolour = 2;
        end
        if colour == [0,1,0,0] %dark green
            finalcolour = 3;
        end
        if colour == [1,0,0,0] %blue
            finalcolour = 4;
        end
        if colour == [0,0,1,1] %light yellow
            finalcolour = 5;
        end
        if colour == [0,1,0,1] %dark yellow
            finalcolour = 6;
        end
        if colour == [1,0,0,1] %purple
            finalcolour = 7;
        end
        if colour == [0,1,1,0] %green
            finalcolour = 8;
        end
        if colour == [0,1,2,0] %slightly light green
            finalcolour = 9;
        end
        if colour == [0,2,1,0] %slightly dark green
            finalcolour = 10;
        end
        if colour == [0,1,1,1] %yellow
            finalcolour = 11;
        end
        Icolours = [Icolours, finalcolour];
        
    end
    %disp(Icolours)
    colours=[colours;Icolours];
end

colours(1,:) = [];
disp(colours)


customColorMap = [
    0, 0, 0        %0 = black
    1, 0, 0;        %1 = red
    0.6, 1, 0;      %2 = light green
    0.12, 0.58, 0;  %3 = dark green
    0, 0, 1;        %4 = blue
    0.97, 1, 0.55;  %5 = light yellow
    0.84, 0.90, 0;  %6 = dark yellow
    0.60, 0, 1     %7 = purple 
    %0, 1, 0;        %8 = green
    %0.12, 0.89, 0;  %9 = slightly light green
    %0.15, 0.72, 0;   %10 = slightly dark green
    %0.93, 1, 0;      %11 = yellow
    ];

%disp(length(customColorMap))

Dspace = linspace(ymin,ymax,num_points);
Ispace = linspace(xmin,xmax,num_points);

%disp(Dspace)
%disp(Ispace)

colormap(customColorMap)
set(gca, 'CLim', [1 7]);
imagesc(Ispace,Dspace,colours)
axis xy
colorbar;