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
    F(1) = D+P1*(sech(M-5)-0.2)/A11 - I/M;
end

%M2
function G = m2solve(x,A22,I,D,P2)
    M = x(1);
    G(1) = D+P2*(sech(M-6.5)-0.2)/A22 - I/M;
end

%Coexistence
function H = m3solve(x,A11,A12,A21,A22,I,D,P1,P2)
    M = x(1);
    H(1) = D+(sech(M-5)-0.2)*(P1*A22-P2*A21)+(sech(M-6.5)-0.2)*(P2*A11-P1*A12) - I/M;
end

% CHOOSE YOUR PARAMETERS %
% HIDE THE TWO YOU WANT TO PLOT %
A11 = 5;
A12 = 0;
A21 = 0;
A22 =  1;
%I = 21.6;
%D = 2.5;
P1 = 10;
P2 = 0.5;

% CHOOSE THE LENGTH AND WIDTH OF THE PLOTS %
% CHOOSE THE PARAMETERS YOU WANT TO PLOT%
ymin = 1;
ymax = 5;
xmin = 1;
xmax = 40;
num_points = 200;

colours = zeros(1,num_points);

for D = linspace(ymin,ymax,num_points) %y-axis
    Icolours = [];
    for I = linspace(xmin,xmax,num_points) %x-axis
        
        equilibria = [0,0,0]; %Create a list for equilibria%
        
        
        %Extinction%
        ext = [0,0,I/D];
        J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {ext(1), ext(2), ext(3), A11,A12,A21,A22,I,D,P1,P2});
        J_num=double(J_num); %Jacobian$
        eigenvalues = real(eig(J_num)); %Real parts of the eigenvalues
        %disp(eigenvalues);
        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check they are all negative
            equilibria = [equilibria; ext]; %If so, add to the equilibria list
        end

        %M1 survival
        M1range = linspace(2,8,5);
        for a = 1:length(M1range)
            x0 = M1range(a);
            solved = fsolve(@(x) m1solve(x,A11,I,D,P1), x0); %Solve for M1
            %disp(solved)
            equil = [(sech(solved-5)-0.2)/A11,0,solved]; %Equilibrium point
            %disp(equil)
            equil = round(equil, 3); %Avoids rounding errors%
                    if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 && ismember(equil, equilibria, 'rows')==0 
                        %If equilibrium unfeasible or already in list, ignore it$
                        J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                        J_num=double(J_num); %Jacobian$
                        eigenvalues = real(eig(J_num)); %Eigenvalues
                        %disp(eigenvalues)
                        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check they are all negative
                            equilibria = [equilibria; equil];
                        end
                    end
        end

        %M2 survival
        M2range = linspace(4,10,5);
        for a = 1:length(M2range)
            x0 = M2range(a);
            solved = fsolve(@(x) m2solve(x,A22,I,D,P2), x0); %Solve for M2
            %disp(solved)
            equil = [0,(sech(solved-6.5)-0.2)/A22,solved]; %Equilibrium point
            %disp(equil)
            equil = round(equil, 3); %Avoids rounding errors%
                    if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 && ismember(equil, equilibria, 'rows')==0 
                        %If equilibrium unfeasible or already in list, ignore it$
                        J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                        J_num=double(J_num); %Jacobian$
                        eigenvalues = real(eig(J_num)); %Eigenvalues
                        %disp(eigenvalues)
                        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check they are all negative
                            equilibria = [equilibria; equil];
                        end
                    end
        end

        %Coexistence
        M3range = linspace(0.1,10,10);
        for a = 1:length(M3range)
            x0 = M3range(a);
            solved = fsolve(@(x) m3solve(x,A11,A12,A21,A22,I,D,P1,P2), x0); %Solve for M3
            %disp(solved)
            equil = [(A22*(sech(solved-5)-0.2)-A12*(sech(solved-6.5)-0.2))/(A11*A22-A12*A21),(A11*(sech(solved-6.5)-0.2)-A21*(sech(solved-5)-0.2))/(A11*A22-A12*A21),solved];
            %disp(equil)
            equil = round(equil, 3); %Avoids rounding errors%
                    if equil(1) >= 0 && equil(2) >=0 && equil(3) >=0 && ismember(equil, equilibria, 'rows')==0 
                        %If equilibrium unfeasible or already in list, ignore it$
                        J_num = subs(J, {C1,C2,M,a11,a12,a21,a22,i,d,p1,p2}, {equil(1), equil(2), equil(3), A11,A12,A21,A22,I,D,P1,P2});
                        J_num=double(J_num); %Jacobian$
                        eigenvalues = real(eig(J_num)); %Eigenvalues
                        %disp(eigenvalues)
                        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0 %Check they are all negative
                            equilibria = [equilibria; equil];
                        end
                    end
        end

        
        equilibria(1,:) = []; %Gets rid of initial [0,0,0]
        %disp(equilibria) %Shows the feasible, stable equilibria for this set of parameters%

        %Each digit represents the number of stable equilibria of that type
        colour = [0,0,0,0]; %[co,m2,m1,ext]
        for k = 1:length(equilibria(:,1))
            if equilibria(k,1) == 0 && equilibria(k,2)==0 %extinction
                colour(4) = colour(4) + 1;
            end
            if equilibria(k,1) > 0 && equilibria(k,2)==0 %m1 survival
                colour(3) = colour(3) + 1;
            end
            if equilibria(k,1) ==0 && equilibria(k,2) > 0 %m2 survival
                colour(2) = colour(2) + 1;
            end
            if equilibria(k,1) > 0 && equilibria(k,2) > 0 %coexistence
                colour(1) = colour(1) + 1;
            end
        end
        %disp(colour)
        
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
%disp(length(colours))

customColorMap = [
    0, 0, 0;        %0 = black
    1, 0, 0;        %1 = red
    0.6, 1, 0;      %2 = light green
    0.12, 0.58, 0;  %3 = dark green
    0, 0, 1;        %4 = blue
    0.97, 1, 0.55;  %5 = light yellow
    0.84, 0.90, 0;  %6 = dark yellow
    0.60, 0, 1;     %7 = purple 
    0, 1, 0;        %8 = green
    0.12, 0.89, 0;  %9 = slightly light green
    0.15, 0.72, 0   %10 = slightly dark green
    %0.93, 1, 0      %11 = yellow
    ];

%disp(length(customColorMap))

Dspace = linspace(ymin,ymax,num_points);
Ispace = linspace(xmin,xmax,num_points);

%disp(Dspace)
%disp(Ispace)

colormap(customColorMap)
set(gca, 'CLim', [0 8]);
imagesc(Ispace,Dspace,colours)
axis xy
colorbar;