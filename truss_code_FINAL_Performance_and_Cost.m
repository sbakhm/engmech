% EK301, Section C1, Summer 2018
% S Bakhmatova
% Truss Performance and Cost Calculation

clear; clc;
 
fid = fopen('Results.txt','w');

% Joint-to-Member Connection matrix (Joints in rows; Members in columns)
C = [
1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
1 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0;
0 0 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0;
0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 1;
0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0;
0 0 0 0 0 1 1 0 0 0 1 1 0 0 0 0 0;
0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0;
0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0;
0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1];

Sx = [1 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0];

Sy = [0 1 0;
0 0 0;
0 0 0;
0 0 0;
0 0 0;
0 0 1;
0 0 0;
0 0 0;
0 0 0;
0 0 0];

L = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];

X = [0.0 10.0 20.0 30.0 40.0 53.5  11.25          25.0         35.618     46.75];
Y = [0.0  0.0  0.0  0.0  0.0  0.0  -9.921567416  -14.14213562  -8.988775000  -7.378177282];

% Setup:
nJoints = length(X);
nMembers = 2*nJoints - 3;
Load = L(L~=0);

% Compilation of coefficient matrix A:
A = zeros(2*nJoints,nMembers);
TotLength = 0; % Combined length of all straws
Lengths = zeros(1,nMembers); % Individual straw lengths

for m = 1:nMembers
    Joints = transpose(find(C(:,m))); % Joints connected by member m
    Joint1 = Joints(1);
    Joint2 = Joints(2);
    XChange = X(Joint2) - X(Joint1);
    YChange = Y(Joint2) - Y(Joint1);
    Dist = sqrt(XChange*XChange + YChange*YChange);
    
    A(Joint1,m) = XChange / Dist;
    A(Joint1+nJoints,m) = YChange / Dist;
    A(Joint2,m) = -XChange / Dist;
    A(Joint2+nJoints,m) = -YChange / Dist;
    
    TotLength = TotLength + Dist;
    Lengths(1,m) = Dist;
end

Sxy = [Sx;Sy];
A = [A Sxy];

% Solving for member forces:
T = inv(A)*transpose(L);

% Determining truss cost, buckling strengths and uncertainty, critical
% member, max theoretical load, and member forces at theoretical max load:
Cost = 10*nJoints + TotLength;
Lengths2 = Lengths.*Lengths;
Lengths3 = Lengths2 .* Lengths;
Fbuckling = 1465.7282 ./ Lengths2;
U = 444.6187 ./ Lengths3;
[minSR,CritMember] = min(T(1:nMembers)./ transpose(Fbuckling)); % Scaling ratio
MaxLoad = Load/abs(minSR);
MaxLoadToCost = MaxLoad/Cost;
T_atMaxLoad = abs(MaxLoad/Load .* T); % Magnitudes of member forces at 
% theoretical max truss load

% Display of analysis results
fprintf(fid, 'EK301, Section C1, Summer 2018: Svetlana Bakhmatova, Sameer Hiremath, Edmund Yu\n');
fprintf(fid, 'Load in N: %2.0d\n',Load);
fprintf(fid, 'Member forces in N:\n');
for m = 1:nMembers
    Tm = round(T(m),3);
    if Tm < 0
        fprintf(fid, 'm%2.2d:  %6.3f  (C)\n', m,abs(Tm));
    elseif Tm == 0.000
        fprintf(fid, 'm%2.2d:  %6.3f  (Z)\n', m,abs(Tm)); % Z = zero-force
    else
        fprintf(fid, 'm%2.2d:  %6.3f  (T)\n', m,abs(Tm));
    end
end

fprintf(fid, 'Reaction forces in N:\n');
fprintf(fid, 'Sx1:  %6.3f\n', T(end-2));
fprintf(fid, 'Sy1:  %6.3f\n', T(end-1));
fprintf(fid, 'Sy2:  %6.3f\n', T(end));

fprintf(fid, 'Cost of truss:  $%6.3f\n',Cost);
fprintf(fid, 'Theoretical max load/cost ratio in N/$: %6.3f\n',MaxLoadToCost);
fprintf(fid, 'Theoretical max load in N:  %6.3f\n',MaxLoad);
fprintf(fid, 'Critical member:  m%2.2d\n',CritMember);

% Display of truss design parameters:
fprintf(fid, 'Truss design parameters:\n');
fprintf(fid, 'mNo   Length   Member     Buckling     Uncertainty  ForceAt\n');
fprintf(fid, '      (cm)     Force(N)   Strength(N)  (N)          MaxLoad(N)\n');
fprintf(fid, '=============================================================\n');
for m = 1:nMembers
    Tm = round(T(m),3);
    if Tm < 0
        fprintf(fid, 'm%2.2d:  %6.3f  %6.3f (C)  %6.3f      %6.3f        %6.3f\n', ...
            m,Lengths(m),abs(Tm),Fbuckling(m),U(m),T_atMaxLoad(m));
    elseif Tm == 0.000
        fprintf(fid, 'm%2.2d:  %6.3f  %6.3f (Z)      -            -        %6.3f\n',...
            m,Lengths(m),abs(Tm),T_atMaxLoad(m));
    else
        fprintf(fid, 'm%2.2d:  %6.3f  %6.3f (T)      -            -        %6.3f\n',...
            m,Lengths(m),abs(Tm),T_atMaxLoad(m));
    end
end

fclose(fid);
type Results.txt;