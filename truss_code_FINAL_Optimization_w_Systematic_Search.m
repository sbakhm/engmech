% EK301, Section C1, Summer 2018
% S Bakhmatova
% Truss Design Optimization (Systematic Search)

clear; clc;

% Joint-to-Member Connection matrix (Joints in rows; Members in columns)
% m1 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
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
Load = L(L~=0);

BestLoad = 0;
WorstLoad = 0;
BestCost = 0;
XBest = 0;
YBest = 0;
BestLengths = zeros(1,17);

% Sets of points:
Radius = 10.0:0.1:15.0;
Angle = transpose(0:1:360);
nR = length(Radius);
nA = length(Angle);
RadiusRep = repmat(Radius,nA,1);
RadiusMtx = transpose(RadiusRep(1:nR*nA));
AngleMtx = repmat(Angle,nR,1);

XY = [RadiusMtx .* cosd(AngleMtx), RadiusMtx .* sind(AngleMtx) ];

Cir1 = round(XY,1);
Cir2 = round([10 0.0] + XY,1);
Cir3 = round([20 0.0] + XY,1);
Cir4 = round([30 0.0] + XY,1);
Cir5 = round([40 0.0] + XY,1);
Cir6 = round([53.5 0.0] + XY,1);
Cir8 = round([25.0 -14.14213562] + XY,1);

Cir1(Cir1(:,2)>0,:) = [];
Cir2(Cir2(:,2)>0,:) = [];
Cir3(Cir3(:,2)>0,:) = [];
Cir4(Cir4(:,2)>0,:) = [];
Cir5(Cir5(:,2)>0,:) = [];
Cir6(Cir6(:,2)>0,:) = [];
Cir8(Cir8(:,2)>0,:) = [];

Coord7 = intersect(intersect(Cir1,Cir2,'rows'),intersect(Cir3,Cir8,'rows'),'rows');
Coord9 = intersect(Cir4,intersect(Cir5,Cir8,'rows'),'rows');
Coord10 = intersect(Cir5,Cir6,'rows');

for i7 = 1:size(Coord7,1)
    J7 = Coord7(i7,:);
    
    for i9 = 1:size(Coord9,1)
        J9 = Coord9(i9,:);
        
      for i10 = 1:size(Coord10,1)
          J10 = Coord10(i10,:);
          
          disp([J7(1,1) J7(1,2) J9(1,1) J9(1,2) J10(1,1) J10(1,2)]);
          
          X = [0.0 10.0 20.0 30.0 40.0 53.5  J7(1,1)  25    J9(1,1)  J10(1,1)];
          Y = [0.0  0.0  0.0  0.0  0.0  0.0  J7(1,2) -14.14213562  J9(1,2)  J10(1,2)];
          
          nJoints = length(X);
          nMembers = 2*nJoints - 3;

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
             Dist = round(sqrt(XChange*XChange + YChange*YChange),5);

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

              if MaxLoad > BestLoad
                 BestCost = Cost;
                 BestLoad = MaxLoad;
                 BestLengths = Lengths;
                 XBest = X;
                 YBest = Y;  
              end
              disp(MaxLoad);
      end        
    end   
end

disp(['Best max load: ',num2str(BestLoad)]);
disp(['Worst max load: ',num2str(WorstLoad)]);
disp(['Best cost: ',num2str(BestCost)]);
disp(['Best X: ',num2str(XBest)]);
disp(['Best Y: ',num2str(YBest)]);
disp('Best Lengths: ');
disp(BestLengths);

% 
