% Img=imread('polygen2.png');
% Img=rgb2gray(Img);
clear all
close all
load 'polygen.txt'
Img=polygen;
BlackPoint=[];
[BlackPoint(:,1), BlackPoint(:,2)]=find(Img==1);
figure
subplot(1,2,1)
scatter(BlackPoint(:,1), BlackPoint(:,2),'filled')
xmin=min(BlackPoint(:,1));
xmax=max(BlackPoint(:,1));
ymin=min(BlackPoint(:,2));
ymax=max(BlackPoint(:,2));
axis([xmin-1 xmax+1 ymin-1 ymax+1])
title('Before Fitting')
Img=Img(xmin:xmax,ymin:ymax);
clear BlackPoint
[BlackPoint(:,1), BlackPoint(:,2)]=find(Img==1);
[row, col]=size(Img);
Num_Points=row*col;
PointIndex=[1:Num_Points]';
Num_BP=size(BlackPoint,1);
idx=(BlackPoint(:,1)-1)*col+BlackPoint(:,2);
PointIndex(idx,2)=1;
clear idx
%top
PointIndex(:,3)=PointIndex(:,1)-col; PointIndex(1:col,3)=0; 
%botom
PointIndex(:, 4) = PointIndex(:, 1)+col; PointIndex(Num_Points-col+1:end, 4) = 0;
%left
PointIndex(:,5)=PointIndex(:,1)-1; PointIndex(1:col:row*col,5)=0;
%right
PointIndex(:,6)=PointIndex(:,1)+1; PointIndex(col:col:row*col,6)=0;
W=2;% weight of corners,1~n for points,n+1~2n for corner (top, bottom),2n+1~3n for corner(left, right)
f=[ones(Num_Points,1); W*ones(2*Num_Points,1)];
%constrains for existing black points
Aeq=[PointIndex(:,2)' zeros(1,2*Num_Points)];
beq=Num_BP;

%constrains for one circle
M=6;
A21=[M*eye(Num_Points);M*eye(Num_Points)];
for i=1:Num_Points
   idx=find(PointIndex(i,3:6)~=0)+2;
   A21(i,PointIndex(i,idx))=1; 
   A21(i+Num_Points,PointIndex(i,idx))=-1; 
end
b2=[(M+2)*ones(Num_Points,1);(M-2)*ones(Num_Points,1)];
A2=[A21 zeros(2*Num_Points, 2*Num_Points)];
clear idx

%constrains for top-bottom corner

b3=M*ones(2*Num_Points,1);
A31=[M*eye(Num_Points);M*eye(Num_Points)];
A32=[-eye(Num_Points);-eye(Num_Points)];
for i=1: Num_Points
    if PointIndex(i,3)~=0
        idxtop = PointIndex(i,3);
        A31(i,idxtop)=1;
        A31(i+Num_Points,idxtop)=-1;
    end
    if PointIndex(i,4)~=0
        idxbot=PointIndex(i,4);
        A31(i,idxbot)=-1;
        A31(i+Num_Points,idxbot)=1;    
    end
end
A3=[A31 A32 zeros(2*Num_Points, Num_Points)];

%constrains for left-right corner
b4=M*ones(2*Num_Points,1);
A41=[M*eye(Num_Points);M*eye(Num_Points)];
A42=[-eye(Num_Points);-eye(Num_Points)];
for i=1: Num_Points
    if PointIndex(i,5)~=0
        idxle = PointIndex(i,5);
        A41(i,idxle)=1;
        A41(i+Num_Points,idxle)=-1;
    end
    if PointIndex(i,6)~=0
        idxr=PointIndex(i,6);
        A41(i,idxr)=-1;
        A41(i+Num_Points,idxr)=1;    
    end
end
A4=[A41 zeros(2*Num_Points, Num_Points) A42];

% solve the integer linear programming problem
intcon=1:3*Num_Points;
A=[A2;A3;A4];
b=[b2;b3;b4];
lb=zeros(1,3*Num_Points);
ub=ones(1,3*Num_Points);
x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
FittingPoint=find(x(1:Num_Points)~=0);
col1=mod(FittingPoint,col);
row1=(FittingPoint-col1)./col+1;
idx=find(col1==0);
col1(idx)=17;
row1(idx)=row1(idx)-1;
subplot(1,2,2)
scatter(row1, col1,'filled')
axis([min(row1)-1 max(row1)+1 min(col1)-1 max(col1)+1])
title('After Fitting')
