%GOPH 547
%Ye Sun
%Safian Omar Qureshi
%ID: 10086638
%Collaborated with Tian Yu, Younghwan Ryan Ok, Ahtisham Mohummad 

clear;
clc;

%%
load('goph_547_w2017_lab4_data.mat')

%Q1
figure
contourf(X,Y,Fz_raw); %using variables supplied in raw data
h_c=colorbar; %simply creating a 2D contour plot 
title('2D Contour Plot of Raw Data')
xlabel('X (m)');
ylabel('Y (m)');
ylabel(h_c,'Magnetic Effect (nT)');

figure
subplot(2,1,1); %subplot function for multiple plots on same figure
plot(X,Fz_raw,'.'); %using linear type of plot for X only
title('Raw Fz vs X')
xlabel('X (m)');
ylabel('Raw Fz (nT)');

a_x=polyfit(X,Fz_raw,1);
linear_x=min(X(:)):1:max(X(:))+1;
linear_y=a_x(1)*linear_x+a_x(2);
hold on;
plot(linear_x,linear_y);

subplot(2,1,2);
plot(Y,Fz_raw,'.'); %using linear type of plot for Y only
title('Raw Fz vs Y')
xlabel('Y (m)');
ylabel('Raw Fz (nT)');

%%
%Q2
a_y=polyfit(Y,Fz_raw,1);
linear_x=min(Y(:)):1:max(Y(:))+1;
linear_y=a_y(1)*linear_x+a_y(2);
hold on;
plot(linear_x,linear_y);


%% 
%Q3 
Fz=Fz_raw-a_y(1)*Y; %here we remove linear component of regional variation in y direction

figure
contourf(X,Y,Fz); %simply creating a 2D contour plot 
h_c=colorbar;
title('2D Contour Plot - x-component regional variations removed')
xlabel('X (m)');
ylabel('Y (m)');
ylabel(h_c,'Magnetic Effect (nT)');

figure %simply creating figures 
subplot(2,1,1);
plot(X,Fz,'.');
title('Fz y-component regional variations removed plotted against X')
xlabel('X (m)');
ylabel('Fz (nT)');

a_x=polyfit(X,Fz,1);
linear_x=min(X(:)):1:max(X(:))+1;
linear_y=a_x(1)*linear_x+a_x(2);
hold on;
plot(linear_x,linear_y);

subplot(2,1,2); %simply creating figures 
plot(Y,Fz,'.');
title('Fz y-component of regional variations removed plotted against Y')
xlabel('Y (m)');
ylabel('Fz (nT)');

a_y=polyfit(Y,Fz,1);
linear_x=min(Y(:)):1:max(Y(:))+1;
linear_y=a_y(1)*linear_x+a_y(2);
hold on;
plot(linear_x,linear_y);

%%
%Q4: 
Fz=Fz-a_x(1)*X; %here we remove linear component of regional variation in x direction, mostly same as question 3

figure
contourf(X,Y,Fz); %simply creating figures 
h_c=colorbar;
title('2D Contour Plot - x-component regional variations removed')
xlabel('X (m)');
ylabel('Y (m)');
ylabel(h_c,'Magnetic Effect (nT)');

figure
subplot(2,1,1);
plot(X,Fz,'.');
title('Fz x-component of regional variations removed plotted against X')
xlabel('X (m)');
ylabel('Fz (nT)');

a_x=polyfit(X,Fz,1);
linear_x=min(X(:)):1:max(X(:))+1;
linear_y=a_x(1)*linear_x+a_x(2);
hold on;
plot(linear_x,linear_y);

subplot(2,1,2);
plot(Y,Fz,'.');
title('Fz x-component of regional variations removed plotted against Y')
xlabel('Y (m)');
ylabel('Fz (nT)');

a_y=polyfit(Y,Fz,1);
linear_x=min(Y(:)):1:max(Y(:))+1;
linear_y=a_y(1)*linear_x+a_y(2);
hold on;
plot(linear_x,linear_y);

%%
%Q5:
Fz=Fz-min(Fz(:)); %removing constant component of regional, subtracting min Fz

figure
contourf(X,Y,Fz); %simply creating a 2D contour plot 
h_c=colorbar;
title('Contour Plot after removal of minimum Fz')
xlabel('X (m)');
ylabel('Y (m)');
ylabel(h_c,'Magnetic Effect (nT)');

%Q6:
dx=30;
dy=dx;
dA=dx*dy;
h=30;

Fz_u=zeros(size(Fz)); %here we have an algorithm to implement the double integral for upward continuation
for i=1:length(Fz(:,1)) %using for loops for every column/row
    for j=1:length(Fz(1,:))
        x_30=[X(i,j),Y(i,j),-h];  %specifying the x,y points in our data set
        for n=1:length(Fz(:,1))
            for m=1:length(Fz(1,:)) %specifying to run loop for all Fz data values
                x_o=[X(n,m),Y(n,m),0]; %the integration itself
                r=sqrt((x_o(1)-x_30(1))^2+(x_o(2)-x_30(2))^2+h^2); %the integration itself
                up_Fz=(h/(2*pi))* (Fz(n,m)/(r^3))*dA; %the integration itself
                Fz_u(i,j)=up_Fz+Fz_u(i,j); %applying upward continuatiion    
            end
        end
    end
end

figure
contourf(X,Y,Fz_u);  %simply creating a 2D contour plot 
h_c=colorbar;
title('Contour Plot - Upward Continuation')
xlabel('X (m)');
ylabel('Y (m)');
ylabel(h_c,'Magnetic Effect (nT)');

%%
%Q7:
Fz_d=Fz;

for j3=2:size(Fz_d,1)-1 %here we have an algorithm to implement downward continuation
    for i3=2:size(Fz_d,2)-1 %specifying the column and of what Fz values to 'hit'
        Fz_d(j3,i3)=6*Fz(j3,i3)-(Fz(j3-1,i3)+Fz(j3+1,i3)+Fz(j3,i3-1)+Fz(j3,i3+1)+Fz_u(j3,i3)); %applying given formula 
    end
end

figure
contourf(X,Y,Fz_d);  %simply creating a 2D contour plot 
h_c=colorbar;
title('Contour Plot - Downward Continuation')
xlabel('X (m)');
ylabel('Y (m)');
ylabel(h_c,'Magnetic Effect (nT)');