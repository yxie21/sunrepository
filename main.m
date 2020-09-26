close all;clear;clc;

%% w dynamic

%% calculate
N=3; k=-0.75; eta=0.1; n=1;
odeFunc=@(t,y)func(y,N,n,k,eta);
ts = linspace(0,100,1e4);
y0 = [0; 1]; % y = [\omega; \zeta];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
[ts,y] = ode113(odeFunc,ts,y0,opts);

%% plot figure
fig = figure;
set(fig, 'position', get(0,'ScreenSize')); % Fullscreen
plotHsv(real(y(:,1)),imag(y(:,1)),angle(y(:,2))/pi);
hold on;
scatter3(real(y(1,1)),imag(y(1,1)),angle(y(1,2))/pi,100,'ro','filled'); % start point
scatter3(real(y(end,1)),imag(y(end,1)),angle(y(end,2))/pi,100,'yo','filled'); % end point
set(gca,'FontSize',30);
xlabel('Re[\omega]');
ylabel('Im[\omega]');
axis equal
zlabel('¡Ï\zeta/\pi');
grid on;
xlim([-1 1]);
ylim([-1 1]);
title(['n = ',num2str(n)],'FontSize',36);
plotcircle();
% view(0,0); % Re[\omega] vs ¡Ï\zeta
% view(90,0) % Im[\omega] vs ¡Ï\zeta
% view(0,90); % Re[\omega] vs Im[\omega] 
view(30,15);

%% main ode function
function dy = func(y,N,n,k,eta)
persistent an tmpn theta0 M;
if isempty(an) || isempty(tmpn) || tmpn ~= n
    an = 2*pi / integral(@(t)((1-cos(t)).^n),0,2*pi);
    tmpn = n; theta0 = linspace(0,2*pi*(N-1)/N,N)';
    M= @(z,w,zeta) zeta.*(z-w)./(1-conj(w).*z);
end
dy = zeros(size(y));
w = y(1); z = y(2); theta = angle(M(exp(1i*theta0),y(1),y(2)));
I = an/N*sum((1-cos(theta)).^n);
dy(1) = -0.5*1i*(1-abs(w)^2)*conj(z)*(eta+k*I-1);
dy(2) = 1i*(1+eta+k*I)*z-0.5*1i*(conj(w)+w*z^2)*(eta+k*I-1);
end

function plotcircle()
t = linspace(0,2*pi,360);
x = cos(t); y = sin(t);
plot(x,y,'color','k','linewidth',5);
end

function plotHsv(x,y,z)
cmap = colormap('hsv');
% change c into an index into the colormap
% min(c) -> 1, max(c) -> number of colors
c = round(1+(size(cmap,1)-1)*(z-min(z))/2);
% make a blank plot
plot3(x,y,z,'LineStyle','None');
% add line segments
for k = 1:(length(x)-1)
    line(x(k:k+1),y(k:k+1),z(k:k+1),'color',cmap(c(k),:),'LineWidth',3);
end
colorbar;
caxis([-1 1]);
end