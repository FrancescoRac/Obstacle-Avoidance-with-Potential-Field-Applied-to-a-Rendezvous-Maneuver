%% Manovra Completa
clc
clear
%% Condizioni operative
tf = 2500; %Tempo
N = 80000; %Passi
h = tf/N; %Passo numerico

%% Condizioni iniziali
m = 600; %Massa
w = 0.00113; %Velocità angolare orbitale
ro = 10E-42; %Densità dell'atmosfera
Cd = 2.20; %Coefficiente di attrito spacecraft
S = 1.2^2; %Area sezione frontale spacecraft
n = 2; %Costante
Tmax = 1; %Impulso massimo propulsori
Cx = 0.005; %Costante positiva

%% Valori modificabili
ka = 1E-10; %Costante campo attrattivo

ki = 1E12; %Costante campo repulsivo

eta_0i1 = 350; %Safety radius
eta_0i2 = 100;

gamma = 2; %Costante campo iperbolico

%% Posizione e velocità iniziali
P = [-16100;
    0;
    3000;];

V = [5;
    0;
    0;];

Z = [       2*w*V(3,1);
    -w*P(2,1);
    -2*w*V(1,1) + 3*(w^2)*P(3,1)];

%% Posizione desiderata del target e definizione ostacoli
Pd = [-3500;
    0;
    0;]; %Posizione desiderta totale

obs1 = [-10000;
    0;
    1500;] ; % Vettore posizione dell'ostacolo

rad1=650;% raggio ostacolo
safety1=rad1+eta_0i1;

obs2= [-2800;
    0;
    0;];

rad2=150;% raggio ostacolo
safety2=rad2+eta_0i2;

obs3= [-2100;
    0;
    200;];

rad3=50;
safety3=rad3+eta_0i2;

Vdmax = 6; %Velocità massima desiderata

%% Forza attrattiva
Fa = -(ka*(P(:,1) - Pd(:,1))); %Forza attrattiva verso il target

%% Forze repulsive
%% Ostacolo 1
eta_i1=norm(P(:,1)-obs1(:,1));

Fr1 = ( (ki/(eta_i1.^2)) * ((1/eta_i1(:,1) - 1/safety1)) ) * ( (P(:,1) - obs1(:,1)) / norm(P(:,1) - obs1(:,1)));
if eta_i1 > safety1
    Fr1(:, 1) = 0;
end

%% Ostacolo 2
eta_i2=norm(P(:,1)-obs2(:,1));

Fr2 = ( (ki/(eta_i2.^2)) * ((1/eta_i2(:,1) - 1/safety2)) ) * ( (P(:,1) - obs2(:,1)) / norm (P(:,1) - obs2(:,1)));
if eta_i2 > safety2
    Fr2(:, 1) = 0;
end

%% Ostacolo 3
eta_i3=norm(P(:,1)-obs3(:,1));

Fr3 = ( (ki/(eta_i3.^2)) * ((1/eta_i3(:,1) - 1/safety3)) ) * ( (P(:,1) - obs3(:,1)) / norm (P(:,1) - obs3(:,1)));
if eta_i3 > safety3
    Fr3(:, 1) = 0;
end

%% Campo totale
Etot = Fa + Fr1 + Fr2 + Fr3;
Vd=Vdmax*Etot;

%% Forze
Fext = [-(1/2)*ro*V(1,1)^2*S*Cd;
    0;
    0]; %Forza di attrito esterna

sigma = V(:,1) - Vd(:,1) + Cx*(P(:,1) - Pd(:,1));

Fthr = -m*eye(3)*n*Tmax*sign(sigma); %Forza propulsori

Ftot = Fext + Fthr;%Forza totale

%% Simulazione numerica
for t = 1:N-1

    Fa(:, t+1) = -(ka.*(P(:,t) - Pd(:,1))); %Forza attrattiva verso il target

    %% Forze repulsive
    %% Ostacolo 1
    eta_i1(:,t+1) = norm(P(:,t)-obs1(:,1));

    Fr1(:,t+1) = ( (ki/(eta_i1(:,t).^2)) * ((1/eta_i1(:,t) - 1/safety1)) ) * ( (P(:,t) - obs1(:,1)) / ( norm(P(:,t)-obs1) ));
    if eta_i1(:,t) > safety1
        Fr1(:,t+1) = 0;
    end

    %% Ostacolo 2
    eta_i2(:,t+1)=norm(P(:,t)-obs2(:,1));

    Fr2(:,t+1) = ( (ki/(eta_i2(:,t).^2)) * ((1/eta_i2(:,t) - 1/safety2)) ) * ( (P(:,t) - obs2(:,1)) / ( norm(P(:,t)-obs2) ));
    if eta_i2(:,t) > safety2
        Fr2(:,t+1) = 0;
    end

    %% Ostacolo 3
    eta_i3(:,t+1)=norm(P(:,t)-obs3(:,1));

    Fr3(:,t+1) = ( (ki/(eta_i3(:,t).^2)) * ((1/eta_i3(:,t) - 1/safety3)) ) * ( (P(:,t) - obs3(:,1)) / ( norm(P(:,t)-obs3) ));
    if eta_i3(:,t) > safety3
        Fr3(:,t+1) = 0;
    end

    %% Campo totale
    
    Etot(:,t+1) = Fa(:, t) + Fr1(:,t) + Fr2(:,t) + Fr3(:,t);

    Vd(:,t+1)=Vdmax*Etot(:,t);

    Fext(1, t+1) = -(1/2)*ro*V(1,t)^2*S*Cd;

    sigma(:, t+1) = V(:,t) - Vd(:,t) + Cx*(P(:,t) - Pd(:,1));

    Fthr(:, t+1) = -m*eye(3)*n*Tmax*sign(sigma(:,t)); %Forza propulsori

    Ftot(:, t+1) = Fext(:, t) + Fthr(:, t);

    if norm( P(:, t) - Pd) < 10
        Pd = [-500;0;0;];
    end 
    
    Z(:,t+1) = [     2*w*V(3,t) ;
        -w*P(2,t);
        -2*w*V(1,t) + 3*(w^2)*P(3,t)];

    P(:, t+1) = P(:, t) + h*V(:,t);

    V(:, t+1) = V(:, t) + h*Z(:, t) + h*(Ftot(:,t)/m);
    %
    %     if norm(P(:,t+1)-obs1)<rad1 || norm(P(:,t+1)-obs2)<rad2 || norm(P(:,t+1)-obs3)<rad3
    %
    %         disp('l ostacolo è stato colpito');
    %         break;
    %     end
%     if norm(P(:,t)-Pd)<25
%         disp('Target raggiunto');
%         break;
%     end

end

t = linspace(0, tf, t+1);

% figure(1)
% 
% %% Plot safety radius dei tre ostacoli
% %%obs1
% plot(obs1(1,1), obs1(3,1), 'bo')
% hold on
% %%radius1
% a=rad1; % horizontal radius
% b=safety1; % vertical radius
% xRadius=obs1(1,1); % x0,y0 ellipse centre coordinates
% yRadius=obs1(3,1);
% g = -pi:0.01:pi;
% x=xRadius+a*cos(g);
% y=yRadius+b*sin(g);
% plot(x,y, 'r--')
% hold on
% 
% %%obs2
% plot(obs2(1,1), obs2(3,1), 'bo')
% hold on
% 
% a=rad2; % horizontal radius
% b=safety2; % vertical radius
% xRadius=obs2(1,1); % x0,y0 ellipse centre coordinates
% yRadius=obs2(3,1);
% g = -pi:0.01:pi;
% x=xRadius+a*cos(g);
% y=yRadius+b*sin(g);
% plot(x,y, 'r--')
% hold on
% 
% %%obs3
% plot(obs3(1,1), obs3(3,1), 'bo')
% hold on
% 
% a=rad3; % horizontal radius
% b=safety3; % vertical radius
% xRadius=obs3(1,1); % x0,y0 ellipse centre coordinates
% yRadius=obs3(3,1);
% g = -pi:0.01:pi;
% x=xRadius+a*cos(g);
% y=yRadius+b*sin(g);
% plot(x,y, 'r--')
% grid on
% 
% %Posizione desiderata
% plot(Pd(1,1), Pd(3,1), 'ro')
% 
% for k = 1:length(t)
%     plot(P(1,1:k), P(3,1:k), 'k')
%     axis([-17000, 10000, -4000, 4000])
%     pause(0.01)
% end

figure(2)

plot(P(1,:), P(3,:), 'k')
hold on
plot(Pd(1,1), Pd(3, 1), 'ro')
hold on
xlabel('Vbar')
ylabel('Rbar')
%% Plot safety radius dei tre ostacoli
%obs1
plot(obs1(1,1), obs1(3,1), 'bo')
hold on

%radius1
a=650; % horizontal radius
b=safety1; % vertical radius
xRadius=obs1(1,1); % x0,y0 ellipse centre coordinates
yRadius=obs1(3,1);
g = -pi:0.01:pi;
x=xRadius+a*cos(g);
y=yRadius+b*sin(g);
plot(x,y, 'r--')
hold on

%obs2
plot(obs2(1,1), obs2(3,1), 'bo')
hold on

%radius2
a=rad2; % horizontal radius
b=eta_0i2; % vertical radius
xRadius=obs2(1,1); % x0,y0 ellipse centre coordinates
yRadius=obs2(3,1);
g = -pi:0.01:pi;
x=xRadius+a*cos(g);
y=yRadius+b*sin(g);
plot(x,y, 'r--')
hold on

%obs3
plot(obs3(1,1), obs3(3,1), 'bo')
hold on

%%radius3
a=rad3; % horizontal radius
b=eta_0i2; % vertical radius
xRadius=obs3(1,1); % x0,y0 ellipse centre coordinates
yRadius=obs3(3,1);
g = -pi:0.01:pi;
x=xRadius+a*cos(g);
y=yRadius+b*sin(g);
plot(x,y, 'r--')
grid on

figure(3)

tiledlayout(2,1)

nexttile
plot(t, V(1,:))
xlabel('Time[s]')
ylabel('Vx[m/s]')
grid on
axis([0 2600 -20 60])

nexttile
plot(t, V(3,:))
xlabel('Time[s]')
ylabel('Vz[m/s]')
grid on
axis([0 2600 -20 50])

figure(4)
tiledlayout(2,1)

nexttile
plot(t, Fthr(1,:))
xlabel('Time[s]')
ylabel('Fx[N]')
grid on
axis([0 2600 -1300 1300])

nexttile
plot(t, Fthr(3,:))
xlabel('Time[s]')
ylabel('Fz[N]')
grid on
axis([0 2600 -1300 1300])