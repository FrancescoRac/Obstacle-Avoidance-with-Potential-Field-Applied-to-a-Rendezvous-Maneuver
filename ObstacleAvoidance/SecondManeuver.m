clear;
clc;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultAxesFontSize',14)
%% Condizioni Operative
w = 0.00113; %Velocità angolare
m = 600; %Massa dello spacecraft

ro = 10E-42; %Densità dell'atmosfera
Cd = 2.20; %Coefficiente di attrito dello spacecraft
S = 1.2^2; %Sezione frontale dello spacecraft
Cx = 0.005; %Costante positiva
Tmax = 1; %Forza propulsori
n = 2; %Propulsori accesi simultaneamente

%% Definizione del passo numerico
tf = 2500; %Tempo di esecuzione del programma
N = 80000; %Passo
h = tf/N; %Passo numerico

%% Posizione e velocità desiderata
Pd = [-500;
    0;
    0]; %Vettore posizione deisderata

Vdmax = 6;%Velocità desiderata massima

obs1 = [-2800;
    0;
    0]; %Vettore posizione del primo ostacolo

rad1 = 150; %Radius del primo ostacolo

obs2 = [-2100;
    0;
    200;]; %Vettore posizione del secondo ostacolo

rad2 = 50; %Radius del secondo ostacolo

%% Matrici posizione e velocità iniziali
P = [-3500;
    0;
    0]; %Vettore posizione, che introduce la posizione iniziale dello spacecraft lungo gli assi.

V = [0;
    0;
    0.5]; %Vettore velocità dello spacecraft, che introduce le velocità iniziali.

Z = [       2*w*V(3,1);
    -w*P(2,1);
    -2*w*V(1,1) + 3*(w^2)*P(3,1)];

%% Condizioni campi potenziali
ka = 1E-10; %Costante campo attrattivo
ki = 1E12; %Costante campo repulsivo

gamma = 2;
eta_0i = 100; %Distanza da mantenere rispetto all'ostacolo

eta_i1 = norm(P(:,1) - obs1(:,1)); %Distanza dall'ostacolo numero 1
safety1 = eta_0i + rad1; %Safety radius del primo ostacolo

eta_i2 = norm(P(:,1) - obs2(:,1)); %Distanza dall'ostacolo numero 2
safety2 = eta_0i + rad2; %Safety radius del secondo ostacolo

%% Potenziale attrattivo verso il target
Ua = (1/2)*ka*((P(1,1) - Pd(1,1))^2 + (P(2,1) - Pd(2,1))^2 + (P(3,1) - Pd(3,1))^2); %Campo potenziale attrattivo

Fa = -(ka*(P(:,1) - Pd(:,1))); %Forza attrattiva verso il target

Eu = Fa(:, 1) / (ka*norm(P(:,1) - Pd(1,1)));

%% Potenziale repulsivo ostacoli
Ur1 = (ki/gamma).*((1/eta_i1(:,1) - 1/safety1).^gamma); % Campo potenziale repulsivo primo ostacolo
if (eta_i1 > safety1)
    Ur1 = 0;
end

Ur2 = (ki/gamma).*((1/eta_i2(:,1) - 1/safety2).^gamma); % Campo potenziale repulsivo secondo ostacolo
if (eta_i2 > safety2)
    Ur2 = 0;
end

Fr1 = (ki/(eta_i1^2))*(( 1 / (eta_i1(:, 1)) - 1/(safety1))^(gamma-1))*((P(:,1) - obs1(:,1)) / norm(P(:,1) - obs1(:,1))); %Forza repulsiva primo ostacolo
if (eta_i1 > safety1)
    Fr1(:, 1) = 0;
end

Fr2 = (ki/(eta_i2^2))*(( 1 / (eta_i2(:, 1)) - 1/(safety2))^(gamma-1))*((P(:,1) - obs2(:,1)) / norm(P(:,1) - obs2(:,1))); %Forza repulsivia secondo ostacolo
if (eta_i2 > safety2)
    Fr2(:, 1) = 0;
end

%% Potenziale totale
Ut = (1/2)*ka*((P(1,1) - Pd(1,1))^2 + (P(2,1) - Pd(2,1))^2 + (P(3,1) + Pd(3,1))^2) + (ki/gamma).*((1/eta_i1(:,1) - 1/eta_0i).^gamma) + (ki/gamma).*((1/eta_i2(:,1) - 1/eta_0i).^gamma);

E = Fa + Fr1 + Fr2;

Vd = Vdmax*E(:,1); %Velocità desiderata

%% Forze in gioco
Fext = [-(1/2)*ro*(V(1,1)^2)*S*Cd;
    0;
    0]; %Attrito dell'atmosfera applicato solo lungo l'asse x.

sigma = (V(:,1) - Vd(:,1)) + Cx*(P(:,1) - Pd(:,1)); %%Problema causato dal versore perchè restituisce un NaN come valore lungo l'asse y
Fthr = -m*eye(3)*n*Tmax*sign(sigma);
F = Fext + Fthr; %Vettore delle forze lungo gli assi

%% Simulazione numerica
for t = 1:N-1
    %% Potenziale attrattivo e velocità desiderata
    Ua(:, t+1) = (1/2)*ka*((P(1,t) - Pd(1,1))^2 + (P(2,t) - Pd(2,1))^2 + (P(3,t) - Pd(3,1))^2);

    Fa(:, t+1) = -(ka*(P(:,t) - Pd(:,1)));

    Eu(:, t+1) = Fa(:, t) / (ka*norm(P(:,t) - Pd(:,1)));

    %% Potenziale repulsivo
    
    %Ostacolo numero 1 
    eta_i1(:, t+1) = norm(P(:,t) - obs1(:,1)); %Distanza dall'ostacolo

    Ur1(:, t+1) = (ki/gamma)*((1/(eta_i1(:,t)) - 1/(safety1)).^gamma);
    if (eta_i1(:, t) > safety1) 
        Ur1(:, t+1) = 0;        
    end

    Fr1(:, t+1) = ((ki/(eta_i1(:,t)^2))*((1/eta_i1(:,t) - 1/safety1)^(gamma-1)))*((P(:,t) - obs1(:,1)) / norm(P(:,t) - obs1(:,1)));
    if eta_i1(:, t) > safety1
        Fr1(:, t+1) = 0;
    end

    %Ostacolo numero 2
    eta_i2(:, t+1) = norm(P(:,t) - obs2(:,1)); %Distanza dall'ostacolo

    Ur2(:, t+1) = (ki/gamma)*((1/(eta_i2(:,t)) - 1/(safety2)).^gamma);
    if (eta_i2(:, t) > safety2) 
        Ur2(:, t+1) = 0;        
    end

    Fr2(:, t+1) = ((ki/(eta_i2(:,t)^2))*((1/eta_i2(:,t) - 1/safety1)^(gamma-1)))*((P(:,t) - obs2(:,1)) / norm(P(:,t) - obs2(:,1)));
    if eta_i2(:, t) > safety2
        Fr2(:, t+1) = 0;
    end

    %% Potenziale totale e versore per la direzione
    Ut(:, t+1) = (1/2)*ka*((P(1,t) - Pd(1,1))^2 + (P(2,t) - Pd(2,1))^2 + (P(3,t) - Pd(3,1))^2) + (ki/gamma)*((1/eta_i1(:,t) - 1/eta_0i).^gamma) + (ki/gamma)*((1/eta_i2(:,t) - 1/eta_0i).^gamma);
    
    E(:, t+1) =  Fa(:, t) + Fr1(:, t) + Fr2(:, t);
    
    Vd(:, t+1) = Vdmax*E(:,t);
    %% Forze
    Fext(1, t+1) = -(1/2)*ro*(V(1,t)^2)*S*Cd;

    F(:,t+1) = Fthr(:, t) + Fext(:, t); %Forza totale

    sigma(:, t+1) = V(:, t) - Vd(:, t) + Cx*(P(:, t) - Pd(:, 1)); %Propulsori accesi/spenti

    Fthr(:,t+1) = -m*eye(3)*n*Tmax*sign(sigma(:,t)); %Forza propulsori

    Z(:,t+1) = [ 2*w*V(3,t) ;
        -w*P(2,t);
        -2*w*V(1,t) + 3*(w^2)*P(3,t)];

    %% Posizione e velocità
    P(:,t+1) = P(:,t) + h*V(:,t); %Vettore posizione

    V(:,t+1) = V(:,t) + h*Z(:,t) + h*(F(:,t)/m); %Vettore velocità

    if norm(P(:,t)-obs1)<rad1 || norm(P(:,t)-obs2)<rad2 
        disp('L ostacolo è stato colpito');
        break;
    end
    
end

t = linspace(0, tf, t+1);

figure(1)

plot (P(1,:), P(3,:), 'k')
hold on
plot(obs1(1,1), obs1(3,1), 'ro')
hold on
plot(obs2(1,1), obs2(3,1), 'ro')
hold on
plot(Pd(1,1), Pd(3,1), 'bo')
xlabel('Vbar')
ylabel('Rbar')
hold on

%% Plot safety radius dall'ostacolo
a=rad1; % horizontal radius
b=eta_0i; % vertical radius
xRadius=obs1(1,1); % x0,y0 ellipse centre coordinates
yRadius=obs1(3,1);
g = -pi:0.01:pi;
x=xRadius+a*cos(g);
y=yRadius+b*sin(g);
plot(x,y, 'r--')
hold on

a=rad2; % horizontal radius
b=eta_0i; % vertical radius
xRadius=obs2(1,1); % x0,y0 ellipse centre coordinates
yRadius=obs2(3,1);
f = -pi:0.01:pi;
x=xRadius+a*cos(f);
y=yRadius+b*sin(f);
plot(x,y, 'r--')
grid on

figure(2)

tiledlayout(2,1)

nexttile
plot(t, V(1,:))
xlabel('Time[s]')
ylabel('Vx[m/s]')
grid on
axis([0 2600 -10 18])

nexttile
plot(t, V(3,:))
xlabel('Time[s]')
ylabel('Vz[m/s]')
grid on
axis([0 2600 -18 22])

figure(3)
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