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
Pd = [-3000;
    0;
    0]; %Vettore posizione deisderata

Vdmax = 6;%Velocità desiderata massima

obs = [-10000;
    0;
    1500]; %Vettore posizione dell'ostacolo

%% Matrici posizione e velocità iniziali
P = [-16100;
    0;
    3000]; %Vettore posizione, che introduce la posizione iniziale dello spacecraft lungo gli assi.

V = [5;
    0;
    0]; %Vettore velocità dello spacecraft, che introduce le velocità iniziali.

Z = [       2*w*V(3,1);
    -w*P(2,1);
    -2*w*V(1,1) + 3*(w^2)*P(3,1)];

%% Condizioni campi potenziali
ka = 1E-10;
ki = 1E12;

gamma = 2;

eta_0i = 350; %Distanza da mantenere rispetto all'ostacolo
rad = 650;
safety = eta_0i + rad;

eta_i = norm(P(:,1) - obs(:,1)); %Distanza dall'ostacolo

%% Potenziale attrattivo verso il target
Ua = (1/2)*ka*((P(1,1) - Pd(1,1))^2 + (P(2,1) - Pd(2,1))^2 + (P(3,1) - Pd(3,1))^2); %Campo potenziale attrattivo

Fa = -(ka*(P(:,1) - Pd(:,1))); %Forza attrattiva verso il target

Eu = Fa(:, 1) / (ka*norm(P(:,1) - Pd(:,1)));

%% Potenziale repulsivo ostacoli
Ur = (ki/gamma).*((1/eta_i(:,1) - 1/safety).^gamma); % Campo potenziale repulsivo
if (eta_i(:, 1) > eta_0i)
    Ur(:, 1) = 0;
end

Fr = (ki/(eta_i^2))*((1/eta_i(:,1) - 1/(safety)).^(gamma-1))*((P(:,1) - obs(:,1))/norm(P(:,1) - obs(:,1)));
if (eta_i(:, 1) > eta_0i)
    Fr(:, 1) = 0;
end

%% Potenziale totale
Ut = (1/2)*ka*((P(1,1) - Pd(1,1))^2 + (P(2,1) - Pd(2,1))^2 + (P(3,1) + Pd(3,1))^2) + (ki/gamma).*((1/eta_i(:,1) - 1/eta_0i).^gamma);

E = Fa + Fr;

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
    eta_i(:, t+1) = norm(P(:,t) - obs(:,1));

    Ur(:, t+1) = (ki/gamma)*((1/(eta_i(:,t)) - 1/(safety)).^gamma);
    if eta_i(:, t+1) > safety
        Ur(:, t) = 0;
    end

    Fr(:, t+1) = ((ki/(eta_i(:,t)^2))*((1/eta_i(:,t) - 1/(safety)).^(gamma-1)))*((P(:,t) - obs(:,1)) / norm(P(:,t) - obs(:,1)));
    if eta_i(:, t) > safety
        Fr(:, t+1) = 0;
    end

    %% Vettore per la direzione e velocità desiderata 
    E(:, t+1) = Fa(:, t) + Fr(:, t);
    
    Vd(:, t+1) = Vdmax*E(:,t);

    %% Forza esterna e propulsori
    Fext(1, t+1) = -(1/2)*ro*(V(1,t)^2)*S*Cd;

    F(:,t+1) = Fthr(:,t) + Fext(:, t); %Forza totale

    sigma(:,t+1) = V(:,t) - Vd(:,t) + Cx*(P(:,t) - Pd(:,1)); %Propulsori accesi/spenti

    Fthr(:,t+1) = -m*eye(3)*n*Tmax*sign(sigma(:,t)); %Forza propulsori

    Z(:,t+1) = [ 2*w*V(3,t) ;
        -w*P(2,t);
        -2*w*V(1,t) + 3*(w^2)*P(3,t)];

    %% Posizione e velocità che evolvono secondo il metodo di Eulero in avanti 
        
    P(:,t+1) = P(:,t) + h*V(:,t); %Vettore posizione

    V(:,t+1) = V(:,t) + h*Z(:,t) + h*(F(:,t)/m); %Vettore velocità
   
    if norm(P(:,t)-obs)<rad 
        disp('L ostacolo è stato colpito');
        break;
    end

end

t = linspace(0, tf, t+1);

figure(1)
plot (P(1,:), P(3,:), 'k')
hold on
plot(P(1,1), P(3,1), 'bo')
hold on
plot(obs(1,1), obs(3,1), 'ro')
hold on
plot(Pd(1,1), Pd(3,1), 'mo')
hold on

%% Plot safety radius dall'ostacolo
a=rad; % horizontal radius
b=safety; % vertical radius
xRadius=obs(1,1); % x0,y0 ellipse centre coordinates
yRadius=obs(3,1);
g = -pi:0.01:pi;
x=xRadius+a*cos(g);
y=yRadius+b*sin(g);
plot(x,y, 'r--')
grid on

xlabel('Vbar')
ylabel('Rbar')

figure(2)

tiledlayout(2,1)

nexttile
plot(t, V(1,:))
xlabel('Time[s]')
ylabel('Vx[m/s]')
axis([0 2600 -20 60])

nexttile
plot(t, V(3,:))
xlabel('Time[s]')
ylabel('Vz[m/s]')
axis([0 2600 -20 50])

figure(3)
tiledlayout(2,1)

nexttile
plot(t, Fthr(1,:))
xlabel('Time[s]')
ylabel('Fx[N]')
axis([0 2600 -1300 1300])

nexttile
plot(t, Fthr(3,:))
xlabel('Time[s]')
ylabel('Fz[N]')
axis([0 2600 -1300 1300])