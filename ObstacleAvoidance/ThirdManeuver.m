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
tf = 4400; %Tempo di esecuzione del programma
N = 140400; %Passo
h = tf/N; %Passo numerico

%% Posizione e velocità desiderata
Pd = [0;
    0;
    0]; %Vettore posizione deisderata

Vdmax = 6;%Velocità desiderata massima

%% Matrici posizione e velocità iniziali
P = [-500;
    0;
    0]; %Vettore posizione, che introduce la posizione iniziale dello spacecraft lungo gli assi.

V = [0.15;
    0;
    0]; %Vettore velocità dello spacecraft, che introduce le velocità iniziali.

Z = [       2*w*V(3,1);
    -w*P(2,1);
    -2*w*V(1,1) + 3*(w^2)*P(3,1)];

y=-1:1; %Estremi del cono
a=-(abs(y)/0.002) + 25; %Funzione per plot di un cono

%% Condizioni campi potenziali
ka = 1E-10;

%% Potenziale attrattivo verso il target
Ua = (1/2)*ka*((P(1,1) - Pd(1,1))^2 + (P(2,1) - Pd(2,1))^2 + (P(3,1) - Pd(3,1))^2); %Campo potenziale attrattivo

Fa = -(ka*(P(:,1) - Pd(:,1))); %Forza attrattiva verso il target

Eu = Fa(:, 1) / (ka*norm(P(:,1) - Pd(:,1)));

%% Potenziale totale
Ut = (1/2)*ka*((P(1,1) - Pd(1,1))^2 + (P(2,1) - Pd(2,1))^2 + (P(3,1) + Pd(3,1))^2);

E(:, 1) = Fa(:, 1);

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

    %% Vettore per la direzione e velocità desiderata 
    E(:, t+1) = Fa(:, t);
    
    Vd(:, t+1) = Vdmax*Eu(:,t);
    %% Forza esterna e propulsori
    Fext(1, t+1) = -(1/2)*ro*(V(1,t)^2)*S*Cd;

    F(:,t+1) = Fthr(:,t) + Fext(:, t); %Forza totale

    sigma(:,t+1) = V(:,t) - Vd(:,t) + Cx*(P(:,t) - Pd(:,1)); %Propulsori accesi/spenti

    Fthr(:,t+1) = -m*eye(3)*n*Tmax*sign(sigma(:,t)); %Forza propulsori

    Z(:,t+1) = [ 2*w*V(3,t) ;
        -w*P(2,t);
        -2*w*V(1,t) + 3*(w^2)*P(3,t)];

    %% Posizione e velocità che evolvono secondo il metodo di Eulero in avanti 
    P(:,t+1) = (P(:,t) + h*V(:,t)); %Vettore posizione

    V(:,t+1) = V(:,t) + h*Z(:,t) + h*(F(:,t)/m); %Vettore velocità
   
     if norm(P(:,t) - Pd(:,1)) < 0.5
        disp('Target raggiunto');
        break;
     end
    
end

t = linspace(0, tf, t+1);

% figure(1)
% 
% for k = 1:length(t)
%     plot(P(1,1:k), P(3,1:k))
%     axis([-500, 0, -1, 1])
%     pause(0.0001)
% end

figure(2)

nexttile
plot (P(1,:), P(3,:), 'k')
hold on
plot(P(1,1), P(3,1), 'bo')
hold on
plot(Pd(1,1), Pd(3,1), 'ro')
xlabel('Vbar')
ylabel('Rbar')

axis([-500 0  -1 1  0 2*pi])
hold on

plot(a,y, 'r--')
grid on

figure(3)

tiledlayout(2,1)

nexttile
plot(t, V(1,:))
xlabel('t[s]')
ylabel('Vx[m/s]')
grid on
axis([0 4500 0 15])

nexttile
plot(t, V(3,:))
xlabel('t[s]')
ylabel('Vz[m/s]')
grid on
axis([0 4500 -1 1])

figure(4)
tiledlayout(2,1)

nexttile
plot(t, Fthr(1,:))
xlabel('Time[s]')
ylabel('Fx[N]')
grid on
axis([0 4500 -1300 1300])

nexttile
plot(t, Fthr(3,:))
xlabel('Time[s]')
ylabel('Fz[N]')
grid on
axis([0 4500 -1300 1300])