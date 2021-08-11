% Simulador de telas inextensibles 
% Franco Coltraro, IRI, CSIC-UPC, Barcelona
% 2019/10/8
close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETROS SIMULACION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tiempo
tf = 1000; dt = 0.015;
times = dt*(0:tf);
%parametros de la superficie
rho = 0.2; %densidad de masa
theta = 0.0; %resistencia  doblarse: rigidez;
alfa = 1.5;  %amortiguamiento oscilaciones lentas (~rozamiento con el aire)
beta = 1.5; %amortiguamiento oscilaciones rapidas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MALLADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%se crea a partir de una parametrizacion de la forma (f1(x,y),f2(x,y),f3(x,y))
%donde (x,y) estan en el rectangulo [ax,bx]x[ay,by]
ax = -0.15; bx = 0.15; ay = -0.15; by = 0.15; 
f1 = @(x,y) x; f2 = @(x,y) 0.75*sqrt(1 - x.^2) -0.1; f3 = @(x,y) y;
npx = 10; npy = 10; n_nodos = npx*npy; %numero de puntos
[X,T] = CreateMesh([ax,bx,ay,by],npx,npy,f1,f2,f3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOUNDARY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deteccion de los puntos del borde para poder imponer isometria tambien
%alli
[Xb,Tb,nodos_borde] = GetBoundary(X,[ax,bx,ay,by],npx,npy,f1,f2,f3); 
n_nodos_bnd = size(Xb.Xb,1);
%interior nodes
nodes_int = setdiff(1:n_nodos,nodos_borde.nodes_bnd); 
%coordenadas de las esquinas
esquinas = [1, npx, npx*(npy-1)+1, npx*npy]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%VARIABLES DE ESTADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Superficie inicial: ver CreateMesh
phi0 = X; 
%Velocidad inicial
dphi0 = [zeros([n_nodos,1]),zeros([n_nodos,1]),zeros([n_nodos,1])]; 
%%%%%%%%%%%%%%%%%%%%%%%%%VARIABLES DE ESTADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%NODOS A CONTROLAR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%damos las coordenadas de los nodos de los que vamos a fijar su trayectoria
controlados = [npx*(npy-1)+1, npx*npy]; %dos esquinas
coord_controlados = [controlados, controlados+n_nodos, controlados+(2*n_nodos)];
%matriz para imponer las condiciones en el contorno
A_b = spalloc(length(coord_controlados),3*n_nodos,length(coord_controlados));
A_b(:,coord_controlados) = eye(length(coord_controlados));

%posicion de los nodos para el control
u = zeros(length(coord_controlados),tf+1);
u(:,1) = reshape(X(controlados,:),[length(coord_controlados),1]);
f = 2;
for tt=2:(tf+1)
    if (tt)>=0 && (tt) <=50
        u(:,tt) = u(:,tt-1);
    elseif (tt)>50 && (tt) <= 350
        u(:,tt) = u(:,tt-1);
        u(3,tt) = u(3,50) + 0.01*sin(2*pi*f*(tt-50)/400);
        u(4,tt) = u(4,50) + 0.01*sin(2*pi*f*(tt-50)/400);
%         u(5,tt) = u(5,tt-1) + 0.0005;
%         u(6,tt) = u(6,tt-1) + 0.0005;
    elseif (tt)>500 && (tt) <= 950
        r1 = 0.1;
        r2 = 0.4;
        angle = (pi/4)*(tt-500)/450;
        u(:,tt) = u(:,tt-1);
        u(1,tt) = u(1,500) + r2-r2*cos(angle);
        u(2,tt) = u(2,500) + r1-r1*cos(angle);
        u(3,tt) = u(3,500) - r2*sin(angle);
        u(4,tt) = u(4,500) - r1*sin(angle);
        u(5,tt) = u(5,tt-1) - 0.0001;
        u(6,tt) = u(6,tt-1) - 0.0001;
%     elseif (tt)>750 && (tt) <= 950
%         r1 = 0.35;
%         r2 = 0.05;
%         angle = (pi/6)*(tt-750)/200;
%         u(:,tt) = u(:,tt-1);
%         u(1,tt) = u(1,750) + r2-r2*cos(angle);
%         u(2,tt) = u(2,750) + r1-r1*cos(angle);
%         u(3,tt) = u(3,750) - r2*sin(angle);
%         u(4,tt) = u(4,750) - r1*sin(angle);
%         u(5,tt) = u(5,tt-1) + 0.0001;
%         u(6,tt) = u(6,tt-1) + 0.0001;
    else
        u(:,tt) = u(:,tt-1);
    end
    
    store_du(:,tt) = abs(u(:,tt) - u(:,tt-1));
%     tt
%     abs(sqrt( (u(1,tt)-u(2,tt))^2 + (u(3,tt)-u(4,tt))^2 +(u(5,tt)-u(6,tt))^2) )
end
max(store_du(:))
figure
plot(u')
legend('1','2','3','4','5','6')
%%
% u=load('uLaw.csv');
%u=u';

% nodos1 = T(1,:);
% nodos2 = T(T(:,2)==npx,:);
% %matriz para imponer las condiciones en el contorno
% A_b = spalloc(6,3*n_nodos,24);
% unos = 0.25*ones([1 length(nodos1)]);
% A_b(1:3,[nodos1,nodos1+n_nodos,nodos1+(2*n_nodos)]) = blkdiag(unos,unos,unos);
% A_b(4:6,[nodos2,nodos2+n_nodos,nodos2+(2*n_nodos)]) = blkdiag(unos,unos,unos);
% %posicion de los nodos para el control
% u = zeros(6,tf+1);
% u(:,1) = A_b*X(:);
% v = [0;-1;0;0;-1;0]; freq = 1;
% for tt=2:(tf+1)
%     if times(tt) < 4
%         u(:,tt) = u(:,tt-1);
%     elseif times(tt) >= 4 && times(tt) <= 6
%         u(:,tt) = u(:,tt-1) + 0.015*cos(2*pi*freq*times(tt)+0.1)*v;% + 0.008*randn([6,1]);
%     else
%         u(:,tt) = u(:,tt-1);
%     end 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%ELEMENTOS DE REFERENCIA%%%%%%%%%%%%%%%%%%%%%%%%%
%cuadrilateros bilineales
theReferenceElement = createReferenceElement();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%SYSTEM OF EQUATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculo de las matrices del sistema: masa Mlum Minv, rigidez K, 
%amortigualmiento D e inextensibilidad C (3-tensor)
[C,n_conds,Mlum,Minv,Mcons,D,K] = computeMatrices(X,T,...
                               nodes_int,nodos_borde,esquinas,...
                               [theta,alfa,beta],...
                               theReferenceElement);
[Cphi0,~] = fun_C(X,C,Mcons,A_b,n_nodos,n_conds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
%%%%%%%%%%%%%%%%%%%%GRAVEDAD Y FUERZAS EXTERNAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fg = Mlum*reshape([zeros([n_nodos,1]),...
                   zeros([n_nodos,1]),...
          -9.8*rho*ones([n_nodos,1])],[3*n_nodos,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%SUPER BUCLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%almacenamiento de las trajectorias
phiPositions = {};  phiPositions{1} = phi0;
phiVelocities = {}; phiVelocities{1} = dphi0;

%para iniciar el integrador de segundo orden
[phi1,dphi1] = explicit_grispun(dt,phi0,dphi0,...
                Fg,rho,A_b,Mlum,Minv,Mcons,K,D,C,...
                n_conds,n_nodos,Cphi0,u(:,2));
            
phiPositions{2} = reshape(phi1,[n_nodos,3]);
phiVelocities{2} = reshape(dphi1,[n_nodos,3]);

cond_initial = phi1;
% csvwrite('xini_nou.csv',cond_initial);
%%
tic
for tt=1:(tf-1)  
%     disp(tt)
    %resolvemos las EDO
    [phi,dphi] = imex_franco(dt,phi0,dphi0,phi1,dphi1,...
                           Fg,rho,A_b,Minv,Mcons,K,D,C,...
                           n_conds,n_nodos,Cphi0,u(:,tt+2));
    store_ref(:,tt) = phi([1 10 101 110 201 210]);
    phiPositions{tt+2} = reshape(phi,[n_nodos,3]);
    phiVelocities{tt+2} = reshape(dphi,[n_nodos,3]);   
    %actualizamos
    phi0 = phi1; dphi0 = dphi1;
    phi1 = phi;  dphi1 = dphi;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(store_ref')
%csvwrite('phi1_3d.csv',store_ref([1 3 5],:)')
%csvwrite('phi4_3d.csv',store_ref([2 4 6],:)')


%%%%%%%%%%%%%%%%%%%%%%%POSTPROCESO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pelicula
figure
peli = ''; vel_vid = 1; 
hacerPeli(phiPositions,T,nodos_borde,peli,vel_vid)
% [Xtri,Ttri] =  traingulateQuadMesh(phiPositions{400},T); 
% [super,curve] = plotMesh(Xtri,Ttri,nodos_borde);

%Area
areaError = computeAreaError(phiPositions,T,theReferenceElement);
figure
plot(times,areaError)
hold on
mu = mean(areaError);
hline = refline([0 mu]);
hline.Color = 'r';
hold off
xlabel('Time')
ylabel('Error')
title('Relative total area error')
disp('Error area medio:')
disp(mu)
ylim([0 max(1, 2*mu)])

%aristas
% aristasError = computeEdgesError(phiPositions,Tb.Tb);
% figure(3)
% plot(times,aristasError)
% hold on
% mu = mean(aristasError);
% hline = refline([0 mu]);
% hline.Color = 'r';
% hold off
% xlabel('Time')
% ylabel('Error')
% title('Relative boundary length error')
% ylim([0 10])

% %energia
% energias = computeEnergies(phiPositions,phiVelocities,rho,Mlum,Fg,K);
% figure(4)
% plot(times,energias)
% hold on
% mu = mean(energias);
% disp('Energia media:')
% disp(mu)
% hline = refline([0 mu]);
% hline.Color = 'r';
% xlabel('Time')
% ylabel('Energy')
% ylim([min(energias)-1 max(energias)+1])

% %guardamos info
% matrices = struct('Mlum',Mlum,'Minv',Minv,'K',(1/theta)*K,'C',C);
% save('manta.mat','phi','T','nodos_borde','matrices','theReferenceElement')

