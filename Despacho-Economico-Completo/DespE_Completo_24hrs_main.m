%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Despacho Economico Completo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objetivo: Realizar un despacho economico optimo en un sistema considerando la red de transmisiÃƒÂ³n y 
%           sus perdidas para un periodo de 24 horas con diferentes curvas de demanda

clear all
clc

%Se carga el caso de estudio
caso14n5g;

%Datos de tolerancia e iteraciones maximas para los metodos de Flujos de potencia y despacho economico.
tolerancia   = 1e-3;  %tolerancia para la convergencia
toleranciaNR = 1e-6;  %tolerancia para Newton Raphson 
iter_max     = 1000  ;  %iteraciones maximas antes de detener el calculo para ambos ciclos

%Se convierten los datos a en por unidad.
Nodos(:,5)       = Nodos(:,5)./Base_MVA;
Nodos(:,6)       = Nodos(:,6)./Base_MVA;
Generadores(:,2) = Generadores(:,2)./Base_MVA;
Generadores(:,3) = Generadores(:,3)./Base_MVA;

%Se hacen respaldos
Nodos_old  = Nodos;
Costos_old = Costos;
Nodos_In = Nodos;

%Se crean vectores con posiciones e indices necesarios
Slk_pos  = find(Nodos(1:size(Nodos,1),2) == 1); %Posicion del nodo Slack
Pgen_pos = find(round(Nodos(:,2)) < 3)        ; %Posicion de los nodos con generadores
PVQ_pos  = find(round(Nodos(:,2)) >=2)        ; %Posicion de los nodos PV y PQ
Num_Nod  = size(Nodos,1)                      ; %Numero de nodos
Num_Dem  = size(Demanda,1)                    ; %Numero de cargas
CT = zeros(1, 24);                              %Costos totales
P_marg = zeros(Num_Nod,24);                     %Precios marginales

for t=1:1:24
    Nodos = Nodos_In;  
    %Se reinicia la matriz de nodos a sus valores iniciales
    for k=1:Num_Dem
        Nodos(Demanda(k,1),5) = Nodos(Demanda(k,1),5)*Demanda(k,t+1);   %Se ajusta la potencia activa de acuerdo a la curva de demanda
        Nodos(Demanda(k,1),6) = Nodos(Demanda(k,1),6)*Demanda(k,t+1);   %Se ajusta la potencia reactiva de acuerdo a la curva de demanda
    end
    Nodos_old = Nodos;
    %Despacho economico simplificado considerando a la red para optener potencia generadas iniciales
    [Pgen] = DespE_LT(Lineas,Nodos,Base_MVA,Costos);

    %Se actualiza la potencia generada en la matriz Nodos
    Generadores(:,2) = Pgen(Pgen_pos); 

    %Se efectua un estudio de flujos de potencia 
    [V,P_per,Qgen,P_te,Pslk_te] = load_flow_NR(Lineas,Nodos,Generadores,toleranciaNR,iter_max);

    %Se actualizan las magnitudes de voltajes y angulos de los nodos en la matriz Nodos
    Nodos(:,4) = angle(V);
    Nodos(:,3) = abs(V);

    %Inicia proceso iterativo para el depsacho economico considerando perdidas
    conv_flag = 0;
    iter      = 0;
    while conv_flag==0 && iter <= iter_max
        P_old        = Pgen; %Se guarda la Pgen incial en una variable
        iter         = iter+1;
        Beta         = Pslk_te/P_te; %Se calculan las Beta

        %Se calculan los ITL
        ITL          = zeros(1,Num_Nod);
        ITL(PVQ_pos) = Beta+1; 

        %Se agregan las perdidas a las potencias activas demandadas
        PL           = zeros(Num_Nod,1)   ;
        PL(Slk_pos)  = P_per              ;
        Nodos(:,5)   = Nodos_old(:,5) + PL;

        %Se obtienen los nuevos valores de "b" y "d" para obtner nuevas curvas de costos
        Costos(:,2)=Costos(:,2)./(1-ITL(Pgen_pos)).';
        Costos(:,3)=Costos(:,3)./(1-ITL(Pgen_pos)).';

        %Despacho simplificado considerando a la red y las nuevas cuvas de costos
        [Pgen,la] = DespE_LT(Lineas,Nodos,Base_MVA,Costos);

        %Se actualiza la potencia generada en la matriz Nodos
        Generadores(:,2) = Pgen(Pgen_pos); 

        %Se efectua un estudio de flujos de potencia 
        [V,P_per,Qgen,P_te,Pslk_te] = load_flow_NR(Lineas,Nodos,Generadores,toleranciaNR,iter_max); 

        %Se actualizan las magnitudes de voltajes y angulos de los nodos en la matriz Nodos
        Nodos(:,4) = angle(V);
        Nodos(:,3) = abs(V)  ;
        %Se calcula la diferencia entre la nueva potencia y  la potencia de la iteracion pasada
        dif = max(Pgen-P_old); 

        if dif > tolerancia 
            conv_flag = 0;
          else
            conv_flag = 1;
        end
    end

    %Se calcula el costo total de generacion
    %Costo total
    CT(1,t) = sum(Costos_old(:,1)+Costos_old(:,2).*((Pgen(Pgen_pos))*Base_MVA)+Costos_old(:,3).*(((Pgen(Pgen_pos))*Base_MVA).^2));
    P_marg(:,t) = la;
end

%Proceso para mostrar resultados en pantalla  
%Mostrar datos generales del sistema 
disp('                             Despacho Eonomico Optimo')
disp(date)
fprintf("    Periodo estudiado    : %i horas\n",t)
fprintf("    Nodo Slack           : %i\n",Slk_pos)
[M1,I1] = max(CT);
fprintf("    Costo maximo         : %.2f [$] en la hora %i\n",M1,I1)
[M2,I2] = min(CT);
fprintf("    Costo minimo         : %.2f [$] en la hora %i\n",M2,I2)
fprintf("\n")
fprintf("    Desplegando graficas...\n")
disp('')

% %Mostrar datos de los nodos del sistema asÃƒÂ­ como los costos incrementales 
% res1=[Nodos(:,1),abs(V),angle(V)*(180/pi) la Pgen Qgen Nodos_old(:,5) Nodos(:,6)];
% disp('                (Pu y Grados)    [$/MWh]     GeneraciÃƒÂ³n(Pu)        Demanda(Pu)')
% fprintf("     Nodo     Voltaje   Angulo      %s        Real    Reactiva    Real   Reactiva \n",char(955));
% disp('    ----------------------------------------------------------------------------- ')
% disp(res1(:,1:8))
% 
% %Se calculan flujos de potencia
% buses    = Nodos(:,1)                                                                           ;
% De       = buses(round(Lineas(:,1)))                                                            ; 
% Hacia    = buses(round(Lineas(:,2)))                                                            ; 
% %Vector de admitancias de los nodos
% YU       = ones(length(Lineas(:,1)),1)./(complex(Lineas(:,3),Lineas(:,4)))                      ; 
% %Flujos de potencia del nodo x al y
% flujosxy = V(De).*conj((V(De)-V(Hacia)).*YU+V(De).*((complex(Lineas(:,5),Lineas(:,6)))/2))      ; 
% %Flujos de potencia del nodo y al x
% flujosyx = V(Hacia).*conj((V(Hacia)-V(De)).*YU+V(Hacia).*((complex(Lineas(:,5),Lineas(:,6)))/2)); 
% Pxy      = real(flujosxy);Qxy=imag(flujosxy)                                                    ;
% Pyx      = real(flujosyx);Qyx=imag(flujosyx)                                                    ;
% 
% %Mostrar datos de los flujos de potencia en las lineas del sistema
% res2=[(1:length(Lineas(:,1))).' Lineas(:,1) Lineas(:,2) Pxy Qxy Pyx Qyx];
% disp('                                      Sxy(Pu)              Syx(Pu)')
% disp('            Flujos de Linea (pu)')
% disp('    Linea    Del Nodo   Al Nodo    Real    Reactiva    Real    Reactiva ')
% disp('    -------------------------------------------------------------------- ')
% disp(res2(:,1:7))
%fprintf("    El costo total de producciÃƒÂ³n es: %f [$/h]\n\n",CT)

%GrÃƒÂ¡ficas
t = 1:1:24;

figure
subplot(1,2,1)
hold on
grid on
title('Costos totales');
plot(t,CT, 'LineWidth', 2)
xlabel('Hora')
ylabel('Costo Total')

subplot(1,2,2)
hold on
grid on
for k=1:Num_Nod
    plot(t,P_marg(k,:), 'LineWidth', 2)
    hold on
end
title('Precios marginales');
legend('N1','N2','N3','N4','N5','N6','N7','N8','N9','Location','northeast')
xlabel('Hora')
ylabel('Precio Marginal')
