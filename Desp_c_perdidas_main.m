%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Despacho Economico Optimo considerando perdidas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

%Se carga el caso de estudio
exa_5_buses;

%Datos de tolerancia e iteraciones máximas para los metodos de Flujos de potencia y despacho ecónomico.
tolerancia   = 1e-3;  %tolerancia para la convergencia
toleranciaNR = 1e-6;  %tolerancia para Newton Raphson 
iter_max     = 50  ;  %iteraciones maximas antes de detener el calculo para ambos ciclos

%Se hacen respaldos
Nodos_old  = Nodos;
Costos_old = Costos;

%Se crean vectores con posiciones e indices necesarios
Slk_pos  = find(Nodos(1:size(Nodos,1),2) == 1); %Posición del nodo Slack
Pgen_pos = find(round(Nodos(:,2)) < 3)        ; %Posición de los nodos con generadores
PVQ_pos  = find(round(Nodos(:,2)) >=2)        ; %Posicion de los nodos PV y PQ
Num_Nod  = size(Nodos,1)                      ; %Numéro de nodos

%Despacho economico simplificado considerando a la red para optener potencia generadas iniciales
[Pgen] = DespE_LT(Lineas,Nodos,Base_MVA,Costos);

%Se actualiza la potencia generada en la matriz Nodos
Nodos(:,5) = Pgen; 

%Se efectua un estudio de flujos de potencia 
[V,P_per,Qgen,P_te,Pslk_te] = load_flow_NR(Lineas,Nodos,toleranciaNR,iter_max);

%Se actualizan las magnitudes de voltajes y ángulos de los nodos en la matriz Nodos
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
    Nodos(:,7)   = Nodos_old(:,7) + PL;
    
    %Se obtienen los nuevos valores de "b" y "d" para obtner nuevas curvas de costos
    Costos(:,2)=Costos(:,2)./(1-ITL(Pgen_pos)).';
    Costos(:,3)=Costos(:,3)./(1-ITL(Pgen_pos)).';
    
    %Despacho simplificado considerando a la red y las nuevas cuvas de costos
    [Pgen,la] = DespE_LT(Lineas,Nodos,Base_MVA,Costos);
    
    %Se actualiza la potencia generada en la matriz Nodos
    Nodos(:,5) = Pgen;
    
    %Se efectua un estudio de flujos de potencia 
    [V,P_per,Qgen,P_te,Pslk_te] = load_flow_NR(Lineas,Nodos,toleranciaNR,iter_max); 
    
    %Se actualizan las magnitudes de voltajes y ángulos de los nodos en la matriz Nodos
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculos Finales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Costo total
CT=sum(Costos_old(:,1)+Costos_old(:,2).*((Pgen(Pgen_pos))*Base_MVA)+Costos_old(:,3).*(((Pgen(Pgen_pos))*Base_MVA).^2));
%Datos de nodos del sistema
res1=[Nodos(:, 1:1),abs(V),angle(V)*(180/pi) la Pgen Qgen Nodos_old(:,7:7) Nodos(:,8:8)];
%Datos de flujos del sistema
buses  = Nodos(:,1);
De     = buses(round(Lineas(:,1))); 
Hacia  = buses(round(Lineas(:,2))); 
YU = ones(length(Lineas(:,1)),1)./(complex(Lineas(:,3),Lineas(:,4)));
flujosxy= V(De).*conj((V(De)-V(Hacia)).*YU+V(De).*((complex(Lineas(:,5),Lineas(:,6)))/2));
flujosyx= V(Hacia).*conj((V(Hacia)-V(De)).*YU+V(Hacia).*((complex(Lineas(:,5),Lineas(:,6)))/2));
Pxy=real(flujosxy);Qxy=imag(flujosxy);
Pyx=real(flujosyx);Qyx=imag(flujosyx);
res2=[(1:length(Lineas(:,1))).' Lineas(:,1:1) Lineas(:,2:2) Pxy Qxy Pyx Qyx];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Despliegue de Resultados
disp('                             Despacho Eonómico Óptimo')
disp(date)
fprintf("    Número de iteraciones: %f\n",iter)
fprintf("    Nodo Slack           : %f\n",Slk_pos)
disp('')
disp('                (Pu y Grados)    [$/MWh]     Generación(Pu)        Demanda(Pu)')
fprintf("     Nodo     Voltaje   Angulo      %s        Real    Reactiva    Real   Reactiva \n",char(955));
%disp('     Nodo     Voltaje   Angulo  Lambda   char(955)          Real    Reactiva     Real   Reactiva ')
disp('    ----------------------------------------------------------------------------- ')
disp(res1(:,1:8))
disp('            Flujos de Linea (pu)')
disp('                                      Sxy(Pu)              Syx(Pu)')
disp('    Linea    Del Nodo   Al Nodo    Real    Reactiva    Real    Reactiva ')
disp('    -------------------------------------------------------------------- ')
disp(res2(:,1:7))
fprintf("    El costo total de producción es: %f [$/h]\n\n",CT)