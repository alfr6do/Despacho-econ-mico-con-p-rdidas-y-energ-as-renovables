function [V,P_per,Qgen,P_te,Pslk_te,Pgen] = load_flow_NR(Lineas,Nodos,Generadores,tolerancia,iter_max)
%
% Sintaxis: 1.- [V,P_per,Qgen,P_te,Pslk_te,Pgen] = load_flow_NR(Lineas,Nodos,tolerancia,iter_max
%           2.- [V,P_per,Qgen,~,~,Pgen] = load_flow_NR(Lineas,Nodos,tolerancia,iter_max
%           3.- [V,P_per,Qgen] = load_flow_NR(Lineas,Nodos,tolerancia,iter_max
%
% Objetivo: Resolver un estudio de flujos de potencia por medio del metodo de Newton Raphson.
% Entradas: 
%           Lineas     - Matriz que contiene los datos:
%                              -r,x,g,b de las lineas del sistema
%           Nodos      - Matriz que contiene los datos:
%                              -Tipo,V,Ang,Pgen,Qgen,Pdem,Qdem,Vmax,Vmin de los nodos del sistema
%           tolerancia - tolerancia para la convergencia del metodo
%           iter_max   - Iteraciones mÃ¡ximas para las convergencia del metodo
% Salidas:  
%           V      - Voltajes en todos los nodos del sistema
%           P_per  - Perdidas totales de potencia activa en el sistema
%           Q_gen  - Potencia reactiva generada
%           P_gen  - Potencia activa generada
%
%           ***Si se esta realizando dentro de un despacho economico***
%           P_slk_te  - Vector del jacobiano de de potencia reactiva respecto al angulo del nodo Slack
%           P_te  - Vector del jacobiano de de potencia reactiva respecto al angulo del sistema


Num_Lin = size(Lineas,1)                  ; %numero de lineas
Num_Nod = size(Nodos,1)                   ; %numero de buses
De      = Lineas(:,1)                     ; 
Hacia   = Lineas(:,2)                     ;
Slk_pos = find(Nodos(1:Num_Nod, 2:2) == 1); %Posicion del nodo Slack
Pgen_pos = find(round(Nodos(:,2)) < 3)        ; %Posicion de los nodos con generadores

%Se calcula la matriz de admitancias nodales
[Y] = Ybus_graf(Lineas,Num_Nod,Num_Lin,De,Hacia,Slk_pos);

%Se extraen los datos necesario de la matriz Nodos
V        = (Nodos(:,3)) ; %magnitud de voltajes
ang      = (Nodos(:,4)) ; %angulos en radianes
Pgen     = zeros(Num_Nod,1);
Qgen     = zeros(Num_Nod,1);
Pgen(Pgen_pos) = (Generadores(:,2)) ; %Potencias Activas Generadas
Qgen(Pgen_pos) = (Generadores(:,3)) ; %Potencias Reactivas Generadas
Pdem     = (Nodos(:,5)) ; %Potencias Activas Demandadas
Qdem     = (Nodos(:,6)) ; %Potencias Reactivas Demandas

%Se crean vectores con posiciones e indices necesarios
tipo_bus                   = round(Nodos(:,2))   ; %columa que indica tipo de buses en orden
PVQ_pos                    = find(tipo_bus >=2)  ; %posicion de los nodos PV Y PQ
PQ_pos                     = find(tipo_bus==3)   ; %Posiciion de los nodos PQ
PV_pos                     = find(tipo_bus==2)   ; %Posicion de los nodos PV
PQ_num                     = length(PQ_pos)      ; %Numero de nodos PQ
PVQ_num                    = length(PVQ_pos)     ; %Numero de nodos PV y PQ
sw_bno                     = ones(Num_Nod,1)     ; %vector de unos de dimension el numero de nodos
g_bno                      = sw_bno              ; %vector de unos de dimension el numero de nodos
bus_zeros                  = zeros(Num_Nod,1)    ; %vector de ceros de dimension el numero de nodos
sw_bno(Slk_pos)            = bus_zeros(Slk_pos)  ; %El vector sw_bno pero con un cero en el nodo Slack
g_bno(PV_pos)              = bus_zeros(PV_pos)   ; %el vector g_bno pero con ceros en los nodos PV
ang_re                     = eye(Num_Nod,Num_Nod); %Matriz diagonal de unos de dimension el numero de nodos
ang_re(Slk_pos,:)          = []                  ; %La matriz ang_re sin la fila del nodo slack
vol_re                     = eye(Num_Nod,Num_Nod); %Matriz diagonal de unos de dimension el numero de nodos
vol_re([Slk_pos;PV_pos],:) = []                  ; %La matriz vol_re sin las fila del nodo slak y nodos PV


%Se calculan los desbalances de potencia incial
S          = (V.*exp(1i*ang)).*conj(Y*(V.*exp(1i*ang))); %Se calculan potencias calculadas
Delta_P    = Pgen-Pdem-real(S)                         ; %Desbalance de potencia activa
Delta_Q    = Qgen-Qdem-imag(S)                         ; %Desbalance de potencia activa
Delta_P    = Delta_P.*sw_bno                           ; %Cero en la posicion del nodo Slack en el desbalance de P
Delta_Q    = Delta_Q.*sw_bno                           ; %Cero en la posicion del nodo Slack en el desbalance de Q
Delta_Q    = Delta_Q.*g_bno                            ; %Cero en la posicion de los nodos PV en el desbalance de Q
Pdes       = max(abs(Delta_P))                         ; %Desbalance maximo de P
Qdes       = max(abs(Delta_Q))                         ; %Desbalance maximo de Q
desbalance = Pdes+Qdes                                 ; %Se obtiene el desbalance total

%Con lo calculado previamente se verifica convergencia
if desbalance > tolerancia 
    conv_flag = 0;
  else
    conv_flag = 1;
end

%Comienza el proceso iterativo
iter=0; %iter incia
while conv_flag==0 && iter <= iter_max
    iter=iter+1;
    
    %CreaciÃ³n del Jacobiano
    P_te= zeros(PVQ_num,PVQ_num) ; %Jacobiano -- Potencia Activa respecto al angulo
    P_ve= zeros(PVQ_num,PQ_num)  ; %Jacobiano -- Potencia Activa respecto al voltaje
    Q_te= zeros(PQ_num,PVQ_num)  ; %Jacobiano -- Potencia Reactiva respecto al angulo
    Q_ve= zeros(PQ_num,PQ_num)   ; %Jacobiano -- Potencia Reactiva respecto al voltaje
    
    %creaciÃ³n P cal respecto a teta
    for s= 1:PVQ_num
        for f = 1:PVQ_num
            if s == f
                P_te(s,f) = -1*(imag(S(PVQ_pos(f))))-imag(Y(PVQ_pos(f),PVQ_pos(f)))*((abs(V(PVQ_pos(f))))^2);            
            else
                P_te(s,f) = (abs(V(PVQ_pos(s))))*abs(V(PVQ_pos(f)))*(((real(Y(PVQ_pos(s),PVQ_pos(f)))*sin((ang(s+1))-(ang(PVQ_pos(f))))))...
                          -(imag(Y(PVQ_pos(s),PVQ_pos(f)))*cos((ang(PVQ_pos(s)))-(ang(PVQ_pos(f))))));
            end

        end
    end
    
    %creaciÃ³n P cal respecto a ve
    for s= 1:PVQ_num
        for f = 1:PQ_num
            if (PVQ_pos(s)) == (PQ_pos(f))
                P_ve(s,f) = ((real(S(PQ_pos(f))))+real(Y(PQ_pos(f),PQ_pos(f)))*((abs(V(PQ_pos(f))))^2))/(abs(V(PQ_pos(f))));            
            else
                P_ve(s,f) = ((abs(V(PVQ_pos(s)))*abs(V(PQ_pos(f))))*((real(Y(PVQ_pos(s),PQ_pos(f)))*cos((ang(PVQ_pos(s)))-(ang(PQ_pos(f))))...
                            +(imag(Y(PVQ_pos(s),PQ_pos(f)))*sin((ang(PVQ_pos(s)))-(ang(PQ_pos(f)))))))/abs(V(PQ_pos(f))));
            end

        end
    end
    
    %creaciÃ³n Q cal respecto a teta
    for s= 1:PQ_num
        for f = 1:PVQ_num
            if (PQ_pos(s)) == (PVQ_pos(f))
                Q_te(s,f) = (real(S(PVQ_pos(f))))-real(Y(PVQ_pos(f),PVQ_pos(f)))*((abs(V(PVQ_pos(f))))^2);            
            else
                Q_te(s,f) = (-1*(abs(V(PQ_pos(s)))*abs(V(PVQ_pos(f)))))*((real(Y(PQ_pos(s),PVQ_pos(f)))*cos((ang(PQ_pos(s)))-(ang(PVQ_pos(f)))))...
                            +(imag(Y(PQ_pos(s),PVQ_pos(f)))*sin((ang(PQ_pos(s)))-(ang(PVQ_pos(f))))));
            end

        end
    end
    
    %creaciÃ³n Q cal respecto a ve
    for s= 1:PQ_num
        for f = 1:PQ_num
            if (PQ_pos(s)) == (PQ_pos(f))
                Q_ve(s,f) = ((imag(S(PQ_pos(s))))-imag(Y(PQ_pos(s),PQ_pos(s)))*((abs(V(PQ_pos(s))))^2))/(abs(V(PQ_pos(s))));            
            else
                Q_ve(s,f) = ((abs(V(PQ_pos(s)))*abs(V(PQ_pos(f))))*((real(Y(PQ_pos(s),PQ_pos(f)))*sin((ang(PQ_pos(s)))-(ang(PQ_pos(f)))))...
                            -(imag(Y(PQ_pos(s),PQ_pos(f)))*cos((ang(PQ_pos(s)))-(ang(PQ_pos(f)))))))/(abs(V(PQ_pos(f))));
            end
        end
    end
    
    %Se concatenan todas las partes para formar al Jacobiano
    Jac=[P_te P_ve; Q_te Q_ve];
    
    %Visualizando a los aumentos como Ax=b donde x son los aumentos.
    %Se obtiene a el vector b
    DeltaP_re = ang_re*Delta_P        ; %Se elimna la fila del vector el Slack del desbalance P
    DeltaQ_re = vol_re*Delta_Q        ; %Se eliminan las filas del nodo Slack y nodos PV del desbalance Q
    b         = [DeltaP_re; DeltaQ_re]; %Se concatenan los desbalances
    
    %Se soluciona para x tomando en cuenta que A = Jacobiano
    temp = Jac\b;
    
    %Se obtienen los aumentos de magnitud de voltaje y Ã¡ngulo 
    Ang_au = ang_re'*temp(1:length(PVQ_pos),:)                               ;
    V_au   = vol_re'*temp(length(PVQ_pos)+1:length(PVQ_pos)+length(PQ_pos),:);
    
    %Se actualizan los valores de magnitud de voltaje y Ã¡ngulo
    V   = abs(V) + V_au;
    ang = ang + Ang_au ;
   
    %Se calculan los desbalances de potencia nuevamente
    V          = V.*exp(1i*ang)   ;
    S          = V.*conj(Y*V)     ; %Se calculan potencias calculadas
    Delta_P    = Pgen-Pdem-real(S); %Desbalance de potencia activa 
    Delta_Q    = Qgen-Qdem-imag(S); %Desbalance de potencia activa
    Delta_P    = Delta_P.*sw_bno  ; %Cero en la posicion del nodo Slack en el desbalance de P
    Delta_Q    = Delta_Q.*sw_bno  ; %Cero en la posicipn del nodo Slack en el desbalance de Q
    Delta_Q    = Delta_Q.*g_bno   ; %Cero en la posicion de los nodos PV en el desbalance de Q
    Pdes       = max(abs(Delta_P)); %Desbalance maximo de P
    Qdes       = max(abs(Delta_Q)); %Desbalance maximo de Q
    desbalance = Pdes+Qdes        ; %Se obtiene el desbalance total
    
    %Con lo calculado previamente se verifica convergencia
    if desbalance > tolerancia 
        conv_flag = 0;
      else
        conv_flag = 1;
    end
end

%Se calculan Potencias generadas y demandadas
S              = V.*conj(Y*V)                    ; %Potencia Calculada
Pgen_pos       = find(tipo_bus < 3)              ;
Pgen(Pgen_pos) = real(S(Pgen_pos))+Pdem(Pgen_pos); %Potencia Activa Generada 
Qgen(Pgen_pos) = imag(S(Pgen_pos))+Qdem(Pgen_pos); %Potencia Reactiva Generada

%Se obtienen las Perdidas
P_per=sum(Pgen)-sum(Pdem);

%Elemento necesario cuando se utlizan flujos de potencia para un despacho economico con perdidas
Pslk_te= zeros(1,PVQ_num);
for f = 1:PVQ_num
    Pslk_te(1,f) = (abs(V(Slk_pos)))*abs(V(PVQ_pos(f)))*(((real(Y(Slk_pos,PVQ_pos(f)))*sin((ang(Slk_pos))-(ang(PVQ_pos(f))))))...
                    -(imag(Y(Slk_pos,PVQ_pos(f)))*cos((ang(Slk_pos))-(ang(PVQ_pos(f))))));
end
end
