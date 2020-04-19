function [sol] = DE_DC(Lineas,Nodos,Generadores,Base_MVA,Costos)
%
% Sintaxis: 1.-  [Pgen_n,la,ang] = DespE_LT(Lineas,Nodos,Base_MVA,Costos)
%           2.-  [Pgen_n,la]     = DespE_LT(Lineas,Nodos,Base_MVA,Costos)
% Objetivo: Resolver un despacho economico considerando a la red.
% Entradas: 
%           Lineas    - Matriz que contiene los datos:
%                             -r,x,g,b de las lineas del sistema
%           Nodos     - Matriz que contiene los datos:
%                             -Tipo,V,Ang,Pgen,Qgen,Pdem,Qdem,Vmax,Vmin de los nodos del sistema
%           Base_MVA  - Base del sistema
%           Costos    - Informaci√≥n de la curva de costos de los generadores
% Salidas:  
%           Pgen      - Vector con las potencias generadas luego del despacho economico
%           la        - Costos incrementales luego del despacho economico
%           ang       - Angulo calculado para los nodos de este sistema luego de este despacho

%Se crean vectores con posiciones e indices necesarios
Num_Nod  = size(Nodos,1)                        ; %Numero de Nodos
Num_Lin  = size(Lineas,1)                       ; %Numero de Lineas
Num_Gen  = size(Generadores,1)                  ; %Numero de generadores
De       = Lineas(:,1)                          ; %Vector De de las Lineas
Hacia    = Lineas(:,2)                          ; %Vector Hacia de las Lineas
Slk_pos  = find(Nodos(1:Num_Nod, 2:2) == 1)     ; %Posicion del nodo Slack
Pgen_pos = find(round(Nodos(:,2))  < 3)         ; %Posicion de generadores
PVQ_pos  = find(round(Nodos(:,2)) >=2)          ; %Posicion de los nodos PV y PQ
Num_Gen  = length(find(round(Nodos(:,2))==2))+1 ; %Numero de generadores
Lim_Pg   = zeros(Num_Gen,1)                     ; %Limite de generadores

%Se crean variables necesarias en para verificar limites
PgLim = 0;
exp   = 0;
x     = 0;

%Se calcula la matriz B
[~,B] = Ybus_graf(Lineas,Num_Nod,Num_Lin,De,Hacia,Slk_pos);

b = (Costos(:,2)); %Se obtiene el parametro "b" de las curvas de costo
d = (Costos(:,3)); %Se obtiene el parametro "d" de las curvas de costos

%Se crea la matriz A de la operacion la operacion Ax=b para resolver este despacho

%Primer Fila llamada A
A1 = diag(2*d)                ; 
A2 = -1*eye(Num_Gen,Num_Nod)  ;
A3 = zeros(Num_Gen,Num_Nod-1) ;

%Segundia Fila llamada B
B1           = A2.'                   ;
B2           = zeros(Num_Nod,Num_Nod) ;
B(:,Slk_pos) = []                     ; %Se elimna de la matriz B la columna del nodo Slack
B3           = B*Base_MVA             ;

%Tercera Fila llamada C
C1 = A3.'                       ;
C2 = B3.'                       ;
C3 = zeros(Num_Nod-1,Num_Nod-1) ;

%Se concatenan todas las columnas para formar a las filas 1 2 y 3 de la matriz
lmat_c1 = horzcat(A1,A2,A3);
lmat_c2 = horzcat(B1,B2,B3);
lmat_c3 = horzcat(C1,C2,C3);

%Se concatenana las filas para fomar a la matriz
A = vertcat(lmat_c1,lmat_c2,lmat_c3);

%Se crea al vector b de la operacion Ax=b
Pdem = (Nodos(:,7))                              ; %Potencias Activas Demandadas
b    = [-1*b;-1*Pdem*Base_MVA;zeros(Num_Nod-1,1)];


while PgLim==0
    
    %Se agregan filas y columnas de los generadores que rompen sus lÌmites
    if exp==1
        x=0;
        for k=1:Num_Gen
            if Lim_Pg(k,1)==1
                Fila = zeros(1,Num_Gen+Num_Nod+Num_Nod-1);
                Fila(1,Pgen_pos(k,1)) = -1;
                Columna = zeros(Num_Gen+Num_Nod+Num_Nod-1,1);
                Columna(Pgen_pos(k,1),1) = -1;
                A((Num_Gen+Num_Nod+Num_Nod),1:(Num_Gen+Num_Nod+Num_Nod-1)) = Fila;
                A(1:(Num_Gen+Num_Nod+Num_Nod-1),(Num_Gen+Num_Nod+Num_Nod)) = Columna;
                b = [b;-Generadores(k,4)*Base_MVA];
                x=1;
            end
            if Lim_Pg(k,1)==-1
                Fila = zeros(1,Num_Gen+Num_Nod+Num_Nod-1);
                Fila(1,Pgen_pos(k,1)) = -1;
                Columna = zeros(Num_Gen+Num_Nod+Num_Nod-1,1);
                Columna(Pgen_pos(k,1),1) = -1;
                A((Num_Gen+Num_Nod+Num_Nod),1:(Num_Gen+Num_Nod+Num_Nod-1)) = Fila;
                A(1:(Num_Gen+Num_Nod+Num_Nod-1),(Num_Gen+Num_Nod+Num_Nod)) = Columna;
                b = [b;-Generadores(k,5)*Base_MVA];
                x=1;
            end
        end
%         mu = zeros(x,x);
%         A((Num_Gen+Num_Nod+Num_Nod):end,(Num_Gen+Num_Nod+Num_Nod):end)=mu;
    end

    %Se obtiene la solucion de la operacion Ax=b
    sol  =A\b;

    %Se extraen a las potencias generadas del vector Sol
    Pgen           = zeros(Num_Nod,1)        ;
    Pgen(Pgen_pos) = sol(1:Num_Gen)./Base_MVA;
    
    %Se extraen el valor de los costos incrementales (Lambda) del vector Sol
    la = sol(length(Pgen_pos)+1:length(Pgen_pos)+Num_Nod);

    %Se extraen los angulos de los nodos calculados en este despacho del vector Sol
    ang          = zeros(Num_Nod,1)                   ;
    ang(PVQ_pos) = sol(length(Pgen_pos)+Num_Nod+1:end-x);
    
    %Se verifican los lÌmites de generaciÛn de potencia
    Lim_Pg  = zeros(Num_Gen,1);
    for k=1:Num_Gen
        if Pgen(Pgen_pos(k,1),1)>Generadores(k,4)
            Lim_Pg(k,1) = 1;
            break
        end
        if Pgen(Pgen_pos(k,1),1)<Generadores(k,5)
            Lim_Pg(k,1) = -1;
            break
        end
    end
    if sum(abs(Lim_Pg))==0
        PgLim = 1;
    else
        exp = 1;
    end
    
end

end
