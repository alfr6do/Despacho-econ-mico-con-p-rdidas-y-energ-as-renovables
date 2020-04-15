function [Pgen_n,la,ang_n,b,d] = Desp_line(Dat,V_nod,Base_MVA,Pdem,b,d,ITL)
%
% Sintaxis: 1.- [Pgen_n,~,ang_n]     = Desp_line(Dat,V_nod,Base_MVA,Pdem,b,d,0)
%          2.-  [Pgen_n,~,ang_n,b,d] = Desp_line(Dat,V_nod,Base_MVA,Pdem,b,d,ITL)
%
% Objetivo: Resolver un despacho ecónomico considerando a la red.
    Entrada:
                Dat - Matriz que contiene los datos de lassiguientes parámetros 
                                        MVA Base    - Base del sistema
                                        Nodos       - Datos de los nodos del sistema
                                        Líneas      - Datos de las líneas de transmisión y transformadores
                                        Gen         - Información necesaria para una simulación dinámica
                                        Simu........- Información necesaria para una simulación dinámica
                Tolerancia  - Tolerancia para la convergencia (1e-6 por defecto)
                Iter_max    - Iteraciones máximas (100 por defecto)
                Método      - Método iterativo que se emplea 
                                        "NR" = 'Newton Raphson' (Por defecto)
                                        "DR" = 'Desacoplado Rapido'
    Salida:
                Resultados  - Diccionario que incluye los siguientes parámetros 
                                        V           - Voltaje en todos los nodos
                                        Y           - Ybus del sistema
                                        Pgen        - Datos de los nodos del sistema
                                        Qgen        - Datos de las líneas de transmisión y transformadores
                                        Pper        - Información necesaria para una simulación dinámica
                                        Qper        - Información necesaria para una simulación dinámica
                                        Pxy         - Flujos de potencia activa del nodos x al y
                                        Qxy         - Flujos de potencia reactiva del nodos x al y
                                        Pyx         - Flujos de potencia activa del nodos y al x
                                        Qyx         - Flujos de potencia reactiva del nodos y al a
                                        Iter        - Iteraciones necesarias para que el metodo convergiera
                                        Metodo      - Método de Solución empleado
                                        tconv       - Tiempo en que convergio el método de solución
                                        ttot        - Tiempo total de solución
                                        conv_flag   - Bandera que indica si el método convergio
    """



%ITL introducir unicamente si se esta actualizando la curva de costos dentro de un ciclo iterativo de despacho con perdidas
%P_per introducir unicamente si se esta actualizando la curva de costos dentro de un ciclo iterativo de despacho con perdidas
% DE LO CONTRARIO, COLOCAR ITL Y P_PER en "0"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Despacho Economico Simplificado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_Nod = size(V_nod,1);
Num_Lin = size(Dat,1); 
De     = Dat(:,1); 
Hacia  = Dat(:,2);
Slk_pos = find(V_nod(1:Num_Nod, 2:2) == 1); %Slack
tipo_bus = round(V_nod(:,2)); %columa que indica tipo de buses en orden
Pgen_pos=find(tipo_bus < 3);
Num_Gen=length(find(tipo_bus==2))+1; %Numero de generadores
[~,B] = Ybus_graf(Num_Nod,Num_Lin,Dat,De,Hacia,Slk_pos);
%Dependiendo de si en la funcion se le mete un valor ITL o no, se actualizara la d y la b o seran las inciales
if sum(ITL)~=0
    for i=1:Num_Gen
        b(i)=b(i)/(1-ITL(Pgen_pos(i)));
        d(i)=d(i)/(1-ITL(Pgen_pos(i)));
    end
end
A1=diag(2*d);
A2=-1*eye(Num_Gen,Num_Nod);
A3=zeros(Num_Gen,Num_Nod-1);
B1=A2.';
B2=zeros(Num_Nod,Num_Nod);
B3=B;
B3(:,Slk_pos)=[];
B3=B3*Base_MVA;
C1=A3.';
C2=B3.';
C3=zeros(Num_Nod-1,Num_Nod-1);
lmat_c1=horzcat(A1,A2,A3);
lmat_c2=horzcat(B1,B2,B3);
lmat_c3=horzcat(C1,C2,C3);
lagmat=vertcat(lmat_c1,lmat_c2,lmat_c3);
sol=lagmat\[-1*b;-1*Pdem*Base_MVA;zeros(Num_Nod-1,1)];
temp=(sol(1:Num_Gen)./Base_MVA);
Pgen_n=zeros(Num_Nod,1);
Pgen_n(Pgen_pos)=temp;
la=sol(length(Pgen_pos)+1:length(Pgen_pos)+Num_Nod);
temp1=sol(length(Pgen_pos)+Num_Nod+1:end);
ang_n=[0 ; temp1];
end
