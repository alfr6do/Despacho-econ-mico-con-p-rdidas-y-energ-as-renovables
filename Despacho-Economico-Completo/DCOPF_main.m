%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Despacho Economico Completo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objetivo: Realizar un despacho economico optimo en un sistema considerando la red de transmisiÃƒÂ³n y 
%           sus perdidas. Ademas, las restricciones que se consideran en este programa son:
%               *RestricciÃƒÂ³n de generaciÃƒÂ³n.
%               *Restriccion de capacidad de lineas de tranmisiÃƒÂ³n.
%               *Restriccion de voltaje en los nodos.

clear all
clc

%Se carga el caso de estudio
caso9n3g;

%Estudio de flujos de potencia optimo en DC
[Pgen,la,ang,Pxy,Gen_Lm,Lin_Lm,mu] = DE_DC(Lineas,Nodos,Generadores,Base_MVA,Costos);

%Se calculan los costos totales de producción
Pgen_pos = find(round(Nodos(:,2)) < 3)        ; %Posicion de los nodos con generadores
CT       = sum(Costos(:,1)+Costos(:,2).*((Pgen(Pgen_pos)))+Costos(:,3).*(((Pgen(Pgen_pos))).^2));

%Despliegue de resultados. 
%Mostrar datos generales del sistema 
disp('                   Flujos de potencia optimos en DC')
disp(date)
Num_Nod  = size(Nodos,1)                        ; %Numero de Nodos
Slk_pos  = find(Nodos(1:Num_Nod,2) == 1)     ; %Posicion del nodo Slack
fprintf("Nodo Slack                : %f\n",Slk_pos)
fprintf("Costo total de produccion : %f [$/h]\n",CT)
fprintf("\n")
%Mostrar datos de los nodos del sistema asi como los costos incrementales 
res1=[Nodos(:,1), Nodos(:,3),ang*(180/pi), la ,Pgen, Nodos(:,5)];
disp('    Informacion de nodos')
disp('                (Pu y Grados)    [$/MWh]  Pgen(MW)   Pdem(MW)')
fprintf("     Nodo     Voltaje   Angulo      %s      Activa    Activa  \n",char(955));
disp('    ----------------------------------------------------------')
disp(res1(:,1:6))

%Mostrar datos de los flujos de potencia en las lineas del sistema
res2=[(1:length(Lineas(:,1))).' Lineas(:,1) Lineas(:,2) Pxy -1 * Pxy ];
fprintf("\n")
disp('    Flujos de Linea (MW)')
disp('    Linea    Del Nodo   Al Nodo    Pxy      Pxy')
disp('    ------------------------------------------------')
disp(res2(:,1:5))

%Se muestra si es que se rompio algun limite

%Limite en generadores
if Gen_Lm ~= 0
    mugen = mu(1:length(Gen_Lm));
    res3  = [Gen_Lm.', Generadores(Gen_Lm).', Generadores(Gen_Lm,4)...
            ,Generadores(Gen_Lm,5), Pgen(Gen_Lm), abs(mugen)];
    fprintf("\n")
    disp('    Restricciones de potencia generada')
    disp('      Gen      Nodo     Pmax       Pmin      Pgen      mu')
    disp('    ---------------------------------------------------------')
    disp(res3(:,1:6))    
end
if Lin_Lm ~= 0 
    mulin = mu(length(Gen_Lm)+1:end);
    res3  = [Lin_Lm.', Lineas(Lin_Lm,1), Lineas(Lin_Lm,2), Lineas(Lin_Lm,7)...
            ,-1 * Lineas(Lin_Lm,7), Pxy(Lin_Lm), abs(mulin)];
    fprintf("\n")
    disp('    Restricciones de flujo de potencia en las linesas')
    disp('    Linea    Del Nodo   Al Nodo    Pmax     Pmin      Pxy         mu')
    disp('    ------------------------------------------------------------------')
    disp(res3(:,1:7))  
end
