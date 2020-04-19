%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Despacho Economico Completo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objetivo: Realizar un despacho economico optimo en un sistema considerando la red de transmisiÃ³n y 
%           sus perdidas. Ademas, las restricciones que se consideran en este programa son:
%               *RestricciÃ³n de generaciÃ³n.
%               *Restriccion de capacidad de lineas de tranmisiÃ³n.
%               *Restriccion de voltaje en los nodos.

clear all
clc

%Se carga el caso de estudio
caso9n3g;

%Despacho simplificado considerando a la red y las nuevas cuvas de costos
[sol] = DE_DC(Lineas,Nodos,Generadores,Base_MVA,Costos);
sol