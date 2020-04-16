function [Y,B] = Ybus_graf(Lineas,Num_Nod,Num_Lin,De,Hacia,Slk_pos)
%
% Sintaxis: 1.-  [Y,B] = Ybus_graf(Lineas,Num_Nod,Num_Lin,De,Hacia,Slk_Pos)
%           2.-  [Y]   = Ybus_graf(Lineas,Num_Nod,Num_Lin,De,Hacia,Slk_Pos)
%           3.-  [~,B] = Ybus_graf(Lineas,Num_Nod,Num_Lin,De,Hacia,Slk_Pos)
%
% Objetivo: Optener la matriz de admitancias nodales Y y las matriz B por el metodo de grafos.
% Entradas: 
%           Lineas    - Matriz que contiene los datos:
%                             -r,x,g,b de las lineas del sistema
%           Num_Nod   - Numéro de Nods
%           Num_Lin   - Numéro de Lineas
%           De        - Vector De de las Lineas
%           Hacia     - Vector Hacia de las Lineas
%           Slk_pos   - Posición del nodo Slack
%
% Salidas:  
%           Y      - Matriz de admitancias nodales Y
%           B      - Matriz B

%Se crean las matrices de incidencias y de matriz primitiva
Ybr_1  = zeros(Num_Lin+Num_Nod,Num_Lin+Num_Nod);
I_A    = zeros(size(Lineas,1)+Num_Nod,Num_Nod) ;

%Ciclo for para llenar matriz de incidencias
for m = 1:Num_Nod
    for n = 1:Num_Lin
        if Lineas(n,1) == m
            I_A(n,m) = 1;
        end
        if Lineas(n,2) == m
            I_A(n,m) = -1;
        end
        Ybr_1(n,n) = (1/(complex(Lineas(n,3),Lineas(n,4))));
    end
    if Lineas(m,6) ~= 0
        I_A(Num_Lin+m,m) = 1;
    end
end

%Ciclo para llenar la matriz primitiva
for k = 1:Num_Nod
    for m = 1:Num_Lin
        if De(m) == k || Hacia(m) == k
            Ybr_1(Num_Lin+k,Num_Lin+k) = Ybr_1(Num_Lin+k,Num_Lin+k)+((complex(Lineas(m,5),Lineas(m,6)))/2);
        end
    end
end

%Se obtiene la matriz Ybus
Y = (I_A')*(Ybr_1)*(I_A);

%Proceso para obtener la matriz B

%Para los componentes diagonales de la matriz B
B  = zeros(Num_Nod,Num_Nod);
for k = 1:Num_Nod
    for m = 1:Num_Lin
        if De(m) == k || Hacia(m) == k
            B(k,k) = B(k,k)+(1/Lineas(m,4));
        end
    end
end

%Para los componentes no diagonales de la matriz B
for n = 1:Num_Lin
    B(De(n),Hacia(n))  = -1*(1/Lineas(n,4));
    B(Hacia(n),De(n))  = B(De(n),Hacia(n)) ;
end
end
