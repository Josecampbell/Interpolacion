function[]=Interpolacion_Hermite()
clc
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Entrada de Datos--------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=[1.3 1.6 1.9];%nodos de el eje x
yn=[0.6200860 0.4554022 0.2818186];%nodos de el eje y
dyn=[-0.5220232 -0.5698959 -0.5811571];%imagenes de las derivadas
x0=1.5;%punto a evaluar en la interpolacion de Hermite
Q=Hermite(xn,yn,dyn);%llamada al metodo de de Hermite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Parametros de salidad-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q es una matrix digonal inferior cuyas entradas son las iteraciones del
% metodo Hermite
P=Polinomio(Q,xn);%Polinomio de Hermite
P=matlabFunction(P);
y0=P(x0);%punto evaludo en el polinomio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Grafico de los Resultados-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
alfa=xn(length(xn))-xn(length(xn)-1);%extremo izquierdo del intervalo donde se quiere visualizar el estudio
alfa1=xn(length(xn))+xn(length(xn)-1);%extremo derecho del intervalo donde se quiere visualizar el estudio
hold on;fplot(P,[alfa,alfa1]);scatter(xn,yn,'MarkerEdgeColor','g');scatter(x0,y0,'MarkerEdgeColor','r');legend('Polinomio de Hermite','nodos interpolantes','Punto evaludo en polinomio');hold off;
xlabel('eje x') %coamando xlabel para colocar nombre a el eje x
ylabel('eje y') %coamando label para colocar nombre a el eje y
title('Polinomio de Hermite')%el comando title para colocar nombre a las graficas
save('Hermite.mat','P') %guardando los resultados de Neville


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Metodo --------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q]=Hermite(xn,yn,dyn)
n=length(xn);%numero de iteraciones del metodo de Interpolacion Hermite

for i=0:n-1
     z(2*i+1)=xn(i+1);
    z(2*i+2)=xn(i+1);
    Q(2*i+1,1)=yn(i+1);
    Q(2*i+2,1)=yn(i+1);
    Q(2*i+1,2)=dyn(i+1);
    if(i~=0)
        Q(2*i,2)=(Q(2*i+1,1)-Q(2*i,1))/( z(2*i+1)-z(2*i));
    end
end
aux(1)=0;
aux(2:length(Q(:,2)))=Q(1:length(Q(:,2))-1,2);
Q(:,2)=aux;    
for i=3:2*n
    for j=3:i
        Q(i,j)=(Q(i,j-1)-Q(i-1,j-1))/(z(i)-z(i-j+1));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Polinomio Interpolante ------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[P]=Polinomio(Q,xn)%funcion para construir el polinomio de Hermite
syms x %comando para definir parametros
n=length(xn);%Grado del polinomio
P=Q(1,1)+Q(2,2)*(x-xn(1));%primeros terminos del polinomio
for i=2:n
    if(mod(i,2)==0)
        variable=kron(x,ones(1,i));% vector [x,x...,x] de dimension i
         variable=variable-xn(1:i);% vector [x-xn(1),x-xn(2),...,x-xn(i)]
         variable=variable.^2;% vector [(x-xn(1))^2,(x-xn(2))^2,...,(x-xn(i-1))^2,(x-xn(i))^2]
         variable=prod(variable);% producto de los terminos del vector variable (x-xn(1))^2*(x-xn(2))^2*...(x-xn(i))^2
         P=P+Q(i,i)*variable;%polinomio de Hermite
    else
        variable=kron(x,ones(1,i));% vector [x,x...,x] de dimension i
         variable=variable-xn(1:i);% vector [x-xn(1),x-xn(2),...,x-xn(i)]
         variable(1:length(variable)-1)=variable(1:length(variable)-1).^2;% vector [(x-xn(1))^2,(x-xn(2))^2,...,(x-xn(i-1))^2,(x-xn(i))]
         variable=prod(variable);% producto de los terminos del vector variable (x-xn(1))^2*(x-xn(2))^2*...(x-xn(i-1))^2*(x-xn(i))
         P=P+Q(i,i)*variable;%polinomio de Hermite
    end
end

end