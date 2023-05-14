clc
clear
clf

n = 8; // panels
phi(1)= %pi / 2; // 90°
for i = 1 : n
    alpha = (2 * %pi) / n; // internal angle
    X(i) = cos(%pi + (alpha / 2) - (i - 1) * alpha);
    Y(i) = sin(%pi + (alpha / 2) - (i - 1) * alpha);
end

X(n + 1) = X(1);
Y(n + 1) = Y(1);

for i = 1 : n;
    x(i) = (X(i+1) + X(i)) / 2;
    y(i)=(Y(i+1) + Y(i)) / 2;
    phi(i)= phi(1) - (i - 1) * alpha; // phi - angle each panel
    bet(i)= phi(i) + (%pi / 2); // beta - normal angle from each panel
    U(i)= -2 * %pi * cos(bet(i)); // vector U
end

//plot(X,Y,'bo') 
//plot(x,y,'ro')

S = sqrt(((X(2) - X(1))^2) + ((Y(2) - Y(1))^2)); // panel length

for i = 1 : n
    for j = 1 : n
     if i == j then
            Iijn(i,j) = %pi;
        else
            A(i,j) = -((x(i) - X(j)) * cos(phi(j))) - ((y(i) - Y(j)) * sin(phi(j)));
            B(i,j) = (x(i) - X(j))^2 + (y(i) - Y(j))^2;
            C(i,j) = sin(phi(i) - phi(j));
            D(i,j) =(y(i) - Y(j)) * cos(phi(i))-(x(i)-X(j))*sin(phi(i));
            E(i,j) = sqrt(B(i,j)-A(i,j)^2);
            // Matrix Iijn divided by three terms: I1, I2, I3
            I1(i,j) = (C(i,j) / 2) * log((S^2 + (2 * A(i,j) * S) + B(i,j)) / B(i,j));
            I2(i,j) = atan((S + A(i,j)) / E(i,j));
            I3(i,j) = atan(A(i,j) / E(i,j));
            Iijn(i,j)= I1(i,j) + ((D(i,j) - (A(i,j) * C(i,j))) / E(i,j)) * (I2(i,j) - I3(i,j));
        end
    end
end

lambida = inv(Iijn) * U; // Iijn^-1 * vector U

for i = 1 : n
    Vs(i) = 0;
    for j = 1 : n
        if i == j then
            IS(i,j) = 0
        else
            IS1(i,j) = (D(i,j) - A(i,j) * C(i,j)) / (2 * E(i,j));
            IS2(i,j) = (S^2 + 2 * A(i,j) * S + B(i,j)) / B(i,j);
            IS3(i,j) = atan((S + A(i,j)) / (E(i,j))) - atan(A(i,j) / E(i,j));
            IS(i,j) = IS1(i,j) * log(IS2(i,j)) - (C(i,j) * IS3(i,j));
        end
        Vs(i) = Vs(i) + lambida(j )*IS(i,j) / (2 * %pi);
    end

    Vs(i) = Vs(i) + sin(bet(i));
    cp(i) = 1 - (Vs(i)^2);
end

teta = linspace(-%pi, %pi, 100);
cpt = 1 - 4 * (sin(teta))^2; // theoretical values of Cp

// graph plot theoretical values of Cp versus panel methods Cp
plot(teta, cpt, ':b', "thickness",3);
xlabel("teta", "fontsize",4);
plot(bet, cp, 'o', "thickness", 3), title("Panel method Cp x theoretical Cp",
"fontsize",5);
ylabel("Cp", "fontsize",4);
legend(['Cp theoretical';'Cp each panel']);