#include <iostream>
#include <cmath>
#include <fstream>


#define w 2.6614e-6
#define G 6.67e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define dTL 3.844e8
#define RT 6.637816e6
#define RL 1.7374e6
#define PI 3.14159265



double f1(double pr);
double f2(double pphi, double r);
double f3(double pphi, double r, double phi, double t, double delta, double mu);
double f4(double delta, double mu, double r, double phi, double t);



using namespace std;

int main()
{   
    double h, t, taux, tmax; //paso y variable temporal
    double y[4], k1[4], k2[4], k3[4], k4[4]; //4 ecuaciones diferenciales y las funciones ki para aplicar el algoritmo
    double aux1[4], aux2[4], aux3[4]; //vectores auxiliares para facilitar el cálculo de ki
    double pr, r, pphi, phi;
    double v_inicial;

    double x_L, y_L; //Posición de la luna
    double x_T, y_T; //Posición de la Tierra (origen)
    double pos_x, pos_y; //para transformar las coordenadas polares en cartesianas para representarlo

    double delta, mu, r_prima;

    double m, H, rL;
    
    int i, j, iter;
    

    ofstream fich_posiciones;
    ofstream Hamiltoniano;


    //Establezco la posición inicial de la luna

    x_L=1.0;
    y_L=0.0;
    x_T=0.0;
    y_T=0.0;


    //Condiciones iniciales (reescalando)

    t=0.0;
    taux=0.0;
    h=60.0*2.0;
    tmax=60000.0*6.0;
    iter=tmax/h;

    
    r=0.017268;  // Radio de la Tierra reescalado en términos de la distancia Tierra-Luna (RT/dTL)
    phi=PI/8.0; //que salga entre el ecuador y el hemisferio norte
    v_inicial=11000.0/dTL; //la velocidad inicial del cohete debe ser la velocidad de escape (también en términos de la distancia Tierra-Luna)

    
    pphi=0.0;
    pr=v_inicial;


    m=20000.0; //masa del cohete
    H=0.0;
    rL=0;


    //Establezco constantes para facilitar los cálculos

    delta=7.014744e-12;  //Constante para simplificar los términos
    mu=0.0123025;                //Constante para simplificar los términos


    //Almaceno en el vector y las condiciones iniciales antes de comenzar el proceso iterativo

    y[0]=r;
    y[1]=phi;
    y[2]=pr;
    y[3]=pphi;


    pos_x=y[0]*cos(y[1]);
    pos_y=y[0]*sin(y[1]);


    //Empiezo el proceso iterativo

    fich_posiciones.open("planets_data.dat");
    Hamiltoniano.open("Hamiltoniano.txt");

    for (i=0; i<=iter; i++)
    {

        taux=0.0;

        //Calculo el Hamiltoniano
        rL=0.0;
        rL=sqrt(pow(dTL*y[0],2)+dTL*dTL-2*y[0]*dTL*dTL*cos(y[1]-w*t));
        H=0.0;
        H=pow(m*dTL*y[2],2)/(2*m)+pow(m*dTL*dTL*y[3],2)/(2*m*pow(dTL*y[0],2))-G*m*MT/(dTL*y[0])-G*m*ML/rL-w*(m*dTL*dTL*y[3]);



        if (y[0]>0.85 && y[0]<1.05)
        {
            h=100.0;
            //reiniciamos las iteraciones 
            i=0;
            iter=tmax/h;
        }

        //Escribo en un fichero de texto los datos de las posiciones para representarlas después

        fich_posiciones << x_T << ", " << y_T << endl;
        fich_posiciones << x_L << ", " << y_L << endl;
        fich_posiciones << pos_x << ", " << pos_y << endl;
        fich_posiciones << endl;

        Hamiltoniano << H << endl;



        //Calculo k1 para cada variable
        k1[0]=h*f1(y[2]);
        k1[1]=h*f2(y[3], y[0]);
        k1[2]=h*f3(y[3], y[0], y[1], t, delta, mu);
        k1[3]=h*f4(delta, mu, y[0], y[1], t);



        //Calculo valores auxiliares para facilitar la expresión de k2
        for (j=0; j<=3; j++)
        {
            aux1[j]=y[j]+k1[j]*0.5;
        }
        taux=t+h*0.5;



        //Calculo k2 para cada variable
        k2[0]=h*f1(aux1[2]);
        k2[1]=h*f2(aux1[3], aux1[0]);
        k2[2]=h*f3(aux1[3], aux1[0], aux1[1], taux, delta, mu);
        k2[3]=h*f4(delta, mu, aux1[0], aux1[1], taux);


        //Calculo valores auxiliares para facilitar la expresión de k3
        for (j=0; j<=3; j++)
        {
            aux2[j]=y[j]+k2[j]*0.5;
        }
        taux=t+h*0.5;



        //Calculo k3 para cada variable
        k3[0]=h*f1(aux2[2]);
        k3[1]=h*f2(aux2[3], aux2[0]);
        k3[2]=h*f3(aux2[3], aux2[0], aux2[1], taux, delta, mu);
        k3[3]=h*f4(delta, mu, aux2[0], aux2[1], taux);


        //Calculo valores auxiliares para facilitar la expresión de k4
        for (j=0; j<=3; j++)
        {
            aux3[j]=y[j]+k3[j];
        }
        taux=t+h;


        
        //Calculo k4 para cada variable
        k4[0]=h*f1(aux3[2]);
        k4[1]=h*f2(aux3[3], aux3[0]);
        k4[2]=h*f3(aux3[3], aux3[0], aux3[1], taux, delta, mu);
        k4[3]=h*f4(delta, mu, aux3[0], aux3[1], taux);


        //Calculo los nuevos valores de cada y
        for (j=0; j<=3; j++)
        {
            y[j]=y[j]+1.0/6.0*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j]);
        }


        //paso de coordenadas polares a cartesianas (posición del cohete)

        pos_x=y[0]*cos(y[1]);
        pos_y=y[0]*sin(y[1]);

        

        //establezco el tiempo para la siguiente iteración

        t=t+h;

        
        //Calculo la posición de la luna para que se esciba en el fichero en la siguiente iteración

        x_L=cos(w*t);
        y_L=sin(w*t);



    }

    fich_posiciones.close();
    Hamiltoniano.close();




    return 0;
}








//-----------------------------------------------FUNCIONES DEL PROGRAMA--------------------------------------//

double f1(double pr)
{
    double f1;

    f1=pr;

    return f1;

}



double f2(double pphi, double r)
{
    double f2;

    f2=pphi/pow(r,2);

    return f2;
}



double f3(double pphi, double r, double phi, double t, double delta, double mu)
{
    double f3, rprima;

    rprima=sqrt(1.0+r*r-2.0*r*cos(phi-w*t));

    f3=pow(pphi,2)/pow(r,3)-delta*(1/pow(r,2)+mu*(r-cos(phi-w*t))/pow(rprima,3));

    return f3;
}



double f4(double delta, double mu, double r, double phi, double t)
{
    double f4, rprima;

    rprima=sqrt(1.0+r*r-2.0*r*cos(phi-w*t));

    f4=-1.0*delta*mu*r/pow(rprima,3)*sin(phi-w*t);

    return f4;
}






 
