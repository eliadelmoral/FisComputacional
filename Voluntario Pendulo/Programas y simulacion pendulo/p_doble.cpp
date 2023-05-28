#include <iostream>
#include <cmath>
#include <fstream>


#define g 9.807
#define PI 3.14159265



double f1(double O1, double O2, double pO1, double pO2);
double f2(double O1, double O2, double pO1, double pO2);
double f3(double O1, double O2, double pO1, double pO2);
double f4(double O1, double O2, double pO1, double pO2);




using namespace std;

int main()
{   
    double h, tmax; //paso y variable temporal
    double y[4], k1[4], k2[4], k3[4], k4[4]; //4 ecuaciones diferenciales y las funciones ki para aplicar el algoritmo
    double aux1[4], aux2[4], aux3[4]; //vectores auxiliares para facilitar el cálculo de ki
    double O1, O2, pO1, pO2;

    
    double pos_x1, pos_y1, pos_x2, pos_y2; //para transformar las coordenadas polares en cartesianas para representarlo

    
    int i, j, iter;
    

    ofstream fich_posiciones;
    ofstream fich_O1_vO1;
    ofstream fich_O2_vO2;


    //Condiciones iniciales

    h=0.1;
    tmax=100.0;
    iter=tmax/h;

    O1=PI/2.0;
    O2=PI/2.0;
    pO1=0.0;
    pO2=0.0;


    //Almaceno en el vector y las condiciones iniciales antes de comenzar el proceso iterativo

    y[0]=O1;
    y[1]=O2;
    y[2]=pO1;
    y[3]=pO2;


    //Posiciones de las masas en coordenadas cartesianas
    pos_x1=sin(y[0]);
    pos_y1=-1.0*cos(y[0]);
    pos_x2=sin(y[1])+sin(y[0]);
    pos_y2=-1.0*cos(y[1])-1.0*cos(y[0]);



    //Empiezo el proceso iterativo

    fich_posiciones.open("planets_data.dat");
    

    for (i=0; i<=iter; i++)
    {

        //Escribo en un fichero de texto los datos de las posiciones para representarlas después

        fich_posiciones << pos_x1 << ", " << pos_y1 << endl;
        fich_posiciones << pos_x2 << ", " << pos_y2 << endl;
        fich_posiciones << endl;



        //Calculo k1 para cada variable
        k1[0]=h*f1(y[0], y[1], y[2], y[3]);
        k1[1]=h*f2(y[0], y[1], y[2], y[3]);
        k1[2]=h*f3(y[0], y[1], y[2], y[3]);
        k1[3]=h*f4(y[0], y[1], y[2], y[3]);



        //Calculo valores auxiliares para facilitar la expresión de k2
        for (j=0; j<=3; j++)
        {
            aux1[j]=y[j]+k1[j]*0.5;
        }



        //Calculo k2 para cada variable
        k2[0]=h*f1(aux1[0], aux1[1], aux1[2], aux1[3]);
        k2[1]=h*f2(aux1[0], aux1[1], aux1[2], aux1[3]);
        k2[2]=h*f3(aux1[0], aux1[1], aux1[2], aux1[3]);
        k2[3]=h*f4(aux1[0], aux1[1], aux1[2], aux1[3]);


        //Calculo valores auxiliares para facilitar la expresión de k3
        for (j=0; j<=3; j++)
        {
            aux2[j]=y[j]+k2[j]*0.5;
        }



        //Calculo k3 para cada variable
        k3[0]=h*f1(aux2[0], aux2[1], aux2[2], aux2[3]);
        k3[1]=h*f2(aux2[0], aux2[1], aux2[2], aux2[3]);
        k3[2]=h*f3(aux2[0], aux2[1], aux2[2], aux2[3]);
        k3[3]=h*f4(aux2[0], aux2[1], aux2[2], aux2[3]);


        //Calculo valores auxiliares para facilitar la expresión de k4
        for (j=0; j<=3; j++)
        {
            aux3[j]=y[j]+k3[j];
        }


        
        //Calculo k4 para cada variable
        k4[0]=h*f1(aux3[0], aux3[1], aux3[2], aux3[3]);
        k4[1]=h*f2(aux3[0], aux3[1], aux3[2], aux3[3]);
        k4[2]=h*f3(aux3[0], aux3[1], aux3[2], aux3[3]);
        k4[3]=h*f4(aux3[0], aux3[1], aux3[2], aux3[3]);


        //Calculo los nuevos valores de cada y
        for (j=0; j<=3; j++)
        {
            y[j]=y[j]+1.0/6.0*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j]);
        }



        //paso de coordenadas polares a cartesianas

        pos_x1=sin(y[0]);
        pos_y1=-1.0*cos(y[0]);
        pos_x2=sin(y[1])+sin(y[0]);
        pos_y2=-1.0*cos(y[1])-1.0*cos(y[0]);
        

    }

    fich_posiciones.close();




    return 0;
}








//-----------------------------------------------FUNCIONES DEL PROGRAMA--------------------------------------//

double f1(double O1, double O2, double pO1, double pO2)
{
    double f1;

    f1=(pO1-1.0*pO2*cos(O1-O2))/(1.0+pow(sin(O1-O2),2));

    return f1;

}



double f2(double O1, double O2, double pO1, double pO2)
{
    double f2;

    f2=(-1.0*pO1*cos(O1-O2)+2.0*pO2)/(1.0+pow(sin(O1-O2),2));

    return f2;
}



double f3(double O1, double O2, double pO1, double pO2)
{
    double f3,h1,h2;

    h1=pO1*pO2*sin(O1-O2)/(1.0+pow(sin(O1-O2),2));
    h2=(pO1*pO1+2.0*pO2*pO2-2.0*pO1*pO2*cos(O1-O2))/(2*pow(1.0+pow(sin(O1-O2),2),2));
    f3=-2.0*g*sin(O1)-h1+h2*sin(2*(O1-O2));

    return f3;
}


double f4(double O1, double O2, double pO1, double pO2)
{
    double f4,h1,h2;

    h1=pO1*pO2*sin(O1-O2)/(1.0+pow(sin(O1-O2),2));
    h2=(pO1*pO1+2.0*pO2*pO2-2.0*pO1*pO2*cos(O1-O2))/(2*pow(1.0+pow(sin(O1-O2),2),2));
    f4=-1.0*g*sin(O2)+h1-h2*sin(2*(O1-O2));

    return f4;
}