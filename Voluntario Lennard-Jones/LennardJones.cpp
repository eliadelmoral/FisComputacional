#include <iostream>
#include <fstream>
#include <cmath>
#include "gsl_rng.h"

#define MAX 100


double modulo (double r1, double r2);
void aceleracion (double m[MAX], double r1[MAX], double r2[MAX], double a[MAX], int atomos);
void posicion (double h,double r[MAX], double v[MAX], double a[MAX], int atomos);
void velocidad (double h, double w[MAX], double a[MAX], double v[MAX], int atomos);
void funcion_w (double h, double v[MAX], double a[MAX], double w[MAX], int atomos);
double potencial (double r1[MAX], double r2[MAX], double vx[MAX], double vy[MAX], int atomos);
double cinetica (double  v1[MAX], double v2[MAX], int atomos);
double energia (double cin, double pot);



using namespace std;

int main ()
{
    int atomos;
    int i, iter;
    double tmax, t;
    double h=0.002;
    double caja[10][10];

    t=0;

    //Pido al usuario que indique el número de partículas
    cout << "Numero de particulas: " << endl;
    cin >> atomos;


    //Pido al usuario que establezca un tiempo máximo para limitar las iteraciones, además del paso entre iteraciones
    cout << "Establezca un tiempo máximo" << endl;
    cin >> tmax;


    //Establezco el número de iteraciones
    iter=(int)tmax/h;




    //Declaro e inicializo los vectores a 0

    double x[atomos-1];
    double y[atomos-1];      
    double vx[atomos-1];     
    double vy[atomos-1];
    double angulo[atomos-1];
    double ax[atomos-1]={0};
    double ay[atomos-1]={0};
    double m[atomos-1]={1.0};
    double wx[atomos-1]={0};
    double wy[atomos-1]={0};
    double pot=0.0;
    double cin=0.0;
    double E=0.0;
 

    //Declaro ficheros donde guardar los datos de las posiciones, de la energía y del momento angular
    ofstream fich_posi;
    ofstream fich_energia_potencial;
    ofstream fich_energia_cinetica;
    ofstream fich_energia;

    
    //Hago que en cada iteración escriba la posición en x e y de todos los planetas, y la energía en función del tiempo

    fich_posi.open("planets_data.dat");
    fich_energia_potencial.open("energia_potencial.txt");
    fich_energia_cinetica.open("energia_cinetica.txt");
    fich_energia.open("energia.txt");
   
    for (i=0; i<=iter; i++)
    {
        E=0;
        pot=0;
        cin=0;

        

        //Calculo las energías
        pot=potencial(x, y, vx, vy, atomos);
        cin=cinetica( vx, vy, atomos);
        E=energia(cin, pot);


        //Completo el fichero de las energía en función del tiempo
        fich_energia << t << "  " << E << endl;
        fich_energia_cinetica << t << " " << cin << endl;
        fich_energia_potencial << t << " " << pot << endl;



        //Completo el fichero de las posiciones

        for (i=0; i<atomos; i++)
        {
            fich_posi << x[i] << " , " << y[i] << endl;
        }

        

        //Calculo la posicion
        posicion (h, x, vx, ax, atomos);
        posicion (h, y, vy, ay, atomos);


        //Calculo la funcion w
        funcion_w (h, vx, ax, wx, atomos);
        funcion_w (h, vy, ay, wy, atomos);

        for (i=0; i<atomos; i++)
        {
            ax[i]=0;
            ay[i]=0;
        }
        
        //Calculo la aceleración 
        aceleracion (m, x, y, ax, atomos);
        aceleracion (m, y, x, ay, atomos);



        //Calculo la velocidad
        velocidad (h, wx, ax, vx, atomos);
        velocidad (h, wy, ay, vy, atomos);


        t=t+h;


    }

    fich_posi.close();
    fich_energia.close();
    fich_energia_cinetica.close();
    fich_energia_potencial.close();


    return;
}






//--------------------------FUNCIONES DEL PROGRAMA--------------------------------------------//

//Cálculo de la temperatura a partir del teorema de equipartición
//void teorema_equipartición


//Función para calcular la energía potencial

double potencial (double r1[MAX], double r2[MAX], double vx[MAX], double vy[MAX], int atomos)
{
    int i, j;
    double pot;

    pot=0.0;

        //Potencial de Lennard-Jones, cada partícula interacciona con todas las demás. Se va superponiendo
        for (i=1; i<atomos; i++)

        {
            for (j=0; j<atomos; j++)
             {
                if(i!=j)
                {
                    pot=pot+4*(1/pow(modulo(r1[i]-r1[j],r2[i]-r2[j]), 12)-pow(modulo(r1[i]-r1[j],r2[i]-r2[j]), 6)) ;
                }

             }
        }
    
        
    return pot;  
}

//Función para calcular la energía cinética

double cinetica (double  v1[MAX], double v2[MAX], int atomos)
{
    int i;
    double cin;

    cin=0.0;
    for (i=0; i<atomos; i++)
    {
        cin=cin+0.5*(v1[i]*v1[i]+v2[i]*v2[i]);
    }

    return cin;
}


//Función para calcular la energía total

double energia (double cin, double pot)
{
    double E=0.0;
    E=cin+pot;

    return E;
}


//-------------------------FUNCIONES PARA EL ALGORITMO DE VERLET-------------------------------//


//Función para obtener el módulo de la posición, facilitando el cálculo de la aceleración

double modulo (double r1, double r2)
{
    double modulo;
    modulo=0.0;

    modulo=sqrt(r1*r1+r2*r2);
    
    return modulo;
    
}


//Función para obtener la aceleración

void aceleracion (double m[MAX], double r1[MAX], double r2[MAX], double a[MAX], int atomos)
{
    //Declaro dos enteros para hacer dos bucles 
    int i,j;

    //Hago un bucle que recorra todo el vector y otro bucle en su interior que opere con todos los elementos del vector de índice distinto

    for (i=0; i<atomos; i++)
    {
        for (j=0; j<atomos-1; j++)
        {
            if (j!=i)
            {
                a[i]=a[i]-(m[j])*(r1[i]-r1[j])/pow(modulo(r1[i]-r1[j], r2[i]-r2[j]),3);
            }
        
        }  
   
    }

    return;  
}



//Función para calcular la posición en tiempo h

void posicion (double h,double r[MAX], double v[MAX], double a[MAX], int atomos)
{
    int i;
    
    for (i=0; i<atomos; i++)
    {
        r[i]=r[i]+h*v[i]+pow(h,2)*a[i]/2.0;
    }

    return;
}


//Función para calcular la velocidad en tiempo h

void velocidad (double h, double w[MAX], double a[MAX], double v[MAX], int atomos)
{
    int i;

    for (i=0; i<atomos; i++)
    {
        v[i]=w[i]+h*a[i]/2.0;
    }

    return;
}




//Función para definir la función w

void funcion_w (double h, double v[MAX], double a[MAX], double w[MAX], int atomos)
{
    int i;

    for (i=0; i<atomos; i++)
    {
        w[i]=v[i]+h*a[i]/2.0;
    
    }

    return;
}
