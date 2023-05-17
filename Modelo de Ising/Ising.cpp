#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#include "gsl_rng.h"

#define MAX 1000
//gsl_rng *tau;


double energia (double s[][MAX], int n, int m, int N);
double minimo (double E, double T);
void orden_matriz (int orden, double s[][MAX],int N);


using namespace std;


int main()
{
    int i,j,k,l; //contadores para los bucles
    int N;   //dimensión de la caja
    int iter; //para obtener el número de iteraciones
    double T; //temperatura
    double s[MAX][MAX]; //matriz de spines
    double delta_E=0.0; //valores de la energía de la iteración anterior y actual, y su diferencia
    int n,m; //posición inicial aleatoria
    double p, e; //valor de la probabilidad y un aleatorio para cambiar el signo del espin
    int orden; //variable para saber si la matriz se inicia de forma ordenada o desordenada

    //Fichero para escribir la matriz de espines en cada iteración
    ofstream fich_s;


    //----------------------GENERADOR DE NÚMEROS ALEATORIOS----------------------------

    //extern gsl_rng *tau;
    //int semilla;

    //semilla=1482673;
    //tau=gsl_rng_alloc(gsl_rng_taus);
    //gsl_rng_set(tau,semilla);

    srand(time(NULL));


    //---------------------------------------------------------------------------------
    //Pregunto al usuario si quiere que la matriz empiece ordenada o desordenada
    cout << "Indique si la matriz comienza ordenada (0) o desordenada (1)" << endl;
    cin >> orden;

    //Pido el valor de la temperatura, la dimensión de la red y los pasos del sistema
    cout << "Valor de la temperatura: " << endl;
    cin >> T ;
    

    cout << "Dimensión de la red: " << endl;
    cin >> N;

    cout << "Numero de iteraciones: " << endl;
    cin >> iter;



    //Inicio todos los espines (configuración ordenada (0) o desordenada(1) )
    orden_matriz (orden, s, N);


    //Inicio el bucle del modelo

    fich_s.open("ising_data.dat");

    for (i=0; i<=iter; i++)
    {
        // Escribo la matriz de espines en un fichero de texto
        for(k=0; k<N; k++)
        {
            for(l=0;l<N-1;l++)
            {
                fich_s << s[k][l] << ", ";
            }
            fich_s << s[k][N-1] << endl;
        }
        fich_s << endl;


        for(j=0; j<N*N; j++)
        {

            //Me "coloco" en una posición aleatoria de la caja 
            n=rand()%(N);
            m=rand()%(N);

        
            //Calculo la variación de la energía
            delta_E=energia(s, n, m, N);

            //Obtengo el valor de p
            p=minimo(delta_E, T);

            //Genero un número aleatorio entre 0 y 1
            e=(double)rand()/RAND_MAX;

            if (e<p)
            {
                s[n][m]=-s[n][m];
            }

        }


    }
    fich_s.close();

    

    return 0;
}









//---------------------------------FUNCIONES DEL SISTEMA----------------------------------
double energia (double s[][MAX], int n, int m, int N)
{
    int i;
    double E=0.0;

    if(n==N-1 && m!=N-1 && m!=0)
    {
        E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    }

    else if(n==0 && m!=N-1 && m!=0)
    {
        E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][m-1]);
    }

    else if (n==N-1 && m==N-1)
    {
        E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
    }

    else if (n==0 && m==N-1)
    {
        E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][0]+s[n][m-1]);
    }

    else if (n==N-1 && m==0)
    {
        E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
    }

    else if (n==0 && m==0)
    {
        E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][N-1]);
    }

    else if (n!=N-1 && n!=0 && m==N-1)
    {
        E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
    }

    else if (n!=N-1 && n!=0 && m==0)
    {
        E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
    }

    else 
    {
        E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    }


    return E;

}


double minimo (double E, double T)
{
    if (exp(-E/T)<1.0)
    { 
        return exp(-E/T);
    }

    else 
    {
        return 1.0;
    }
}


void orden_matriz (int orden, double s[][MAX],int N)
{
    int pot_aleatoria;
    int i,j;

    if(orden==0)
    {
        pot_aleatoria=rand();

        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                s[i][j]=pow(-1,pot_aleatoria);
            }
        }

    }

    else if(orden==1)
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                s[i][j]=pow(-1,rand());
            }
        }
    }

    return;
}