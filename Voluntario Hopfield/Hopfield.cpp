#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#define MAX 1000

void orden_matriz (int orden, double s[][MAX],int N);
double minimo (double E, double T);


using namespace std;



int main ()
{
    int i,j,k,l,z,x;
    double s[MAX][MAX], patrones[MAX][MAX][MAX];
    double w[MAX][MAX][MAX][MAX], umbral[MAX][MAX];
    double a[MAX];
    int recuerda;
    double T, H, e;
    int N, P, mu, p;
    int orden, iter;
    int n,m;
    double solapamiento[MAX];

    ofstream fich_s;
    ofstream fich_solapamiento;
    ofstream fich_serecuerda;

    srand(time(NULL));

    //Los patrones valdrán 0 o 1 aleatoriamente

    for (mu=0; mu<P; mu++)
    {
        for (i=0; i<N; i++)
        {
            for (j=0; j<N; j++)
            {
                patrones[i][j][mu]=rand()%(2);
            }
        }

    }
    

    //Calculo a (proporcion de píxeles blancos/negros)

    for(i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            a[mu]=a[mu]+patrones[i][j][mu];
        }
    }
    
    a[mu]=a[mu]/pow(N,2);



    //Calculo w (pesos sinápticos)

    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            for (k=0; k<N; k++)
            {
                for (l=0; l<N; l++)
                {
                    if (i==j && k==l)
                    {
                        w[i][j][k][l]=0.0;
                    }
                    else 
                    {
                        for (mu=0; mu<P; mu++)
                        {
                            w[i][j][k][l]=w[i][j][k][l]+(patrones[i][j][mu]-a[mu])*(patrones[k][l][mu]-a[mu]);
                        }
                        w[i][j][k][l]=w[i][j][k][l]/pow(N,2);
                    }
                    
                }
            }
        }
    }



    //Calculo theta (umbral de disparo)

    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            for (k=0; k<N; k++)
            {
                for (l=0; l<N; l++)
                {
                    umbral[i][j]=umbral[i][j]+0.5*w[i][j][k][l];
                }
            }
        }
    }






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



    //Inicio todos los espines (configuración ordenada o desordenada)
    orden_matriz (orden, s, N);


    //Inicio el bucle del modelo

    fich_s.open("ising_data.dat");
    fich_solapamiento.open("solapamiento.txt");
    fich_serecuerda.open("recuerda.txt");


    recuerda=0;


    for (z=0; z<=iter; z++)
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


        //Escribo el solapamiento en un fichero de texto
        for(mu=0; mu<P; mu++)
        {
            fich_solapamiento << solapamiento[mu] << endl;
        }


        //Escribo el número de patrones recordados en un fichero de texto
        for(mu=0; mu<P; mu++)
        {
            fich_serecuerda << recuerda;
        }

       

        for(x=0; x<N*N; x++)
        {

            //Me "coloco" en una posición aleatoria de la caja (condición inicial aleatoria)
            n=rand()%(N);
            m=rand()%(N);

        
            //Calculo la variación de la energía
            H=0.0;

            for(k=0; k<N; k++)
            {
                for (j=0; j<N; j++)
                {
                    if(k!=n && j!=m)
                    {
                        H=H+0.5*(w[n][m][k][l]*(double)s[n][m]);
                    }
                }
            }

            H=umbral[n][m]*(double)s[n][m]-H;
            
            


            //Obtengo el valor de p
            p=minimo(H, T);

            //Genero un número aleatorio entre 0 y 1
            e=(double)rand()/RAND_MAX;

            if (e<p)
            {
                s[n][m]=1-s[n][m];
            }



            //Calculo el solapamiento
            
            for(mu=0; mu<P; mu++)
            {
                for(i=0; i<N; i++)
                {   
                    for (j=0; j<N; j++)
                    {
                        solapamiento[mu]=solapamiento[mu]+(patrones[i][j][mu]-a[mu])*((double)s[i][j]-a[mu]);
                    }   
                }
                solapamiento[mu]=solapamiento[mu]/((double)N*(double)N*a[mu]*(1.0-a[mu]));
            }

            //¿Se recuerda el patrón?

            for(mu=0; mu<P; mu++)
            {
                if((double)abs(solapamiento[mu]>0.75))
                {
                    recuerda=recuerda+1;
                }
            }


        }


    }
    fich_s.close();
    fich_solapamiento.close();
    fich_serecuerda.close();

    

    return 0;
}




//----------------------FUNCIONES DEL PROGRAMA---------------------------//


//Función orden, para que el usuario indique si quiere empezar con un patrón ordenado o desordenado
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




//Función mínimo, que devuelve el mínimo entre dos valores escogidos
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

