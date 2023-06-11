#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <time.h>
#include "gsl_rng.h"



#define MAX 2500
#define PI 3.14159265

gsl_rng *tau;


using namespace std;


int main()
{
    int n, j;
    int n_D;
    double m_T;
    double K, K_teo;
    double p, P_d, P_i, P_d_max;


    double S;
    double V[MAX];
    int N;
    float n_ciclos;
    double  K_0; 
    float lambda;
    complex<double> phi[MAX], beta[MAX], alpha[MAX], chi[MAX], b[MAX];
    double norma; 
    double sigma, x_0; //Centro y anchura de la gaussiana


    
    int Nexp, t; //número de experimentos y variable temporal
    double aux; //variable auxiliar para determinar si se ha detectado la partícula



    //generador de números aleatorios
    extern gsl_rng *tau;
    int semilla=1297534;

    tau=gsl_rng_alloc(gsl_rng_taus); //puntero
    gsl_rng_set(tau,semilla);
    





    //ficheros de texto para guardar los resultados

    ofstream fich_K;
    ofstream fich_probabilidad;



    complex<double> im_puro (0.0,1.0);
    complex<double> A0[MAX];
 

    //--------------------CONDICIONES INICIALES Y CONDICIONES DE CONTORNO-----------------//

    //Establezco el número de iteraciones, la longitud de onda, el número de ciclos, s y h

    N=1000; //tamaño del sistema 
    lambda=0.5;
    n_ciclos=(double)N/16.0; //Si n_ciclos=N/4 la frecuencia es demasiado alta y no se distinguen bien la parte real e imaginaria
    t=0; 
    Nexp=1000;  //número de experimentos
    


    //Calculo el coeficiente de transmisión teórico
    
    if(lambda<1)
    {
        K_teo=4*(1-lambda)/(4*(1-lambda)+pow(lambda,2)*pow(sin(2*PI/5*n_ciclos*sqrt(1-lambda)),2));
    }
    
    else if(lambda>1)
    {
        K_teo=4*(lambda-1)/(4*(lambda-1)+pow(lambda,2)*pow(sinh(2*PI/5*n_ciclos*sqrt(lambda-1)),2));
    }

    else if(lambda==1)
    {
        K=0.0;
    }
    
    


    //Calculo los parámetros iniciales

    K_0=(2.0*PI*n_ciclos)/(double)N;
    S=1.0/(4.0*K_0*K_0);


   //Potencial
    
    for (j=0; j<N; j++)
    {
        if (j>0.4*(double)N && j<0.6*(double)N)
        {
            V[j]=lambda*K_0*K_0;
        }

        else V[j]=0.0;
        
    }


    //Variable compleja A0

    for (j=0; j<N; j++)
    {
        A0[j]=-2.0+2.0*im_puro/S-V[j];
    }


    //Obtengo alpha (iteraciones hacia atrás)

    
    alpha[N-1]=0.0;
    

    for (j=0; j<N-1; j++)
    {
        alpha[(N-2)-j]=-1.0/(A0[(N-1)-j]+alpha[(N-1)-j]);
    }

    

    sigma=(double)N/16.0; //Anchura de la gaussiana
    x_0=(double)N/4.0; //Centro de la gaussiana
    



    //------------------PROCESO ITERATIVO---------------------//



    fich_K.open("coef_trans.txt");
    fich_probabilidad.open("P_D.txt");

    m_T=0.0;
    n_D=500; //establecewmos el máximo aproximadamente en N/2
    P_d=0.0;
    P_i=0.0;
    P_d_max=0.0;


    for (n=0; n<Nexp; n++)
    {

        //Función de onda 
        //---------Condiciones de contorno----------//

        phi[0]=0.0;
        phi[N]=0.0;

        //---------Función de onda ----------//

        norma=0.0;

        for (j=1; j<N-1; j++)
        {
            phi[j]=exp(im_puro*(double)j*K_0)*exp(-1.0*(pow((double)j-x_0,2)/(2*pow(sigma,2))));

            norma=norma+pow(abs(phi[j]),2);

        }

        //Normalizo la función de onda
        
        for (j=0; j<N; j++)
        {
            phi[j]=phi[j]/sqrt(norma);
        }



        //--------------------------CÁLCULOS----------------------------//

        t=0;
        aux=0;

        while (aux==0) //mientras la partícula todavía se encuentre a la izquierda y no se haya alcanzado el tiempo máximo de iteraciones
        {
            t=t+1;

            //Calculo b, beta

            for (j=0; j<N; j++)
            {
                b[j]=4.0*im_puro*phi[j]/S;
            }


            beta[N-1]=0.0;

            for (j=0; j<N-1; j++)
            {
                beta[(N-2)-j]=(b[(N-1)-j]-beta[(N-1)-j])/(A0[(N-1)-j]+alpha[(N-1)-j]);
            }


            //Calculo chi (iteraciones hacia alante)

            chi[0]=0.0;
            chi[N-1]=0.0;

            for (j=0; j<N-1; j++)
            {
                chi[j+1]=alpha[j]*chi[j]+beta[j];
            }


            //Calculamos la función de onda y su norma (para comprobar que se conserva)

            norma=0.0;


            for(j=0; j<N; j++)
            {
                phi[j]=chi[j]-phi[j];

            }



            //----------------COEFICIENTE DE TRANSMISIÓN-----------------------//

                    

            //Calculo la probabilidad de encontrar la partícula a la derecha 

            if(t%n_D==0) //mientras no se encuentre el máximo local
            {
                P_d=0.0;


                for(j=4.0*N/5.0; j<N; j++)
                {
                    P_d=P_d+pow(abs(phi[j]),2);
                }
                P_d_max=P_d_max+P_d;
                

                //Genero un número aleatorio entre 0 y 1 para ver si hemos detectado la partícula
                p=gsl_rng_uniform(tau);

                if (p<P_d)
                {
                    m_T=m_T+1.0;
                    aux=1; //¡aquí se acabaría el bucle!
                }  
                    

                //Si no la detectamos a la derecha, la probabilidad a la derecha colapsa a  y la función de onda se redistribuye en el resto del espacio
                else
                {
                    for(j=4.0*N/5.0; j<N; j++)
                    {
                        phi[j]=0.0;
                    }

                    norma=0.0;

                    //Calculamos la función de onda en el resto del espacio porque cambia el estado físico del sistema

                    for(j=0; j<N; j++)
                    {
                        norma=norma+pow(abs(phi[j]),2);
                    }

                    for(j=0; j<N; j++)
                    {
                        phi[j]=phi[j]/sqrt(norma);
                    }


                    //Veamos ahora si se detecta la partícula en el detector izquierdo
                    P_i=0.0;

                    for(j=0; j<1.0*N/5.0; j++)
                    {
                        P_i=P_i+pow(abs(phi[j]),2);
                    }

                    //Genero un número aleatorio entre 0 y 1 para ver si hemos detectado la partícula
                    p=gsl_rng_uniform(tau);

                    if(p<P_i)
                    {
                        aux=1; //también se acabaría el bucle
                    }

                    //Si no se ha encontrado, vuelve a cambiar el estado físico del sistema

                    else
                    {
                        for(j=0; j<1.0*N/5.0; j++)
                        {
                            phi[j]=0.0;
                        }

                        norma=0.0;

                        //Calculamos la función de onda en el resto del espacio porque cambia el estado físico del sistema

                        for(j=0; j<N; j++)
                        {
                            norma=norma+pow(abs(phi[j]),2);
                        }

                        for(j=0; j<N; j++)
                        {
                            phi[j]=phi[j]/sqrt(norma);
                        }
                    }


                }
            }
        }

            

            
    }

    //Calculo el coeficiente de transmisión 
    K=(double)m_T/Nexp;

    fich_K << K << endl;
    fich_K << K_teo << endl;
    fich_K << endl;


    P_d_max=P_d_max/Nexp;

    fich_probabilidad << P_d_max << endl;



    fich_K.close();
    fich_probabilidad.close();




    return 0;
}



