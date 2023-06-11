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


    complex<double> posi;
    complex<double> momento;
    complex<double> Ecinetica;
    complex<double> Etotal;
    complex<double> error_posi;
    complex<double> error_momento;
    complex<double> error_Ecinetica;
    complex<double> error_Etotal;
    complex<double> phi1[MAX], phi2[MAX];



    
    int Nexp, t, tmax; //número de experimentos, tiempo que se realizan los experimentos y variable temporal
    double aux; //variable auxiliar para determinar si se ha detectado la partícula



    //generador de números aleatorios
    extern gsl_rng *tau;
    int semilla=1297534;

    tau=gsl_rng_alloc(gsl_rng_taus); //puntero
    gsl_rng_set(tau,semilla);
    





    //ficheros de texto para guardar los resultados


    ofstream fich_posi;
    ofstream fich_Ecinetica;
    ofstream fich_Etotal;
    ofstream fich_momento;




    complex<double> im_puro (0.0,1.0);
    complex<double> A0[MAX];
 

    //--------------------CONDICIONES INICIALES Y CONDICIONES DE CONTORNO-----------------//

    //Establezco el número de iteraciones, la longitud de onda, el número de ciclos, s y h

    N=1000; //tamaño del sistema 
    lambda=0.5;
    n_ciclos=(double)N/16.0; //Si n_ciclos=N/4 la frecuencia es demasiado alta y no se distinguen bien la parte real e imaginaria
    t=0; 
    Nexp=1000;  //número de experimentos
    tmax=1000; //tiempo máximo de iteraciones
    


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




    fich_posi.open("posicion.txt");
    fich_Ecinetica.open("Ecinetica.txt");
    fich_Etotal.open("Etotal.txt");
    fich_momento.open("momento.txt");



    n_D=500; //establecewmos el máximo aproximadamente en N/2
    P_d=0.0;
    P_i=0.0;



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

        for(t=0; t<tmax; t++)
        {

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




            //-----------------------VALORES ESPERADOS---------------------------//


            //Calculo los valores esperados de la posición, el momento, la energía cinética y la energía total
            //Para ello necesito calcular los valores de la primera y la segunda derivada de la función de onda 
            //Estas derivadas las aproximo como variaciones de la función de onda entre iteraciones

            //Conduciones de contorno
            phi1[0]=phi[1]-phi[0];
            phi1[N]=phi[N]-phi[N-1];

            for (j=1; j<N-1; j++)
            {
                phi1[j]=(phi[j+1]-phi[j-1])/2.0;
            }



            phi2[0]=phi1[1]-phi1[0];
            phi2[N]=phi1[N]-phi1[N-1];

            for (j=1; j<N-1; j++)
            {
                phi2[j]=(phi1[j+1]-phi1[j-1])/2.0;
            }


            //Valores esperados

            posi=0.0;
            momento=0.0;
            Ecinetica=0.0;
            Etotal=0.0;

            error_posi=0.0;
            error_momento=0.0;
            error_Ecinetica=0.0;
            error_Etotal=0.0;
            


            for(j=0; j<N; j++)
            {
                //Valores esperados
                posi=posi+j*pow(abs(phi[j]),2);
                momento=momento+conj(phi[j])*phi1[j];
                Ecinetica=Ecinetica-conj(phi[j])*phi2[j];
                Etotal=Etotal-conj(phi[j])*phi2[j]+V[j]*pow(abs(phi[j]),2);


                //Incertidumbres
                error_posi=error_posi+j*j*pow(abs(phi[j]),2);
                error_momento=error_momento+pow(conj(phi[j])*phi1[j],2);
                error_Ecinetica=error_Ecinetica-pow(conj(phi[j])*phi2[j],2);
                error_Etotal=error_Etotal+pow(-conj(phi[j])*phi2[j]+V[j]*pow(abs(phi[j]),2),2);

            }
            momento=momento*(-1.0)*im_puro;

            error_posi=sqrt(error_posi-posi*posi);
            error_momento=sqrt(error_momento-momento*momento);
            error_Ecinetica=sqrt(error_Ecinetica-Ecinetica*Ecinetica);
            error_Etotal=sqrt(error_Etotal-Etotal*Etotal);


        
            //-------------------------FICHEROS DE TEXTO-----------------------------//

            fich_posi << t << " " << real(posi) << "    " << real(error_posi) << endl;
            fich_momento << t << "  " << real(momento) << "   " << real(error_momento) << endl;
            fich_Ecinetica << t << "    " << sqrt(pow(real(Ecinetica),2)+pow(imag(Ecinetica),2)) << "    " << sqrt(pow(real(error_Ecinetica),2)+pow(imag(error_Ecinetica),2)) << endl;
            fich_Etotal << t << "   " << sqrt(pow(real(Etotal),2)+pow(imag(Etotal),2)) << "    " << sqrt(pow(real(error_Etotal),2)+pow(imag(error_Etotal),2)) << endl;


            //Hacemos la corrección en nD, donde se encuentra el máximo de probabilidad, puesto que cambia el comportamiento del sistema

            if(t==n_D)
            {
                P_d=0.0;


                for(j=4.0*N/5.0; j<N; j++)
                {
                    P_d=P_d+pow(abs(phi[j]),2);
                }
                    

                //Si no la detectamos a la derecha, la probabilidad a la derecha colapsa a  y la función de onda se redistribuye en el resto del espacio
                if (p>P_d)
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


        fich_posi.close();
        fich_Ecinetica.close();
        fich_Etotal.close();
        fich_momento.close();    


    return 0;
}