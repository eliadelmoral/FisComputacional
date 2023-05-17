#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>


#define MAX 1000
#define PI 3.14159265




using namespace std;

int main()
{
    int j, n, t;
    double S;
    double V[MAX];
    int N;
    float n_ciclos;
    double  K_0; 
    float lambda;
    complex<double> phi[MAX], beta[MAX], alpha[MAX], chi[MAX], b[MAX];
    double norma; 
    double sigma, x_0; //Centro y anchura de la gaussiana




    //ficheros de texto para guardar los resultados

    ofstream fich_funcion;
    ofstream fich_norma;


    complex<double> im_puro (0.0,1.0);
    complex<double> A0[MAX];
 

    //--------------------CONDICIONES INICIALES Y CONDICIONES DE CONTORNO-----------------//

    //Establezco el número de iteraciones, la longitud de onda, el número de ciclos, s y h

    N=1000;
    lambda=5.0;
    n_ciclos=(double)N/16.0; //Si n_ciclos=N/4 la frecuencia es demasiado alta y no se distinguen bien la parte real e imaginaria
    t=1000; 



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





    //Función de onda inicial
    //---------Condiciones de contorno----------//

    phi[0]=0.0;
    phi[N]=0.0;

    //---------Función de onda inicial----------//

    norma=0.0;

    sigma=(double)N/16.0; //Anchura de la gaussiana
    x_0=(double)N/4.0; //Centro de la gaussiana
    


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



    //------------------PROCESO ITERATIVO---------------------//

    fich_funcion.open("schrodinger_data.dat");
    fich_norma.open("norma.txt");



    for (n=0; n<t; n++)
    {

        //Escribo cada uno de los ficheros de texto

        for (j=0; j<N; j++)
        {
            fich_funcion << j << ", " << real(phi[j]) << ", " << imag(phi[j]) << ", " << pow(abs(phi[j]),2) << endl;
        }
        fich_funcion << endl;
        


        //Calculo b, beta (iteraciones hacia atrás al igual que alpha)


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
            norma=norma+pow(abs(phi[j]),2);

        }
        fich_norma << norma << endl;


    }

    fich_funcion.close();
    fich_norma.close();




    return 0;
}










