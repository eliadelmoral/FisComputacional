#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <time.h>


#define MAX 1000
#define PI 3.14159265


using namespace std;


int main()
{
    int j, n, t;
    int m_T, n_D;
    int K;
    double p, P_d, P_d_antes, P_d_max;


    double S;
    double V[MAX];
    int N;
    float n_ciclos;
    double  K_0; 
    float lambda;
    complex<double> phi[MAX], beta[MAX], alpha[MAX], chi[MAX], b[MAX];
    double norma; 
    double sigma, x_0; //Centro y anchura de la gaussiana


    double posi;
    complex<double> momento;
    double Epotencial;
    complex<double> Ecinetica;
    complex<double> Etotal;
    complex<double> phi1[MAX], phi2[MAX];
    



    srand(time(NULL));




    //ficheros de texto para guardar los resultados

    ofstream fich_funcion;
    ofstream fich_norma;
    ofstream fich_K;
    ofstream fich_Pdmax;
    ofstream fich_posi;
    ofstream fich_energias;
    ofstream fich_momento;
    ofstream fich_probabilidad;



    complex<double> im_puro (0.0,1.0);
    complex<double> A0[MAX];
 

    //--------------------CONDICIONES INICIALES Y CONDICIONES DE CONTORNO-----------------//

    //Establezco el número de iteraciones, la longitud de onda, el número de ciclos, s y h

    N=500;
    lambda=0.1;
    n_ciclos=(double)N/16.0; //Si n_ciclos=N/4 la frecuencia es demasiado alta y no se distinguen bien la parte real e imaginaria
    t=500; 



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
    fich_K.open("coef_trans.txt");
    fich_Pdmax.open("Pdmax.txt");
    fich_posi.open("posicion.txt");
    fich_energias.open("energias.txt");
    fich_momento.open("momento.txt");
    fich_probabilidad.open("P_D.txt");

    m_T=0;
    n_D=0;


    for (n=0; n<t; n++)
    {

        //-------------------------FICHEROS DE TEXTO-----------------------------//


        for (j=0; j<N; j++)
        {
            fich_funcion << j << ", " << real(phi[j]) << ", " << imag(phi[j]) << ", " << pow(abs(phi[j]),2) << endl;
            fich_posi << j << ", " << posi << endl;
            fich_energias << j << ", " << Epotencial << ", " << Ecinetica << ", " << Etotal << endl;
            fich_momento << j << ", " << momento << endl;
        }
        fich_funcion << endl;
        fich_posi << endl;
        fich_energias << endl;
        fich_momento << endl;




        //--------------------------CÁLCULOS----------------------------//

        
        //Volvemos a generar la función de onda en cada paso

        phi[0]=(0.0,0.0);
        phi[N]=(0.0,0.0);

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

        
        //Guardo P_d de la iteración anterior

        P_d_antes=P_d;



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


        //Calculo la probabilidad de encontrar la partícula a la derecha y la copio en un fichero

        P_d=0.0;

        for(j=round(4.0*n/5.0); j=N; j++)
        {
            P_d=P_d+norma;
        }

        fich_probabilidad << n << ", " << P_d << endl;

        //Compruebo si el valor de P_d obtenido es el máximo

        if(P_d>P_d_antes) 
        {
            n_D=t;
            P_d_max=P_d;

        }


        //Genero un número aleatorio entre 0 y 1 para ver si hemos detectado la partícula

        p=rand()%2;

        if (p>P_d)     m_T=m_T+1;
        else             m_T=m_T;


        //Calculo el coeficiente de transmisión
        K=m_T/t;


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
        Epotencial=0.0;
        Etotal=0.0;


        for(j=0; j<N; j++)
        {
            posi=posi+j*pow(abs(phi[j]),2);
            momento=momento+conj(phi[j])*phi1[j];
            Ecinetica=Ecinetica-conj(phi[j])*phi2[j];
            Epotencial=Epotencial+V[j]*pow(abs(phi[j]),2);
            Etotal=Ecinetica+Epotencial;
        }
        momento=momento*(-1.0)*im_puro;


    }

    fich_funcion.close();
    fich_norma.close();
    fich_K.close();
    fich_Pdmax.close();
    fich_posi.close();
    fich_energias.close();
    fich_momento.close();
    fich_probabilidad.close();




    return 0;
}



