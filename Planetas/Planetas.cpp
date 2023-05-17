#include <iostream>
#include <fstream>
#include <cmath>

#define G 6.67e-20  //km³/(kgs)
#define c 1.5e8     //km
#define Ms 1.99e30  //kg


using namespace std;



//Declaro las funciones del programa
void copiar_fichero (double vector[9], string fichero);
void aceleracion (double m[9],double r1[9], double r2[9], double a[9]);
void posicion (double h, double r[9], double v[9], double a[9]);
void velocidad (double h, double w[9], double a[9], double v[9]);
void funcion_w (double h, double v[9], double a[9], double w[9]);
double modulo (double r1, double r2);
double energia (double r1[9], double r2[9], double m[9], double vx[9], double vy[9]);
double energia (double r1[9], double r2[9], double m[9], double vx[9], double vy[9]);

//----------------------------------------------FUNCIÓN PRINCIPAL----------------------------------------------

int main ()
{
    //Declaro los contadores para los bucles y el real para el número de iteraciones
    int j=0;
    int k=0;
    double i=0;
    double contador[9]={0};

    //Declaro e inicio los valores del tiempo y del paso
    double h;
    double t;
    t=0;
    double tmax;


    //Declaro e inicializo los vectores a 0

    double x[9];
    double y[9]={0};      //Obligo a que los planetas se situen alineados sobre el eje x y empecen a desplazarse desde ahí 
    double vx[9]={0};     //Obligo a que los planetas empiecen desplazandose en la direccion y
    double vy[9];
    double ax[9]={0};
    double ay[9]={0};
    double m[9];
    double wx[9]={0};
    double wy[9]={0};
    double E=0;
    double y_ant[9]={0};
 

    //Declaro ficheros donde guardar los datos de las posiciones, de la energía y del momento angular
    ofstream fich_posi;
    ofstream fich_energia;
    ofstream fich_periodo;


    //Pido al usuario que establezca un tiempo máximo para limitar las iteraciones, además del paso entre iteraciones

    cout << "Establezca un tiempo máximo" << endl;
    cin >> tmax;

    cout << "Establezca un paso" << endl;
    cin >> h;

    //Para calcular el número de iteraciones, divido el tiempo máximo entre h

    i= (int)tmax/h;



//-------------------------------VALORES INICIALES---------------------------


    //Copio los valores iniciales de posiciones y velocidades en un vector
    copiar_fichero (m, "masas.txt");
    copiar_fichero (x, "posiciones_iniciales.txt");
    copiar_fichero (vy, "velocidades_iniciales.txt");


    //Reescalo las velocidades
    for (j=0; j<9; j++)
    {
        vy[j]=vy[j]/sqrt(G*Ms/c);
    }


    //Calculo la aceleración inicial
    aceleracion (m, x, y, ax);
    aceleracion (m, y, x, ay);


    


//-------------------------------ALGORITMO DE VERLET------------------------------

    //Mediante una iteración, aplico el algoritmo de Verlet mientras no se alcance el tiempo máximo



    //Hago que en cada iteración escriba la posición en x e y de todos los planetas, y la energía en función del tiempo

    fich_posi.open("planets_data.dat");
    fich_energia.open("energia.txt");
    fich_periodo.open("periodo.txt");

    for (j=0; j<=i; j++)
    {
        E=0;

        //Calculo la posición de los planetas en el eje y en la iteración anterior para calcular el periodo
        for (k=0; k<9; k++)
        {
            y_ant[k]=y[k];
        }

        

        //Calculo la energía
        E=energia (x, y, m, vx, vy);


        //Completo el fichero de la energía en función del tiempo
        fich_energia << t << "  " << E << endl;



        //Completo el fichero de las posiciones

        fich_posi << x[0] << "," << y[0] << endl;
        fich_posi << x[1] << "," << y[1] << endl;
        fich_posi << x[2] << "," << y[2] << endl;
        fich_posi << x[3] << "," << y[3] << endl;
        fich_posi << x[4] << "," << y[4] << endl;
        fich_posi << x[5] << "," << y[5] << endl;
        fich_posi << x[6] << "," << y[6] << endl;
        fich_posi << x[7] << "," << y[7] << endl;
        fich_posi << x[8] << "," << y[8] << endl;
        fich_posi << endl;

        

        //Calculo la posicion
        posicion (h, x, vx, ax);
        posicion (h, y, vy, ay);


        //Calculo la funcion w
        funcion_w (h, vx, ax, wx);
        funcion_w (h, vy, ay, wy);

        for (k=0;k<9; k++)
        {
            ax[k]=0;
            ay[k]=0;
        }
        
        //Calculo la aceleración 
        aceleracion (m, x, y, ax);
        aceleracion (m, y, x, ay);



        //Calculo la velocidad
        velocidad (h, wx, ax, vx);
        velocidad (h, wy, ay, vy);


        //Calculo el periodo y lo paso a un fichero
        for (k=1; k<9; k++)
        {
            
            if (contador[k]==0)
            {
                if ((y_ant[k]<0) && y[k]>0)
                {
                    contador[k]=1;
                    fich_periodo << k << " , " << t << " , " << t*58.1 << endl;

                }
            }
        }




        t=t+h;


    }

    fich_posi.close();
    fich_energia.close();
    fich_periodo.close();



    

    return 0;

}




//--------------------------------------------FUNCIONES DEL PROGRAMA--------------------------------------------

//Función para pasar los datos de un fichero de texto a un vector

void copiar_fichero (double vector[9], string fichero)
{
    int i;
    ifstream fich;

    fich.open (fichero);

   for (i=0; i<9; i++)
   {
    fich >> vector[i];
   }

    fich.close();

    return;
}


//Función para obtener el módulo de la posición, facilitando el cálculo de la aceleración

double modulo (double r1, double r2)
{
    double modulo;
    modulo=0.0;

    modulo=sqrt(r1*r1+r2*r2);
    
    return modulo;
    
}


//Función para obtener la aceleración

void aceleracion (double m[9],double r1[9], double r2[9], double a[9])
{
    //Declaro dos enteros para hacer dos bucles 
    int i,j;

    //Hago un bucle que recorra todo el vector y otro bucle en su interior que opere con todos los elementos del vector de índice distinto

    for (i=0; i<9; i++)
    {
        for (j=0; j<9; j++)
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

void posicion (double h,double r[9], double v[9], double a[9])
{
    int i;
    
    for (i=0; i<9; i++)
    {
        r[i]=r[i]+h*v[i]+pow(h,2)*a[i]/2.0;
    }

    return;
}


//Función para calcular la velocidad en tiempo h

void velocidad (double h, double w[9], double a[9], double v[9])
{
    int i;

    for (i=0; i<9; i++)
    {
        v[i]=w[i]+h*a[i]/2.0;
    }

    return;
}




//Función para definir la función w

void funcion_w (double h, double v[9], double a[9], double w[9])
{
    int i;

    for (i=0; i<9; i++)
    {
        w[i]=v[i]+h*a[i]/2.0;
    
    }

    return;
}



double energia (double r1[9], double r2[9], double m[9], double vx[9], double vy[9])
{
    int i, j;
    double E;

    E=0;

        //En este caso no tendremos en cuenta el Sol
        for (i=1; i<9; i++)

        {
            //Término de la energía cinética. Cada planeta tiene la suya propia
            E=E+0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]);

            //Término de la energía potencial gravitatoria. Cada planeta con todos los demás. Exijo que no la calcule sobre sí mismo
            for (j=0; j<9; j++)
             {
                if(i!=j)
                {
                    E=E-m[i]*m[j]/modulo(r1[i]-r1[j],r2[i]-r2[j]);
                }

             }
        }
    
        
    return E;  
}


