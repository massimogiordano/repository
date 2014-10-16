#include <iostream>
#include <solarsystem.h>
#include <planet.h>
#include <cmath>
#include <armadillo>

#define RUDDAKUTA 1;
using namespace arma;
using namespace std;



int main()
{
    int RK4 = 0; //set zero for not use this method
    solarsystem mysystem;

    //add a planet    Mass, x,y,z, vx,vy,vz;
   planet sun(   1.98e30, 0, 0,     0,0,0,   0 );
   planet earth( 5.24e24, 1.5e11, 0 ,0,0,30000,0);

    //au
  //  planet sun(   1, 0, 0, 0,0,0, 0 );
  // planet earth( 1, 1, 0 ,0,0,0.1,0);

    planet marth( 6.41e23, 2.2e11 ,0,0,0,24000,0);
    planet giove(1.89e27,7.78e11, 0,0,0,13000,0);
    planet saturno(5.7e26,1.4e12,0,0,0,9000,0);
    planet k( 6.41e23, 0, 2.2e1,0,0,24000,2);
    planet satellite(5.24e26, -1.5e11,0,0,0,30000,0);
    planet luna(7.3e2, 1.5e11, 3.84e5,0,-1000, 30000,0);
    planet prova(1.98e30, 0, 0 ,0,0,-10,0);
    planet prova2(1.98e30, 1e11, 0 ,0,0,30000,0);


    mysystem.add(sun);
    mysystem.add(earth);
   // mysystem.add(marth);
    //mysystem.add(giove);
    //mysystem.add(saturno);
    //mysystem.add(satellite);
    //mysystem.add(luna);


    //mysystem.add(prova);
    //mysystem.add(prova2);

    mysystem.print_position(mysystem.all_planets,3);

    int elements = mysystem.number_planets;
    cout << "number of element" << elements<< endl;
    //vettore con tutti i pianeti, h , t_max



    if(RK4) {
    mysystem.solverRK4(mysystem.all_planets, 500, 6e7 );
    cout << "rudda";
   }else{
    mysystem.solverVERLET(mysystem.all_planets, 5000, 9e7 );
    cout << "vertel";
    }

}

