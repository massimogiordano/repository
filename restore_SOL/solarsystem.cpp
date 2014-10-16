#include "solarsystem.h"
#include "planet.h"
#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;


solarsystem::solarsystem()
{
}
void solarsystem::add(planet n){
     number_planets++;
     all_planets.push_back(n);  //questa funzione inserisce un elementp nella push_bac
 }
void solarsystem::print_position(vector<planet> vec){
    print_position(vec, 3);
}

void solarsystem::print_position(vector<planet> vec, int n){
    if(n>3 || n<=0) n=3;
    for(int i=0; i<vec.size(); i++){
        planet &questo = vec[i];
        std::cout << std::scientific;
        for(int j=0; j<n;j++){
        std::cout << questo.position[j] << "   ";
        }
        std::cout << "       ";
       }
    std::cout << std::endl;
}

void solarsystem::insert_data(vector<planet> vec, mat &ma){
    for(int i=0; i<vec.size(); i++){
        planet &questo = vec[i];
        ma(i,6)=questo.mass;

        for(int k=0; k<3;k++){
            ma(i,k)=questo.position[k];
            ma(i,k+3)=questo.velocity[k];
        }
    }
}

void solarsystem::synctroniz(vector<planet> vec, mat &ma){
    int n = vec.size();

    for(int j=0; j<n;j++){
        planet &questo = vec[j];
    for (int i = 0; i < 3; ++i){
       questo.position[i] =  ma(j,i);
       questo.velocity[i] = ma(j,i+3);
}
}

}

void solarsystem::solverRK4(vector<planet> vec, double h, double tmax){

    mat y_i(number_planets,7);
    mat y_i_temp(number_planets,7);
    mat k1(number_planets,7);
    mat k2(number_planets,7);
    mat k3(number_planets,7);
    mat k4(number_planets,7);

    insert_data(vec , y_i);
    double t=0;

    while(t<tmax){

        derivate(y_i, k1, number_planets);

        sum_matrix(y_i_temp, 1, y_i, 0.5*h, k1, number_planets);
        derivate(y_i_temp, k2, number_planets);

        sum_matrix( y_i_temp, 1,  y_i, 0.5*h,  k2, number_planets);

        derivate( y_i_temp,  k3, number_planets);

        sum_matrix( y_i_temp, 1,  y_i, h,  k3, number_planets);

        derivate( y_i_temp,  k4, number_planets);

        for(int j=0; j<number_planets; j++){

             for(int i=0; i<6; i++){
                 y_i(j,i) = y_i(j,i) + h*(k1(j,i) + 2*k2(j,i) + 2*k3(j,i) + k4(j,i))/6;
             }

             //aggiorno posizione  e velocità.
             planet &questo = vec[j];
             for(int i=0; i<3; i++){
             questo.position[i] = y_i(j,i);
             questo.velocity[i] = y_i(j,i+3);
             }

        }


print_position(vec,3);


 t+=h;
}



    }

void solarsystem::solverVERLET(vector<planet> vec, double h, double tmax){

    mat y_i(number_planets,7);
    mat r_i_dt(number_planets,7);
    mat a_dt(number_planets,7);
    mat v_dt(number_planets,7);
    mat v_dt_2(number_planets,7);
    //mat a(number_planets,7);

    insert_data(vec , y_i);

    double t=0;

    while(t<tmax){

        derivate(y_i, a_dt, number_planets);

        for(int j=0; j<number_planets; j++){

             for(int i=0; i<3; i++){

                 y_i(j,i) = y_i(j,i) + h*y_i(j,i+3) + 0.5*h*h*a_dt(j,i+3);

                 v_dt_2(j,i+3) = y_i(j,i+3) + 0.5*h*a_dt(j,i+3);
             }

        }

        derivate(y_i, a_dt, number_planets);

        for(int j=0; j<number_planets; j++){

             for(int i=3; i<6; i++){

                 y_i(j,i) = v_dt_2(j,i) + 0.5*h*a_dt(j,i);

             }


             //aggiorno posizione  e velocità.
             planet &questo = vec[j];
             for(int i=0; i<3; i++){
             questo.position[i] = y_i(j,i);
             questo.velocity[i] = y_i(j,i+3);
             }
        }

//cout << sqrt(y_i(1,0)*y_i(1,0) + y_i(1,1)*y_i(1,1)) << "  ";
print_position(vec,2);


 t+=h;
}



}


//______________________ FUNZIONI PRESE DAL MAIN ______________________________
void solarsystem::sum_matrix(mat &result, int coeff_one, mat &first,int coeff_two, mat &second, int n){
    for(int i=0; i<7; i++){
         for(int j=0; j<n; j++){
            result(j,i) = coeff_one*first(j,i) + coeff_two*second(j,i);
         }
    }
}
void solarsystem::printmat(mat &ma, int n){
    cout << endl;
    for(int i=0; i<7; i++){

        for(int k=0; k<n;k++){
            cout <<  ma(k,i)<<" " ;
        } cout << endl;}
}
double solarsystem::force(double x, double y, double z, double Mothers){
    double G=6.67e-11;
    double force=0;
    double distance=0;

    distance = x*x + y*y + z*z;

    force = G*Mothers/pow(distance, 1.5);

    return force;
}
void solarsystem::derivate(mat &dat, mat &de, int n){

    double accelleration_x=0,accelleration_y=0,accelleration_z=0, mod_force;
    for(int i=0; i<n; i++){

        accelleration_x=0,accelleration_y=0,accelleration_z=0;
        for(int j=0; j<n; j++){
            if(i!=j){

mod_force = force(dat(j,0)-dat(i,0),dat(j,1)-dat(i,1) ,dat(j,2)-dat(i,2),  dat(j,6));


accelleration_x += mod_force*(dat(j,0)-dat(i,0));
accelleration_y += mod_force*(dat(j,1)-dat(i,1));
accelleration_z += mod_force*(dat(j,2)-dat(i,2));
}
}
        de(i,3) = accelleration_x; //velx
        de(i,4) = accelleration_y; //vely
        de(i,5) = accelleration_z; //velz

    }


    for(int i=0; i<n; i++){
        de(i,0) = dat(i,3); //velx
        de(i,1) = dat(i,4); //vely
        de(i,2) = dat(i,5); //velz
    }
}










