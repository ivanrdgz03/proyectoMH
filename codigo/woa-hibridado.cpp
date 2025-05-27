extern "C"
{
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <cstring>
#include "random.hpp"
#include <numeric>
using namespace std;

using Random = effolkronium::random_static;

//Funciones auxiliares
vector<double> abs_vector(const vector<double> &vec)
{
    vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = abs(vec[i]);
    }
    return result;
}

vector<double> multiply_scalar_vector(const double &scalar, const vector<double> &vec)
{
    vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = scalar * vec[i];
    }
    return result;
}

vector<double> add_vectors(const vector<double> &vec1, const vector<double> &vec2)
{
    vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i)
    {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

vector<double> subtract_vectors(const vector<double> &vec1, const vector<double> &vec2)
{
    vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i)
    {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

double generate_random_0_1()
{
    return Random::get<double>(0.0, 1.0);
}

int generate_random_int(int min, int max)
{
    return Random::get<int>(min, max);
}

//
double WOA_para_CEC(const int &func_id, const double &mu, const int& n_ballenas, const int &phi, const int& M,  const int &dim, const vector<double> &lim_inferior, const vector<double> &lim_superior, const double &b, const unsigned int &maxevals)
{
    unsigned int n_evals = 0;

    vector<vector<double>> poblacion_ballenas(n_ballenas, vector<double>(dim));

    for (size_t i = 0; i < n_ballenas; ++i)
    {
        for (size_t j = 0; j < dim; ++j)
        {
            poblacion_ballenas[i][j] = lim_inferior[j] + (lim_superior[j] - lim_inferior[j]) * generate_random_0_1();
        }
    }

    double mejor_aptitud_global = numeric_limits<double>::max();
    vector<double> mejor_posicion_global(dim);

    for (size_t i = 0; i < n_ballenas; ++i)
    {
        double aptitud_actual = cec17_fitness(&poblacion_ballenas[i][0]);
        n_evals++;

        if (aptitud_actual < mejor_aptitud_global)
        {
            mejor_aptitud_global = aptitud_actual;
            mejor_posicion_global = poblacion_ballenas[i];
        }
    }

    //Parámetros para el enfriamiento simulado
    //double mu = 0.2;
    //double phi = 0.3;
    double T0;
    double Tf = 1e-6;
    //int M = maxevals/n_ballenas;

    if (mejor_aptitud_global <= 0) { // Evitar logaritmo de cero o negativo
        T0 = 1.0; // Valor por defecto
    } else {
        T0 = (mu * mejor_aptitud_global) / (-log(phi));
    }
    double current_T = T0;

    double beta = (T0 - Tf) / (static_cast<double>(M) * T0 * Tf);

    int current_iteration = 0;
    while (n_evals < maxevals)
    {
        if (current_iteration > 0) { // No enfriar en la primera iteración
             current_T = current_T / (1.0 + beta * current_T);
        }

        if (current_T < Tf) {
            current_T = Tf;
        }

        // a es proporcional a la temperatura normalizada entre T0 y Tf.
        double a = 2.0 * (current_T / T0);
        if (a < 0.0) a = 0.0;
        if (a > 2.0) a = 2.0;


        for (int i = 0; i < n_ballenas; ++i)
        {
            if (n_evals >= maxevals)
            {
                break;
            }

            double r1 = generate_random_0_1();
            double r2 = generate_random_0_1();
            double A = 2.0 * a * r1 - a;
            double C = 2.0 * r2;
            //double b = 0.1;
            double l = 2.0 * generate_random_0_1() - 1.0;

            double p = generate_random_0_1();

            vector<double> nueva_posicion(dim);

            for (int j = 0; j < dim; ++j)
            {
                if (p < 0.5)
                {
                    if (abs(A) >= 1.0)
                    {
                        int randomLeaderIndex = generate_random_int(0, n_ballenas - 1);
                        const vector<double> &randomPosition = poblacion_ballenas[randomLeaderIndex];

                        double D_X_rand = abs(C * randomPosition[j] - poblacion_ballenas[i][j]);
                        nueva_posicion[j] = randomPosition[j] - A * D_X_rand;
                    }
                    else
                    {
                        double D_Leader = abs(C * mejor_posicion_global[j] - poblacion_ballenas[i][j]);
                        nueva_posicion[j] = mejor_posicion_global[j] - A * D_Leader;
                    }
                }
                else
                {
                    double distanceToLeader = abs(mejor_posicion_global[j] - poblacion_ballenas[i][j]);
                    nueva_posicion[j] = distanceToLeader * exp(b * l) * cos(l * 2.0 * M_PI) + mejor_posicion_global[j];
                }
            }

            // Aplicar límites volviendo a un valor dentro de los límites si es necesario
            for (int z = 0; z < dim; ++z)
            {
                if (nueva_posicion[z] < lim_inferior[z] || nueva_posicion[z] > lim_superior[z])
                    nueva_posicion[z] = lim_inferior[z] + generate_random_0_1() * (lim_superior[z] - lim_inferior[z]);
            }

            poblacion_ballenas[i] = nueva_posicion;

            double aptitud_actual = cec17_fitness(&poblacion_ballenas[i][0]);
            n_evals++;

            if (aptitud_actual < mejor_aptitud_global)
            {
                mejor_aptitud_global = aptitud_actual;
                mejor_posicion_global = poblacion_ballenas[i];
            }
        }

        current_iteration++;
        if (n_evals >= maxevals)
        {
            break;
        }
    }
    return mejor_aptitud_global;
}

int main(int argc, char *argv[])
{
    if(argc != 2){
        cerr<<"Uso: " << argv[0] << " <carpeta_resultados>" << endl;
        exit(1);
    }
    int dim = 30; // Dimensión del problema
    double mu = 4.65;
    int num_ballenas = 50;
    double phi = 0.05;
    int M = 1100;
    double b = 0.13;

    unsigned int max_evals = 10000 * dim;
    //Random::seed(3);
    vector<double> lower_bound(dim, -100.0);
    vector<double> upper_bound(dim, 100.0);
    vector<double> resultados;

    for (int func_id = 1; func_id <= 30; ++func_id)
    {
        cec17_init(argv[1], func_id, dim);
        double best_fitness_found = WOA_para_CEC(func_id, mu, num_ballenas, phi, M, dim, lower_bound, upper_bound, b, max_evals);
        
        cout << "WOA [F" << func_id << "]: " << scientific << cec17_error(best_fitness_found) << endl;
    }

    return 0;
}