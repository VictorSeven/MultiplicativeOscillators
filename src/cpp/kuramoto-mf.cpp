#include<iostream>
#include<vector>
#include<random>
#include<fstream>
#include<complex>
#include<algorithm>

//Compiler definitions: true or false are shortcuts

#define FALSE 0
#define TRUE 1

//Order of the Kuramoto-Daido parameters
#ifndef ORDER
#define ORDER 10
#endif // ORDER

//Define integration type:
//Series: print temporal series to file
//Diagram: useful to print a phase diagram
#define SERIES 0
#define DIAGRAM 1

//Select current mode
#ifndef MODE
#define MODE SERIES 
#endif // MODE


using namespace std;





//Auxiliary stuff
static const complex<double> I(0.0, 1.0); //Imaginary unit, always useful
mt19937 gen;
uniform_real_distribution<double> ran_u(0.0, 1.0);
normal_distribution<double> ran_g(0.0, 1.0);

//Generate random initial conditions
void initial_conditions(const int N, vector<double> &phi, complex<double> &kuramoto, double &r, double &psi)
{
    int i;
    double x,y;
    phi = vector<double>(N);

    kuramoto = complex<double>(0.0, 0.0);

    x = y = 0.0;
    for (i=0; i < N; i++)
    {
        phi[i] = ran_u(gen) * 6.2838;
        x += cos(phi[i]);
        y += sin(phi[i]);
    }

    //Initial conditions only need the first moment to start.
    kuramoto = complex<double>(x,y) / (1.0 * N);
    r = abs(kuramoto);
    psi = arg(kuramoto);
}

//Step after relaxation
void step(const int N, const double dt, const double sqdt, const double w, const double q, const double s, vector<double> &phi, complex<double> &kuramoto, double &r, double &psi, vector<double> &xy)
{
    int i,j;

    xy = vector<double>(1+2*ORDER, 0.0);
    double auxc, auxs, prevc, prevs, nextc, nexts;

    for (i=0; i < N; i++)
    {
        //Simulation step using MF dynamics
        phi[i] += dt * (w + q * r * sin(psi - phi[i])) + sqdt * ran_g(gen) * s;
        phi[i] = fmod(phi[i], 2*M_PI);

        auxc = cos(phi[i]);
        auxs = sin(phi[i]);

        prevc = auxc;
        prevs = auxs;

        xy[1] += auxc;
        xy[1+ORDER] += auxs;

        //Then compute all requested Kuramoto-Daido parameters
        for (j=2; j <= ORDER; j++)
        {
            //Use trigonometry to obtain sin(k*x) as a function of sin((k-1)x). Same with cosine
            nextc = auxc * prevc - auxs * prevs;
            nexts = auxs * prevc + auxc * prevs;
            xy[j] +=  nextc;
            xy[j+ORDER] += nexts;

            prevc = nextc;
            prevs = nexts;
        }
    }

    //Get the Kuramoto parameter
    kuramoto = complex<double>(xy[1], xy[1+ORDER]) / (1.0 * N);

    //Update r and psi which are needed for MF
    r = abs(kuramoto);
    psi = arg(kuramoto);
}

//Step during relaxation: no computation of KD nor output, to make it fast
void step_relax(const int N, const double dt, const double sqdt, const double w, const double q, const double s, vector<double> &phi, complex<double> &kuramoto, double &r, double &psi)
{
    int i;
    double auxc, auxs;

    auxc = auxs = 0.0;
    for (i=0; i < N; i++)
    {
        //Simulation step using MF dynamics
        phi[i] += dt * (w + q * r * sin(psi - phi[i])) + sqdt * ran_g(gen) * s;
        phi[i] = fmod(phi[i], 2*M_PI);

        auxc += cos(phi[i]);
        auxs += sin(phi[i]);
    }

    //Get the Kuramoto parameter
    kuramoto = complex<double>(auxc, auxs) / (1.0 * N);

    //Update r and psi which are needed for MF
    r = abs(kuramoto);
    psi = arg(kuramoto);
}

//Main code
int main(int argc, char* argv[])
{
    //Main simulation variables
    int N = 5000; //Number of oscillators
    vector<double> phi; //Phase
    double w, s, q; //Parameters
    complex<double> kuramoto;
    double r, psi;  //Kuramoto parameters
    vector<double> xy; //Kuramoto Daido pairs

    //Simulation parameters 
    double tf, trelax;
    double t;
    const double dt = 0.01;
    const double sqdt = sqrt(dt);
        
    //Output
    ofstream output;
    string filename = "prueba";

    #if MODE==DIAGRAM //The program executes several times, to get several values of sigma

        //Counters
        int i;

        //Diagram variables
        double q0,qf,dq;
        int nq;

        //Averages
        int seed;
        complex<double> kd;
        double avr, avr2;
        double avpsi, avpsi2, psi2pi;
        vector< complex<double>  > avkd; //Kuramoto-Daido order parameters
        vector< complex<double>  > avkd2; 
        double nits;

        if (argc == 11)
        {
            N = stoi(argv[1]);
            w = stod(argv[2]);
            s = stod(argv[3]);
            q0 = stod(argv[4]);
            qf = stod(argv[5]);
            nq = stoi(argv[6]);
            dq = (qf - q0) / nq;
            trelax = stod(argv[7]);
            tf = stod(argv[8]);
            filename = string(argv[9]);
            seed = stoi(argv[10]);
            gen.seed(seed);
        }
        else
        {
            cout << "WRONG NUMBER OF PARAMETERS FED" << endl;
            return EXIT_SUCCESS;
        }

        avkd = vector< complex<double> >(ORDER+2, 0.0);

        output.open(filename);
        for (q=q0; q < qf; q += dq)
        {

            //Generate the initial conditions and relax the system
            initial_conditions(N, phi, kuramoto, r, psi);
            for (t=0; t < trelax; t += dt) step_relax(N, dt, sqdt, w, q, s, phi, kuramoto, r, psi);

            //Make measurements of our observables
            avr = avr2 = 0.0;
            avpsi = avpsi2 = 0.0;
            for (t=0.0; t < tf; t += dt)
            {
                step(N, dt, sqdt, w, q, s, phi, kuramoto, r, psi, xy);

                //Average of Kuramoto parameter
                avr += r;
                avr2 += r*r;

                //Average of global phase
                psi2pi = fmod(psi, M_2_PI);
                avpsi += psi2pi; 
                avpsi2 += psi2pi*psi2pi;

                //Get the average of KD parameters
                for (i=2; i <= ORDER; i++)
                {
                    kd = complex<double>(xy[i], xy[i+ORDER]) / (1.0*N);
                    avkd[i] += kd; 
                }
            }

            //Finish averages for Kuramoto
            nits = tf / dt;
            avr  /= nits;
            avr2 /= nits;
            avpsi /= nits;
            avpsi2 /= nits;

            //Write mean and variance to disk
            output << q << " " << avr << " " << avr2 - avr*avr << " " << avpsi2 - avpsi * avpsi << " ";

            //Next do the same with KD parameters
            for (i=2; i <= ORDER; i++) 
            {
                avkd[i] /= nits;
                output << avkd[i].real() << " " << avkd[i].imag() << " ";
            }
            output << endl;
        }
        output.close();


    //Time series mode!
    #elif MODE==SERIES

        int i,j;
        int seed;

        //Get arguments
        if (argc == 9)
        {
            N = stod(argv[1]);
            w = stod(argv[2]);
            s = stod(argv[3]);
            q = stod(argv[4]);
            trelax = stod(argv[5]);
            tf = stod(argv[6]);
            filename = string(argv[7]);
            seed = stoi(argv[8]);
            gen.seed(seed);
        }
        else
        {
            cout << "WRONG NUMBER OF PARAMETERS FED" << endl;
            return EXIT_SUCCESS;
        }

        //Generate the initial conditions and relax the system
        initial_conditions(N, phi, kuramoto, r, psi);
        for (t=0; t < trelax; t += dt) step_relax(N, dt, sqdt, w, q, s, phi, kuramoto, r, psi);


        //Write temporal series while its done!
        output.open(filename);
        i = 0;
        for (t=0.0; t < tf; t += dt)
        {
            step(N, dt, sqdt, w, q, s, phi, kuramoto, r, psi, xy);

            //Write all Kuramoto-Daido from time to time
            //Format: each row has the KD parameters at selected time. Even columns 0,2,4... have the real part, odd the imaginary
            if (i % 100 == 0) 
            {
                output << t << " ";
                for (j=1; j <= ORDER; j++)
                {
                    output << xy[j] / (1.0*N) << " " << xy[j+ORDER] / (1.0*N) << " "; 
                }
                output << endl;
            }
            i++;
        }

        output << t << " ";
        for (j=1; j <= ORDER; j++)
        {
            output << xy[j] / (1.0*N) << " " << xy[j+ORDER] / (1.0*N) << " "; 
        }
        output << endl;
        output.close();

        output.open("phases");
        for (i=0; i<N-1;i++) output << phi[i] << " ";
        output << phi[N-1] << endl;
        output.close();

    #endif

    return 0; 
}


