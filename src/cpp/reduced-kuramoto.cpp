#include<iostream>
#include<vector>
#include<random>
#include<fstream>
#include<complex>
#include<algorithm>

#define DIAGRAM 0
#define TIMESERIES 1

#ifndef MODE
#define MODE TIMESERIES 
#endif

using namespace std;

//RNG
mt19937 gen(563153215);
normal_distribution<double> ran_g(0,1);
uniform_real_distribution<double> ran_u(0,1);

//Generation of multivariable normal distribution
void cholesky(const int nharm, const vector<vector<double>> &corr, vector<vector<double>> &b);
void get_gaussian_vars(const int nharm, const vector<vector<double>> &b, vector<double> &rg);
void get_corrs(const int nharm, const double nsize, const vector<double> &r, vector<vector<double>> &corr);

void initial_conditions(const int nharm, vector<double> &r);

void step_tl(const int nharm, const vector<double> &old_r, vector<double> &r, 
vector<vector<double>> &corr, vector<vector<double>> &b, 
const double q, const double s2, const double nsize, const double dt, const double sqdt, bool flag);

void step_tl_measures(const int nharm, const vector<double> &old_r, vector<double> &r, 
vector<vector<double>> &corr, vector<vector<double>> &b, 
const double q, const double s2, const double nsize, const double dt, const double sqdt,
vector<double> &avr, vector<double> &avr2);

void step_fs(const int nharm, const vector<double> &old_r, vector<double> &r, 
const double q, const double s2, const double nsize,
const double dt, const double sqdt);

int main(int argc, char* argv[])
{
    int i,j,k;
    const int nharm = 10;
    vector<double> r(nharm);
    vector<double> old_r(nharm);

    vector<vector<double>> b = vector<vector<double>>(nharm, vector<double>(nharm));
    vector<vector<double>> corr = vector<vector<double>>(nharm, vector<double>(nharm));

    ofstream output;
    string filepath;

    double t;
    const double tf = 400.0;
    const double t_thermal = 0.0;
    double dt = 0.01;
    double sqdt = sqrt(dt);
    int nits = tf/dt;

    #if MODE==DIAGRAM
        vector<double> avr(nharm, 0.0);
        vector<double> avr2(nharm, 0.0);

        double q;
        double q0, qf, dq;
        int nq;

        double sys_size;
        double s2;
        double nsize;

        if (argc==7)
        {
            q0 = stod(argv[1]);
            qf = stod(argv[2]);
            nq = stoi(argv[3]);
            dq = (qf - q0) / nq;

            sys_size = stod(argv[4]);
            s2 = stod(argv[5]);
            s2 *= s2;
            nsize = 0.5*s2/sys_size;

            filepath = string(argv[6]);
        }


        output.open(filepath);
        for (q=q0; q < qf; q += dq)
        {
            initial_conditions(nharm, old_r);
            r[0] = 1.0; //To avoid problems with swapping later
            avr[0] = 1.0;
            avr2[0] = 1.0;

            for (t=0; t<t_thermal; t += dt)
            {
                step_tl(nharm, old_r, r, corr, b, q, s2, nsize, dt, sqdt);
                //if (t < 2*dt) cout << old_r[0] << " " << old_r[5] << endl;
                swap(old_r, r);
                //if (t < 2*dt) cout << old_r[0] << " " << old_r[5] << endl;
            }

            avr = vector<double>(nharm, 0.0);
            avr2 = vector<double>(nharm, 0.0);
            for (t=0; t<tf; t += dt)
            {
                step_tl_measures(nharm, old_r, r, corr, b, q, s2, nsize, dt, sqdt, avr, avr2);
                swap(old_r, r);
            }
            
            output << q << " "; 
            for (k=0; k < nharm-1; k++) output << avr[k]/(1.0*nits) << " " << avr2[k]/(1.0*nits) - avr[k]*avr[k]/(1.0*nits*nits) << " ";
            output << endl;
        }
        output.close();

    #elif MODE==TIMESERIES
        double q;
        double sys_size;
        double s2;
        double nsize;

        if (argc==5)
        {
            q = stod(argv[1]);
            sys_size = stod(argv[2]);
            s2 = stod(argv[3]);
            s2 *= s2;
            nsize = 0.5*s2/sys_size;

            filepath = string(argv[4]);
        }


        initial_conditions(nharm, old_r);
        r[0] = 1.0; //To avoid problems with swapping later

        get_corrs(nharm, nsize, old_r, corr);
    for (i=0; i < nharm; i++)
    {
        for (j=0; j <= i; j++)
        {
            cout << corr[i][j] << " "; 
        }
        cout << endl;
    }
    cout << endl << "r:"<< endl; 
    for (i=0; i < nharm; i++)
    {
        cout << old_r[i] << " "; 
    }
    cout << endl;

        

        output.open(filepath);
        for (t=0; t<t_thermal; t += dt)
        {
            step_tl(nharm, old_r, r, corr, b, q, s2, nsize, dt, sqdt, false);
            output << t << " ";
            for (k=0; k < nharm; k++) output << old_r[k] << " ";
            output << endl;
            swap(old_r, r);
        }
        //dt = 1e-4;
        //sqdt = sqrt(dt);
        for (t=0; t<tf; t += dt)
        {
            step_tl(nharm, old_r, r, corr, b, q, s2, nsize, dt, sqdt, true);

            output << t+t_thermal << " ";
            for (k=0; k < nharm; k++) output << old_r[k] << " ";
            output << endl;

            swap(old_r, r);
        }
        output.close();
        

    #endif

    return EXIT_SUCCESS;
}

void initial_conditions(const int nharm, vector<double> &r)
{
    int k;

    r[0] = 1.0; //Always this!
    r[1] = 0.01 * ran_u(gen);
    for (k=2; k < nharm; k++)
    {
        r[k] = pow(r[1], k); 
    }
}

void step_tl(const int nharm, const vector<double> &old_r, vector<double> &r, 
vector<vector<double>> &corr, vector<vector<double>> &b, 
const double q, const double s2, const double nsize, const double dt, const double sqdt, bool flag)
{
    int k, j;
    double det, sto;

    vector<double> z(nharm);

    //Get the matrix b from the correlation matrix
    get_corrs(nharm, nsize, old_r, corr);
    /*for (int i=0; i < nharm; i++)
    {
        for (j=0; j <= i; j++)
        {
            cout << corr[i][j] << " "; 
        }
        cout << endl;
    }*/
    //cout << endl << "b:"<< endl; 
    cholesky(nharm, corr, b);
    /*for (int i=0; i < nharm; i++)
    {
        for (j=0; j <= i; j++)
        {
            cout << b[i][j] << " "; 
        }
        cout << endl;
    }*/
    //cout << endl << "r:"<< endl; 


    //Index k=0 corresponds to r[0]=1 always
    for (k=1; k < nharm-1; k++)
    {
        //Deterministic part
        det = 0.5*k* (q * old_r[1] * (old_r[k-1] - old_r[k+1]) - k*s2*old_r[k]); 
        //cout << old_r[k-1] << " " << det << endl;

        //Stochastic part
        z[k] = ran_g(gen);
        sto = 0.0;
        for (j=1; j <= k; j++)
        {
            sto += b[k][j] * z[j]; 
        }

        //Update harmonic
        if (flag) r[k] = old_r[k] + dt * det + sqdt * sto;
        else r[k] = old_r[k] + dt * det;
        //r[k] = old_r[k] + dt * det;
        r[k] = r[k] < 0 ? 1e-50 : r[k];
    }


    //For last one, make r[k+1]=0
    //Note that here k=nharm-1 after the loop
    det = 0.5*k* (q * old_r[1] * old_r[k-1] - k*s2*old_r[k]); 

    //Stochastic part
    z[k] = ran_g(gen);
    sto = 0.0;
    for (j=1; j <= k; j++)
    {
        sto += b[k][j] * z[j]; 
    }

    //Update last
    if (flag) r[k] = old_r[k] + dt * det + sqdt * sto;
    else r[k] = old_r[k] + dt * det;
    //r[k] = old_r[k] + dt * det;
    r[k] = r[k] < 0 ? 1e-50 : r[k];

    /*
    for (int i=0; i < nharm; i++)
    {
        cout << r[i] << " "; 
    }
    cout << endl;*/

    return;
}

void step_tl_measures(const int nharm, const vector<double> &old_r, vector<double> &r, 
vector<vector<double>> &corr, vector<vector<double>> &b, 
const double q, const double s2, const double nsize, const double dt, const double sqdt,
vector<double> &avr, vector<double> &avr2)
{
    int k, j;
    double det, sto;

    vector<double> z(nharm);

    //Get the matrix b from the correlation matrix
    get_corrs(nharm, nsize, old_r, corr);
    cholesky(nharm, corr, b);

    //Index k=0 corresponds to r[0]=1 always
    for (k=1; k < nharm-1; k++)
    {
        //Deterministic part
        det = (0.5*k* (q * old_r[1] * (old_r[k-1] - old_r[k+1]) - k*s2*old_r[k])); 
        //cout << k << " " << det << " " << old_r[k] << " " << old_r[k-1] << " " << old_r[k+1] << endl;

        //Stochastic part
        z[k] = ran_g(gen);
        sto = 0.0;
        for (j=1; j <= k; j++)
        {
            sto += b[k][j] * z[j]; 
        }

        //Update harmonic
        r[k] = old_r[k] + dt * det + sqdt * sto;
        r[k] = max(r[k], 0.0);

        avr[k] += r[k];
        avr2[k] += r[k]*r[k];
    }


    //For last one, make r[k+1]=0
    //Note that here k=nharm-1 after the loop
    det = (0.5*k* (q * old_r[1] * old_r[k-1] - k*s2*old_r[k])); 

    //Stochastic part
    z[k] = ran_g(gen);
    sto = 0.0;
    for (j=1; j <= k; j++)
    {
        sto += b[k][j] * z[j]; 
    }

    //Update last
    r[k] = old_r[k] + dt * det + sqdt * sto;
    r[k] = max(r[k], 0.0);
    avr[k] += r[k];
    avr2[k] += r[k]*r[k];

    return;
}

void step_fs(const int nharm, const vector<double> &old_r, vector<double> &r, 
const double q, const double s2, const double nsize,
const double dt, const double sqdt)
{
    int k;

    //Index k=0 corresponds to r[0]=1 always
    for (k=1; k < nharm-1; k++)
    {
        r[k] += dt * (0.5*k* (q * old_r[1] * (old_r[k-1] - old_r[k+1]) - k*s2*old_r[k]) + nsize * (1+old_r[k]*old_r[k])/old_r[k]); 
    }
    //For last one, make r[k+1]=0
    //Note that here k=nharm-1 after the loop
    r[k] += dt * (0.5*k* (q * old_r[1] * old_r[k-1] - k*s2*old_r[k]) + nsize * (1+old_r[k]*old_r[k])/old_r[k]); 

    return;
}

void get_corrs(const int nharm, const double nsize, const vector<double> &r, vector<vector<double>> &corr)
{
    int i,j;
    double ipj;

    corr[0][0] = 10000; 

    for (i=1; i < nharm; i++)
    {
        corr[0][i] = 0.0;
        corr[i][0] = 0.0; 
        for (j=1; j <= i; j++)
        {
            ipj = (i+j < nharm) ? r[i+j] : 0.0; 
            //cout << i << " " << j << " "  << i+j << " " << nharm << " " << ipj << endl;
            corr[i][j] = i*j*nsize*(r[i-j] - ipj);
            corr[j][i] = corr[i][j]; 
        }
    }

    return;
}

//Simple iterative Cholesky-like computation of square root of correlation matrix
void cholesky(const int nharm, const vector<vector<double>> &corr, vector<vector<double>> &b)
{
    double sum;
    int i,j,k;
    for (i = 0; i < nharm; i++) {
        for (j = 0; j <= i; j++) {
            sum = 0.;
            for (k = 0; k < j; k++)
                sum += b[i][k] * b[j][k];

            if (i == j)
                b[i][j] = sqrt(corr[i][i] - sum);
            else
                b[i][j] = (1.0 / b[j][j] * (corr[i][j] - sum));
        }
    }
}

//Sample Gaussian variables from matrix b
void get_gaussian_vars(const int nharm, const vector<vector<double>> &b, vector<double> &rg)
{
    int i,j;
    vector<double> z(nharm);

    rg = vector<double>(nharm, 0.0);

    for (i=0; i < nharm; i++)
    {
        //Generate a new random variable and get the sums
        z[i] = ran_g(gen);
        for (j=0; j <= i; j++)
        {
            rg[i] += b[i][j] * z[j]; 
        }
    }
}








