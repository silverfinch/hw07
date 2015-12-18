#include <stdio.h>
#include <math.h>
#include <string.h>
#include "adjust.c"
#include "timer.c"
#include <gsl/gsl_integration.h>

double fun (double x, void *params);
void integral_recur (int nmin, int nmax, double vals[]);
void integral_gen (int nmin, int nmax, double vals[]);

void integral_recur (int nmin, int nmax, double vals[])
{
    vals[nmax] = 0.00367843;
    for (int n = nmax - 1; n >= nmin; n--)      //is nmax-1 because vals[nmax] has already been provided
    {
        vals[n] = (vals[n + 1] + 1 / M_E) / (n + 1);
    }
}

void integral_gen (int nmin, int nmax, double *vals)
{
    double a = 0., b = 1.;      // limits of integration
    double abserr = 1.e-10, relerr = 1.e-10; // requested errors
    double result;              // the integral value
    double error;
    size_t np = 1000;           // work area size
    int n;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc (np);

    gsl_function F;

    F.function = &fun;
    F.params = &n;

    for (n = nmin; n <= nmax; n++)
    {
        gsl_integration_qag (&F, a, b, abserr, relerr, np,
                 GSL_INTEG_GAUSS15, w, &result, &error);
        vals[n] = result;
    }

    gsl_integration_workspace_free (w);
}

double fun (double x, void *params)
{
    double n = (double) *(int *) params;//converts to a double, but treats the point initially as an integer since n is defined as an integer
    double f = exp (-x) * pow (x, n);//this is the function to integrate

    return f;
}

#define NMAX 100 //100 terms

int main (void)
{
    int nmax = NMAX, nmin = 0;
    int n;
    double vals1[NMAX + 1], vals2[NMAX + 1];

    integral_gen(nmin, nmax, vals2);
    //performs the integrals with each method
    integral_recur (nmin, nmax, vals1);
    printf ("Index   Abs diff b/t recurrence and general algorithms\n");
    for (n = nmax; n >= nmin; n--)
    {
	printf ("%5d   %e\n", n, fabs(vals1[n] - vals2[n]));//abs difference in results between recurrence and general integrators for each index
    }
    timer_start();
    integral_recur (nmin, nmax, vals1);
    double t_recur = timer_stop();//measuring the time of one recurrence integrator call
    timer_start();
    integral_gen (nmin, nmax, vals2);
    double t_gen = timer_stop();//measuring the time of one general integrator call
    int count_recur = adjust_rep_count(1, t_recur, 1., 2.);
    int count_gen = adjust_rep_count(5, t_gen, 1., 2.);//like in midterm 1, using adjust_rep_count to get a number of total calls that will take 1-2 seconds all together
    for (int k = 1; k <= count_recur; k++)//running each function for the given number of calls, and dividing by the total time it takes to get the ratios printed at the bottom.
    {
        timer_start();
        integral_recur(nmin, nmax, vals1);
        t_recur += timer_stop();
    }
    for (int k = 1; k <= count_gen; k++)
    {
        timer_start();
        integral_gen(nmin, nmax, vals2);
        t_gen += timer_stop();
    }
    double ratio_recur = t_recur / count_recur;
    double ratio_gen = t_gen / count_gen;
    printf("\nGen call time|Recur call time|General/recurrence");
    printf("\n%e   %e   %e\n", ratio_gen, ratio_recur, ratio_gen / ratio_recur);
}
