#include <stdio.h>
#include <math.h>

#include "ising_datastructure.h"
#include "ising.h"
#include "metropolis.h"
#include "progressbar.h"
#include "ising-demo-metropolis.h"
#include "plot_fun.h"

int main (int argc, char **argv)
{
    struct ising model;

    int step_final = 20;
    double start_temp = 3.0;
    double final_temp = 4.0;

    double temp_array[step_final];
    double data_array[step_final];

    for (int step = 0; step <= step_final; step++)
    {

        double t = start_temp + (final_temp - start_temp) / (double) step_final * (double) step;                  // temperature (kT)
        double h = 0.;                   // magnetic field
        double j = 1.;                   // ferromagnetic coupling

        int lx = 32;
        int ly = lx;

        int mcsteps = 1024 * 128;        // 2^{17}

        if( ising_init (&model, lx, ly, j, h, t) != 0)
        {
        fprintf(stderr, "%s - cannot allocate spin array\n", argv[0]);
        argc = 1;
        exit(1);
        }

        print_header (&model, stderr);

        int therm_steps = mcsteps / 4;
        int report = therm_steps / 32;

        for (int s = 1; s <= therm_steps; s++)
        {
        (void) one_metropolis_step_per_spin (&model);
        if (s % report == 0)
        {
            text_progress_bar (s, therm_steps, 32, "Thermalization", stderr);
        }
        }
        fprintf (stderr, " Done\n");

        report = mcsteps / 32;
        double mav = 0., m2av = 0., eav = 0., e2av = 0.;
        double cv, xi;

        for (int s = 1; s <= mcsteps; s++)
        {
        (void) one_metropolis_step_per_spin (&model);
        double m = magnetization_per_spin (&model);
        double e = energy_per_spin (&model);

        mav += m;
        m2av += m * m;
        eav += e;
        e2av += e * e;

        if (s % report == 0)
        {
            text_progress_bar (s, mcsteps, 32, "Production", stderr);
        }

        }
        fprintf (stderr, " Done\n\n");

        mav /= mcsteps;
        eav /= mcsteps;
        m2av /= mcsteps;
        e2av /= mcsteps;

        // specific heat dE/dt
        cv = (e2av - eav*eav)/(t*t); 

        // magnetic susceptibility
        xi = (m2av - mav*mav)/t;

        printf ("%.3f   %.3f     % f  %f     % f  %f\n", t, h, mav, 
        xi, eav, cv);

        temp_array[step] = t;
        data_array[step] = cv;

        printf("t_before: %f\n", t);
        t = (final_temp - start_temp) / (double) step_final + t;
        printf("t_after: %f\n", t);

        ising_reinit(&model, h, t);
    }

    for (int step = 0; step < step_final; step++)
    {

        printf("%f, %f\n", temp_array[step], data_array[step]);

    }

    plot(temp_array, data_array, step_final);

    ising_free (model); 

    return 0;
}

void print_header (struct ising *model, FILE *fp)
{
    fprintf (fp, "\n");
    fprintf (fp, " Two-dimensional Ising model - Metropolis simulation\n");
    fprintf (fp, " ---------------------------------------------------\n");
    fprintf (fp, " T/J = %.3f  H/J = %.3f  %dx%d\n\n", model->T/model->J, 
        model->H/model->J, model->Lx, model->Ly);
}

