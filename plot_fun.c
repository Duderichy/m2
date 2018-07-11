#include <stdio.h>
#include <stdlib.h>

//function which pipes data straight to gnuplot so that the program creates the plot automatically
void plot (double x[], double y[], int points);
void plot (double x[], double y[], int points)
{
    char *commandsForGnuplot[] =
        { "set title \"Specific Heat Vs. Temperature\"", "plot 'res' w l" };
    FILE *temp = fopen ("res", "w");    // write coordinates here.
    FILE *gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;

    for (int q = 0; q < points; q++)
    {
        fprintf (temp, "%lf %lf \n", x[q], y[q]);       //Write the data to a temporary file
    }
    for (i = 0; i < 2; i++)
    {
        fprintf (gnuplotPipe, "%s \n", commandsForGnuplot[i]);  //Send commands to gnuplot one by one.
    }
}
