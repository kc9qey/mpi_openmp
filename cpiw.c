#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <omp.h>

int main(int argc, char *argv[])
{
	int n, myid, numprocs, i;
	double PI25DT = 3.141592653589793238462643;
	double mypi, pi, h, sum, x;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// a parallel on the while loop will run and give results (perhaps correct) but 
// duplicates the work done by threads and is wasteful, the correct way to mark the
// outter loop is to just declare an omp region
#pragma omp 
	while (1) {
// without the #prama omp single to restrict MPI related activities
// to a single thread, code crashed when while loop was marked "parallel"
#pragma omp single
		if (myid == 0) {
			printf("Enter the number of intervals: (0 quits) ");
			scanf("%d",&n);
		}
#pragma omp single
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (n == 0)
			break;
		else {
			h   = 1.0 / (double) n;
			sum = 0.0;
#pragma omp parallel for reduction(+:sum)
			for (i = myid + 1; i <= n; i += numprocs) {
				x = h * ((double)i - 0.5);
				sum += (4.0 / (1.0 + x*x));
			        if (omp_get_thread_num() %2 == 0)
                                        putchar('.');
                                else
                                        putchar('-');
			}
			mypi = h * sum;
#pragma omp single
			MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0,
					MPI_COMM_WORLD);
#pragma omp single
			if (myid == 0)
				printf("\npi is approximately %.16f, Error is %.16f\n",
						pi, fabs(pi - PI25DT));
		}
	}
	MPI_Finalize();
	return 0;
}
