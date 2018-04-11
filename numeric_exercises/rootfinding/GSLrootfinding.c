#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

int rosenbrock_f(const gsl_vector* x, void* params, gsl_vector* f)
{
   const double x0=gsl_vector_get(x,0);
   const double x1=gsl_vector_get(x,1);
   const double y0=2*x0-2-400*x0*(x1-x0*x0);
   const double y1=200*(x1-x0*x0);
   gsl_vector_set(f,0,y0);
   gsl_vector_set(f,1,y1);
   return GSL_SUCCESS;
}

int rosenbrock_fun_df(const gsl_vector* x, void* params, gsl_matrix* J){
   const double x0=gsl_vector_get(x,0); 
   const double y=gsl_vector_get(x,1);
   gsl_matrix_set(J,0,0,1200*x0*x0-400*y+2); gsl_matrix_set(J,0,1,-400*x0);
   gsl_matrix_set(J,1,0,-400*x0); gsl_matrix_set(J,1,1,200);
   return GSL_SUCCESS;
}

int
rosenbrock_fun_fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J) {
   const double x0 = gsl_vector_get(x,0);
   const double x1 = gsl_vector_get(x,1);
   const double y0=2*x0-2-400*x0*(x1-x0*x0); 
   const double y1=200*(x1-x0*x0); 
   gsl_vector_set(f,0,y0);
   gsl_vector_set(f,1,y1);
   gsl_matrix_set(J,0,0,1200*x0*x0-400*x1+2); gsl_matrix_set(J,0,1,-400*x0);
   gsl_matrix_set(J,1,0,-400*x0); gsl_matrix_set(J,1,1,200);
   return GSL_SUCCESS;
}

int main()
{
   printf("Finding the root of the Rosenbrock gradient using GSL routines. The solver type is the discrete Newton algorithm using finite differences.\n");
   int vectorDim = 2;
   const gsl_multiroot_fsolver_type* solvertype = gsl_multiroot_fsolver_dnewton; 
   gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(solvertype,vectorDim);

   gsl_multiroot_function ros_fun;
   ros_fun.f=rosenbrock_f;
   ros_fun.n=vectorDim;
   ros_fun.params=NULL;

   double startx, starty;
   startx=3;
   starty=5;

   gsl_vector* startpoint = gsl_vector_alloc(vectorDim);
   gsl_vector_set(startpoint,0,startx);
   gsl_vector_set(startpoint,1,starty);
   gsl_multiroot_fsolver_set(solver,&ros_fun,startpoint);

   printf("Start guess (x,y) = (%g,%g)\n", startx,starty);

   double before_iterx,before_itery,after_iterx,after_itery,deltax,deltay;
   int flag, numberOfSteps;
   numberOfSteps=0;
   do
   {
      before_iterx = gsl_vector_get(solver->f, 0);
      before_itery = gsl_vector_get(solver->f, 1);
      gsl_multiroot_fsolver_iterate(solver);
      after_iterx = gsl_vector_get(solver->f, 0);
      after_itery = gsl_vector_get(solver->f, 1);
      deltax=after_iterx-before_iterx;
      deltay=after_itery-before_itery;

      if(deltax != 0 || deltay != 0)
      {
         numberOfSteps++;
      }

      flag = gsl_multiroot_test_residual(solver->f,1e-6);
      
   }
   while(flag==GSL_CONTINUE);
   double resultx = gsl_vector_get(solver->x,0);
   double resulty = gsl_vector_get(solver->x,1);

   printf("Number of steps taken: %i \n", numberOfSteps);
   printf("Root found at: (x,y)=(%g,%g)\n",resultx,resulty);

   printf("Now using the Newton method with supplied Jacobian using the same start guess:\n");

   const gsl_multiroot_fdfsolver_type* solvertype1 = gsl_multiroot_fdfsolver_newton; 
   gsl_multiroot_fdfsolver* solver1 = gsl_multiroot_fdfsolver_alloc(solvertype1,vectorDim);
   gsl_multiroot_function_fdf ros_funfdf = {&rosenbrock_f, &rosenbrock_fun_df, &rosenbrock_fun_fdf, vectorDim, NULL};
   gsl_multiroot_fdfsolver_set(solver1, &ros_funfdf, startpoint);

   numberOfSteps=0;
   do
   {
      before_iterx = gsl_vector_get(solver1->f, 0);
      before_itery = gsl_vector_get(solver1->f, 1);
      gsl_multiroot_fdfsolver_iterate(solver1);
      after_iterx = gsl_vector_get(solver1->f, 0);
      after_itery = gsl_vector_get(solver1->f, 1);
      deltax=after_iterx-before_iterx;
      deltay=after_itery-before_itery;

      if(deltax != 0 || deltay != 0)
      {
         numberOfSteps++;
      }

      flag = gsl_multiroot_test_residual(solver1->f,1e-6);
      
   }
   while(flag==GSL_CONTINUE);

   resultx = gsl_vector_get(solver1->x,0);
   resulty = gsl_vector_get(solver1->x,1);

   printf("Number of steps taken: %i \n", numberOfSteps);
   printf("Root found at: (x,y)=(%g,%g)\n",resultx,resulty);

   gsl_vector_free(startpoint); gsl_multiroot_fsolver_free(solver); gsl_multiroot_fdfsolver_free(solver1);
   return 0;
}