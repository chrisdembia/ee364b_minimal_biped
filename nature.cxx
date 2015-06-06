
#include "psopt.h"



//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

double m = 50.0;
double g = 9.81;

double Lmax = 0.5;
double D = 0.5;
double V = 0.5;

double d = D * Lmax;
double v = V * sqrt(g * Lmax);

adouble endpoint_cost(adouble* initial_states, adouble* final_states, 
                      adouble* parameters,adouble& t0, adouble& tf, 
                      adouble* xad, int iphase)
{
   return 0;
} 

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters, 
                     adouble& time, adouble* xad, int iphase)
{
    if (iphase == 1) {

        adouble x = states[0];
        adouble y = states[1];

        adouble xd = states[2];
        adouble yd = states[3];

        adouble Ldot = (x * xd + y * yd) / sqrt(x * x + y * y);
        
        adouble F = controls[0];

        //return Ldot * Ldot / (m * g * d);
        //return F * F / (m * g * d);
        //return F / (m * g * d);
        //return F * Ldot / (m * g * d);
        /*
        if ((F * Ldot) > 0) {
            return F * Ldot / (m * g * d);
        } else {
            return 0.0;
        }
        */
        /*
        adouble W = F * Ldot;
        adouble cost;
        if (W > 0) {
            cost = 1.0 / 0.25 * W / (m * g * d) ;
        } else {
            cost = 1.0 / 0.25 * -W / (m * g * d);
        }
        adouble k = 5;
        return 2.0 / k * log(1 + exp(k * cost)) - cost - 2.0 / k * log(2);
        */
        return (F * Ldot) * (F * Ldot) / (m * g * d);
    }
    else {
        return 0.0;
    }
} 


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states, 
         adouble* controls, adouble* parameters, adouble& time, 
         adouble* xad, int iphase)
{

    adouble x = states[0];
    adouble y = states[1];

    adouble xd = states[2];
    adouble yd = states[3];

    derivatives[0] = xd;
    derivatives[1] = yd;

    derivatives[2] = 0;
    derivatives[3] = -g;

    if (iphase == 1) {
        adouble L = sqrt(x * x + y * y);
        adouble F = controls[0];
        derivatives[2] += F * x / L / m;
        derivatives[3] += F * y / L / m;

        path[0] = L;
    }
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states, 
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad, 
            int iphase) 

{
    e[0] = 0;
}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad)
{
    adouble xi_1[4];
    get_initial_states(xi_1, xad, 1);
    adouble xf_1[4];
    get_final_states(xf_1, xad, 1);
    adouble xi_2[4];
    get_initial_states(xi_2, xad, 2);
    adouble xf_2[4];
    get_final_states(xf_2, xad, 2);

    adouble tf_1 = get_final_time(xad, 1);
    adouble ti_2 = get_initial_time(xad, 2);

    linkages[0] = xf_1[0] - xi_2[0];
    linkages[1] = xf_1[1] - xi_2[1];
    linkages[2] = xf_1[2] - xi_2[2];
    linkages[3] = xf_1[3] - xi_2[3];
    linkages[4] = ti_2 - tf_1;

    // cyclic
    adouble params[1];
    get_parameters(params, xad, 1);
    linkages[5] = (xf_2[0] - xi_1[0]) - params[0];
    linkages[6] = xi_1[1] - xf_2[1];
    linkages[7] = xi_1[2] - xf_2[2];
    linkages[8] = xi_1[3] - xf_2[3];
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        		= "nature";
    problem.outfilename                 = "nature.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 2;
    problem.nlinkages                   = 9;
    //problem.nphases   			= 1;
    //problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 4;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents = 1;
    problem.phases(1).npath     		= 1;
    problem.phases(1).nparameters = 1;
    problem.phases(1).nodes             = 100;

    problem.phases(2).nstates   		= 4;
    problem.phases(2).ncontrols 		= 0;
    problem.phases(2).nevents = 1;
    problem.phases(2).npath     		= 0;
    problem.phases(2).nodes             = 100;
    /*
    */

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(1) = -2;
    problem.phases(1).bounds.upper.states(1) = 2;
    problem.phases(1).bounds.lower.states(2) = 0;
    problem.phases(1).bounds.upper.states(2) = 2;

    problem.phases(1).bounds.lower.states(3) = -5;
    problem.phases(1).bounds.upper.states(3) = 5;
    problem.phases(1).bounds.lower.states(4) = -5;
    problem.phases(1).bounds.upper.states(4) = 5;

    problem.phases(2).bounds.lower.states = problem.phases(1).bounds.lower.states;
    problem.phases(2).bounds.upper.states = problem.phases(1).bounds.upper.states;

    problem.phases(1).bounds.lower.controls(1) = 0.000;
    problem.phases(1).bounds.upper.controls(1) = 1000.0;

    problem.phases(1).bounds.lower.path(1) = 0;
    problem.phases(1).bounds.upper.path(1) = Lmax;

    problem.phases(1).bounds.lower.events(1) = 0;
    problem.phases(1).bounds.lower.events(1) = 0;
    problem.phases(2).bounds.upper.events(1) = 0;
    problem.phases(2).bounds.upper.events(1) = 0;

    problem.phases(1).bounds.lower.parameters(1) = d;
    problem.phases(1).bounds.upper.parameters(1) = d;

    double tstep = d / v;
    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.01;
    problem.phases(1).bounds.upper.EndTime      = tstep;

    problem.phases(2).bounds.lower.StartTime    = 0.01;
    problem.phases(2).bounds.upper.StartTime    = tstep;

    problem.phases(2).bounds.lower.EndTime      = tstep;
    problem.phases(2).bounds.upper.EndTime      = tstep;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    /*
    problem.phases(1).guess.parameters = zeros(4, 1);
    problem.phases(1).guess.parameters(1) = -0.1;
    */

    /*
    problem.phases(1).guess.states = zeros(4, 20);
    problem.phases(1).guess.states(1, colon()) = linspace(-1, 1, 20);
    */

    /*
    problem.phases(2).guess.time           = linspace(0.0, 3.0, 20); 
    DMatrix x0 = zeros(6, 20);
    x0(1, colon()) = linspace(0.5, 0.0, 20);
    x0(2, colon()) = linspace(0.0, -0.5, 20);
    x0(3, colon()) = linspace(0.0, -1.57, 20);
    problem.phases(2).guess.states         = x0;
    */

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-5;
    algorithm.nlp_method                  = "SNOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "numerical";
    //algorithm.mesh_refinement = "automatic";
    //TODO trapezoidal seems worse in this case.
    algorithm.collocation_method = "trapezoidal";
    //algorithm.collocation_method = "Hermite-Simpson";
    //algorithm.collocation_method = "Chebyshev";
    //algorithm.defect_scaling              = "jacobian-based";
    algorithm.switch_order = 0;
    algorithm.diff_matrix                 = "reduced-roundoff";
    

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, t, x1, t1, x2, u1, t2;

    //x 		= solution.get_states_in_phase(1);
    //t 		= solution.get_time_in_phase(1);
    x1 		= solution.get_states_in_phase(1);
    u1      = solution.get_controls_in_phase(1);
    t1 		= solution.get_time_in_phase(1);
    t1.Print("time1");

    x2 		= solution.get_states_in_phase(2);
    t2 		= solution.get_time_in_phase(2);
    t2.Print("time2");

    x = (x1 || x2);
    t = (t1 || t2);
    /*
    */

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x1.Save("nature_x1.dat");
    u1.Save("nature_u1.dat");
    t1.Save("nature_t1.dat");

    x2.Save("nature_x2.dat");
    t2.Save("nature_t2.dat");

    x.Save("nature_x.dat");
    t.Save("nature_t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    plot(t1,u1,problem.name+": control","time (s)", "control 1", "F");
    /*
    */

    plot(t,x(colon(1, 2), colon()),problem.name+": state","time (s)", "state", "x y");
    plot(t,x(colon(3, 4), colon()),problem.name+": state","time (s)", "state", "xd yd");
    //plot(t,x(2, colon()),problem.name+": state","time (s)", "state", "y");
    
    
    /*
    plot(t,u,problem.name+": control","time (s)", "control 1", "u1",
	      "pdf", "nature_controls.pdf");

    plot(t,x,problem.name+": states","time (s)", "states", "x1 x2",
	      "pdf", "nature_states.pdf");
          */

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
