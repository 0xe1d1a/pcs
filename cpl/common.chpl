
const FPOPS_PER_POINT_PER_ITERATION = 12;

record params {
    const N : int;
    const M : int;
    const maxiter : int;
    const period : int;
    const threshold : real;
    const io_tmin : real;
    const io_tmax : real;
    const nthreads : int;
    const tinit;
    const tcond;
    const printreports;
}

record results {
    var niter : int;
    var tmin : real;
    var tmax : real;
    var maxdiff : real;
    var tavg : real;
    var time : real;
}

proc print_header(){ 
    writeln("Output from heat-1.15 (pcs2014@list.uva.nl):\n\n",
               "   Iterations",
               "        T(min)", 
               "        T(max)", 
               "       T(diff)", 
               "        T(avg)", 
               "          Time",
               "        FLOP/s");
}

proc report_results(p : params, r : results) 
{
    writeln(format("%-13zu ", r.niter), 
            format("% .6e " , r.tmin),
            format("% .6e " , r.tmax),
            format("% .6e " , r.maxdiff),
            format("% .6e " , r.tavg),
            format("% .6e " , r.time),
            format("% .6e" , p.N:real * p.M:real * 
                             (r.niter:real * FPOPS_PER_POINT_PER_ITERATION +
                              r.niter:real / p.period:real) / r.time));
}
