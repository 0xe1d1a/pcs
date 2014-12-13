use common;
use Time;

proc do_compute(p : params) 
{
    /* Results structure */
    var r : results;

    /* Timer */
    var t : Timer;

    /* Start timer. */
    t.start();

    writeln(p.tinit);

    /* Intermediate timing results. */
    r.time = t.elapsed();

    /* Stop timer, set final timing */
    t.stop();
    r.time = t.elapsed();

    return r;
}
