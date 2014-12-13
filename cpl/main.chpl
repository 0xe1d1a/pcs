use common;
use compute;

config const Gm = 150;
config const Gn = 100;
config const Gi = 42;
config const Gk = 1000;
config const Ge = 0.1;
config const GL = -100.0;
config const GH = 100.0;
config const Gp = 1;
config const Gc = "pattern_100x150.pgm";
config const Gt = "pattern_100x150.pgm";
config const Gh = 0;

proc usage()
{
    writeln("Usage:  [OPTION]...\n",
           "\n",
           "  -sGn=NUM     Set cylinder height to ROWS.\n",
           "  -sGm=NUM     Set cylinder width to COLUMNS.\n",
           "  -sGi=NUM     Set the maximum number of iterations to NUM.\n",
           "  -sGk=NUM     Set the reporting period to NUM.\n",
           "  -sGe=NUM     The the convergence threshold to NUM.\n",
           "  -sGc=FILE    Read conductivity values from FILE.\n",
           "  -sGt=FILE    Read initial temperature values from FILE.\n",
           "  -sGL=NUM     Coldest temperature in input/output images.\n",
           "  -sGH=NUM     Warmest temperature in input/output images.\n",
           "  -sGp=NUM     Number of threads to use (when applicable).\n",
           "  -sGh=1       Print this help.\n"
           );
    exit(0);
}

proc readpgm(fname, height, width, min, max)
{
    var f = open(fname, iomode.r).reader();

    var header = f.readln(string);
    var (w,h) = f.readln(int, int);
    var maxv = f.readln(int);

    var im : [1..w, 1..w] real;

    for i in 1..h {
        for j in 1..w {
            im[i,j] = min + f.read(int):real * (max - min) / maxv:real;
        }
    }

    return im;
}


proc read_parameters()
{

    writeln("Parameters:\n",
           "  -sGn=", Gn, " # number of rows\n",
           "  -sGm=", Gm, " # number of columns\n",
           "  -sGi=", Gi, " # maximum number of iterations\n",
           "  -sGk=", Gk, " # reporting period\n",
           "  -sGe=", Ge, " # convergence threshold\n",
           "  -sGc=", Gc, " # input file for conductivity\n",
           "  -sGt=", Gt, " # input file for initial temperatures\n",
           "  -sGL=", GL, " # coolest temperature in input/output\n",
           "  -sGH=", GH, " # highest temperature in input/output\n",
           "  -sGp=", Gp, " # number of threads (if applicable)");

    var tinit = readpgm(Gt, Gn, Gm, GL, GH);
    var tcond = readpgm(Gc, Gn, Gm, 0.0, 1.0);

    return (tinit, tcond);
}

proc fill_record(tinit, tcond)
{
    var p = new params(N = Gn, M = Gm, maxiter = Gi, period = Gk, threshold = Ge, nthreads = Gp, 
                       io_tmin = GL, io_tmax = GH, tinit = tinit, tcond = tcond);

    return p;
}

proc main()
{
    if(Gh == 1) then {
      usage();
    }

    var (tinit, tcond) = read_parameters();
    var p = fill_record(tinit, tcond);

    print_header();

    var r = do_compute(p);

    report_results(p, r);
}


