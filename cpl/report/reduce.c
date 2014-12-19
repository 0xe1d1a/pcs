class heatReduction : ReduceScanOp
{
    type eltType;
    var tmin: eltType = max(eltType);
    var tmax: eltType = min(eltType);
    var tsum: chpl__sumType(eltType);
    proc accumulate(val: eltType) {
        if (val < tmin) then tmin = val;
        if (val > tmax) then tmax = val;
        tsum += val;
    }
    proc combine(other: heatReduction) {
        if (other.tmin < tmin) then tmin = other.tmin;
        if (other.tmax > tmax) then tmax = other.tmax;
        tsum += other.tsum;
    }
    proc generate() {
        return new heatReductionResults(eltType, tmin, tmax, tsum);
    }
}

/*  code between  */

r.maxdiff = max reduce [ij in ProblemSpace] abs(dst[ij] - src[ij]);
var reduction = heatReduction reduce dst[ProblemSpace];
r.tmin = reduction.tmin;
r.tmax = reduction.tmax;
r.tavg = reduction.tsum / (p.N * p.M);

/* code continues */