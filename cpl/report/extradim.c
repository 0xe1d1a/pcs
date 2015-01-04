const TwoBigDomains = {0..p.N+1, 0..p.M+1, 0..1};
const BigDomain = TwoBigDomains[0..p.N+1, 0..p.M+1, 0];
var data: [TwoBigDomains] real;
var src = 0, dst = 1;

// buffer swapping
dst <=> src;
var cursrc: [BigDomain] => data[..,..,src];
var curdst: [BigDomain] => data[..,..,dst];

