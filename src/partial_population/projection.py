"""Certified outer enclosures of the fixed-count joint e-value projection.

Floating-point conic optimization only proposes certificates. FLINT/Arb
validates continuum kernel bounds; entropy duality then certifies endpoints.
All calibration choices are fixed by design, never fitted to the fingerprint.
"""
from dataclasses import dataclass, asdict
from fractions import Fraction
import hashlib
import json
import math
import warnings

import numpy as np
from scipy.special import gammaln, logsumexp
from flint import arb, ctx

SCHEMA = 'partial-population-projection-v1'


def outward(v, upper=True):
    """A binary64 bound enclosing an Arb ball, with a checked conversion."""
    x = float(v.upper() if upper else v.lower())
    if not math.isfinite(x):
        raise ArithmeticError('Nonfinite certified quantity')
    x = math.nextafter(x, math.inf if upper else -math.inf)
    assert arb(x) >= v if upper else arb(x) <= v
    return x


@dataclass(frozen=True)
class Design:
    n: int
    target: str = 'missing_mass'
    m: int = 0
    alpha: str = '1/20'
    # Fixed scales for the singleton/mass or truncated Good--Toulmin score.
    scales: tuple = ()
    # Optional predeclared linear fingerprint score. Empty uses the classical
    # singleton or three-term Good--Toulmin score. Adaptive fitting needs its
    # own calibration; supplying data-fitted weights does not provide that.
    score_weights: tuple = ()

    def __post_init__(self):
        if type(self.n) is not int or self.n < 1:
            raise ValueError('n must be a positive integer')
        if self.target not in ('missing_mass', 'discoveries'):
            raise ValueError('unknown target')
        if type(self.m) is not int or self.m < 0:
            raise ValueError('m must be a nonnegative integer')
        if self.target == 'discoveries' and self.m < 1:
            raise ValueError('future sample size must be positive')
        if self.target == 'missing_mass' and self.m != 0:
            raise ValueError('missing mass uses m=0')
        if not 0 < Fraction(self.alpha) < 1:
            raise ValueError('alpha must lie in (0,1)')
        if not self.scales:
            base = Fraction(1,max(4,math.ceil(math.sqrt(self.n/(math.log(self.n+1)+3)))))
            chosen = [base*v for v in (Fraction(-2),Fraction(-1),Fraction(-1,2),Fraction(-1,4),Fraction(1,4),Fraction(1,2),Fraction(1),Fraction(2))]
            object.__setattr__(self,'scales',tuple(str(s) for s in chosen))
        object.__setattr__(self, 'scales', tuple(self.scales))
        object.__setattr__(self, 'score_weights', tuple(str(Fraction(s)) for s in self.score_weights))
        if not self.scales or any(abs(Fraction(s)) > Fraction(1,2) or not Fraction(s) for s in self.scales):
            raise ValueError('nonzero fixed scales must have magnitude <= 1/2')
        if not (any(Fraction(s)<0 for s in self.scales) and any(Fraction(s)>0 for s in self.scales)):
            raise ValueError('two-sided projection requires both scale signs')

    @property
    def target_scale(self):
        return 1 if self.target == 'missing_mass' else self.m

    def directions(self):
        ans = []
        for raw in self.scales:
            s = Fraction(raw)
            if self.target == 'missing_mass':
                ans.append(([-s*Fraction(w) for w in self.score_weights] if self.score_weights else [-s], s*self.n))
            else:
                ratio = Fraction(self.m, self.n)
                # Three preassigned coordinates; truncation is part of the
                # score, not a claim that the predictor is unbiased.
                ans.append(([-s*Fraction(w) for w in self.score_weights] if self.score_weights else [-s*((-1)**(r+1))*ratio**r for r in range(1,4)], s))
        return ans


def B(q):
    return arb(str(q)) if isinstance(q, Fraction) else arb(q)


def factor_ball(x, design, direction, derivative=False):
    lam, eta = direction
    t = arb(design.n)
    e = (-t*x).exp()
    f, df = arb(1), arb(0)
    for r, lr in enumerate(lam, 1):
        c = B(lr).expm1()*t**r/math.factorial(r)
        f += c*e*x**r
        df += c*e*(r*x**(r-1)-t*x**r)
    if design.target == 'missing_mass':
        # Difference of exponentials avoids exp(eta*x) overflow.
        a = (B(eta)-t)
        f += ((a*x).exp()*(-(-B(eta)*x).expm1())
              if eta > 0 else e*(B(eta)*x).expm1())
        df += a*(a*x).exp()+t*e
    else:
        u, c = arb(design.m), B(eta).expm1()
        f += c*e*(-(-u*x).expm1())
        df += c*(-t*e+(t+u)*(-(t+u)*x).exp())
    return (f, df) if derivative else f


def kernel_ball(lo, hi, design, direction):
    """Enclose kappa on [lo,hi], including its removable singularity at 0."""
    x = B(lo).union(B(hi))
    if lo == 0:
        f, df = factor_ball(x, design, direction, True)
        # log f(x)/x = integral_0^1 f'(sx)/f(sx) ds, f(0)=1.
        return df/f
    return factor_ball(x, design, direction).log()/x


def kernel_float(x, design):
    x = np.asarray(x, dtype=float)
    rows = []
    t = design.n
    for lam, eta in design.directions():
        z = t*x
        delta = sum(np.exp(-z)*z**r/math.factorial(r)*math.expm1(float(lr))
                    for r, lr in enumerate(lam, 1))
        if design.target == 'missing_mass':
            delta += (np.exp((float(eta)-t)*x)*(-np.expm1(-float(eta)*x))
                      if eta > 0 else np.exp(-t*x)*np.expm1(float(eta)*x))
            zero = t*math.expm1(float(lam[0]))+float(eta)
        else:
            delta += np.exp(-z)*(-np.expm1(-design.m*x))*math.expm1(float(eta))
            zero = t*math.expm1(float(lam[0]))+design.m*math.expm1(float(eta))
        with np.errstate(divide='ignore', invalid='ignore'):
            rows.append(np.where(x == 0, zero, np.log1p(delta)/x))
    return np.array(rows)


def penalty_ball(design):
    def one(n):
        return arb(n)-n*arb(n).log()+arb(n+1).lgamma()
    return one(design.n)+(one(design.m) if design.target == 'discoveries' else 0)


def penalty_float(design):
    def one(n):
        return n-n*math.log(n)+gammaln(n+1)
    return one(design.n)+(one(design.m) if design.target == 'discoveries' else 0)


def design_key(design):
    fields=asdict(design)
    if not fields['score_weights']:fields.pop('score_weights')
    return hashlib.sha256(json.dumps(fields, sort_keys=True).encode()).hexdigest()[:16]


def make_kernel_table(design, resolution=None):
    """Shared continuum enclosures; a grid alone is never a certificate.

    Endpoints are exact binary floats. Subdivision retries any nonfinite
    interval calculation. No probability floor is imposed.
    """
    if resolution is None:
        resolution = max(1400,min(14000,int(16*math.sqrt(design.n))))
    with ctx.workprec(128):
        smallest = min(1e-10, 1e-6/design.n)
        edges = np.unique(np.r_[0., np.geomspace(smallest, 1., resolution)])
        cells, bounds = [], []
        pending = [(float(a),float(b),0) for a,b in zip(edges[:-1],edges[1:])]
        directions = design.directions()
        while pending:
            lo, hi, depth = pending.pop()
            vals = [kernel_ball(lo, hi, design, d) for d in directions]
            if any(not v.is_finite() for v in vals):
                if depth >= 30:
                    raise ArithmeticError('Kernel interval cannot be bounded')
                mid = (lo+hi)/2
                pending.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
                continue
            cells.append([lo,hi])
            bounds.append([outward(v) for v in vals])
        order = np.argsort([c[0] for c in cells])
        return {'schema':SCHEMA, 'design':asdict(design), 'precision_bits':128,
                'cells':[cells[i] for i in order],
                'upper':[bounds[i] for i in order]}


def verify_kernel_table(table):
    design = Design(**table['design'])
    if table['schema'] != SCHEMA or len(table['cells']) != len(table['upper']):
        raise ValueError('invalid kernel table')
    end = 0.
    with ctx.workprec(128):
        for (lo,hi), row in zip(table['cells'],table['upper']):
            if lo != end or not lo < hi <= 1 or len(row) != len(design.scales):
                raise ValueError('kernel cells must partition [0,1]')
            for bound, direction in zip(row,design.directions()):
                if not math.isfinite(bound) or not B(bound) >= kernel_ball(lo,hi,design,direction):
                    raise ValueError('false continuum bound')
            end = hi
    if end != 1:
        raise ValueError('incomplete probability domain')
    return True


def rational_weights(weights, denominator=10**12):
    w = np.maximum(np.asarray(weights,float),0)
    if not np.isfinite(w).all() or w.sum() == 0:
        raise ValueError('invalid proposed weights')
    ints = np.floor(w/w.sum()*denominator).astype(np.int64)
    ints[int(np.argmax(w))] += denominator-int(ints.sum())
    return [int(x) for x in ints]


def fingerprint(counts, dimension):
    return [sum(int(c == r) for c in counts) for r in range(1,dimension+1)]


def validate_fingerprint(phi, design):
    d = len(design.directions()[0][0])
    if len(phi) != d or any(type(v) is not int or v < 0 for v in phi):
        raise ValueError('invalid fingerprint')
    if sum((r+1)*v for r,v in enumerate(phi)) > design.n:
        raise ValueError('fingerprint exceeds known sample size')


def score_ball(phi, lam):
    return sum((v*B(l) for v,l in zip(phi,lam)),arb(0))


def dual_endpoint(design, phi, qints, table):
    """Entropy lower bound on log mixture; returns a rigorous endpoint ball."""
    validate_fingerprint(phi,design)
    if len(qints) != len(design.scales) or any(type(v) is not int or v < 0 for v in qints) or sum(qints) != 10**12:
        raise ValueError('dual weights must be nonnegative and sum exactly to 10^12')
    with ctx.workprec(128):
        qs = [arb(q)/10**12 for q in qints]
        # Upper bounds on all kernels, jointly valid on each probability cell.
        support = max((sum((q*B(u) for q,u in zip(qs,row)),arb(0)).upper()
                       for row in table['upper']))
        intercept, slope = arb(0), arb(0)
        J = len(qs)
        penalty = penalty_ball(design)
        for q,qi,(lam,eta) in zip(qs,qints,design.directions()):
            if qi:
                intercept += q*(score_ball(phi,lam)-penalty-arb(J).log()-q.log())
                slope += q*B(eta)*design.target_scale
        if not (slope > 0 or slope < 0):
            raise ValueError('dual slope has no certified sign')
        value = (-B(Fraction(design.alpha)).log()-intercept+support)/slope
        return value, slope > 0


def feasible_witness(design, phi, nodes, mass, z):
    """Verify one actual point in the relaxed convex region (not integer-realizability)."""
    if len(nodes) != len(mass) or not nodes or any(not 0 <= x <= 1 for x in nodes):
        return False
    if any(type(w) is not int or w < 0 for w in mass) or sum(mass) != 10**12 or not 0 <= z <= 1:
        return False
    with ctx.workprec(128):
        value = arb(0)
        for lameta in design.directions():
            k = sum((arb(w)/10**12*kernel_ball(x,x,design,lameta)
                     for x,w in zip(nodes,mass) if w),arb(0))
            lam,eta = lameta
            value += (score_ball(phi,lam)+B(eta)*design.target_scale*B(z)-penalty_ball(design)-k).exp()/len(design.scales)
        return bool(value <= 1/B(Fraction(design.alpha)))


class Projector:
    """Reusable fixed-design optimizer; certificates carry their own inputs."""
    def __init__(self, design, table=None, grid_size=100):
        import cvxpy as cp
        self.design = design
        self.table = make_kernel_table(design) if table is None else table
        if self.table['design'] != json.loads(json.dumps(asdict(design))):
            # tuples in freshly constructed tables are harmless
            if Design(**self.table['design']) != design:
                raise ValueError('table/design mismatch')
        self.nodes = np.unique(np.r_[0.,np.geomspace(1e-8/design.n,1.,grid_size)])
        self.K = kernel_float(self.nodes,design)
        self.eta = np.array([float(eta)*design.target_scale for _,eta in design.directions()])
        self.lam = np.array([[float(x) for x in lam] for lam,_ in design.directions()])
        self.mass = cp.Variable(len(self.nodes),nonneg=True)
        self.z = cp.Variable()
        self.score = cp.Parameter(len(design.scales))
        self.sign = cp.Parameter()
        expr = self.score+self.eta*self.z-self.K@self.mass-penalty_float(design)-math.log(len(design.scales))
        self.constraints = [cp.sum(self.mass)==1,self.z>=0,self.z<=1,
                            cp.log_sum_exp(expr) <= -math.log(float(Fraction(design.alpha)))]
        self.problem = cp.Problem(cp.Minimize(self.sign*self.z),self.constraints)

    def interval(self, phi):
        from cvxpy.error import SolverError
        validate_fingerprint(phi,self.design)
        self.score.value = self.lam@np.asarray(phi)
        bounds, certs = [], []
        for side in ('lower','upper'):
            self.sign.value = 1 if side == 'lower' else -1
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    self.problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,
                                       tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=200)
                z = float(self.z.value)
                mass = np.maximum(np.asarray(self.mass.value),0)
                mass /= mass.sum()
                logs = self.score.value+self.eta*z-self.K@mass
                qints = rational_weights(np.exp(logs-logsumexp(logs)))
                value, positive = dual_endpoint(self.design,phi,qints,self.table)
                if positive != (side == 'upper'):
                    raise ValueError('wrong dual sign')
                bound = outward(value, upper=positive)
                bound = min(1.,max(0.,bound))
                # Move into the numerical region and verify a witness independently.
                witness = None
                wm = rational_weights(mass)
                for step in (1e-7,1e-5,1e-3,.01,.05):
                    wz = min(1.,max(0.,z+(step if side=='lower' else -step)))
                    if feasible_witness(self.design,phi,self.nodes.tolist(),wm,wz):
                        witness = {'nodes':self.nodes.tolist(),'mass':wm,'z':wz}
                        break
                cert = {'side':side,'q':qints,'bound_unit':bound,
                        'proposed_unit':z,'witness':witness,'solver_status':self.problem.status}
            except (ValueError,ArithmeticError,TypeError,OverflowError,SolverError) as exc:
                candidates=[]
                for j,eta in enumerate(self.eta):
                    if (eta>0) != (side=='upper'): continue
                    qi=[0]*len(self.eta);qi[j]=10**12
                    v,pos=dual_endpoint(self.design,phi,qi,self.table)
                    candidates.append((min(1.,max(0.,outward(v,upper=pos))),qi))
                bound,qints=(min(candidates,key=lambda x:x[0]) if side=='upper' else max(candidates,key=lambda x:x[0]))
                cert={'side':side,'q':qints,'bound_unit':bound,'witness':None,
                      'solver_status':'failed; certified one-direction bound','reason':str(exc)}
            bounds.append(bound)
            certs.append(cert)
        # Empty projections always report full support, preserving unconditional coverage.
        if bounds[0] > bounds[1]:
            bounds = [0.,1.]
            fallback = 'empty-certified-projection'
        else:
            fallback = None
        scale = self.design.target_scale
        reported = [bounds[0]*scale,bounds[1]*scale]
        # Outward rounding also covers the final unit conversion.
        reported = [max(0.,math.nextafter(reported[0],-math.inf)),
                    min(float(scale),math.nextafter(reported[1],math.inf))]
        if self.design.target == 'discoveries':
            reported = [math.ceil(reported[0]),math.floor(reported[1])]
            if reported[0] > reported[1]:
                reported = [0,scale]
                fallback = 'empty-integer-projection'
        return {'schema':SCHEMA,'design':asdict(self.design),'fingerprint':phi,
                'interval':reported,'endpoints':certs,'fallback':fallback}


def verify_interval(cert,table):
    """Independent of CVXPY; call verify_kernel_table once per shared table."""
    design = Design(**cert['design'])
    if design != Design(**table['design']) or cert['schema'] != SCHEMA:
        raise ValueError('certificate/table mismatch')
    phi = cert['fingerprint']
    validate_fingerprint(phi,design)
    expected = []
    for side,endpoint in zip(('lower','upper'),cert['endpoints']):
        if endpoint['side'] != side:
            raise ValueError('endpoint order')
        if 'fallback' in endpoint:
            bound = 0. if side=='lower' else 1.
        else:
            value,positive = dual_endpoint(design,phi,endpoint['q'],table)
            if positive != (side=='upper'):
                raise ValueError('wrong certificate sign')
            bound = min(1.,max(0.,outward(value,upper=positive)))
            if endpoint.get('witness'):
                w = endpoint['witness']
                if not feasible_witness(design,phi,w['nodes'],w['mass'],w['z']):
                    raise ValueError('invalid primal witness')
        if endpoint['bound_unit'] != bound:
            raise ValueError('false endpoint')
        expected.append(bound)
    if len(cert['endpoints']) != 2:
        raise ValueError('two endpoint certificates required')
    if expected[0]>expected[1]: expected=[0.,1.]
    scale=design.target_scale
    expected=[max(0.,math.nextafter(expected[0]*scale,-math.inf)),
              min(float(scale),math.nextafter(expected[1]*scale,math.inf))]
    if design.target=='discoveries':
        expected=[math.ceil(expected[0]),math.floor(expected[1])]
        if expected[0]>expected[1]: expected=[0,scale]
    if expected != cert['interval']:
        raise ValueError('reported interval differs from certificate')
    return True
