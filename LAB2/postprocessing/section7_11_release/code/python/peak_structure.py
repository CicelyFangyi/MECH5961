"""8-connected superlevel persistence on original MATLAB power maps."""
import numpy as np
from scipy import ndimage
from numba import njit

TIE = 64 * np.finfo(float).eps
LOCATION_TOL_M = 0.04  # identity tolerance, not a resolution distance

@njit(cache=False)
def merge_tree(v, ny, nx):
    order = np.argsort(-v, kind='mergesort')
    uf = np.full(v.size, -1, np.int64)
    elder = np.arange(v.size)
    parent = np.full(v.size, -1, np.int64)
    death = np.full(v.size, np.nan)
    for p in order:
        uf[p] = p
        iy, ix = p // nx, p % nx
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                yy, xx = iy+dy, ix+dx
                if (dy == 0 and dx == 0) or yy < 0 or yy >= ny or xx < 0 or xx >= nx:
                    continue
                q = yy*nx+xx
                if uf[q] < 0: continue
                a, b = p, q
                while uf[a] != a:
                    uf[a] = uf[uf[a]]; a = uf[a]
                while uf[b] != b:
                    uf[b] = uf[uf[b]]; b = uf[b]
                if a == b: continue
                ea, eb = elder[a], elder[b]
                if v[ea] < v[eb] or (v[ea] == v[eb] and ea > eb):
                    a,b,ea,eb = b,a,eb,ea
                parent[eb] = ea; death[eb] = v[p]
                uf[b] = a; elder[a] = ea
    return parent, death

def analyse(power, axis, reference=None):
    p = np.asarray(power, float)
    assert p.ndim == 2 and p.shape == (len(axis), len(axis))
    assert np.all(np.isfinite(p)) and p.min() >= 0 and p.max() > 0
    # Quantization only groups floating ties for topology; reported values use raw B.
    v = np.rint(p/p.max()/TIE)*TIE
    candidate = v == ndimage.maximum_filter(v, size=3, mode='nearest')
    labels,n = ndimage.label(candidate, np.ones((3,3)))
    parent,death = merge_tree(v.ravel(),*v.shape)
    peaks=[]
    for k in range(1,n+1):
        yy,xx=np.where(labels==k)
        index=int(np.min(yy*len(axis)+xx))
        y,x=divmod(index,len(axis))
        boundary=bool(np.any((xx==0)|(xx==len(axis)-1)|(yy==0)|(yy==len(axis)-1)))
        saddle=death[index]
        if np.isfinite(saddle) and v[y,x]-saddle<=TIE:
            continue  # Flat shoulders are not separate maxima.
        peaks.append(dict(index=index,x_m=float(axis[x]),y_m=float(axis[y]),
                          level_dB=float(10*np.log10(p[y,x]/p.max())),
                          prominence_dB=float(10*np.log10(v[y,x]/saddle)) if saddle>0 else None,
                          boundary=boundary,plateau_cells=len(xx)))
    peaks.sort(key=lambda z:(-z['level_dB'],z['index']))
    interior=[z for z in peaks if not z['boundary']]
    pair=None; depth=None; saddle=None
    if reference is not None and len(interior)>=2:
        # Primary pair is the strongest two interior maxima, chosen BEFORE matching.
        a,b=sorted(interior[:2],key=lambda z:z['x_m'])
        ref=np.array([[reference[0],0],[reference[1],0]])
        coords=np.array([[a['x_m'],a['y_m']],[b['x_m'],b['y_m']]])
        if np.all(np.linalg.norm(coords-ref,axis=1)<(reference[1]-reference[0])/2):
            pair=[a,b]
            def chain(index):
                d={index:np.inf}; bottleneck=np.inf
                while parent[index]>=0:
                    bottleneck=min(bottleneck,death[index]);index=int(parent[index]);d[index]=bottleneck
                return d
            ca,cb=chain(a['index']),chain(b['index'])
            saddle=max(min(ca[z],cb[z]) for z in ca.keys() & cb.keys())
            depth=float(min(a['level_dB'],b['level_dB'])-10*np.log10(saddle))
    return dict(peaks=peaks,pair=pair,saddle_depth_dB=depth,
                saddle_relative_power=None if saddle is None else float(saddle))

def same_pair(a,b,tol=LOCATION_TOL_M):
    if a['pair'] is None or b['pair'] is None:return False
    return all(np.hypot(x['x_m']-y['x_m'],x['y_m']-y['y_m'])<=tol
               for x,y in zip(a['pair'],b['pair']))

def sidelobe_conflict(result,control,grid_step):
    """Control flags matching secondary displacement AND no stronger persistence.
    Different control geometry prevents treating this alone as proof of an artifact.
    """
    if result['pair'] is None:return False
    p=sorted(result['pair'],key=lambda z:-z['level_dB'])
    c=[z for z in control['peaks'] if not z['boundary']]
    if len(c)<2:return False
    delta=np.array([p[1]['x_m']-p[0]['x_m'],p[1]['y_m']-p[0]['y_m']])
    for q in c[1:]:
        d=np.array([q['x_m']-c[0]['x_m'],q['y_m']-c[0]['y_m']])
        if np.linalg.norm(delta-d)<=max(LOCATION_TOL_M,2*grid_step):
            if q['prominence_dB'] is not None and result['saddle_depth_dB']<=q['prominence_dB']:
                return True
    return False

def minima_span(axis,profile,index):
    ids=np.flatnonzero((profile[1:-1]<profile[:-2])&(profile[1:-1]<profile[2:]))+1
    left=ids[ids<index];right=ids[ids>index]
    if not len(left) or not len(right):return None
    return float(axis[right[0]]-axis[left[-1]])
