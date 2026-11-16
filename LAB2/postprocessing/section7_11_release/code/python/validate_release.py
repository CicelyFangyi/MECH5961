"""Numerical interchange, exported artifact and source-integrity checks."""
from pathlib import Path
import hashlib,json,csv,sys
import numpy as np
from scipy.io import loadmat
from PIL import Image
import pymupdf
import xml.etree.ElementTree as ET
from analyse_release import csvwrite

ROOT=Path(__file__).resolve().parents[2]
def validate():
    G=loadmat(ROOT/'data_frozen'/'geometry.mat',squeeze_me=True);rows=[];datasets=[]
    units={'xy':'m','x':'m','xb':'m','xf':'m','frequencies':'Hz','fs':'Hz','L':'samples','c0':'m/s','h':'m',
           'f':'Hz','fineF':'Hz','C':'Pa^2 per FFT bin','macroCSM':'Pa^2 per FFT bin',
           'meanMaps':'Pa^2 m^2 under 1/(4 pi r) Green function','macroMaps':'Pa^2 m^2 under 1/(4 pi r) Green function',
           'fineMaps':'Pa^2 m^2 under 1/(4 pi r) Green function','spectrumFrequency':'Hz',
           'spectrumPower':'Pa^2 per FFT bin','run':'identifier','referenceDeviation':'dimensionless'}
    for file in sorted((ROOT/'data_frozen').glob('*.mat')):
        for key,value in loadmat(file,squeeze_me=True).items():
            if key.startswith('_'):continue
            a=np.asarray(value)
            datasets.append(dict(file=file.name,field=key,shape=list(a.shape),dtype=str(a.dtype),units=units[key]))
    (ROOT/'docs'/'dataset_schema.json').write_text(json.dumps(datasets,indent=2))
    for run in range(4):
        D=loadmat(ROOT/'data_frozen'/f'run{run}.mat',squeeze_me=True)
        meanMacro=D['macroMaps'].mean(axis=3)
        e=np.max(abs(meanMacro-D['meanMaps'][::2,::2,:]))/np.max(D['meanMaps'])
        assert e<1e-12;rows.append(dict(check=f'run{run}_macro_mean_vs_full',max_relative_error=float(e),result='PASS'))
        for j,f in enumerate(D['fineF']):
            k=int(np.flatnonzero(D['f']==f)[0]);p=D['fineMaps'][:,:,j]
            e=np.max(abs(p[::2,::2]-D['meanMaps'][:,:,k]))/p.max();assert e<1e-12
            coords=[(0,0),(131,277),(300,300),(600,600)]
            es=[]
            for iy,ix in coords:
                r=np.sqrt(np.sum((G['xy']-np.array([G['xf'][ix],G['xf'][iy]]))**2,axis=1)+1)
                t=np.exp(-1j*2*np.pi*float(f)*r/343)/(4*np.pi*r);w=t/np.sum(abs(t)**2)
                scalar=np.real(np.vdot(w,D['C'][:,:,k]@w));es.append(abs(scalar-p[iy,ix])/p.max())
            assert max(es)<1e-12
            rows.append(dict(check=f'run{run}_f{f}_grid_and_scalar',max_relative_error=float(max(e,max(es))),result='PASS'))
    with open(ROOT/'figures_tables'/'FIGURE_TABLE_MANIFEST.csv',encoding='utf8') as f:manifest=list(csv.DictReader(f))
    for item in manifest[:7]:
        name=item['filename'].split(';')[0][:-4];base=ROOT/'figures_tables'/name
        im=Image.open(base.with_suffix('.png'));doc=pymupdf.open(base.with_suffix('.pdf'))
        w,h=float(item['width_mm']),float(item['height_mm'])
        assert abs(im.width-w/25.4*600)<=1 and abs(im.height-h/25.4*600)<=1
        assert abs(doc[0].rect.width-w/25.4*72)<.01 and abs(doc[0].rect.height-h/25.4*72)<.01
        tree=ET.parse(base.with_suffix('.svg'));assert len(list(tree.iter('{http://www.w3.org/2000/svg}text')))>0
        sizes=[span['size'] for block in doc[0].get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']]
        assert min(sizes)>=5
        rows.append(dict(check=name+'_exports',result='PASS',width_mm=w,height_mm=h,png_width_px=im.width,png_height_px=im.height,min_pdf_font_pt=min(sizes)))
    assert not any(p.name.startswith('F8') for p in ROOT.rglob('*'))
    csvwrite(ROOT/'qa'/'acceptance.csv',rows)
    hashes=[dict(file=str(p.relative_to(ROOT)).replace('\\','/'),bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((ROOT/'data_frozen').glob('*.mat'))]
    csvwrite(ROOT/'docs'/'frozen_manifest.csv',hashes)
    print(json.dumps(dict(numerical_and_export_checks=len(rows),maximum_relative_error=max(r.get('max_relative_error',0) for r in rows),status='PASS')))

if __name__=='__main__':validate()
