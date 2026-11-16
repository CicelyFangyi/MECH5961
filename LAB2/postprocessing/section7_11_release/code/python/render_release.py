"""Independent release reproduction: python code/python/render_release.py --from-frozen.
All scientific values are recalculated from the bundled MATLAB-derived arrays.
No MATLAB, TDMS, skill installation, private path or network is required.
"""
from pathlib import Path
import argparse,json,os,sys,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from analyse_release import process,csvwrite

ROOT=Path(__file__).resolve().parents[2]
BLUE='#1B619A';OCHRE='#D9A42B';GRAY='#555555'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
    'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.8,'lines.linewidth':1.4,
    'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,
    'mathtext.fontset':'dejavusans','savefig.facecolor':'white'})

def layout(rows,cols,width=183,height=165):
    fig,axs=plt.subplots(rows,cols,figsize=(width/25.4,height/25.4),squeeze=False)
    fig.subplots_adjust(left=.11,right=.98,bottom=.16,top=.97,wspace=.38,hspace=.65)
    return fig,axs

def condition(ax,text):
    ax.text(.5,-.32,text,transform=ax.transAxes,ha='center',va='top',fontsize=9)

def alignment(fig,name):
    fig.canvas.draw();boxes=np.array([a.get_position().bounds for a in fig.axes])
    pt=boxes*np.array([fig.get_figwidth()*72,fig.get_figheight()*72]*2)
    max_width=float(np.ptp(pt[:,2]));max_height=float(np.ptp(pt[:,3]))
    assert max_width<=1.5 and max_height<=1.5
    report=dict(verdict='PASS' if len(boxes)>1 else 'NOT APPLICABLE',
                plot_areas_pt=pt.tolist(),width_spread_pt=max_width,height_spread_pt=max_height)
    (ROOT/'qa'/f'{name}.alignment.json').write_text(json.dumps(report,indent=2))
    # Optional local skill audit. Scientific reproduction does not depend on it.
    if os.environ.get('NATURE_FIGURE_SCRIPTS'):
        sys.path.insert(0,os.environ['NATURE_FIGURE_SCRIPTS'])
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(fig,json_out=str(ROOT/'qa'/f'{name}.skill-alignment.json'),
                                           tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)

def save(fig,name,manifest,title,section,source,method,limits):
    alignment(fig,name)
    for ext in ['png','svg','pdf']:
        kw={'dpi':600} if ext=='png' else {}
        if ext=='pdf':kw['metadata']={'CreationDate':None,'ModDate':None}
        fig.savefig(ROOT/'figures_tables'/f'{name}.{ext}',**kw)
    manifest.append(dict(filename=f'{name}.png;{name}.svg;{name}.pdf',english_title=title,
        placement=section,source=source,method_parameters_units=method,limitations_review=limits,
        width_mm=fig.get_figwidth()*25.4,height_mm=fig.get_figheight()*25.4,dpi=600,
        python_source='code/python/render_release.py'))
    plt.close(fig)

def level(p):return 10*np.log10(np.maximum(p,np.finfo(float).tiny)/p.max())

def render():
    G,runs,fine,tables=process(ROOT);x=G['xf'];xy=G['xy'];manifest=[]
    fig,axs=layout(1,1,100,95);ax=axs[0,0]
    fig.subplots_adjust(left=.2,bottom=.18,right=.96,top=.97)
    ax.scatter(xy[:,0],xy[:,1],s=12,color=BLUE);ax.plot(0,0,'+',color=GRAY,ms=6)
    ax.set(xlabel='x (m)',ylabel='y (m)',xlim=(-.4,.4),ylim=(-.4,.4),xticks=[-.3,0,.3],yticks=[-.3,0,.3]);ax.set_aspect('equal')
    save(fig,'F7A_array_geometry',manifest,'Actual 56-channel microphone geometry','Section 7 instrumentation',
         'geometry.mat:xy','Untransformed teacher coordinate rows; m; equal aspect','Physical channel wiring unverified; Dmax and Dx distinct in T7')
    fig,axs=layout(2,2,183,145)
    for i,ax in enumerate(axs.flat):
        quantitative=i in [1,2]
        ax.axhline(0,color=GRAY,lw=1.2);ax.axhline(1,color='.8',lw=.8)
        if quantitative:
            ref=[-.2,0] if i==1 else [-.2,.2]
            ax.plot(ref[0],1,'o',color=BLUE,ms=6);ax.plot(ref[1],1,'s',color=OCHRE,ms=6)
            ax.set(xlabel='x (m)',ylabel='z (m)',xticks=[-.2,0,.2],yticks=[0,1])
        else:
            ax.plot(0,1,'D' if i==0 else 'x',color=GRAY,ms=8,mew=1.5)
            ax.set(xticks=[],yticks=[])
            for s in ax.spines.values():s.set_visible(False)
        ax.set(xlim=(-.4,.4),ylim=(-.12,1.2))
        condition(ax,['(a) Single loudspeaker','(b) Δx = 0.20 m','(c) Δx = 0.40 m','(d) Two-bladed propeller'][i])
    save(fig,'F7B_four_source_geometries',manifest,'Four experimental geometries','Section 7 setup',
        'User/teammate confirmed runs 1,2; qualitative runs 0,3','Quantitative x-z plane, y=0; h=1 m',
        'Panels a,d are qualitative schematics: marker placement is not a measured lateral coordinate. Left source common to b,c.')
    fig,axs=layout(2,3)
    for row,run in enumerate([2,1]):
        for col,f in enumerate([1000,3000,6000]):
            ax=axs[row,col];D=fine[(run,f)];L=level(D['power']);profile=L[len(x)//2,:]
            ax.plot(x,profile,color=BLUE)
            for z in D['ref']:ax.axvline(z,ymin=.4,color='.75',lw=.8,ls=':')
            if D['supported']:
                a=D['analysis']
                for p in a['pair']:
                    ix=np.argmin(abs(x-p['x_m']));ax.plot(x[ix],profile[ix],'o',color=OCHRE,ms=4)
                # Sparse annotations in the lower clear region; exact values in CSV.
                ax.text(.04,.07,f"$\\hat x_1$ = {a['pair'][0]['x_m']:.3f} m\n$\\hat x_2$ = {a['pair'][1]['x_m']:.3f} m",
                        transform=ax.transAxes,fontsize=8)
                ax.text(.04,.26,f"$\\Delta L_s$ = {a['saddle_depth_dB']:.2f} dB",transform=ax.transAxes,fontsize=8)
            ax.set(xlim=(-.5,.5),ylim=(-40,1),xticks=[-.4,0,.4],yticks=[-40,-20,0],xlabel='x (m)')
            if col==0:ax.set_ylabel('Relative beamforming level (dB)')
            condition(ax,f"({chr(97+row*3+col)}) Δx = {D['ref'][1]+.2:.2f} m, f = {f/1000:g} kHz")
    save(fig,'F11A_two_source_profiles',manifest,'Horizontal two-source beamforming profiles','Appendix Q1',
        'run1.mat;run2.mat:fineMaps','Unmodified y=0 slice; 601 grid; normalized 10 log10 B/max B',
        'Markers project supported 2D peaks onto y=0; actual peak y is nonzero. ΔLs is a 2D saddle descriptor, not linecut depth.')
    fig,axs=layout(1,2,183,90)
    fig.subplots_adjust(left=.17,right=.98,bottom=.28,top=.96,wspace=.4)
    for ax,run,letter in zip(axs.flat,[2,1],'ab'):
        rows=[r for r in tables['T11_frequency_evidence'] if r['run']==run]
        ax.plot([r['frequency_Hz']/1000 for r in rows],[r['evidence_stage'] for r in rows],color=BLUE,marker='o',ms=3)
        ax.set(xlim=(1,6),ylim=(-.2,2.2),xticks=[1,2,3,4,5,6],yticks=[0,1,2],
               yticklabels=['No pair','Mean only','Stable blocks'],xlabel='f (kHz)')
        condition(ax,f'({letter}) Δx = {0.4 if run==2 else 0.2:.2f} m')
    save(fig,'F11B_robust_resolvability',manifest,'Deterministic two-dimensional peak stability','Appendix Q1',
         'run0..2.mat:meanMaps,macroMaps','301 full-record; 151 five 2 s blocks, two splits and five leave-one-out aggregates; 250 Hz step',
         'Ordered evidence categories, not probabilities. Connecting lines guide reading, not interpolation. f_onset indeterminate; no claim between sampled frequencies.')
    P=runs[3];sf=P['spectrumFrequency'].astype(float)
    assert np.all(P['spectrumPower']>0)
    sp=10*np.log10(P['spectrumPower']/(20e-6)**2)
    fig,axs=layout(1,2,183,90);fig.subplots_adjust(left=.11,right=.98,bottom=.27,top=.97,wspace=.38)
    for j,ax in enumerate(axs.flat):
        mask=(sf>= (50 if j==0 else 100))&(sf<=(10000 if j==0 else 1000));scale=1000 if j==0 else 1
        ax.plot(sf[mask]/scale,sp[mask],color=BLUE,lw=1.2)
        guides=P['fineF'] if j==0 else [200,400]
        for k,f in enumerate(guides):
            ax.axvline(f/scale,color=OCHRE if k%2==0 else GRAY,ls=':' if k%2==0 else '--',lw=.8)
            if j==0:ax.plot(f/scale,sp[np.argmin(abs(sf-f))],'o',color=OCHRE,ms=3)
        ax.set(xlabel='f (kHz)' if j==0 else 'f (Hz)',ylabel='SPL (dB re 20 µPa)',ylim=(0,65))
        if j==0:
            ax.set(xlim=(.05,10),xticks=[0,2,4,6,8,10])
        else:
            ax.set(xlim=(100,1000),xticks=[200,400,600,800,1000])
            ax.text(.13,.91,'fBPF',transform=ax.transAxes,fontsize=8)
            ax.text(.36,.91,'2 fBPF',transform=ax.transAxes,fontsize=8)
        condition(ax,['(a) 50 Hz–10 kHz','(b) 100–1000 Hz'][j])
    save(fig,'F11C_propeller_spectrum',manifest,'Propeller microphone 1 spectrum and blade-passing references','Appendix Q2',
        'run3.mat:spectrumPower,spectrumFrequency','Mic1 Hann energy-normalized power per 5 Hz bin; 50 blocks; 20 µPa; ENBW 7.5 Hz',
        'Pa-tagged estimate; calibration chain unverified. Guides 200/400 Hz are theory, not proof of components. Selected 2000/3000/6185 Hz; last mechanism unassigned.')
    fig,axs=layout(2,3)
    for col,f in enumerate(P['fineF']):
        L=level(P['fineMaps'][:,:,col]);iy,ix=np.unravel_index(np.argmax(L),L.shape)
        for row in [0,1]:
            ax=axs[row,col];profile=L[iy,:] if row==0 else L[:,ix];idx=ix if row==0 else iy
            ax.plot(x,profile,color=BLUE);ax.plot(x[idx],profile[idx],'o',ms=4,color=OCHRE)
            ax.set(xlim=(-.5,.5),ylim=(-25,1),xticks=[-.4,0,.4],yticks=[-20,-10,0],xlabel=('x (m)' if row==0 else 'y (m)'))
            if col==0:ax.set_ylabel('Relative beamforming level (dB)')
            condition(ax,f"({chr(97+row*3+col)}) f = {f/1000:.3f} kHz, {'x' if row==0 else 'y'}-profile")
    save(fig,'F11D_propeller_profiles',manifest,'Propeller dominant-peak cross sections','Appendix Q2',
         'run3.mat:fineMaps','x at peak y; y at peak x; 601 grid; actual coordinates preserved',
         'Main-lobe spans in T11 use flanking local minima; N/A when truncated. Map width is not propeller diameter.')
    fig,axs=layout(1,1,120,95);ax=axs[0,0];fig.subplots_adjust(left=.2,bottom=.2,right=.96,top=.94)
    pos=tables['T11_position_accuracy']
    for src,color,marker in [(1,BLUE,'o'),(2,OCHRE,'s')]:
        z=[p for p in pos if p['source']==src]
        ax.plot(np.arange(len(z)),[p['error_m']*1000 for p in z],marker=marker,color=color,ls='none',label=f'Source {src}')
    ax.set(xticks=[0,1,2],xticklabels=['0.40 m\n3 kHz','0.40 m\n6 kHz','0.20 m\n6 kHz'],ylabel='Nominal position error (mm)',ylim=(0,55),xlim=(-.4,2.4))
    ax.legend(loc='upper right',fontsize=8,ncol=2)
    save(fig,'F11E_nominal_error',manifest,'Nominal source-position error','Appendix Q3 — human-review candidate',
         'T11_position_accuracy.csv','Euclidean distance to user-confirmed nominal source coordinates',
         'Optional inclusion only; nominal coordinates and physical channel wiring limit accuracy interpretation; no errors for unresolved maps.')
    for name in tables:
        manifest.append(dict(filename=name+'.csv',english_title=name.replace('_',' '),placement='Section 7' if name.startswith('T7') else 'Appendix Q1–Q3',
            source='data_frozen/*.mat',method_parameters_units='Recomputed by analyse_release.py; units in column names',
            limitations_review='Single recording; unverified calibration/wiring; empty cells indicate scientifically undefined metrics'))
    csvwrite(ROOT/'figures_tables'/'FIGURE_TABLE_MANIFEST.csv',manifest)
    from review_overview import make_overview
    make_overview(ROOT,manifest[:7])
    print('Recreated seven Section 7/11 figures and all tables from frozen data.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--from-frozen',action='store_true',required=True);parser.parse_args()
    render()
