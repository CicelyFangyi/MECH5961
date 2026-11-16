"""Recompute all report-facing measurements from frozen MATLAB arrays."""
from pathlib import Path
import csv,json
import numpy as np
from scipy.io import loadmat
from peak_structure import analyse,same_pair,sidelobe_conflict,minima_span

def csvwrite(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def process(root):
    data=root/'data_frozen';out=root/'figures_tables';qa=root/'qa'
    G=loadmat(data/'geometry.mat',squeeze_me=True)
    runs=[loadmat(data/f'run{r}.mat',squeeze_me=True) for r in range(4)]
    x,xb,xf=G['x'],G['xb'],G['xf'];freq=G['frequencies']
    control=[analyse(runs[0]['meanMaps'][:,:,k],x) for k in range(len(freq))]
    allrows=[];summaries=[];positions=[];details={};fine={};onsets=[];candidates=[]
    for run,ref in [(2,[-.2,.2]),(1,[-.2,0])]:
        D=runs[run];evidence=[]
        for k,f in enumerate(freq):
            A=analyse(D['meanMaps'][:,:,k],x,ref)
            blocks=[analyse(D['macroMaps'][:,:,k,b],xb,ref) for b in range(5)]
            loo=[analyse(np.mean(np.delete(D['macroMaps'][:,:,k,:],b,axis=2),axis=2),xb,ref) for b in range(5)]
            splits=[analyse(np.mean(D['macroMaps'][:,:,k,ids],axis=2),xb,ref) for ids in [slice(0,2),slice(2,5)]]
            control_conflict=sidelobe_conflict(A,control[k],x[1]-x[0])
            block_ok=[same_pair(A,b) for b in blocks];loo_ok=[same_pair(A,b) for b in loo];split_ok=[same_pair(A,b) for b in splits]
            # Independent identity-tolerance sensitivity, no dB threshold.
            sensitivity={str(t):all(same_pair(A,b,t) for b in blocks+loo+splits) for t in [.03,.04,.05]}
            temporal=all(block_ok+loo_ok+split_ok)
            robust=A['pair'] is not None and temporal and not control_conflict
            state=2 if robust else (1 if A['pair'] is not None else 0)
            row=dict(run=run,separation_m=ref[1]-ref[0],frequency_Hz=int(f),
                     detected_interior_maxima=sum(not z['boundary'] for z in A['peaks']),
                     mean_pair=A['pair'] is not None,saddle_depth_dB=A['saddle_depth_dB'],
                     matched_macroblocks=sum(block_ok),macroblocks=5,leave_one_out_pass=all(loo_ok),
                     split_pass=all(split_ok),single_source_conflict=control_conflict,
                     temporal_support=temporal,evidence_stage=state,
                     identity_30mm_pass=sensitivity['0.03'],identity_50mm_pass=sensitivity['0.05'])
            allrows.append(row);evidence.append(row)
            details[f'{run}_{int(f)}']=dict(mean=A,macroblocks=blocks,leave_one_out=loo,splits=splits,control=control[k])
            for rank,p in enumerate(A['peaks']):
                candidates.append(dict(run=run,frequency_Hz=int(f),rank=rank+1,**p,
                    association='PRIMARY_PAIR' if A['pair'] and p in A['pair'] else 'NOT_ASSOCIATED',
                    accepted=bool(robust and A['pair'] and p in A['pair'])))
            if f in [1000,3000,6000]:
                j=int(np.flatnonzero(D['fineF']==f)[0]);B=analyse(D['fineMaps'][:,:,j],xf,ref)
                grid=same_pair(A,B) if A['pair'] is not None else B['pair'] is None
                supported=robust and grid
                fine[(run,int(f))]=dict(analysis=B,supported=supported,power=D['fineMaps'][:,:,j],ref=ref)
                s=dict(run=run,separation_m=ref[1]-ref[0],frequency_Hz=int(f),
                       status='ROBUST_DUAL_PEAK' if supported else 'UNRESOLVED',
                       saddle_depth_dB=B['saddle_depth_dB'] if supported else None,
                       grid_301_601_consistent=grid,matched_macroblocks=sum(block_ok),
                       reason='All deterministic tests pass' if supported else 'No pair or insufficient temporal/control/grid evidence')
                if supported:
                    pp=D['fineMaps'][:,:,j]
                    s['peak_separation_m']=float(np.hypot(B['pair'][1]['x_m']-B['pair'][0]['x_m'],B['pair'][1]['y_m']-B['pair'][0]['y_m']))
                    for q,p in enumerate(B['pair']):
                        ix=int(np.argmin(abs(xf-p['x_m'])));iy=int(np.argmin(abs(xf-p['y_m'])))
                        s[f'x{q+1}_minima_span_m']=minima_span(xf,pp[iy,:],ix)
                        s[f'y{q+1}_minima_span_m']=minima_span(xf,pp[:,ix],iy)
                        err=np.hypot(p['x_m']-ref[q],p['y_m'])
                        s.update({f'x{q+1}_m':p['x_m'],f'y{q+1}_m':p['y_m'],f'error{q+1}_m':err})
                        positions.append(dict(run=run,separation_m=ref[1]-ref[0],frequency_Hz=int(f),source=q+1,
                            nominal_x_m=ref[q],nominal_y_m=0,x_m=p['x_m'],y_m=p['y_m'],error_m=err,
                            interpretation='Nominal geometry; physical channel wiring unverified'))
                summaries.append(s)
        sustained=[k for k in range(len(evidence)-2) if all(z['evidence_stage']==2 for z in evidence[k:k+3])]
        onset=None if not sustained else int(freq[sustained[0]])
        # Full-grid and finer frequency-neighborhood checks intentionally required before onset claims.
        onsets.append(dict(run=run,separation_m=ref[1]-ref[0],f_onset_Hz='indeterminate',
            first_coarse_sustained_candidate_Hz=onset,scan_step_Hz=250,FFT_bin_Hz=5,
            consecutive_frequencies_required=3,
            later_coarse_losses_Hz=';'.join(str(int(z['frequency_Hz'])) for z in evidence if onset and z['frequency_Hz']>onset and z['evidence_stage']<2),
            reason='Coarse candidate only: dense 601-grid neighborhood and time-block independence not established; no physical onset claimed'))
    prop=[];P=runs[3];sp=10*np.log10(P['spectrumPower']/(20e-6)**2)
    # Reconfirm the existing predeclared selection, not a new acoustic resolution cutoff.
    sf=P['spectrumFrequency'].astype(float);peaks=np.flatnonzero((sp[1:-1]>sp[:-2])&(sp[1:-1]>=sp[2:]))+1
    prom=np.array([sp[i]-np.median(sp[(abs(sf-sf[i])<=50)&(abs(sf-sf[i])>10)]) for i in peaks])
    selected=[]
    for lo,hi in [(1000,2000),(2500,4000),(4500,6500)]:
        ids=peaks[(sf[peaks]>=lo)&(sf[peaks]<=hi)&(prom>=6)]
        if not len(ids):ids=np.flatnonzero((sf>=lo)&(sf<=hi))
        selected.append(sf[ids[np.argmax(sp[ids])]])
    assert np.array_equal(selected,P['fineF'])
    for k,f in enumerate(P['fineF']):
        p=P['fineMaps'][:,:,k];iy,ix=np.unravel_index(np.argmax(p),p.shape)
        prop.append(dict(frequency_Hz=int(f),SPL_estimate_dB_re20uPa=float(sp[np.argmin(abs(sf-f))]),
            peak_x_m=xf[ix],peak_y_m=xf[iy],x_minima_span_m=minima_span(xf,p[iy,:],ix),
            y_minima_span_m=minima_span(xf,p[:,ix],iy),
            selection='Prominent spectral feature' if k<2 else 'Strongest band bin; mechanism unassigned',
            calibration='Pa metadata; sensitivity calibration chain unverified'))
    tables={'T11_resolution_summary':summaries,'T11_onset_summary':onsets,'T11_position_accuracy':positions,
            'T11_propeller_frequency_summary':prop,'T11_frequency_evidence':allrows,'T11_peak_candidates':candidates}
    xy=G['xy'];diam=np.max(np.linalg.norm(xy[:,None]-xy[None,:],axis=2))
    tables['T7_acquisition_geometry']=[dict(run=r,condition=['Single loudspeaker','Two loudspeakers','Two loudspeakers','Two-bladed propeller'][r],
        channels=56,sample_count=500000,nominal_fs_Hz=50000,metadata_fs_Hz=50000.0000101864,duration_s=10,
        source_x_m=['unknown','-0.20;0.00','-0.20;0.20','unknown'][r],h_m=1 if r in [1,2] else 'not independently recorded',
        array_Dx_m=np.ptp(xy[:,0]),array_Dy_m=np.ptp(xy[:,1]),array_Dmax_m=diam,
        blades=2 if r==3 else '',rotation_rev_s=100 if r==3 else '',BPF_Hz=200 if r==3 else '',
        calibration='Unverified chain; stored Pascals',channel_order='Sound Pressure_0..55 to coordinate rows 1..56; physical wiring unverified') for r in range(4)]
    for name,rows in tables.items():csvwrite(out/(name+'.csv'),rows)
    (qa/'peak_evidence.json').write_text(json.dumps(details,indent=2,allow_nan=False),encoding='utf8')
    return G,runs,fine,tables
