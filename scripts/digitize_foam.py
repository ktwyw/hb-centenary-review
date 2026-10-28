"""Extract plotted marker coordinates, not raw observations, from Ovarlez et al.

Source: arXiv:1007.4133v2, PDF page 5, Fig. 3b (SLES foam).
Run with an explicit PDF path. The published PDF is not redistributed.
Reproduction of the review's figures only requires the supplied derived CSV.
"""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def extract(pdf_path):
    import fitz
    digest = hashlib.sha256(Path(pdf_path).read_bytes()).hexdigest()
    expected = '8bf8cdb7d69100aacb1490bc989c09eb5b91dd93dc813bed2d75cef0d0d6138a'
    if digest != expected:
        raise ValueError('Source hash differs from the calibrated arXiv v2 PDF; extraction stopped.')
    document = fitz.open(pdf_path)
    page = document[4]
    drawings = page.get_drawings()
    # PDF-point coordinates calibrated to major tick lines: g=0.01 and 10 /s;
    # tau=0 and 100 Pa. These are page coordinates, not raster-pixel guesses.
    x0, x3 = 362.3312683105469, 475.0485534667969
    y0, y100 = 326.6993408203125, 224.02581787109375
    rows=[]
    for index,d in enumerate(drawings):
        r=d['rect']
        if d['type']!='f' or d['fill']!=(1.,1.,1.): continue
        if not (360<r.x0<501 and 225<r.y0<310 and 3<r.width<6 and 3<r.height<6): continue
        items=d['items']; center=np.array([(r.x0+r.x1)/2,(r.y0+r.y1)/2])
        if items[0][0]=='re': group='square'
        elif any(t[0]=='c' for t in items): group='circle'
        else:
            vertices=[]
            for t in items:
                if t[0]=='l':
                    for pt in t[1:3]:
                        a=np.array([pt.x,pt.y])
                        if not any(np.linalg.norm(a-b)<0.01 for b in vertices): vertices.append(a)
            if len(vertices)==3:
                center=np.mean(vertices,axis=0)
                far=max(vertices,key=lambda a:np.linalg.norm(a-center))
                if r.width>r.height:
                    # Two vertices share the base y coordinate; the remaining
                    # vertex fixes orientation even for rounded PDF values.
                    ys=sorted(v[1] for v in vertices)
                    group='triangle_up' if ys[1]-ys[0]>ys[2]-ys[1] else 'triangle_down'
                else:
                    xs=sorted(v[0] for v in vertices)
                    group='triangle_left' if xs[1]-xs[0]>xs[2]-xs[1] else 'triangle_right'
            elif len(vertices)==4: group='diamond'
            else: raise ValueError(f'Unknown marker {index}: {len(vertices)} vertices')
        g=10**(-2+3*(center[0]-x0)/(x3-x0))
        tau=100*(y0-center[1])/(y0-y100)
        rows.append(dict(path_index=index,group=group,x_pt=center[0],y_pt=center[1],shear_rate_s_inv=g,stress_Pa=tau))
    rows.sort(key=lambda q:(q['group'],q['shear_rate_s_inv']))
    assert len(rows)==152, 'Unexpected marker count'
    assert len({r['path_index'] for r in rows})==152
    dest=ROOT/'data/foam_sles_fig3b_digitized.csv'
    with dest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    meta={'source':'Ovarlez, Krishan and Cohen-Addad (2010), EPL 91, 68005',
          'doi':'10.1209/0295-5075/91/68005','inspected_version':'arXiv:1007.4133v2',
          'source_url':'https://arxiv.org/pdf/1007.4133v2','pdf_page':5,'figure':'3b',
          'sha256':hashlib.sha256(Path(pdf_path).read_bytes()).hexdigest(),
          'count':len(rows),'group_counts':{g:sum(r['group']==g for r in rows) for g in sorted({r['group'] for r in rows})},
          'axes':{'x_at_0.01':x0,'x_at_10':x3,'y_at_0_Pa':y0,'y_at_100_Pa':y100},
          'status':'Coordinates of published vector markers; not original raw MRI/torque data.',
          'group_meaning':'Marker shape identifies an acquisition series; speeds are not assigned because Fig. 3 has no speed legend.',
          'uncertainty':'No experimental covariance or error bars recoverable. Coordinate perturbation is only a digitisation sensitivity.'}
    (ROOT/'data/foam_digitization_provenance.json').write_text(json.dumps(meta,indent=2))
    broad=[]
    for index,d in enumerate(drawings):
        r=d['rect']
        if d['type']=='f' and d['fill']==(1.,1.,1.) and 350<r.x0<510 and 210<r.y0<335 and 3<r.width<6 and 3<r.height<6:
            broad.append(index)
    selected={row['path_index'] for row in rows}
    xy=np.array([[row['x_pt'],row['y_pt']] for row in rows])
    delta=np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=-1)
    np.fill_diagonal(delta,np.inf)
    first,second=np.unravel_index(np.argmin(delta),delta.shape)
    audit={'source_sha256':digest,'panel_bounds_PDF_pt':[360,225,501,310],
      'broader_bounds_PDF_pt':[350,210,510,335],'marker_size_range_PDF_pt':[3,6],
      'selected':len(rows),'broader_candidates':len(broad),
      'excluded_broader_candidate_paths':sorted(set(broad)-selected),
      'unique_path_ids':len(selected),'unique_centres':len(np.unique(xy,axis=0)),
      'minimum_centre_distance_PDF_pt':float(delta[first,second]),
      'closest_paths':[rows[first]['path_index'],rows[second]['path_index']],
      'interpretation':'No points removed on proximity alone; different marker groups can overlap. Raw-data independence is unresolved.'}
    (ROOT/'data/foam_selection_audit.json').write_text(json.dumps(audit,indent=2))
    for row in rows:
        page.draw_circle(fitz.Point(row['x_pt'],row['y_pt']),radius=1.2,color=(1,0,0),width=.35,overlay=True)
    page.get_pixmap(matrix=fitz.Matrix(3,3),clip=fitz.Rect(325,200,530,347),alpha=False).save(ROOT/'data/foam_extraction_overlay.png')
    document.close()
    print(json.dumps({'provenance':meta,'selection_audit':audit},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('pdf',type=Path);extract(p.parse_args().pdf)
