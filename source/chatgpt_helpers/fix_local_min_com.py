from pathlib import Path
p=Path(r"E:\chatgpt\analyze_dinp_40_50_local.py")
s=p.read_text(encoding='utf-8')
old="commind=float(bestmin[0])"
new="rotmin,tranmin=k(ca.positions.copy().astype(float),rca_min); bestmin=None\n for sh in itertools.product([-1,0,1],repeat=3):\n  image_m=lu+np.asarray(sh,dtype=float)*box; fit_m=image_m@rotmin+tranmin; score_m=float(np.linalg.norm(fit_m.mean(0)-rcom_min))\n  if bestmin is None or score_m<bestmin[0]: bestmin=(score_m,fit_m)\n commind=float(bestmin[0])"
if old not in s: raise SystemExit('old not found')
# Remove the prior bestmin loop, leaving the corrected one.
prior="bestmin=None\n for sh in itertools.product([-1,0,1],repeat=3):\n  image_m=lu+np.asarray(sh,dtype=float)*box; fit_m=image_m@rotp+tranp; score_m=float(np.linalg.norm(fit_m.mean(0)-rcom_min))\n  if bestmin is None or score_m<bestmin[0]: bestmin=(score_m,fit_m)\n "
s=s.replace(prior,'')
s=s.replace(old,new)
p.write_text(s,encoding='utf-8')
