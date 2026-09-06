from pathlib import Path
p=Path(r"E:\chatgpt\analyze_dinp_40_50_local.py")
s=p.read_text(encoding='utf-8')
old="image,fit=best[1],best[2]; core=image[core_i]; b1=image[b1_i]; b2=image[b2_i]; core_pos=rms(fit[core_i],core40); rc,tr=k(core,core40); core_int=rms(core@rc+tr,core40); b1r=rms(b1@rc+tr,b140); b2r=rms(b2@rc+tr,b240); com40d=float(np.linalg.norm(fit.mean(0)-com40)); commind=float(np.linalg.norm(fit.mean(0)-rcom_min))"
new="image,fit=best[1],best[2]; core=image[core_i]; b1=image[b1_i]; b2=image[b2_i]; core_pos=rms(fit[core_i],core40); rc,tr=k(core,core40); core_int=rms(core@rc+tr,core40); b1r=rms(b1@rc+tr,b140); b2r=rms(b2@rc+tr,b240); com40d=float(np.linalg.norm(fit.mean(0)-com40)); bestmin=None\n for sh in itertools.product([-1,0,1],repeat=3):\n  image_m=lu+np.asarray(sh,dtype=float)*box; fit_m=image_m@rotp+tranp; score_m=float(np.linalg.norm(fit_m.mean(0)-rcom_min))\n  if bestmin is None or score_m<bestmin[0]: bestmin=(score_m,fit_m)\n commind=float(bestmin[0])"
if old not in s: raise SystemExit('old not found')
p.write_text(s.replace(old,new),encoding='utf-8')
