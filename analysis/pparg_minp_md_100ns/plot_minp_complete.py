from pathlib import Path
import csv, math
from PIL import Image, ImageDraw, ImageFont

RUN = Path(r"E:\chatgpt\pparg_minp_md\run_20260915_seed20260917")
FONT = ImageFont.load_default()

def read_csv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

def panel(draw, box, title, series, xmin, xmax, ymin=None, ymax=None, hline=None):
    x0, y0, x1, y1 = box
    draw.rectangle(box, outline=(75, 75, 75), width=1)
    draw.text((x0 + 8, y0 + 6), title, fill=(20, 20, 20), font=FONT)
    all_y = [float(v) for _, rows, key, _ in series for row in rows for v in [row[key]] if row.get(key) not in (None, "")]
    if not all_y:
        return
    lo = min(all_y) if ymin is None else ymin
    hi = max(all_y) if ymax is None else ymax
    if hi <= lo: hi = lo + 1
    pad = 0.05 * (hi - lo)
    lo -= pad; hi += pad
    def xy(x, y):
        px = x0 + 48 + (x - xmin) / max(1e-9, xmax - xmin) * (x1 - x0 - 62)
        py = y1 - 28 - (y - lo) / (hi - lo) * (y1 - y0 - 55)
        return int(px), int(py)
    draw.text((x0 + 5, y1 - 21), f"{xmin:g}", fill=(90, 90, 90), font=FONT)
    draw.text((x1 - 35, y1 - 21), f"{xmax:g}", fill=(90, 90, 90), font=FONT)
    draw.text((x0 + 5, y0 + 24), f"{hi:.3g}", fill=(90, 90, 90), font=FONT)
    draw.text((x0 + 5, y1 - 43), f"{lo:.3g}", fill=(90, 90, 90), font=FONT)
    if hline is not None:
        yy = xy(xmin, hline)[1]
        draw.line((x0 + 48, yy, x1 - 14, yy), fill=(190, 80, 80), width=1)
    lx = x0 + 65
    for label, rows, key, color in series:
        step = max(1, len(rows) // 900)
        pts = [xy(float(r["time_ns"] if "time_ns" in r else r["bin_start_ns"]), float(r[key])) for r in rows[::step] if r.get(key) not in (None, "")]
        if len(pts) > 1: draw.line(pts, fill=color, width=2)
        for p in pts[::max(1, len(pts)//80)]: draw.ellipse((p[0]-2, p[1]-2, p[0]+2, p[1]+2), fill=color)
        draw.line((lx, y0 + 25, lx + 18, y0 + 25), fill=color, width=3)
        draw.text((lx + 22, y0 + 18), label, fill=(25, 25, 25), font=FONT)
        lx += 95 + len(label) * 2

def main():
    metrics = read_csv(RUN / "trajectory_metrics_0_100ns_complete.csv")
    contacts = read_csv(RUN / "contacts_5ns_0_100ns.csv")
    image = Image.new("RGB", (1400, 920), "white")
    draw = ImageDraw.Draw(image)
    draw.text((30, 18), "MiNP-PPARG MD complete 0-100 ns", fill=(10, 10, 10), font=FONT)
    panel(draw, (30, 55, 680, 455), "Ligand RMSD (A)", [("Whole", metrics, "whole_ligand_rmsd_A", (70,70,70)), ("Core", metrics, "core_rmsd_A", (190,50,45))], 0, 100, 0, 12)
    panel(draw, (720, 55, 1370, 455), "COM displacement (A)", [("Whole COM", metrics, "whole_com_displacement_A", (110,110,110)), ("Core COM", metrics, "core_com_displacement_A", (140,60,160))], 0, 100, 0, 12)
    panel(draw, (30, 485, 680, 885), "Minimum distance to original pocket (A)", [("Core-pocket", metrics, "pocket_min_distance_A", (35,125,80))], 0, 100, 0, 8, hline=4.5)
    panel(draw, (720, 485, 1370, 885), "Protein / pocket RMSD (A)", [("Protein backbone", metrics, "protein_backbone_rmsd_A", (40,70,130)), ("Pocket CA", metrics, "pocket_ca_rmsd_A", (210,110,20))], 0, 100, 0, 4)
    image.save(RUN / "minp_0_100ns_curves.png")

    image = Image.new("RGB", (1400, 520), "white")
    draw = ImageDraw.Draw(image)
    draw.text((30, 18), "MiNP-PPARG contacts, 5 ns bins", fill=(10,10,10), font=FONT)
    panel(draw, (30, 55, 680, 485), "Contact residues / atom pairs", [("Residues", contacts, "total_contact_residues_mean", (50,80,150)), ("Hydrophobic residues", contacts, "hydrophobic_contact_residues_mean", (35,145,80)), ("Atom pairs", contacts, "total_contact_pairs_mean", (110,110,110))], 0, 100, 0, 130)
    panel(draw, (720, 55, 1370, 485), "Hydrogen bonds / initial retention", [("H-bonds", contacts, "hbond_count_mean", (170,60,150)), ("Initial retention", contacts, "initial_contact_retention_mean", (210,100,0))], 0, 100, 0, 1)
    image.save(RUN / "minp_contacts_5ns.png")

    # 2D projection of the four representative ligand poses from the overlay PDB.
    models = {}
    current = None
    with (RUN / "representative_initial_rearranged_final_cluster_overlay.pdb").open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("MODEL"):
                current = line[10:].strip()
                models[current] = []
            elif current and line.startswith("HETATM") and line[17:20].strip() == "UNK" and line[76:78].strip().upper() != "H":
                models[current].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    colors = [(30,30,30), (210,60,50), (45,110,185), (35,150,85)]
    image = Image.new("RGB", (1400, 650), "white")
    draw = ImageDraw.Draw(image)
    draw.text((30, 18), "MiNP representative pose overlay (heavy atoms)", fill=(10,10,10), font=FONT)
    for panel_idx, (title, axis) in enumerate([("XY projection", (0,1)), ("XZ projection", (0,2))]):
        x0, y0, x1, y1 = (40, 70, 680, 610) if panel_idx == 0 else (720, 70, 1360, 610)
        draw.rectangle((x0,y0,x1,y1), outline=(80,80,80), width=1)
        pts_all = [p for coords in models.values() for p in coords]
        xs = [p[axis[0]] for p in pts_all]; ys = [p[axis[1]] for p in pts_all]
        lo_x, hi_x = min(xs)-2, max(xs)+2; lo_y, hi_y = min(ys)-2, max(ys)+2
        def xy(p):
            return (int(x0+35+(p[axis[0]]-lo_x)/(hi_x-lo_x)*(x1-x0-55)), int(y1-35-(p[axis[1]]-lo_y)/(hi_y-lo_y)*(y1-y0-55)))
        draw.text((x0+8,y0+8), title, fill=(20,20,20), font=FONT)
        lx=x0+80
        for (name, coords), color in zip(models.items(), colors):
            for p in coords: 
                q=xy(p); draw.ellipse((q[0]-3,q[1]-3,q[0]+3,q[1]+3), fill=color)
            draw.line((lx,y0+28,lx+18,y0+28), fill=color, width=3); draw.text((lx+22,y0+21),name[:18],fill=(20,20,20),font=FONT); lx += 125
    image.save(RUN / "minp_representative_overlay_projection.png")
    print("wrote plots")

if __name__ == "__main__":
    main()
